use crate::faest::{Prover as FaestProver, Verifier as FaestVerifier};
use crate::{PublicKey, SecretKey};
use blake3;
use core::marker::PhantomData;
use digest::Digest;

pub trait Signer {
    type Signature;
    fn new(secret_key: SecretKey, public_key: PublicKey) -> Self;
    fn sign(self, message: &[u8]) -> Self::Signature;
}

pub trait SignatureVerifier {
    type Signature;
    fn new(public_key: PublicKey) -> Self;
    fn verify(self, signature: &Self::Signature, message: &[u8]) -> bool;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FsSignature<P>
where
    P: FaestProver,
    P::Commitment: bincode::Encode,
    P::Response: bincode::Encode,
    P::Proof: bincode::Encode,
    P::Decommitment: bincode::Encode,
{
    commitment: P::Commitment,
    response: P::Response,
    proof: P::Proof,
    decommitment: P::Decommitment,
}

impl<P> bincode::Encode for FsSignature<P>
where
    P: FaestProver,
    P::Commitment: bincode::Encode,
    P::Response: bincode::Encode,
    P::Proof: bincode::Encode,
    P::Decommitment: bincode::Encode,
{
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.commitment, encoder)?;
        bincode::Encode::encode(&self.response, encoder)?;
        bincode::Encode::encode(&self.proof, encoder)?;
        bincode::Encode::encode(&self.decommitment, encoder)?;
        Ok(())
    }
}

pub struct FsSigner<P: FaestProver, V: FaestVerifier> {
    public_key: PublicKey,
    prover: P,
    _phantom_v: PhantomData<V>,
}

impl<P, V> Signer for FsSigner<P, V>
where
    P: FaestProver,
    V: FaestVerifier<
        Commitment = P::Commitment,
        Challenge1 = P::Challenge1,
        Response = P::Response,
        Challenge2 = P::Challenge2,
        Proof = P::Proof,
        Choice = P::Choice,
        Decommitment = P::Decommitment,
    >,
    P::Commitment: bincode::Encode,
    P::Response: bincode::Encode,
    P::Proof: bincode::Encode,
    P::Decommitment: bincode::Encode,
{
    type Signature = FsSignature<P>;
    fn new(secret_key: SecretKey, public_key: PublicKey) -> Self {
        Self {
            public_key,
            prover: P::new(secret_key, public_key),
            _phantom_v: PhantomData,
        }
    }
    fn sign(mut self, message: &[u8]) -> Self::Signature {
        let bincode_cfg = bincode::config::standard()
            .with_little_endian()
            .with_fixed_int_encoding()
            .skip_fixed_array_length();

        let commitment = self.prover.commit();
        let h1 = {
            let mut hasher = blake3::Hasher::new_derive_key("FS-FAEST-1");
            bincode::encode_into_std_write(self.public_key, &mut hasher, bincode_cfg).unwrap();
            hasher.update(&(message.len() as u64).to_le_bytes());
            hasher.update(message);
            bincode::encode_into_std_write(&commitment, &mut hasher, bincode_cfg).unwrap();
            hasher.finalize()
        };

        let challenge1 = V::generate_commit_challenge_from_seed(h1.into());
        let response = self.prover.commit_prove_consistency(challenge1);
        let h2 = {
            let mut hasher = blake3::Hasher::new_derive_key("FS-FAEST-2").chain_update(&h1);
            bincode::encode_into_std_write(&response, &mut hasher, bincode_cfg).unwrap();
            hasher.finalize()
        };

        let challenge2 = V::generate_challenge_from_seed(h2.into());
        let proof = self.prover.prove(challenge2);
        let h3 = {
            let mut hasher = blake3::Hasher::new_derive_key("FS-FAEST-3").chain_update(&h2);
            bincode::encode_into_std_write(&proof, &mut hasher, bincode_cfg).unwrap();
            hasher.finalize()
        };

        let choice = V::generate_choice_from_seed(h3.into());
        let decommitment = self.prover.transfer(choice);

        FsSignature {
            commitment,
            response,
            proof,
            decommitment,
        }
    }
}

pub struct FsVerifier<P: FaestProver, V: FaestVerifier> {
    public_key: PublicKey,
    verifier: V,
    _phantom_p: PhantomData<P>,
}

impl<P, V> SignatureVerifier for FsVerifier<P, V>
where
    P: FaestProver,
    V: FaestVerifier<
        Commitment = P::Commitment,
        Challenge1 = P::Challenge1,
        Response = P::Response,
        Challenge2 = P::Challenge2,
        Proof = P::Proof,
        Choice = P::Choice,
        Decommitment = P::Decommitment,
    >,
    P::Commitment: bincode::Encode,
    P::Response: bincode::Encode,
    P::Proof: bincode::Encode,
    P::Decommitment: bincode::Encode,
{
    type Signature = FsSignature<P>;
    fn new(public_key: PublicKey) -> Self {
        Self {
            public_key,
            verifier: V::new(public_key),
            _phantom_p: PhantomData,
        }
    }
    fn verify(mut self, signature: &Self::Signature, message: &[u8]) -> bool {
        let bincode_cfg = bincode::config::standard()
            .with_little_endian()
            .with_fixed_int_encoding()
            .skip_fixed_array_length();

        let Self::Signature {
            commitment,
            response,
            proof,
            decommitment,
        } = signature;

        let h1 = {
            let mut hasher = blake3::Hasher::new_derive_key("FS-FAEST-1");
            bincode::encode_into_std_write(self.public_key, &mut hasher, bincode_cfg).unwrap();
            hasher.update(&(message.len() as u64).to_le_bytes());
            hasher.update(message);
            bincode::encode_into_std_write(&commitment, &mut hasher, bincode_cfg).unwrap();
            hasher.finalize()
        };
        let _challenge1 = self
            .verifier
            .commit_send_challenge_from_seed(commitment.clone(), h1.into());

        self.verifier.commit_receive_response(response.clone());
        let h2 = {
            let mut hasher = blake3::Hasher::new_derive_key("FS-FAEST-2").chain_update(&h1);
            bincode::encode_into_std_write(&response, &mut hasher, bincode_cfg).unwrap();
            hasher.finalize()
        };
        let _challenge2 = self.verifier.send_challenge_from_seed(h2.into());

        let h3 = {
            let mut hasher = blake3::Hasher::new_derive_key("FS-FAEST-3").chain_update(&h2);
            bincode::encode_into_std_write(&proof, &mut hasher, bincode_cfg).unwrap();
            hasher.finalize()
        };
        let _choice = self.verifier.choose_from_seed(proof.clone(), h3.into());

        self.verifier.verify(decommitment.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Block128;
    use crate::faest::keygen;
    use crate::faest::{FaestProverFromHC, FaestVerifierFromHC};
    use crate::homcom::{HomCom128ReceiverFromVitH, HomCom128SenderFromVitH};
    use crate::primitives::{Aes128CtrLdPrg, Blake3LE};
    use crate::veccom::GgmVecCom;
    use crate::voleith::{VoleInTheHeadReceiverFromVC, VoleInTheHeadSenderFromVC};

    type VC = GgmVecCom<Block128, Aes128CtrLdPrg, blake3::Hasher, Blake3LE<Block128>>;
    type FP = FaestProverFromHC<HomCom128SenderFromVitH<VoleInTheHeadSenderFromVC<VC>>>;
    type FV = FaestVerifierFromHC<HomCom128ReceiverFromVitH<VoleInTheHeadReceiverFromVC<VC>>>;
    type FaestSigner = FsSigner<FP, FV>;
    type FaestVerifier = FsVerifier<FP, FV>;

    #[test]
    fn test_correctness() {
        let (sk, pk) = keygen();
        let message = "Am I happy or in misery?";
        let signer = FaestSigner::new(sk, pk);
        let signature = signer.sign(message.as_bytes());
        let verifier = FaestVerifier::new(pk);
        let result = verifier.verify(&signature, message.as_bytes());
        assert!(result);
        let bincode_cfg = bincode::config::standard()
            .with_little_endian()
            .with_fixed_int_encoding()
            .skip_fixed_array_length();
        let encoded_signature = bincode::encode_to_vec(signature, bincode_cfg).unwrap();
        assert_eq!(encoded_signature.len(), 7191);
    }
}
