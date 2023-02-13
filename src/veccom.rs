use crate::common::prefix_decompose;
use crate::traits::{Block, LeafExpander, LengthDoubling};
use core::marker::PhantomData;
use digest::{Digest, Output as DigestOutput, XofReader};
use rand::thread_rng;

pub trait VecCom {
    type Commitment: Clone + AsRef<[u8]>;
    type DecommitmentKey: Copy;
    type Decommitment: Clone;
    type XofReader: XofReader;
    fn commit(
        log_num_messages: u32,
    ) -> (
        Self::Commitment,
        Self::DecommitmentKey,
        Vec<Self::XofReader>,
    );
    fn decommit(
        log_num_messages: u32,
        decommitment_key: Self::DecommitmentKey,
        index: usize,
    ) -> Self::Decommitment;
    fn verify(
        log_num_messages: u32,
        commitment: &Self::Commitment,
        decommitment: &Self::Decommitment,
        index: usize,
    ) -> Option<Vec<Self::XofReader>>;
}

#[derive(Debug, PartialEq, Eq)]
pub struct Commitment<H: Digest> {
    hash: DigestOutput<H>,
}

impl<H: Digest> From<DigestOutput<H>> for Commitment<H> {
    fn from(hash: DigestOutput<H>) -> Self {
        Self { hash }
    }
}

impl<H: Digest> Into<DigestOutput<H>> for &Commitment<H> {
    fn into(self) -> DigestOutput<H> {
        self.hash.clone()
    }
}

impl<H: Digest> Into<DigestOutput<H>> for Commitment<H> {
    fn into(self) -> DigestOutput<H> {
        self.hash
    }
}

impl<H: Digest> Clone for Commitment<H> {
    fn clone(&self) -> Self {
        Self {
            hash: self.hash.clone(),
        }
    }
}

impl<H: Digest> AsRef<[u8]> for Commitment<H> {
    fn as_ref(&self) -> &[u8] {
        &self.hash
    }
}

impl<H: Digest> bincode::Encode for Commitment<H> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.hash.as_slice(), encoder)?;
        Ok(())
    }
}

pub struct GgmVecCom<B, Prg, H, LE>
where
    B: Block,
    Prg: LengthDoubling<Block = B>,
    H: Digest,
    LE: LeafExpander,
{
    _phantom_b: PhantomData<B>,
    _phantom_prg: PhantomData<Prg>,
    _phantom_h: PhantomData<H>,
    _phantom_le: PhantomData<LE>,
}

impl<B, Prg, H, LE> VecCom for GgmVecCom<B, Prg, H, LE>
where
    B: Block + bincode::Encode,
    Prg: LengthDoubling<Block = B>,
    H: Digest,
    LE: LeafExpander<Com = B>,
{
    type Commitment = Commitment<H>;
    type DecommitmentKey = B;
    type Decommitment = (Vec<B>, B);
    type XofReader = LE::XofR;

    fn commit(
        log_num_messages: u32,
    ) -> (
        Self::Commitment,
        Self::DecommitmentKey,
        Vec<Self::XofReader>,
    ) {
        let num_messages = 1usize << log_num_messages;
        let tree_height = log_num_messages;

        let decommitment_key = B::random(&mut thread_rng());
        let mut seeds = vec![B::default(); num_messages];
        seeds[0] = decommitment_key;
        for level_l in 0..tree_height {
            for j in (0..1 << level_l).rev() {
                let (k0, k1) = Prg::expand(seeds[j]);
                seeds[2 * j] = k0;
                seeds[2 * j + 1] = k1;
            }
        }
        let mut xofs = Vec::<LE::XofR>::with_capacity(num_messages);
        let commitment = {
            let mut hasher = H::new();
            seeds.into_iter().for_each(|s_j| {
                let (com_j, xof_j) = LE::expand_leaf(s_j);
                hasher.update(com_j);
                xofs.push(xof_j);
            });
            hasher.finalize()
        }
        .into();
        (commitment, decommitment_key, xofs)
    }

    fn decommit(log_num_messages: u32, decommitment_key: B, index: usize) -> Self::Decommitment {
        let num_messages = 1usize << log_num_messages;
        let tree_height = log_num_messages as usize;
        assert!(index < num_messages);
        let index_prefixes = prefix_decompose(index, log_num_messages);

        let mut seeds = vec![B::default(); num_messages];
        let mut decommitment = Vec::with_capacity(tree_height);
        seeds[0] = decommitment_key;
        (seeds[0], seeds[1]) = Prg::expand(seeds[0]);
        for level_l in 1..tree_height {
            decommitment.push(seeds[index_prefixes[level_l - 1] ^ 1]);
            let i = index_prefixes[level_l - 1];
            let (k0, k1) = Prg::expand(seeds[i]);
            seeds[2 * i] = k0;
            seeds[2 * i + 1] = k1;
        }
        decommitment.push(seeds[index_prefixes[tree_height - 1] ^ 1]);
        let (com_i, _) = LE::expand_leaf(seeds[index]);
        (decommitment, com_i)
    }

    fn verify(
        log_num_messages: u32,
        commitment: &Self::Commitment,
        decommitment: &Self::Decommitment,
        index: usize,
    ) -> Option<Vec<Self::XofReader>> {
        let num_messages = 1usize << log_num_messages;
        let tree_height = log_num_messages as usize;
        assert!(index < num_messages);
        let index_prefixes = prefix_decompose(index, log_num_messages);

        let (decommitment, com_i) = decommitment;

        let mut seeds = vec![B::default(); num_messages];
        for level_l in 1..tree_height {
            let i = index_prefixes[level_l - 1];
            seeds[i ^ 1] = decommitment[level_l - 1];
            for j in (0..1 << level_l).rev() {
                if j == i {
                    continue;
                }
                let (k0, k1) = Prg::expand(seeds[j]);
                seeds[2 * j] = k0;
                seeds[2 * j + 1] = k1;
            }
        }
        seeds[index ^ 1] = decommitment[tree_height - 1];
        let mut xofs = Vec::<LE::XofR>::with_capacity(num_messages);
        let recomputed_commitment = {
            let mut hasher = H::new();
            seeds.into_iter().enumerate().for_each(|(j, s_j)| {
                let (com_j, xof_j) = LE::expand_leaf(s_j);
                if j == index {
                    hasher.update(com_i);
                    // this XofReader is not supposed to be used, but we need to
                    // put something here
                    xofs.push(xof_j);
                } else {
                    hasher.update(com_j);
                    xofs.push(xof_j);
                }
            });
            hasher.finalize()
        };

        if recomputed_commitment == commitment.into() {
            Some(xofs)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Block128;
    use crate::primitives::{Aes128CtrLdPrg, Blake3LE};
    use blake3;

    #[test]
    fn test_veccom_correctness() {
        let log_num_messages = 6;
        let num_messages = 1 << log_num_messages;
        type Msg = [u8; 32];
        type VC = GgmVecCom<Block128, Aes128CtrLdPrg, blake3::Hasher, Blake3LE<Block128>>;

        let (commitment, decommitment_key, xofs_sender) = VC::commit(log_num_messages);
        let messages_sender: Vec<_> = xofs_sender
            .into_iter()
            .map(|mut xof| {
                let mut m = Msg::default();
                xof.read(&mut m);
                m
            })
            .collect();
        for index in 0..num_messages {
            let decommitment = VC::decommit(log_num_messages, decommitment_key, index);
            let xofs_receiver = VC::verify(log_num_messages, &commitment, &decommitment, index);
            assert!(xofs_receiver.is_some());
            let mut xofs_receiver = xofs_receiver.unwrap();
            for i in 0..num_messages {
                if i == index {
                    continue;
                }
                let mut message_receiver_i = Msg::default();
                xofs_receiver[i].read(&mut message_receiver_i);
                assert_eq!(message_receiver_i, messages_sender[i]);
            }
        }
    }
}
