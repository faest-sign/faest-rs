use crate::arithmetic::{bit_xor_assign, hash_bitvector_and_matrix, hash_matrix};
use crate::field::BitMulAccumulate;
use crate::gf2::{GF2Vector, GF2View};
use crate::gf2psmall::SmallGF;
use crate::veccom::VecCom;
use core::marker::PhantomData;
use core::mem;
use digest::XofReader;
use digest::{Digest, Output as DigestOutput};
use itertools::izip;
use ndarray::{s, Array1, Array2, ArrayView1, ArrayView2, Axis};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaChaRng;

type Vector<T> = Array1<T>;
type VectorView<'a, T> = ArrayView1<'a, T>;
type Matrix<T> = Array2<T>;
type MatrixView<'a, T> = ArrayView2<'a, T>;

pub trait VoleInTheHeadSender {
    type Commitment: Clone;
    type Challenge: Clone;
    type Response: Clone;
    type Decommitment: Clone;
    type Field: Clone;

    const FIELD_SIZE: usize;

    fn new(vole_length: usize, num_repetitions: usize) -> Self;
    fn commit_message(&mut self, message: GF2Vector) -> Self::Commitment;
    fn commit_random(&mut self) -> Self::Commitment;
    fn consistency_check_respond(&mut self, random_points: Self::Challenge) -> Self::Response;
    #[allow(non_snake_case)]
    fn decommit(&mut self, Deltas: Vector<Self::Field>) -> Self::Decommitment;
    fn get_output(&self) -> (&GF2View, MatrixView<'_, Self::Field>);
}

pub trait VoleInTheHeadReceiver {
    type Commitment: Clone;
    type Challenge: Clone;
    type Response: Clone;
    type Decommitment: Clone;
    type Field: Clone;

    const FIELD_SIZE: usize;

    fn new(vole_length: usize, num_repetitions: usize) -> Self;
    fn receive_commitment(&mut self, commitment: Self::Commitment);
    fn receive_random_commitment(&mut self, commitment: Self::Commitment);
    fn consistency_challenge_from_seed(num_repetitions: usize, seed: [u8; 32]) -> Self::Challenge;
    fn generate_consistency_challenge_from_seed(&mut self, seed: [u8; 32]) -> Self::Challenge;
    fn generate_consistency_challenge(&mut self) -> Self::Challenge;
    fn store_consistency_response(&mut self, consistency_response: Self::Response);
    fn final_challenge_from_seed(num_repetitions: usize, seed: [u8; 32]) -> Vector<Self::Field>;
    fn generate_final_challenge_from_seed(&mut self, seed: [u8; 32]) -> Vector<Self::Field>;
    fn generate_final_challenge(&mut self) -> Vector<Self::Field>;
    fn receive_decommitment(&mut self, decommitments: &Self::Decommitment) -> bool;
    fn get_output(&self) -> (VectorView<'_, Self::Field>, MatrixView<'_, Self::Field>);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VoleInTheHeadSenderState {
    New,
    Committed,
    RespondedToConsistencyChallenge,
    Ready,
}

#[allow(non_snake_case)]
pub struct VoleInTheHeadSenderFromVC<F: SmallGF, VC: VecCom, H: Digest = blake3::Hasher> {
    vole_length: usize,
    num_repetitions: usize,
    state: VoleInTheHeadSenderState,
    u: GF2Vector,
    V: Matrix<F>,
    decommitment_keys: Vec<VC::DecommitmentKey>,
    _phantom_vc: PhantomData<VC>,
    _phantom_h: PhantomData<H>,
}

impl<F: SmallGF, VC: VecCom, H: Digest> VoleInTheHeadSenderFromVC<F, VC, H> {
    fn commit_impl(
        &mut self,
        message: Option<GF2Vector>,
    ) -> <Self as VoleInTheHeadSender>::Commitment {
        assert_eq!(self.state, VoleInTheHeadSenderState::New);

        let log_q = F::LOG_ORDER;
        let tau = self.num_repetitions;
        let ell_hat = self.vole_length + self.num_repetitions;

        let mut hasher = H::new();
        let mut correction_values =
            Vec::with_capacity(if message.is_some() { tau } else { tau - 1 });

        // iteration i = 0
        {
            let (commitment, decommitment_key, mut xofs) = VC::commit(log_q);
            self.decommitment_keys.push(decommitment_key);
            hasher.update(commitment.as_ref());
            // iteration x = 0
            let u_0 = &mut self.u;
            xofs[0].read(u_0.bits.as_raw_mut_slice());
            let mut v_0 = self.V.row_mut(0);
            debug_assert_eq!(xofs.len(), F::ORDER);
            for (x, xof_x) in xofs.iter_mut().enumerate().skip(1) {
                let mut r_x_0 = GF2Vector::with_capacity(ell_hat);
                r_x_0.resize(ell_hat, false);
                xof_x.read(r_x_0.as_raw_mut_slice());
                bit_xor_assign(u_0, &r_x_0);
                BitMulAccumulate::bitmul_accumulate(
                    v_0.as_slice_mut().unwrap(),
                    F::from(x),
                    r_x_0.as_raw_slice(),
                );
            }
            if let Some(msg) = message {
                let mut msg_correction = msg.clone();
                bit_xor_assign(&mut msg_correction, u_0);
                msg_correction.bits.set_uninitialized(false);
                debug_assert_eq!(msg_correction.bits.len(), self.vole_length);
                correction_values.push(msg_correction);
                u_0.bits[0..self.vole_length].copy_from_bitslice(msg.as_ref());
            }
        }

        // other iterations
        for i in 1..tau {
            let (commitment, decommitment_key, mut xofs) = VC::commit(log_q);
            self.decommitment_keys.push(decommitment_key);
            hasher.update(commitment.as_ref());
            // iteration x = 0
            let mut u_i = {
                let mut r_0_1 = GF2Vector::with_capacity(ell_hat);
                r_0_1.resize(ell_hat, false);
                xofs[0].read(r_0_1.as_raw_mut_slice());
                r_0_1
            };
            let mut v_i = self.V.row_mut(i);
            debug_assert_eq!(xofs.len(), F::ORDER);
            for (x, xof_x) in xofs.iter_mut().enumerate().skip(1) {
                let mut r_x_i = GF2Vector::with_capacity(ell_hat);
                r_x_i.resize(ell_hat, false);
                xof_x.read(r_x_i.as_raw_mut_slice());
                bit_xor_assign(&mut u_i, &r_x_i);
                BitMulAccumulate::bitmul_accumulate(
                    v_i.as_slice_mut().unwrap(),
                    F::from(x),
                    r_x_i.as_raw_slice(),
                );
            }
            bit_xor_assign(&mut u_i, &self.u);
            correction_values.push(u_i);
        }

        let vc_commitment_hash = hasher.finalize();

        // transpose V such that we can now access it row-wise
        {
            // reversed_axes consumes self, but we cannot move out of the &mut ref.
            // Hence, we do the swapping with a dummy value.
            let mut tmp = Matrix::<F>::default((1, 1));
            mem::swap(&mut tmp, &mut self.V);
            tmp = tmp.reversed_axes();
            mem::swap(&mut tmp, &mut self.V);
            self.V = self.V.as_standard_layout().into_owned();
        }

        self.state = VoleInTheHeadSenderState::Committed;
        Commitment {
            vc_commitment_hash,
            correction_values,
        }
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Commitment<H: Digest> {
    pub vc_commitment_hash: DigestOutput<H>,
    pub correction_values: Vec<GF2Vector>,
}

impl<H: Digest> Clone for Commitment<H> {
    fn clone(&self) -> Self {
        Self {
            vc_commitment_hash: self.vc_commitment_hash.clone(),
            correction_values: self.correction_values.clone(),
        }
    }
}

impl<H: Digest> bincode::Encode for Commitment<H> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.vc_commitment_hash.as_slice(), encoder)?;
        bincode::Encode::encode(&self.correction_values, encoder)?;
        Ok(())
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Response<F, H: Digest> {
    pub vector: Vector<F>,
    pub hash: DigestOutput<H>,
}

impl<F: Clone, H: Digest> Clone for Response<F, H> {
    fn clone(&self) -> Self {
        Self {
            vector: self.vector.clone(),
            hash: self.hash.clone(),
        }
    }
}

impl<F: bincode::Encode, H: Digest> bincode::Encode for Response<F, H> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.hash.as_slice(), encoder)?;
        bincode::Encode::encode(&self.vector.as_slice(), encoder)?;
        Ok(())
    }
}

#[allow(non_snake_case)]
impl<F: SmallGF, VC: VecCom, H: Digest> VoleInTheHeadSender
    for VoleInTheHeadSenderFromVC<F, VC, H>
{
    type Commitment = Commitment<H>;
    type Challenge = Vector<Self::Field>;
    type Response = Response<Self::Field, H>;
    type Decommitment = Vec<VC::Decommitment>;
    type Field = F;

    const FIELD_SIZE: usize = 8;

    fn new(vole_length: usize, num_repetitions: usize) -> Self {
        let ell_hat = vole_length + num_repetitions;
        let mut output = Self {
            vole_length,
            num_repetitions,
            state: VoleInTheHeadSenderState::New,
            u: GF2Vector::new(),
            V: Matrix::<F>::zeros((num_repetitions, ell_hat)),
            decommitment_keys: Vec::with_capacity(num_repetitions),
            _phantom_vc: PhantomData,
            _phantom_h: PhantomData,
        };
        output.u.resize(ell_hat, false);
        output
    }

    fn commit_message(&mut self, message: GF2Vector) -> Self::Commitment {
        self.commit_impl(Some(message))
    }

    fn commit_random(&mut self) -> Self::Commitment {
        self.commit_impl(None)
    }

    fn consistency_check_respond(&mut self, random_points: Vector<Self::Field>) -> Self::Response {
        assert_eq!(self.state, VoleInTheHeadSenderState::Committed);
        assert_eq!(random_points.len(), self.num_repetitions);

        let (u_tilde, V_tilde) =
            hash_bitvector_and_matrix((&random_points).into(), self.u.as_ref(), (&self.V).into());
        let V_tilde_hash = H::digest(
            V_tilde
                .as_slice()
                .expect("not in standard memory order")
                .iter()
                .flat_map(|x| x.to_repr())
                .collect::<Vec<u8>>(),
        );

        // ignore the extra entries
        self.u.resize(self.vole_length, false);
        // cannot resize the matrix V, so keep in mind to ignore the last rows

        self.state = VoleInTheHeadSenderState::RespondedToConsistencyChallenge;
        Self::Response {
            vector: u_tilde,
            hash: V_tilde_hash,
        }
    }

    #[allow(non_snake_case)]
    fn decommit(&mut self, Deltas: Vector<F>) -> Self::Decommitment {
        assert_eq!(
            self.state,
            VoleInTheHeadSenderState::RespondedToConsistencyChallenge
        );
        assert_eq!(Deltas.len(), self.num_repetitions);
        let log_q = F::LOG_ORDER;
        let decommitments = Deltas
            .iter()
            .zip(self.decommitment_keys.iter())
            .map(|(Delta_i, &decommitment_key)| {
                VC::decommit(log_q, decommitment_key, (*Delta_i).into())
            })
            .collect();
        self.state = VoleInTheHeadSenderState::Ready;
        decommitments
    }

    fn get_output(&self) -> (&GF2View, MatrixView<'_, Self::Field>) {
        assert_ne!(self.state, VoleInTheHeadSenderState::New);
        let tau = self.num_repetitions;
        let ell = self.vole_length;
        assert_eq!(self.V.shape(), &[ell + tau, tau]);
        (&self.u.as_ref()[..ell], self.V.slice(s![..ell, ..]))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VoleInTheHeadReceiverState {
    New,
    CommitmentReceived,
    ConsistencyChallengeGenerated,
    ConsistencyChallengeResponseReceived,
    FinalChallengeGenerated,
    Failed,
    Ready,
}

#[allow(non_snake_case)]
pub struct VoleInTheHeadReceiverFromVC<F: SmallGF, VC: VecCom, H: Digest = blake3::Hasher> {
    vole_length: usize,
    num_repetitions: usize,
    state: VoleInTheHeadReceiverState,
    // commitment: Option<Blake3Hash>,
    random_commitment: bool,
    vc_commitment_hash: DigestOutput<H>,
    correction_values: Vec<GF2Vector>,
    consistency_challenges: Vector<F>,
    u_tilde: Vector<F>,
    V_tilde_hash: DigestOutput<H>,
    final_challenges: Vector<F>,
    W: Matrix<F>,
    _phantom_vc: PhantomData<VC>,
    _phantom_h: PhantomData<H>,
}

#[allow(non_snake_case)]
impl<F: SmallGF, VC: VecCom, H: Digest> VoleInTheHeadReceiver
    for VoleInTheHeadReceiverFromVC<F, VC, H>
{
    type Commitment = Commitment<H>;
    type Challenge = Vector<Self::Field>;
    type Response = Response<Self::Field, H>;
    type Decommitment = Vec<VC::Decommitment>;
    type Field = F;

    const FIELD_SIZE: usize = 8;

    fn new(vole_length: usize, num_repetitions: usize) -> Self {
        Self {
            vole_length,
            num_repetitions,
            state: VoleInTheHeadReceiverState::New,
            // commitment: None,
            random_commitment: false,
            vc_commitment_hash: Default::default(),
            correction_values: Vec::new(),
            consistency_challenges: Default::default(),
            u_tilde: Default::default(),
            V_tilde_hash: Default::default(),
            final_challenges: Default::default(),
            W: Default::default(),
            _phantom_vc: PhantomData,
            _phantom_h: PhantomData,
        }
    }

    fn receive_commitment(&mut self, commitment: Self::Commitment) {
        assert_eq!(self.state, VoleInTheHeadReceiverState::New);
        let Commitment {
            vc_commitment_hash,
            correction_values,
        } = commitment;
        assert_eq!(correction_values.len(), self.num_repetitions);
        self.vc_commitment_hash = vc_commitment_hash;
        self.correction_values = correction_values;
        self.random_commitment = false;
        self.state = VoleInTheHeadReceiverState::CommitmentReceived;
    }

    fn receive_random_commitment(&mut self, commitment: Self::Commitment) {
        assert_eq!(self.state, VoleInTheHeadReceiverState::New);
        let Commitment {
            vc_commitment_hash,
            correction_values,
        } = commitment;
        assert_eq!(correction_values.len(), self.num_repetitions - 1);
        self.vc_commitment_hash = vc_commitment_hash;
        self.correction_values = correction_values;
        self.random_commitment = true;
        self.state = VoleInTheHeadReceiverState::CommitmentReceived;
    }

    fn consistency_challenge_from_seed(num_repetitions: usize, seed: [u8; 32]) -> Vector<F> {
        let mut rng = ChaChaRng::from_seed(seed);
        (0..num_repetitions).map(|_| F::random(&mut rng)).collect()
    }

    fn generate_consistency_challenge_from_seed(&mut self, seed: [u8; 32]) -> Vector<F> {
        assert_eq!(self.state, VoleInTheHeadReceiverState::CommitmentReceived);
        self.consistency_challenges =
            Self::consistency_challenge_from_seed(self.num_repetitions, seed);
        self.state = VoleInTheHeadReceiverState::ConsistencyChallengeGenerated;
        self.consistency_challenges.clone()
    }

    fn generate_consistency_challenge(&mut self) -> Vector<F> {
        self.generate_consistency_challenge_from_seed(thread_rng().gen())
    }

    fn store_consistency_response(&mut self, consistency_response: Self::Response) {
        assert_eq!(
            self.state,
            VoleInTheHeadReceiverState::ConsistencyChallengeGenerated
        );
        self.u_tilde = consistency_response.vector;
        self.V_tilde_hash = consistency_response.hash;
        self.state = VoleInTheHeadReceiverState::ConsistencyChallengeResponseReceived;
    }

    fn final_challenge_from_seed(num_repetitions: usize, seed: [u8; 32]) -> Vector<F> {
        let mut rng = ChaChaRng::from_seed(seed);
        (0..num_repetitions).map(|_| F::random(&mut rng)).collect()
    }

    fn generate_final_challenge_from_seed(&mut self, seed: [u8; 32]) -> Vector<F> {
        assert_eq!(
            self.state,
            VoleInTheHeadReceiverState::ConsistencyChallengeResponseReceived
        );
        self.final_challenges = Self::final_challenge_from_seed(self.num_repetitions, seed);
        self.state = VoleInTheHeadReceiverState::FinalChallengeGenerated;
        self.final_challenges.clone()
    }

    fn generate_final_challenge(&mut self) -> Vector<F> {
        self.generate_final_challenge_from_seed(thread_rng().gen())
    }

    fn receive_decommitment(&mut self, decommitments: &Vec<VC::Decommitment>) -> bool {
        assert_eq!(
            self.state,
            VoleInTheHeadReceiverState::FinalChallengeGenerated
        );
        let tau = self.num_repetitions;
        let ell_hat = self.vole_length + self.num_repetitions;
        let log_q = F::LOG_ORDER;

        let mut W = Matrix::<F>::zeros((tau, ell_hat));

        let mut r_x_i = GF2Vector::with_capacity(ell_hat);
        r_x_i.resize(ell_hat, false);

        let mut hasher = H::new();

        // iteration i = 0
        {
            let Delta_0 = self.final_challenges[0];
            let (recomputed_commitment, mut xofs) =
                VC::recompute_commitment(log_q, &decommitments[0], Delta_0.into());
            hasher.update(recomputed_commitment.as_ref());
            let mut w_0 = W.row_mut(0);
            debug_assert_eq!(xofs.len(), F::ORDER);
            for (x, xof_x) in xofs.iter_mut().enumerate() {
                xof_x.read(r_x_i.as_raw_mut_slice());
                let Delta_0_minus_x = Delta_0 - F::from(x);
                BitMulAccumulate::bitmul_accumulate(
                    w_0.as_slice_mut().unwrap(),
                    Delta_0_minus_x,
                    r_x_i.as_raw_slice(),
                );
            }

            if !self.random_commitment {
                debug_assert_eq!(self.correction_values.len(), self.num_repetitions);
                debug_assert_eq!(self.correction_values[0].len(), self.vole_length);
                BitMulAccumulate::bitmul_accumulate(
                    &mut w_0.as_slice_mut().unwrap()[..self.vole_length],
                    Delta_0,
                    self.correction_values[0].as_raw_slice(),
                );
            }
        }

        // other iterations
        debug_assert_eq!(self.final_challenges.len(), tau);
        debug_assert_eq!(decommitments.len(), tau);
        for ((i, (&Delta_i, decommitment_i)), correction_value_i) in
            izip!(self.final_challenges.iter(), decommitments.iter())
                .enumerate()
                .skip(1)
                .zip(
                    self.correction_values
                        .iter()
                        .skip(if self.random_commitment { 0 } else { 1 }),
                )
        {
            // for i in 1..tau {
            let (recomputed_commitment, mut xofs) =
                VC::recompute_commitment(log_q, decommitment_i, Delta_i.into());
            hasher.update(recomputed_commitment.as_ref());
            let mut w_i = W.row_mut(i);
            debug_assert_eq!(xofs.len(), F::ORDER);
            for (x, xof_x) in xofs.iter_mut().enumerate() {
                xof_x.read(r_x_i.as_raw_mut_slice());
                let Delta_i_minus_x = Delta_i - F::from(x);
                BitMulAccumulate::bitmul_accumulate(
                    w_i.as_slice_mut().unwrap(),
                    Delta_i_minus_x,
                    r_x_i.as_raw_slice(),
                );
            }

            BitMulAccumulate::bitmul_accumulate(
                w_i.as_slice_mut().unwrap(),
                Delta_i,
                correction_value_i.as_raw_slice(),
            );
        }

        if hasher.finalize() != self.vc_commitment_hash {
            self.state = VoleInTheHeadReceiverState::Failed;
            return false;
        }

        // transpose W s.t. we can access it row-wise
        self.W = W.reversed_axes();
        self.W = self.W.as_standard_layout().to_owned();
        let Deltas = &self.final_challenges;
        let u_tilde = &self.u_tilde;

        // compute rhs = RW + [\tilde{u} * \Delta_1, ..., \tilde{u} * \Delta_\tau]
        let mut rhs = hash_matrix((&self.consistency_challenges).into(), (&self.W).into());
        for (mut rw_i, &u_tilde_i) in rhs.axis_iter_mut(Axis(0)).zip(u_tilde) {
            rw_i.scaled_add(u_tilde_i, Deltas);
        }
        let rw_hash = H::digest(
            rhs.as_slice()
                .expect("not in standard memory order")
                .iter()
                .flat_map(|x| x.to_repr())
                .collect::<Vec<u8>>(),
        );
        assert_eq!(self.V_tilde_hash, rw_hash);
        if self.V_tilde_hash == rw_hash {
            self.state = VoleInTheHeadReceiverState::Ready;
            true
        } else {
            self.state = VoleInTheHeadReceiverState::Failed;
            false
        }
    }

    fn get_output(&self) -> (VectorView<F>, MatrixView<F>) {
        assert_eq!(self.state, VoleInTheHeadReceiverState::Ready);
        let tau = self.num_repetitions;
        let ell = self.vole_length;
        assert_eq!(self.W.shape(), &[ell + tau, tau]);
        (
            (&self.final_challenges).into(),
            self.W.slice(s![..self.vole_length, ..]),
        )
    }
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use crate::common::Block128;
    use crate::gf2psmall::{GF2p10, GF2p8};
    use crate::primitives::{Aes128CtrLdPrg, Blake3LE};
    use crate::veccom::GgmVecCom;
    use blake3;
    use rand::{thread_rng, Rng};

    fn test_voleith_correctness_with_param<F: SmallGF, const RANDOM: bool>() {
        let vole_length = 16;
        let num_repetitions = 5;
        type VC = GgmVecCom<Block128, Aes128CtrLdPrg, blake3::Hasher, Blake3LE<Block128>>;

        let mut sender = VoleInTheHeadSenderFromVC::<F, VC>::new(vole_length, num_repetitions);
        let mut receiver = VoleInTheHeadReceiverFromVC::<F, VC>::new(vole_length, num_repetitions);

        let messages = if RANDOM {
            Default::default()
        } else {
            let mut messages = GF2Vector::with_capacity(vole_length);
            messages.resize(vole_length, false);
            thread_rng().fill(messages.as_raw_mut_slice());
            messages
        };

        if RANDOM {
            let commitment = sender.commit_random();
            receiver.receive_random_commitment(commitment);
        } else {
            let commitment = sender.commit_message(messages.clone());
            receiver.receive_commitment(commitment);
        }
        let consistency_challenge = receiver.generate_consistency_challenge();
        let consistency_response = sender.consistency_check_respond(consistency_challenge);
        receiver.store_consistency_response(consistency_response);

        let final_challenge = receiver.generate_final_challenge();
        let decommitment = sender.decommit(final_challenge);
        let receiver_accepts = receiver.receive_decommitment(&decommitment);
        assert!(receiver_accepts);

        let (u, V) = sender.get_output();
        if !RANDOM {
            assert_eq!(u, messages.as_ref());
        }
        let (Deltas, W) = receiver.get_output();
        assert_eq!(u.len(), vole_length);
        assert_eq!(Deltas.len(), num_repetitions);
        assert_eq!(V.shape(), &[vole_length, num_repetitions]);
        assert_eq!(W.shape(), &[vole_length, num_repetitions]);
        for i in 0..vole_length {
            if u[i] {
                assert_eq!(&V.row(i) + &Deltas, W.row(i));
            } else {
                assert_eq!(V.row(i), W.row(i));
            }
        }
    }

    #[test]
    fn test_voleith_correctness_with_msg_gf2p8() {
        test_voleith_correctness_with_param::<GF2p8, false>();
    }

    #[test]
    fn test_voleith_correctness_random_gf2p8() {
        test_voleith_correctness_with_param::<GF2p8, true>();
    }

    #[test]
    fn test_voleith_correctness_with_msg_gf2p10() {
        test_voleith_correctness_with_param::<GF2p10, false>();
    }

    #[test]
    fn test_voleith_correctness_random_gf2p10() {
        test_voleith_correctness_with_param::<GF2p10, true>();
    }
}
