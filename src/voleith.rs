use crate::arithmetic::{hash_bitvector_and_matrix, hash_matrix};
use crate::field::GF2p8;
use crate::veccom::VecCom;
use bitvec::{slice::BitSlice, vec::BitVec};
use core::marker::PhantomData;
use core::mem;
use digest::XofReader;
use ff::Field;
use ndarray::{s, Array1, Array2, ArrayView1, ArrayView2, Axis};
use rand::thread_rng;
// use blake3::{Hash as Blake3Hash, Hasher as Blake3Hasher};

type GF2Vector = BitVec<u8>;
type GF2View = BitSlice<u8>;
type GF2p8Vector = Array1<GF2p8>;
type GF2p8VectorView<'a> = ArrayView1<'a, GF2p8>;
type GF2p8Matrix = Array2<GF2p8>;
type GF2p8MatrixView<'a> = ArrayView2<'a, GF2p8>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VoleInTheHeadSenderState {
    New,
    Committed,
    RespondedToConsistencyChallenge,
    Ready,
}

#[allow(non_snake_case)]
pub struct VoleInTheHeadSender<VC: VecCom> {
    vole_length: usize,
    num_repetitions: usize,
    state: VoleInTheHeadSenderState,
    u: GF2Vector,
    V: GF2p8Matrix,
    decommitment_keys: Vec<VC::DecommitmentKey>,
    _phantom_hc: PhantomData<VC>,
}

#[allow(non_snake_case)]
impl<'a, VC: VecCom> VoleInTheHeadSender<VC> {
    pub fn new(vole_length: usize, num_repetitions: usize) -> Self {
        let ell_hat = vole_length + num_repetitions;
        let mut output = Self {
            vole_length,
            num_repetitions,
            state: VoleInTheHeadSenderState::New,
            u: GF2Vector::new(),
            V: GF2p8Matrix::zeros((num_repetitions, ell_hat)),
            decommitment_keys: Vec::with_capacity(num_repetitions),
            _phantom_hc: PhantomData,
        };
        output.u.resize(ell_hat, false);
        output
    }

    pub fn commit(&mut self) -> (Vec<VC::Commitment>, Vec<GF2Vector>) {
        assert_eq!(self.state, VoleInTheHeadSenderState::New);

        let log_q = GF2p8::LOG_ORDER;
        let tau = self.num_repetitions;
        let ell_hat = self.vole_length + self.num_repetitions;

        // let mut hasher = Blake3Hasher::new();
        let mut commitments = Vec::with_capacity(tau);
        let mut correction_values = Vec::with_capacity(tau - 1);

        // iteration i = 0
        {
            let (commitment, decommitment_key, mut xofs) = VC::commit(log_q);
            self.decommitment_keys.push(decommitment_key);
            // hasher.update(commitment.as_ref());
            commitments.push(commitment);
            // iteration x = 0
            let u_0 = &mut self.u;
            xofs[0].read(u_0.as_raw_mut_slice());
            let mut v_0 = self.V.row_mut(0);
            for x in 1..=255 {
                let mut r_x_0 = GF2Vector::with_capacity(ell_hat);
                r_x_0.resize(ell_hat, false);
                xofs[x].read(r_x_0.as_raw_mut_slice());
                *u_0 ^= &r_x_0;
                for (i, b) in r_x_0.iter().enumerate() {
                    if *b {
                        v_0[i] -= GF2p8(x as u8);
                    }
                }
            }
        }

        // other iterations
        for i in 1..tau {
            let (commitment, decommitment_key, mut xofs) = VC::commit(log_q);
            self.decommitment_keys.push(decommitment_key);
            // hasher.update(commitment.as_ref());
            commitments.push(commitment);
            // iteration x = 0
            let mut u_i = {
                let mut r_0_1 = GF2Vector::with_capacity(ell_hat);
                r_0_1.resize(ell_hat, false);
                xofs[0].read(r_0_1.as_raw_mut_slice());
                r_0_1
            };
            let mut v_i = self.V.row_mut(i);
            for x in 1..=255 {
                let mut r_x_i = GF2Vector::with_capacity(ell_hat);
                r_x_i.resize(ell_hat, false);
                xofs[x].read(r_x_i.as_raw_mut_slice());
                u_i ^= &r_x_i;
                for (i, b) in r_x_i.iter().enumerate() {
                    if *b {
                        v_i[i] -= GF2p8(x as u8);
                    }
                }
            }
            u_i ^= &self.u;
            correction_values.push(u_i);
        }

        // let commitment = hasher.finalize();

        // transpose V such that we can now access it row-wise
        {
            // reversed_axes consumes self, but we cannot move out of the &mut ref.
            // Hence, we do the swapping with a dummy value.
            let mut tmp = GF2p8Matrix::default((1, 1));
            mem::swap(&mut tmp, &mut self.V);
            tmp = tmp.reversed_axes();
            mem::swap(&mut tmp, &mut self.V);
        }

        self.state = VoleInTheHeadSenderState::Committed;
        // (commitment, correction_values)
        (commitments, correction_values)
    }

    pub fn consistency_check_respond(
        &mut self,
        random_points: GF2p8Vector,
    ) -> (GF2p8Vector, GF2p8Matrix) {
        assert_eq!(self.state, VoleInTheHeadSenderState::Committed);
        assert_eq!(random_points.len(), self.num_repetitions);

        let (u_tilde, V_tilde) =
            hash_bitvector_and_matrix((&random_points).into(), &self.u, (&self.V).into());

        // ignore the extra entries
        self.u.resize(self.vole_length, false);
        // cannot resize the matrix V, so keep in mind to ignore the last rows

        self.state = VoleInTheHeadSenderState::RespondedToConsistencyChallenge;
        (u_tilde, V_tilde)
    }

    pub fn decommit(&mut self, Deltas: GF2p8Vector) -> Vec<VC::Decommitment> {
        assert_eq!(
            self.state,
            VoleInTheHeadSenderState::RespondedToConsistencyChallenge
        );
        assert_eq!(Deltas.len(), self.num_repetitions);
        let log_q = GF2p8::LOG_ORDER;
        let decommitments = Deltas
            .iter()
            .zip(self.decommitment_keys.iter())
            .map(|(Delta_i, &decommitment_key)| {
                VC::decommit(log_q, decommitment_key, Delta_i.0 as usize)
            })
            .collect();
        self.state = VoleInTheHeadSenderState::Ready;
        decommitments
    }

    pub fn get_output(&'a self) -> (&'a GF2View, GF2p8MatrixView<'a>) {
        assert_eq!(self.state, VoleInTheHeadSenderState::Ready);
        let tau = self.num_repetitions;
        let ell = self.vole_length;
        assert_eq!(self.V.shape(), &[ell + tau, tau]);
        (&self.u[..ell], self.V.slice(s![..ell, ..]))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VoleInTheHeadReceiverState {
    New,
    ConsistencyChallengeGenerated,
    ConsistencyChallengeResponseReceived,
    FinalChallengeGenerated,
    Failed,
    Ready,
}

#[allow(non_snake_case)]
pub struct VoleInTheHeadReceiver<VC: VecCom> {
    vole_length: usize,
    num_repetitions: usize,
    state: VoleInTheHeadReceiverState,
    // commitment: Option<Blake3Hash>,
    commitments: Vec<VC::Commitment>,
    correction_values: Vec<GF2Vector>,
    consistency_challenges: GF2p8Vector,
    u_tilde: GF2p8Vector,
    V_tilde: GF2p8Matrix,
    final_challenges: GF2p8Vector,
    W: GF2p8Matrix,
    // decommitment_keys: Vec<VC::DecommitmentKey>,
    _phantom_hc: PhantomData<VC>,
}

#[allow(non_snake_case)]
impl<'a, VC: VecCom> VoleInTheHeadReceiver<VC> {
    pub fn new(vole_length: usize, num_repetitions: usize) -> Self {
        let output = Self {
            vole_length,
            num_repetitions,
            state: VoleInTheHeadReceiverState::New,
            // commitment: None,
            commitments: Vec::new(),
            correction_values: Vec::new(),
            consistency_challenges: Default::default(),
            u_tilde: Default::default(),
            V_tilde: Default::default(),
            final_challenges: Default::default(),
            W: Default::default(),
            _phantom_hc: PhantomData,
        };
        output
    }

    pub fn generate_consistency_challenge(
        &mut self,
        // commitment: Blake3Hash,
        commitments: Vec<VC::Commitment>,
        correction_values: Vec<GF2Vector>,
    ) -> GF2p8Vector {
        assert_eq!(self.state, VoleInTheHeadReceiverState::New);
        // self.commitment = Some(commitment);
        self.commitments = commitments;
        self.correction_values = correction_values;
        self.consistency_challenges = (0..self.num_repetitions)
            .map(|_| GF2p8::random(thread_rng()))
            .collect();
        self.state = VoleInTheHeadReceiverState::ConsistencyChallengeGenerated;
        self.consistency_challenges.clone()
    }

    pub fn store_consistency_response(&mut self, consistency_response: (GF2p8Vector, GF2p8Matrix)) {
        assert_eq!(
            self.state,
            VoleInTheHeadReceiverState::ConsistencyChallengeGenerated
        );
        (self.u_tilde, self.V_tilde) = consistency_response;
        self.state = VoleInTheHeadReceiverState::ConsistencyChallengeResponseReceived;
    }

    pub fn generate_final_challenge(&mut self) -> GF2p8Vector {
        assert_eq!(
            self.state,
            VoleInTheHeadReceiverState::ConsistencyChallengeResponseReceived
        );
        self.final_challenges = (0..self.num_repetitions)
            .map(|_| GF2p8::random(thread_rng()))
            .collect();
        self.state = VoleInTheHeadReceiverState::FinalChallengeGenerated;
        self.final_challenges.clone()
    }

    pub fn receive_decommitment(&mut self, decommitments: &[VC::Decommitment]) -> bool {
        assert_eq!(
            self.state,
            VoleInTheHeadReceiverState::FinalChallengeGenerated
        );
        let tau = self.num_repetitions;
        let ell_hat = self.vole_length + self.num_repetitions;
        let log_q = GF2p8::LOG_ORDER;

        let mut W = GF2p8Matrix::zeros((tau, ell_hat));

        let mut r_x_i = GF2Vector::with_capacity(ell_hat);
        r_x_i.resize(ell_hat, false);

        // iteration i = 0
        {
            let Delta_i = self.final_challenges[0];
            let xofs = VC::verify(
                log_q,
                &self.commitments[0],
                &decommitments[0],
                Delta_i.0 as usize,
            );
            if xofs.is_none() {
                self.state = VoleInTheHeadReceiverState::Failed;
                return false;
            }
            let mut xofs = xofs.unwrap();
            let mut w_0 = W.row_mut(0);
            for x in 0..=255 {
                xofs[x].read(r_x_i.as_raw_mut_slice());
                let Delta_i_minus_x = Delta_i - GF2p8(x as u8);
                for (j, b) in r_x_i.iter().enumerate() {
                    if *b {
                        w_0[j] += Delta_i_minus_x;
                    }
                }
            }
        }

        // other iterations
        for i in 1..tau {
            let Delta_i = self.final_challenges[i];
            let xofs = VC::verify(
                log_q,
                &self.commitments[i],
                &decommitments[i],
                Delta_i.0 as usize,
            );
            if xofs.is_none() {
                self.state = VoleInTheHeadReceiverState::Failed;
                return false;
            }
            let mut xofs = xofs.unwrap();
            let mut w_i = W.row_mut(i);
            for x in 0..=255 {
                xofs[x].read(r_x_i.as_raw_mut_slice());
                let Delta_i_minus_x = Delta_i - GF2p8(x as u8);
                for (j, b) in r_x_i.iter().enumerate() {
                    if *b {
                        w_i[j] += Delta_i_minus_x;
                    }
                }
            }

            for (j, b) in self.correction_values[i - 1].iter().enumerate() {
                if *b {
                    w_i[j] += Delta_i;
                }
            }
        }

        // transpose W s.t. we can access it row-wise
        self.W = W.reversed_axes();
        let mut lhs = Default::default();
        let Deltas = &self.final_challenges;
        let u_tilde = &self.u_tilde;
        mem::swap(&mut lhs, &mut self.V_tilde);

        for (mut v_tilde_i, &u_tilde_i) in lhs.axis_iter_mut(Axis(0)).zip(u_tilde) {
            v_tilde_i.scaled_add(u_tilde_i, Deltas);
        }
        // lhs is now \tilde{V] + [\tilde{u} * \Delta_1, ..., \tilde{u} * \Delta_\tau]

        // compute rhs = RW
        let rhs = hash_matrix((&self.consistency_challenges).into(), (&self.W).into());
        assert_eq!(lhs, rhs);
        if lhs == rhs {
            self.state = VoleInTheHeadReceiverState::Ready;
            true
        } else {
            self.state = VoleInTheHeadReceiverState::Failed;
            false
        }
    }

    pub fn get_output(&'a self) -> (GF2p8VectorView, GF2p8MatrixView<'a>) {
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
    use crate::veccom::GgmVecCom;
    use crate::primitives::{Aes128CtrLdPrg, Blake3LE};
    use blake3;

    #[test]
    fn test_voleith_correctness() {
        let vole_length = 16;
        let num_repetitions = 5;
        type VC = GgmVecCom<Block128, Aes128CtrLdPrg, blake3::Hasher, Blake3LE<Block128>>;

        let mut sender = VoleInTheHeadSender::<VC>::new(vole_length, num_repetitions);
        let mut receiver = VoleInTheHeadReceiver::<VC>::new(vole_length, num_repetitions);

        let (commitment, correction_values) = sender.commit();
        let consistency_challenge =
            receiver.generate_consistency_challenge(commitment, correction_values);
        let consistency_response = sender.consistency_check_respond(consistency_challenge);
        receiver.store_consistency_response(consistency_response);

        let final_challenge = receiver.generate_final_challenge();
        let decommitment = sender.decommit(final_challenge);
        let receiver_accepts = receiver.receive_decommitment(&decommitment);
        assert!(receiver_accepts);

        let (u, V) = sender.get_output();
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
}
