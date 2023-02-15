use crate::field::{GF2Vector, GF2p128, VecToGF2p128};
use crate::voleith::{VoleInTheHeadReceiver, VoleInTheHeadSender};
use core::fmt;
use ndarray::{Array1, Axis};

type Vector<T> = Array1<T>;

pub trait HomComSender {
    type Choice: Clone;
    type Tag;
    type Commitment: Clone;
    type Challenge: Clone;
    type Response: Clone;
    type Decommitment: Clone;

    fn new(num_bit_commitments: usize) -> Self;
    fn commit_to_bits(&mut self, bits: GF2Vector) -> (Vec<Self::Tag>, Self::Commitment);
    fn commit_prove_consistency(&mut self, cons_challenge: Self::Challenge) -> Self::Response;
    fn transfer(&mut self, choice: Self::Choice) -> Self::Decommitment;
}

pub trait HomComReceiver {
    type Choice: Clone;
    type GlobalKey;
    type Key;
    type Commitment: Clone;
    type Challenge: Clone;
    type Response: Clone;
    type Decommitment: Clone;

    fn new(num_bit_commitments: usize) -> Self;
    fn generate_challenge_from_seed(seed: [u8; 32]) -> Self::Challenge;
    fn commit_send_challenge_from_seed(
        &mut self,
        cons_commitment: Self::Commitment,
        seed: [u8; 32],
    ) -> Self::Challenge;
    fn commit_send_challenge(&mut self, cons_commitment: Self::Commitment) -> Self::Challenge;
    fn commit_receive_response(&mut self, cons_response: Self::Response);
    fn generate_choice_from_seed(seed: [u8; 32]) -> Self::Choice;
    fn choose_from_seed(&mut self, seed: [u8; 32]) -> Self::Choice;
    fn choose(&mut self) -> Self::Choice;
    fn receive(
        &mut self,
        decommitment: Self::Decommitment,
    ) -> Option<(Self::GlobalKey, Vec<Self::Key>)>;
}

pub struct HomCom128SenderFromVitH<VitHS> {
    vith_sender: VitHS,
}

impl<VitHS> HomComSender for HomCom128SenderFromVitH<VitHS>
where
    VitHS: VoleInTheHeadSender,
    VitHS::Field: VecToGF2p128<GF2p128> + fmt::Debug,
{
    type Choice = Vector<VitHS::Field>;
    type Tag = GF2p128;
    type Commitment = VitHS::Commitment;
    type Challenge = VitHS::Challenge;
    type Response = VitHS::Response;
    type Decommitment = VitHS::Decommitment;

    fn new(num_bit_commitments: usize) -> Self {
        Self {
            vith_sender: VitHS::new(num_bit_commitments, VitHS::Field::VECTOR_SIZE),
        }
    }

    #[allow(non_snake_case)]
    fn commit_to_bits(&mut self, bits: GF2Vector) -> (Vec<Self::Tag>, Self::Commitment) {
        let commitment = self.vith_sender.commit_message(bits);
        let (_, V) = self.vith_sender.get_output();
        debug_assert!(V.is_standard_layout());
        let tags: Vec<GF2p128> = V
            .axis_iter(Axis(0))
            .map(|r| VecToGF2p128::convert(r.as_slice().expect("row as_slice failed")))
            .collect();
        (tags, commitment)
    }

    fn commit_prove_consistency(&mut self, cons_challenge: Self::Challenge) -> Self::Response {
        self.vith_sender.consistency_check_respond(cons_challenge)
    }

    fn transfer(&mut self, choice: Self::Choice) -> Self::Decommitment {
        self.vith_sender.decommit(choice)
    }

    //     fn lift_gf2p8(tags: &[Self::Tag]) -> Self::Tag {
    //         assert_eq!(tags.len(), 8);
    //         InnerProduct::inner_product(tags.iter(), GF2p128::GF2P8_EMBEDDING_POX.iter())
    //     }
}

pub struct HomCom128ReceiverFromVitH<VitHR> {
    vith_receiver: VitHR,
}

impl<VitHR> HomComReceiver for HomCom128ReceiverFromVitH<VitHR>
where
    VitHR: VoleInTheHeadReceiver,
    VitHR::Field: VecToGF2p128<GF2p128>,
{
    type Choice = Vector<VitHR::Field>;
    type GlobalKey = GF2p128;
    type Key = GF2p128;
    type Commitment = VitHR::Commitment;
    type Challenge = VitHR::Challenge;
    type Response = VitHR::Response;
    type Decommitment = VitHR::Decommitment;

    fn new(num_bit_commitments: usize) -> Self {
        Self {
            vith_receiver: VitHR::new(num_bit_commitments, VitHR::Field::VECTOR_SIZE),
        }
    }

    fn generate_challenge_from_seed(seed: [u8; 32]) -> Self::Challenge {
        VitHR::consistency_challenge_from_seed(VitHR::Field::VECTOR_SIZE, seed)
    }

    fn commit_send_challenge_from_seed(
        &mut self,
        cons_commitment: Self::Commitment,
        seed: [u8; 32],
    ) -> Self::Challenge {
        self.vith_receiver.receive_commitment(cons_commitment);
        self.vith_receiver
            .generate_consistency_challenge_from_seed(seed)
    }

    fn commit_send_challenge(&mut self, cons_commitment: Self::Commitment) -> Self::Challenge {
        self.vith_receiver.receive_commitment(cons_commitment);
        self.vith_receiver.generate_consistency_challenge()
    }

    fn commit_receive_response(&mut self, cons_response: Self::Response) {
        self.vith_receiver.store_consistency_response(cons_response)
    }

    fn generate_choice_from_seed(seed: [u8; 32]) -> Self::Choice {
        VitHR::final_challenge_from_seed(VitHR::Field::VECTOR_SIZE, seed)
    }

    fn choose_from_seed(&mut self, seed: [u8; 32]) -> Self::Choice {
        self.vith_receiver.generate_final_challenge_from_seed(seed)
    }

    fn choose(&mut self) -> Self::Choice {
        self.vith_receiver.generate_final_challenge()
    }

    #[allow(non_snake_case)]
    fn receive(
        &mut self,
        decommitment: Self::Decommitment,
    ) -> Option<(Self::GlobalKey, Vec<Self::Key>)> {
        if !self.vith_receiver.receive_decommitment(&decommitment) {
            return None;
        }
        let (Deltas, W) = self.vith_receiver.get_output();
        debug_assert!(W.is_standard_layout());
        let global_key = VecToGF2p128::convert(Deltas.as_slice().unwrap());
        let keys: Vec<GF2p128> = W
            .axis_iter(Axis(0))
            .map(|r| VecToGF2p128::convert(r.as_slice().expect("row as_slice failed")))
            .collect();
        Some((global_key, keys))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Block128;
    use crate::primitives::{Aes128CtrLdPrg, Blake3LE};
    use crate::veccom::GgmVecCom;
    use crate::voleith::{VoleInTheHeadReceiverFromVC, VoleInTheHeadSenderFromVC};
    use blake3;
    use itertools::izip;
    use rand::{thread_rng, Rng};

    type VC = GgmVecCom<Block128, Aes128CtrLdPrg, blake3::Hasher, Blake3LE<Block128>>;
    type VitHSender = VoleInTheHeadSenderFromVC<VC>;
    type VitHReceiver = VoleInTheHeadReceiverFromVC<VC>;

    #[test]
    fn test_correctness() {
        let num_bits = 42;
        let mut hc_sender = HomCom128SenderFromVitH::<VitHSender>::new(num_bits);
        let mut hc_receiver = HomCom128ReceiverFromVitH::<VitHReceiver>::new(num_bits);

        let bits = {
            let mut bits = GF2Vector::with_capacity(num_bits);
            bits.resize(num_bits, false);
            thread_rng().fill(bits.as_raw_mut_slice());
            bits
        };

        let (tags, commitment) = hc_sender.commit_to_bits(bits.clone());
        let challenge = hc_receiver.commit_send_challenge(commitment);
        let response = hc_sender.commit_prove_consistency(challenge);
        hc_receiver.commit_receive_response(response);

        let choice = hc_receiver.choose();
        let decommitment = hc_sender.transfer(choice);

        let receiver_output = hc_receiver.receive(decommitment);
        assert!(receiver_output.is_some());
        let (global_key, keys) = receiver_output.unwrap();

        assert_eq!(tags.len(), num_bits);
        assert_eq!(keys.len(), num_bits);
        for (b_i, key_i, tag_i) in izip!(
            bits.bits.iter().by_vals(),
            keys.iter().copied(),
            tags.iter().copied()
        ) {
            if b_i {
                assert_eq!(global_key + key_i, tag_i);
            } else {
                assert_eq!(key_i, tag_i);
            }
        }
    }
}
