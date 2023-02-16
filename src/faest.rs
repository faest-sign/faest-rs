use crate::aes::{Aes, AesState, KeyExpansion, EXTENDED_WITNESS_BYTE_SIZE, ROUND_CONSTANTS};
use crate::field::{BytesRepr, InnerProduct};
use crate::gf2::GF2Vector;
use crate::gf2p128::GF2p128;
use crate::gf2psmall::GF2p8;
use crate::homcom::{HomComReceiver, HomComSender};
use ff::Field;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaChaRng;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SecretKey {
    Aes128Key { key: [u8; 16] },
    Aes192Key { key: [u8; 24] },
    Aes256Key { key: [u8; 32] },
}

impl SecretKey {
    pub fn as_slice(&self) -> &[u8] {
        match self {
            Self::Aes128Key { key } => key,
            Self::Aes192Key { key } => key,
            Self::Aes256Key { key } => key,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, bincode::Encode)]
pub struct PublicKey {
    pub input: [u8; 16],
    pub output: [u8; 16],
}

pub fn keygen(variant: Aes) -> (SecretKey, PublicKey) {
    let mut rng = thread_rng();

    loop {
        let key = {
            let mut key = vec![0u8; variant.get_key_size()];
            rng.fill(key.as_mut_slice());
            key
        };
        let key_as_gf2p8: Vec<_> = key.iter().copied().map(GF2p8).collect();
        let (round_keys, inv_outs) =
            KeyExpansion::expand_and_collect_inv_outputs(variant, &key_as_gf2p8);
        if inv_outs.iter().any(|x| x.iter().any(|y| y.0 == 0)) {
            continue;
        }
        let input: [u8; 16] = rng.gen();
        let (output, inv_outs) =
            AesState::from(input).encrypt_and_collect_inv_outputs(variant, &round_keys);
        if inv_outs.iter().any(|x| x.0.iter().any(|y| y.0 == 0)) {
            continue;
        }
        return (
            match variant {
                Aes::Aes128 => SecretKey::Aes128Key {
                    key: key.try_into().unwrap(),
                },
                Aes::Aes192 => SecretKey::Aes192Key {
                    key: key.try_into().unwrap(),
                },
                Aes::Aes256 => SecretKey::Aes256Key {
                    key: key.try_into().unwrap(),
                },
            },
            PublicKey {
                input,
                output: output.0.map(Into::into),
            },
        );
    }
}

pub trait Prover {
    type Commitment: Clone;
    type Challenge1: Clone;
    type Response: Clone;
    type Challenge2: Clone;
    type Proof: Clone;
    type Choice: Clone;
    type Decommitment: Clone;

    fn new(secret_key: SecretKey, public_key: PublicKey) -> Self;
    fn commit(&mut self) -> Self::Commitment;
    fn commit_prove_consistency(&mut self, challenge: Self::Challenge1) -> Self::Response;
    fn prove(&mut self, challenge: Self::Challenge2) -> Self::Proof;
    fn transfer(&mut self, choice: Self::Choice) -> Self::Decommitment;
}

pub trait Verifier {
    type Commitment: Clone;
    type Challenge1: Clone;
    type Response: Clone;
    type Challenge2: Clone;
    type Proof: Clone;
    type Choice: Clone;
    type Decommitment: Clone;
    fn new(public_key: PublicKey) -> Self;
    fn generate_commit_challenge_from_seed(seed: [u8; 32]) -> Self::Challenge1;
    fn commit_send_challenge_from_seed(
        &mut self,
        commitment: Self::Commitment,
        seed: [u8; 32],
    ) -> Self::Challenge1;
    fn commit_send_challenge(&mut self, commitment: Self::Commitment) -> Self::Challenge1;
    fn commit_receive_response(&mut self, response: Self::Response);
    fn generate_challenge_from_seed(seed: [u8; 32]) -> Self::Challenge2;
    fn send_challenge_from_seed(&mut self, seed: [u8; 32]) -> Self::Challenge2;
    fn send_challenge(&mut self) -> Self::Challenge2;
    fn generate_choice_from_seed(seed: [u8; 32]) -> Self::Choice;
    fn choose_from_seed(&mut self, proof: Self::Proof, seed: [u8; 32]) -> Self::Choice;
    fn choose(&mut self, proof: Self::Proof) -> Self::Choice;
    fn verify(&mut self, decommitment: Self::Decommitment) -> bool;
}

fn rotate_gf2p8(x: GF2p8) -> GF2p8 {
    x + x.rotate_left(1) + x.rotate_left(2) + x.rotate_left(3) + x.rotate_left(4)
}

fn lift_into_gf2p8_subfield(tags: &[GF2p128]) -> GF2p128 {
    assert_eq!(tags.len(), 8);
    InnerProduct::inner_product(tags.iter(), GF2p128::GF2P8_EMBEDDING_POX.iter())
}

fn lift_into_rotated_gf2p8_subfield(tags: &[GF2p128]) -> GF2p128 {
    assert_eq!(tags.len(), 8);
    let mut rotated_tags: [GF2p128; 8] = tags.try_into().unwrap();
    for i in 0..8 {
        rotated_tags[i] +=
            tags[(7 + i) % 8] + tags[(6 + i) % 8] + tags[(5 + i) % 8] + tags[(4 + i) % 8];
    }
    InnerProduct::inner_product(rotated_tags.iter(), GF2p128::GF2P8_EMBEDDING_POX.iter())
}

pub struct FaestProverFromHC<HCS>
where
    HCS: HomComSender,
{
    secret_key: SecretKey,
    public_key: PublicKey,
    hc_sender: HCS,
    witness: GF2Vector,
    tags: Vec<GF2p128>,
    mask: GF2p128,
    mask_tag: GF2p128,
    key_schedule_intermediate_states: Vec<[GF2p8; 4]>,
    key_schedule_intermediate_state_tags: Vec<[GF2p128; 4]>,
    round_keys: Vec<[GF2p8; 16]>,
    round_key_tags: Vec<[GF2p128; 16]>,
    intermediate_states: Vec<[GF2p8; 16]>,
    intermediate_state_tags: Vec<[GF2p128; 16]>,
    rotated_intermediate_states: Vec<[GF2p8; 16]>,
    rotated_intermediate_state_tags: Vec<[GF2p128; 16]>,
}

impl<HCS> Prover for FaestProverFromHC<HCS>
where
    HCS: HomComSender<Tag = GF2p128>,
{
    type Commitment = HCS::Commitment;
    type Challenge1 = HCS::Challenge;
    type Response = HCS::Response;
    type Challenge2 = [u8; 32];
    type Proof = (GF2p128, GF2p128);
    type Choice = HCS::Choice;
    type Decommitment = HCS::Decommitment;

    fn new(secret_key: SecretKey, public_key: PublicKey) -> Self {
        Self {
            secret_key,
            public_key,
            hc_sender: HCS::new(EXTENDED_WITNESS_BYTE_SIZE * 8 + 128),
            witness: Default::default(),
            tags: Default::default(),
            mask: Default::default(),
            mask_tag: Default::default(),
            key_schedule_intermediate_states: Default::default(),
            key_schedule_intermediate_state_tags: Default::default(),
            round_keys: Default::default(),
            round_key_tags: Default::default(),
            intermediate_states: Default::default(),
            intermediate_state_tags: Default::default(),
            rotated_intermediate_states: Default::default(),
            rotated_intermediate_state_tags: Default::default(),
        }
    }

    fn commit(&mut self) -> HCS::Commitment {
        let witness = Aes::compute_extended_witness(
            Aes::Aes128,
            self.secret_key.as_slice(),
            &self.public_key.input,
        );

        assert!(witness.is_some());
        let (mut witness, output) = witness.unwrap();
        assert_eq!(output, self.public_key.output);
        witness.extend((0..16).map(|_| thread_rng().gen::<u8>())); // add 128 more random bit
                                                                   // for the mask
        self.witness = GF2Vector::from_vec(witness);
        assert_eq!(self.witness.len(), 8 * EXTENDED_WITNESS_BYTE_SIZE + 128);
        let (tags, commitment) = self.hc_sender.commit_to_bits(self.witness.clone());
        self.tags = tags;
        commitment
    }

    fn commit_prove_consistency(&mut self, challenge: HCS::Challenge) -> HCS::Response {
        self.hc_sender.commit_prove_consistency(challenge)
    }

    fn prove(&mut self, challenge_seed: [u8; 32]) -> (GF2p128, GF2p128) {
        self.lift_commitments();
        self.aggregates_conditions(challenge_seed)
    }

    fn transfer(&mut self, choice: HCS::Choice) -> HCS::Decommitment {
        self.hc_sender.transfer(choice)
    }
}

impl<HCS> FaestProverFromHC<HCS>
where
    HCS: HomComSender<Tag = GF2p128>,
{
    fn lift_commitments(&mut self) {
        let witness_bytes = self.witness.as_raw_slice();
        assert_eq!(witness_bytes.len(), EXTENDED_WITNESS_BYTE_SIZE + 16);

        // prepare uniform mask
        (self.mask, self.mask_tag) = {
            let byte_offset = EXTENDED_WITNESS_BYTE_SIZE;
            let bit_offset = 8 * byte_offset;
            let mask = GF2p128::from_repr(witness_bytes[byte_offset..].try_into().unwrap());
            let mask_tag = InnerProduct::<GF2p128>::inner_product(
                self.tags[bit_offset..].iter().copied(),
                (0..128).into_iter().map(|i| GF2p128::from_u128(1 << i)),
            );
            (mask, mask_tag)
        };

        self.key_schedule_intermediate_states = Vec::with_capacity(10);
        self.key_schedule_intermediate_state_tags = Vec::with_capacity(10);
        self.round_keys = Vec::with_capacity(11);
        self.round_key_tags = Vec::with_capacity(11);

        // handle the first round key
        {
            let witness_key_chunk: [u8; 16] = witness_bytes[0..16].try_into().unwrap();
            let tag_key_chunk = &self.tags[0..128];
            let round_key = witness_key_chunk.map(GF2p8);
            self.round_keys.push(round_key);
            let mut round_key_tags: [GF2p128; 16] = Default::default();
            for j in 0..16 {
                round_key_tags[j] = lift_into_gf2p8_subfield(&tag_key_chunk[8 * j..8 * (j + 1)]);
            }
            self.round_key_tags.push(round_key_tags);
        }
        // handle the other round keys
        {
            // iterate through the witness and the tags in 32 bit chunks
            let mut witness_word_chunks_it = witness_bytes[16..16 + 4 * 10].chunks_exact(4);
            let mut tag_word_chunks_it = self.tags[128..128 + 32 * 10].chunks_exact(32);
            assert_eq!(witness_word_chunks_it.len(), 10);
            assert_eq!(tag_word_chunks_it.len(), 10);

            for round_key_i in 1..11 {
                {
                    let wit_chunk = <&[u8] as TryInto<[u8; 4]>>::try_into(
                        witness_word_chunks_it.next().unwrap(),
                    )
                    .unwrap()
                    .map(GF2p8);
                    self.key_schedule_intermediate_states.push(wit_chunk);
                    let mut round_key = self.round_keys[round_key_i - 1];
                    let after_sbox = wit_chunk.map(rotate_gf2p8).map(|x| x + GF2p8(0x63));
                    round_key[0] += after_sbox[0] + ROUND_CONSTANTS[round_key_i - 1];
                    round_key[1] += after_sbox[1];
                    round_key[2] += after_sbox[2];
                    round_key[3] += after_sbox[3];
                    for j in 1..4 {
                        round_key[j * 4] += round_key[(j - 1) * 4];
                        round_key[j * 4 + 1] += round_key[(j - 1) * 4 + 1];
                        round_key[j * 4 + 2] += round_key[(j - 1) * 4 + 2];
                        round_key[j * 4 + 3] += round_key[(j - 1) * 4 + 3];
                    }
                    self.round_keys.push(round_key);
                }
                {
                    let tag_chunk = tag_word_chunks_it.next().unwrap();
                    {
                        let mut inv_out = [GF2p128::ZERO; 4];
                        for j in 0..4 {
                            inv_out[j] = lift_into_gf2p8_subfield(&tag_chunk[8 * j..8 * (j + 1)]);
                        }
                        self.key_schedule_intermediate_state_tags.push(inv_out);
                    }
                    let mut after_sbox = [GF2p128::ZERO; 4];
                    for j in 0..4 {
                        after_sbox[j] =
                            lift_into_rotated_gf2p8_subfield(&tag_chunk[8 * j..8 * (j + 1)]);
                    }
                    let mut round_key_tag = self.round_key_tags[round_key_i - 1];
                    round_key_tag[0] += after_sbox[0];
                    round_key_tag[1] += after_sbox[1];
                    round_key_tag[2] += after_sbox[2];
                    round_key_tag[3] += after_sbox[3];
                    for i in 1..4 {
                        round_key_tag[i * 4] += round_key_tag[(i - 1) * 4];
                        round_key_tag[i * 4 + 1] += round_key_tag[(i - 1) * 4 + 1];
                        round_key_tag[i * 4 + 2] += round_key_tag[(i - 1) * 4 + 2];
                        round_key_tag[i * 4 + 3] += round_key_tag[(i - 1) * 4 + 3];
                    }
                    self.round_key_tags.push(round_key_tag);
                }
            }
        }

        // prepare intermediate states
        (
            self.intermediate_states,
            self.intermediate_state_tags,
            self.rotated_intermediate_states,
            self.rotated_intermediate_state_tags,
        ) = {
            let mut inv_outputs = Vec::with_capacity(10);
            let mut tags = Vec::with_capacity(10);
            let mut l_inv_outputs = Vec::with_capacity(10);
            let mut l_tags = Vec::with_capacity(10);

            let mut witness_state_chunks_it = witness_bytes[16 + 4 * 10..].chunks_exact(16);
            let mut tag_state_chunks_it = self.tags[128 + 32 * 10..].chunks_exact(128);
            assert_eq!(witness_state_chunks_it.len(), 10 + 1);
            assert_eq!(tag_state_chunks_it.len(), 10 + 1);

            for _ in 0..10 {
                let wit_chunk =
                    <&[u8] as TryInto<[u8; 16]>>::try_into(witness_state_chunks_it.next().unwrap())
                        .unwrap();
                let inv_out_state = wit_chunk.map(GF2p8);
                l_inv_outputs.push(inv_out_state.map(rotate_gf2p8));
                inv_outputs.push(inv_out_state);

                let tag_chunk = tag_state_chunks_it.next().unwrap();
                let mut inv_out_state_tag: [GF2p128; 16] = Default::default();
                let mut l_inv_out_state_tag: [GF2p128; 16] = Default::default();
                for j in 0..16 {
                    inv_out_state_tag[j] = lift_into_gf2p8_subfield(&tag_chunk[8 * j..8 * (j + 1)]);
                    l_inv_out_state_tag[j] =
                        lift_into_rotated_gf2p8_subfield(&tag_chunk[8 * j..8 * (j + 1)]);
                }
                tags.push(inv_out_state_tag);
                l_tags.push(l_inv_out_state_tag);
            }

            (inv_outputs, tags, l_inv_outputs, l_tags)
        };
    }

    #[allow(clippy::type_complexity)]
    pub fn get_lifted_commitments(
        &mut self,
    ) -> (
        (GF2p128, GF2p128),
        (&[[GF2p8; 4]], &[[GF2p128; 4]]),
        (&[[GF2p8; 16]], &[[GF2p128; 16]]),
        (&[[GF2p8; 16]], &[[GF2p128; 16]]),
        (&[[GF2p8; 16]], &[[GF2p128; 16]]),
    ) {
        (
            (self.mask, self.mask_tag),
            (
                &self.key_schedule_intermediate_states,
                &self.key_schedule_intermediate_state_tags,
            ),
            (&self.round_keys, &self.round_key_tags),
            (&self.intermediate_states, &self.intermediate_state_tags),
            (
                &self.rotated_intermediate_states,
                &self.rotated_intermediate_state_tags,
            ),
        )
    }

    #[allow(non_snake_case)]
    fn aggregates_conditions(&self, challenge_seed: [u8; 32]) -> (GF2p128, GF2p128) {
        let mut chi_rng = ChaChaRng::from_seed(challenge_seed);
        let mut A0 = self.mask_tag;
        let mut A1 = self.mask;

        // key schedule conditions
        for key_schedule_round_i in 0..10 {
            let inv_input = {
                let round_key = self.round_keys[key_schedule_round_i];
                [round_key[13], round_key[14], round_key[15], round_key[12]]
            };
            let inv_input_tag = {
                let round_key_tag = self.round_key_tags[key_schedule_round_i];
                [
                    round_key_tag[13],
                    round_key_tag[14],
                    round_key_tag[15],
                    round_key_tag[12],
                ]
            };
            let inv_output = self.key_schedule_intermediate_states[key_schedule_round_i];
            let inv_output_tag = self.key_schedule_intermediate_state_tags[key_schedule_round_i];
            for j in 0..4 {
                let chi_k = GF2p128::random(&mut chi_rng);
                A0 += chi_k * (inv_input_tag[j] * inv_output_tag[j]);
                A1 += chi_k
                    * (inv_input_tag[j] * GF2p128::embed_gf2p8(inv_output[j])
                        + inv_output_tag[j] * GF2p128::embed_gf2p8(inv_input[j]));
                assert_eq!(
                    inv_input[j] * inv_output[j],
                    GF2p8::ONE,
                    "condition failed in key schedule round {key_schedule_round_i}"
                );
            }
        }

        // first round
        {
            let inv_input = {
                let mut input = self.public_key.input.map(GF2p8);
                input
                    .iter_mut()
                    .zip(self.round_keys[0].iter())
                    .for_each(|(i, k)| *i += k);
                input
            };
            let inv_input_tag = self.round_key_tags[0];
            let inv_output = self.intermediate_states[0];
            let inv_output_tag = self.intermediate_state_tags[0];
            for j in 0..16 {
                let chi_k = GF2p128::random(&mut chi_rng);
                A0 += chi_k * (inv_input_tag[j] * inv_output_tag[j]);
                A1 += chi_k
                    * (inv_input_tag[j] * GF2p128::embed_gf2p8(inv_output[j])
                        + inv_output_tag[j] * GF2p128::embed_gf2p8(inv_input[j]));
                assert_eq!(inv_input[j] * inv_output[j], GF2p8::ONE);
            }
        }

        let Gamma = GF2p8(0b01100011);
        let emb_2 = GF2p128::embed_gf2p8(GF2p8(2));
        let emb_3 = GF2p128::embed_gf2p8(GF2p8(3));

        // middle rounds
        for sbox_layer_i in 1..10 {
            let inv_input = {
                let mut input = self.round_keys[sbox_layer_i];
                let Ls = self.rotated_intermediate_states[sbox_layer_i - 1];

                input[0] += GF2p8(2) * Ls[0] + GF2p8(3) * Ls[5] + Ls[10] + Ls[15] + Gamma;
                input[1] += Ls[0] + GF2p8(2) * Ls[5] + GF2p8(3) * Ls[10] + Ls[15] + Gamma;
                input[2] += Ls[0] + Ls[5] + GF2p8(2) * Ls[10] + GF2p8(3) * Ls[15] + Gamma;
                input[3] += GF2p8(3) * Ls[0] + Ls[5] + Ls[10] + GF2p8(2) * Ls[15] + Gamma;

                input[4] += GF2p8(2) * Ls[4] + GF2p8(3) * Ls[9] + Ls[14] + Ls[3] + Gamma;
                input[5] += Ls[4] + GF2p8(2) * Ls[9] + GF2p8(3) * Ls[14] + Ls[3] + Gamma;
                input[6] += Ls[4] + Ls[9] + GF2p8(2) * Ls[14] + GF2p8(3) * Ls[3] + Gamma;
                input[7] += GF2p8(3) * Ls[4] + Ls[9] + Ls[14] + GF2p8(2) * Ls[3] + Gamma;

                input[8] += GF2p8(2) * Ls[8] + GF2p8(3) * Ls[13] + Ls[2] + Ls[7] + Gamma;
                input[9] += Ls[8] + GF2p8(2) * Ls[13] + GF2p8(3) * Ls[2] + Ls[7] + Gamma;
                input[10] += Ls[8] + Ls[13] + GF2p8(2) * Ls[2] + GF2p8(3) * Ls[7] + Gamma;
                input[11] += GF2p8(3) * Ls[8] + Ls[13] + Ls[2] + GF2p8(2) * Ls[7] + Gamma;

                input[12] += GF2p8(2) * Ls[12] + GF2p8(3) * Ls[1] + Ls[6] + Ls[11] + Gamma;
                input[13] += Ls[12] + GF2p8(2) * Ls[1] + GF2p8(3) * Ls[6] + Ls[11] + Gamma;
                input[14] += Ls[12] + Ls[1] + GF2p8(2) * Ls[6] + GF2p8(3) * Ls[11] + Gamma;
                input[15] += GF2p8(3) * Ls[12] + Ls[1] + Ls[6] + GF2p8(2) * Ls[11] + Gamma;

                input
            };
            let inv_input_tag = {
                let mut tags = self.round_key_tags[sbox_layer_i];
                let Ls = self.rotated_intermediate_state_tags[sbox_layer_i - 1];

                tags[0] += emb_2 * Ls[0] + emb_3 * Ls[5] + Ls[10] + Ls[15];
                tags[1] += Ls[0] + emb_2 * Ls[5] + emb_3 * Ls[10] + Ls[15];
                tags[2] += Ls[0] + Ls[5] + emb_2 * Ls[10] + emb_3 * Ls[15];
                tags[3] += emb_3 * Ls[0] + Ls[5] + Ls[10] + emb_2 * Ls[15];

                tags[4] += emb_2 * Ls[4] + emb_3 * Ls[9] + Ls[14] + Ls[3];
                tags[5] += Ls[4] + emb_2 * Ls[9] + emb_3 * Ls[14] + Ls[3];
                tags[6] += Ls[4] + Ls[9] + emb_2 * Ls[14] + emb_3 * Ls[3];
                tags[7] += emb_3 * Ls[4] + Ls[9] + Ls[14] + emb_2 * Ls[3];

                tags[8] += emb_2 * Ls[8] + emb_3 * Ls[13] + Ls[2] + Ls[7];
                tags[9] += Ls[8] + emb_2 * Ls[13] + emb_3 * Ls[2] + Ls[7];
                tags[10] += Ls[8] + Ls[13] + emb_2 * Ls[2] + emb_3 * Ls[7];
                tags[11] += emb_3 * Ls[8] + Ls[13] + Ls[2] + emb_2 * Ls[7];

                tags[12] += emb_2 * Ls[12] + emb_3 * Ls[1] + Ls[6] + Ls[11];
                tags[13] += Ls[12] + emb_2 * Ls[1] + emb_3 * Ls[6] + Ls[11];
                tags[14] += Ls[12] + Ls[1] + emb_2 * Ls[6] + emb_3 * Ls[11];
                tags[15] += emb_3 * Ls[12] + Ls[1] + Ls[6] + emb_2 * Ls[11];

                tags
            };
            let inv_output = self.intermediate_states[sbox_layer_i];
            let inv_output_tag = self.intermediate_state_tags[sbox_layer_i];
            for j in 0..16 {
                let chi_k = GF2p128::random(&mut chi_rng);
                A0 += chi_k * (inv_input_tag[j] * inv_output_tag[j]);
                A1 += chi_k
                    * (inv_input_tag[j] * GF2p128::embed_gf2p8(inv_output[j])
                        + inv_output_tag[j] * GF2p128::embed_gf2p8(inv_input[j]));
                assert_eq!(
                    inv_input[j] * inv_output[j],
                    GF2p8::ONE,
                    "wrong values in sbox layer {sbox_layer_i}, byte {j}"
                );
            }
        }

        // last round
        {
            let input = {
                let mut input = self.round_keys[10];
                let Ls = self.rotated_intermediate_states[9];

                input[0] += Ls[0] + Gamma;
                input[1] += Ls[5] + Gamma;
                input[2] += Ls[10] + Gamma;
                input[3] += Ls[15] + Gamma;

                input[4] += Ls[4] + Gamma;
                input[5] += Ls[9] + Gamma;
                input[6] += Ls[14] + Gamma;
                input[7] += Ls[3] + Gamma;

                input[8] += Ls[8] + Gamma;
                input[9] += Ls[13] + Gamma;
                input[10] += Ls[2] + Gamma;
                input[11] += Ls[7] + Gamma;

                input[12] += Ls[12] + Gamma;
                input[13] += Ls[1] + Gamma;
                input[14] += Ls[6] + Gamma;
                input[15] += Ls[11] + Gamma;

                input
            };
            let output = self.public_key.output.map(GF2p8);
            assert_eq!(input, output);
            let input_tag = {
                let mut tags = self.round_key_tags[10];
                let Ls = self.rotated_intermediate_state_tags[9];

                tags[0] += Ls[0];
                tags[1] += Ls[5];
                tags[2] += Ls[10];
                tags[3] += Ls[15];

                tags[4] += Ls[4];
                tags[5] += Ls[9];
                tags[6] += Ls[14];
                tags[7] += Ls[3];

                tags[8] += Ls[8];
                tags[9] += Ls[13];
                tags[10] += Ls[2];
                tags[11] += Ls[7];

                tags[12] += Ls[12];
                tags[13] += Ls[1];
                tags[14] += Ls[6];
                tags[15] += Ls[11];

                tags
            };
            for input_tag_j in input_tag {
                let chi_k = GF2p128::random(&mut chi_rng);
                A1 += chi_k * input_tag_j;
            }
        }

        (A0, A1)
    }
}

pub struct FaestVerifierFromHC<HCR>
where
    HCR: HomComReceiver,
{
    public_key: PublicKey,
    hc_receiver: HCR,
    chi_seed: [u8; 32],
    proof: (GF2p128, GF2p128),
    global_key: GF2p128,
    keys: Vec<GF2p128>,
    mask_key: GF2p128,
    key_schedule_intermediate_state_keys: Vec<[GF2p128; 4]>,
    round_key_keys: Vec<[GF2p128; 16]>,
    intermediate_state_keys: Vec<[GF2p128; 16]>,
    rotated_intermediate_state_keys: Vec<[GF2p128; 16]>,
}

impl<HCR> Verifier for FaestVerifierFromHC<HCR>
where
    HCR: HomComReceiver<GlobalKey = GF2p128, Key = GF2p128>,
{
    type Commitment = HCR::Commitment;
    type Challenge1 = HCR::Challenge;
    type Response = HCR::Response;
    type Challenge2 = [u8; 32];
    type Proof = (GF2p128, GF2p128);
    type Choice = HCR::Choice;
    type Decommitment = HCR::Decommitment;

    fn new(public_key: PublicKey) -> Self {
        Self {
            public_key,
            hc_receiver: HCR::new(EXTENDED_WITNESS_BYTE_SIZE * 8 + 128),
            chi_seed: Default::default(),
            proof: Default::default(),
            key_schedule_intermediate_state_keys: Default::default(),
            global_key: Default::default(),
            keys: Default::default(),
            mask_key: Default::default(),
            round_key_keys: Default::default(),
            intermediate_state_keys: Default::default(),
            rotated_intermediate_state_keys: Default::default(),
        }
    }

    fn generate_commit_challenge_from_seed(seed: [u8; 32]) -> HCR::Challenge {
        HCR::generate_challenge_from_seed(seed)
    }

    fn commit_send_challenge_from_seed(
        &mut self,
        commitment: HCR::Commitment,
        seed: [u8; 32],
    ) -> HCR::Challenge {
        self.hc_receiver
            .commit_send_challenge_from_seed(commitment, seed)
    }

    fn commit_send_challenge(&mut self, commitment: HCR::Commitment) -> HCR::Challenge {
        self.hc_receiver.commit_send_challenge(commitment)
    }

    fn commit_receive_response(&mut self, response: HCR::Response) {
        self.hc_receiver.commit_receive_response(response);
    }

    fn generate_challenge_from_seed(seed: [u8; 32]) -> [u8; 32] {
        seed
    }

    fn send_challenge_from_seed(&mut self, seed: [u8; 32]) -> [u8; 32] {
        self.chi_seed = seed;
        seed
    }

    fn send_challenge(&mut self) -> [u8; 32] {
        self.chi_seed = thread_rng().gen();
        self.chi_seed
    }

    fn generate_choice_from_seed(seed: [u8; 32]) -> HCR::Choice {
        HCR::generate_choice_from_seed(seed)
    }

    fn choose_from_seed(&mut self, proof: Self::Proof, seed: [u8; 32]) -> HCR::Choice {
        self.proof = proof;
        self.hc_receiver.choose_from_seed(seed)
    }

    fn choose(&mut self, proof: Self::Proof) -> HCR::Choice {
        self.proof = proof;
        self.hc_receiver.choose()
    }

    #[allow(non_snake_case)]
    fn verify(&mut self, decommitment: HCR::Decommitment) -> bool {
        if let Some((global_key, keys)) = self.hc_receiver.receive(decommitment) {
            self.global_key = global_key;
            self.keys = keys;
        } else {
            return false;
        }
        self.lift_commitments();
        let B = self.aggregates_conditions();
        let (A0, A1) = self.proof;
        B == A0 + A1 * self.global_key
    }
}

impl<HCR> FaestVerifierFromHC<HCR>
where
    HCR: HomComReceiver<GlobalKey = GF2p128, Key = GF2p128>,
{
    #[allow(non_snake_case)]
    fn lift_commitments(&mut self) {
        // prepare uniform mask
        self.mask_key = {
            let bit_offset = 8 * EXTENDED_WITNESS_BYTE_SIZE;
            InnerProduct::<GF2p128>::inner_product(
                self.keys[bit_offset..].iter().copied(),
                (0..128).into_iter().map(|i| GF2p128::from_u128(1 << i)),
            )
        };

        self.round_key_keys = Vec::with_capacity(11);

        // handle the first round key
        {
            let key_chunk = &self.keys[0..128];
            let mut round_key_keys: [GF2p128; 16] = Default::default();
            for j in 0..16 {
                round_key_keys[j] = lift_into_gf2p8_subfield(&key_chunk[8 * j..8 * (j + 1)]);
            }
            self.round_key_keys.push(round_key_keys);
        }
        // handle the other round keys
        {
            // iterate through the witness and the tags in 32 bit chunks
            let mut key_word_chunks_it = self.keys[128..128 + 32 * 10].chunks_exact(32);
            assert_eq!(key_word_chunks_it.len(), 10);

            let Delta = self.global_key;
            let Delta_x_emb_Gamma = Delta * GF2p128::embed_gf2p8(GF2p8(0x63));

            for round_key_i in 1..11 {
                let key_chunk = key_word_chunks_it.next().unwrap();
                {
                    let mut inv_out = [GF2p128::ZERO; 4];
                    for j in 0..4 {
                        inv_out[j] = lift_into_gf2p8_subfield(&key_chunk[8 * j..8 * (j + 1)]);
                    }
                    self.key_schedule_intermediate_state_keys.push(inv_out);
                }
                let mut after_sbox = [GF2p128::ZERO; 16];
                for j in 0..4 {
                    after_sbox[j] =
                        lift_into_rotated_gf2p8_subfield(&key_chunk[8 * j..8 * (j + 1)])
                            - Delta_x_emb_Gamma;
                }
                let mut round_key_key = self.round_key_keys[round_key_i - 1];
                round_key_key[0] +=
                    after_sbox[0] - Delta * GF2p128::embed_gf2p8(ROUND_CONSTANTS[round_key_i - 1]);
                round_key_key[1] += after_sbox[1];
                round_key_key[2] += after_sbox[2];
                round_key_key[3] += after_sbox[3];
                for i in 1..4 {
                    round_key_key[i * 4] += round_key_key[(i - 1) * 4];
                    round_key_key[i * 4 + 1] += round_key_key[(i - 1) * 4 + 1];
                    round_key_key[i * 4 + 2] += round_key_key[(i - 1) * 4 + 2];
                    round_key_key[i * 4 + 3] += round_key_key[(i - 1) * 4 + 3];
                }
                self.round_key_keys.push(round_key_key);
            }
        }

        // prepare intermediate states
        (
            self.intermediate_state_keys,
            self.rotated_intermediate_state_keys,
        ) = {
            let mut keys = Vec::with_capacity(10);
            let mut l_keys = Vec::with_capacity(10);

            // iterate through the keys in 128 bit chunks
            let mut key_state_chunks_it = self.keys[128 + 32 * 10..].chunks_exact(128);
            assert_eq!(key_state_chunks_it.len(), 10 + 1);

            for _ in 0..10 {
                let key_chunk = key_state_chunks_it.next().unwrap();
                let mut inv_out_state_key: [GF2p128; 16] = Default::default();
                let mut l_inv_out_state_key: [GF2p128; 16] = Default::default();
                for j in 0..16 {
                    inv_out_state_key[j] = lift_into_gf2p8_subfield(&key_chunk[8 * j..8 * (j + 1)]);
                    l_inv_out_state_key[j] =
                        lift_into_rotated_gf2p8_subfield(&key_chunk[8 * j..8 * (j + 1)]);
                }
                keys.push(inv_out_state_key);
                l_keys.push(l_inv_out_state_key);
            }

            (keys, l_keys)
        };
    }

    #[allow(clippy::type_complexity)]
    pub fn get_lifted_commitments(
        &mut self,
    ) -> (
        GF2p128,
        GF2p128,
        &[[GF2p128; 4]],
        &[[GF2p128; 16]],
        &[[GF2p128; 16]],
        &[[GF2p128; 16]],
    ) {
        (
            self.global_key,
            self.mask_key,
            &self.key_schedule_intermediate_state_keys,
            &self.round_key_keys,
            &self.intermediate_state_keys,
            &self.rotated_intermediate_state_keys,
        )
    }

    #[allow(non_snake_case)]
    fn aggregates_conditions(&self) -> GF2p128 {
        let mut chi_rng = ChaChaRng::from_seed(self.chi_seed);
        let mut B = self.mask_key;

        let Delta = self.global_key;
        let Delta_squared = Delta.square();

        // key schedule conditions
        for key_schedule_round_i in 0..10 {
            let inv_input_key = {
                let round_key_key = self.round_key_keys[key_schedule_round_i];
                [
                    round_key_key[13],
                    round_key_key[14],
                    round_key_key[15],
                    round_key_key[12],
                ]
            };
            let inv_output_key = self.key_schedule_intermediate_state_keys[key_schedule_round_i];
            for j in 0..4 {
                let chi_k = GF2p128::random(&mut chi_rng);
                B += chi_k * ((inv_input_key[j] * inv_output_key[j]) - Delta_squared);
            }
        }

        // first round
        {
            let inv_input_key = {
                let input = &self.public_key.input;
                let mut key = self.round_key_keys[0];
                for j in 0..16 {
                    key[j] -= self.global_key * GF2p128::embed_gf2p8(GF2p8(input[j]));
                }
                key
            };
            let inv_output_key = self.intermediate_state_keys[0];
            for j in 0..16 {
                let chi_k = GF2p128::random(&mut chi_rng);
                B += chi_k * ((inv_input_key[j] * inv_output_key[j]) - Delta_squared);
            }
        }

        let Gamma = GF2p8(0b01100011);
        let emb_Gamma = GF2p128::embed_gf2p8(Gamma);
        let Delta_x_emb_Gamma = Delta * emb_Gamma;
        let emb_2 = GF2p128::embed_gf2p8(GF2p8(2));
        let emb_3 = GF2p128::embed_gf2p8(GF2p8(3));

        // middle rounds
        for sbox_layer_i in 1..10 {
            let inv_input_key = {
                let mut keys = self.round_key_keys[sbox_layer_i];
                let Ls = self.rotated_intermediate_state_keys[sbox_layer_i - 1];

                keys[0] += emb_2 * Ls[0] + emb_3 * Ls[5] + Ls[10] + Ls[15] - Delta_x_emb_Gamma;
                keys[1] += Ls[0] + emb_2 * Ls[5] + emb_3 * Ls[10] + Ls[15] - Delta_x_emb_Gamma;
                keys[2] += Ls[0] + Ls[5] + emb_2 * Ls[10] + emb_3 * Ls[15] - Delta_x_emb_Gamma;
                keys[3] += emb_3 * Ls[0] + Ls[5] + Ls[10] + emb_2 * Ls[15] - Delta_x_emb_Gamma;

                keys[4] += emb_2 * Ls[4] + emb_3 * Ls[9] + Ls[14] + Ls[3] - Delta_x_emb_Gamma;
                keys[5] += Ls[4] + emb_2 * Ls[9] + emb_3 * Ls[14] + Ls[3] - Delta_x_emb_Gamma;
                keys[6] += Ls[4] + Ls[9] + emb_2 * Ls[14] + emb_3 * Ls[3] - Delta_x_emb_Gamma;
                keys[7] += emb_3 * Ls[4] + Ls[9] + Ls[14] + emb_2 * Ls[3] - Delta_x_emb_Gamma;

                keys[8] += emb_2 * Ls[8] + emb_3 * Ls[13] + Ls[2] + Ls[7] - Delta_x_emb_Gamma;
                keys[9] += Ls[8] + emb_2 * Ls[13] + emb_3 * Ls[2] + Ls[7] - Delta_x_emb_Gamma;
                keys[10] += Ls[8] + Ls[13] + emb_2 * Ls[2] + emb_3 * Ls[7] - Delta_x_emb_Gamma;
                keys[11] += emb_3 * Ls[8] + Ls[13] + Ls[2] + emb_2 * Ls[7] - Delta_x_emb_Gamma;

                keys[12] += emb_2 * Ls[12] + emb_3 * Ls[1] + Ls[6] + Ls[11] - Delta_x_emb_Gamma;
                keys[13] += Ls[12] + emb_2 * Ls[1] + emb_3 * Ls[6] + Ls[11] - Delta_x_emb_Gamma;
                keys[14] += Ls[12] + Ls[1] + emb_2 * Ls[6] + emb_3 * Ls[11] - Delta_x_emb_Gamma;
                keys[15] += emb_3 * Ls[12] + Ls[1] + Ls[6] + emb_2 * Ls[11] - Delta_x_emb_Gamma;

                keys
            };
            let inv_output_key = self.intermediate_state_keys[sbox_layer_i];
            for j in 0..16 {
                let chi_k = GF2p128::random(&mut chi_rng);
                B += chi_k * ((inv_input_key[j] * inv_output_key[j]) - Delta_squared);
            }
        }

        // last round
        {
            let input_key = {
                let mut keys = self.round_key_keys[10];
                let Ls = self.rotated_intermediate_state_keys[9];

                keys[0] += Ls[0] - Delta_x_emb_Gamma;
                keys[1] += Ls[5] - Delta_x_emb_Gamma;
                keys[2] += Ls[10] - Delta_x_emb_Gamma;
                keys[3] += Ls[15] - Delta_x_emb_Gamma;

                keys[4] += Ls[4] - Delta_x_emb_Gamma;
                keys[5] += Ls[9] - Delta_x_emb_Gamma;
                keys[6] += Ls[14] - Delta_x_emb_Gamma;
                keys[7] += Ls[3] - Delta_x_emb_Gamma;

                keys[8] += Ls[8] - Delta_x_emb_Gamma;
                keys[9] += Ls[13] - Delta_x_emb_Gamma;
                keys[10] += Ls[2] - Delta_x_emb_Gamma;
                keys[11] += Ls[7] - Delta_x_emb_Gamma;

                keys[12] += Ls[12] - Delta_x_emb_Gamma;
                keys[13] += Ls[1] - Delta_x_emb_Gamma;
                keys[14] += Ls[6] - Delta_x_emb_Gamma;
                keys[15] += Ls[11] - Delta_x_emb_Gamma;

                keys
            };
            let output = self.public_key.output.map(GF2p8);
            for j in 0..16 {
                let chi_k = GF2p128::random(&mut chi_rng);
                B += chi_k
                    * (input_key[j] * Delta - GF2p128::embed_gf2p8(output[j]) * Delta_squared);
            }
        }

        B
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Block128;
    use crate::gf2psmall::{GF2p10, GF2p8, SmallGF};
    use crate::homcom::{HomCom128ReceiverFromVitH, HomCom128SenderFromVitH};
    use crate::primitives::{Aes128CtrLdPrg, Blake3LE};
    use crate::veccom::GgmVecCom;
    use crate::voleith::{VoleInTheHeadReceiverFromVC, VoleInTheHeadSenderFromVC};
    use itertools::izip;

    type VC = GgmVecCom<Block128, Aes128CtrLdPrg, blake3::Hasher, Blake3LE<Block128>>;
    type FaestProver<F> =
        FaestProverFromHC<HomCom128SenderFromVitH<VoleInTheHeadSenderFromVC<F, VC>>>;
    type FaestVerifier<F> =
        FaestVerifierFromHC<HomCom128ReceiverFromVitH<VoleInTheHeadReceiverFromVC<F, VC>>>;

    const SECRET_KEY_128: SecretKey = SecretKey::Aes128Key {
        key: [
            0x42, 0x13, 0x4f, 0x71, 0x34, 0x89, 0x1b, 0x16, 0x82, 0xa8, 0xab, 0x56, 0x76, 0x27,
            0x30, 0x0c,
        ],
    };
    const PUBLIC_KEY: PublicKey = PublicKey {
        input: [
            0x7b, 0x60, 0x66, 0xd6, 0x5a, 0x73, 0xd6, 0x00, 0xa0, 0xae, 0xf0, 0x1c, 0x8b, 0x19,
            0x17, 0x40,
        ],
        output: [
            0x89, 0x13, 0xaf, 0x07, 0x31, 0xb5, 0x81, 0xdf, 0x8a, 0x0a, 0xb5, 0x6e, 0x08, 0x3c,
            0x6b, 0x7c,
        ],
    };

    #[test]
    #[allow(non_snake_case)]
    fn test_lifting() {
        let n = 1;
        let Delta = GF2p128::random(thread_rng());
        let u_bytes: Vec<u8> = (0..n).map(|_| thread_rng().gen()).collect();
        let u_bits = GF2Vector::from_vec(u_bytes);
        let u_as_gf2p128: Vec<GF2p128> = u_bits
            .bits
            .iter()
            .by_vals()
            .map(|b| if b { GF2p128::ONE } else { GF2p128::ZERO })
            .collect();
        let v: Vec<GF2p128> = (0..8 * n).map(|_| GF2p128::random(thread_rng())).collect();
        let w: Vec<GF2p128> = u_as_gf2p128
            .iter()
            .zip(v.iter())
            .map(|(&u_i, &v_i)| u_i * Delta + v_i)
            .collect();

        let lifted_u: Vec<GF2p8> = u_bits.as_raw_slice().iter().copied().map(GF2p8).collect();
        let lifted_v: Vec<GF2p128> = v.chunks_exact(8).map(lift_into_gf2p8_subfield).collect();
        let lifted_w: Vec<GF2p128> = w.chunks_exact(8).map(lift_into_gf2p8_subfield).collect();
        assert_eq!(lifted_u.len(), n);
        assert_eq!(lifted_v.len(), n);
        assert_eq!(lifted_w.len(), n);
        for i in 0..n {
            assert_eq!(
                GF2p128::embed_gf2p8(lifted_u[i]) * Delta + lifted_v[i],
                lifted_w[i]
            );
        }
    }

    fn test_correctness<F: SmallGF>() {
        let mut prover = FaestProver::<F>::new(SECRET_KEY_128, PUBLIC_KEY);
        let mut verifier = FaestVerifier::<F>::new(PUBLIC_KEY);

        let commitment = prover.commit();
        let challenge = verifier.commit_send_challenge(commitment);
        let response = prover.commit_prove_consistency(challenge);
        verifier.commit_receive_response(response);

        let challenge = verifier.send_challenge();
        let proof = prover.prove(challenge);
        let choice = verifier.choose(proof);
        let decommitment = prover.transfer(choice);
        let result = verifier.verify(decommitment);

        // check correctness of commitments here
        {
            let (
                (mask, mask_tag),
                (key_schedule_intermediate_states, key_schedule_intermediate_state_tags),
                (round_keys, round_key_tags),
                (intermediate_states, intermediate_state_tags),
                (rotated_intermediate_states, rotated_intermediate_state_tags),
            ) = prover.get_lifted_commitments();
            let (
                global_key,
                mask_key,
                key_schedule_intermediate_state_keys,
                round_key_keys,
                intermediate_state_keys,
                rotated_intermediate_state_keys,
            ) = verifier.get_lifted_commitments();
            assert_eq!(mask_tag, mask * global_key + mask_key);
            for (is_i, is_t_i, is_k_i) in izip!(
                key_schedule_intermediate_states.iter().copied(),
                key_schedule_intermediate_state_tags.iter().copied(),
                key_schedule_intermediate_state_keys.iter().copied(),
            ) {
                for j in 0..4 {
                    assert_eq!(
                        is_t_i[j],
                        GF2p128::embed_gf2p8(is_i[j]) * global_key + is_k_i[j]
                    );
                }
            }
            for (rk_i, rk_t_i, rk_k_i) in izip!(
                round_keys.iter().copied(),
                round_key_tags.iter().copied(),
                round_key_keys.iter().copied(),
            ) {
                for j in 0..16 {
                    assert_eq!(
                        rk_t_i[j],
                        GF2p128::embed_gf2p8(rk_i[j]) * global_key + rk_k_i[j]
                    );
                }
            }
            for (is_i, is_t_i, is_k_i) in izip!(
                intermediate_states.iter().copied(),
                intermediate_state_tags.iter().copied(),
                intermediate_state_keys.iter().copied(),
            ) {
                for j in 0..16 {
                    assert_eq!(
                        is_t_i[j],
                        GF2p128::embed_gf2p8(is_i[j]) * global_key + is_k_i[j]
                    );
                }
            }
            for (ris_i, ris_t_i, ris_k_i) in izip!(
                rotated_intermediate_states.iter().copied(),
                rotated_intermediate_state_tags.iter().copied(),
                rotated_intermediate_state_keys.iter().copied(),
            ) {
                for j in 0..16 {
                    assert_eq!(
                        ris_t_i[j],
                        GF2p128::embed_gf2p8(ris_i[j]) * global_key + ris_k_i[j]
                    );
                }
            }
        }

        assert!(result, "verifier does not accept");
    }

    #[test]
    fn test_correctness_with_gf2p8() {
        test_correctness::<GF2p8>();
    }

    #[test]
    fn test_correctness_with_gf2p10() {
        test_correctness::<GF2p10>();
    }
}
