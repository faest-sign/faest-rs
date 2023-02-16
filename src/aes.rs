use crate::gf2psmall::GF2p8;
use ff::Field;
use itertools::izip;
use std::fmt;

pub const EXTENDED_WITNESS_BYTE_SIZE: usize = 16 + 10 * 4 + 16 * 10;

pub const ROUND_CONSTANTS: [GF2p8; 10] = [
    F(0x01),
    F(0x02),
    F(0x04),
    F(0x08),
    F(0x10),
    F(0x20),
    F(0x40),
    F(0x80),
    F(0x1b),
    F(0x36),
];

#[allow(non_snake_case)]
const fn F(x: u8) -> GF2p8 {
    GF2p8(x)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Aes {
    Aes128,
    Aes192,
    Aes256,
}

impl Aes {
    pub const fn get_key_size(self) -> usize {
        match self {
            Aes::Aes128 => 16,
            Aes::Aes192 => 24,
            Aes::Aes256 => 32,
        }
    }

    pub const fn get_words_per_key(self) -> usize {
        match self {
            Aes::Aes128 => 4,
            Aes::Aes192 => 6,
            Aes::Aes256 => 8,
        }
    }

    pub const fn get_number_rounds(self) -> usize {
        match self {
            Aes::Aes128 => 11,
            Aes::Aes192 => 13,
            Aes::Aes256 => 15,
        }
    }

    pub const fn get_number_sbox_layers_in_key_schedule(self) -> usize {
        match self {
            Aes::Aes128 => 10,
            Aes::Aes192 => 8,
            Aes::Aes256 => 13,
        }
    }

    pub const fn get_extended_witness_size(self) -> usize {
        self.get_key_size()
            + 4 * (self.get_number_sbox_layers_in_key_schedule())
            + 16 * (self.get_number_rounds() - 1)
    }

    pub fn compute_extended_witness(
        self,
        key: &[u8],
        input: &[u8; 16],
    ) -> Option<(Vec<u8>, [u8; 16])> {
        assert_eq!(key.len(), self.get_key_size());

        let mut buf = Vec::with_capacity(self.get_extended_witness_size());

        let key_as_gf2p8: Vec<_> = key.iter().copied().map(GF2p8).collect();

        let (round_keys, inv_outputs) =
            KeyExpansion::expand_and_collect_inv_outputs(self, &key_as_gf2p8);
        if inv_outputs
            .iter()
            .any(|x| x.iter().any(|y| *y == GF2p8::ZERO))
        {
            return None;
        }
        buf.extend(key);
        assert_eq!(buf.len(), self.get_key_size());
        buf.extend(inv_outputs.iter().flatten().copied().map(Into::<u8>::into));
        assert_eq!(
            buf.len(),
            self.get_key_size() + 4 * self.get_number_sbox_layers_in_key_schedule()
        );

        let input_state = AesState(input.map(GF2p8));
        let (output_state, inv_outputs) =
            input_state.encrypt_and_collect_inv_outputs(self, &round_keys);
        if inv_outputs
            .iter()
            .any(|x| x.0.iter().any(|y| *y == GF2p8::ZERO))
        {
            return None;
        }

        buf.extend(inv_outputs.iter().flat_map(|x| x.0).map(Into::<u8>::into));
        let output = output_state.0.map(Into::into);

        assert_eq!(buf.len(), self.get_extended_witness_size());

        Some((buf, output))
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub struct AesState(pub [GF2p8; 16]);

impl AesState {
    pub fn add_round_key(&self, round_key: [GF2p8; 16]) -> Self {
        let mut new = *self;
        new.0
            .iter_mut()
            .zip(round_key.iter().copied())
            .for_each(|(b, k)| *b += k);
        new
    }

    pub fn shift_rows(&self) -> Self {
        let old = &self.0;
        let new = [
            old[0], old[5], old[10], old[15], old[4], old[9], old[14], old[3], old[8], old[13],
            old[2], old[7], old[12], old[1], old[6], old[11],
        ];
        Self(new)
    }

    pub fn mix_columns(&self) -> Self {
        let old = &self.0;
        let mut new = [GF2p8::ZERO; 16];
        for i in 0..4 {
            new[4 * i] =
                F(2) * old[4 * i] + F(3) * old[4 * i + 1] + old[4 * i + 2] + old[4 * i + 3];
            new[4 * i + 1] =
                old[4 * i] + F(2) * old[4 * i + 1] + F(3) * old[4 * i + 2] + old[4 * i + 3];
            new[4 * i + 2] =
                old[4 * i] + old[4 * i + 1] + F(2) * old[4 * i + 2] + F(3) * old[4 * i + 3];
            new[4 * i + 3] =
                F(3) * old[4 * i] + old[4 * i + 1] + old[4 * i + 2] + F(2) * old[4 * i + 3];
        }
        Self(new)
    }

    pub fn sbox_invert_single(b: GF2p8) -> GF2p8 {
        if b == GF2p8::ZERO {
            b
        } else {
            b.invert().unwrap()
        }
    }

    pub fn sbox_rotations_single(b: GF2p8) -> GF2p8 {
        b + b.rotate_left(1) + b.rotate_left(2) + b.rotate_left(3) + b.rotate_left(4)
    }

    pub fn sbox_add_single(b: GF2p8) -> GF2p8 {
        b + GF2p8(0x63)
    }

    pub fn sbox_single(b: GF2p8) -> GF2p8 {
        Self::sbox_add_single(Self::sbox_rotations_single(Self::sbox_invert_single(b)))
    }

    pub fn sbox_invert(&self) -> Self {
        let mut new = *self;
        new.0
            .iter_mut()
            .for_each(|b| *b = Self::sbox_invert_single(*b));
        new
    }

    pub fn sbox_rotations(&self) -> Self {
        let mut new = *self;
        new.0
            .iter_mut()
            .for_each(|b| *b = Self::sbox_rotations_single(*b));
        new
    }

    pub fn sbox_add(&self) -> Self {
        let mut new = *self;
        new.0
            .iter_mut()
            .for_each(|b| *b = Self::sbox_add_single(*b));
        new
    }

    pub fn sub_bytes(&self) -> (Self, Self) {
        let inv_out = self.sbox_invert();
        (inv_out.sbox_rotations().sbox_add(), inv_out)
    }

    pub fn round(&self, round_key: [GF2p8; 16]) -> (Self, Self) {
        let (tmp, inv_out) = self.sub_bytes();
        (
            tmp.shift_rows().mix_columns().add_round_key(round_key),
            inv_out,
        )
    }
    pub fn last_round(&self, round_key: [GF2p8; 16]) -> (Self, Self) {
        let (tmp, inv_out) = self.sub_bytes();
        (tmp.shift_rows().add_round_key(round_key), inv_out)
    }

    pub fn encrypt(&self, variant: Aes, round_keys: &[[GF2p8; 16]]) -> Self {
        self.encrypt_and_collect_inv_outputs(variant, round_keys).0
    }

    #[allow(non_snake_case)]
    pub fn encrypt_and_collect_inv_outputs(
        &self,
        variant: Aes,
        round_keys: &[[GF2p8; 16]],
    ) -> (Self, Vec<Self>) {
        let R = variant.get_number_rounds();
        assert_eq!(round_keys.len(), R);
        let mut intermediate_states = vec![Self::default(); R - 1];
        let mut tmp = self.add_round_key(round_keys[0]);
        let mut inv_out;
        for (i, round_key_i) in round_keys.iter().take(R - 1).skip(1).enumerate() {
            (tmp, inv_out) = tmp.round(*round_key_i);
            intermediate_states[i] = inv_out;
        }
        (tmp, inv_out) = tmp.last_round(round_keys[10]);
        *intermediate_states.last_mut().unwrap() = inv_out;
        (tmp, intermediate_states)
    }
}

impl From<[GF2p8; 16]> for AesState {
    fn from(x: [GF2p8; 16]) -> Self {
        Self(x)
    }
}

impl From<[u8; 16]> for AesState {
    fn from(x: [u8; 16]) -> Self {
        let mut new = Self::default();
        for (i, b) in x.into_iter().enumerate() {
            new.0[i] = GF2p8(b);
        }
        new
    }
}

impl From<u128> for AesState {
    fn from(x: u128) -> Self {
        Self::from(x.to_be_bytes())
    }
}

impl From<AesState> for [GF2p8; 16] {
    fn from(val: AesState) -> Self {
        val.0
    }
}

impl From<AesState> for [u8; 16] {
    fn from(val: AesState) -> Self {
        let mut new = [0u8; 16];
        for (i, b) in val.0.into_iter().enumerate() {
            new[i] = b.0;
        }
        new
    }
}

impl From<AesState> for u128 {
    fn from(val: AesState) -> Self {
        u128::from_be_bytes(val.into())
    }
}

impl fmt::Debug for AesState {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        f.write_str("[")?;
        for i in (8..16).rev() {
            f.write_str(&format!("{:02x}", self.0[i].0))?;
        }
        f.write_str(" ")?;
        for i in (0..8).rev() {
            f.write_str(&format!("{:02x}", self.0[i].0))?;
        }
        // f.write_str(&format!("0x{:02x}", self.0[0].0))?;
        // for i in 1..16 {
        //     f.write_str(&format!(", 0x{:02x}", self.0[i].0))?;
        // }
        f.write_str("]")
    }
}

pub struct KeyExpansion {}

impl KeyExpansion {
    fn rot_word(w: &[GF2p8]) -> [GF2p8; 4] {
        assert_eq!(w.len(), 4);
        [w[1], w[2], w[3], w[0]]
    }

    fn sub_word(w: &[GF2p8]) -> ([GF2p8; 4], [GF2p8; 4]) {
        assert_eq!(w.len(), 4);
        let inv_out = [
            AesState::sbox_invert_single(w[0]),
            AesState::sbox_invert_single(w[1]),
            AesState::sbox_invert_single(w[2]),
            AesState::sbox_invert_single(w[3]),
        ];
        let out = [
            AesState::sbox_add_single(AesState::sbox_rotations_single(inv_out[0])),
            AesState::sbox_add_single(AesState::sbox_rotations_single(inv_out[1])),
            AesState::sbox_add_single(AesState::sbox_rotations_single(inv_out[2])),
            AesState::sbox_add_single(AesState::sbox_rotations_single(inv_out[3])),
        ];
        (out, inv_out)
    }

    #[allow(non_snake_case)]
    pub fn expand_and_collect_inv_outputs(
        variant: Aes,
        key: &[GF2p8],
    ) -> (Vec<[GF2p8; 16]>, Vec<[GF2p8; 4]>) {
        let R = variant.get_number_rounds();
        let N = variant.get_words_per_key();

        assert_eq!(key.len(), variant.get_key_size());

        let mut inv_outputs =
            Vec::<[GF2p8; 4]>::with_capacity(variant.get_number_sbox_layers_in_key_schedule());
        let mut key_words = [[GF2p8::ZERO; 4]; 4 * Aes::Aes256.get_number_rounds()];

        for i in 0..4 * R {
            if i < N {
                key_words[i].copy_from_slice(&key[i * 4..(i + 1) * 4]);
            } else if i >= N && (i % N) == 0 {
                let (mut new_word, inv_out) = Self::sub_word(&Self::rot_word(&key_words[i - 1]));
                inv_outputs.push(inv_out);
                new_word[0] += ROUND_CONSTANTS[(i - 1) / N];
                for j in 0..4 {
                    key_words[i][j] = key_words[i - N][j] + new_word[j]
                }
            } else if i >= N && N > 6 && (i % N) == 4 {
                let (new_word, inv_out) = Self::sub_word(&key_words[i - 1]);
                inv_outputs.push(inv_out);
                for j in 0..4 {
                    key_words[i][j] = key_words[i - N][j] + new_word[j]
                }
            } else {
                for j in 0..4 {
                    key_words[i][j] = key_words[i - N][j] + key_words[i - 1][j];
                }
            }
        }

        let mut round_keys = vec![[GF2p8::default(); 16]; R];
        for (rk, kws) in izip!(round_keys.iter_mut(), key_words.chunks_exact(4)) {
            rk[0..4].copy_from_slice(&kws[0]);
            rk[4..8].copy_from_slice(&kws[1]);
            rk[8..12].copy_from_slice(&kws[2]);
            rk[12..16].copy_from_slice(&kws[3]);
        }

        (round_keys, inv_outputs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sub_bytes() {
        let input = AesState::from(0x00102030405060708090a0b0c0d0e0f0u128);
        let expected_output = AesState::from(0x63cab7040953d051cd60e0e7ba70e18cu128);
        let (output, _) = input.sub_bytes();
        assert_eq!(output, expected_output);
    }

    #[test]
    fn test_mix_columns() {
        let input = AesState::from(0x6353e08c0960e104cd70b751bacad0e7u128);
        let expected_output = AesState::from(0x5f72641557f5bc92f7be3b291db9f91au128);
        let output = input.mix_columns();
        assert_eq!(output, expected_output);
    }

    #[test]
    fn test_shift_rows() {
        let input = AesState::from(0x63cab7040953d051cd60e0e7ba70e18cu128);
        let expected_output = AesState::from(0x6353e08c0960e104cd70b751bacad0e7u128);
        let output = input.shift_rows();
        assert_eq!(output, expected_output);
    }

    #[test]
    fn test_aes128_key_expand() {
        let key: [GF2p8; 16] = [
            F(0x2b),
            F(0x7e),
            F(0x15),
            F(0x16),
            F(0x28),
            F(0xae),
            F(0xd2),
            F(0xa6),
            F(0xab),
            F(0xf7),
            F(0x15),
            F(0x88),
            F(0x09),
            F(0xcf),
            F(0x4f),
            F(0x3c),
        ];
        let expected_round_keys: [[GF2p8; 16]; 11] = [
            [
                F(0x2b),
                F(0x7e),
                F(0x15),
                F(0x16),
                F(0x28),
                F(0xae),
                F(0xd2),
                F(0xa6),
                F(0xab),
                F(0xf7),
                F(0x15),
                F(0x88),
                F(0x09),
                F(0xcf),
                F(0x4f),
                F(0x3c),
            ],
            [
                F(0xa0),
                F(0xfa),
                F(0xfe),
                F(0x17),
                F(0x88),
                F(0x54),
                F(0x2c),
                F(0xb1),
                F(0x23),
                F(0xa3),
                F(0x39),
                F(0x39),
                F(0x2a),
                F(0x6c),
                F(0x76),
                F(0x05),
            ],
            [
                F(0xf2),
                F(0xc2),
                F(0x95),
                F(0xf2),
                F(0x7a),
                F(0x96),
                F(0xb9),
                F(0x43),
                F(0x59),
                F(0x35),
                F(0x80),
                F(0x7a),
                F(0x73),
                F(0x59),
                F(0xf6),
                F(0x7f),
            ],
            [
                F(0x3d),
                F(0x80),
                F(0x47),
                F(0x7d),
                F(0x47),
                F(0x16),
                F(0xfe),
                F(0x3e),
                F(0x1e),
                F(0x23),
                F(0x7e),
                F(0x44),
                F(0x6d),
                F(0x7a),
                F(0x88),
                F(0x3b),
            ],
            [
                F(0xef),
                F(0x44),
                F(0xa5),
                F(0x41),
                F(0xa8),
                F(0x52),
                F(0x5b),
                F(0x7f),
                F(0xb6),
                F(0x71),
                F(0x25),
                F(0x3b),
                F(0xdb),
                F(0x0b),
                F(0xad),
                F(0x00),
            ],
            [
                F(0xd4),
                F(0xd1),
                F(0xc6),
                F(0xf8),
                F(0x7c),
                F(0x83),
                F(0x9d),
                F(0x87),
                F(0xca),
                F(0xf2),
                F(0xb8),
                F(0xbc),
                F(0x11),
                F(0xf9),
                F(0x15),
                F(0xbc),
            ],
            [
                F(0x6d),
                F(0x88),
                F(0xa3),
                F(0x7a),
                F(0x11),
                F(0x0b),
                F(0x3e),
                F(0xfd),
                F(0xdb),
                F(0xf9),
                F(0x86),
                F(0x41),
                F(0xca),
                F(0x00),
                F(0x93),
                F(0xfd),
            ],
            [
                F(0x4e),
                F(0x54),
                F(0xf7),
                F(0x0e),
                F(0x5f),
                F(0x5f),
                F(0xc9),
                F(0xf3),
                F(0x84),
                F(0xa6),
                F(0x4f),
                F(0xb2),
                F(0x4e),
                F(0xa6),
                F(0xdc),
                F(0x4f),
            ],
            [
                F(0xea),
                F(0xd2),
                F(0x73),
                F(0x21),
                F(0xb5),
                F(0x8d),
                F(0xba),
                F(0xd2),
                F(0x31),
                F(0x2b),
                F(0xf5),
                F(0x60),
                F(0x7f),
                F(0x8d),
                F(0x29),
                F(0x2f),
            ],
            [
                F(0xac),
                F(0x77),
                F(0x66),
                F(0xf3),
                F(0x19),
                F(0xfa),
                F(0xdc),
                F(0x21),
                F(0x28),
                F(0xd1),
                F(0x29),
                F(0x41),
                F(0x57),
                F(0x5c),
                F(0x00),
                F(0x6e),
            ],
            [
                F(0xd0),
                F(0x14),
                F(0xf9),
                F(0xa8),
                F(0xc9),
                F(0xee),
                F(0x25),
                F(0x89),
                F(0xe1),
                F(0x3f),
                F(0x0c),
                F(0xc8),
                F(0xb6),
                F(0x63),
                F(0x0c),
                F(0xa6),
            ],
        ];
        let (expanded_key, _) = KeyExpansion::expand_and_collect_inv_outputs(Aes::Aes128, &key);
        for i in 0..11 {
            assert_eq!(
                expanded_key[i], expected_round_keys[i],
                "round key {i} incorrect"
            );
        }
    }

    #[test]
    fn test_aes192_key_expand() {
        let key: [GF2p8; Aes::Aes192.get_key_size()] = [
            F(0x8e),
            F(0x73),
            F(0xb0),
            F(0xf7),
            F(0xda),
            F(0x0e),
            F(0x64),
            F(0x52),
            F(0xc8),
            F(0x10),
            F(0xf3),
            F(0x2b),
            F(0x80),
            F(0x90),
            F(0x79),
            F(0xe5),
            F(0x62),
            F(0xf8),
            F(0xea),
            F(0xd2),
            F(0x52),
            F(0x2c),
            F(0x6b),
            F(0x7b),
        ];
        let expected_round_keys: [[GF2p8; 16]; Aes::Aes192.get_number_rounds()] = [
            [
                F(0x8e),
                F(0x73),
                F(0xb0),
                F(0xf7),
                F(0xda),
                F(0x0e),
                F(0x64),
                F(0x52),
                F(0xc8),
                F(0x10),
                F(0xf3),
                F(0x2b),
                F(0x80),
                F(0x90),
                F(0x79),
                F(0xe5),
            ],
            [
                F(0x62),
                F(0xf8),
                F(0xea),
                F(0xd2),
                F(0x52),
                F(0x2c),
                F(0x6b),
                F(0x7b),
                F(0xfe),
                F(0x0c),
                F(0x91),
                F(0xf7),
                F(0x24),
                F(0x02),
                F(0xf5),
                F(0xa5),
            ],
            [
                F(0xec),
                F(0x12),
                F(0x06),
                F(0x8e),
                F(0x6c),
                F(0x82),
                F(0x7f),
                F(0x6b),
                F(0x0e),
                F(0x7a),
                F(0x95),
                F(0xb9),
                F(0x5c),
                F(0x56),
                F(0xfe),
                F(0xc2),
            ],
            [
                F(0x4d),
                F(0xb7),
                F(0xb4),
                F(0xbd),
                F(0x69),
                F(0xb5),
                F(0x41),
                F(0x18),
                F(0x85),
                F(0xa7),
                F(0x47),
                F(0x96),
                F(0xe9),
                F(0x25),
                F(0x38),
                F(0xfd),
            ],
            [
                F(0xe7),
                F(0x5f),
                F(0xad),
                F(0x44),
                F(0xbb),
                F(0x09),
                F(0x53),
                F(0x86),
                F(0x48),
                F(0x5a),
                F(0xf0),
                F(0x57),
                F(0x21),
                F(0xef),
                F(0xb1),
                F(0x4f),
            ],
            [
                F(0xa4),
                F(0x48),
                F(0xf6),
                F(0xd9),
                F(0x4d),
                F(0x6d),
                F(0xce),
                F(0x24),
                F(0xaa),
                F(0x32),
                F(0x63),
                F(0x60),
                F(0x11),
                F(0x3b),
                F(0x30),
                F(0xe6),
            ],
            [
                F(0xa2),
                F(0x5e),
                F(0x7e),
                F(0xd5),
                F(0x83),
                F(0xb1),
                F(0xcf),
                F(0x9a),
                F(0x27),
                F(0xf9),
                F(0x39),
                F(0x43),
                F(0x6a),
                F(0x94),
                F(0xf7),
                F(0x67),
            ],
            [
                F(0xc0),
                F(0xa6),
                F(0x94),
                F(0x07),
                F(0xd1),
                F(0x9d),
                F(0xa4),
                F(0xe1),
                F(0xec),
                F(0x17),
                F(0x86),
                F(0xeb),
                F(0x6f),
                F(0xa6),
                F(0x49),
                F(0x71),
            ],
            [
                F(0x48),
                F(0x5f),
                F(0x70),
                F(0x32),
                F(0x22),
                F(0xcb),
                F(0x87),
                F(0x55),
                F(0xe2),
                F(0x6d),
                F(0x13),
                F(0x52),
                F(0x33),
                F(0xf0),
                F(0xb7),
                F(0xb3),
            ],
            [
                F(0x40),
                F(0xbe),
                F(0xeb),
                F(0x28),
                F(0x2f),
                F(0x18),
                F(0xa2),
                F(0x59),
                F(0x67),
                F(0x47),
                F(0xd2),
                F(0x6b),
                F(0x45),
                F(0x8c),
                F(0x55),
                F(0x3e),
            ],
            [
                F(0xa7),
                F(0xe1),
                F(0x46),
                F(0x6c),
                F(0x94),
                F(0x11),
                F(0xf1),
                F(0xdf),
                F(0x82),
                F(0x1f),
                F(0x75),
                F(0x0a),
                F(0xad),
                F(0x07),
                F(0xd7),
                F(0x53),
            ],
            [
                F(0xca),
                F(0x40),
                F(0x05),
                F(0x38),
                F(0x8f),
                F(0xcc),
                F(0x50),
                F(0x06),
                F(0x28),
                F(0x2d),
                F(0x16),
                F(0x6a),
                F(0xbc),
                F(0x3c),
                F(0xe7),
                F(0xb5),
            ],
            [
                F(0xe9),
                F(0x8b),
                F(0xa0),
                F(0x6f),
                F(0x44),
                F(0x8c),
                F(0x77),
                F(0x3c),
                F(0x8e),
                F(0xcc),
                F(0x72),
                F(0x04),
                F(0x01),
                F(0x00),
                F(0x22),
                F(0x02),
            ],
        ];
        let (expanded_key, _) = KeyExpansion::expand_and_collect_inv_outputs(Aes::Aes192, &key);
        for i in 0..Aes::Aes192.get_number_rounds() {
            assert_eq!(
                expanded_key[i], expected_round_keys[i],
                "round key {i} incorrect"
            );
        }
    }

    #[test]
    fn test_aes256_key_expand() {
        let key: [GF2p8; Aes::Aes256.get_key_size()] = [
            F(0x60),
            F(0x3d),
            F(0xeb),
            F(0x10),
            F(0x15),
            F(0xca),
            F(0x71),
            F(0xbe),
            F(0x2b),
            F(0x73),
            F(0xae),
            F(0xf0),
            F(0x85),
            F(0x7d),
            F(0x77),
            F(0x81),
            F(0x1f),
            F(0x35),
            F(0x2c),
            F(0x07),
            F(0x3b),
            F(0x61),
            F(0x08),
            F(0xd7),
            F(0x2d),
            F(0x98),
            F(0x10),
            F(0xa3),
            F(0x09),
            F(0x14),
            F(0xdf),
            F(0xf4),
        ];
        let expected_round_keys: [[GF2p8; 16]; Aes::Aes256.get_number_rounds()] = [
            [
                F(0x60),
                F(0x3d),
                F(0xeb),
                F(0x10),
                F(0x15),
                F(0xca),
                F(0x71),
                F(0xbe),
                F(0x2b),
                F(0x73),
                F(0xae),
                F(0xf0),
                F(0x85),
                F(0x7d),
                F(0x77),
                F(0x81),
            ],
            [
                F(0x1f),
                F(0x35),
                F(0x2c),
                F(0x07),
                F(0x3b),
                F(0x61),
                F(0x08),
                F(0xd7),
                F(0x2d),
                F(0x98),
                F(0x10),
                F(0xa3),
                F(0x09),
                F(0x14),
                F(0xdf),
                F(0xf4),
            ],
            [
                F(0x9b),
                F(0xa3),
                F(0x54),
                F(0x11),
                F(0x8e),
                F(0x69),
                F(0x25),
                F(0xaf),
                F(0xa5),
                F(0x1a),
                F(0x8b),
                F(0x5f),
                F(0x20),
                F(0x67),
                F(0xfc),
                F(0xde),
            ],
            [
                F(0xa8),
                F(0xb0),
                F(0x9c),
                F(0x1a),
                F(0x93),
                F(0xd1),
                F(0x94),
                F(0xcd),
                F(0xbe),
                F(0x49),
                F(0x84),
                F(0x6e),
                F(0xb7),
                F(0x5d),
                F(0x5b),
                F(0x9a),
            ],
            [
                F(0xd5),
                F(0x9a),
                F(0xec),
                F(0xb8),
                F(0x5b),
                F(0xf3),
                F(0xc9),
                F(0x17),
                F(0xfe),
                F(0xe9),
                F(0x42),
                F(0x48),
                F(0xde),
                F(0x8e),
                F(0xbe),
                F(0x96),
            ],
            [
                F(0xb5),
                F(0xa9),
                F(0x32),
                F(0x8a),
                F(0x26),
                F(0x78),
                F(0xa6),
                F(0x47),
                F(0x98),
                F(0x31),
                F(0x22),
                F(0x29),
                F(0x2f),
                F(0x6c),
                F(0x79),
                F(0xb3),
            ],
            [
                F(0x81),
                F(0x2c),
                F(0x81),
                F(0xad),
                F(0xda),
                F(0xdf),
                F(0x48),
                F(0xba),
                F(0x24),
                F(0x36),
                F(0x0a),
                F(0xf2),
                F(0xfa),
                F(0xb8),
                F(0xb4),
                F(0x64),
            ],
            [
                F(0x98),
                F(0xc5),
                F(0xbf),
                F(0xc9),
                F(0xbe),
                F(0xbd),
                F(0x19),
                F(0x8e),
                F(0x26),
                F(0x8c),
                F(0x3b),
                F(0xa7),
                F(0x09),
                F(0xe0),
                F(0x42),
                F(0x14),
            ],
            [
                F(0x68),
                F(0x00),
                F(0x7b),
                F(0xac),
                F(0xb2),
                F(0xdf),
                F(0x33),
                F(0x16),
                F(0x96),
                F(0xe9),
                F(0x39),
                F(0xe4),
                F(0x6c),
                F(0x51),
                F(0x8d),
                F(0x80),
            ],
            [
                F(0xc8),
                F(0x14),
                F(0xe2),
                F(0x04),
                F(0x76),
                F(0xa9),
                F(0xfb),
                F(0x8a),
                F(0x50),
                F(0x25),
                F(0xc0),
                F(0x2d),
                F(0x59),
                F(0xc5),
                F(0x82),
                F(0x39),
            ],
            [
                F(0xde),
                F(0x13),
                F(0x69),
                F(0x67),
                F(0x6c),
                F(0xcc),
                F(0x5a),
                F(0x71),
                F(0xfa),
                F(0x25),
                F(0x63),
                F(0x95),
                F(0x96),
                F(0x74),
                F(0xee),
                F(0x15),
            ],
            [
                F(0x58),
                F(0x86),
                F(0xca),
                F(0x5d),
                F(0x2e),
                F(0x2f),
                F(0x31),
                F(0xd7),
                F(0x7e),
                F(0x0a),
                F(0xf1),
                F(0xfa),
                F(0x27),
                F(0xcf),
                F(0x73),
                F(0xc3),
            ],
            [
                F(0x74),
                F(0x9c),
                F(0x47),
                F(0xab),
                F(0x18),
                F(0x50),
                F(0x1d),
                F(0xda),
                F(0xe2),
                F(0x75),
                F(0x7e),
                F(0x4f),
                F(0x74),
                F(0x01),
                F(0x90),
                F(0x5a),
            ],
            [
                F(0xca),
                F(0xfa),
                F(0xaa),
                F(0xe3),
                F(0xe4),
                F(0xd5),
                F(0x9b),
                F(0x34),
                F(0x9a),
                F(0xdf),
                F(0x6a),
                F(0xce),
                F(0xbd),
                F(0x10),
                F(0x19),
                F(0x0d),
            ],
            [
                F(0xfe),
                F(0x48),
                F(0x90),
                F(0xd1),
                F(0xe6),
                F(0x18),
                F(0x8d),
                F(0x0b),
                F(0x04),
                F(0x6d),
                F(0xf3),
                F(0x44),
                F(0x70),
                F(0x6c),
                F(0x63),
                F(0x1e),
            ],
        ];
        let (expanded_key, _) = KeyExpansion::expand_and_collect_inv_outputs(Aes::Aes256, &key);
        for i in 0..Aes::Aes256.get_number_rounds() {
            assert_eq!(
                expanded_key[i], expected_round_keys[i],
                "round key {i} incorrect"
            );
        }
    }

    #[test]
    fn test_aes128_encrypt_block() {
        let key = [
            F(0x00),
            F(0x01),
            F(0x02),
            F(0x03),
            F(0x04),
            F(0x05),
            F(0x06),
            F(0x07),
            F(0x08),
            F(0x09),
            F(0x0a),
            F(0x0b),
            F(0x0c),
            F(0x0d),
            F(0x0e),
            F(0x0f),
        ];
        let input_block = [
            F(0x00),
            F(0x11),
            F(0x22),
            F(0x33),
            F(0x44),
            F(0x55),
            F(0x66),
            F(0x77),
            F(0x88),
            F(0x99),
            F(0xaa),
            F(0xbb),
            F(0xcc),
            F(0xdd),
            F(0xee),
            F(0xff),
        ];
        let expected_output_block = [
            F(0x69),
            F(0xc4),
            F(0xe0),
            F(0xd8),
            F(0x6a),
            F(0x7b),
            F(0x04),
            F(0x30),
            F(0xd8),
            F(0xcd),
            F(0xb7),
            F(0x80),
            F(0x70),
            F(0xb4),
            F(0xc5),
            F(0x5a),
        ];
        let state = AesState::from(input_block);
        let (round_keys, _) = KeyExpansion::expand_and_collect_inv_outputs(Aes::Aes128, &key);
        let output_block = state.encrypt(Aes::Aes128, &round_keys);
        assert_eq!(output_block.0, expected_output_block);
    }

    #[test]
    fn test_aes128_compute_extended_witness() {
        let key = [
            0x42, 0x13, 0x4f, 0x71, 0x34, 0x89, 0x1b, 0x16, 0x82, 0xa8, 0xab, 0x56, 0x76, 0x27,
            0x30, 0x0c,
        ];
        let input = [
            0x7b, 0x60, 0x66, 0xd6, 0x5a, 0x73, 0xd6, 0x00, 0xa0, 0xae, 0xf0, 0x1c, 0x8b, 0x19,
            0x17, 0x40,
        ];
        let output = [
            0x89, 0x13, 0xaf, 0x07, 0x31, 0xb5, 0x81, 0xdf, 0x8a, 0x0a, 0xb5, 0x6e, 0x08, 0x3c,
            0x6b, 0x7c,
        ];
        let result = Aes::Aes128.compute_extended_witness(&key, &input);
        assert!(result.is_some());
        let (_witness, computed_output) = result.unwrap();
        assert_eq!(computed_output, output);
    }
}
