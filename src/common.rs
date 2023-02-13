use crate::traits::Block;
use rand::Rng;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, bincode::Encode)]
pub struct Block128(pub [u8; 16]);
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Block256(pub [u8; 32]);

impl AsRef<[u8]> for Block128 {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for Block128 {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl Block for Block128 {
    const BIT_LENGTH: usize = 128;
    fn random(rng: &mut impl Rng) -> Self {
        Self(rng.gen())
    }
}

impl AsRef<[u8]> for Block256 {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for Block256 {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl Block for Block256 {
    const BIT_LENGTH: usize = 256;
    fn random(rng: &mut impl Rng) -> Self {
        Self(rng.gen())
    }
}

pub fn prefix_decompose(x: usize, n_bits: u32) -> Vec<usize> {
    assert!(n_bits <= usize::BITS);
    assert!(n_bits == usize::BITS || x < 1 << n_bits);
    (0..n_bits).rev().map(|i| x >> i).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prefix_decompose() {
        assert_eq!(
            prefix_decompose(0b1100110100, 10),
            vec![
                0b1,
                0b11,
                0b110,
                0b1100,
                0b11001,
                0b110011,
                0b1100110,
                0b11001101,
                0b110011010,
                0b1100110100,
            ]
        );
        assert_eq!(
            prefix_decompose(0b0011101100, 10),
            vec![
                0b0,
                0b00,
                0b001,
                0b0011,
                0b00111,
                0b001110,
                0b0011101,
                0b00111011,
                0b001110110,
                0b0011101100,
            ]
        );
        assert_eq!(
            prefix_decompose(usize::MIN, usize::BITS),
            (0..usize::BITS).map(|_| 0).collect::<Vec<usize>>()
        );
        let mut prefixes = Vec::<usize>::with_capacity(usize::BITS as usize);
        prefixes.push(1);
        for i in 0..usize::BITS - 1 {
            prefixes.push((prefixes[i as usize] << 1) | 1)
        }
        assert_eq!(prefix_decompose(usize::MAX, usize::BITS), prefixes);
    }
}
