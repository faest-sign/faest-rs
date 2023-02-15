use crate::arithmetic::{clmul_u8, clmul_u8_x86};
use crate::field::{BytesRepr, InnerProduct, LazyField, LazyMul, UnreducedField};
use core::arch::x86_64::*;
use core::iter::{Product, Sum};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use ff::Field;
use num_traits::identities::{One, Zero};
use rand::{Rng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, bincode::Encode)]
pub struct GF2p8(pub u8);

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct UnreducedGF2p8(pub u16);

impl GF2p8 {
    // x^8 + x^4 + x^3 + x + 1
    const POLYNOMIAL: u16 = 0x11b;
    pub const ORDER: usize = 256;
    pub const LOG_ORDER: u32 = 8;

    const LOG_TABLE: [u8; 256] = [
        0x00, 0x00, 0x19, 0x01, 0x32, 0x02, 0x1a, 0xc6, 0x4b, 0xc7, 0x1b, 0x68, 0x33, 0xee, 0xdf,
        0x03, 0x64, 0x04, 0xe0, 0x0e, 0x34, 0x8d, 0x81, 0xef, 0x4c, 0x71, 0x08, 0xc8, 0xf8, 0x69,
        0x1c, 0xc1, 0x7d, 0xc2, 0x1d, 0xb5, 0xf9, 0xb9, 0x27, 0x6a, 0x4d, 0xe4, 0xa6, 0x72, 0x9a,
        0xc9, 0x09, 0x78, 0x65, 0x2f, 0x8a, 0x05, 0x21, 0x0f, 0xe1, 0x24, 0x12, 0xf0, 0x82, 0x45,
        0x35, 0x93, 0xda, 0x8e, 0x96, 0x8f, 0xdb, 0xbd, 0x36, 0xd0, 0xce, 0x94, 0x13, 0x5c, 0xd2,
        0xf1, 0x40, 0x46, 0x83, 0x38, 0x66, 0xdd, 0xfd, 0x30, 0xbf, 0x06, 0x8b, 0x62, 0xb3, 0x25,
        0xe2, 0x98, 0x22, 0x88, 0x91, 0x10, 0x7e, 0x6e, 0x48, 0xc3, 0xa3, 0xb6, 0x1e, 0x42, 0x3a,
        0x6b, 0x28, 0x54, 0xfa, 0x85, 0x3d, 0xba, 0x2b, 0x79, 0x0a, 0x15, 0x9b, 0x9f, 0x5e, 0xca,
        0x4e, 0xd4, 0xac, 0xe5, 0xf3, 0x73, 0xa7, 0x57, 0xaf, 0x58, 0xa8, 0x50, 0xf4, 0xea, 0xd6,
        0x74, 0x4f, 0xae, 0xe9, 0xd5, 0xe7, 0xe6, 0xad, 0xe8, 0x2c, 0xd7, 0x75, 0x7a, 0xeb, 0x16,
        0x0b, 0xf5, 0x59, 0xcb, 0x5f, 0xb0, 0x9c, 0xa9, 0x51, 0xa0, 0x7f, 0x0c, 0xf6, 0x6f, 0x17,
        0xc4, 0x49, 0xec, 0xd8, 0x43, 0x1f, 0x2d, 0xa4, 0x76, 0x7b, 0xb7, 0xcc, 0xbb, 0x3e, 0x5a,
        0xfb, 0x60, 0xb1, 0x86, 0x3b, 0x52, 0xa1, 0x6c, 0xaa, 0x55, 0x29, 0x9d, 0x97, 0xb2, 0x87,
        0x90, 0x61, 0xbe, 0xdc, 0xfc, 0xbc, 0x95, 0xcf, 0xcd, 0x37, 0x3f, 0x5b, 0xd1, 0x53, 0x39,
        0x84, 0x3c, 0x41, 0xa2, 0x6d, 0x47, 0x14, 0x2a, 0x9e, 0x5d, 0x56, 0xf2, 0xd3, 0xab, 0x44,
        0x11, 0x92, 0xd9, 0x23, 0x20, 0x2e, 0x89, 0xb4, 0x7c, 0xb8, 0x26, 0x77, 0x99, 0xe3, 0xa5,
        0x67, 0x4a, 0xed, 0xde, 0xc5, 0x31, 0xfe, 0x18, 0x0d, 0x63, 0x8c, 0x80, 0xc0, 0xf7, 0x70,
        0x07,
    ];
    const ANTI_LOG_TABLE: [Self; 255] = [
        Self(0x01),
        Self(0x03),
        Self(0x05),
        Self(0x0f),
        Self(0x11),
        Self(0x33),
        Self(0x55),
        Self(0xff),
        Self(0x1a),
        Self(0x2e),
        Self(0x72),
        Self(0x96),
        Self(0xa1),
        Self(0xf8),
        Self(0x13),
        Self(0x35),
        Self(0x5f),
        Self(0xe1),
        Self(0x38),
        Self(0x48),
        Self(0xd8),
        Self(0x73),
        Self(0x95),
        Self(0xa4),
        Self(0xf7),
        Self(0x02),
        Self(0x06),
        Self(0x0a),
        Self(0x1e),
        Self(0x22),
        Self(0x66),
        Self(0xaa),
        Self(0xe5),
        Self(0x34),
        Self(0x5c),
        Self(0xe4),
        Self(0x37),
        Self(0x59),
        Self(0xeb),
        Self(0x26),
        Self(0x6a),
        Self(0xbe),
        Self(0xd9),
        Self(0x70),
        Self(0x90),
        Self(0xab),
        Self(0xe6),
        Self(0x31),
        Self(0x53),
        Self(0xf5),
        Self(0x04),
        Self(0x0c),
        Self(0x14),
        Self(0x3c),
        Self(0x44),
        Self(0xcc),
        Self(0x4f),
        Self(0xd1),
        Self(0x68),
        Self(0xb8),
        Self(0xd3),
        Self(0x6e),
        Self(0xb2),
        Self(0xcd),
        Self(0x4c),
        Self(0xd4),
        Self(0x67),
        Self(0xa9),
        Self(0xe0),
        Self(0x3b),
        Self(0x4d),
        Self(0xd7),
        Self(0x62),
        Self(0xa6),
        Self(0xf1),
        Self(0x08),
        Self(0x18),
        Self(0x28),
        Self(0x78),
        Self(0x88),
        Self(0x83),
        Self(0x9e),
        Self(0xb9),
        Self(0xd0),
        Self(0x6b),
        Self(0xbd),
        Self(0xdc),
        Self(0x7f),
        Self(0x81),
        Self(0x98),
        Self(0xb3),
        Self(0xce),
        Self(0x49),
        Self(0xdb),
        Self(0x76),
        Self(0x9a),
        Self(0xb5),
        Self(0xc4),
        Self(0x57),
        Self(0xf9),
        Self(0x10),
        Self(0x30),
        Self(0x50),
        Self(0xf0),
        Self(0x0b),
        Self(0x1d),
        Self(0x27),
        Self(0x69),
        Self(0xbb),
        Self(0xd6),
        Self(0x61),
        Self(0xa3),
        Self(0xfe),
        Self(0x19),
        Self(0x2b),
        Self(0x7d),
        Self(0x87),
        Self(0x92),
        Self(0xad),
        Self(0xec),
        Self(0x2f),
        Self(0x71),
        Self(0x93),
        Self(0xae),
        Self(0xe9),
        Self(0x20),
        Self(0x60),
        Self(0xa0),
        Self(0xfb),
        Self(0x16),
        Self(0x3a),
        Self(0x4e),
        Self(0xd2),
        Self(0x6d),
        Self(0xb7),
        Self(0xc2),
        Self(0x5d),
        Self(0xe7),
        Self(0x32),
        Self(0x56),
        Self(0xfa),
        Self(0x15),
        Self(0x3f),
        Self(0x41),
        Self(0xc3),
        Self(0x5e),
        Self(0xe2),
        Self(0x3d),
        Self(0x47),
        Self(0xc9),
        Self(0x40),
        Self(0xc0),
        Self(0x5b),
        Self(0xed),
        Self(0x2c),
        Self(0x74),
        Self(0x9c),
        Self(0xbf),
        Self(0xda),
        Self(0x75),
        Self(0x9f),
        Self(0xba),
        Self(0xd5),
        Self(0x64),
        Self(0xac),
        Self(0xef),
        Self(0x2a),
        Self(0x7e),
        Self(0x82),
        Self(0x9d),
        Self(0xbc),
        Self(0xdf),
        Self(0x7a),
        Self(0x8e),
        Self(0x89),
        Self(0x80),
        Self(0x9b),
        Self(0xb6),
        Self(0xc1),
        Self(0x58),
        Self(0xe8),
        Self(0x23),
        Self(0x65),
        Self(0xaf),
        Self(0xea),
        Self(0x25),
        Self(0x6f),
        Self(0xb1),
        Self(0xc8),
        Self(0x43),
        Self(0xc5),
        Self(0x54),
        Self(0xfc),
        Self(0x1f),
        Self(0x21),
        Self(0x63),
        Self(0xa5),
        Self(0xf4),
        Self(0x07),
        Self(0x09),
        Self(0x1b),
        Self(0x2d),
        Self(0x77),
        Self(0x99),
        Self(0xb0),
        Self(0xcb),
        Self(0x46),
        Self(0xca),
        Self(0x45),
        Self(0xcf),
        Self(0x4a),
        Self(0xde),
        Self(0x79),
        Self(0x8b),
        Self(0x86),
        Self(0x91),
        Self(0xa8),
        Self(0xe3),
        Self(0x3e),
        Self(0x42),
        Self(0xc6),
        Self(0x51),
        Self(0xf3),
        Self(0x0e),
        Self(0x12),
        Self(0x36),
        Self(0x5a),
        Self(0xee),
        Self(0x29),
        Self(0x7b),
        Self(0x8d),
        Self(0x8c),
        Self(0x8f),
        Self(0x8a),
        Self(0x85),
        Self(0x94),
        Self(0xa7),
        Self(0xf2),
        Self(0x0d),
        Self(0x17),
        Self(0x39),
        Self(0x4b),
        Self(0xdd),
        Self(0x7c),
        Self(0x84),
        Self(0x97),
        Self(0xa2),
        Self(0xfd),
        Self(0x1c),
        Self(0x24),
        Self(0x6c),
        Self(0xb4),
        Self(0xc7),
        Self(0x52),
        Self(0xf6),
    ];

    #[inline(always)]
    fn mul_reduce(x: u8, y: u8) -> u8 {
        unsafe {
            let x = _mm_set_epi64x(x as i64, y as i64);
            let x = _mm_clmulepi64_si128(x, x, 0b01);
            let c = _mm_srli_epi16(x, 8);
            let z = _mm_xor_si128(
                _mm_xor_si128(c, _mm_srli_epi16(c, 4)),
                _mm_xor_si128(_mm_srli_epi16(c, 5), _mm_srli_epi16(c, 7)),
            );
            let z_xor_x = _mm_xor_si128(z, x);
            let p_xor_x = _mm_xor_si128(
                _mm_xor_si128(_mm_slli_epi16(z, 4), _mm_slli_epi16(z, 3)),
                _mm_xor_si128(_mm_slli_epi16(z, 1), z_xor_x),
            );
            _mm_extract_epi32(p_xor_x, 0) as u8
        }
    }

    #[inline(always)]
    // fn reduce(x: u16) -> u8 {
    fn reduce(x: u32) -> u8 {
        // let lo = x as u8;
        // let c = (x >> 8) as u8;
        let lo = x;
        let c = x >> 8;
        let z = c ^ (c >> 4) ^ (c >> 5) ^ (c >> 7);
        let p = (z << 4) ^ (z << 3) ^ (z << 1) ^ z;
        (lo ^ p) as u8
    }

    pub fn rotate_left(self, n: u32) -> Self {
        Self(self.0.rotate_left(n))
    }
}

impl BytesRepr for GF2p8 {
    type Repr = [u8; 1];
    #[inline(always)]
    fn to_repr(self) -> Self::Repr {
        [self.0]
    }
    #[inline(always)]
    fn from_repr(repr: Self::Repr) -> Self {
        Self(repr[0])
    }
}

impl LazyField for GF2p8 {
    type UF = UnreducedGF2p8;
}

impl UnreducedField for UnreducedGF2p8 {
    type ReducedField = GF2p8;
    const ZERO: Self = Self(0);
    #[inline(always)]
    fn reduce(self) -> Self::ReducedField {
        // GF2p8(GF2p8::reduce(self.0))
        GF2p8(GF2p8::reduce(self.0 as u32))
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 ^ other.0)
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for UnreducedGF2p8 {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 ^ other.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add<<GF2p8 as LazyField>::UF> for GF2p8 {
    type Output = <GF2p8 as LazyField>::UF;
    #[inline(always)]
    fn add(self, other: <GF2p8 as LazyField>::UF) -> Self::Output {
        UnreducedGF2p8(self.0 as u16 ^ other.0)
    }
}

impl LazyMul for GF2p8 {
    type Output = UnreducedGF2p8;
    #[inline(always)]
    fn lazy_mul(self, other: Self) -> Self::Output {
        // UnreducedGF2p8(clmul_u8(self.0, other.0))
        // UnreducedGF2p8(clmul_u8_x86(self.0, other.0))
        UnreducedGF2p8(clmul_u8_x86(self.0, other.0) as u16)
    }
}

impl Mul for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: Self) -> Self::Output {
        // 1) not that great in micro benchmark, but ok overall
        // Self(Self::reduce(clmul_u8(self.0, other.0) as u32))
        // 2) good in micro benchmark, but bad overall?
        // Self(Self::reduce(clmul_u8_x86(self.0, other.0)))
        // 3) good in micro benchmark, but bad overall?
        // Self(Self::mul_reduce(self.0, other.0))
        // 4) not constant time, but decent (and currently best) overall
        if self.0 == 0 || other.0 == 0 {
            Self(0)
        } else {
            Self::ANTI_LOG_TABLE[((Self::LOG_TABLE[self.0 as usize] as u32
                + Self::LOG_TABLE[other.0 as usize] as u32)
                % 255) as usize]
        }
    }
}
impl Mul<bool> for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: bool) -> Self::Output {
        Self(self.0 * other as u8)
    }
}
impl MulAssign<bool> for GF2p8 {
    #[inline(always)]
    fn mul_assign(&mut self, other: bool) {
        self.0 *= other as u8;
    }
}

impl ConstantTimeEq for GF2p8 {
    #[inline(always)]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}
impl ConditionallySelectable for GF2p8 {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self(u8::conditional_select(&a.0, &b.0, choice))
    }
}

impl_additional_field_arithmetic!(GF2p8, UnreducedGF2p8);

impl Field for GF2p8 {
    const ZERO: Self = Self(0u8);
    const ONE: Self = Self(1u8);

    #[inline(always)]
    fn random(mut rng: impl RngCore) -> Self {
        Self(rng.gen())
    }

    #[inline(always)]
    fn square(&self) -> Self {
        *self * self
    }

    #[inline(always)]
    fn double(&self) -> Self {
        *self + self
    }

    #[inline(always)]
    fn invert(&self) -> CtOption<Self> {
        let mut x = *self;
        let mut y = Self::ONE;
        x = x.square();
        for _ in 0..6 {
            y = x * y;
            x = x.square();
        }
        CtOption::new(x * y, !self.ct_eq(&Self::ZERO))
    }

    fn sqrt_ratio(_num: &Self, _div: &Self) -> (Choice, Self) {
        panic!("not implemented")
    }
}

impl From<u8> for GF2p8 {
    fn from(val: u8) -> Self {
        Self(val)
    }
}

impl From<GF2p8> for u8 {
    fn from(val: GF2p8) -> u8 {
        val.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    const GF2P8_VALUES: [GF2p8; 4] = [GF2p8(0x4b), GF2p8(0x27), GF2p8(0xb1), GF2p8(0xfc)];

    #[test]
    fn test_gf2p8_add() {
        for _ in 0..1000 {
            let val1 = GF2p8::random(thread_rng());
            assert_eq!(val1 + GF2p8::ZERO, val1);
            assert_eq!(val1 + val1, GF2p8::ZERO);
            let val2 = GF2p8::random(thread_rng());
            assert_eq!(val1 + val2, val2 + val1);
        }

        let values = &GF2P8_VALUES;
        assert_eq!(values[0] + values[0], GF2p8(0x00));
        assert_eq!(values[0] + values[1], GF2p8(0x6c));
        assert_eq!(values[0] + values[2], GF2p8(0xfa));
        assert_eq!(values[0] + values[3], GF2p8(0xb7));
        assert_eq!(values[1] + values[0], GF2p8(0x6c));
        assert_eq!(values[1] + values[1], GF2p8(0x00));
        assert_eq!(values[1] + values[2], GF2p8(0x96));
        assert_eq!(values[1] + values[3], GF2p8(0xdb));
        assert_eq!(values[2] + values[0], GF2p8(0xfa));
        assert_eq!(values[2] + values[1], GF2p8(0x96));
        assert_eq!(values[2] + values[2], GF2p8(0x00));
        assert_eq!(values[2] + values[3], GF2p8(0x4d));
        assert_eq!(values[3] + values[0], GF2p8(0xb7));
        assert_eq!(values[3] + values[1], GF2p8(0xdb));
        assert_eq!(values[3] + values[2], GF2p8(0x4d));
        assert_eq!(values[3] + values[3], GF2p8(0x00));
    }

    #[test]
    fn test_gf2p8_mul() {
        for _ in 0..1000 {
            let val1 = GF2p8::random(thread_rng());
            assert_eq!(val1 * GF2p8::ZERO, GF2p8::ZERO);
            assert_eq!(val1 * GF2p8::ONE, val1);
            let val2 = GF2p8::random(thread_rng());
            assert_eq!(val1 * val2, val2 * val1);
        }

        let values = &GF2P8_VALUES;
        assert_eq!(values[0] * values[0], GF2p8(0xee));
        assert_eq!(values[0] * values[1], GF2p8(0x49));
        assert_eq!(values[0] * values[2], GF2p8(0x8e));
        assert_eq!(values[0] * values[3], GF2p8(0xc1));
        assert_eq!(values[1] * values[0], GF2p8(0x49));
        assert_eq!(values[1] * values[1], GF2p8(0x79));
        assert_eq!(values[1] * values[2], GF2p8(0xeb));
        assert_eq!(values[1] * values[3], GF2p8(0x70));
        assert_eq!(values[2] * values[0], GF2p8(0x8e));
        assert_eq!(values[2] * values[1], GF2p8(0xeb));
        assert_eq!(values[2] * values[2], GF2p8(0xec));
        assert_eq!(values[2] * values[3], GF2p8(0xe9));
        assert_eq!(values[3] * values[0], GF2p8(0xc1));
        assert_eq!(values[3] * values[1], GF2p8(0x70));
        assert_eq!(values[3] * values[2], GF2p8(0xe9));
        assert_eq!(values[3] * values[3], GF2p8(0x16));
    }

    #[test]
    fn test_gf2p8_neg() {
        for _ in 0..1000 {
            let val1 = GF2p8::random(thread_rng());
            assert_eq!(-val1, val1);
        }
    }

    #[test]
    fn test_gf2p8_sub() {
        for _ in 0..1000 {
            let val1 = GF2p8::random(thread_rng());
            let val2 = GF2p8::random(thread_rng());
            assert_eq!(val1 + val2, val1 - val2);
        }
    }

    #[test]
    fn test_gf2p8_inv() {
        assert!(<Choice as Into<bool>>::into(GF2p8::ZERO.invert().is_none()));
        assert_eq!(GF2p8::ONE.invert().unwrap(), GF2p8::ONE);

        for _ in 0..1000 {
            let val1 = GF2p8::random(thread_rng());
            if val1 == GF2p8::ZERO {
                continue;
            }
            assert_eq!(val1 * val1.invert().unwrap(), GF2p8::ONE);
        }

        let values = &GF2P8_VALUES;
        assert_eq!(values[0].invert().unwrap(), GF2p8(0x13));
        assert_eq!(values[1].invert().unwrap(), GF2p8(0xc9));
        assert_eq!(values[2].invert().unwrap(), GF2p8(0xe0));
        assert_eq!(values[3].invert().unwrap(), GF2p8(0xcd));
    }

    #[test]
    fn test_gf2p8_inner_product() {
        for _ in 0..100 {
            let xs: Vec<_> = (0..3).map(|_| GF2p8::random(thread_rng())).collect();
            let ys: Vec<_> = (0..3).map(|_| GF2p8::random(thread_rng())).collect();
            let expected = xs[0] * ys[0] + xs[1] * ys[1] + xs[2] * ys[2];
            let result: GF2p8 = InnerProduct::inner_product(xs.iter(), ys.iter());
            assert_eq!(result, expected);
        }
    }
}
