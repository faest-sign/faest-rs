use crate::arithmetic::clmul_u64;
use crate::field::{BytesRepr, InnerProduct, LazyField, LazyMul, UnreducedField, VecToGF2p128};
use crate::gf2psmall::GF2p8;
use core::arch::x86_64::*;
use core::fmt;
use core::iter::{Product, Sum};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use ff::Field;
use num_traits::identities::{One, Zero};
use rand::{Rng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// pub type GF2p128 = GF2p128Naive;
pub type GF2p128 = GF2p128Fast;

#[derive(Default, Clone, Copy, PartialEq, Eq, bincode::Encode)]
pub struct GF2p128Naive(pub [u64; 2]);

union XmmInitHelper {
    a: __m128i,
    b: u128,
}

#[derive(Clone, Copy)]
pub struct GF2p128Fast(pub __m128i);

impl Default for GF2p128Fast {
    fn default() -> Self {
        Self(unsafe { _mm_setzero_si128() })
    }
}
impl PartialEq for GF2p128Fast {
    fn eq(&self, other: &Self) -> bool {
        unsafe {
            let tmp = _mm_xor_si128(self.0, other.0);
            _mm_test_all_zeros(tmp, tmp) == 1
        }
    }
}
impl Eq for GF2p128Fast {}

impl bincode::Encode for GF2p128Fast {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.to_u128(), encoder)?;
        Ok(())
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct UnreducedGF2p128Naive(pub [u64; 4]);

#[derive(Debug, Clone, Copy)]
pub struct UnreducedGF2p128Fast(pub [__m128i; 2]);

impl Default for UnreducedGF2p128Fast {
    fn default() -> Self {
        Self(unsafe { [_mm_setzero_si128(), _mm_setzero_si128()] })
    }
}
impl PartialEq for UnreducedGF2p128Fast {
    fn eq(&self, other: &Self) -> bool {
        unsafe {
            let lhs = _mm256_set_m128i(self.0[1], self.0[0]);
            let rhs = _mm256_set_m128i(other.0[1], other.0[0]);
            let tmp = _mm256_xor_si256(lhs, rhs);
            _mm256_testz_si256(tmp, tmp) == 0
        }
    }
}
impl Eq for UnreducedGF2p128Fast {}

impl fmt::Debug for GF2p128Fast {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), std::fmt::Error> {
        write!(f, "GF2p128F(0x{:x?})", self.to_u128())
    }
}

impl fmt::Debug for GF2p128Naive {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), std::fmt::Error> {
        write!(f, "GF2p128(0x{:x?})", self.to_u128())
    }
}

impl GF2p128Naive {
    // x^128 + x^7 + x^2 + x + 1
    const POLYNOMIAL: u64 = 0b1000_0111;
    pub const LOG_ORDER: u32 = 128;

    pub const GF2P8_EMBEDDING_POX: [Self; 8] = [
        Self::ONE,
        Self::from_u128(0x053d8555a9979a1ca13fe8ac5560ce0d),
        Self::from_u128(0x4cf4b7439cbfbb84ec7759ca3488aee1),
        Self::from_u128(0x35ad604f7d51d2c6bfcf02ae363946a8),
        Self::from_u128(0x0dcb364640a222fe6b8330483c2e9849),
        Self::from_u128(0x549810e11a88dea5252b49277b1b82b4),
        Self::from_u128(0xd681a5686c0c1f75c72bf2ef2521ff22),
        Self::from_u128(0x0950311a4fb78fe07a7a8e94e136f9bc),
    ];

    pub fn embed_gf2p8(x: GF2p8) -> Self {
        let mut y = Self::ZERO;
        for i in 0..8 {
            if x.0 & (1 << i) != 0 {
                y += Self::GF2P8_EMBEDDING_POX[i];
            }
        }
        y
    }

    fn reduce(mut x: [u64; 4]) -> [u64; 2] {
        debug_assert_eq!(x[3] & (1 << 63), 0);
        for i in (0..63).rev() {
            let t = (x[3] >> i) & 1;
            if i > 46 {
                x[2] ^= t * (Self::POLYNOMIAL >> (64 - i));
            }
            x[1] ^= t * (Self::POLYNOMIAL << i);
        }
        for i in (0..64).rev() {
            let t = (x[2] >> i) & 1;
            if i > 46 {
                x[1] ^= t * (Self::POLYNOMIAL >> (64 - i));
            }
            x[0] ^= t * (Self::POLYNOMIAL << i);
        }
        [x[0], x[1]]
    }

    pub const fn from_u128(val: u128) -> Self {
        Self([val as u64, (val >> 64) as u64])
    }

    pub const fn to_u128(self) -> u128 {
        (self.0[1] as u128) << 64 | self.0[0] as u128
    }
}

impl GF2p128Fast {
    // x^128 + x^17 + x^9 + x^7 + 1)
    // const POLYNOMIAL: u64 = 0b0010_0000_0010_1000_0001;
    // alternative: x^128 + x^7 + x^2 + x + 1
    //              const POLYNOMIAL: u64 = 0b1000_0111;
    pub const LOG_ORDER: u32 = 128;

    pub const GF2P8_EMBEDDING_POX: [Self; 8] = [
        Self::ONE,
        Self::from_u128(0x053d8555a9979a1ca13fe8ac5560ce0d),
        Self::from_u128(0x4cf4b7439cbfbb84ec7759ca3488aee1),
        Self::from_u128(0x35ad604f7d51d2c6bfcf02ae363946a8),
        Self::from_u128(0x0dcb364640a222fe6b8330483c2e9849),
        Self::from_u128(0x549810e11a88dea5252b49277b1b82b4),
        Self::from_u128(0xd681a5686c0c1f75c72bf2ef2521ff22),
        Self::from_u128(0x0950311a4fb78fe07a7a8e94e136f9bc),
    ];

    pub fn embed_gf2p8(x: GF2p8) -> Self {
        let mut y = Self::ZERO;
        for i in 0..8 {
            if x.0 & (1 << i) != 0 {
                y += Self::GF2P8_EMBEDDING_POX[i];
            }
        }
        y
    }

    #[inline(always)]
    fn gfmul(a: __m128i, b: __m128i) -> [__m128i; 2] {
        unsafe {
            let c = _mm_clmulepi64_si128(a, b, 0x11);
            let d = _mm_clmulepi64_si128(a, b, 0x00);
            let a_with_swapped_words = _mm_shuffle_epi32(a, 0b_01_00_11_10);
            let b_with_swapped_words = _mm_shuffle_epi32(b, 0b_01_00_11_10);
            let e = _mm_clmulepi64_si128(
                _mm_xor_si128(a, a_with_swapped_words),
                _mm_xor_si128(b, b_with_swapped_words),
                0x00,
            );
            let tmp = _mm_xor_si128(e, _mm_xor_si128(c, d));
            let res_lo = _mm_xor_si128(d, _mm_slli_si128(tmp, 8));
            let res_hi = _mm_xor_si128(c, _mm_srli_si128(tmp, 8));
            [res_lo, res_hi]
        }
    }

    #[inline(always)]
    fn reduce(x: [__m128i; 2]) -> __m128i {
        unsafe {
            let [lo, hi] = x;
            let xmmmask = _mm_setr_epi32(i32::MAX, 0x0, 0x0, 0x0);
            let tmp7 = _mm_srli_epi32(hi, 31);
            let tmp8 = _mm_srli_epi32(hi, 30);
            let tmp9 = _mm_srli_epi32(hi, 25);
            let tmp7 = _mm_xor_si128(tmp7, _mm_xor_si128(tmp8, tmp9));
            let tmp8 = _mm_shuffle_epi32(tmp7, 0b_10_01_00_11);
            let tmp7 = _mm_and_si128(xmmmask, tmp8);
            let tmp8 = _mm_andnot_si128(xmmmask, tmp8);
            let tmp3 = _mm_xor_si128(lo, tmp8);
            let tmp6 = _mm_xor_si128(hi, tmp7);
            let tmp10 = _mm_slli_epi32(tmp6, 1);
            let tmp3 = _mm_xor_si128(tmp3, tmp10);
            let tmp11 = _mm_slli_epi32(tmp6, 2);
            let tmp3 = _mm_xor_si128(tmp3, tmp11);
            let tmp12 = _mm_slli_epi32(tmp6, 7);
            let tmp3 = _mm_xor_si128(tmp3, tmp12);
            _mm_xor_si128(tmp3, tmp6)
        }
    }

    pub const fn from_u128(val: u128) -> Self {
        Self(unsafe { XmmInitHelper { b: val }.a })
    }

    pub const fn to_u128(self) -> u128 {
        unsafe { XmmInitHelper { a: self.0 }.b }
    }
}

impl From<u128> for GF2p128Naive {
    fn from(val: u128) -> Self {
        Self::from_u128(val)
    }
}

impl BytesRepr for GF2p128Naive {
    type Repr = [u8; 16];
    #[inline(always)]
    fn to_repr(self) -> Self::Repr {
        let mut z = [0u8; 16];
        z[..8].copy_from_slice(&self.0[0].to_le_bytes());
        z[8..].copy_from_slice(&self.0[1].to_le_bytes());
        z
    }
    #[inline(always)]
    fn from_repr(repr: Self::Repr) -> Self {
        Self([
            u64::from_le_bytes(repr[..8].try_into().unwrap()),
            u64::from_le_bytes(repr[8..].try_into().unwrap()),
        ])
    }
}

impl BytesRepr for GF2p128Fast {
    type Repr = [u8; 16];
    #[inline(always)]
    fn to_repr(self) -> Self::Repr {
        let mut z = [0u8; 16];
        unsafe {
            _mm_storeu_si128(z.as_mut_ptr() as *mut __m128i, self.0);
        }
        z
    }
    #[inline(always)]
    fn from_repr(repr: Self::Repr) -> Self {
        Self(unsafe { _mm_loadu_si128(repr.as_ptr() as *const __m128i) })
    }
}

impl LazyField for GF2p128Naive {
    type UF = UnreducedGF2p128Naive;
}
impl LazyField for GF2p128Fast {
    type UF = UnreducedGF2p128Fast;
}

impl UnreducedField for UnreducedGF2p128Naive {
    type ReducedField = GF2p128Naive;
    const ZERO: Self = Self([0, 0, 0, 0]);
    #[inline(always)]
    fn reduce(self) -> Self::ReducedField {
        GF2p128Naive(GF2p128Naive::reduce(self.0))
    }
}
impl UnreducedField for UnreducedGF2p128Fast {
    type ReducedField = GF2p128Fast;
    const ZERO: Self = Self(unsafe { [XmmInitHelper { b: 0 }.a, XmmInitHelper { b: 0 }.a] });
    #[inline(always)]
    fn reduce(self) -> Self::ReducedField {
        GF2p128Fast(GF2p128Fast::reduce(self.0))
    }
}

impl Add for GF2p128Naive {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: Self) -> Self::Output {
        Self([self.0[0] ^ other.0[0], self.0[1] ^ other.0[1]])
    }
}
impl Add for UnreducedGF2p128Naive {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: Self) -> Self::Output {
        Self([
            self.0[0] ^ other.0[0],
            self.0[1] ^ other.0[1],
            self.0[2] ^ other.0[2],
            self.0[3] ^ other.0[3],
        ])
    }
}

impl Add for GF2p128Fast {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: Self) -> Self::Output {
        Self(unsafe { _mm_xor_si128(self.0, other.0) })
    }
}
impl Add for UnreducedGF2p128Fast {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: Self) -> Self::Output {
        Self(unsafe {
            [
                _mm_xor_si128(self.0[0], other.0[0]),
                _mm_xor_si128(self.0[1], other.0[1]),
            ]
        })
    }
}
impl Add<<GF2p128Naive as LazyField>::UF> for GF2p128Naive {
    type Output = <GF2p128Naive as LazyField>::UF;
    #[inline(always)]
    fn add(self, other: <GF2p128Naive as LazyField>::UF) -> Self::Output {
        let mut z = other.0;
        z[0] ^= self.0[0];
        z[1] ^= self.0[1];
        UnreducedGF2p128Naive(z)
    }
}
impl Add<<GF2p128Fast as LazyField>::UF> for GF2p128Fast {
    type Output = <GF2p128Fast as LazyField>::UF;
    #[inline(always)]
    fn add(self, other: <GF2p128Fast as LazyField>::UF) -> Self::Output {
        UnreducedGF2p128Fast(unsafe { [_mm_xor_si128(self.0, other.0[0]), other.0[1]] })
    }
}

impl LazyMul for GF2p128Naive {
    type Output = UnreducedGF2p128Naive;
    #[inline(always)]
    fn lazy_mul(self, other: Self) -> Self::Output {
        let x = self.0;
        let y = other.0;
        let xy_00 = clmul_u64(x[0], y[0]);
        let xy_10 = clmul_u64(x[1], y[0]);
        let xy_01 = clmul_u64(x[0], y[1]);
        let xy_11 = clmul_u64(x[1], y[1]);
        let mut z = [0u64; 4];
        z[0] = xy_00 as u64;
        z[1] = (xy_00 >> 64) as u64 ^ xy_10 as u64 ^ xy_01 as u64;
        z[2] = (xy_10 >> 64) as u64 ^ (xy_01 >> 64) as u64 ^ xy_11 as u64;
        z[3] = (xy_11 >> 64) as u64;
        UnreducedGF2p128Naive(z)
    }
}
impl LazyMul for GF2p128Fast {
    type Output = UnreducedGF2p128Fast;
    #[inline(always)]
    fn lazy_mul(self, other: Self) -> Self::Output {
        UnreducedGF2p128Fast(Self::gfmul(self.0, other.0))
    }
}

impl Mul for GF2p128Naive {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: Self) -> Self::Output {
        let x = self.0;
        let y = other.0;
        let xy_00 = clmul_u64(x[0], y[0]);
        let xy_10 = clmul_u64(x[1], y[0]);
        let xy_01 = clmul_u64(x[0], y[1]);
        let xy_11 = clmul_u64(x[1], y[1]);
        let mut z = [0u64; 4];
        z[0] = xy_00 as u64;
        z[1] = (xy_00 >> 64) as u64 ^ xy_10 as u64 ^ xy_01 as u64;
        z[2] = (xy_10 >> 64) as u64 ^ (xy_01 >> 64) as u64 ^ xy_11 as u64;
        z[3] = (xy_11 >> 64) as u64;
        Self(Self::reduce(z))
    }
}
impl Mul for GF2p128Fast {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: Self) -> Self::Output {
        Self(Self::reduce(Self::gfmul(self.0, other.0)))
    }
}
impl ConstantTimeEq for GF2p128Naive {
    #[inline(always)]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}
impl ConditionallySelectable for GF2p128Naive {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
        ])
    }
}
impl ConstantTimeEq for GF2p128Fast {
    #[inline(always)]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.to_u128().ct_eq(&other.to_u128())
    }
}
impl ConditionallySelectable for GF2p128Fast {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self::from_u128(u128::conditional_select(&a.to_u128(), &b.to_u128(), choice))
    }
}

impl Field for GF2p128Naive {
    const ZERO: Self = Self([0, 0]);
    const ONE: Self = Self([1, 0]);

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
        for _ in 0..126 {
            y = x * y;
            x = x.square();
        }
        CtOption::new(x * y, !self.ct_eq(&Self::ZERO))
    }

    fn sqrt_ratio(_num: &Self, _div: &Self) -> (Choice, Self) {
        panic!("not implemented")
    }
}

impl Field for GF2p128Fast {
    const ZERO: Self = Self::from_u128(0);
    const ONE: Self = Self::from_u128(1);

    #[inline(always)]
    fn random(mut rng: impl RngCore) -> Self {
        let bytes: [u8; 16] = rng.gen();
        Self(unsafe { _mm_loadu_si128(bytes.as_ptr() as *const __m128i) })
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
        for _ in 0..126 {
            y = x * y;
            x = x.square();
        }
        CtOption::new(x * y, !self.ct_eq(&Self::ZERO))
    }

    fn sqrt_ratio(_num: &Self, _div: &Self) -> (Choice, Self) {
        panic!("not implemented")
    }
}

impl VecToGF2p128<GF2p128Naive> for GF2p8 {
    const VECTOR_SIZE: usize = 16;
    fn convert(vec: &[Self]) -> GF2p128Naive {
        debug_assert_eq!(vec.len(), <Self as VecToGF2p128<GF2p128Naive>>::VECTOR_SIZE);
        let mut bytes = [0u8; <Self as VecToGF2p128<GF2p128Naive>>::VECTOR_SIZE];
        for (b, x) in bytes.iter_mut().zip(vec.iter()) {
            *b = x.0;
        }
        GF2p128Naive::from_repr(bytes)
    }
}
impl VecToGF2p128<GF2p128Fast> for GF2p8 {
    const VECTOR_SIZE: usize = 16;
    fn convert(vec: &[Self]) -> GF2p128Fast {
        debug_assert_eq!(vec.len(), <Self as VecToGF2p128<GF2p128Fast>>::VECTOR_SIZE);
        let mut bytes = [0u8; <Self as VecToGF2p128<GF2p128Fast>>::VECTOR_SIZE];
        for (b, x) in bytes.iter_mut().zip(vec.iter()) {
            *b = x.0;
        }
        GF2p128Fast::from_repr(bytes)
    }
}

impl_additional_field_arithmetic!(GF2p128Naive, UnreducedGF2p128Naive);
impl_additional_field_arithmetic!(GF2p128Fast, UnreducedGF2p128Fast);

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    const GF2P128_VALUES: [u128; 4] = [
        0x616ce329b8aee6b6c752890eaec5fdca,
        0xc9d404e824cf265e31f38f85087a788f,
        0x7eceb30086f94c398a6881ea95c275b2,
        0x41a376dbe3eebf96a8f49dff52f12e1c,
    ];
    const GF2P128_SUMS: [u128; 16] = [
        0x00000000000000000000000000000000,
        0xa8b8e7c19c61c0e8f6a1068ba6bf8545,
        0x1fa250293e57aa8f4d3a08e43b078878,
        0x20cf95f25b4059206fa614f1fc34d3d6,
        0xa8b8e7c19c61c0e8f6a1068ba6bf8545,
        0x00000000000000000000000000000000,
        0xb71ab7e8a2366a67bb9b0e6f9db80d3d,
        0x88777233c72199c89907127a5a8b5693,
        0x1fa250293e57aa8f4d3a08e43b078878,
        0xb71ab7e8a2366a67bb9b0e6f9db80d3d,
        0x00000000000000000000000000000000,
        0x3f6dc5db6517f3af229c1c15c7335bae,
        0x20cf95f25b4059206fa614f1fc34d3d6,
        0x88777233c72199c89907127a5a8b5693,
        0x3f6dc5db6517f3af229c1c15c7335bae,
        0x00000000000000000000000000000000,
    ];
    const GF2P128_PRODS: [u128; 16] = [
        0x3c98549feed83d303eb7a796f31e041e,
        0x1a5c768b84c413168688f04271930c5b,
        0x5c05635ff5c2e1a60d827ab85b1cfb41,
        0x73e65c51b95ada0467c8bc6f148c64fa,
        0x1a5c768b84c413168688f04271930c5b,
        0x946f6a75480ecdd314198e6d032488a6,
        0x9b4072097abf4b1a997e59768adb412e,
        0xe89a5742eb99f471742434defe7a4910,
        0x5c05635ff5c2e1a60d827ab85b1cfb41,
        0x9b4072097abf4b1a997e59768adb412e,
        0x80c18fcf199a54658b565a4b4ca3fa75,
        0xd47377c7c8349377b23204a9f433768e,
        0x73e65c51b95ada0467c8bc6f148c64fa,
        0xe89a5742eb99f471742434defe7a4910,
        0xd47377c7c8349377b23204a9f433768e,
        0x34e48b81a0144125bfb4d28e745e4804,
    ];
    const GF2P128_INVS: [u128; 4] = [
        0xf5d6d8bbcc5ea1a34385c82096abca34,
        0xd2153fdf6c54531d35c78c543208cd55,
        0x6261060d33482531fab6c55e5dfe7155,
        0xb27283763acb647fd8600be2a72b7d66,
    ];

    const GF2P8_VALUES: [GF2p8; 4] = [GF2p8(0x4b), GF2p8(0x27), GF2p8(0xb1), GF2p8(0xfc)];
    const GF2P8_IN_GF2P128_EMBEDDINGS: [u128; 4] = [
        0xe6114072b8ca57afd9db18ed46787786,
        0x1d5122f72fa0ff3d6863f8411af3e259,
        0x500317bd159d73bb34d2f7fba603e340,
        0xffdb65d9987f058ca0415e708193f42a,
    ];

    macro_rules! make_tests_gf2p128 {
        (
            $field_name:ident,
            $test_add_name:ident,
            $test_mul_name:ident,
            $test_neg_name:ident,
            $test_sub_name:ident,
            $test_inv_name:ident,
            $test_gf2p8_embedding_name:ident,
            $test_inner_product_name:ident,
            $test_to_from_u128_name:ident
        ) => {
            #[test]
            fn $test_add_name() {
                type F = $field_name;
                for _ in 0..1000 {
                    let val1 = F::random(thread_rng());
                    assert_eq!(val1 + F::ZERO, val1);
                    assert_eq!(val1 + val1, F::ZERO);
                    let val2 = F::random(thread_rng());
                    assert_eq!(val1 + val2, val2 + val1);
                }

                let values = GF2P128_VALUES.map(F::from_u128);
                let sums = GF2P128_SUMS.map(F::from_u128);
                assert_eq!(values[0] + values[0], sums[0]);
                assert_eq!(values[0] + values[1], sums[1]);
                assert_eq!(values[0] + values[2], sums[2]);
                assert_eq!(values[0] + values[3], sums[3]);
                assert_eq!(values[1] + values[0], sums[4]);
                assert_eq!(values[1] + values[1], sums[5]);
                assert_eq!(values[1] + values[2], sums[6]);
                assert_eq!(values[1] + values[3], sums[7]);
                assert_eq!(values[2] + values[0], sums[8]);
                assert_eq!(values[2] + values[1], sums[9]);
                assert_eq!(values[2] + values[2], sums[10]);
                assert_eq!(values[2] + values[3], sums[11]);
                assert_eq!(values[3] + values[0], sums[12]);
                assert_eq!(values[3] + values[1], sums[13]);
                assert_eq!(values[3] + values[2], sums[14]);
                assert_eq!(values[3] + values[3], sums[15]);
            }

            #[test]
            fn $test_mul_name() {
                type F = $field_name;
                for _ in 0..1000 {
                    let val1 = F::random(thread_rng());
                    assert_eq!(
                        val1 * F::ZERO,
                        F::ZERO,
                        "multiplication with zero is not zero"
                    );
                    assert_eq!(val1 * F::ONE, val1, "one is not the identity");
                    let val2 = F::random(thread_rng());
                    assert_eq!(
                        val1 * val2,
                        val2 * val1,
                        "multiplication is not commutative"
                    );
                }

                let values = GF2P128_VALUES.map(F::from_u128);
                let prods = GF2P128_PRODS.map(F::from_u128);
                assert_eq!(values[0] * values[0], prods[0]);
                assert_eq!(values[0] * values[1], prods[1]);
                assert_eq!(values[0] * values[2], prods[2]);
                assert_eq!(values[0] * values[3], prods[3]);
                assert_eq!(values[1] * values[0], prods[4]);
                assert_eq!(values[1] * values[1], prods[5]);
                assert_eq!(values[1] * values[2], prods[6]);
                assert_eq!(values[1] * values[3], prods[7]);
                assert_eq!(values[2] * values[0], prods[8]);
                assert_eq!(values[2] * values[1], prods[9]);
                assert_eq!(values[2] * values[2], prods[10]);
                assert_eq!(values[2] * values[3], prods[11]);
                assert_eq!(values[3] * values[0], prods[12]);
                assert_eq!(values[3] * values[1], prods[13]);
                assert_eq!(values[3] * values[2], prods[14]);
                assert_eq!(values[3] * values[3], prods[15]);
            }
            #[test]
            fn $test_neg_name() {
                type F = $field_name;
                for _ in 0..1000 {
                    let val1 = F::random(thread_rng());
                    assert_eq!(-val1, val1);
                }
            }

            #[test]
            fn $test_sub_name() {
                type F = $field_name;
                for _ in 0..1000 {
                    let val1 = F::random(thread_rng());
                    let val2 = F::random(thread_rng());
                    assert_eq!(val1 + val2, val1 - val2);
                }
            }

            #[test]
            fn $test_inv_name() {
                type F = $field_name;
                assert!(<Choice as Into<bool>>::into(F::ZERO.invert().is_none()));
                assert_eq!(F::ONE.invert().unwrap(), F::ONE);

                for _ in 0..100 {
                    let val1 = F::random(thread_rng());
                    if val1 == F::ZERO {
                        continue;
                    }
                    assert_eq!(val1 * val1.invert().unwrap(), F::ONE);
                }

                let values = GF2P128_VALUES.map(F::from_u128);
                let invs = GF2P128_INVS.map(F::from_u128);
                assert_eq!(values[0].invert().unwrap(), invs[0]);
                assert_eq!(values[1].invert().unwrap(), invs[1]);
                assert_eq!(values[2].invert().unwrap(), invs[2]);
                assert_eq!(values[3].invert().unwrap(), invs[3]);
            }

            #[test]
            fn $test_gf2p8_embedding_name() {
                type F = $field_name;
                for _ in 0..100 {
                    let val1 = GF2p8::random(thread_rng());
                    let val2 = GF2p8::random(thread_rng());
                    assert_eq!(
                        F::embed_gf2p8(val1) + F::embed_gf2p8(val2),
                        F::embed_gf2p8(val1 + val2)
                    );
                    assert_eq!(
                        F::embed_gf2p8(val1) * F::embed_gf2p8(val2),
                        F::embed_gf2p8(val1 * val2)
                    );
                    if val1 != GF2p8::ZERO {
                        assert_eq!(
                            F::embed_gf2p8(val1).invert().unwrap(),
                            F::embed_gf2p8(val1.invert().unwrap())
                        );
                    }
                }

                let values = &GF2P8_VALUES;
                let embedded_values = GF2P8_IN_GF2P128_EMBEDDINGS.map(F::from_u128);
                assert_eq!(F::embed_gf2p8(GF2p8::ZERO), F::ZERO);
                assert_eq!(F::embed_gf2p8(GF2p8::ONE), F::ONE);
                assert_eq!(F::embed_gf2p8(values[0]), embedded_values[0]);
                assert_eq!(F::embed_gf2p8(values[1]), embedded_values[1]);
                assert_eq!(F::embed_gf2p8(values[2]), embedded_values[2]);
                assert_eq!(F::embed_gf2p8(values[3]), embedded_values[3]);
            }

            #[test]
            fn $test_inner_product_name() {
                type F = $field_name;
                for _ in 0..100 {
                    let xs: Vec<_> = (0..3).map(|_| F::random(thread_rng())).collect();
                    let ys: Vec<_> = (0..3).map(|_| F::random(thread_rng())).collect();
                    let expected = xs[0] * ys[0] + xs[1] * ys[1] + xs[2] * ys[2];
                    let result: F = InnerProduct::inner_product(xs.iter(), ys.iter());
                    assert_eq!(result, expected);
                }
            }

            #[test]
            fn $test_to_from_u128_name() {
                type F = $field_name;
                for _ in 0..100 {
                    let x = thread_rng().gen();
                    let y = F::from_u128(x).to_u128();
                    assert_eq!(x, y);
                }
            }
        };
    }

    make_tests_gf2p128!(
        GF2p128Naive,
        test_gf2p128_naive_add,
        test_gf2p128_naive_mul,
        test_gf2p128_naive_neg,
        test_gf2p128_naive_sub,
        test_gf2p128_naive_inv,
        test_gf2p128_naive_gf2p8_embedding,
        test_gf2p128_naive_inner_product,
        test_gf2p128_naive_to_from_u128
    );
    make_tests_gf2p128!(
        GF2p128Fast,
        test_gf2p128_fast_add,
        test_gf2p128_fast_mul,
        test_gf2p128_fast_neg,
        test_gf2p128_fast_sub,
        test_gf2p128_fast_inv,
        test_gf2p128_fast_gf2p8_embedding,
        test_gf2p128_fast_inner_product,
        test_gf2p128_fast_to_from_u128
    );
}
