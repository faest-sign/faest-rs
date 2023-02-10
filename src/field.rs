use core::fmt;
use core::iter::{Product, Sum};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use ff::Field;
use num_traits::identities::{One, Zero};
use rand::{Rng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

macro_rules! make_clmul {
    ($function_name:ident, $small_uint:ty, $large_uint:ty) => {
        fn $function_name(x: $small_uint, y: $small_uint) -> $large_uint {
            let x = x as $large_uint;
            let y = y as $large_uint;
            let mut result = 0;
            for i in 0..<$small_uint>::BITS {
                result ^= ((x >> i) & 1) * (y << i);
            }
            result
        }
    };
}

make_clmul!(clmul_u8, u8, u16);
make_clmul!(clmul_u64, u64, u128);

macro_rules! impl_arithmetic {
    ($field_name:ident, $unreduced_field_name:ident) => {
        impl Add<&Self> for $field_name {
            type Output = Self;
            #[inline(always)]
            fn add(self, other: &Self) -> Self::Output {
                self + *other
            }
        }
        impl AddAssign for $field_name {
            #[inline(always)]
            fn add_assign(&mut self, other: Self) {
                *self = *self + other;
            }
        }
        impl AddAssign<&Self> for $field_name {
            #[inline(always)]
            fn add_assign(&mut self, other: &Self) {
                *self = *self + other;
            }
        }
        impl Div for $field_name {
            type Output = Self;
            #[inline(always)]
            fn div(self, other: Self) -> Self::Output {
                if other == Self::ZERO {
                    panic!("division by zero");
                }
                self * other.invert().unwrap()
            }
        }
        impl Div<&Self> for $field_name {
            type Output = Self;
            #[inline(always)]
            fn div(self, other: &Self) -> Self::Output {
                self / *other
            }
        }
        impl DivAssign for $field_name {
            #[inline(always)]
            fn div_assign(&mut self, other: Self) {
                *self = *self / other;
            }
        }
        impl DivAssign<&Self> for $field_name {
            #[inline(always)]
            fn div_assign(&mut self, other: &Self) {
                *self = *self / other;
            }
        }
        impl Mul<&Self> for $field_name {
            type Output = Self;
            #[inline(always)]
            fn mul(self, other: &Self) -> Self::Output {
                self * *other
            }
        }
        impl MulAssign for $field_name {
            #[inline(always)]
            fn mul_assign(&mut self, other: Self) {
                *self = *self * other;
            }
        }
        impl MulAssign<&Self> for $field_name {
            #[inline(always)]
            fn mul_assign(&mut self, other: &Self) {
                *self = *self * other;
            }
        }
        impl Neg for $field_name {
            type Output = Self;
            #[inline(always)]
            fn neg(self) -> Self::Output {
                self
            }
        }
        #[allow(clippy::suspicious_arithmetic_impl)]
        impl Sub for $field_name {
            type Output = Self;
            #[inline(always)]
            fn sub(self, other: Self) -> Self::Output {
                self + other
            }
        }
        #[allow(clippy::suspicious_arithmetic_impl)]
        impl Sub<&Self> for $field_name {
            type Output = Self;
            #[inline(always)]
            fn sub(self, other: &Self) -> Self::Output {
                self + other
            }
        }
        #[allow(clippy::suspicious_op_assign_impl)]
        impl SubAssign for $field_name {
            #[inline(always)]
            fn sub_assign(&mut self, other: Self) {
                *self += other
            }
        }
        #[allow(clippy::suspicious_op_assign_impl)]
        impl SubAssign<&Self> for $field_name {
            #[inline(always)]
            fn sub_assign(&mut self, other: &Self) {
                *self += other
            }
        }
        #[allow(clippy::suspicious_arithmetic_impl)]
        impl Sub<<$field_name as LazyField>::UF> for $field_name {
            type Output = <$field_name as LazyField>::UF;
            #[inline(always)]
            fn sub(self, other: <$field_name as LazyField>::UF) -> Self::Output {
                self + other
            }
        }
        impl Sum for $field_name {
            #[inline(always)]
            fn sum<I>(iter: I) -> Self
            where
                I: Iterator<Item = Self>,
            {
                iter.fold(Self::ZERO, Add::add)
            }
        }
        impl<'a> Sum<&'a Self> for $field_name {
            #[inline(always)]
            fn sum<I>(iter: I) -> Self
            where
                I: Iterator<Item = &'a Self>,
            {
                iter.copied().sum()
            }
        }
        impl Sum for $unreduced_field_name {
            #[inline(always)]
            fn sum<I>(iter: I) -> Self
            where
                I: Iterator<Item = Self>,
            {
                iter.fold(Self::ZERO, Add::add)
            }
        }
        impl<'a> Sum<&'a Self> for $unreduced_field_name {
            #[inline(always)]
            fn sum<I>(iter: I) -> Self
            where
                I: Iterator<Item = &'a Self>,
            {
                iter.copied().sum()
            }
        }

        impl Product for $field_name {
            #[inline(always)]
            fn product<I>(iter: I) -> Self
            where
                I: Iterator<Item = Self>,
            {
                iter.fold(Self::ONE, Mul::mul)
            }
        }
        impl<'a> Product<&'a Self> for $field_name {
            #[inline(always)]
            fn product<I>(iter: I) -> Self
            where
                I: Iterator<Item = &'a Self>,
            {
                iter.copied().product()
            }
        }
        impl InnerProduct for $field_name {
            #[inline(always)]
            fn inner_product<I, J>(iter1: I, iter2: J) -> Self
            where
                I: Iterator<Item = Self>,
                J: Iterator<Item = Self>,
            {
                iter1
                    .zip(iter2)
                    .map(|(x, y)| LazyMul::lazy_mul(x, y))
                    .sum::<<$field_name as LazyMul>::Output>()
                    .reduce()
            }
        }
        impl<'a> InnerProduct<&'a Self> for $field_name {
            #[inline(always)]
            fn inner_product<I, J>(iter1: I, iter2: J) -> Self
            where
                I: Iterator<Item = &'a Self>,
                J: Iterator<Item = &'a Self>,
            {
                InnerProduct::inner_product(iter1.copied(), iter2.copied())
            }
        }
        impl ConstantTimeEq for $field_name {
            #[inline(always)]
            fn ct_eq(&self, other: &Self) -> Choice {
                self.0.ct_eq(&other.0)
            }
        }
        impl Zero for $field_name {
            fn zero() -> Self {
                Self::ZERO
            }
            fn is_zero(&self) -> bool {
                *self == Self::ZERO
            }
        }

        impl One for $field_name {
            fn one() -> Self {
                Self::ONE
            }
            fn is_one(&self) -> bool {
                *self == Self::ONE
            }
        }
    };
}

pub trait UnreducedField {
    type ReducedField: Field;
    const ZERO: Self;
    fn reduce(self) -> Self::ReducedField;
}

pub trait LazyMul {
    type Output: UnreducedField;
    fn lazy_mul(self, other: Self) -> Self::Output;
}

pub trait LazyField: Add<Self::UF> + Sub<Self::UF> + LazyMul<Output = Self::UF> {
    type UF: UnreducedField;
}

pub trait InnerProduct<A = Self> {
    fn inner_product<I, J>(iter1: I, iter2: J) -> Self
    where
        I: Iterator<Item = A>,
        J: Iterator<Item = A>;
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct GF2p8(pub u8);

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct UnreducedGF2p8(pub u16);

#[derive(Default, Clone, Copy, PartialEq, Eq)]
pub struct GF2p128(pub [u64; 2]);

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct UnreducedGF2p128(pub [u64; 4]);

impl fmt::Debug for GF2p128 {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), std::fmt::Error> {
        write!(f, "GF2p128(0x{:x?})", self.to_u128())
    }
}

pub trait BytesRepr {
    type Repr: AsRef<[u8]> + AsMut<[u8]>;
    fn to_repr(self) -> Self::Repr;
    fn from_repr(repr: Self::Repr) -> Self;
}

impl GF2p8 {
    // x^8 + x^4 + x^3 + x + 1
    const POLYNOMIAL: u16 = 0x11b;
    pub const ORDER: usize = 256;
    pub const LOG_ORDER: u32 = 8;

    fn reduce(mut x: u16) -> u8 {
        for i in (u8::BITS..(u16::BITS - 1)).rev() {
            x ^= ((x >> i) & 1) * (Self::POLYNOMIAL << (i - u8::BITS));
        }
        x as u8
    }

    pub fn rotate_left(self, n: u32) -> Self {
        Self(self.0.rotate_left(n))
    }
}

impl GF2p128 {
    // x^128 + x^17 + x^9 + x^7 + 1)
    const POLYNOMIAL: u64 = 0b0010_0000_0010_1000_0001;
    // alternative: x^128 + x^7 + x^2 + x + 1
    //              const POLYNOMIAL: u64 = 0b1000_0111;
    pub const LOG_ORDER: u32 = 128;

    pub const GF2P8_EMBEDDING_POX: [Self; 8] = [
        Self::ONE,
        Self::from_u128(0x9ccca334694bcdda1d9c309153a6e2eb),
        Self::from_u128(0xe82911a043adc4b187665cfa2996a803),
        Self::from_u128(0x8e7992c810cd11f7dcab7d771ad26064),
        Self::from_u128(0xbcddb0d51a7b316694cc2f6ce7fa6a86),
        Self::from_u128(0x23c67f98ba9ec974c7f80e43d7ba4cae),
        Self::from_u128(0xb9a6f12e7e5de757206e89d43c7b758d),
        Self::from_u128(0xda3096a65c7ae8953643e598c84ff4ce),
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

impl From<u128> for GF2p128 {
    fn from(val: u128) -> Self {
        Self::from_u128(val)
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

impl BytesRepr for GF2p128 {
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

impl LazyField for GF2p8 {
    type UF = UnreducedGF2p8;
}
impl LazyField for GF2p128 {
    type UF = UnreducedGF2p128;
}

impl UnreducedField for UnreducedGF2p8 {
    type ReducedField = GF2p8;
    const ZERO: Self = Self(0);
    #[inline(always)]
    fn reduce(self) -> Self::ReducedField {
        GF2p8(GF2p8::reduce(self.0))
    }
}
impl UnreducedField for UnreducedGF2p128 {
    type ReducedField = GF2p128;
    const ZERO: Self = Self([0, 0, 0, 0]);
    #[inline(always)]
    fn reduce(self) -> Self::ReducedField {
        GF2p128(GF2p128::reduce(self.0))
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

impl Add for GF2p128 {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: Self) -> Self::Output {
        Self([self.0[0] ^ other.0[0], self.0[1] ^ other.0[1]])
    }
}
impl Add for UnreducedGF2p128 {
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

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add<<GF2p8 as LazyField>::UF> for GF2p8 {
    type Output = <GF2p8 as LazyField>::UF;
    #[inline(always)]
    fn add(self, other: <GF2p8 as LazyField>::UF) -> Self::Output {
        UnreducedGF2p8(self.0 as u16 ^ other.0)
    }
}

// #[allow(clippy::suspicious_arithmetic_impl)]
impl Add<<GF2p128 as LazyField>::UF> for GF2p128 {
    type Output = <GF2p128 as LazyField>::UF;
    #[inline(always)]
    fn add(self, other: <GF2p128 as LazyField>::UF) -> Self::Output {
        let mut z = other.0;
        z[0] ^= self.0[0];
        z[1] ^= self.0[1];
        UnreducedGF2p128(z)
    }
}

impl LazyMul for GF2p8 {
    type Output = UnreducedGF2p8;
    #[inline(always)]
    fn lazy_mul(self, other: Self) -> Self::Output {
        UnreducedGF2p8(clmul_u8(self.0, other.0))
    }
}
impl LazyMul for GF2p128 {
    type Output = UnreducedGF2p128;
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
        UnreducedGF2p128(z)
    }
}

impl Mul for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: Self) -> Self::Output {
        Self(Self::reduce(clmul_u8(self.0, other.0)))
    }
}
impl Mul for GF2p128 {
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

impl ConditionallySelectable for GF2p8 {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self(u8::conditional_select(&a.0, &b.0, choice))
    }
}
impl ConditionallySelectable for GF2p128 {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
        ])
    }
}

impl_arithmetic!(GF2p8, UnreducedGF2p8);
impl_arithmetic!(GF2p128, UnreducedGF2p128);

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

impl Field for GF2p128 {
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
        for i in 0..126 {
            y = x * y;
            x = x.square();
            eprintln!("{i:3}: y = {y:?}");
            eprintln!("{i:3}: x = {x:?}");
        }
        CtOption::new(x * y, !self.ct_eq(&Self::ZERO))
    }

    fn sqrt_ratio(_num: &Self, _div: &Self) -> (Choice, Self) {
        panic!("not implemented")
    }
}

pub trait VecToGF2p128: Sized {
    const VECTOR_SIZE: usize;
    fn convert(vec: &[Self]) -> GF2p128;
}

impl VecToGF2p128 for GF2p8 {
    const VECTOR_SIZE: usize = 16;
    fn convert(vec: &[Self]) -> GF2p128 {
        debug_assert_eq!(vec.len(), Self::VECTOR_SIZE);
        let mut bytes = [0u8; Self::VECTOR_SIZE];
        for (b, x) in bytes.iter_mut().zip(vec.iter()) {
            *b = x.0;
        }
        GF2p128::from_repr(bytes)
    }
}

impl From<u8> for GF2p8 {
    fn from(val: u8) -> Self {
        Self(val)
    }
}

impl Into<u8> for GF2p8 {
    fn into(self) -> u8 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn test_clmul_u8() {
        let a = 0xe6;
        let b = 0xcd;
        let c = 0x4ece;
        assert_eq!(clmul_u8(a, b), c);
    }

    #[test]
    fn test_clmul_u64() {
        let a = 0x8cedd8695a31a3c0;
        let b = 0x89eb6fd18fbf3991;
        let c = 0x42ec785c6a2a23f8276b5a8ce243bfc0;
        assert_eq!(clmul_u64(a, b), c);
    }

    const GF2P8_VALUES: [GF2p8; 4] = [GF2p8(0x4b), GF2p8(0x27), GF2p8(0xb1), GF2p8(0xfc)];
    const GF2P128_VALUES: [GF2p128; 4] = [
        GF2p128::from_u128(0x616ce329b8aee6b6c752890eaec5fdca),
        GF2p128::from_u128(0xc9d404e824cf265e31f38f85087a788f),
        GF2p128::from_u128(0x7eceb30086f94c398a6881ea95c275b2),
        GF2p128::from_u128(0x41a376dbe3eebf96a8f49dff52f12e1c),
    ];

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
    fn test_gf2p128_add() {
        for _ in 0..1000 {
            let val1 = GF2p128::random(thread_rng());
            assert_eq!(val1 + GF2p128::ZERO, val1);
            assert_eq!(val1 + val1, GF2p128::ZERO);
            let val2 = GF2p128::random(thread_rng());
            assert_eq!(val1 + val2, val2 + val1);
        }

        let values = &GF2P128_VALUES;
        assert_eq!(values[0] + values[0], GF2p128::from_u128(0x00));
        assert_eq!(
            values[0] + values[1],
            GF2p128::from_u128(0xa8b8e7c19c61c0e8f6a1068ba6bf8545)
        );
        assert_eq!(
            values[0] + values[2],
            GF2p128::from_u128(0x1fa250293e57aa8f4d3a08e43b078878)
        );
        assert_eq!(
            values[0] + values[3],
            GF2p128::from_u128(0x20cf95f25b4059206fa614f1fc34d3d6)
        );
        assert_eq!(
            values[1] + values[0],
            GF2p128::from_u128(0xa8b8e7c19c61c0e8f6a1068ba6bf8545)
        );
        assert_eq!(values[1] + values[1], GF2p128::from_u128(0x00));
        assert_eq!(
            values[1] + values[2],
            GF2p128::from_u128(0xb71ab7e8a2366a67bb9b0e6f9db80d3d)
        );
        assert_eq!(
            values[1] + values[3],
            GF2p128::from_u128(0x88777233c72199c89907127a5a8b5693)
        );
        assert_eq!(
            values[2] + values[0],
            GF2p128::from_u128(0x1fa250293e57aa8f4d3a08e43b078878)
        );
        assert_eq!(
            values[2] + values[1],
            GF2p128::from_u128(0xb71ab7e8a2366a67bb9b0e6f9db80d3d)
        );
        assert_eq!(values[2] + values[2], GF2p128::from_u128(0x00));
        assert_eq!(
            values[2] + values[3],
            GF2p128::from_u128(0x3f6dc5db6517f3af229c1c15c7335bae)
        );
        assert_eq!(
            values[3] + values[0],
            GF2p128::from_u128(0x20cf95f25b4059206fa614f1fc34d3d6)
        );
        assert_eq!(
            values[3] + values[1],
            GF2p128::from_u128(0x88777233c72199c89907127a5a8b5693)
        );
        assert_eq!(
            values[3] + values[2],
            GF2p128::from_u128(0x3f6dc5db6517f3af229c1c15c7335bae)
        );
        assert_eq!(values[3] + values[3], GF2p128::from_u128(0x00));
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
    fn test_gf2p128_mul() {
        for _ in 0..1000 {
            let val1 = GF2p128::random(thread_rng());
            assert_eq!(val1 * GF2p128::ZERO, GF2p128::ZERO);
            assert_eq!(val1 * GF2p128::ONE, val1);
            let val2 = GF2p128::random(thread_rng());
            assert_eq!(val1 * val2, val2 * val1);
        }

        let values = &GF2P128_VALUES;
        assert_eq!(
            values[0] * values[0],
            GF2p128::from_u128(0x6e1625dc144c2cbda9163eeff9c1cf70)
        );
        assert_eq!(
            values[0] * values[1],
            GF2p128::from_u128(0xc9b21ac4a8e59b89f625fac9c2fe74de)
        );
        assert_eq!(
            values[0] * values[2],
            GF2p128::from_u128(0xc8a38099a47f698ed723268533f955be)
        );
        assert_eq!(
            values[0] * values[3],
            GF2p128::from_u128(0xd1fa610a3764b149244bdf5b76d73432)
        );
        assert_eq!(
            values[1] * values[0],
            GF2p128::from_u128(0xc9b21ac4a8e59b89f625fac9c2fe74de)
        );
        assert_eq!(
            values[1] * values[1],
            GF2p128::from_u128(0x556aac35c047bc7b8c72cdb250d2e20b)
        );
        assert_eq!(
            values[1] * values[2],
            GF2p128::from_u128(0xbc764270e1f03836f279ca879a014c5d)
        );
        assert_eq!(
            values[1] * values[3],
            GF2p128::from_u128(0x331876535842a5e5149919956bbf7f5d)
        );
        assert_eq!(
            values[2] * values[0],
            GF2p128::from_u128(0xc8a38099a47f698ed723268533f955be)
        );
        assert_eq!(
            values[2] * values[1],
            GF2p128::from_u128(0xbc764270e1f03836f279ca879a014c5d)
        );
        assert_eq!(
            values[2] * values[2],
            GF2p128::from_u128(0xf7304cb68d84d4cc8907074dd29b1c4d)
        );
        assert_eq!(
            values[2] * values[3],
            GF2p128::from_u128(0xfcf856eba6540a75f8b10ed6c197b877)
        );
        assert_eq!(
            values[3] * values[0],
            GF2p128::from_u128(0xd1fa610a3764b149244bdf5b76d73432)
        );
        assert_eq!(
            values[3] * values[1],
            GF2p128::from_u128(0x331876535842a5e5149919956bbf7f5d)
        );
        assert_eq!(
            values[3] * values[2],
            GF2p128::from_u128(0xfcf856eba6540a75f8b10ed6c197b877)
        );
        assert_eq!(
            values[3] * values[3],
            GF2p128::from_u128(0xde61339d54458418e5ab0957830f836e)
        );

        assert_eq!(
            GF2p128::from_u128(0x87c0fd359cea490832952ab4744fb570)
                * GF2p128::from_u128(0x1035fddaca9982ebb85b3e30054be76b),
            GF2p128::from_u128(0xf044f80c9fdbf7ec297ccd3d28817f57)
        );
    }

    #[test]
    fn test_gf2p8_neg() {
        for _ in 0..1000 {
            let val1 = GF2p8::random(thread_rng());
            assert_eq!(-val1, val1);
        }
    }

    #[test]
    fn test_gf2p128_neg() {
        for _ in 0..1000 {
            let val1 = GF2p128::random(thread_rng());
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
    fn test_gf2p128_sub() {
        for _ in 0..1000 {
            let val1 = GF2p128::random(thread_rng());
            let val2 = GF2p128::random(thread_rng());
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
    fn test_gf2p128_inv() {
        assert!(<Choice as Into<bool>>::into(
            GF2p128::ZERO.invert().is_none()
        ));
        assert_eq!(GF2p128::ONE.invert().unwrap(), GF2p128::ONE);

        for _ in 0..100 {
            let val1 = GF2p128::random(thread_rng());
            if val1 == GF2p128::ZERO {
                continue;
            }
            assert_eq!(val1 * val1.invert().unwrap(), GF2p128::ONE);
        }

        let values = &GF2P128_VALUES;
        assert_eq!(
            values[0].invert().unwrap(),
            GF2p128::from_u128(0xc876170a464ab346687c36e1b1864093)
        );
        assert_eq!(
            values[1].invert().unwrap(),
            GF2p128::from_u128(0xbcf64c0569e68857b0bfe6c00a23c877)
        );
        assert_eq!(
            values[2].invert().unwrap(),
            GF2p128::from_u128(0xb5f8fca4aa516a279b530e6da097717a)
        );
        assert_eq!(
            values[3].invert().unwrap(),
            GF2p128::from_u128(0x78ccdee85419f2d3a7e9df6e082fcf8a)
        );
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

    #[test]
    fn test_gf2p128_inner_product() {
        for _ in 0..100 {
            let xs: Vec<_> = (0..3).map(|_| GF2p128::random(thread_rng())).collect();
            let ys: Vec<_> = (0..3).map(|_| GF2p128::random(thread_rng())).collect();
            let expected = xs[0] * ys[0] + xs[1] * ys[1] + xs[2] * ys[2];
            let result: GF2p128 = InnerProduct::inner_product(xs.iter(), ys.iter());
            assert_eq!(result, expected);
        }
    }

    #[test]
    fn test_gf2p128_gf2p8_embedding() {
        for _ in 0..100 {
            let val1 = GF2p8::random(thread_rng());
            let val2 = GF2p8::random(thread_rng());
            assert_eq!(
                GF2p128::embed_gf2p8(val1) + GF2p128::embed_gf2p8(val2),
                GF2p128::embed_gf2p8(val1 + val2)
            );
            assert_eq!(
                GF2p128::embed_gf2p8(val1) * GF2p128::embed_gf2p8(val2),
                GF2p128::embed_gf2p8(val1 * val2)
            );
            if val1 != GF2p8::ZERO {
                assert_eq!(
                    GF2p128::embed_gf2p8(val1).invert().unwrap(),
                    GF2p128::embed_gf2p8(val1.invert().unwrap())
                );
            }
        }

        let values = &GF2P8_VALUES;
        assert_eq!(GF2p128::embed_gf2p8(GF2p8::ZERO), GF2p128::ZERO);
        assert_eq!(GF2p128::embed_gf2p8(GF2p8::ONE), GF2p128::ONE);
        assert_eq!(
            GF2p128::embed_gf2p8(values[0]),
            GF2p128::from_u128(0xab13c0d207db3b7ae159c432750ff703)
        );
        assert_eq!(
            GF2p128::embed_gf2p8(values[1]),
            GF2p128::from_u128(0x5723cd0c9078c01f5d026228ad8a0647)
        );
        assert_eq!(
            GF2p128::embed_gf2p8(values[2]),
            GF2p128::from_u128(0x452b59ebfc9f10876577c4b7f80fd2e7)
        );
        assert_eq!(
            GF2p128::embed_gf2p8(values[3]),
            GF2p128::from_u128(0x9add2badd1a222961ed46ceef7306f0c)
        );
    }
}
