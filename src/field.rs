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

pub trait UnreducedField {
    type ReducedField: Field;
    fn reduce(self) -> Self::ReducedField;
}

trait LazyMul {
    type Output: UnreducedField;
    fn lazy_mul(self, other: Self) -> Self::Output;
}

trait LazyField: Add<Self::UF> + Sub<Self::UF> + LazyMul<Output = Self::UF> {
    type UF: UnreducedField;
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct GF2p8(pub u8);

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct UnreducedGF2p8(pub u16);

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
    #[inline(always)]
    fn reduce(self) -> Self::ReducedField {
        GF2p8(GF2p8::reduce(self.0))
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
impl Add<&Self> for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn add(self, other: &Self) -> Self::Output {
        self + *other
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Add<UnreducedGF2p8> for GF2p8 {
    type Output = UnreducedGF2p8;
    #[inline(always)]
    fn add(self, other: UnreducedGF2p8) -> Self::Output {
        UnreducedGF2p8(self.0 as u16 ^ other.0)
    }
}
impl AddAssign for GF2p8 {
    #[inline(always)]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}
impl AddAssign<&Self> for GF2p8 {
    #[inline(always)]
    fn add_assign(&mut self, other: &Self) {
        *self = *self + other;
    }
}

impl Div for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn div(self, other: Self) -> Self::Output {
        if other == Self::ZERO {
            panic!("division by zero");
        }
        self * other.invert().unwrap()
    }
}
impl Div<&Self> for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn div(self, other: &Self) -> Self::Output {
        self / *other
    }
}
impl DivAssign for GF2p8 {
    #[inline(always)]
    fn div_assign(&mut self, other: Self) {
        *self = *self / other;
    }
}
impl DivAssign<&Self> for GF2p8 {
    #[inline(always)]
    fn div_assign(&mut self, other: &Self) {
        *self = *self / other;
    }
}

impl LazyMul for GF2p8 {
    type Output = UnreducedGF2p8;
    #[inline(always)]
    fn lazy_mul(self, other: Self) -> Self::Output {
        UnreducedGF2p8(clmul_u8(self.0, other.0))
    }
}

impl Mul for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: Self) -> Self::Output {
        Self(Self::reduce(clmul_u8(self.0, other.0)))
    }
}
impl Mul<&Self> for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: &Self) -> Self::Output {
        self * *other
    }
}
impl Mul<bool> for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, other: bool) -> Self::Output {
        Self(self.0 * other as u8)
    }
}
impl MulAssign for GF2p8 {
    #[inline(always)]
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}
impl MulAssign<&Self> for GF2p8 {
    #[inline(always)]
    fn mul_assign(&mut self, other: &Self) {
        *self = *self * other;
    }
}
impl MulAssign<bool> for GF2p8 {
    #[inline(always)]
    fn mul_assign(&mut self, other: bool) {
        self.0 *= other as u8;
    }
}

impl Neg for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        self
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, other: Self) -> Self::Output {
        self + other
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub<&Self> for GF2p8 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, other: &Self) -> Self::Output {
        self + other
    }
}
#[allow(clippy::suspicious_arithmetic_impl)]
impl Sub<UnreducedGF2p8> for GF2p8 {
    type Output = UnreducedGF2p8;
    #[inline(always)]
    fn sub(self, other: UnreducedGF2p8) -> Self::Output {
        self + other
    }
}
#[allow(clippy::suspicious_op_assign_impl)]
impl SubAssign for GF2p8 {
    #[inline(always)]
    fn sub_assign(&mut self, other: Self) {
        *self += other
    }
}
#[allow(clippy::suspicious_op_assign_impl)]
impl SubAssign<&Self> for GF2p8 {
    #[inline(always)]
    fn sub_assign(&mut self, other: &Self) {
        *self += other
    }
}

impl Sum for GF2p8 {
    #[inline(always)]
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self::ZERO, Add::add)
    }
}
impl<'a> Sum<&'a Self> for GF2p8 {
    #[inline(always)]
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.copied().sum()
    }
}

impl Product for GF2p8 {
    #[inline(always)]
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self::ONE, Mul::mul)
    }
}
impl<'a> Product<&'a Self> for GF2p8 {
    #[inline(always)]
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.copied().product()
    }
}

impl ConditionallySelectable for GF2p8 {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self(u8::conditional_select(&a.0, &b.0, choice))
    }
}

impl ConstantTimeEq for GF2p8 {
    #[inline(always)]
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

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

impl Zero for GF2p8 {
    fn zero() -> Self {
        Self::ZERO
    }
    fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }
}

impl One for GF2p8 {
    fn one() -> Self {
        Self::ONE
    }
    fn is_one(&self) -> bool {
        *self == Self::ONE
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
}
