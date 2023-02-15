use core::ops::{Add, Sub};
use ff::Field;

// pub use crate::gf2::{GF2Vector, GF2View};
// pub use crate::gf2p128::{GF2p128Fast, GF2p128Naive};
// pub use crate::gf2psmall::GF2p8;

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

pub trait BytesRepr {
    type Repr: AsRef<[u8]> + AsMut<[u8]> + IntoIterator<Item = u8>;
    fn to_repr(self) -> Self::Repr;
    fn from_repr(repr: Self::Repr) -> Self;
}

pub trait VecToGF2p128<F>: Sized {
    const VECTOR_SIZE: usize;
    fn convert(vec: &[Self]) -> F;
}

pub trait BitMulAccumulate: Sized {
    fn bitmul_accumulate(ys: &mut [Self], x: Self, bs: &[u8]);
}
