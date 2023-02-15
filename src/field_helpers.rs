macro_rules! impl_additional_field_arithmetic {
    ($field_name:ident) => {
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
                    .map(|(x, y)| Mul::mul(x, y))
                    .sum::<$field_name>()
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
