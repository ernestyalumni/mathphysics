use crate::Half;
use qft_symbolics::{
    polynomial::Poly,
    scalar::{C, Q},
};
use std::{
    fmt::{Debug, Display},
    ops::{Add, Div, Mul, Neg, Sub},
};
/// Commutative coefficient ring contract. Generic code assumes these laws;
/// Rust's trait system checks operations, not their mathematical axioms.
pub trait Scalar:
    Clone
    + PartialEq
    + Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
{
    fn from_i64(value: i64) -> Self;
    fn zero() -> Self {
        Self::from_i64(0)
    }
    fn one() -> Self {
        Self::from_i64(1)
    }
}
pub trait RationalScalar: Scalar {
    fn from_q(value: Q) -> Self;
}
pub trait Conjugate: Scalar {
    fn conjugate(&self) -> Self;
}
/// Monomorphized real numerical backend. Add an in-crate implementation or use
/// a newtype around an external half type (Rust orphan rules apply). No virtual
/// dispatch or runtime scalar selection is used.
pub trait Real: RationalScalar + Conjugate + Copy + PartialOrd + Div<Output = Self> {
    fn from_f64(value: f64) -> Self;
    fn to_f64(self) -> f64;
    fn epsilon() -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn hypot(self, b: Self) -> Self;
    fn abs(self) -> Self;
    fn is_finite(self) -> bool;
}
macro_rules! integer{($($t:ty),*)=>{$(
    impl Scalar for $t{fn from_i64(v:i64)->Self{v.try_into().expect("integer scalar range")}}
    impl Conjugate for $t{fn conjugate(&self)->Self{*self}}
)*};}
integer!(i32, i64);
macro_rules! real{($($t:ty),*)=>{$(
    impl Scalar for $t{fn from_i64(v:i64)->Self{v as Self}}
    impl RationalScalar for $t{fn from_q(q:Q)->Self{q.numerator() as Self/q.denominator() as Self}}
    impl Conjugate for $t{fn conjugate(&self)->Self{*self}}
    impl Real for $t{
        fn from_f64(v:f64)->Self{v as Self} fn to_f64(self)->f64{self as f64}
        fn epsilon()->Self{Self::EPSILON} fn sin(self)->Self{self.sin()}
        fn cos(self)->Self{self.cos()} fn hypot(self,b:Self)->Self{self.hypot(b)}
        fn abs(self)->Self{self.abs()} fn is_finite(self)->bool{self.is_finite()}
    }
)*};}
real!(f32, f64);
impl Scalar for Half {
    fn from_i64(v: i64) -> Self {
        Self::from_f64(v as f64)
    }
}
impl RationalScalar for Half {
    fn from_q(q: Q) -> Self {
        Self::from_f64(q.numerator() as f64 / q.denominator() as f64)
    }
}
impl Conjugate for Half {
    fn conjugate(&self) -> Self {
        *self
    }
}
impl Real for Half {
    fn from_f64(v: f64) -> Self {
        Self::from_f64(v)
    }
    fn to_f64(self) -> f64 {
        self.to_f64()
    }
    fn epsilon() -> Self {
        Self::from_bits(0x1400)
    }
    fn sin(self) -> Self {
        Self::from_f64(f64::from(self.to_f32().sin()))
    }
    fn cos(self) -> Self {
        Self::from_f64(f64::from(self.to_f32().cos()))
    }
    fn hypot(self, b: Self) -> Self {
        Self::from_f64(f64::from(self.to_f32().hypot(b.to_f32())))
    }
    fn abs(self) -> Self {
        Self::from_bits(self.to_bits() & 0x7fff)
    }
    fn is_finite(self) -> bool {
        self.to_f32().is_finite()
    }
}
impl Scalar for Q {
    fn from_i64(v: i64) -> Self {
        v.into()
    }
}
impl RationalScalar for Q {
    fn from_q(v: Q) -> Self {
        v
    }
}
impl Conjugate for Q {
    fn conjugate(&self) -> Self {
        *self
    }
}
impl Scalar for C {
    fn from_i64(v: i64) -> Self {
        v.into()
    }
}
impl RationalScalar for C {
    fn from_q(v: Q) -> Self {
        v.into()
    }
}
impl Conjugate for C {
    fn conjugate(&self) -> Self {
        self.conj()
    }
}
impl Scalar for Poly {
    fn from_i64(v: i64) -> Self {
        Self::constant(v)
    }
}
impl RationalScalar for Poly {
    fn from_q(v: Q) -> Self {
        Self::constant(v)
    }
}
/// Coefficient conjugation, fixing every formal variable. For the Lie group
/// interpretation these variables are real; this is not conjugation of
/// arbitrary complex-valued variable assignments.
impl Conjugate for Poly {
    fn conjugate(&self) -> Self {
        self.conjugate_coefficients()
    }
}
/// A real coefficient ring extended by a central i. Exact Q, Trig and numerical
/// real backends work alike; do not interpret complex-coefficient T as real.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Complex<T> {
    pub re: T,
    pub im: T,
}
impl<T: Scalar> Complex<T> {
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }
    pub fn real(re: T) -> Self {
        Self::new(re, T::zero())
    }
    pub fn conj(&self) -> Self
    where
        T: Conjugate,
    {
        Self::new(self.re.conjugate(), -self.im.conjugate())
    }
    /// re^2+im^2: squared modulus when the component coefficient ring is real.
    pub fn norm_squared(&self) -> T {
        self.re.clone() * self.re.clone() + self.im.clone() * self.im.clone()
    }
}
impl<T: Display> Display for Complex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({})+I*({})", self.re, self.im)
    }
}
impl<T: Scalar> Add for Complex<T> {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self::new(self.re + b.re, self.im + b.im)
    }
}
impl<T: Scalar> Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}
impl<T: Scalar> Sub for Complex<T> {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl<T: Scalar> Mul for Complex<T> {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self::new(
            self.re.clone() * b.re.clone() - self.im.clone() * b.im.clone(),
            self.re * b.im + self.im * b.re,
        )
    }
}
impl<T: Scalar> Scalar for Complex<T> {
    fn from_i64(v: i64) -> Self {
        Self::real(T::from_i64(v))
    }
}
impl<T: RationalScalar> RationalScalar for Complex<T> {
    fn from_q(v: Q) -> Self {
        Self::real(T::from_q(v))
    }
}
impl<T: Conjugate> Conjugate for Complex<T> {
    fn conjugate(&self) -> Self {
        self.conj()
    }
}
