//! Exact rational and Gaussian rational coefficients.
//! Preflights conservatively bound component bits and bit-weighted arithmetic
//! work, not allocator usage or elapsed time. Cancellation may be rejected.
use crate::{
    coefficient::Coefficient,
    computation::{Budget, Error, Resource, Result},
    scalar::{C, Q},
};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{ToPrimitive, Zero};
use std::fmt;
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct BigQ(BigRational);
impl BigQ {
    pub fn new(n: BigInt, d: BigInt, b: &mut Budget) -> Result<Self> {
        if d.is_zero() {
            return Err(Error::ZeroDenominator);
        }
        let bits = n.bits().max(d.bits());
        preflight(b, bits, bits.saturating_mul(bits).max(1))?;
        Ok(Self(BigRational::new(n, d)))
    }
    pub fn parse(n: &str, d: &str, b: &mut Budget) -> Result<Self> {
        b.check(
            Resource::InputBytes,
            (n.len() as u64).saturating_add(d.len() as u64),
        )?;
        fn integer(s: &str, b: &mut Budget) -> Result<BigInt> {
            let digits = s
                .strip_prefix('-')
                .or_else(|| s.strip_prefix('+'))
                .unwrap_or(s);
            if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
                return Err(Error::InvalidInteger);
            }
            let bits = (digits.len() as u64).saturating_mul(4);
            preflight(b, bits, bits.saturating_mul(bits).max(1))?;
            BigInt::parse_bytes(s.as_bytes(), 10).ok_or(Error::InvalidInteger)
        }
        let n = integer(n, b)?;
        let d = integer(d, b)?;
        Self::new(n, d, b)
    }
    pub fn numerator(&self) -> &BigInt {
        self.0.numer()
    }
    pub fn denominator(&self) -> &BigInt {
        self.0.denom()
    }
    pub fn bits(&self) -> u64 {
        self.numerator().bits().max(self.denominator().bits())
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    pub fn try_to_small(&self) -> Result<Q> {
        let fail = || Error::Overflow {
            operation: "exact to small",
        };
        let n = self.numerator().to_i128().ok_or_else(fail)?;
        let d = self.denominator().to_i128().ok_or_else(fail)?;
        if n == i128::MIN || d == i128::MIN {
            return Err(fail());
        }
        Ok(Q::new(n, d))
    }
    fn binary(&self, r: &Self, b: &mut Budget, op: u8) -> Result<Self> {
        if op == 3 && r.is_zero() {
            return Err(Error::DivisionByZero);
        }
        let bits = self
            .bits()
            .saturating_add(r.bits())
            .saturating_add(if op < 2 { 1 } else { 0 });
        preflight(b, bits, bits.saturating_mul(bits).saturating_mul(8).max(1))?;
        Ok(Self(match op {
            0 => &self.0 + &r.0,
            1 => &self.0 - &r.0,
            2 => &self.0 * &r.0,
            _ => &self.0 / &r.0,
        }))
    }
    pub fn checked_add(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.binary(r, b, 0)
    }
    pub fn checked_sub(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.binary(r, b, 1)
    }
    pub fn checked_mul(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.binary(r, b, 2)
    }
    pub fn checked_div(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.binary(r, b, 3)
    }
    pub fn checked_neg(&self, b: &mut Budget) -> Result<Self> {
        preflight(b, self.bits(), self.bits().max(1))?;
        Ok(Self(-&self.0))
    }
    pub fn checked_pow(&self, mut n: usize, b: &mut Budget) -> Result<Self> {
        preflight(b, self.bits(), 1)?;
        let mut a = Self::from(1);
        let mut x = self.clone();
        while n > 0 {
            if n & 1 == 1 {
                a = a.checked_mul(&x, b)?
            }
            n >>= 1;
            if n > 0 {
                x = x.checked_mul(&x, b)?
            }
        }
        Ok(a)
    }
}
fn preflight(b: &mut Budget, bits: u64, work: u64) -> Result<()> {
    b.check(Resource::CoefficientBits, bits)?;
    b.charge(work)
}
impl From<i64> for BigQ {
    fn from(n: i64) -> Self {
        Self(BigRational::from_integer(n.into()))
    }
}
impl From<Q> for BigQ {
    fn from(q: Q) -> Self {
        Self(BigRational::new(
            q.numerator().into(),
            q.denominator().into(),
        ))
    }
}
impl fmt::Display for BigQ {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct BigC {
    pub re: BigQ,
    pub im: BigQ,
}
impl BigC {
    pub fn new(re: BigQ, im: BigQ) -> Self {
        Self { re, im }
    }
    pub fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
    pub fn is_real(&self) -> bool {
        self.im.is_zero()
    }
    pub fn conjugate(&self) -> Self {
        Self::new(self.re.clone(), BigQ(-&self.im.0))
    }
    pub fn bits(&self) -> u64 {
        self.re.bits().max(self.im.bits())
    }
    pub fn try_to_small(&self) -> Result<C> {
        Ok(C::new(self.re.try_to_small()?, self.im.try_to_small()?))
    }
    pub fn checked_add(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        Ok(Self::new(
            self.re.checked_add(&r.re, b)?,
            self.im.checked_add(&r.im, b)?,
        ))
    }
    pub fn checked_sub(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        Ok(Self::new(
            self.re.checked_sub(&r.re, b)?,
            self.im.checked_sub(&r.im, b)?,
        ))
    }
    pub fn checked_neg(&self, b: &mut Budget) -> Result<Self> {
        Ok(Self::new(self.re.checked_neg(b)?, self.im.checked_neg(b)?))
    }
    pub fn checked_mul(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        let re = self
            .re
            .checked_mul(&r.re, b)?
            .checked_sub(&self.im.checked_mul(&r.im, b)?, b)?;
        let im = self
            .re
            .checked_mul(&r.im, b)?
            .checked_add(&self.im.checked_mul(&r.re, b)?, b)?;
        Ok(Self::new(re, im))
    }
    pub fn checked_div(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        if r.is_zero() {
            return Err(Error::DivisionByZero);
        }
        let d =
            r.re.checked_mul(&r.re, b)?
                .checked_add(&r.im.checked_mul(&r.im, b)?, b)?;
        let conjugate = Self::new(r.re.clone(), r.im.checked_neg(b)?);
        let n = self.checked_mul(&conjugate, b)?;
        Ok(Self::new(
            n.re.checked_div(&d, b)?,
            n.im.checked_div(&d, b)?,
        ))
    }
    pub fn checked_pow(&self, mut n: usize, b: &mut Budget) -> Result<Self> {
        preflight(b, self.bits(), 1)?;
        let mut a = Self::from(1);
        let mut x = self.clone();
        while n > 0 {
            if n & 1 == 1 {
                a = a.checked_mul(&x, b)?
            }
            n >>= 1;
            if n > 0 {
                x = x.checked_mul(&x, b)?
            }
        }
        Ok(a)
    }
}
impl From<i64> for BigC {
    fn from(n: i64) -> Self {
        Self::new(n.into(), 0.into())
    }
}
impl From<Q> for BigC {
    fn from(q: Q) -> Self {
        Self::new(q.into(), 0.into())
    }
}
impl From<C> for BigC {
    fn from(c: C) -> Self {
        Self::new(c.re.into(), c.im.into())
    }
}
impl fmt::Display for BigC {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_real() {
            write!(f, "{}", self.re)
        } else {
            write!(f, "({}+{}i)", self.re, self.im)
        }
    }
}
impl Coefficient for BigC {
    fn from_i64(n: i64) -> Self {
        n.into()
    }
    fn is_zero(&self) -> bool {
        BigC::is_zero(self)
    }
    fn is_real(&self) -> bool {
        BigC::is_real(self)
    }
    fn conjugate(&self) -> Self {
        BigC::conjugate(self)
    }
    fn bits(&self) -> u64 {
        BigC::bits(self)
    }
    fn checked_add(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        BigC::checked_add(self, r, b)
    }
    fn checked_sub(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        BigC::checked_sub(self, r, b)
    }
    fn checked_mul(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        BigC::checked_mul(self, r, b)
    }
    fn checked_div(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        BigC::checked_div(self, r, b)
    }
    fn checked_neg(&self, b: &mut Budget) -> Result<Self> {
        BigC::checked_neg(self, b)
    }
}
