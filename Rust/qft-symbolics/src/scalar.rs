//! Reduced rational arithmetic with checked i128 intermediates, and Q(i).
use std::{
    fmt,
    ops::{Add, Div, Mul, Neg, Sub},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Q {
    n: i128,
    d: i128,
}
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
fn mul(a: i128, b: i128) -> i128 {
    a.checked_mul(b).expect("rational overflow")
}
fn add(a: i128, b: i128) -> i128 {
    a.checked_add(b).expect("rational overflow")
}
impl Q {
    pub const ZERO: Self = Self { n: 0, d: 1 };
    pub const ONE: Self = Self { n: 1, d: 1 };
    pub fn new(n: i128, d: i128) -> Self {
        assert!(d != 0, "zero rational denominator");
        // Exclude MIN so abs/negation is representable in all subsequent steps.
        assert!(n != i128::MIN && d != i128::MIN, "rational range exceeded");
        let g = gcd(n.abs(), d.abs());
        let (n, d) = (n / g, d / g);
        if d < 0 {
            Self { n: -n, d: -d }
        } else {
            Self { n, d }
        }
    }
    pub fn sign(self) -> i8 {
        self.n.signum() as i8
    }
    pub fn is_zero(self) -> bool {
        self.n == 0
    }
    pub fn pow(self, n: usize) -> Self {
        (0..n).fold(Self::ONE, |a, _| a * self)
    }
    pub fn numerator(self) -> i128 {
        self.n
    }
    pub fn denominator(self) -> i128 {
        self.d
    }
}
impl From<i64> for Q {
    fn from(n: i64) -> Self {
        Self { n: n as i128, d: 1 }
    }
}
impl Add for Q {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        let g = gcd(self.d, b.d);
        Self::new(
            add(mul(self.n, b.d / g), mul(b.n, self.d / g)),
            mul(self.d, b.d / g),
        )
    }
}
impl Neg for Q {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.n, self.d)
    }
}
impl Sub for Q {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl Mul for Q {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let g = gcd(self.n.abs(), b.d);
        let h = gcd(b.n.abs(), self.d);
        Self::new(mul(self.n / g, b.n / h), mul(self.d / h, b.d / g))
    }
}
impl Div for Q {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        assert!(!b.is_zero(), "division by zero");
        self * Self::new(b.d, b.n)
    }
}
impl fmt::Display for Q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.d == 1 {
            write!(f, "{}", self.n)
        } else {
            write!(f, "{}/{}", self.n, self.d)
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct C {
    pub re: Q,
    pub im: Q,
}
impl C {
    pub const ZERO: Self = Self {
        re: Q::ZERO,
        im: Q::ZERO,
    };
    pub const ONE: Self = Self {
        re: Q::ONE,
        im: Q::ZERO,
    };
    pub const I: Self = Self {
        re: Q::ZERO,
        im: Q::ONE,
    };
    pub fn new(re: Q, im: Q) -> Self {
        Self { re, im }
    }
    pub fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
    pub fn is_zero(self) -> bool {
        self == Self::ZERO
    }
}
impl From<i64> for C {
    fn from(n: i64) -> Self {
        Self {
            re: n.into(),
            im: Q::ZERO,
        }
    }
}
impl From<Q> for C {
    fn from(re: Q) -> Self {
        Self { re, im: Q::ZERO }
    }
}
impl Add for C {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            re: self.re + b.re,
            im: self.im + b.im,
        }
    }
}
impl Neg for C {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}
impl Sub for C {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl Mul for C {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            re: self.re * b.re - self.im * b.im,
            im: self.re * b.im + self.im * b.re,
        }
    }
}
impl Div for C {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        let d = b.re * b.re + b.im * b.im;
        let n = self * b.conj();
        Self {
            re: n.re / d,
            im: n.im / d,
        }
    }
}
impl fmt::Display for C {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.im.is_zero() {
            write!(f, "{}", self.re)
        } else {
            write!(f, "({}+{}i)", self.re, self.im)
        }
    }
}

// Fallible compatibility API. Native checked intermediates allocate no heap.
// Like legacy operators, intermediate i128 overflow is an error even if a later
// cancellation could fit. Bit preflights and bit-weighted charges match BigQ;
// complex operations charge each rational step. Limits are not OS memory bounds.
use crate::computation::{Budget, Error, Resource, Result};
impl Default for Q {
    fn default() -> Self {
        Self::ZERO
    }
}
impl Default for C {
    fn default() -> Self {
        Self::ZERO
    }
}
fn scalar_preflight(b: &mut Budget, bits: u64, work: u64) -> Result<()> {
    b.check(Resource::CoefficientBits, bits)?;
    b.charge(work)
}
fn scalar_overflow() -> Error {
    Error::Overflow {
        operation: "checked rational intermediate",
    }
}
impl Q {
    pub fn bits(&self) -> u64 {
        (128 - self.n.unsigned_abs().leading_zeros())
            .max(128 - self.d.unsigned_abs().leading_zeros()) as u64
    }
    pub fn try_new(n: i128, d: i128) -> Result<Self> {
        if d == 0 {
            return Err(Error::ZeroDenominator);
        }
        // Unsigned magnitudes include abs(MIN); reduce before checking the
        // legacy symmetric range so reducible MIN inputs remain supported.
        let (mut a, mut b) = (n.unsigned_abs(), d.unsigned_abs());
        while b != 0 {
            (a, b) = (b, a % b);
        }
        let numerator = n.unsigned_abs() / a;
        let denominator = d.unsigned_abs() / a;
        if numerator > i128::MAX as u128 || denominator > i128::MAX as u128 {
            return Err(scalar_overflow());
        }
        let numerator = numerator as i128;
        Ok(Self {
            n: if (n < 0) != (d < 0) {
                -numerator
            } else {
                numerator
            },
            d: denominator as i128,
        })
    }
    fn checked_binary(&self, r: &Self, b: &mut Budget, op: u8) -> Result<Self> {
        if op == 3 && r.is_zero() {
            return Err(Error::DivisionByZero);
        }
        let bits = self.bits() + r.bits() + u64::from(op < 2);
        scalar_preflight(b, bits, (bits * bits * 8).max(1))?;
        let (n, d) = match op {
            0 | 1 => {
                let g = gcd(self.d, r.d);
                let a = self.n.checked_mul(r.d / g).ok_or_else(scalar_overflow)?;
                let c = r.n.checked_mul(self.d / g).ok_or_else(scalar_overflow)?;
                let n = if op == 0 {
                    a.checked_add(c)
                } else {
                    a.checked_sub(c)
                }
                .ok_or_else(scalar_overflow)?;
                (n, self.d.checked_mul(r.d / g).ok_or_else(scalar_overflow)?)
            }
            2 => {
                let g = gcd(self.n.abs(), r.d);
                let h = gcd(r.n.abs(), self.d);
                (
                    (self.n / g)
                        .checked_mul(r.n / h)
                        .ok_or_else(scalar_overflow)?,
                    (self.d / h)
                        .checked_mul(r.d / g)
                        .ok_or_else(scalar_overflow)?,
                )
            }
            _ => {
                let g = gcd(self.n.abs(), r.n.abs());
                let h = gcd(self.d, r.d);
                (
                    (self.n / g)
                        .checked_mul(r.d / h)
                        .ok_or_else(scalar_overflow)?,
                    (self.d / h)
                        .checked_mul(r.n / g)
                        .ok_or_else(scalar_overflow)?,
                )
            }
        };
        Self::try_new(n, d)
    }
    pub fn checked_add(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.checked_binary(r, b, 0)
    }
    pub fn checked_sub(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.checked_binary(r, b, 1)
    }
    pub fn checked_mul(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.checked_binary(r, b, 2)
    }
    pub fn checked_div(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.checked_binary(r, b, 3)
    }
    pub fn checked_neg(&self, b: &mut Budget) -> Result<Self> {
        scalar_preflight(b, self.bits(), self.bits().max(1))?;
        Ok(Self {
            n: self.n.checked_neg().ok_or_else(scalar_overflow)?,
            d: self.d,
        })
    }
}
impl C {
    pub fn bits(&self) -> u64 {
        self.re.bits().max(self.im.bits())
    }
    pub fn is_real(&self) -> bool {
        self.im.is_zero()
    }
    pub fn conjugate(&self) -> Self {
        self.conj()
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
        let conjugate = Self::new(r.re, r.im.checked_neg(b)?);
        let n = self.checked_mul(&conjugate, b)?;
        Ok(Self::new(
            n.re.checked_div(&d, b)?,
            n.im.checked_div(&d, b)?,
        ))
    }
}
impl Q {
    pub fn checked_pow(&self, mut n: usize, b: &mut Budget) -> Result<Self> {
        scalar_preflight(b, self.bits(), 1)?;
        let mut a = Self::ONE;
        let mut x = *self;
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
impl C {
    pub fn checked_pow(&self, mut n: usize, b: &mut Budget) -> Result<Self> {
        scalar_preflight(b, self.bits(), 1)?;
        let mut a = Self::ONE;
        let mut x = *self;
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
impl crate::coefficient::Coefficient for C {
    fn from_i64(n: i64) -> Self {
        n.into()
    }
    fn is_zero(&self) -> bool {
        C::is_zero(*self)
    }
    fn is_real(&self) -> bool {
        C::is_real(self)
    }
    fn conjugate(&self) -> Self {
        self.conj()
    }
    fn bits(&self) -> u64 {
        C::bits(self)
    }
    fn checked_add(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        C::checked_add(self, r, b)
    }
    fn checked_sub(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        C::checked_sub(self, r, b)
    }
    fn checked_mul(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        C::checked_mul(self, r, b)
    }
    fn checked_div(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        C::checked_div(self, r, b)
    }
    fn checked_neg(&self, b: &mut Budget) -> Result<Self> {
        C::checked_neg(self, b)
    }
}
