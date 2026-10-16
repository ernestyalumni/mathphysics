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
