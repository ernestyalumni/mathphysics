//! Sparse commuting polynomials over Q(i). Variables are opaque strings.
//! No quotient ideals, dummy-index semantics, or general factorisation.
use crate::scalar::C;
use std::{
    collections::BTreeMap,
    fmt,
    ops::{Add, Mul, Neg, Sub},
};
const MAX_TERMS: usize = 8192;
const MAX_DEGREE: usize = 64;
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Poly(BTreeMap<Vec<String>, C>);
impl Poly {
    pub fn constant(c: impl Into<C>) -> Self {
        let mut p = Self::default();
        p.insert(vec![], c.into());
        p
    }
    pub fn var(s: impl Into<String>) -> Self {
        Self(BTreeMap::from([(vec![s.into()], C::ONE)]))
    }
    fn insert(&mut self, key: Vec<String>, c: C) {
        assert!(key.len() <= MAX_DEGREE, "polynomial degree budget");
        let c = self.0.get(&key).copied().unwrap_or(C::ZERO) + c;
        if c.is_zero() {
            self.0.remove(&key);
        } else {
            self.0.insert(key, c);
        }
        assert!(self.0.len() <= MAX_TERMS, "polynomial term budget");
    }
    /// Conjugate coefficients only; variables are left unchanged. This equals
    /// complex conjugation of a polynomial value only when variables are real.
    pub fn conjugate_coefficients(&self) -> Self {
        Self(self.0.iter().map(|(m, c)| (m.clone(), c.conj())).collect())
    }
    /// Variables remain formal; this does not assert that their values are real.
    pub fn has_real_coefficients(&self) -> bool {
        self.0.values().all(|c| c.im.is_zero())
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    pub fn term_count(&self) -> usize {
        self.0.len()
    }
    pub fn pow(&self, n: usize) -> Self {
        (0..n).fold(Self::constant(1), |a, _| a * self.clone())
    }
    pub fn eval(&self, values: &BTreeMap<String, C>) -> Result<C, String> {
        let mut out = C::ZERO;
        for (vars, c) in &self.0 {
            let mut t = *c;
            for v in vars {
                t = t * values
                    .get(v)
                    .copied()
                    .ok_or_else(|| format!("missing variable {v}"))?;
            }
            out = out + t;
        }
        Ok(out)
    }
}
impl Add for Poly {
    type Output = Self;
    fn add(mut self, b: Self) -> Self {
        for (m, c) in b.0 {
            self.insert(m, c);
        }
        self
    }
}
impl Neg for Poly {
    type Output = Self;
    fn neg(self) -> Self {
        Self(self.0.into_iter().map(|(m, c)| (m, -c)).collect())
    }
}
impl Sub for Poly {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl Mul for Poly {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let mut out = Self::default();
        for (a, c) in &self.0 {
            for (b, d) in &b.0 {
                let mut m = a.clone();
                m.extend(b.iter().cloned());
                m.sort();
                out.insert(m, *c * *d);
            }
        }
        out
    }
}
impl fmt::Display for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        for (i, (vars, c)) in self.0.iter().enumerate() {
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "{c}")?;
            for v in vars {
                write!(f, "*{v}")?;
            }
        }
        Ok(())
    }
}

/// Rational functions retain explicit denominator polynomials.
/// Equality by cross multiplication is valid only on the common nonzero domain.
/// The nonzero-denominator invariant cannot be bypassed by callers:
/// ```compile_fail
/// use qft_symbolics::polynomial::{Poly, Ratio};
/// let mut r = Ratio::new(Poly::constant(0), Poly::constant(1)).unwrap();
/// r.denominator = Poly::default();
/// ```
/// ```compile_fail
/// use qft_symbolics::polynomial::{Poly, Ratio};
/// let r = Ratio { numerator: Poly::default(), denominator: Poly::default() };
/// ```
#[derive(Clone, Debug)]
pub struct Ratio {
    numerator: Poly,
    denominator: Poly,
}
impl Ratio {
    pub fn new(numerator: Poly, denominator: Poly) -> Result<Self, String> {
        if denominator.is_zero() {
            return Err("identically zero denominator".into());
        }
        Ok(Self {
            numerator,
            denominator,
        })
    }
    pub fn numerator(&self) -> &Poly {
        &self.numerator
    }
    pub fn denominator(&self) -> &Poly {
        &self.denominator
    }
    pub fn residual(&self, b: &Self) -> Poly {
        self.numerator.clone() * b.denominator.clone()
            - b.numerator.clone() * self.denominator.clone()
    }
    pub fn eval(&self, values: &BTreeMap<String, C>) -> Result<C, String> {
        let d = self.denominator.eval(values)?;
        if d.is_zero() {
            return Err("pole: denominator vanishes".into());
        }
        Ok(self.numerator.eval(values)? / d)
    }
}
