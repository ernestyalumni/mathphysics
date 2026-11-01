//! Sparse commuting polynomials with checked shared coefficient arithmetic.
//! Legacy operators/constructors use `expect` and may panic for compatibility.
//! Budgets bound algebraic work, not allocator OOM or elapsed time.
//! InputBytes bounds aggregate UTF-8 symbol bytes in each polynomial separately,
//! counting every variable occurrence, and aggregate names in an evaluation map.
//! Inputs and intermediate/output polynomials each obey that cap; it does not
//! bound combined live heap size or cumulative allocation. Operations count
//! logical algebra/validation steps, not individual string bytes.
use crate::{
    coefficient::Coefficient,
    computation::{Budget, Error, Resource, Result as Checked},
    exact::BigC,
    scalar::C,
};
use std::{
    collections::BTreeMap,
    fmt,
    ops::{Add, Mul, Neg, Sub},
};
/// Compatibility wrappers retain historical term/degree limits without adding
/// new global work, string-byte, or coefficient-bit limits. Use try_* to bound work.
pub(crate) fn legacy_budget() -> Budget {
    Budget::new(crate::computation::Limits {
        max_operations: u64::MAX,
        max_coefficient_bits: u64::MAX,
        max_input_bytes: usize::MAX,
        ..crate::computation::Limits::default()
    })
}
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Polynomial<A: Coefficient>(BTreeMap<Vec<String>, A>);
pub type Poly = Polynomial<C>;
pub type BigPoly = Polynomial<BigC>;
impl<A: Coefficient> Polynomial<A> {
    pub fn terms(&self) -> &BTreeMap<Vec<String>, A> {
        &self.0
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    pub fn term_count(&self) -> usize {
        self.0.len()
    }
    fn validate(&self, b: &mut Budget) -> Checked<u64> {
        b.charge(1)?;
        b.check(Resource::Terms, self.0.len() as u64)?;
        let mut bytes = 0u64;
        for (m, c) in &self.0 {
            b.charge(1)?;
            b.check(Resource::Degree, m.len() as u64)?;
            b.check(Resource::CoefficientBits, c.bits())?;
            for v in m {
                b.charge(1)?;
                bytes = bytes.saturating_add(v.len() as u64);
            }
            b.check(Resource::InputBytes, bytes)?;
        }
        Ok(bytes)
    }
    pub fn try_constant(c: impl Into<A>, b: &mut Budget) -> Checked<Self> {
        let c = c.into();
        b.charge(1)?;
        b.check(Resource::CoefficientBits, c.bits())?;
        if c.is_zero() {
            return Ok(Self::default());
        }
        b.check(Resource::Terms, 1)?;
        Ok(Self(BTreeMap::from([(vec![], c)])))
    }
    /// Checks the borrowed name before allocating its owned polynomial key.
    pub fn try_var(s: impl AsRef<str>, b: &mut Budget) -> Checked<Self> {
        let s = s.as_ref();
        Self::try_var_built(s.len(), || s.to_owned(), b)
    }
    /// All constructor resource checks precede the name-building closure.
    /// The crate-owned builder must return precisely `bytes` UTF-8 bytes.
    pub(crate) fn try_var_built(
        bytes: usize,
        build: impl FnOnce() -> String,
        b: &mut Budget,
    ) -> Checked<Self> {
        b.charge(1)?;
        b.check(Resource::InputBytes, bytes as u64)?;
        b.check(Resource::Degree, 1)?;
        b.check(Resource::Terms, 1)?;
        let c = A::from_i64(1);
        b.check(Resource::CoefficientBits, c.bits())?;
        let name = build();
        debug_assert_eq!(name.len(), bytes);
        Ok(Self(BTreeMap::from([(vec![name], c)])))
    }
    fn insert(&mut self, m: &[String], c: A, bytes: &mut u64, b: &mut Budget) -> Checked<()> {
        b.charge((m.len() as u64).saturating_add(1))?;
        self.insert_precharged(m, c, bytes, b)
    }
    fn insert_precharged(
        &mut self,
        m: &[String],
        c: A,
        bytes: &mut u64,
        b: &mut Budget,
    ) -> Checked<()> {
        let key_bytes = m
            .iter()
            .fold(0u64, |sum, v| sum.saturating_add(v.len() as u64));
        let present = self.0.contains_key(m);
        b.check(Resource::Degree, m.len() as u64)?;
        b.check(Resource::CoefficientBits, c.bits())?;
        let c = if let Some(old) = self.0.get(m) {
            old.checked_add(&c, b)?
        } else {
            c
        };
        b.check(Resource::CoefficientBits, c.bits())?;
        if c.is_zero() {
            if present {
                *bytes = bytes.saturating_sub(key_bytes);
            }
            self.0.remove(m);
        } else {
            let required = if present {
                *bytes
            } else {
                bytes.saturating_add(key_bytes)
            };
            b.check(Resource::InputBytes, required)?;
            b.check(
                Resource::Terms,
                (self.0.len() + usize::from(!self.0.contains_key(m))) as u64,
            )?;
            self.0.insert(m.to_vec(), c);
            *bytes = required;
        }
        Ok(())
    }
    pub fn try_add(&self, rhs: &Self, b: &mut Budget) -> Checked<Self> {
        self.validate(b)?;
        rhs.validate(b)?;
        let mut bytes = 0u64;
        let mut out = Self::default();
        for (m, c) in self.0.iter().chain(&rhs.0) {
            if !out.0.contains_key(m) {
                b.check(Resource::Terms, (out.0.len() + 1) as u64)?;
                let key_bytes = m
                    .iter()
                    .fold(0u64, |sum, v| sum.saturating_add(v.len() as u64));
                b.charge(m.len() as u64)?;
                b.check(Resource::InputBytes, bytes.saturating_add(key_bytes))?;
            }
            out.insert(m, c.clone(), &mut bytes, b)?;
        }
        Ok(out)
    }
    pub fn try_neg(&self, b: &mut Budget) -> Checked<Self> {
        self.validate(b)?;
        let mut out = Self::default();
        let mut bytes = 0u64;
        for (m, c) in &self.0 {
            let negated = c.checked_neg(b)?;
            out.insert(m, negated, &mut bytes, b)?;
        }
        Ok(out)
    }
    pub fn try_sub(&self, rhs: &Self, b: &mut Budget) -> Checked<Self> {
        self.validate(b)?;
        rhs.validate(b)?;
        self.try_add(&rhs.try_neg(b)?, b)
    }
    pub fn try_mul(&self, rhs: &Self, b: &mut Budget) -> Checked<Self> {
        self.validate(b)?;
        rhs.validate(b)?;
        let mut out = Self::default();
        let mut bytes = 0u64;
        for (m, c) in &self.0 {
            for (n, d) in &rhs.0 {
                b.charge(1)?;
                b.check(
                    Resource::Degree,
                    (m.len() as u64).saturating_add(n.len() as u64),
                )?;
                b.charge((m.len() + n.len()) as u64)?;
                let key_bytes = m
                    .iter()
                    .chain(n)
                    .fold(0u64, |sum, v| sum.saturating_add(v.len() as u64));
                b.check(Resource::InputBytes, key_bytes)?;
                // At a resource boundary, compare a borrowed sorted merge before
                // allocating the key, so merging existing terms remains possible.
                if out.0.len() == b.limits().max_terms
                    || bytes.saturating_add(key_bytes) > b.limits().max_input_bytes as u64
                {
                    let mut present = false;
                    for key in out.0.keys() {
                        b.charge(1)?;
                        let (mut left, mut right) = (m.iter().peekable(), n.iter().peekable());
                        let merged = std::iter::from_fn(|| match (left.peek(), right.peek()) {
                            (Some(a), Some(c)) if a <= c => left.next(),
                            (Some(_), Some(_)) => right.next(),
                            (Some(_), None) => left.next(),
                            (None, Some(_)) => right.next(),
                            (None, None) => None,
                        });
                        if key.iter().eq(merged) {
                            present = true;
                            break;
                        }
                    }
                    if !present {
                        b.check(Resource::Terms, (out.0.len() + 1) as u64)?;
                        b.check(Resource::InputBytes, bytes.saturating_add(key_bytes))?;
                    }
                }
                let value = c.checked_mul(d, b)?;
                // Reserve insertion work before allocating the merged owned key.
                b.charge(
                    (m.len() as u64)
                        .saturating_add(n.len() as u64)
                        .saturating_add(1),
                )?;
                let mut key = m.clone();
                key.extend(n.iter().cloned());
                key.sort();
                out.insert_precharged(&key, value, &mut bytes, b)?;
            }
        }
        Ok(out)
    }
    pub fn try_pow(&self, mut n: usize, b: &mut Budget) -> Checked<Self> {
        self.validate(b)?;
        if let Some(degree) = self.0.keys().map(Vec::len).max() {
            b.check(Resource::Degree, (degree as u64).saturating_mul(n as u64))?;
        }
        let mut out = Self::try_constant(A::from_i64(1), b)?;
        let mut base = self.clone();
        while n > 0 {
            b.charge(1)?;
            if n % 2 == 1 {
                out = out.try_mul(&base, b)?;
            }
            n /= 2;
            if n > 0 {
                base = base.try_mul(&base, b)?;
            }
        }
        Ok(out)
    }
    pub fn try_eval(&self, values: &BTreeMap<String, A>, b: &mut Budget) -> Checked<A> {
        self.validate(b)?;
        let mut name_bytes = 0u64;
        for (name, c) in values {
            b.charge(1)?;
            name_bytes = name_bytes.saturating_add(name.len() as u64);
            b.check(Resource::InputBytes, name_bytes)?;
            b.check(Resource::CoefficientBits, c.bits())?;
        }
        let mut out = A::from_i64(0);
        b.check(Resource::CoefficientBits, out.bits())?;
        for (m, c) in &self.0 {
            let mut t = c.clone();
            for v in m {
                b.charge(1)?;
                let x = values
                    .get(v)
                    .ok_or_else(|| Error::MissingVariable { name: v.clone() })?;
                t = t.checked_mul(x, b)?;
            }
            out = out.checked_add(&t, b)?;
        }
        Ok(out)
    }
    pub fn constant(c: impl Into<A>) -> Self {
        Self::try_constant(c, &mut legacy_budget()).expect("polynomial constant")
    }
    pub fn var(s: impl Into<String>) -> Self {
        Self::try_var(s.into(), &mut legacy_budget()).expect("polynomial variable")
    }
    pub fn pow(&self, n: usize) -> Self {
        self.try_pow(n, &mut legacy_budget())
            .expect("polynomial power")
    }
    pub fn eval(&self, v: &BTreeMap<String, A>) -> Result<A, String> {
        self.try_eval(v, &mut legacy_budget())
            .map_err(|e| e.to_string())
    }
    /// Conjugates coefficients only; variables remain formal.
    pub fn conjugate_coefficients(&self) -> Self {
        Self(
            self.0
                .iter()
                .map(|(m, c)| (m.clone(), c.conjugate()))
                .collect(),
        )
    }
    pub fn has_real_coefficients(&self) -> bool {
        self.0.values().all(Coefficient::is_real)
    }
}
impl<A: Coefficient> Add for Polynomial<A> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.try_add(&rhs, &mut legacy_budget())
            .expect("polynomial addition")
    }
}
impl<A: Coefficient> Sub for Polynomial<A> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.try_sub(&rhs, &mut legacy_budget())
            .expect("polynomial subtraction")
    }
}
impl<A: Coefficient> Mul for Polynomial<A> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.try_mul(&rhs, &mut legacy_budget())
            .expect("polynomial multiplication")
    }
}
impl<A: Coefficient> Neg for Polynomial<A> {
    type Output = Self;
    fn neg(self) -> Self {
        self.try_neg(&mut legacy_budget())
            .expect("polynomial negation")
    }
}
impl<A: Coefficient> fmt::Display for Polynomial<A> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        for (i, (m, c)) in self.0.iter().enumerate() {
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "{c}")?;
            for v in m {
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
    pub fn try_residual(&self, rhs: &Self, budget: &mut Budget) -> Checked<Poly> {
        let left = self.numerator.try_mul(&rhs.denominator, budget)?;
        let right = rhs.numerator.try_mul(&self.denominator, budget)?;
        left.try_sub(&right, budget)
    }
    pub fn residual(&self, b: &Self) -> Poly {
        self.try_residual(b, &mut legacy_budget())
            .expect("rational residual")
    }
    pub fn try_eval(&self, values: &BTreeMap<String, C>, budget: &mut Budget) -> Checked<C> {
        let d = self.denominator.try_eval(values, budget)?;
        let n = self.numerator.try_eval(values, budget)?;
        if d.is_zero() {
            return Err(Error::DivisionByZero);
        }
        n.checked_div(&d, budget)
    }
    pub fn eval(&self, values: &BTreeMap<String, C>) -> Result<C, String> {
        let d = self.denominator.eval(values)?;
        if d.is_zero() {
            return Err("pole: denominator vanishes".into());
        }
        Ok(self.numerator.eval(values)? / d)
    }
}
