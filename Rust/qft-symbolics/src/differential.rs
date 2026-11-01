//! Bounded exact polynomial differential operators with commuting right partials.
use crate::{
    computation::{Budget, Error, Resource, Result},
    exact::{BigC, BigQ},
    polynomial::BigPoly,
};
use std::collections::BTreeMap;

pub const MAX_DIMENSION: usize = 16;
pub const MAX_DERIVATIVE_ORDER: usize = 32;
fn invalid(message: &'static str) -> Error {
    Error::UnsupportedDomain { message }
}

fn poly_size(p: &BigPoly, b: &mut Budget) -> Result<(u64, u64)> {
    b.check(Resource::Terms, p.term_count() as u64)?;
    let mut bytes = 0u64;
    for (m, c) in p.terms() {
        b.charge(1)?;
        b.check(Resource::Degree, m.len() as u64)?;
        b.check(Resource::CoefficientBits, c.bits())?;
        for v in m {
            b.charge(1)?;
            bytes = bytes.saturating_add(v.len() as u64);
            b.check(Resource::InputBytes, bytes)?;
        }
    }
    Ok((p.term_count() as u64, bytes))
}

/// Exact first partial. Symbols other than `variable` are constant parameters.
pub fn differentiate(p: &BigPoly, variable: &str, b: &mut Budget) -> Result<BigPoly> {
    b.charge(1)?;
    b.check(Resource::InputBytes, variable.len() as u64)?;
    poly_size(p, b)?;
    let mut out = BigPoly::default();
    for (m, c) in p.terms() {
        b.charge(m.len() as u64 + 1)?;
        let count = m.iter().filter(|v| v.as_str() == variable).count();
        if count == 0 {
            continue;
        }
        let count =
            i64::try_from(count).map_err(|_| invalid("monomial multiplicity exceeds i64"))?;
        let coefficient = c.checked_mul(&BigC::from(count), b)?;
        let mut term = BigPoly::try_constant(coefficient, b)?;
        let mut skipped = false;
        for v in m {
            if !skipped && v == variable {
                skipped = true;
                continue;
            }
            term = term.try_mul(&BigPoly::try_var(v, b)?, b)?;
        }
        out = out.try_add(&term, b)?;
    }
    Ok(out)
}

/// A sum a_alpha(x) D^alpha. `compose(rhs)` acts on rhs first.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DifferentialOperator<const N: usize> {
    coords: [String; N],
    terms: BTreeMap<[u16; N], BigPoly>,
}
impl<const N: usize> DifferentialOperator<N> {
    fn context(coords: &[String; N], b: &mut Budget) -> Result<u64> {
        if N == 0 || N > MAX_DIMENSION {
            return Err(invalid("unsupported differential dimension"));
        }
        b.charge(1)?;
        let mut bytes = 0u64;
        for (i, c) in coords.iter().enumerate() {
            b.charge((i + 1) as u64)?;
            if c.is_empty() || coords[..i].contains(c) {
                return Err(invalid("empty or colliding coordinates"));
            }
            bytes = bytes.saturating_add(c.len() as u64);
            b.check(Resource::InputBytes, bytes)?;
        }
        Ok(bytes)
    }
    fn order(a: &[u16; N], b: &mut Budget) -> Result<()> {
        b.charge(N as u64 + 1)?;
        let n = a.iter().map(|v| *v as usize).sum::<usize>();
        if n > MAX_DERIVATIVE_ORDER {
            return Err(invalid("hard derivative order exceeded"));
        }
        b.check(Resource::Degree, n as u64)?;
        b.check(Resource::WordLength, n as u64)
    }
    fn validate(&self, b: &mut Budget) -> Result<()> {
        let mut bytes = Self::context(&self.coords, b)?;
        let mut nodes = self.terms.len() as u64;
        b.check(Resource::Terms, nodes)?;
        for (a, p) in &self.terms {
            Self::order(a, b)?;
            let (n, s) = poly_size(p, b)?;
            nodes = nodes.saturating_add(n);
            bytes = bytes.saturating_add(s);
            b.check(Resource::Terms, nodes)?;
            b.check(Resource::InputBytes, bytes)?;
        }
        Ok(())
    }
    fn compatible(&self, rhs: &Self, b: &mut Budget) -> Result<()> {
        b.charge(1)?;
        if self.coords != rhs.coords {
            return Err(invalid("differential coordinate contexts differ"));
        }
        self.validate(b)?;
        rhs.validate(b)
    }
    fn empty_like(&self, b: &mut Budget) -> Result<Self> {
        // Caller has validated before any context clone.
        b.charge(N as u64 + 1)?;
        Ok(Self {
            coords: self.coords.clone(),
            terms: BTreeMap::new(),
        })
    }
    fn insert(&mut self, a: [u16; N], p: BigPoly, b: &mut Budget) -> Result<()> {
        Self::order(&a, b)?;
        let p = if let Some(old) = self.terms.get(&a) {
            old.try_add(&p, b)?
        } else {
            p
        };
        self.check_slot(&a, &p, b)?;
        b.charge(1)?;
        if p.is_zero() {
            self.terms.remove(&a);
        } else {
            self.terms.insert(a, p);
        }
        Ok(())
    }
    // Validate the aggregate result before inserting a key or cloning a new coefficient.
    fn check_slot(&self, a: &[u16; N], p: &BigPoly, b: &mut Budget) -> Result<()> {
        let (pn, pb) = poly_size(p, b)?;
        let mut nodes = 0u64;
        let mut bytes = Self::context(&self.coords, b)?;
        for (key, value) in &self.terms {
            b.charge(1)?;
            if key == a {
                continue;
            }
            let (n, s) = poly_size(value, b)?;
            nodes = nodes.saturating_add(n + 1);
            bytes = bytes.saturating_add(s);
        }
        if !p.is_zero() {
            nodes = nodes.saturating_add(pn + 1);
            bytes = bytes.saturating_add(pb);
        }
        b.check(Resource::Terms, nodes)?;
        b.check(Resource::InputBytes, bytes)?;
        Ok(())
    }
    pub fn from_terms(
        coords: [String; N],
        terms: BTreeMap<[u16; N], BigPoly>,
        b: &mut Budget,
    ) -> Result<Self> {
        let mut out = Self { coords, terms };
        out.validate(b)?;
        out.terms.retain(|_, p| !p.is_zero());
        Ok(out)
    }
    pub fn coords(&self) -> &[String; N] {
        &self.coords
    }
    pub fn terms(&self) -> &BTreeMap<[u16; N], BigPoly> {
        &self.terms
    }
    pub fn partial(coords: [String; N], axis: usize, b: &mut Budget) -> Result<Self> {
        Self::context(&coords, b)?;
        if axis >= N {
            return Err(invalid("invalid differential axis"));
        }
        let mut a = [0; N];
        a[axis] = 1;
        Self::order(&a, b)?;
        // One derivative slot plus its constant monomial, before allocating either.
        b.check(Resource::Terms, 2)?;
        let mut out = Self {
            coords,
            terms: BTreeMap::new(),
        };
        let p = BigPoly::try_constant(1, b)?;
        out.insert(a, p, b)?;
        Ok(out)
    }
    pub fn multiplication(coords: [String; N], p: BigPoly, b: &mut Budget) -> Result<Self> {
        Self::context(&coords, b)?;
        let mut out = Self {
            coords,
            terms: BTreeMap::new(),
        };
        out.insert([0; N], p, b)?;
        Ok(out)
    }
    pub fn try_add(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.compatible(rhs, b)?;
        let mut out = self.empty_like(b)?;
        for (a, p) in self.terms.iter().chain(&rhs.terms) {
            if !out.terms.contains_key(a) {
                out.check_slot(a, p, b)?;
            }
            b.charge(p.term_count() as u64 + 1)?;
            out.insert(*a, p.clone(), b)?;
        }
        Ok(out)
    }
    pub fn try_sub(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.compatible(rhs, b)?;
        let minus = BigPoly::try_constant(-1, b)?;
        self.try_add(&rhs.try_scale(&minus, b)?, b)
    }
    /// Left multiplication of every coefficient by a polynomial.
    pub fn try_scale(&self, factor: &BigPoly, b: &mut Budget) -> Result<Self> {
        self.validate(b)?;
        poly_size(factor, b)?;
        let mut out = self.empty_like(b)?;
        for (a, p) in &self.terms {
            out.insert(*a, p.try_mul(factor, b)?, b)?;
        }
        Ok(out)
    }
    pub fn apply(&self, p: &BigPoly, b: &mut Budget) -> Result<BigPoly> {
        self.validate(b)?;
        poly_size(p, b)?;
        let mut out = BigPoly::default();
        for (a, c) in &self.terms {
            b.charge(p.term_count() as u64 + 1)?;
            let mut value = p.clone();
            for (axis, n) in a.iter().enumerate() {
                for _ in 0..*n {
                    value = differentiate(&value, &self.coords[axis], b)?;
                }
            }
            out = out.try_add(&c.try_mul(&value, b)?, b)?;
        }
        Ok(out)
    }
    pub fn compose(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.compatible(rhs, b)?;
        let mut out = self.empty_like(b)?;
        for (alpha, a) in &self.terms {
            for (beta, c) in &rhs.terms {
                // All intermediate derivative orders are at most this sum.
                let mut highest = [0; N];
                for i in 0..N {
                    highest[i] = alpha[i]
                        .checked_add(beta[i])
                        .ok_or_else(|| invalid("derivative index overflow"))?;
                }
                Self::order(&highest, b)?;
                let mut gamma = [0u16; N];
                loop {
                    b.charge(N as u64 + 1)?;
                    let mut index = [0; N];
                    b.charge(c.term_count() as u64 + 1)?;
                    let mut differentiated = c.clone();
                    let mut weight = BigQ::from(1);
                    for i in 0..N {
                        index[i] = alpha[i] - gamma[i] + beta[i];
                        for _ in 0..gamma[i] {
                            differentiated = differentiate(&differentiated, &self.coords[i], b)?;
                        }
                        // n<=32: binomial(n,k) <= 601080390, intermediates fit u64.
                        let mut choose = 1u64;
                        for j in 1..=gamma[i] {
                            b.charge(1)?;
                            choose = choose
                                .checked_mul((alpha[i] - j + 1) as u64)
                                .ok_or_else(|| invalid("binomial overflow"))?
                                / j as u64;
                        }
                        weight = weight.checked_mul(&BigQ::from(choose as i64), b)?;
                    }
                    if !differentiated.is_zero() {
                        let w = BigPoly::try_constant(BigC::new(weight, BigQ::from(0)), b)?;
                        let p = a.try_mul(&differentiated, b)?.try_mul(&w, b)?;
                        out.insert(index, p, b)?;
                    }
                    let mut i = 0;
                    while i < N && gamma[i] == alpha[i] {
                        gamma[i] = 0;
                        i += 1;
                    }
                    if i == N {
                        break;
                    }
                    gamma[i] += 1;
                }
            }
        }
        Ok(out)
    }
    pub fn commutator(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.compatible(rhs, b)?;
        self.compose(rhs, b)?.try_sub(&rhs.compose(self, b)?, b)
    }
}
impl DifferentialOperator<4> {
    pub fn momentum(coords: [String; 4], mu: usize, hbar: BigQ, b: &mut Budget) -> Result<Self> {
        Self::context(&coords, b)?;
        if mu >= 4 {
            return Err(invalid("invalid momentum index"));
        }
        b.check(Resource::CoefficientBits, hbar.bits())?;
        let im = if mu == 0 { hbar } else { hbar.checked_neg(b)? };
        let factor = BigPoly::try_constant(BigC::new(BigQ::from(0), im), b)?;
        Self::partial(coords, mu, b)?.try_scale(&factor, b)
    }
    pub fn orbital_lorentz(
        coords: [String; 4],
        mu: usize,
        nu: usize,
        hbar: BigQ,
        b: &mut Budget,
    ) -> Result<Self> {
        Self::context(&coords, b)?;
        if mu >= 4 || nu >= 4 {
            return Err(invalid("invalid Lorentz index"));
        }
        b.check(Resource::CoefficientBits, hbar.bits())?;
        b.charge(8)?;
        let pm = Self::momentum(coords.clone(), mu, hbar.clone(), b)?;
        let pn = Self::momentum(coords.clone(), nu, hbar, b)?;
        let xm = BigPoly::try_var(&coords[mu], b)?;
        let xn = BigPoly::try_var(&coords[nu], b)?;
        pn.try_scale(&xm, b)?.try_sub(&pm.try_scale(&xn, b)?, b)
    }
}
