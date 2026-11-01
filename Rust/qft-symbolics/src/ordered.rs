//! Bounded free associative algebra over commuting arbitrary precision polynomials.
use crate::{
    computation::{Budget, Error, Resource, Result},
    polynomial::BigPoly,
};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OrderedPolynomial(BTreeMap<Vec<String>, BigPoly>);
const HARD_LENGTH: u64 = 256;
fn bounded(b: &Budget, r: Resource, n: u64, hard: u64) -> Result<()> {
    b.check(r, n)?;
    if n > hard {
        return Err(Error::ResourceLimit {
            resource: r,
            limit: hard,
            required: n,
        });
    }
    Ok(())
}
fn name(s: &str, b: &mut Budget) -> Result<()> {
    b.charge((s.len() as u64).saturating_add(1))?;
    b.check(Resource::InputBytes, s.len() as u64)?;
    if s.is_empty() || s.chars().any(|c| c.is_whitespace() || c.is_control()) {
        return Err(Error::UnsupportedDomain {
            message: "names must be nonempty without whitespace or control characters",
        });
    }
    Ok(())
}
// Validate aggregate storage, including every coefficient monomial and symbol occurrence.
fn stats<'a>(
    items: impl Iterator<Item = (&'a [String], &'a BigPoly)>,
    b: &mut Budget,
) -> Result<()> {
    let (mut terms, mut bytes, mut words, mut bits) = (0u64, 0u64, 0u64, 0u64);
    for (w, p) in items {
        b.charge(1)?;
        bounded(b, Resource::WordLength, w.len() as u64, HARD_LENGTH)?;
        words = words.saturating_add(w.len() as u64);
        for s in w {
            name(s, b)?;
            bytes = bytes.saturating_add(s.len() as u64);
        }
        for (m, c) in p.terms() {
            b.charge(1)?;
            terms = terms.saturating_add(1);
            bounded(b, Resource::Degree, m.len() as u64, HARD_LENGTH)?;
            bits = bits.saturating_add(c.bits());
            b.check(Resource::CoefficientBits, bits)?;
            for s in m {
                name(s, b)?;
                bytes = bytes.saturating_add(s.len() as u64);
            }
        }
        b.check(Resource::Terms, terms)?;
        b.check(Resource::InputBytes, bytes)?;
    }
    // Aggregate word storage is bounded by both the term and per-word caps.
    b.check(
        Resource::Terms,
        words.div_ceil(b.limits().max_word_length.max(1) as u64),
    )?;
    Ok(())
}
impl OrderedPolynomial {
    pub fn terms(&self) -> &BTreeMap<Vec<String>, BigPoly> {
        &self.0
    }
    fn validate(&self, b: &mut Budget) -> Result<()> {
        stats(self.0.iter().map(|(w, p)| (w.as_slice(), p)), b)
    }
    fn insert(&mut self, w: &[String], p: BigPoly, b: &mut Budget) -> Result<()> {
        // Conservative temporary storage admission before coefficient merging or key cloning.
        stats(
            self.0
                .iter()
                .map(|(w, p)| (w.as_slice(), p))
                .chain(std::iter::once((w, &p))),
            b,
        )?;
        b.charge((w.len() as u64).saturating_add(1))?;
        let p = if let Some(old) = self.0.get(w) {
            old.try_add(&p, b)?
        } else {
            p
        };
        stats(
            self.0
                .iter()
                .filter(|(key, _)| key.as_slice() != w)
                .map(|(key, c)| (key.as_slice(), c))
                .chain(std::iter::once((w, &p))),
            b,
        )?;
        if p.is_zero() {
            self.0.remove(w);
        } else {
            self.0.insert(w.to_vec(), p);
        }
        Ok(())
    }
    pub fn scalar(p: BigPoly, b: &mut Budget) -> Result<Self> {
        stats(std::iter::once((&[][..], &p)), b)?;
        let mut out = Self::default();
        out.insert(&[], p, b)?;
        Ok(out)
    }
    pub fn generator(s: &str, b: &mut Budget) -> Result<Self> {
        name(s, b)?;
        bounded(b, Resource::WordLength, 1, HARD_LENGTH)?;
        b.check(Resource::Terms, 1)?;
        b.charge(1)?;
        let p = BigPoly::try_constant(1, b)?;
        Ok(Self(BTreeMap::from([(vec![s.to_owned()], p)])))
    }
    pub fn try_add(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.validate(b)?;
        rhs.validate(b)?;
        let mut out = Self::default();
        for (w, p) in self.0.iter().chain(&rhs.0) {
            stats(
                out.0
                    .iter()
                    .map(|(w, p)| (w.as_slice(), p))
                    .chain(std::iter::once((w.as_slice(), p))),
                b,
            )?;
            out.insert(w, p.try_add(&BigPoly::default(), b)?, b)?;
        }
        Ok(out)
    }
    pub fn try_neg(&self, b: &mut Budget) -> Result<Self> {
        self.validate(b)?;
        let mut out = Self::default();
        for (w, p) in &self.0 {
            out.insert(w, p.try_neg(b)?, b)?;
        }
        Ok(out)
    }
    pub fn try_sub(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.validate(b)?;
        rhs.validate(b)?;
        self.try_add(&rhs.try_neg(b)?, b)
    }
    pub fn try_mul(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.validate(b)?;
        rhs.validate(b)?;
        let mut out = Self::default();
        for (w, p) in &self.0 {
            for (v, q) in &rhs.0 {
                b.charge(1)?;
                let len = (w.len() as u64).saturating_add(v.len() as u64);
                bounded(b, Resource::WordLength, len, HARD_LENGTH)?;
                let bytes = w
                    .iter()
                    .chain(v)
                    .fold(0u64, |n, s| n.saturating_add(s.len() as u64));
                b.check(Resource::InputBytes, bytes)?;
                // Admit a conservative Cartesian coefficient expansion before multiplication.
                let count = (p.term_count() as u64).saturating_mul(q.term_count() as u64);
                b.check(
                    Resource::Terms,
                    out.0
                        .values()
                        .map(|x| x.term_count() as u64)
                        .sum::<u64>()
                        .saturating_add(count),
                )?;
                let old_bytes = out.0.iter().fold(0u64, |n, (k, c)| {
                    n.saturating_add(k.iter().map(|s| s.len() as u64).sum::<u64>())
                        .saturating_add(
                            c.terms()
                                .keys()
                                .flatten()
                                .map(|s| s.len() as u64)
                                .sum::<u64>(),
                        )
                });
                let coeff_bytes = p
                    .terms()
                    .keys()
                    .flatten()
                    .map(|s| s.len() as u64)
                    .sum::<u64>()
                    .saturating_mul(q.term_count() as u64)
                    .saturating_add(
                        q.terms()
                            .keys()
                            .flatten()
                            .map(|s| s.len() as u64)
                            .sum::<u64>()
                            .saturating_mul(p.term_count() as u64),
                    );
                b.check(
                    Resource::InputBytes,
                    old_bytes.saturating_add(bytes).saturating_add(coeff_bytes),
                )?;
                let old_bits = out
                    .0
                    .values()
                    .flat_map(|c| c.terms().values())
                    .fold(0u64, |n, c| n.saturating_add(c.bits()));
                let mut candidate_bits = old_bits;
                for c in p.terms().values() {
                    for d in q.terms().values() {
                        b.charge(1)?;
                        candidate_bits = candidate_bits.saturating_add(
                            c.bits()
                                .saturating_add(d.bits())
                                .saturating_add(1)
                                .saturating_mul(4),
                        );
                        b.check(Resource::CoefficientBits, candidate_bits)?;
                    }
                }
                for m in p.terms().keys() {
                    for n in q.terms().keys() {
                        b.charge(1)?;
                        bounded(
                            b,
                            Resource::Degree,
                            (m.len() as u64).saturating_add(n.len() as u64),
                            HARD_LENGTH,
                        )?;
                    }
                }
                b.charge(len)?;
                let c = p.try_mul(q, b)?;
                let key: Vec<_> = w.iter().chain(v).cloned().collect();
                out.insert(&key, c, b)?;
            }
        }
        Ok(out)
    }
    pub fn try_scale(&self, p: &BigPoly, b: &mut Budget) -> Result<Self> {
        self.validate(b)?;
        stats(std::iter::once((&[][..], p)), b)?;
        self.try_mul(&Self::scalar(p.try_add(&BigPoly::default(), b)?, b)?, b)
    }
    pub fn commutator(&self, rhs: &Self, b: &mut Budget) -> Result<Self> {
        self.try_mul(rhs, b)?.try_sub(&rhs.try_mul(self, b)?, b)
    }
    fn select(
        &self,
        parameters: &[&str],
        degree: usize,
        extract: bool,
        b: &mut Budget,
    ) -> Result<Self> {
        self.validate(b)?;
        let mut bytes = 0u64;
        for s in parameters {
            name(s, b)?;
            bytes = bytes.saturating_add(s.len() as u64);
            b.check(Resource::InputBytes, bytes)?;
        }
        let mut out = Self::default();
        for (w, p) in &self.0 {
            let mut coeff = BigPoly::default();
            for (m, c) in p.terms() {
                b.charge(
                    (m.len() as u64)
                        .saturating_mul(parameters.len() as u64)
                        .saturating_add(1),
                )?;
                let n = m
                    .iter()
                    .filter(|s| parameters.contains(&s.as_str()))
                    .count();
                if (extract && n != degree) || (!extract && n > degree) {
                    continue;
                }
                let mut term = BigPoly::try_constant(c.clone(), b)?;
                for s in m {
                    if !extract || !parameters.contains(&s.as_str()) {
                        term = term.try_mul(&BigPoly::try_var(s, b)?, b)?;
                    }
                }
                coeff = coeff.try_add(&term, b)?;
            }
            out.insert(w, coeff, b)?;
        }
        Ok(out)
    }
    pub fn truncate(&self, parameters: &[&str], max_degree: usize, b: &mut Budget) -> Result<Self> {
        self.select(parameters, max_degree, false, b)
    }
    pub fn coefficient(&self, parameter: &str, power: usize, b: &mut Budget) -> Result<Self> {
        self.select(&[parameter], power, true, b)
    }
    pub fn substitute(
        &self,
        replacements: &BTreeMap<String, Self>,
        b: &mut Budget,
    ) -> Result<Self> {
        self.validate(b)?;
        let mut bytes = 0u64;
        for (s, p) in replacements {
            name(s, b)?;
            bytes = bytes.saturating_add(s.len() as u64);
            b.check(Resource::InputBytes, bytes)?;
            p.validate(b)?;
            for (w, c) in &p.0 {
                bytes = bytes.saturating_add(w.iter().map(|x| x.len() as u64).sum::<u64>());
                bytes = bytes.saturating_add(
                    c.terms()
                        .keys()
                        .flatten()
                        .map(|x| x.len() as u64)
                        .sum::<u64>(),
                );
                b.check(Resource::InputBytes, bytes)?;
            }
        }
        stats(
            replacements
                .values()
                .flat_map(|p| p.0.iter().map(|(w, c)| (w.as_slice(), c))),
            b,
        )?;
        let mut out = Self::default();
        for (w, c) in &self.0 {
            let mut term = Self::scalar(c.try_add(&BigPoly::default(), b)?, b)?;
            for s in w {
                b.charge(1)?;
                term = if let Some(p) = replacements.get(s) {
                    term.try_mul(p, b)?
                } else {
                    term.try_mul(&Self::generator(s, b)?, b)?
                };
            }
            out = out.try_add(&term, b)?;
        }
        Ok(out)
    }
}
