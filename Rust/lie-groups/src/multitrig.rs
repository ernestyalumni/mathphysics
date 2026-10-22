//! Exact tensor product of N unit-circle coordinate rings over Q:
//! Q[s_i,c_i]/(s_i²+c_i²-1). Independent angles only; no addition identities,
//! inverse trig, arbitrary function simplification, or complex-valued arguments.
//! Cosine powers reduce to 0/1. At most 8 angles, 4096 normal-form terms and
//! total sine degree 64; budget/checked-i128 exhaustion panics. This is a
//! focused exact coefficient ring, not a general CAS.
use crate::{scalar::Conjugate, Error, RationalScalar, Scalar};
use qft_symbolics::scalar::Q;
use std::{
    collections::BTreeMap,
    fmt,
    ops::{Add, Mul, Neg, Sub},
};
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiTrig<const N: usize> {
    terms: BTreeMap<([u16; N], u8), Q>,
}
impl<const N: usize> MultiTrig<N> {
    fn add_term(&mut self, key: ([u16; N], u8), value: Q) {
        let next = self.terms.get(&key).copied().unwrap_or(Q::ZERO) + value;
        if next == Q::ZERO {
            self.terms.remove(&key);
        } else {
            self.terms.insert(key, next);
        }
        assert!(self.terms.len() <= 4096, "MultiTrig term budget exceeded");
    }
    pub fn sin(index: usize) -> Result<Self, Error> {
        if N > 8 {
            return Err(Error::Domain(
                "MultiTrig supports at most 8 independent angles",
            ));
        }
        if index >= N {
            return Err(Error::IndexOutOfRange {
                index,
                dimension: N,
            });
        }
        let mut powers = [0; N];
        powers[index] = 1;
        Ok(Self {
            terms: BTreeMap::from([((powers, 0), Q::ONE)]),
        })
    }
    pub fn cos(index: usize) -> Result<Self, Error> {
        if N > 8 {
            return Err(Error::Domain(
                "MultiTrig supports at most 8 independent angles",
            ));
        }
        if index >= N {
            return Err(Error::IndexOutOfRange {
                index,
                dimension: N,
            });
        }
        Ok(Self {
            terms: BTreeMap::from([(([0; N], 1 << index), Q::ONE)]),
        })
    }
    /// Rational evaluation only, rejecting assignments off any unit circle.
    pub fn at(&self, s: [Q; N], c: [Q; N]) -> Result<Q, Error> {
        if (0..N).any(|i| s[i] * s[i] + c[i] * c[i] != Q::ONE) {
            return Err(Error::Domain("not a unit circle"));
        }
        let mut result = Q::ZERO;
        for ((powers, mask), coefficient) in &self.terms {
            let mut term = *coefficient;
            for i in 0..N {
                for _ in 0..powers[i] {
                    term = term * s[i];
                }
                if mask & (1 << i) != 0 {
                    term = term * c[i];
                }
            }
            result = result + term;
        }
        Ok(result)
    }
}
impl<const N: usize> Scalar for MultiTrig<N> {
    fn from_i64(v: i64) -> Self {
        Self::from_q(Q::from(v))
    }
}
impl<const N: usize> RationalScalar for MultiTrig<N> {
    fn from_q(v: Q) -> Self {
        assert!(N <= 8, "MultiTrig supports at most 8 angles");
        Self {
            terms: if v == Q::ZERO {
                BTreeMap::new()
            } else {
                BTreeMap::from([(([0; N], 0), v)])
            },
        }
    }
}
impl<const N: usize> Conjugate for MultiTrig<N> {
    fn conjugate(&self) -> Self {
        self.clone()
    }
}
impl<const N: usize> Add for MultiTrig<N> {
    type Output = Self;
    fn add(mut self, b: Self) -> Self {
        for (k, v) in b.terms {
            self.add_term(k, v);
        }
        self
    }
}
impl<const N: usize> Neg for MultiTrig<N> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            terms: self.terms.into_iter().map(|(k, v)| (k, -v)).collect(),
        }
    }
}
impl<const N: usize> Sub for MultiTrig<N> {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl<const N: usize> Mul for MultiTrig<N> {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let mut result = Self::zero();
        for ((a, am), av) in &self.terms {
            for ((b, bm), bv) in &b.terms {
                let powers = std::array::from_fn(|i| {
                    a[i].checked_add(b[i]).expect("MultiTrig degree overflow")
                });
                assert!(
                    powers.iter().map(|x| usize::from(*x)).sum::<usize>() <= 64,
                    "MultiTrig degree budget exceeded"
                );
                let overlap = am & bm;
                let mut branches = vec![(powers, av.to_owned() * bv.to_owned())];
                for i in 0..N {
                    if overlap & (1 << i) != 0 {
                        let mut extra = Vec::new();
                        for (p, v) in &branches {
                            let mut p2 = *p;
                            p2[i] = p2[i].checked_add(2).expect("MultiTrig degree overflow");
                            extra.push((p2, -*v));
                        }
                        branches.extend(extra);
                    }
                }
                for (p, v) in branches {
                    assert!(
                        p.iter().map(|x| usize::from(*x)).sum::<usize>() <= 64,
                        "MultiTrig degree budget exceeded"
                    );
                    result.add_term((p, am ^ bm), v);
                }
            }
        }
        result
    }
}

/// Fixed-symbol expression export for trusted CAS interchange. No parser is
/// provided; variables are s0,c0,... and powers use Python/SymPy ** syntax.
impl<const N: usize> fmt::Display for MultiTrig<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }
        for (index, ((powers, mask), coefficient)) in self.terms.iter().enumerate() {
            if index > 0 {
                write!(f, "+")?;
            }
            write!(f, "({coefficient})")?;
            for (i, power) in powers.iter().enumerate() {
                if *power > 0 {
                    write!(f, "*s{i}**{power}")?;
                }
                if mask & (1 << i) != 0 {
                    write!(f, "*c{i}")?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type T = MultiTrig<2>;
    #[test]
    fn independent_cosine_square_reductions_expand_all_branches() {
        let a = T::cos(0).unwrap();
        let b = T::cos(1).unwrap();
        let s = T::sin(0).unwrap();
        let t = T::sin(1).unwrap();
        let s2 = s.clone() * s;
        let t2 = t.clone() * t;
        assert_eq!(
            (a.clone() * b.clone()) * (a * b),
            T::one() - s2.clone() - t2.clone() + s2 * t2
        );
    }
    #[test]
    fn normal_form_cancels_terms_and_preserves_independent_angles() {
        let a = T::sin(0).unwrap();
        let b = T::sin(1).unwrap();
        assert_eq!(a.clone() - a.clone(), T::zero());
        assert_ne!(a, b);
        assert_eq!((T::from_i64(3) * a.clone()) - a.clone(), T::from_i64(2) * a);
    }
    #[test]
    #[should_panic(expected = "MultiTrig degree budget exceeded")]
    fn expression_budget_is_enforced() {
        let mut a = T::sin(0).unwrap();
        for _ in 0..7 {
            a = a.clone() * a;
        }
    }
}
