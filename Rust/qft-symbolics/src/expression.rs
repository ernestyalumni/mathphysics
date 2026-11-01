//! Bounded structural AST. No parser, implicit contraction or tensor symmetry reduction.
pub use crate::tensor::{Space, Variance};
use crate::{
    computation::{Budget, Error, Resource, Result},
    exact::BigC,
};
use std::collections::BTreeMap;
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Expr {
    Constant(BigC),
    Symbol(String),
    Sum(Vec<Expr>),
    Product(Vec<Expr>),
    OrderedProduct(Vec<Expr>),
    Quotient(Box<Expr>, Box<Expr>),
    Operator(String),
    Tensor(TensorAtom),
    Contract {
        binder: String,
        space: Space,
        body: Box<Expr>,
    },
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum IndexName {
    Free(String),
    Bound(String),
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Index {
    pub name: IndexName,
    pub space: Space,
    pub variance: Variance,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TensorAtom {
    pub name: String,
    pub slots: Vec<(Space, Variance)>,
    pub indices: Vec<Index>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Domain {
    Scalar,
    Tensor,
    Operator,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Signature {
    pub domain: Domain,
    pub free_indices: Vec<Index>,
}
fn invalid(message: &'static str) -> Error {
    Error::UnsupportedDomain { message }
}
fn name(s: &str) -> String {
    format!("{}:{}", s.len(), s)
}
impl Expr {
    fn preflight(
        &self,
        b: &mut Budget,
        depth: u64,
        nodes: &mut u64,
        bytes: &mut u64,
    ) -> Result<()> {
        b.check(Resource::Degree, depth)?;
        if depth > 256 {
            return Err(invalid("AST depth exceeds hard cap 256"));
        }
        *nodes = nodes.saturating_add(1);
        b.check(Resource::Terms, *nodes)?;
        b.charge(1)?;
        let mut add = 64u64;
        match self {
            Self::Constant(c) => {
                b.check(Resource::CoefficientBits, c.bits())?;
                add = add.saturating_add(c.bits().saturating_mul(4));
            }
            Self::Symbol(s) | Self::Operator(s) => {
                add = add.saturating_add(s.len() as u64);
            }
            Self::Tensor(t) => {
                add = add.saturating_add(t.name.len() as u64);
                b.check(Resource::Terms, t.indices.len() as u64)?;
                b.check(Resource::Terms, t.slots.len() as u64)?;
                for i in &t.indices {
                    let (IndexName::Free(s) | IndexName::Bound(s)) = &i.name;
                    add = add.saturating_add(s.len() as u64).saturating_add(64);
                }
            }
            Self::Contract { binder, body, .. } => {
                add = add.saturating_add(binder.len() as u64);
                body.preflight(b, depth + 1, nodes, bytes)?;
            }
            Self::Quotient(x, y) => {
                x.preflight(b, depth + 1, nodes, bytes)?;
                y.preflight(b, depth + 1, nodes, bytes)?;
            }
            Self::Sum(v) | Self::Product(v) | Self::OrderedProduct(v) => {
                b.check(Resource::Terms, v.len() as u64)?;
                for x in v {
                    x.preflight(b, depth + 1, nodes, bytes)?;
                }
            }
        }
        *bytes = bytes.saturating_add(add);
        b.check(Resource::InputBytes, *bytes)?;
        Ok(())
    }
    pub fn validate(&self, b: &mut Budget) -> Result<Signature> {
        let (mut nodes, mut bytes) = (0, 0);
        self.preflight(b, 1, &mut nodes, &mut bytes)?;
        // Conservative bound for signature copies, name comparisons and sorting.
        b.charge(
            bytes
                .saturating_mul(nodes)
                .saturating_mul(bytes.max(1).ilog2() as u64 + 1),
        )?;
        let s = self.typed(&mut vec![], b)?;
        if s.free_indices
            .iter()
            .any(|i| matches!(i.name, IndexName::Bound(_)))
        {
            return Err(invalid("unbound index"));
        }
        Ok(s)
    }
    fn typed(&self, env: &mut Vec<(String, Space)>, b: &mut Budget) -> Result<Signature> {
        b.charge(1)?;
        let scalar = || Signature {
            domain: Domain::Scalar,
            free_indices: vec![],
        };
        match self {
            Self::Constant(_) => Ok(scalar()),
            Self::Symbol(s) | Self::Operator(s) => {
                if s.is_empty() {
                    return Err(invalid("empty symbol"));
                }
                Ok(Signature {
                    domain: if matches!(self, Self::Operator(_)) {
                        Domain::Operator
                    } else {
                        Domain::Scalar
                    },
                    free_indices: vec![],
                })
            }
            Self::Tensor(t) => {
                if t.name.is_empty() || t.slots.len() != t.indices.len() {
                    return Err(invalid("malformed tensor atom"));
                }
                for (i, slot) in t.indices.iter().zip(&t.slots) {
                    if (i.space, i.variance) != *slot {
                        return Err(invalid("tensor slot space/variance mismatch"));
                    }
                    match &i.name {
                        IndexName::Free(n) if !n.is_empty() => {}
                        IndexName::Bound(n)
                            if !n.is_empty()
                                && env.iter().any(|(s, sp)| s == n && *sp == i.space) => {}
                        _ => return Err(invalid("empty or dangling index")),
                    }
                }
                let mut indices = t.indices.clone();
                indices.sort();
                check_free(&indices)?;
                Ok(Signature {
                    domain: Domain::Tensor,
                    free_indices: indices,
                })
            }
            Self::Sum(v) => {
                let mut sig = None;
                for x in v {
                    let s = x.typed(env, b)?;
                    if sig.as_ref().is_some_and(|a| a != &s) {
                        return Err(invalid("sum signature/domain mismatch"));
                    }
                    sig = Some(s);
                }
                Ok(sig.unwrap_or_else(scalar))
            }
            Self::Product(v) | Self::OrderedProduct(v) => {
                let ordered = matches!(self, Self::OrderedProduct(_));
                let mut out = scalar();
                for x in v {
                    let s = x.typed(env, b)?;
                    if !ordered && s.domain != Domain::Scalar {
                        return Err(invalid("commuting Product requires scalars"));
                    }
                    if s.domain == Domain::Operator {
                        if out.domain == Domain::Tensor {
                            return Err(invalid("mixed tensor/operator product"));
                        }
                        out.domain = Domain::Operator;
                    } else if s.domain == Domain::Tensor {
                        if out.domain == Domain::Operator {
                            return Err(invalid("mixed tensor/operator product"));
                        }
                        out.domain = Domain::Tensor;
                    }
                    if s.domain == Domain::Operator && !out.free_indices.is_empty() {
                        return Err(invalid("mixed tensor/operator product"));
                    }
                    out.free_indices.extend(s.free_indices);
                }
                out.free_indices.sort();
                check_free(&out.free_indices)?;
                Ok(out)
            }
            Self::Quotient(x, y) => {
                let a = x.typed(env, b)?;
                let d = y.typed(env, b)?;
                if d.domain != Domain::Scalar {
                    return Err(invalid("non-scalar denominator"));
                }
                if matches!(y.as_ref(),Self::Constant(c) if c.is_zero()) {
                    return Err(Error::DivisionByZero);
                }
                Ok(a)
            }
            Self::Contract {
                binder,
                space,
                body,
            } => {
                if binder.is_empty() || env.iter().any(|(n, _)| n == binder) {
                    return Err(invalid("empty or shadowed binder"));
                }
                env.push((binder.clone(), *space));
                let result = body.typed(env, b);
                env.pop();
                let mut s = result?;
                let hits: Vec<_> = s
                    .free_indices
                    .iter()
                    .filter(|i| matches!(&i.name,IndexName::Bound(n) if n==binder))
                    .collect();
                if hits.len() != 2
                    || hits[0].space != *space
                    || hits[1].space != *space
                    || hits[0].variance == hits[1].variance
                {
                    return Err(invalid(
                        "binder requires two opposite-variance occurrences per branch",
                    ));
                }
                s.free_indices
                    .retain(|i| !matches!(&i.name,IndexName::Bound(n) if n==binder));
                if s.free_indices.is_empty() && s.domain == Domain::Tensor {
                    s.domain = Domain::Scalar;
                }
                Ok(s)
            }
        }
    }
    pub fn canonical(&self, b: &mut Budget) -> Result<String> {
        self.validate(b)?;
        let rendered = self.render(&mut vec![], b)?;
        b.check(
            Resource::InputBytes,
            (rendered.len() as u64).saturating_add(12),
        )?;
        b.charge(12)?;
        Ok(format!("qft-expr.v1:{rendered}"))
    }
    /// Canonical nonzero predicates, including those beneath zero numerators/products.
    pub fn nonzero_obligations(&self, b: &mut Budget) -> Result<Vec<String>> {
        let (mut nodes, mut bytes) = (0, 0);
        self.preflight(b, 1, &mut nodes, &mut bytes)?;
        // A denominator can appear in several ancestor obligations.
        b.check(Resource::InputBytes, bytes.saturating_mul(nodes))?;
        self.validate(b)?;
        let mut out = vec![];
        self.obligations(&mut vec![], b, &mut out)?;
        out.sort();
        out.dedup();
        Ok(out)
    }
    fn obligations(
        &self,
        env: &mut Vec<String>,
        b: &mut Budget,
        out: &mut Vec<String>,
    ) -> Result<()> {
        b.charge(1)?;
        match self {
            Self::Quotient(x, y) => {
                out.push(format!("nonzero(qft-expr.v1:{})", y.render(env, b)?));
                x.obligations(env, b, out)?;
                y.obligations(env, b, out)?;
            }
            Self::Sum(v) | Self::Product(v) | Self::OrderedProduct(v) => {
                for x in v {
                    x.obligations(env, b, out)?;
                }
            }
            Self::Contract { binder, body, .. } => {
                env.push(binder.clone());
                body.obligations(env, b, out)?;
                env.pop();
            }
            _ => {}
        }
        Ok(())
    }
    fn render(&self, env: &mut Vec<String>, b: &mut Budget) -> Result<String> {
        b.charge(1)?;
        let s = match self {
            Self::Constant(c) => format!("C({},{})", c.re, c.im),
            Self::Symbol(s) => format!("S({})", name(s)),
            Self::Operator(s) => format!("O({})", name(s)),
            Self::Tensor(t) => {
                let mut s = format!("T({};", name(&t.name));
                for i in &t.indices {
                    let n = match &i.name {
                        IndexName::Free(n) => format!("f{}", name(n)),
                        IndexName::Bound(n) => format!(
                            "b{}",
                            env.iter()
                                .rev()
                                .position(|v| v == n)
                                .ok_or(invalid("dangling binder"))?
                        ),
                    };
                    s.push_str(&format!("{:?}:{:?}:{};", i.space, i.variance, n));
                }
                s.push(')');
                s
            }
            Self::Quotient(x, y) => format!("Q({},{})", x.render(env, b)?, y.render(env, b)?),
            Self::Contract {
                binder,
                space,
                body,
            } => {
                env.push(binder.clone());
                let s = body.render(env, b)?;
                env.pop();
                format!("B({space:?},{s})")
            }
            Self::Sum(v) | Self::Product(v) | Self::OrderedProduct(v) => {
                let mut parts = v
                    .iter()
                    .map(|x| x.render(env, b))
                    .collect::<Result<Vec<_>>>()?;
                if !matches!(self, Self::OrderedProduct(_)) {
                    let bytes = parts.iter().map(|s| s.len() as u64).sum::<u64>();
                    b.charge(bytes.saturating_mul((parts.len() as u64).max(1).ilog2() as u64 + 1))?;
                    parts.sort();
                }
                format!(
                    "{}({})",
                    if matches!(self, Self::Sum(_)) {
                        "A"
                    } else if matches!(self, Self::Product(_)) {
                        "M"
                    } else {
                        "P"
                    },
                    parts.join(",")
                )
            }
        };
        b.check(Resource::InputBytes, s.len() as u64)?;
        b.charge(s.len() as u64)?;
        Ok(s)
    }
}
fn check_free(indices: &[Index]) -> Result<()> {
    let mut seen = BTreeMap::new();
    for i in indices {
        if let IndexName::Free(n) = &i.name {
            if seen.insert(n, (i.space, i.variance)).is_some() {
                return Err(invalid("repeated free name requires explicit binder"));
            }
        }
    }
    Ok(())
}
