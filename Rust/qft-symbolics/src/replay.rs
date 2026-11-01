//! Independent bounded algebra replay, sharing only the trusted coefficient core.
use crate::{
    computation::{Budget, Error, Resource, Result},
    exact::BigC,
    expression::Expr,
};
use std::collections::{BTreeMap, BTreeSet};
pub const CERTIFICATE_VERSION: u32 = 1;
pub const SCALAR_DOMAIN: &str = "gaussian-rational-commuting-scalars/explicit-nonzero-poles/v1";
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rule {
    PolynomialIdentity,
    RationalIdentity,
}
#[derive(Clone, Debug)]
pub struct ReplayStep {
    pub before: Expr,
    pub after: Expr,
    pub rule: Rule,
}
#[derive(Clone, Debug)]
pub struct Certificate {
    pub version: u32,
    pub domain: String,
    pub start: Expr,
    pub end: Expr,
    pub steps: Vec<ReplayStep>,
    pub nonzero_assumptions: Vec<Expr>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReplayReport {
    pub steps_checked: usize,
    pub conditional_domain: Vec<String>,
    pub statement: &'static str,
    pub formal_proof: bool,
}
fn bad(message: &'static str) -> Error {
    Error::UnsupportedDomain { message }
}
type Mono = Vec<String>;
type Poly = BTreeMap<Mono, BigC>;
struct Fraction {
    n: Poly,
    d: Poly,
}
fn scan(e: &Expr, depth: u64, nodes: &mut u64, bytes: &mut u64, b: &mut Budget) -> Result<()> {
    if depth > 128 {
        return Err(bad("replay AST hard depth cap 128"));
    }
    b.check(Resource::Degree, depth)?;
    *nodes = nodes.saturating_add(1);
    b.check(Resource::Terms, *nodes)?;
    b.charge(1)?;
    let extra = match e {
        Expr::Constant(c) => {
            b.check(Resource::CoefficientBits, c.bits())?;
            64u64.saturating_add(c.bits().saturating_mul(4))
        }
        Expr::Symbol(s) => 64u64.saturating_add(s.len() as u64),
        Expr::Sum(v) | Expr::Product(v) => {
            b.check(Resource::Terms, v.len() as u64)?;
            64
        }
        Expr::Quotient(_, _) => 64,
        _ => return Err(bad("replay requires commuting scalar AST only")),
    };
    *bytes = bytes.saturating_add(extra);
    b.check(Resource::InputBytes, *bytes)?;
    match e {
        Expr::Sum(v) | Expr::Product(v) => {
            for x in v {
                scan(x, depth + 1, nodes, bytes, b)?;
            }
        }
        Expr::Quotient(x, y) => {
            scan(x, depth + 1, nodes, bytes, b)?;
            scan(y, depth + 1, nodes, bytes, b)?;
        }
        _ => (),
    }
    Ok(())
}
fn key(e: &Expr, b: &mut Budget) -> Result<String> {
    e.validate(b)?;
    e.canonical(b)
}
fn mb(m: &Mono) -> u64 {
    m.iter().fold(64u64, |n, s| {
        n.saturating_add(32).saturating_add(s.len() as u64)
    })
}
fn size(p: &Poly) -> u64 {
    p.iter().fold(0u64, |n, (m, c)| {
        n.saturating_add(mb(m))
            .saturating_add(c.bits().saturating_mul(4))
            .saturating_add(64)
    })
}
fn inspect(p: &Poly, b: &mut Budget) -> Result<()> {
    b.charge(p.len() as u64)?;
    b.check(Resource::Terms, p.len() as u64)?;
    b.check(Resource::InputBytes, size(p))?;
    for (m, c) in p {
        b.check(Resource::Degree, m.len() as u64)?;
        b.check(Resource::CoefficientBits, c.bits())?;
    }
    b.charge(size(p).max(1))
}
fn constant(c: &BigC, b: &mut Budget) -> Result<Poly> {
    b.check(Resource::CoefficientBits, c.bits())?;
    b.check(
        Resource::InputBytes,
        128u64.saturating_add(c.bits().saturating_mul(4)),
    )?;
    b.charge(c.bits().max(1))?;
    let mut p = Poly::new();
    if !c.is_zero() {
        b.check(Resource::Terms, 1)?;
        p.insert(vec![], c.clone());
    }
    Ok(p)
}
fn one(b: &mut Budget) -> Result<Poly> {
    constant(&BigC::from(1), b)
}
fn insert(p: &mut Poly, m: Mono, c: BigC, b: &mut Budget) -> Result<()> {
    b.check(Resource::Degree, m.len() as u64)?;
    if c.is_zero() {
        return Ok(());
    }
    let v = match p.get(&m) {
        Some(old) => old.checked_add(&c, b)?,
        None => c,
    };
    b.check(Resource::Terms, (p.len() as u64).saturating_add(1))?;
    b.charge(p.len() as u64)?;
    b.check(
        Resource::InputBytes,
        size(p)
            .saturating_add(mb(&m))
            .saturating_add(v.bits().saturating_mul(4))
            .saturating_add(64),
    )?;
    b.charge(mb(&m).saturating_mul((p.len() as u64).saturating_add(1).ilog2() as u64 + 1))?;
    if v.is_zero() {
        p.remove(&m);
    } else {
        p.insert(m, v);
    }
    Ok(())
}
fn add(mut a: Poly, c: Poly, b: &mut Budget) -> Result<Poly> {
    b.charge((a.len() as u64).saturating_add(c.len() as u64))?;
    b.check(
        Resource::Terms,
        (a.len() as u64).saturating_add(c.len() as u64),
    )?;
    b.check(Resource::InputBytes, size(&a).saturating_add(size(&c)))?;
    for (m, v) in c {
        insert(&mut a, m, v, b)?;
    }
    Ok(a)
}
fn mul(a: &Poly, c: &Poly, b: &mut Budget) -> Result<Poly> {
    inspect(a, b)?;
    inspect(c, b)?;
    let pairs = (a.len() as u64).saturating_mul(c.len() as u64);
    b.check(Resource::Terms, pairs)?;
    b.charge(pairs)?;
    // Reject an unaffordable dense coefficient product before entering its loops.
    let ab = a.values().map(BigC::bits).max().unwrap_or(0);
    let cb = c.values().map(BigC::bits).max().unwrap_or(0);
    let bits = ab.saturating_add(cb).saturating_add(2);
    let work = pairs.saturating_mul(bits.saturating_mul(bits).saturating_mul(128));
    b.check(Resource::Operations, b.used().saturating_add(work))?;

    let ad = a.keys().map(Vec::len).max().unwrap_or(0) as u64;
    let cd = c.keys().map(Vec::len).max().unwrap_or(0) as u64;
    b.check(Resource::Degree, ad.saturating_add(cd))?;
    let largest = |p: &Poly| {
        p.iter()
            .map(|(m, v)| {
                mb(m)
                    .saturating_add(v.bits().saturating_mul(16))
                    .saturating_add(128)
            })
            .max()
            .unwrap_or(0)
    };
    b.check(
        Resource::InputBytes,
        pairs.saturating_mul(largest(a).saturating_add(largest(c))),
    )?;
    let mut out = Poly::new();
    for (m, x) in a {
        for (n, y) in c {
            let coefficient = x.checked_mul(y, b)?;
            b.charge(
                mb(m)
                    .saturating_add(mb(n))
                    .saturating_mul(ad.saturating_add(cd).saturating_add(1)),
            )?;
            let mut monomial = Vec::with_capacity(m.len() + n.len());
            monomial.extend(m.iter().cloned());
            monomial.extend(n.iter().cloned());
            monomial.sort();
            insert(&mut out, monomial, coefficient, b)?;
        }
    }
    Ok(out)
}
fn equal(a: &Poly, c: &Poly, b: &mut Budget) -> Result<bool> {
    inspect(a, b)?;
    inspect(c, b)?;
    Ok(a == c)
}
fn expand(
    e: &Expr,
    rational: bool,
    assumptions: &BTreeSet<String>,
    b: &mut Budget,
) -> Result<Fraction> {
    b.charge(1)?;
    match e {
        Expr::Constant(c) => Ok(Fraction {
            n: constant(c, b)?,
            d: one(b)?,
        }),
        Expr::Symbol(s) => {
            b.check(Resource::Terms, 1)?;
            b.check(Resource::Degree, 1)?;
            b.check(Resource::InputBytes, 256u64.saturating_add(s.len() as u64))?;
            b.charge(s.len() as u64 + 1)?;
            let mut n = Poly::new();
            n.insert(vec![s.clone()], BigC::from(1));
            Ok(Fraction { n, d: one(b)? })
        }
        Expr::Sum(v) | Expr::Product(v) => {
            let product = matches!(e, Expr::Product(_));
            let mut out = Fraction {
                n: constant(&BigC::from(if product { 1 } else { 0 }), b)?,
                d: one(b)?,
            };
            for x in v {
                let rhs = expand(x, rational, assumptions, b)?;
                let n = if product {
                    mul(&out.n, &rhs.n, b)?
                } else {
                    add(mul(&out.n, &rhs.d, b)?, mul(&rhs.n, &out.d, b)?, b)?
                };
                let d = mul(&out.d, &rhs.d, b)?;
                out = Fraction { n, d };
            }
            Ok(out)
        }
        Expr::Quotient(x, y) => {
            if !rational {
                return Err(bad("polynomial rule forbids quotients"));
            }
            // No short circuit at zero: every nested pole survives as an obligation.
            let a = expand(x, true, assumptions, b)?;
            let d = expand(y, true, assumptions, b)?;
            if d.n.is_empty() {
                return Err(Error::DivisionByZero);
            }
            let predicate = key(y, b)?;
            b.charge((predicate.len() as u64).saturating_mul(assumptions.len() as u64 + 1))?;
            if !assumptions.contains(&predicate) {
                return Err(bad("missing explicit nonzero denominator assumption"));
            }
            Ok(Fraction {
                n: mul(&a.n, &d.d, b)?,
                d: mul(&a.d, &d.n, b)?,
            })
        }
        _ => Err(bad("unsupported replay AST")),
    }
}
/// Replays exact identities on a conditional domain; emits no formal proof.
pub fn replay(c: &Certificate, b: &mut Budget) -> Result<ReplayReport> {
    b.check(Resource::InputBytes, c.domain.len() as u64)?;
    b.charge(c.domain.len() as u64 + 1)?;
    if c.version != CERTIFICATE_VERSION || c.domain != SCALAR_DOMAIN {
        return Err(bad("certificate version/domain convention mismatch"));
    }
    if c.steps.is_empty() {
        return Err(bad("empty replay certificate"));
    }
    b.check(Resource::Terms, c.steps.len() as u64)?;
    b.check(Resource::Terms, c.nonzero_assumptions.len() as u64)?;
    let mut nodes = 0;
    let mut bytes = c.domain.len() as u64;
    for e in std::iter::once(&c.start)
        .chain(std::iter::once(&c.end))
        .chain(c.nonzero_assumptions.iter())
        .chain(c.steps.iter().flat_map(|s| [&s.before, &s.after]))
    {
        scan(e, 1, &mut nodes, &mut bytes, b)?;
        e.validate(b)?;
    }
    let mut assumptions = BTreeSet::new();
    let mut key_bytes = 0u64;
    for e in &c.nonzero_assumptions {
        let k = key(e, b)?;
        key_bytes = key_bytes.saturating_add(k.len() as u64);
        b.check(Resource::InputBytes, key_bytes)?;
        b.charge(key_bytes.max(1))?;
        assumptions.insert(k);
    }
    for e in &c.nonzero_assumptions {
        if expand(e, true, &assumptions, b)?.n.is_empty() {
            return Err(bad("nonzero assumption normalizes to zero"));
        }
    }
    expand(&c.start, true, &assumptions, b)?;
    expand(&c.end, true, &assumptions, b)?;
    let mut previous = key(&c.start, b)?;
    for step in &c.steps {
        let before = key(&step.before, b)?;
        b.charge(previous.len().saturating_add(before.len()) as u64)?;
        if previous != before {
            return Err(bad("noncontiguous replay chain"));
        }
        let rational = step.rule == Rule::RationalIdentity;
        let a = expand(&step.before, rational, &assumptions, b)?;
        let d = expand(&step.after, rational, &assumptions, b)?;
        let valid = if rational {
            equal(&mul(&a.n, &d.d, b)?, &mul(&d.n, &a.d, b)?, b)?
        } else {
            equal(&a.n, &d.n, b)?
        };
        if !valid {
            return Err(bad("replay algebraic identity failed"));
        }
        previous = key(&step.after, b)?;
    }
    let end = key(&c.end, b)?;
    b.charge(previous.len().saturating_add(end.len()) as u64)?;
    if previous != end {
        return Err(bad("replay end binding failed"));
    }
    Ok(ReplayReport{steps_checked:c.steps.len(),conditional_domain:assumptions.into_iter().map(|s|format!("nonzero({s})")).collect(),statement:"Independently replayed bounded algebra over Gaussian-rational commuting scalars, conditional on supplied nonzero assumptions. Not a Lean/formal QFT proof.",formal_proof:false})
}
