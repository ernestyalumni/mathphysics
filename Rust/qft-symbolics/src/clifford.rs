//! Ordered Clifford algebra with {gamma(a),gamma(b)}=2h(a,b), h=-eta.
//! Labels are unsummed slots. No gamma5 prescription. Legacy wrappers panic.
use crate::{
    coefficient::Coefficient,
    computation::{Budget, Error, Resource, Result},
    polynomial::{BigPoly, Poly, Polynomial},
};
use std::collections::BTreeMap;
pub type NormalWord = BTreeMap<Vec<String>, Poly>;
pub type GenericNormalWord<A> = BTreeMap<Vec<String>, Polynomial<A>>;
fn labels_checked(labels: &[&str], cap: usize, b: &mut Budget) -> Result<()> {
    b.charge(1)?;
    if labels.len() > cap {
        return Err(Error::ResourceLimit {
            resource: Resource::WordLength,
            limit: cap as u64,
            required: labels.len() as u64,
        });
    }
    b.check(Resource::WordLength, labels.len() as u64)?;
    let bytes = labels
        .iter()
        .fold(0u64, |n, s| n.saturating_add(s.len() as u64));
    b.check(Resource::InputBytes, bytes)?;
    b.charge(labels.len() as u64)?;
    Ok(())
}
pub fn try_metric_with<A: Coefficient>(a: &str, c: &str, b: &mut Budget) -> Result<Polynomial<A>> {
    labels_checked(&[a, c], 12, b)?;
    let (a, c) = if a <= c { (a, c) } else { (c, a) };
    // Include encoding overhead in the preflight, before format allocates.
    let digits = |n: usize| if n == 0 { 1 } else { n.ilog10() as usize + 1 };
    let bytes = a
        .len()
        .saturating_add(c.len())
        .saturating_add(digits(a.len()))
        .saturating_add(digits(c.len()))
        .saturating_add(6);
    Polynomial::try_var_built(
        bytes,
        || format!("h({}:{};{}:{})", a.len(), a, c.len(), c),
        b,
    )
}
pub fn metric(a: &str, b: &str) -> Poly {
    try_metric_with(a, b, &mut crate::polynomial::legacy_budget()).expect("Clifford metric")
}
fn trace_inner<A: Coefficient>(labels: &[&str], b: &mut Budget) -> Result<Polynomial<A>> {
    b.charge(1)?;
    if labels.is_empty() {
        return Polynomial::try_constant(A::from_i64(4), b);
    }
    if labels.len() % 2 == 1 {
        return Ok(Polynomial::default());
    }
    let mut out = Polynomial::default();
    for j in 1..labels.len() {
        b.charge(1)?;
        let h = try_metric_with::<A>(labels[0], labels[j], b)?;
        let rest: Vec<_> = labels[1..j]
            .iter()
            .chain(&labels[j + 1..])
            .copied()
            .collect();
        let t = h.try_mul(&trace_inner::<A>(&rest, b)?, b)?;
        out = if j % 2 == 1 {
            out.try_add(&t, b)?
        } else {
            out.try_sub(&t, b)?
        };
    }
    Ok(out)
}
/// Four-component Dirac trace, without gamma5; hard algorithm cap 12.
pub fn try_trace_with<A: Coefficient>(labels: &[&str], b: &mut Budget) -> Result<Polynomial<A>> {
    labels_checked(labels, 12, b)?;
    trace_inner(labels, b)
}
pub fn try_trace(labels: &[&str], b: &mut Budget) -> Result<Poly> {
    try_trace_with(labels, b)
}
pub fn try_trace_exact(labels: &[&str], b: &mut Budget) -> Result<BigPoly> {
    try_trace_with(labels, b)
}
pub fn trace(labels: &[&str]) -> Poly {
    try_trace(labels, &mut crate::polynomial::legacy_budget()).expect("Clifford trace")
}
fn add<A: Coefficient>(
    out: &mut GenericNormalWord<A>,
    w: Vec<String>,
    p: Polynomial<A>,
    b: &mut Budget,
) -> Result<()> {
    b.charge(1)?;
    let p = if let Some(old) = out.get(&w) {
        old.try_add(&p, b)?
    } else {
        p
    };
    if p.is_zero() {
        out.remove(&w);
    } else {
        b.check(
            Resource::Terms,
            (out.len() + usize::from(!out.contains_key(&w))) as u64,
        )?;
        out.insert(w, p);
    }
    Ok(())
}
fn reduce<A: Coefficient>(
    w: Vec<String>,
    p: Polynomial<A>,
    out: &mut GenericNormalWord<A>,
    b: &mut Budget,
) -> Result<()> {
    b.charge(1)?;
    for i in 0..w.len().saturating_sub(1) {
        if w[i] >= w[i + 1] {
            let h = try_metric_with::<A>(&w[i], &w[i + 1], b)?;
            let mut contracted = w.clone();
            contracted.drain(i..i + 2);
            if w[i] == w[i + 1] {
                reduce(contracted, p.try_mul(&h, b)?, out, b)?;
            } else {
                let mut swapped = w;
                swapped.swap(i, i + 1);
                reduce(swapped, p.try_neg(b)?, out, b)?;
                let twice = Polynomial::try_constant(A::from_i64(2), b)?;
                reduce(contracted, p.try_mul(&twice, b)?.try_mul(&h, b)?, out, b)?;
            }
            return Ok(());
        }
    }
    add(out, w, p, b)
}
pub fn try_normal_order_with<A: Coefficient>(
    labels: &[&str],
    b: &mut Budget,
) -> Result<GenericNormalWord<A>> {
    labels_checked(labels, 12, b)?;
    let p = Polynomial::try_constant(A::from_i64(1), b)?;
    let mut out = BTreeMap::new();
    reduce(
        labels.iter().map(|s| s.to_string()).collect(),
        p,
        &mut out,
        b,
    )?;
    Ok(out)
}
pub fn try_normal_order(labels: &[&str], b: &mut Budget) -> Result<NormalWord> {
    try_normal_order_with(labels, b)
}
pub fn normal_order(labels: &[&str]) -> NormalWord {
    try_normal_order(labels, &mut crate::polynomial::legacy_budget())
        .expect("Clifford normal ordering")
}
/// Explicit contracted gamma^mu W gamma_mu, gamma^mu gamma_mu=-d.
/// Formal dimension, ordinary labels unsummed; no gamma5 prescription.
pub fn try_sandwich_srednicki_with<A: Coefficient>(
    labels: &[&str],
    dimension: &Polynomial<A>,
    b: &mut Budget,
) -> Result<GenericNormalWord<A>> {
    labels_checked(labels, 8, b)?;
    let sign = Polynomial::try_constant(
        A::from_i64(if labels.len().is_multiple_of(2) {
            -1
        } else {
            1
        }),
        b,
    )?;
    let leading = sign.try_mul(dimension, b)?;
    let mut result = BTreeMap::new();
    for (w, p) in try_normal_order_with::<A>(labels, b)? {
        add(&mut result, w, leading.try_mul(&p, b)?, b)?;
    }
    for j in 0..labels.len() {
        b.charge(1)?;
        let mut rest: Vec<_> = labels
            .iter()
            .enumerate()
            .filter_map(|(i, s)| if i == j { None } else { Some(*s) })
            .collect();
        rest.push(labels[j]);
        let c = Polynomial::try_constant(A::from_i64(if j % 2 == 0 { -2 } else { 2 }), b)?;
        for (w, p) in try_normal_order_with::<A>(&rest, b)? {
            add(&mut result, w, c.try_mul(&p, b)?, b)?;
        }
    }
    Ok(result)
}
pub fn try_sandwich_srednicki(
    labels: &[&str],
    dimension: &Poly,
    b: &mut Budget,
) -> Result<NormalWord> {
    try_sandwich_srednicki_with(labels, dimension, b)
}
pub fn sandwich_srednicki(labels: &[&str], dimension: &Poly) -> NormalWord {
    try_sandwich_srednicki(labels, dimension, &mut crate::polynomial::legacy_budget())
        .expect("Clifford sandwich")
}
