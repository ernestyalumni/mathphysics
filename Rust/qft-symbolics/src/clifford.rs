//! Ordered gamma algebra with {gamma(a),gamma(b)} = 2 h(a,b).
//! h is a symbolic symmetric Clifford form. For Srednicki h=-eta.
//! Labels denote slots/components, NEVER implicitly summed dummy indices.
//! Canonicalisation uses oriented swaps, not general tensor canonicalisation.
use crate::polynomial::Poly;
use std::collections::BTreeMap;
pub fn metric(a: &str, b: &str) -> Poly {
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    // Length prefixes prevent arbitrary input labels from colliding.
    Poly::var(format!("h({}:{};{}:{})", a.len(), a, b.len(), b))
}
pub type NormalWord = BTreeMap<Vec<String>, Poly>;
fn add(out: &mut NormalWord, w: Vec<String>, p: Poly) {
    let p = out.remove(&w).unwrap_or_default() + p;
    if !p.is_zero() {
        out.insert(w, p);
    }
}
fn reduce(w: Vec<String>, p: Poly, out: &mut NormalWord, budget: &mut usize) {
    assert!(*budget > 0, "Clifford expansion budget");
    *budget -= 1;
    for i in 0..w.len().saturating_sub(1) {
        if w[i] >= w[i + 1] {
            let h = metric(&w[i], &w[i + 1]);
            let mut contracted = w.clone();
            contracted.drain(i..i + 2);
            if w[i] == w[i + 1] {
                reduce(contracted, p * h, out, budget);
            } else {
                let mut swapped = w;
                swapped.swap(i, i + 1);
                reduce(swapped, -p.clone(), out, budget);
                reduce(contracted, p * Poly::constant(2) * h, out, budget);
            }
            return;
        }
    }
    add(out, w, p);
}
pub fn normal_order(labels: &[&str]) -> NormalWord {
    assert!(labels.len() <= 12, "gamma word length budget");
    let mut out = NormalWord::new();
    reduce(
        labels.iter().map(|s| s.to_string()).collect(),
        Poly::constant(1),
        &mut out,
        &mut 100_000,
    );
    out
}
/// Four-component Dirac trace without gamma5, exact in symbolic h.
/// Oriented pairing recurrence, independent of normal_order's swap algorithm.
pub fn trace(labels: &[&str]) -> Poly {
    assert!(labels.len() <= 12, "gamma trace length budget");
    if labels.is_empty() {
        return Poly::constant(4);
    }
    if labels.len() % 2 == 1 {
        return Poly::default();
    }
    let mut out = Poly::default();
    for j in 1..labels.len() {
        let rest: Vec<_> = labels[1..j]
            .iter()
            .chain(&labels[j + 1..])
            .copied()
            .collect();
        let s = if j % 2 == 1 { 1 } else { -1 };
        out = out + Poly::constant(s) * metric(labels[0], labels[j]) * trace(&rest);
    }
    out
}
