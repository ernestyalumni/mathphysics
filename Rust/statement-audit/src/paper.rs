//! The informal statements of `formal/informal-statements.json` as executable predicates.
//!
//! Written by hand from the informal text, with ordinary (partial) arithmetic: a division by zero,
//! the log of a nonpositive number, or `Z(0)` is undecided rather than given a conventional value.
//! A spectrum `λ₁, λ₂, …` is stored 0-based: element `i` is `λ_{i+1}`.
use crate::lean::FORALL_NAT;
use crate::num::{self, Iv};
use crate::rng::Rng;
use crate::value::{self as v, arith, nat, real, real_int, Mat, Op, Sem, Tri, V};
use std::collections::BTreeMap;
use std::rc::Rc;

pub type Inst = BTreeMap<String, V>;
type Pred = fn(&Inst) -> Tri;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Real,
    Nat,
    /// `λ₁, λ₂, …` (informal indices start at 1).
    Seq1,
    /// `δ₀, δ₁, …` (informal indices start at 0).
    Seq0,
    /// `N × N` matrix; `N` is the variable named `N`.
    Mat { complex: bool },
}

pub struct PVar {
    pub name: &'static str,
    pub aliases: &'static [&'static str],
    pub kind: Kind,
}

pub struct Paper {
    pub name: &'static str,
    pub vars: Vec<PVar>,
    pub hyps: Vec<(&'static str, Pred)>,
    pub concl: Vec<(&'static str, Pred)>,
    pub sample: fn(&mut Rng) -> Inst,
}

const M: Sem = Sem::Math;

// ---------- helpers ----------

fn seq(i: &Inst, k: &str) -> Vec<V> {
    match i.get(k) {
        Some(V::Seq(s)) => s.as_ref().clone(),
        _ => vec![],
    }
}
fn get(i: &Inst, k: &str) -> V {
    i.get(k).cloned().unwrap_or(V::Undecided)
}
fn at(s: &[V], j: usize) -> V {
    s.get(j).cloned().unwrap_or_else(|| real_int(0))
}
fn add(a: &V, b: &V) -> V {
    arith(Op::Add, a, b, M).unwrap_or(V::Undecided)
}
fn sub(a: &V, b: &V) -> V {
    arith(Op::Sub, a, b, M).unwrap_or(V::Undecided)
}
fn mul(a: &V, b: &V) -> V {
    arith(Op::Mul, a, b, M).unwrap_or(V::Undecided)
}
fn div(a: &V, b: &V) -> V {
    arith(Op::Div, a, b, M).unwrap_or(V::Undecided)
}
fn powv(a: &V, k: usize) -> V {
    v::power(a, &nat(k as u64), M).unwrap_or(V::Undecided)
}
fn to_real(n: &V) -> V {
    match n {
        V::Nat(k) | V::Int(k) => real(num::qbig(k.clone(), num_bigint::BigInt::from(1))),
        other => other.clone(),
    }
}
fn nat_of(i: &Inst, k: &str) -> Option<usize> {
    get(i, k).as_usize()
}
fn all<F: Fn(usize) -> Tri>(n: usize, f: F) -> Tri {
    (0..n).fold(Tri::T, |acc, j| if acc == Tri::F { acc } else { acc.and(f(j)) })
}
fn ln(a: &V) -> V {
    v::log(a, M)
}

/// `Z(m) = Σ_j λ_j^m`; undecided for `m = 0` (the sum diverges).
fn partition(lam: &[V], m: usize) -> V {
    if m == 0 {
        return V::Undecided;
    }
    lam.iter().fold(real_int(0), |acc, x| add(&acc, &powv(x, m)))
}
fn purity(lam: &[V], n: usize) -> V {
    div(&partition(lam, 2 * n), &powv(&partition(lam, n), 2))
}
/// `Σ_{j ≥ 2} (λ_j / λ_1)^k`.
fn excited(lam: &[V], k: usize) -> V {
    let top = at(lam, 0);
    lam.iter().skip(1).fold(real_int(0), |acc, x| add(&acc, &powv(&div(x, &top), k)))
}
fn spectrum_hyps(i: &Inst, need_second: bool) -> Tri {
    let lam = seq(i, "lam");
    let top = at(&lam, 0);
    let mut t = v::lt(&real_int(0), &top);
    if need_second {
        t = t.and(v::lt(&real_int(0), &at(&lam, 1)));
    }
    // every λ_j, including the zeros beyond the support
    let bounds = all(lam.len() + 1, |j| v::le(&real_int(0), &at(&lam, j)).and(v::le(&at(&lam, j), &top)));
    t.and(bounds)
}

// ---------- the seven statements ----------

fn h_summable(_: &Inst) -> Tri {
    Tri::T // finitely supported by construction
}
fn h_top_pos(i: &Inst) -> Tri {
    v::lt(&real_int(0), &at(&seq(i, "lam"), 0))
}
fn h_second_pos(i: &Inst) -> Tri {
    v::lt(&real_int(0), &at(&seq(i, "lam"), 1))
}
fn h_bounds(i: &Inst) -> Tri {
    let lam = seq(i, "lam");
    let top = at(&lam, 0);
    all(lam.len() + 1, |j| v::le(&real_int(0), &at(&lam, j)).and(v::le(&at(&lam, j), &top)))
}
fn h_n_pos(i: &Inst) -> Tri {
    v::le(&nat(1), &get(i, "n"))
}

fn pb_parts(i: &Inst) -> Option<(V, V, V, V)> {
    let lam = seq(i, "lam");
    let n = nat_of(i, "n")?;
    Some((purity(&lam, n), excited(&lam, n), excited(&lam, 2 * n), lam.iter().skip(1).fold(real_int(0), |acc, x| {
        add(&acc, &powv(&powv(&div(x, &at(&lam, 0)), n), 2))
    })))
}
fn c_purity_formula(i: &Inst) -> Tri {
    let Some((p, s, _, q)) = pb_parts(i) else { return Tri::U };
    v::eq(&p, &div(&add(&real_int(1), &q), &powv(&add(&real_int(1), &s), 2)))
}
fn c_purity_le_one(i: &Inst) -> Tri {
    let Some((p, ..)) = pb_parts(i) else { return Tri::U };
    v::le(&p, &real_int(1))
}
fn c_lower(i: &Inst) -> Tri {
    let Some((p, s, ..)) = pb_parts(i) else { return Tri::U };
    v::le(&s, &sub(&div(&real_int(1), &p), &real_int(1)))
}
fn c_upper(i: &Inst) -> Tri {
    let Some((p, s, ..)) = pb_parts(i) else { return Tri::U };
    v::le(&sub(&div(&real_int(1), &p), &real_int(1)), &add(&mul(&real_int(2), &s), &powv(&s, 2)))
}
fn c_doubling(i: &Inst) -> Tri {
    let Some((_, s, s2, _)) = pb_parts(i) else { return Tri::U };
    v::le(&s2, &powv(&s, 2))
}

fn gap(i: &Inst) -> V {
    let lam = seq(i, "lam");
    ln(&div(&at(&lam, 0), &at(&lam, 1)))
}
fn c_exp_gap(i: &Inst) -> Tri {
    let lam = seq(i, "lam");
    let Some(n) = nat_of(i, "n") else { return Tri::U };
    let p = purity(&lam, n);
    let lhs = v::exp(&v::negate(&mul(&to_real(&nat(n as u64)), &gap(i)), M).unwrap_or(V::Undecided));
    v::le(&lhs, &sub(&div(&real_int(1), &p), &real_int(1)))
}
fn c_gap_bound(i: &Inst) -> Tri {
    let lam = seq(i, "lam");
    let Some(n) = nat_of(i, "n") else { return Tri::U };
    let p = purity(&lam, n);
    let premise = v::lt(&real(num::q(1, 2)), &p);
    if premise == Tri::F {
        return Tri::T;
    }
    let rhs = mul(
        &div(&real_int(1), &to_real(&nat(n as u64))),
        &ln(&div(&real_int(1), &sub(&div(&real_int(1), &p), &real_int(1)))),
    );
    premise.implies(v::le(&rhs, &gap(i)))
}

fn zs(i: &Inst) -> Vec<V> {
    ["z1", "z2", "z3", "z4", "z5"].iter().map(|k| get(i, k)).collect()
}
fn ps(z: &[V]) -> Vec<V> {
    (0..4).map(|k| div(&z[k + 1], &powv(&z[k], 2))).collect()
}
fn defect(z: &[V]) -> V {
    sub(&mul(&real_int(16), &ln(&z[0])), &ln(&z[4]))
}
fn h_z_pos(i: &Inst) -> Tri {
    zs(i).iter().fold(Tri::T, |acc, z| acc.and(v::lt(&real_int(0), z)))
}
fn h_p_le_one(i: &Inst) -> Tri {
    ps(&zs(i)).iter().fold(Tri::T, |acc, p| acc.and(v::le(p, &real_int(1))))
}
fn c_telescope(i: &Inst) -> Tri {
    let z = zs(i);
    let p = ps(&z);
    let rhs = mul(&mul(&mul(&powv(&p[0], 8), &powv(&p[1], 4)), &powv(&p[2], 2)), &p[3]);
    v::eq(&v::exp(&v::negate(&defect(&z), M).unwrap_or(V::Undecided)), &rhs)
}
fn c_defect_nonneg(i: &Inst) -> Tri {
    v::le(&real_int(0), &defect(&zs(i)))
}
fn c_last(i: &Inst) -> Tri {
    let z = zs(i);
    v::le(&v::exp(&v::negate(&defect(&z), M).unwrap_or(V::Undecided)), &ps(&z)[3])
}

fn delta_v(i: &Inst) -> V {
    get(i, "Delta")
}
fn h_delta_pos(i: &Inst) -> Tri {
    v::lt(&real_int(0), &delta_v(i))
}
fn h_delta_small(i: &Inst) -> Tri {
    v::lt(&delta_v(i), &ln(&real_int(2)))
}
fn h_purity_vs_defect(i: &Inst) -> Tri {
    let Some(n) = nat_of(i, "n") else { return Tri::U };
    let p = purity(&seq(i, "lam"), n);
    v::le(&v::exp(&v::negate(&delta_v(i), M).unwrap_or(V::Undecided)), &p)
}
fn c_gap_from_defect(i: &Inst) -> Tri {
    let Some(n) = nat_of(i, "n") else { return Tri::U };
    let rhs = mul(
        &div(&real_int(1), &to_real(&nat(n as u64))),
        &ln(&div(&real_int(1), &sub(&v::exp(&delta_v(i)), &real_int(1)))),
    );
    v::le(&rhs, &gap(i))
}

fn c1(i: &Inst) -> V {
    get(i, "c1")
}
fn threshold(i: &Inst) -> V {
    div(&real_int(1), &mul(&real_int(2), &c1(i)))
}
fn bound_k(i: &Inst, k: usize) -> V {
    let Some(n) = nat_of(i, "n") else { return V::Undecided };
    let d = at(&seq(i, "delta"), k);
    let scale = div(&real_int(1), &mul(&powv(&real_int(2), k), &to_real(&nat(n as u64))));
    mul(&scale, &ln(&div(&real_int(1), &sub(&v::exp(&d), &real_int(1)))))
}
fn h_c1(i: &Inst) -> Tri {
    v::le(&real_int(2), &c1(i))
}
fn h_delta_seq_pos(i: &Inst) -> Tri {
    let d = seq(i, "delta");
    all(FORALL_NAT, |k| v::lt(&real_int(0), &at(&d, k)))
}
fn h_start(i: &Inst) -> Tri {
    v::le(&at(&seq(i, "delta"), 0), &threshold(i))
}
fn h_recurrence(i: &Inst) -> Tri {
    let d = seq(i, "delta");
    let t = threshold(i);
    all(FORALL_NAT, |k| {
        v::le(&at(&d, k), &t).implies(v::le(&at(&d, k + 1), &mul(&c1(i), &powv(&at(&d, k), 2))))
    })
}
fn h_from_defect(i: &Inst) -> Tri {
    let m = seq(i, "m");
    all(FORALL_NAT, |k| v::le(&bound_k(i, k), &at(&m, k)))
}
fn c_conditional_gap(i: &Inst) -> Tri {
    let Some(n) = nat_of(i, "n") else { return Tri::U };
    let m = seq(i, "m");
    let rhs = div(&ln(&real_int(2)), &to_real(&nat(n as u64)));
    all(FORALL_NAT, |k| v::le(&rhs, &at(&m, k)))
}

fn mat(i: &Inst, k: &str) -> V {
    get(i, k)
}
/// Instances whose matrices are not `N × N` are ill-typed: every predicate is undecided on them.
fn dims_ok(i: &Inst, names: &[&str]) -> Option<usize> {
    let n = nat_of(i, "N")?;
    for k in names {
        match i.get(*k) {
            Some(V::Mat(m)) if m.rows == n && m.cols == n => {}
            _ => return None,
        }
    }
    Some(n)
}
fn mm(a: &V, b: &V) -> V {
    mul(a, b)
}
fn tr(a: &V) -> V {
    v::trace(a, M).unwrap_or(V::Undecided)
}
fn lin_trace(n: usize, complex: bool, f: impl Fn(V) -> V + 'static) -> V {
    let g = V::Fun(Rc::new(move |z: V| Ok(f(z))));
    v::linear_map_trace(&g, n, complex, M).unwrap_or(V::Undecided)
}
fn commutator(a: &V, b: &V) -> V {
    sub(&mm(a, b), &mm(b, a))
}
fn c_sandwich(i: &Inst) -> Tri {
    let Some(n) = dims_ok(i, &["A", "B"]) else { return Tri::U };
    let (a, b) = (mat(i, "A"), mat(i, "B"));
    let (a2, b2) = (a.clone(), b.clone());
    let t = lin_trace(n, false, move |z| mm(&mm(&a2, &z), &b2));
    v::eq(&t, &mm(&tr(&a), &tr(&b)))
}
fn c_ad_ad(i: &Inst) -> Tri {
    let Some(n) = dims_ok(i, &["X", "Y"]) else { return Tri::U };
    let (x, y) = (mat(i, "X"), mat(i, "Y"));
    let (x2, y2) = (x.clone(), y.clone());
    let t = lin_trace(n, false, move |z| commutator(&x2, &commutator(&y2, &z)));
    let rhs = sub(&mul(&mul(&real_int(2), &to_real(&nat(n as u64))), &tr(&mm(&x, &y))), &mul(&real_int(2), &mm(&tr(&x), &tr(&y))));
    v::eq(&t, &rhs)
}

fn cx(r: crate::num::Iv, i: i64) -> V {
    V::Cx(r.lo, num::qi(i))
}
fn is_su(n: usize, x: &V) -> Tri {
    let V::Mat(m) = x else { return Tri::U };
    if m.rows != n || m.cols != n {
        return Tri::F;
    }
    let skew = v::eq(&v::conj_transpose(x).unwrap_or(V::Undecided), &v::negate(x, M).unwrap_or(V::Undecided));
    let zero = V::Cx(num::qi(0), num::qi(0));
    skew.and(v::eq(&tr(x), &zero))
}
/// An orthogonal basis of su(N) for the inner product `Re tr(A* B)`.
fn su_basis(n: usize) -> Vec<V> {
    let zero = || V::Cx(num::qi(0), num::qi(0));
    let unit = |entries: Vec<(usize, usize, V)>| {
        let mut e = vec![zero(); n * n];
        for (a, b, x) in entries {
            e[a * n + b] = x;
        }
        V::Mat(Rc::new(Mat { rows: n, cols: n, e }))
    };
    let mut out = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            out.push(unit(vec![(a, b, cx(Iv::int(1), 0)), (b, a, cx(Iv::int(-1), 0))]));
            out.push(unit(vec![(a, b, cx(Iv::int(0), 1)), (b, a, cx(Iv::int(0), 1))]));
        }
    }
    for l in 1..n {
        let mut d = Vec::new();
        for a in 0..l {
            d.push((a, a, cx(Iv::int(0), 1)));
        }
        d.push((l, l, cx(Iv::int(0), -(l as i64))));
        out.push(unit(d));
    }
    out
}
fn re_inner(a: &V, b: &V) -> V {
    v::real_part(&tr(&mm(&v::conj_transpose(a).unwrap_or(V::Undecided), b)))
}
/// Killing form of su(N): the trace of the real-linear map `Z ↦ [X,[Y,Z]]` on su(N).
fn killing(n: usize, x: &V, y: &V) -> V {
    let mut acc = real_int(0);
    for b in su_basis(n) {
        let image = commutator(x, &commutator(y, &b));
        acc = add(&acc, &div(&re_inner(&b, &image), &re_inner(&b, &b)));
    }
    acc
}
fn h_n_ge_2(i: &Inst) -> Tri {
    v::le(&nat(2), &get(i, "N"))
}
fn h_x_su(i: &Inst) -> Tri {
    nat_of(i, "N").map(|n| is_su(n, &mat(i, "X"))).unwrap_or(Tri::U)
}
fn h_y_su(i: &Inst) -> Tri {
    nat_of(i, "N").map(|n| is_su(n, &mat(i, "Y"))).unwrap_or(Tri::U)
}
fn c_killing(i: &Inst) -> Tri {
    let Some(n) = dims_ok(i, &["X", "Y"]) else { return Tri::U };
    if h_x_su(i).and(h_y_su(i)) != Tri::T {
        return Tri::U; // K is defined on su(N) only
    }
    let (x, y) = (mat(i, "X"), mat(i, "Y"));
    let rhs = mul(&mul(&real_int(2), &to_real(&nat(n as u64))), &tr(&mm(&x, &y)));
    v::eq(&killing(n, &x, &y), &rhs)
}
fn c_casimir(i: &Inst) -> Tri {
    let Some(n) = dims_ok(i, &["X", "Y"]) else { return Tri::U };
    if h_x_su(i).and(h_y_su(i)) != Tri::T {
        return Tri::U;
    }
    let (x, y) = (mat(i, "X"), mat(i, "Y"));
    let nn = to_real(&nat(n as u64));
    let inner = mul(&v::negate(&div(&real_int(1), &nn), M).unwrap_or(V::Undecided), &v::real_part(&tr(&mm(&x, &y))));
    let c = mul(&real_int(2), &powv(&nn, 2));
    v::eq(&killing(n, &x, &y), &mul(&v::negate(&c, M).unwrap_or(V::Undecided), &inner))
}

// ---------- samplers of instances satisfying the hypotheses ----------

fn spectrum(r: &mut Rng, second_positive: bool, small: bool) -> V {
    let top = r.pos_rational();
    let len = 1 + r.below(6) as usize;
    let mut out = vec![real(top.clone())];
    for j in 1..len {
        let u = if small { num::mul(&r.unit_rational(), &num::q(1, 3)) } else { r.unit_rational() };
        let x = match r.below(6) {
            0 if !(second_positive && j == 1) => num::qi(0),
            1 if !small => top.clone(),
            _ => num::mul(&top, &u),
        };
        let x = if second_positive && j == 1 && x.is_zero() { num::mul(&top, &num::q(1, 7)) } else { x };
        out.push(real(x));
    }
    if second_positive && out.len() < 2 {
        out.push(real(num::mul(&top, &num::q(1, 2))));
    }
    V::Seq(Rc::new(out))
}
fn s_purity(r: &mut Rng) -> Inst {
    let mut i = Inst::new();
    i.insert("lam".into(), spectrum(r, false, false));
    i.insert("n".into(), nat(1 + r.below(4)));
    i
}
fn s_gap_purity(r: &mut Rng) -> Inst {
    let mut i = Inst::new();
    let small = r.below(2) == 0;
    i.insert("lam".into(), spectrum(r, true, small));
    i.insert("n".into(), nat(1 + r.below(4)));
    i
}
fn s_defect(r: &mut Rng) -> Inst {
    let mut i = Inst::new();
    let mut z = r.pos_rational();
    for k in 1..=5 {
        i.insert(format!("z{k}"), real(z.clone()));
        z = num::round_to(&num::mul(&num::mul(&z, &z), &r.unit_rational()), false, 64);
        if z.is_zero() {
            z = num::q(1, 1 << 20);
        }
    }
    i
}
fn s_gap_defect(r: &mut Rng) -> Inst {
    loop {
        let mut i = Inst::new();
        let lam = spectrum(r, true, true);
        let n = 1 + r.below(3) as usize;
        let p = match &lam {
            V::Seq(s) => purity(s, n),
            _ => unreachable!(),
        };
        let V::Real(p) = p else { continue };
        let Some(min_delta) = num::log_iv(&p, false).map(|x| x.neg()) else { continue };
        let (l2, _) = num::log2_bounds();
        let lo = num::round_to(&min_delta.hi, true, 40);
        let lo = if num::sign(&lo) <= 0 { num::q(1, 1 << 20) } else { lo };
        if num::cmp(&lo, &l2) != std::cmp::Ordering::Less {
            continue;
        }
        let d = num::add(&lo, &num::mul(&num::sub(&l2, &lo), &num::mul(&r.unit_rational(), &num::q(9, 10))));
        let d = num::round_to(&d, true, 40);
        if num::cmp(&d, &l2) != std::cmp::Ordering::Less {
            continue;
        }
        i.insert("lam".into(), lam);
        i.insert("n".into(), nat(n as u64));
        i.insert("Delta".into(), real(d));
        return i;
    }
}
fn s_conditional(r: &mut Rng) -> Inst {
    let mut i = Inst::new();
    let c = num::add(&num::qi(2), &num::mul(&r.unit_rational(), &num::qi(3)));
    let n = 1 + r.below(3);
    let mut d = num::mul(&num::div(&num::qi(1), &num::mul(&num::qi(2), &c)), &r.unit_rational_pos());
    let mut deltas = Vec::new();
    for _ in 0..=FORALL_NAT {
        deltas.push(real(d.clone()));
        let next = num::mul(&num::mul(&c, &num::mul(&d, &d)), &r.unit_rational_pos());
        d = num::round_to(&next, false, 64);
    }
    i.insert("c1".into(), real(c));
    i.insert("n".into(), nat(n));
    i.insert("delta".into(), V::Seq(Rc::new(deltas)));
    let mut ms = Vec::new();
    for k in 0..=FORALL_NAT {
        let b = bound_k(&i, k);
        let V::Real(b) = b else { ms.push(real_int(100)); continue };
        let slack = if r.below(3) == 0 { num::qi(0) } else { num::mul(&r.unit_rational(), &num::q(1, 4)) };
        ms.push(real(num::add(&num::round_to(&b.hi, true, 64), &slack)));
    }
    i.insert("m".into(), V::Seq(Rc::new(ms)));
    i
}
fn s_trace(r: &mut Rng) -> Inst {
    let mut i = Inst::new();
    let n = 1 + r.below(3) as usize;
    i.insert("N".into(), nat(n as u64));
    for k in ["A", "B", "X", "Y"] {
        i.insert(k.into(), r.matrix(n, false));
    }
    i
}
pub fn random_su(r: &mut Rng, n: usize) -> V {
    let m = r.matrix(n, true);
    let x = sub(&m, &v::conj_transpose(&m).unwrap());
    let t = tr(&x);
    let shift = div(&t, &to_real(&nat(n as u64)));
    sub(&x, &mul(&shift, &v::identity(n, Some(&V::Cx(num::qi(0), num::qi(0))))))
}
fn s_casimir(r: &mut Rng) -> Inst {
    let mut i = Inst::new();
    let n = 2 + r.below(3) as usize;
    i.insert("N".into(), nat(n as u64));
    i.insert("X".into(), random_su(r, n));
    i.insert("Y".into(), random_su(r, n));
    i
}

pub fn statements() -> Vec<Paper> {
    let lam = || PVar { name: "lam", aliases: &["lam", "lambda", "eig", "l", "ev"], kind: Kind::Seq1 };
    let n = || PVar { name: "n", aliases: &["n"], kind: Kind::Nat };
    let spectrum_h: Vec<(&'static str, Pred)> =
        vec![("summable", h_summable), ("λ₁ > 0", h_top_pos), ("0 ≤ λ_j ≤ λ₁", h_bounds), ("n ≥ 1", h_n_pos)];
    let mut gap_h = spectrum_h.clone();
    gap_h.insert(2, ("λ₂ > 0", h_second_pos));
    let mut defect_h = gap_h.clone();
    defect_h.extend([("Δ > 0", h_delta_pos as Pred), ("Δ < log 2", h_delta_small), ("P(n) ≥ e^{-Δ}", h_purity_vs_defect)]);
    let mats = |names: &[&'static str], complex: bool| -> Vec<PVar> {
        let mut v = vec![PVar { name: "N", aliases: &["N"], kind: Kind::Nat }];
        for k in names {
            let aliases: &'static [&'static str] = match *k {
                "A" => &["A"],
                "B" => &["B"],
                "X" => &["X"],
                _ => &["Y"],
            };
            v.push(PVar { name: k, aliases, kind: Kind::Mat { complex } });
        }
        v
    };
    let _ = spectrum_hyps;
    vec![
        Paper {
            name: "casimir_su",
            vars: mats(&["X", "Y"], true),
            hyps: vec![("N ≥ 2", h_n_ge_2), ("X ∈ su(N)", h_x_su), ("Y ∈ su(N)", h_y_su)],
            concl: vec![("K(X,Y) = 2N tr(XY)", c_killing), ("K = -C⟨X,Y⟩, C = 2N²", c_casimir)],
            sample: s_casimir,
        },
        Paper {
            name: "trace_ad_ad",
            vars: mats(&["A", "B", "X", "Y"], false),
            hyps: vec![],
            concl: vec![("tr(Z ↦ AZB) = tr A tr B", c_sandwich), ("tr(ad X ad Y) = 2N tr(XY) - 2 tr X tr Y", c_ad_ad)],
            sample: s_trace,
        },
        Paper {
            name: "purity_bounds",
            vars: vec![lam(), n()],
            hyps: spectrum_h.clone(),
            concl: vec![
                ("P = (1+q)/(1+s)²", c_purity_formula),
                ("P ≤ 1", c_purity_le_one),
                ("s ≤ 1/P - 1", c_lower),
                ("1/P - 1 ≤ 2s + s²", c_upper),
                ("s_{2n} ≤ s²", c_doubling),
            ],
            sample: s_purity,
        },
        Paper {
            name: "gap_from_purity",
            vars: vec![lam(), n()],
            hyps: gap_h,
            concl: vec![("e^{-n m} ≤ 1/P - 1", c_exp_gap), ("P > 1/2 → m ≥ (1/n) log(1/(1/P - 1))", c_gap_bound)],
            sample: s_gap_purity,
        },
        Paper {
            name: "defect_telescope",
            vars: ["z1", "z2", "z3", "z4", "z5"]
                .iter()
                .map(|k| PVar {
                    name: k,
                    aliases: match *k {
                        "z1" => &["z1"],
                        "z2" => &["z2"],
                        "z3" => &["z3"],
                        "z4" => &["z4"],
                        _ => &["z5"],
                    },
                    kind: Kind::Real,
                })
                .collect(),
            hyps: vec![("z_k > 0", h_z_pos), ("P_k ≤ 1", h_p_le_one)],
            concl: vec![("e^{-Δ} = P₁⁸P₂⁴P₃²P₄", c_telescope), ("Δ ≥ 0", c_defect_nonneg), ("P₄ ≥ e^{-Δ}", c_last)],
            sample: s_defect,
        },
        Paper {
            name: "gap_from_defect",
            vars: vec![lam(), n(), PVar { name: "Delta", aliases: &["Delta", "D"], kind: Kind::Real }],
            hyps: defect_h,
            concl: vec![("log(λ₁/λ₂) ≥ (1/n) log(1/(e^Δ - 1))", c_gap_from_defect)],
            sample: s_gap_defect,
        },
        Paper {
            name: "conditional_gap",
            vars: vec![
                PVar { name: "c1", aliases: &["c1", "c"], kind: Kind::Real },
                n(),
                PVar { name: "delta", aliases: &["delta", "d"], kind: Kind::Seq0 },
                PVar { name: "m", aliases: &["m", "gap", "g"], kind: Kind::Seq0 },
            ],
            hyps: vec![
                ("c₁ ≥ 2", h_c1),
                ("n ≥ 1", h_n_pos),
                ("δ_k > 0", h_delta_seq_pos),
                ("δ₀ ≤ 1/(2c₁)", h_start),
                ("recurrence", h_recurrence),
                ("m_k ≥ bound", h_from_defect),
            ],
            concl: vec![("m_k ≥ log 2 / n", c_conditional_gap)],
            sample: s_conditional,
        },
    ]
}

pub fn by_name(name: &str) -> Option<Paper> {
    statements().into_iter().chain(toys()).find(|p| p.name == name)
}

// ---------- development statements (formal/Audit/Toy.lean); not part of the benchmark ----------

fn t_a_nonneg(i: &Inst) -> Tri {
    v::le(&real_int(0), &get(i, "a"))
}
fn t_a_le_b(i: &Inst) -> Tri {
    v::le(&get(i, "a"), &get(i, "b"))
}
fn t_sq(i: &Inst) -> Tri {
    v::le(&powv(&get(i, "a"), 2), &powv(&get(i, "b"), 2))
}
fn t_x_gt_1(i: &Inst) -> Tri {
    v::lt(&real_int(1), &get(i, "x"))
}
fn t_log_pos(i: &Inst) -> Tri {
    v::lt(&real_int(0), &ln(&get(i, "x")))
}
fn t_exp_log(i: &Inst) -> Tri {
    v::eq(&v::exp(&ln(&get(i, "x"))), &get(i, "x"))
}
fn s_toy_sq(r: &mut Rng) -> Inst {
    let a = r.pos_rational();
    let b = num::add(&a, &r.pos_rational());
    Inst::from([("a".to_string(), real(a)), ("b".to_string(), real(b))])
}
fn s_toy_log(r: &mut Rng) -> Inst {
    Inst::from([("x".to_string(), real(num::add(&num::qi(1), &r.pos_rational())))])
}
pub fn toys() -> Vec<Paper> {
    vec![
        Paper {
            name: "toy_sq",
            vars: vec![PVar { name: "a", aliases: &["a"], kind: Kind::Real }, PVar { name: "b", aliases: &["b"], kind: Kind::Real }],
            hyps: vec![("0 ≤ a", t_a_nonneg), ("a ≤ b", t_a_le_b)],
            concl: vec![("a² ≤ b²", t_sq)],
            sample: s_toy_sq,
        },
        Paper {
            name: "toy_log",
            vars: vec![PVar { name: "x", aliases: &["x"], kind: Kind::Real }],
            hyps: vec![("1 < x", t_x_gt_1)],
            concl: vec![("0 < log x", t_log_pos), ("exp (log x) = x", t_exp_log)],
            sample: s_toy_log,
        },
    ]
}
