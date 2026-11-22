//! Compare a Lean statement with an informal statement on instances.
use crate::lean::{bind, defined_value, defining_hypotheses, Env, Part, Statement};
use crate::paper::{Inst, Kind, Paper};
use crate::rng::Rng;
use crate::value::{self as v, nat, real, Tri, Ty, V};
use crate::num;
use serde_json::{json, Value as J};
use std::collections::HashMap;
use std::rc::Rc;

/// Lean binder names normalised: subscripts to digits, Greek letters spelled out, underscores dropped.
pub fn normalise(name: &str) -> String {
    let mut out = String::new();
    for c in name.chars() {
        match c {
            '₀'..='₉' => out.push(char::from_digit(c as u32 - '₀' as u32, 10).unwrap()),
            'λ' => out.push_str("lam"),
            'δ' => out.push_str("delta"),
            'Δ' => out.push_str("Delta"),
            '_' | '\'' => {}
            c => out.push(c),
        }
    }
    out
}

#[derive(Clone, Debug)]
struct Link {
    paper_var: usize,
    lean_id: String,
    offset: usize,
}

fn candidates(st: &Statement, p: &Paper) -> (Vec<Vec<Link>>, Vec<usize>) {
    let lean: Vec<_> = st.vars().map(|b| (b.id.clone(), normalise(&b.name))).collect();
    let mut options: Vec<Vec<Link>> = vec![vec![]];
    let mut unmatched = Vec::new();
    for (k, pv) in p.vars.iter().enumerate() {
        let found = lean.iter().find(|(_, n)| pv.aliases.iter().any(|a| a == n));
        match found {
            None => unmatched.push(k),
            Some((id, _)) => {
                let offsets: &[usize] = if pv.kind == Kind::Seq1 { &[1, 0] } else { &[0] };
                let mut next = Vec::new();
                for opt in &options {
                    for &o in offsets {
                        let mut o2 = opt.clone();
                        o2.push(Link { paper_var: k, lean_id: id.clone(), offset: o });
                        next.push(o2);
                    }
                }
                options = next;
            }
        }
    }
    (options, unmatched)
}

/// A random value of a Lean type.
fn sample_type(r: &mut Rng, t: &Ty) -> V {
    match t {
        Ty::Nat => nat(r.below(6)),
        Ty::Int => V::Int(num_bigint::BigInt::from(r.small_int())),
        Ty::Real => r.real(),
        Ty::Complex => V::Cx(num::qi(r.small_int()), num::qi(r.small_int())),
        Ty::Fin(n) => nat(r.below(*n as u64)),
        Ty::Fun(a, b) if matches!(**a, Ty::Nat) && matches!(**b, Ty::Real) => {
            let len = r.below(7) as usize;
            { let c = r.coin(); r.sequence(len, c) }
        }
        Ty::Fun(a, b) if matches!(**a, Ty::Nat) && matches!(**b, Ty::Nat) => {
            V::Seq(Rc::new((0..r.below(7)).map(|_| nat(r.below(6))).collect()))
        }
        Ty::Matrix(a, b, ring) => match (a.size(), b.size()) {
            (Some(n), Some(m)) if n == m => r.matrix(n, matches!(**ring, Ty::Complex)),
            _ => V::Undecided,
        },
        Ty::Type => V::Ty(Ty::Real),
        _ => V::Undecided,
    }
}

/// A random instance of the informal variables, ignoring the hypotheses.
fn generic_instance(r: &mut Rng, p: &Paper) -> Inst {
    let mut i = Inst::new();
    for pv in &p.vars {
        let val = match pv.kind {
            Kind::Real => r.real(),
            Kind::Nat => nat(if pv.name == "N" { r.below(5) } else { r.below(6) }),
            Kind::Seq1 => {
                let len = 1 + r.below(6) as usize;
                { let c = r.coin(); r.sequence(len, c) }
            }
            Kind::Seq0 => { let c = r.coin(); r.sequence(1 + crate::lean::FORALL_NAT, c) }
            Kind::Mat { complex } => {
                let n = i.get("N").and_then(|x| x.as_usize()).unwrap_or(2);
                if complex && r.coin() { crate::paper::random_su(r, n) } else { r.matrix(n, complex) }
            }
        };
        i.insert(pv.name.to_string(), val);
    }
    i
}

/// Change one variable of an instance a little.
fn perturb(r: &mut Rng, p: &Paper, mut i: Inst) -> Inst {
    for _ in 0..1 + r.below(2) {
        let pv = &p.vars[r.below(p.vars.len() as u64) as usize];
        let old = i.get(pv.name).cloned().unwrap_or(V::Undecided);
        let new = match (&old, pv.kind) {
            (V::Real(x), _) => match r.below(5) {
                0 => real(num::qi(0)),
                1 => V::Real(x.neg()),
                2 => V::Real(x.mul(&num::Iv::int(2))),
                3 => V::Real(x.mul(&num::Iv::exact(num::q(1, 2)))),
                _ => r.real(),
            },
            (V::Nat(k), _) => {
                let k = k.clone();
                match r.below(3) {
                    0 => nat(0),
                    1 => V::Nat(k + 1),
                    _ => V::Nat(if k > num_bigint::BigInt::from(0) { k - 1 } else { k }),
                }
            }
            (V::Seq(s), _) => {
                let mut s = s.as_ref().clone();
                if s.is_empty() {
                    s.push(r.real());
                } else {
                    let j = r.below(s.len() as u64) as usize;
                    s[j] = match r.below(5) {
                        0 => real(num::qi(0)),
                        1 => v::negate(&s[j], v::Sem::Math).unwrap_or(V::Undecided),
                        2 => s[0].clone(),
                        3 => v::arith(v::Op::Mul, &s[j], &real(num::qi(3)), v::Sem::Math).unwrap_or(V::Undecided),
                        _ => r.real(),
                    };
                }
                V::Seq(Rc::new(s))
            }
            (V::Mat(m), Kind::Mat { complex }) => {
                let mut m2 = m.as_ref().clone();
                if !m2.e.is_empty() {
                    let k = r.below(m2.e.len() as u64) as usize;
                    m2.e[k] = if complex { V::Cx(num::qi(r.small_int()), num::qi(r.small_int())) } else { real(num::qi(r.small_int())) };
                }
                V::Mat(Rc::new(m2))
            }
            _ => old,
        };
        i.insert(pv.name.to_string(), new);
    }
    i
}

fn instance(r: &mut Rng, p: &Paper, t: usize) -> Inst {
    match t % 3 {
        0 => (p.sample)(r),
        1 => {
            let base = (p.sample)(r);
            perturb(r, p, base)
        }
        _ => generic_instance(r, p),
    }
}

fn paper_eval(p: &Paper, i: &Inst) -> (Tri, Vec<Tri>, Vec<Tri>) {
    let hs: Vec<Tri> = p.hyps.iter().map(|(_, h)| h(i)).collect();
    let h = hs.iter().fold(Tri::T, |a, b| a.and(*b));
    let c = p.concl.iter().map(|(_, c)| c(i)).collect();
    (h, hs, c)
}

/// The Lean value of a linked informal variable.
fn lean_value(link: &Link, p: &Paper, i: &Inst, ty: &Ty, slot: &V) -> Option<V> {
    let pv = &p.vars[link.paper_var];
    let val = i.get(pv.name)?.clone();
    match (pv.kind, ty) {
        (Kind::Seq1, _) => {
            let V::Seq(s) = &val else { return None };
            if link.offset == 1 {
                Some(val)
            } else {
                let mut out = vec![slot.clone()];
                out.extend(s.iter().cloned());
                Some(V::Seq(Rc::new(out)))
            }
        }
        (Kind::Nat, Ty::Real) => match val {
            V::Nat(k) => Some(real(num::qbig(k, num_bigint::BigInt::from(1)))),
            other => Some(other),
        },
        (Kind::Real, Ty::Nat) => {
            let k = val.as_usize()?;
            Some(nat(k as u64))
        }
        _ => Some(val),
    }
}

struct Trial {
    hl: Tri,
    cl: Vec<Tri>,
}

/// Build the Lean environment for an instance; the free (Lean-only) variables are drawn from `r`.
fn lean_trial(st: &Statement, p: &Paper, links: &[Link], i: &Inst, r: &mut Rng) -> Result<(Trial, Env), String> {
    let defs = defining_hypotheses(st);
    let defined: Vec<&str> = defs.iter().map(|d| d.var_id.as_str()).collect();
    let mut env: Env = Rc::new(HashMap::new());
    let slot = r.real();
    for part in &st.parts {
        match part {
            Part::Var(b) => {
                let ty = st.ev.ty(&b.ty, &env)?;
                let val = if let Some(link) = links.iter().find(|l| l.lean_id == b.id) {
                    let val = lean_value(link, p, i, &ty, &slot).ok_or("informal value does not fit the Lean type")?;
                    // an instance whose matrices are not the shape their Lean type says is ill-typed: skip it
                    if let (Ty::Matrix(r, c, _), V::Mat(m)) = (&ty, &val) {
                        if r.size() != Some(m.rows) || c.size() != Some(m.cols) {
                            return Err("ill-typed instance: matrix shape differs from its Lean type".into());
                        }
                    }
                    val
                } else if defined.contains(&b.id.as_str()) {
                    V::Undecided
                } else {
                    sample_type(r, &ty)
                };
                env = bind(&env, &b.id, val);
            }
            Part::Hyp(b) | Part::Inst(b) => env = bind(&env, &b.id, V::Opaque),
        }
        // a defining hypothesis sets its variable as soon as everything it uses is bound
        for d in &defs {
            if matches!(env.get(&d.var_id), Some(V::Undecided)) {
                if let Ok(val) = defined_value(st, d, &env) {
                    env = bind(&env, &d.var_id, val);
                }
            }
        }
    }
    let mut hl = Tri::T;
    for b in st.hyps() {
        hl = hl.and(st.ev.prop(&b.ty, &env)?);
    }
    let cl = st.conjuncts(&env)?;
    Ok((Trial { hl, cl }, env))
}

#[derive(Default)]
struct Stats {
    trials: usize,
    errors: usize,
    first_error: Option<String>,
    decided: usize,
    hl_true: usize,
    false_count: usize,
    false_example: Option<String>,
    stronger: usize,
    stronger_example: Option<String>,
    weaker: usize,
    agree: Vec<Vec<usize>>,
    disagree: Vec<Vec<usize>>,
    disagree_example: HashMap<(usize, usize), String>,
}

fn describe(i: &Inst) -> String {
    i.iter().map(|(k, v)| format!("{k}={v:?}")).collect::<Vec<_>>().join(", ")
}

fn run(st: &Statement, p: &Paper, links: &[Link], trials: usize, seed: u64) -> Stats {
    let mut r = Rng::new(seed);
    let nl = st.conjunct_count();
    let np = p.concl.len();
    let mut s = Stats { agree: vec![vec![0; np]; nl], disagree: vec![vec![0; np]; nl], ..Default::default() };
    for t in 0..trials {
        let inst = instance(&mut r, p, t);
        let (hp, _, cp) = paper_eval(p, &inst);
        let mut res = lean_trial(st, p, links, &inst, &mut r);
        if let Ok((tr, _)) = &res {
            if hp == Tri::T && tr.hl != Tri::T {
                // search the Lean-only variables for values that satisfy the Lean hypotheses
                let mut any_undecided = tr.hl == Tri::U;
                let mut found = None;
                for _ in 0..20 {
                    match lean_trial(st, p, links, &inst, &mut r) {
                        Ok((t2, e2)) if t2.hl == Tri::T => {
                            found = Some(Ok((t2, e2)));
                            break;
                        }
                        Ok((t2, _)) => any_undecided |= t2.hl == Tri::U,
                        Err(_) => {}
                    }
                }
                match found {
                    Some(f) => res = f,
                    None => {
                        if !any_undecided {
                            s.stronger += 1;
                            s.stronger_example.get_or_insert_with(|| describe(&inst));
                        }
                    }
                }
            }
        }
        s.trials += 1;
        let (tr, _) = match res {
            Ok(x) => x,
            Err(e) => {
                s.errors += 1;
                s.first_error.get_or_insert(e);
                continue;
            }
        };
        if tr.cl.len() != nl {
            s.errors += 1;
            continue;
        }
        if tr.hl.decided() && hp.decided() {
            s.decided += 1;
        }
        if tr.hl == Tri::T {
            s.hl_true += 1;
            if tr.cl.iter().any(|c| *c == Tri::F) {
                s.false_count += 1;
                s.false_example.get_or_insert_with(|| describe(&inst));
            }
            if hp == Tri::F {
                s.weaker += 1;
            }
        }
        for (a, ca) in tr.cl.iter().enumerate() {
            for (b, cb) in cp.iter().enumerate() {
                if ca.decided() && cb.decided() {
                    if ca == cb {
                        s.agree[a][b] += 1;
                    } else {
                        s.disagree[a][b] += 1;
                        s.disagree_example.entry((a, b)).or_insert_with(|| describe(&inst));
                    }
                }
            }
        }
    }
    s
}

fn matches_of(s: &Stats, a: usize) -> Vec<usize> {
    (0..s.agree.get(a).map(|r| r.len()).unwrap_or(0)).filter(|&b| s.disagree[a][b] == 0 && s.agree[a][b] > 0).collect()
}

pub fn audit(st: &Statement, p: &Paper, formalizer: bool, trials: usize, seed: u64) -> J {
    let (options, unmatched) = candidates(st, p);
    let mut best: Option<(usize, usize, Vec<String>, Stats, Vec<Link>)> = None;
    for links in options {
        let s = run(st, p, &links, trials, seed);
        let mut flags = Vec::new();
        if s.errors == s.trials {
            flags.push("ERROR".to_string());
        }
        if s.false_count > 0 {
            flags.push("FALSE".into());
        }
        if s.hl_true == 0 {
            flags.push("VACUOUS".into());
        }
        if s.stronger > 0 {
            flags.push("STRONGER".into());
        }
        let nl = st.conjunct_count();
        if (0..nl).any(|a| matches_of(&s, a).is_empty()) {
            flags.push("MISMATCH".into());
        }
        let covered: Vec<bool> = (0..p.concl.len()).map(|b| (0..nl).any(|a| matches_of(&s, a).contains(&b))).collect();
        if formalizer && covered.iter().any(|c| !c) {
            flags.push("MISSING".into());
        }
        if formalizer && !unmatched.is_empty() {
            flags.push("UNMATCHED".into());
        }
        if s.decided < 50 {
            flags.push("UNDECIDED".into());
        }
        let disagreements: usize = s.disagree.iter().flatten().sum();
        let better = match &best {
            None => true,
            Some((f, d, ..)) => (flags.len(), disagreements) < (*f, *d),
        };
        if better {
            best = Some((flags.len(), disagreements, flags, s, links));
        }
    }
    let (_, _, flags, s, links) = best.expect("at least one mapping");
    let nl = st.conjunct_count();
    let lean_names: HashMap<String, String> = st.vars().map(|b| (b.id.clone(), b.name.clone())).collect();
    let conj: Vec<J> = (0..nl)
        .map(|a| {
            let m = matches_of(&s, a);
            json!({
                "lean_conjunct": a,
                "matches": m.iter().map(|b| p.concl[*b].0).collect::<Vec<_>>(),
                "disagreements": (0..p.concl.len()).map(|b| json!({
                    "informal": p.concl[b].0, "agree": s.agree[a][b], "disagree": s.disagree[a][b],
                    "example": s.disagree_example.get(&(a, b)),
                })).collect::<Vec<_>>(),
            })
        })
        .collect();
    let covered: Vec<&str> = (0..p.concl.len())
        .filter(|b| (0..nl).any(|a| matches_of(&s, a).contains(b)))
        .map(|b| p.concl[b].0)
        .collect();
    json!({
        "lean": st.name,
        "informal": p.name,
        "verdict": if flags.is_empty() { "PASS" } else { "FLAGGED" },
        "flags": flags,
        "weaker": s.weaker > 0,
        "mapping": links.iter().map(|l| json!({
            "informal": p.vars[l.paper_var].name, "lean": lean_names.get(&l.lean_id), "offset": l.offset})).collect::<Vec<_>>(),
        "unmatched_informal": unmatched.iter().map(|k| p.vars[*k].name).collect::<Vec<_>>(),
        "covered_informal_conjuncts": covered,
        "trials": s.trials, "errors": s.errors, "first_error": s.first_error,
        "decided": s.decided, "lean_hypotheses_true": s.hl_true,
        "false": s.false_count, "false_example": s.false_example,
        "stronger": s.stronger, "stronger_example": s.stronger_example,
        "weaker_count": s.weaker,
        "conjuncts": conj,
    })
}

/// Self-check of an informal statement and the necessity of each of its hypotheses.
pub fn scan(p: &Paper, trials: usize, seed: u64) -> J {
    let mut r = Rng::new(seed);
    let mut self_false = 0;
    let mut self_example = None;
    let mut sat = 0;
    let nh = p.hyps.len();
    let mut needed = vec![0usize; nh];
    let mut needed_example: Vec<Option<String>> = vec![None; nh];
    let mut others_only = vec![0usize; nh];
    for t in 0..trials {
        // perturbed instances dominate: they break one hypothesis while keeping the others
        let inst = match t % 4 {
            0 => (p.sample)(&mut r),
            3 => generic_instance(&mut r, p),
            _ => {
                let b = (p.sample)(&mut r);
                perturb(&mut r, p, b)
            }
        };
        let (h, hs, c) = paper_eval(p, &inst);
        let concl = c.iter().fold(Tri::T, |a, b| a.and(*b));
        if h == Tri::T {
            sat += 1;
            if concl == Tri::F {
                self_false += 1;
                self_example.get_or_insert_with(|| describe(&inst));
            }
        }
        for k in 0..nh {
            let others = hs.iter().enumerate().filter(|(j, _)| *j != k).fold(Tri::T, |a, (_, b)| a.and(*b));
            if others == Tri::T && hs[k] == Tri::F {
                others_only[k] += 1;
                if concl == Tri::F {
                    needed[k] += 1;
                    needed_example[k].get_or_insert_with(|| describe(&inst));
                }
            }
        }
    }
    json!({
        "informal": p.name,
        "instances_satisfying_hypotheses": sat,
        "self_check_failures": self_false, "self_check_example": self_example,
        "hypotheses": (0..nh).map(|k| json!({
            "hypothesis": p.hyps[k].0,
            "instances_violating_only_this": others_only[k],
            "counterexamples_without_it": needed[k],
            "needed": needed[k] > 0,
            "example": needed_example[k],
        })).collect::<Vec<_>>(),
    })
}
