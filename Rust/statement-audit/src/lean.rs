//! Lean terms exported by `formal/Audit/Export.lean`, and their evaluation with Lean semantics.
use crate::value::{self, apply, arith, nat, real, real_int, Mat, Op, Sem, Tri, Ty, V, R};
use crate::num;
use num_bigint::BigInt;
use num_traits::One;
use serde_json::Value as J;
use std::collections::HashMap;
use std::rc::Rc;

/// How many values a universal quantifier over `ℕ` is checked for.
pub const FORALL_NAT: usize = 8;
/// A `tsum` over `ℕ` is decided only if its terms vanish on `[SUPPORT, SUPPORT + 8)`; it then sums `[0, SUPPORT)`.
const SUPPORT: usize = 16;

#[derive(Clone, Debug)]
pub struct Binder {
    pub name: String,
    pub id: String,
    pub bi: String,
    pub prop: bool,
    pub ty: Rc<E>,
}

#[derive(Clone, Debug)]
pub enum E {
    Const(String),
    App(Rc<E>, Vec<Rc<E>>),
    Lam(Binder, Rc<E>),
    Pi(Binder, Rc<E>, bool),
    Let { id: String, name: String, value: Rc<E>, body: Rc<E> },
    Fvar { id: String, name: String },
    Nat(BigInt),
    Sort(bool),
    Other(String),
}

pub fn parse(j: &J) -> R<Rc<E>> {
    let k = j["k"].as_str().ok_or("node without kind")?;
    let s = |f: &str| j[f].as_str().unwrap_or("").to_string();
    Ok(Rc::new(match k {
        "const" => E::Const(s("n")),
        "app" => {
            let args = j["args"].as_array().ok_or("app args")?.iter().map(parse).collect::<R<Vec<_>>>()?;
            E::App(parse(&j["f"])?, args)
        }
        "lam" | "pi" => {
            let b = Binder {
                name: s("name"),
                id: s("id"),
                bi: s("bi"),
                prop: j["prop"].as_bool().unwrap_or(false),
                ty: parse(&j["type"])?,
            };
            let body = parse(&j["body"])?;
            if k == "lam" {
                E::Lam(b, body)
            } else {
                E::Pi(b, body, j["self_prop"].as_bool().unwrap_or(false))
            }
        }
        "let" => E::Let { id: s("id"), name: s("name"), value: parse(&j["value"])?, body: parse(&j["body"])? },
        "fvar" => E::Fvar { id: s("id"), name: s("name") },
        "nat" => E::Nat(s("v").parse::<BigInt>().map_err(|e| e.to_string())?),
        "sort" => E::Sort(j["prop"].as_bool().unwrap_or(false)),
        other => E::Other(other.to_string()),
    }))
}

#[derive(Clone, Debug)]
pub enum Def {
    Value(Rc<E>),
    Structure { params: usize, ctor: Rc<E> },
}

pub type Env = Rc<HashMap<String, V>>;

pub fn bind(env: &Env, id: &str, v: V) -> Env {
    let mut m = (**env).clone();
    m.insert(id.to_string(), v);
    Rc::new(m)
}

/// Evaluator for one exported statement.
#[derive(Clone)]
pub struct Ev {
    pub defs: Rc<HashMap<String, Def>>,
}

fn head(e: &E) -> Option<(&str, &[Rc<E>])> {
    match e {
        E::App(f, args) => match &**f {
            E::Const(c) => Some((c.as_str(), args.as_slice())),
            _ => None,
        },
        E::Const(c) => Some((c.as_str(), &[])),
        _ => None,
    }
}

impl Ev {
    pub fn eval(&self, e: &Rc<E>, env: &Env) -> R<V> {
        match &**e {
            E::Fvar { id, name } => env.get(id).cloned().ok_or_else(|| format!("unbound variable {name}")),
            E::Nat(n) => Ok(V::Nat(n.clone())),
            E::Sort(p) => Ok(V::Ty(if *p { Ty::Prop } else { Ty::Type })),
            E::Const(c) => self.call(c, &[], env),
            E::App(f, args) => {
                if let E::Const(c) = &**f {
                    return self.call(c, args, env);
                }
                let mut v = self.eval(f, env)?;
                for a in args {
                    v = apply(&v, self.eval(a, env)?)?;
                }
                Ok(v)
            }
            E::Lam(b, body) => {
                let (this, env, id, body) = (self.clone(), env.clone(), b.id.clone(), body.clone());
                Ok(V::Fun(Rc::new(move |x: V| this.eval(&body, &bind(&env, &id, x)))))
            }
            E::Pi(b, body, self_prop) => {
                if *self_prop {
                    Ok(V::Prop(self.forall(b, body, env)?))
                } else {
                    let dom = self.ty(&b.ty, env)?;
                    let cod = self.ty(body, &bind(env, &b.id, V::Opaque)).unwrap_or(Ty::Other("dependent".into()));
                    Ok(V::Ty(Ty::Fun(Box::new(dom), Box::new(cod))))
                }
            }
            E::Let { id, value, body, .. } => {
                let v = self.eval(value, env)?;
                self.eval(body, &bind(env, id, v))
            }
            E::Other(k) => Err(format!("unsupported term kind {k}")),
        }
    }

    pub fn ty(&self, e: &Rc<E>, env: &Env) -> R<Ty> {
        match self.eval(e, env)? {
            V::Ty(t) => Ok(t),
            other => Err(format!("expected a type, got {other:?}")),
        }
    }

    pub fn prop(&self, e: &Rc<E>, env: &Env) -> R<Tri> {
        Ok(self.eval(e, env)?.as_tri())
    }

    /// `∀ x : A, B x` or `A → B` as a proposition.
    fn forall(&self, b: &Binder, body: &Rc<E>, env: &Env) -> R<Tri> {
        if b.bi == "inst" {
            return self.prop(body, &bind(env, &b.id, V::Opaque));
        }
        if b.prop {
            let a = self.prop(&b.ty, env)?;
            if a == Tri::F {
                return Ok(Tri::T);
            }
            let c = self.prop(body, &bind(env, &b.id, V::Opaque))?;
            return Ok(a.implies(c));
        }
        let dom = self.ty(&b.ty, env)?;
        let values: Vec<V> = match &dom {
            Ty::Nat => (0..FORALL_NAT as u64).map(nat).collect(),
            Ty::Fin(n) => (0..*n as u64).map(nat).collect(),
            Ty::Subtype(base, p) if matches!(**base, Ty::Nat) => {
                let mut v = Vec::new();
                let mut k = 0u64;
                while v.len() < FORALL_NAT && k < 64 {
                    match p(nat(k))?.as_tri() {
                        Tri::T => v.push(nat(k)),
                        Tri::F => {}
                        Tri::U => return Ok(Tri::U),
                    }
                    k += 1;
                }
                v
            }
            _ => return Ok(Tri::U),
        };
        let mut acc = Tri::T;
        for x in values {
            acc = acc.and(self.prop(body, &bind(env, &b.id, x))?);
            if acc == Tri::F {
                break;
            }
        }
        Ok(acc)
    }

    fn arg(&self, args: &[Rc<E>], i: usize, env: &Env) -> R<V> {
        let a = args.get(i).ok_or_else(|| format!("missing argument {i}"))?;
        self.eval(a, env)
    }
    fn last(&self, args: &[Rc<E>], env: &Env) -> R<V> {
        self.arg(args, args.len().checked_sub(1).ok_or("no arguments")?, env)
    }
    /// Apply remaining arguments beyond the first `n`.
    fn rest(&self, mut v: V, args: &[Rc<E>], n: usize, env: &Env) -> R<V> {
        for a in args.iter().skip(n) {
            v = apply(&v, self.eval(a, env)?)?;
        }
        Ok(v)
    }

    fn numeral(&self, t: &Ty, k: &BigInt) -> R<V> {
        Ok(match t {
            Ty::Nat => V::Nat(k.clone()),
            Ty::Int => V::Int(k.clone()),
            Ty::Real => real(num::qbig(k.clone(), BigInt::one())),
            Ty::Complex => V::Cx(num::qbig(k.clone(), BigInt::one()), num::qi(0)),
            Ty::Fin(_) => V::Nat(k.clone()),
            Ty::Matrix(a, b, r) => {
                let (n, m) = (a.size().ok_or("matrix index type")?, b.size().ok_or("matrix index type")?);
                let c = self.numeral(r, k)?;
                let z = self.numeral(r, &BigInt::from(0))?;
                let e = (0..n * m).map(|p| if p / m == p % m { c.clone() } else { z.clone() }).collect();
                V::Mat(Rc::new(Mat { rows: n, cols: m, e }))
            }
            Ty::Fun(_, cod) => {
                let c = self.numeral(cod, k)?;
                V::Fun(Rc::new(move |_| Ok(c.clone())))
            }
            other => return Err(format!("numeral of type {other:?}")),
        })
    }

    fn cast(&self, t: &Ty, v: V) -> R<V> {
        Ok(match (t, &v) {
            (_, V::Undecided) => V::Undecided,
            (Ty::Real, V::Nat(n) | V::Int(n)) => real(num::qbig(n.clone(), BigInt::one())),
            (Ty::Complex, V::Nat(n) | V::Int(n)) => V::Cx(num::qbig(n.clone(), BigInt::one()), num::qi(0)),
            (Ty::Int, V::Nat(n)) => V::Int(n.clone()),
            (Ty::Complex, V::Real(x)) => match x.exact_value() {
                Some(q) => V::Cx(q.clone(), num::qi(0)),
                None => V::Undecided,
            },
            _ => v,
        })
    }

    /// Sum of `f` over the domain type `dom`; `None` when undecided.
    fn sum_over(&self, dom: &Ty, f: &V, total_if_infinite: bool) -> R<V> {
        let mut acc = real_int(0);
        let mut first = true;
        let mut add = |acc: &mut V, x: V| -> R<()> {
            if first {
                *acc = x;
                first = false;
                Ok(())
            } else {
                *acc = arith(Op::Add, acc, &x, Sem::Lean)?;
                Ok(())
            }
        };
        match dom {
            Ty::Fin(n) => {
                for i in 0..*n as u64 {
                    add(&mut acc, apply(f, nat(i))?)?;
                }
                Ok(acc)
            }
            Ty::Nat | Ty::Subtype(..) => {
                let member = |k: u64| -> R<Tri> {
                    match dom {
                        Ty::Subtype(_, p) => Ok(p(nat(k))?.as_tri()),
                        _ => Ok(Tri::T),
                    }
                };
                for k in SUPPORT as u64..SUPPORT as u64 + 8 {
                    match member(k)? {
                        Tri::F => continue,
                        Tri::U => return Ok(V::Undecided),
                        Tri::T => {}
                    }
                    if value::eq(&apply(f, nat(k))?, &real_int(0)) != Tri::T {
                        return Ok(if total_if_infinite { V::Undecided } else { V::Prop(Tri::U) });
                    }
                }
                for k in 0..SUPPORT as u64 {
                    match member(k)? {
                        Tri::F => continue,
                        Tri::U => return Ok(V::Undecided),
                        Tri::T => add(&mut acc, apply(f, nat(k))?)?,
                    }
                }
                Ok(acc)
            }
            other => Err(format!("sum over {other:?}")),
        }
    }

    fn domain_of(&self, f: &Rc<E>, env: &Env) -> R<Ty> {
        match &**f {
            E::Lam(b, _) => self.ty(&b.ty, env),
            _ => Err("summand is not a lambda".into()),
        }
    }

    fn call(&self, c: &str, args: &[Rc<E>], env: &Env) -> R<V> {
        let a = |i: usize| self.arg(args, i, env);
        match c {
            // types
            "Real" => Ok(V::Ty(Ty::Real)),
            "Nat" => Ok(V::Ty(Ty::Nat)),
            "Int" => Ok(V::Ty(Ty::Int)),
            "Complex" => Ok(V::Ty(Ty::Complex)),
            "Fin" => Ok(V::Ty(Ty::Fin(a(0)?.as_usize().ok_or("Fin size")?))),
            "Matrix" => {
                let (m, n, r) = (self.ty(&args[0], env)?, self.ty(&args[1], env)?, self.ty(&args[2], env)?);
                Ok(V::Ty(Ty::Matrix(Box::new(m), Box::new(n), Box::new(r))))
            }
            "Subtype" => {
                let base = self.ty(&args[0], env)?;
                let p = a(1)?;
                Ok(V::Ty(Ty::Subtype(Box::new(base), Rc::new(move |x| apply(&p, x)))))
            }
            "optParam" | "outParam" | "autoParam" => a(0),
            "LinearMap" | "Set" | "Finset" | "Module.End" => Ok(V::Ty(Ty::Other(c.to_string()))),
            // numbers
            "OfNat.ofNat" => {
                let t = self.ty(&args[0], env)?;
                match &*args[1] {
                    E::Nat(k) => self.numeral(&t, k),
                    _ => match a(1)? {
                        V::Nat(k) => self.numeral(&t, &k),
                        other => Err(format!("numeral literal {other:?}")),
                    },
                }
            }
            "Zero.zero" => self.numeral(&self.ty(&args[0], env)?, &BigInt::from(0)),
            "One.one" => self.numeral(&self.ty(&args[0], env)?, &BigInt::from(1)),
            "Nat.cast" | "NatCast.natCast" | "Int.cast" | "IntCast.intCast" => {
                let t = self.ty(&args[0], env)?;
                let v = self.cast(&t, a(2)?)?;
                self.rest(v, args, 3, env)
            }
            "Complex.ofReal" | "Complex.ofReal'" => self.cast(&Ty::Complex, a(0)?),
            "Complex.re" => Ok(value::real_part(&a(0)?)),
            "Nat.succ" => arith(Op::Add, &a(0)?, &nat(1), Sem::Lean),
            // arithmetic
            "HAdd.hAdd" | "HSub.hSub" | "HMul.hMul" | "HDiv.hDiv" | "HSMul.hSMul" => {
                let op = match c {
                    "HAdd.hAdd" => Op::Add,
                    "HSub.hSub" => Op::Sub,
                    "HDiv.hDiv" => Op::Div,
                    _ => Op::Mul,
                };
                let v = arith(op, &a(4)?, &a(5)?, Sem::Lean)?;
                self.rest(v, args, 6, env)
            }
            "Add.add" | "Sub.sub" | "Mul.mul" | "Div.div" => {
                let op = match c {
                    "Add.add" => Op::Add,
                    "Sub.sub" => Op::Sub,
                    "Mul.mul" => Op::Mul,
                    _ => Op::Div,
                };
                arith(op, &a(2)?, &a(3)?, Sem::Lean)
            }
            "HPow.hPow" => {
                let v = value::power(&a(4)?, &a(5)?, Sem::Lean)?;
                self.rest(v, args, 6, env)
            }
            "Monoid.npow" => value::power(&a(3)?, &a(2)?, Sem::Lean),
            "Neg.neg" => {
                let v = value::negate(&a(2)?, Sem::Lean)?;
                self.rest(v, args, 3, env)
            }
            "Inv.inv" => value::inverse(&a(2)?, Sem::Lean),
            "abs" => {
                let x = self.last(args, env)?;
                Ok(match value::le(&real_int(0), &x) {
                    Tri::T => x,
                    Tri::F => value::negate(&x, Sem::Lean)?,
                    Tri::U => V::Undecided,
                })
            }
            "Max.max" | "Min.min" => {
                let (x, y) = (a(2)?, a(3)?);
                let le = value::le(&x, &y);
                Ok(match (c, le) {
                    (_, Tri::U) => V::Undecided,
                    ("Max.max", Tri::T) | ("Min.min", Tri::F) => y,
                    _ => x,
                })
            }
            "Real.log" => Ok(value::log(&a(0)?, Sem::Lean)),
            "Real.exp" => Ok(value::exp(&a(0)?)),
            "Real.rpow" => value::power(&a(0)?, &a(1)?, Sem::Lean),
            "Real.pi" => Ok(V::Real(num::Iv::new(
                num::q(3_141_592_653_589_793, 1_000_000_000_000_000),
                num::q(3_141_592_653_589_794, 1_000_000_000_000_000),
            ))),
            // logic
            "True" => Ok(V::Prop(Tri::T)),
            "False" => Ok(V::Prop(Tri::F)),
            "And" => Ok(V::Prop(a(0)?.as_tri().and(a(1)?.as_tri()))),
            "Or" => Ok(V::Prop(a(0)?.as_tri().or(a(1)?.as_tri()))),
            "Not" => Ok(V::Prop(a(0)?.as_tri().not())),
            "Iff" => Ok(V::Prop(a(0)?.as_tri().iff(a(1)?.as_tri()))),
            "Eq" => Ok(V::Prop(value::eq(&a(1)?, &a(2)?))),
            "Ne" => Ok(V::Prop(value::eq(&a(1)?, &a(2)?).not())),
            "LE.le" => Ok(V::Prop(value::le(&a(2)?, &a(3)?))),
            "LT.lt" => Ok(V::Prop(value::lt(&a(2)?, &a(3)?))),
            "GE.ge" => Ok(V::Prop(value::le(&a(3)?, &a(2)?))),
            "GT.gt" => Ok(V::Prop(value::lt(&a(3)?, &a(2)?))),
            "ite" => match a(1)?.as_tri() {
                Tri::T => a(3),
                Tri::F => a(4),
                Tri::U => Ok(V::Undecided),
            },
            "Exists" => Ok(V::Prop(Tri::U)),
            // sums
            "tsum" => {
                let dom = self.domain_of(&args[4], env).or_else(|_| self.ty(&args[1], env))?;
                let f = a(4)?;
                self.sum_over(&dom, &f, true)
            }
            "Summable" => {
                let dom = self.domain_of(&args[4], env).or_else(|_| self.ty(&args[1], env))?;
                let f = a(4)?;
                Ok(match self.sum_over(&dom, &f, false)? {
                    V::Prop(t) => V::Prop(t),
                    V::Undecided => V::Prop(Tri::U),
                    _ => V::Prop(Tri::T),
                })
            }
            "Finset.sum" => {
                let f = a(4)?;
                let dom = match head(&args[3]) {
                    Some(("Finset.univ", sargs)) => self.ty(&sargs[0], env)?,
                    Some(("Finset.range", sargs)) => {
                        let n = self.arg(sargs, 0, env)?.as_usize().ok_or("range size")?;
                        Ty::Fin(n)
                    }
                    _ => return Err("Finset.sum over an unsupported set".into()),
                };
                self.sum_over(&dom, &f, true)
            }
            // matrices
            "Matrix.trace" => value::trace(&self.last(args, env)?, Sem::Lean),
            "Matrix.transpose" | "Matrix.conjTranspose" => {
                let m = a(if c == "Matrix.transpose" { 3 } else { 4 })?;
                let v = value::conj_transpose(&m)?;
                self.rest(v, args, if c == "Matrix.transpose" { 4 } else { 5 }, env)
            }
            // linear maps
            "DFunLike.coe" => {
                let f = a(4)?;
                self.rest(f, args, 5, env)
            }
            "LinearMap.mulLeft" | "LinearMap.mulRight" => {
                let x = self.last(args, env)?;
                let left = c == "LinearMap.mulLeft";
                Ok(V::Fun(Rc::new(move |z: V| {
                    if left { arith(Op::Mul, &x, &z, Sem::Lean) } else { arith(Op::Mul, &z, &x, Sem::Lean) }
                })))
            }
            "LinearMap.comp" => {
                let n = args.len();
                let (f, g) = (a(n - 2)?, a(n - 1)?);
                Ok(V::Fun(Rc::new(move |z: V| apply(&f, apply(&g, z)?))))
            }
            "LinearMap.id" => Ok(V::Fun(Rc::new(Ok))),
            "LinearMap.trace" => {
                let ring = self.ty(&args[0], env)?;
                // `LinearMap.trace R [inst] M [insts]`: the module is the first matrix type after the ring
                let m = args[1..]
                    .iter()
                    .find_map(|x| match self.ty(x, env) {
                        Ok(t @ Ty::Matrix(..)) => Some(t),
                        _ => None,
                    })
                    .ok_or("LinearMap.trace without a matrix module")?;
                let Ty::Matrix(r, _, _) = &m else {
                    return Err(format!("LinearMap.trace on {m:?}"));
                };
                let n = r.size().ok_or("matrix size")?;
                let complex = matches!(ring, Ty::Complex);
                Ok(V::Fun(Rc::new(move |f: V| value::linear_map_trace(&f, n, complex, Sem::Lean))))
            }
            "Subtype.val" => {
                let v = a(2)?;
                self.rest(v, args, 3, env)
            }
            "Subtype.mk" => a(2),
            "Fin.val" => a(1),
            _ => self.call_def(c, args, env),
        }
    }

    fn call_def(&self, c: &str, args: &[Rc<E>], env: &Env) -> R<V> {
        match self.defs.get(c) {
            Some(Def::Value(body)) => {
                let mut v = self.eval(body, &Rc::new(HashMap::new()))?;
                for x in args {
                    v = apply(&v, self.eval(x, env)?)?;
                }
                Ok(v)
            }
            Some(Def::Structure { params, ctor }) => {
                // the conjunction of the structure's fields, with its parameters instantiated
                let mut e = ctor.clone();
                let mut local = Rc::new(HashMap::new());
                let mut i = 0;
                let mut acc = Tri::T;
                while let E::Pi(b, body, _) = &*e.clone() {
                    if i < *params {
                        let v = self.eval(args.get(i).ok_or("structure parameter")?, env)?;
                        local = bind(&local, &b.id, v);
                    } else {
                        acc = acc.and(self.prop(&b.ty, &local)?);
                        local = bind(&local, &b.id, V::Opaque);
                    }
                    i += 1;
                    e = body.clone();
                }
                Ok(V::Prop(acc))
            }
            // Instance arguments and other constants with no numerical meaning: an opaque placeholder.
            // Using one as a number, function or proposition is an error or undecided, never a wrong value.
            None => Ok(V::Opaque),
        }
    }
}

/// A statement split into binders, hypotheses and conclusion.
#[derive(Clone, Debug)]
pub enum Part {
    /// A value or type variable.
    Var(Binder),
    /// A hypothesis.
    Hyp(Binder),
    /// An instance argument.
    Inst(Binder),
}

#[derive(Clone)]
pub struct Statement {
    pub name: String,
    pub parts: Vec<Part>,
    pub conclusion: Rc<E>,
    pub ev: Ev,
}

impl Statement {
    pub fn from_export(j: &J) -> R<Statement> {
        let mut defs = HashMap::new();
        for d in j["defs"].as_array().ok_or("defs")? {
            let name = d["name"].as_str().unwrap_or("").to_string();
            match d["kind"].as_str() {
                Some("def") => {
                    defs.insert(name, Def::Value(parse(&d["value"])?));
                }
                Some("structure") => {
                    defs.insert(
                        name,
                        Def::Structure { params: d["params"].as_u64().unwrap_or(0) as usize, ctor: parse(&d["ctor"])? },
                    );
                }
                _ => {}
            }
        }
        let ty = parse(&j["type"])?;
        // telescope up to and including the last hypothesis; the rest is the conclusion
        let mut chain = Vec::new();
        let mut e = ty.clone();
        while let E::Pi(b, body, _) = &*e.clone() {
            chain.push((b.clone(), body.clone()));
            e = body.clone();
        }
        let last_hyp = chain.iter().rposition(|(b, _)| b.prop && b.bi != "inst");
        let cut = match last_hyp {
            Some(i) => i + 1,
            None => chain.iter().take_while(|(b, _)| b.bi == "inst" || !b.prop).count().min(
                // without hypotheses, keep leading binders as variables
                chain.len(),
            ),
        };
        let mut parts = Vec::new();
        for (b, _) in chain.iter().take(cut) {
            parts.push(if b.bi == "inst" {
                Part::Inst(b.clone())
            } else if b.prop {
                Part::Hyp(b.clone())
            } else {
                Part::Var(b.clone())
            });
        }
        let conclusion = if cut == 0 { ty } else { chain[cut - 1].1.clone() };
        Ok(Statement {
            name: j["name"].as_str().unwrap_or("").to_string(),
            parts,
            conclusion,
            ev: Ev { defs: Rc::new(defs) },
        })
    }

    pub fn vars(&self) -> impl Iterator<Item = &Binder> {
        self.parts.iter().filter_map(|p| if let Part::Var(b) = p { Some(b) } else { None })
    }
    pub fn hyps(&self) -> impl Iterator<Item = &Binder> {
        self.parts.iter().filter_map(|p| if let Part::Hyp(b) = p { Some(b) } else { None })
    }

    /// Top-level conjuncts of the conclusion, through `let`s, evaluated in `env`.
    pub fn conjuncts(&self, env: &Env) -> R<Vec<Tri>> {
        let mut out = Vec::new();
        self.split(&self.conclusion, env, &mut out)?;
        Ok(out)
    }
    fn split(&self, e: &Rc<E>, env: &Env, out: &mut Vec<Tri>) -> R<()> {
        match &**e {
            E::Let { id, value, body, .. } => {
                let v = self.ev.eval(value, env)?;
                self.split(body, &bind(env, id, v), out)
            }
            _ => {
                if let Some(("And", args)) = head(e) {
                    self.split(&args[0], env, out)?;
                    self.split(&args[1], env, out)
                } else {
                    out.push(self.ev.prop(e, env)?);
                    Ok(())
                }
            }
        }
    }
    pub fn conjunct_count(&self) -> usize {
        fn count(e: &Rc<E>) -> usize {
            match &**e {
                E::Let { body, .. } => count(body),
                _ => match head(e) {
                    Some(("And", args)) => count(&args[0]) + count(&args[1]),
                    _ => 1,
                },
            }
        }
        count(&self.conclusion)
    }
}

/// A defining hypothesis `v = e` or `∀ x, (guards →) v x = e`, for a variable `v`.
pub struct Definition {
    pub var_id: String,
    pub hyp_index: usize,
}

pub fn defining_hypotheses(st: &Statement) -> Vec<Definition> {
    let mut out = Vec::new();
    for (i, b) in st.hyps().enumerate() {
        let mut e = b.ty.clone();
        let mut bound = Vec::new();
        loop {
            match &*e.clone() {
                E::Pi(pb, body, _) => {
                    if !pb.prop {
                        bound.push(pb.id.clone());
                    }
                    e = body.clone();
                }
                _ => break,
            }
        }
        if let Some(("Eq", args)) = head(&e) {
            let lhs = &args[1];
            let target = match &**lhs {
                E::Fvar { id, .. } if bound.is_empty() => Some(id.clone()),
                E::App(f, xs) => match &**f {
                    E::Fvar { id, .. }
                        if xs.len() == bound.len()
                            && xs.iter().zip(&bound).all(|(x, b)| matches!(&**x, E::Fvar { id, .. } if id == b)) =>
                    {
                        Some(id.clone())
                    }
                    _ => None,
                },
                _ => None,
            };
            if let Some(var_id) = target {
                out.push(Definition { var_id, hyp_index: i });
            }
        }
    }
    out
}

/// Value of the variable defined by hypothesis `def`, as a function of its bound arguments.
pub fn defined_value(st: &Statement, def: &Definition, env: &Env) -> R<V> {
    let b = st.hyps().nth(def.hyp_index).ok_or("hypothesis")?;
    fn build(ev: &Ev, e: &Rc<E>, env: &Env) -> R<V> {
        match &**e {
            E::Pi(pb, body, _) if pb.prop => build(ev, body, &bind(env, &pb.id, V::Opaque)),
            E::Pi(pb, body, _) => {
                let (ev2, env2, id, body2) = (ev.clone(), env.clone(), pb.id.clone(), body.clone());
                Ok(V::Fun(Rc::new(move |x: V| build(&ev2, &body2, &bind(&env2, &id, x)))))
            }
            _ => match head(e) {
                Some(("Eq", args)) => ev.eval(&args[2], env),
                _ => Err("not a definition".into()),
            },
        }
    }
    build(&st.ev, &b.ty, env)
}
