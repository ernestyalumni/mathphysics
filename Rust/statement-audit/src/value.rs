//! Values shared by the Lean evaluator and the informal predicates.
//!
//! Two semantics: `Sem::Lean` is Lean's total arithmetic (`x / 0 = 0`, `log 0 = 0`,
//! `log (-x) = log x`, truncated `ℕ` subtraction); `Sem::Math` is ordinary partial arithmetic, where
//! a division by zero or the log of a nonpositive number is undecided.
use crate::num::{self, compare, Iv};
use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive, Zero};
use qft_symbolics::exact::BigQ;
use std::cmp::Ordering;
use std::fmt;
use std::rc::Rc;

pub type R<T> = Result<T, String>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tri {
    T,
    F,
    U,
}
impl Tri {
    pub fn from_bool(b: bool) -> Tri {
        if b { Tri::T } else { Tri::F }
    }
    pub fn and(self, o: Tri) -> Tri {
        match (self, o) {
            (Tri::F, _) | (_, Tri::F) => Tri::F,
            (Tri::T, Tri::T) => Tri::T,
            _ => Tri::U,
        }
    }
    pub fn or(self, o: Tri) -> Tri {
        match (self, o) {
            (Tri::T, _) | (_, Tri::T) => Tri::T,
            (Tri::F, Tri::F) => Tri::F,
            _ => Tri::U,
        }
    }
    pub fn not(self) -> Tri {
        match self {
            Tri::T => Tri::F,
            Tri::F => Tri::T,
            Tri::U => Tri::U,
        }
    }
    pub fn implies(self, o: Tri) -> Tri {
        self.not().or(o)
    }
    pub fn iff(self, o: Tri) -> Tri {
        match (self, o) {
            (Tri::U, _) | (_, Tri::U) => Tri::U,
            (a, b) => Tri::from_bool(a == b),
        }
    }
    pub fn decided(self) -> bool {
        self != Tri::U
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sem {
    Lean,
    Math,
}

pub type Func = Rc<dyn Fn(V) -> R<V>>;

#[derive(Clone)]
pub enum Ty {
    Nat,
    Int,
    Real,
    Complex,
    Prop,
    Type,
    Fin(usize),
    Matrix(Box<Ty>, Box<Ty>, Box<Ty>),
    Fun(Box<Ty>, Box<Ty>),
    Subtype(Box<Ty>, Func),
    Other(String),
}
impl fmt::Debug for Ty {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Ty::Nat => write!(f, "ℕ"),
            Ty::Int => write!(f, "ℤ"),
            Ty::Real => write!(f, "ℝ"),
            Ty::Complex => write!(f, "ℂ"),
            Ty::Prop => write!(f, "Prop"),
            Ty::Type => write!(f, "Type"),
            Ty::Fin(n) => write!(f, "Fin {n}"),
            Ty::Matrix(a, b, c) => write!(f, "Matrix ({a:?}) ({b:?}) {c:?}"),
            Ty::Fun(a, b) => write!(f, "{a:?} → {b:?}"),
            Ty::Subtype(a, _) => write!(f, "{{x : {a:?} // _}}"),
            Ty::Other(s) => write!(f, "{s}"),
        }
    }
}
impl Ty {
    pub fn size(&self) -> Option<usize> {
        match self {
            Ty::Fin(n) => Some(*n),
            _ => None,
        }
    }
}

#[derive(Clone)]
pub struct Mat {
    pub rows: usize,
    pub cols: usize,
    pub e: Vec<V>,
}
impl Mat {
    pub fn get(&self, i: usize, j: usize) -> &V {
        &self.e[i * self.cols + j]
    }
}

#[derive(Clone)]
pub enum V {
    Nat(BigInt),
    Int(BigInt),
    Real(Iv),
    /// Exact Gaussian rational.
    Cx(BigQ, BigQ),
    Prop(Tri),
    Mat(Rc<Mat>),
    /// A finitely supported sequence indexed from 0 (zero beyond the vector).
    Seq(Rc<Vec<V>>),
    Fun(Func),
    Ty(Ty),
    Opaque,
    Undecided,
}
impl fmt::Debug for V {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            V::Nat(n) => write!(f, "{n}"),
            V::Int(n) => write!(f, "{n}"),
            V::Real(x) => {
                if let Some(v) = x.exact_value() {
                    write!(f, "{v}")
                } else {
                    write!(f, "~{:.6e}", x.to_f64())
                }
            }
            V::Cx(a, b) => write!(f, "({a} + {b}i)"),
            V::Prop(t) => write!(f, "{t:?}"),
            V::Mat(m) => {
                write!(f, "[")?;
                for i in 0..m.rows {
                    if i > 0 {
                        write!(f, "; ")?;
                    }
                    for j in 0..m.cols {
                        if j > 0 {
                            write!(f, " ")?;
                        }
                        write!(f, "{:?}", m.get(i, j))?;
                    }
                }
                write!(f, "]")
            }
            V::Seq(s) => write!(f, "seq{:?}", s.as_slice()),
            V::Fun(_) => write!(f, "<fun>"),
            V::Ty(t) => write!(f, "{t:?}"),
            V::Opaque => write!(f, "_"),
            V::Undecided => write!(f, "?"),
        }
    }
}

pub fn real(q: BigQ) -> V {
    V::Real(Iv::exact(q))
}
pub fn real_int(n: i64) -> V {
    V::Real(Iv::int(n))
}
pub fn nat(n: u64) -> V {
    V::Nat(BigInt::from(n))
}

impl V {
    pub fn as_tri(&self) -> Tri {
        match self {
            V::Prop(t) => *t,
            _ => Tri::U, // includes opaque constants: undecided, never true or false
        }
    }
    pub fn as_usize(&self) -> Option<usize> {
        match self {
            V::Nat(n) | V::Int(n) => n.to_usize(),
            V::Real(x) => {
                let v = x.exact_value()?;
                if v.denominator().is_one() { v.numerator().to_usize() } else { None }
            }
            _ => None,
        }
    }
    pub fn as_ty(&self) -> Option<&Ty> {
        match self {
            V::Ty(t) => Some(t),
            _ => None,
        }
    }
    pub fn is_undecided(&self) -> bool {
        matches!(self, V::Undecided)
    }
}

/// Coerce a number to an interval, if it is real-valued.
fn to_iv(v: &V) -> Option<Iv> {
    match v {
        V::Nat(n) | V::Int(n) => Some(Iv::exact(num::qbig(n.clone(), BigInt::one()))),
        V::Real(x) => Some(x.clone()),
        V::Cx(a, b) if b.is_zero() => Some(Iv::exact(a.clone())),
        _ => None,
    }
}
fn to_cx(v: &V) -> Option<(BigQ, BigQ)> {
    match v {
        V::Cx(a, b) => Some((a.clone(), b.clone())),
        V::Nat(n) | V::Int(n) => Some((num::qbig(n.clone(), BigInt::one()), num::qi(0))),
        V::Real(x) => x.exact_value().map(|q| (q.clone(), num::qi(0))),
        _ => None,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Add,
    Sub,
    Mul,
    Div,
}

/// Apply a function-like value to an argument.
pub fn apply(f: &V, x: V) -> R<V> {
    match f {
        V::Fun(g) => g(x),
        V::Seq(s) => {
            if x.is_undecided() {
                return Ok(V::Undecided);
            }
            let i = x.as_usize().ok_or("sequence index is not a natural number")?;
            Ok(s.get(i).cloned().unwrap_or_else(|| zero_like(s.first())))
        }
        V::Mat(m) => {
            let i = x.as_usize().ok_or("matrix row index")?;
            if i >= m.rows {
                return Err("matrix row index out of range".into());
            }
            let row: Vec<V> = (0..m.cols).map(|j| m.get(i, j).clone()).collect();
            Ok(V::Seq(Rc::new(row)))
        }
        V::Undecided => Ok(V::Undecided),
        V::Opaque => Err("cannot apply an opaque constant".into()),
        other => Err(format!("cannot apply {other:?}")),
    }
}
fn zero_like(v: Option<&V>) -> V {
    match v {
        Some(V::Nat(_)) => nat(0),
        Some(V::Cx(..)) => V::Cx(num::qi(0), num::qi(0)),
        _ => real_int(0),
    }
}

pub fn arith(op: Op, a: &V, b: &V, sem: Sem) -> R<V> {
    use V::*;
    if a.is_undecided() || b.is_undecided() {
        return Ok(Undecided);
    }
    match (a, b) {
        (Nat(x), Nat(y)) => Ok(match op {
            Op::Add => Nat(x + y),
            Op::Sub => {
                if sem == Sem::Lean {
                    if x >= y { Nat(x - y) } else { Nat(BigInt::zero()) }
                } else {
                    Int(x - y)
                }
            }
            Op::Mul => Nat(x * y),
            Op::Div => {
                if y.is_zero() {
                    if sem == Sem::Lean { Nat(BigInt::zero()) } else { Undecided }
                } else {
                    Nat(x / y)
                }
            }
        }),
        (Int(x), Int(y)) | (Int(x), Nat(y)) | (Nat(x), Int(y)) if op != Op::Div => Ok(Int(match op {
            Op::Add => x + y,
            Op::Sub => x - y,
            _ => x * y,
        })),
        (Mat(x), Mat(y)) => mat_arith(op, x, y, sem),
        (Mat(x), s) | (s, Mat(x)) if op == Op::Mul => {
            let e = x.e.iter().map(|v| arith(Op::Mul, s, v, sem)).collect::<R<Vec<_>>>()?;
            Ok(Mat(Rc::new(self::Mat { rows: x.rows, cols: x.cols, e })))
        }
        (Fun(_), _) | (_, Fun(_)) | (Seq(_), _) | (_, Seq(_)) => {
            // pointwise operations on functions (Pi instances, linear maps)
            let (f, g) = (a.clone(), b.clone());
            let fun_a = matches!(f, Fun(_) | Seq(_));
            let fun_b = matches!(g, Fun(_) | Seq(_));
            Ok(Fun(Rc::new(move |x: V| {
                let fa = if fun_a { apply(&f, x.clone())? } else { f.clone() };
                let gb = if fun_b { apply(&g, x)? } else { g.clone() };
                arith(op, &fa, &gb, sem)
            })))
        }
        (Cx(..), _) | (_, Cx(..)) => {
            let (Some((ar, ai)), Some((br, bi))) = (to_cx(a), to_cx(b)) else {
                return Ok(Undecided);
            };
            let (r, i) = match op {
                Op::Add => (num::add(&ar, &br), num::add(&ai, &bi)),
                Op::Sub => (num::sub(&ar, &br), num::sub(&ai, &bi)),
                Op::Mul => (
                    num::sub(&num::mul(&ar, &br), &num::mul(&ai, &bi)),
                    num::add(&num::mul(&ar, &bi), &num::mul(&ai, &br)),
                ),
                Op::Div => {
                    let d = num::add(&num::mul(&br, &br), &num::mul(&bi, &bi));
                    if d.is_zero() {
                        return Ok(if sem == Sem::Lean { Cx(num::qi(0), num::qi(0)) } else { Undecided });
                    }
                    let nr = num::add(&num::mul(&ar, &br), &num::mul(&ai, &bi));
                    let ni = num::sub(&num::mul(&ai, &br), &num::mul(&ar, &bi));
                    (num::div(&nr, &d), num::div(&ni, &d))
                }
            };
            Ok(Cx(r, i))
        }
        _ => {
            let (Some(x), Some(y)) = (to_iv(a), to_iv(b)) else {
                return Err(format!("unsupported arithmetic {op:?} on {a:?} and {b:?}"));
            };
            Ok(match op {
                Op::Add => Real(x.add(&y)),
                Op::Sub => Real(x.sub(&y)),
                Op::Mul => Real(x.mul(&y)),
                Op::Div => {
                    if y.is_exact_zero() {
                        if sem == Sem::Lean { real_int(0) } else { Undecided }
                    } else {
                        match x.div(&y) {
                            Some(z) => Real(z),
                            None => Undecided,
                        }
                    }
                }
            })
        }
    }
}

fn mat_arith(op: Op, x: &Mat, y: &Mat, sem: Sem) -> R<V> {
    match op {
        Op::Add | Op::Sub => {
            if x.rows != y.rows || x.cols != y.cols {
                return Err("matrix shapes differ".into());
            }
            let e = x.e.iter().zip(&y.e).map(|(a, b)| arith(op, a, b, sem)).collect::<R<Vec<_>>>()?;
            Ok(V::Mat(Rc::new(Mat { rows: x.rows, cols: x.cols, e })))
        }
        Op::Mul => {
            if x.cols != y.rows {
                return Err("matrix shapes do not compose".into());
            }
            let mut e = Vec::with_capacity(x.rows * y.cols);
            for i in 0..x.rows {
                for j in 0..y.cols {
                    let mut acc = zero_like(x.e.first());
                    for k in 0..x.cols {
                        acc = arith(Op::Add, &acc, &arith(Op::Mul, x.get(i, k), y.get(k, j), sem)?, sem)?;
                    }
                    e.push(acc);
                }
            }
            Ok(V::Mat(Rc::new(Mat { rows: x.rows, cols: y.cols, e })))
        }
        Op::Div => Err("matrix division".into()),
    }
}

pub fn negate(a: &V, sem: Sem) -> R<V> {
    match a {
        V::Nat(_) => Ok(if sem == Sem::Lean { nat(0) } else { V::Undecided }),
        V::Int(n) => Ok(V::Int(-n)),
        V::Real(x) => Ok(V::Real(x.neg())),
        V::Cx(r, i) => Ok(V::Cx(num::neg(r), num::neg(i))),
        V::Undecided => Ok(V::Undecided),
        other => arith(Op::Sub, &zero_of(other), other, sem),
    }
}
fn zero_of(v: &V) -> V {
    match v {
        V::Mat(m) => V::Mat(Rc::new(Mat { rows: m.rows, cols: m.cols, e: vec![zero_like(m.e.first()); m.e.len()] })),
        _ => real_int(0),
    }
}

pub fn inverse(a: &V, sem: Sem) -> R<V> {
    arith(Op::Div, &real_int(1), a, sem)
}

/// `a ^ b`. Natural exponents are exact; real exponents use `exp (b log a)` for `a > 0`.
pub fn power(a: &V, b: &V, sem: Sem) -> R<V> {
    if a.is_undecided() || b.is_undecided() {
        return Ok(V::Undecided);
    }
    if let V::Nat(k) = b {
        let k = k.to_u64().ok_or("exponent too large")?;
        if k > 4096 {
            return Ok(V::Undecided);
        }
        return match a {
            V::Nat(x) => Ok(V::Nat(num_traits::pow::pow(x.clone(), k as usize))),
            V::Int(x) => Ok(V::Int(num_traits::pow::pow(x.clone(), k as usize))),
            V::Real(x) => Ok(V::Real(x.powi(k))),
            V::Cx(..) | V::Mat(_) => {
                let mut acc = match a {
                    V::Mat(m) => identity(m.rows, m.e.first()),
                    _ => V::Cx(num::qi(1), num::qi(0)),
                };
                for _ in 0..k {
                    acc = arith(Op::Mul, &acc, a, sem)?;
                }
                Ok(acc)
            }
            _ => Err(format!("unsupported power base {a:?}")),
        };
    }
    let (Some(x), Some(y)) = (to_iv(a), to_iv(b)) else {
        return Err(format!("unsupported power {a:?} ^ {b:?}"));
    };
    if x.is_exact_zero() {
        return Ok(if y.is_exact_zero() {
            real_int(1)
        } else if !y.contains_zero() {
            real_int(0)
        } else {
            V::Undecided
        });
    }
    if num::sign(&x.lo) > 0 {
        let l = num::log_iv(&x, false).unwrap();
        return Ok(V::Real(num::exp_iv(&l.mul(&y))));
    }
    Ok(V::Undecided)
}

pub fn identity(n: usize, like: Option<&V>) -> V {
    let one = match like {
        Some(V::Cx(..)) => V::Cx(num::qi(1), num::qi(0)),
        _ => real_int(1),
    };
    let zero = zero_like(like);
    let e = (0..n * n).map(|k| if k / n == k % n { one.clone() } else { zero.clone() }).collect();
    V::Mat(Rc::new(Mat { rows: n, cols: n, e }))
}

pub fn log(a: &V, sem: Sem) -> V {
    match to_iv(a) {
        Some(x) => match num::log_iv(&x, sem == Sem::Lean) {
            Some(y) => V::Real(y),
            None => V::Undecided,
        },
        None => V::Undecided,
    }
}
pub fn exp(a: &V) -> V {
    match to_iv(a) {
        Some(x) => V::Real(num::exp_iv(&x)),
        None => V::Undecided,
    }
}

/// Order comparison of two real-valued numbers.
pub fn order(a: &V, b: &V) -> Option<Ordering> {
    match (a, b) {
        (V::Nat(x), V::Nat(y)) | (V::Int(x), V::Int(y)) | (V::Nat(x), V::Int(y)) | (V::Int(x), V::Nat(y)) => {
            Some(x.cmp(y))
        }
        _ => compare(&to_iv(a)?, &to_iv(b)?),
    }
}
pub fn le(a: &V, b: &V) -> Tri {
    match order(a, b) {
        Some(o) => Tri::from_bool(o != Ordering::Greater),
        None => Tri::U,
    }
}
pub fn lt(a: &V, b: &V) -> Tri {
    match order(a, b) {
        Some(o) => Tri::from_bool(o == Ordering::Less),
        None => Tri::U,
    }
}
pub fn eq(a: &V, b: &V) -> Tri {
    match (a, b) {
        (V::Undecided, _) | (_, V::Undecided) => Tri::U,
        (V::Prop(x), V::Prop(y)) => x.iff(*y),
        (V::Mat(x), V::Mat(y)) => {
            if x.rows != y.rows || x.cols != y.cols {
                return Tri::F;
            }
            x.e.iter().zip(&y.e).fold(Tri::T, |acc, (p, q)| acc.and(eq(p, q)))
        }
        (V::Cx(..), _) | (_, V::Cx(..)) => match (to_cx(a), to_cx(b)) {
            (Some((ar, ai)), Some((br, bi))) => {
                Tri::from_bool(num::cmp(&ar, &br) == Ordering::Equal && num::cmp(&ai, &bi) == Ordering::Equal)
            }
            _ => Tri::U,
        },
        (V::Seq(x), V::Seq(y)) => {
            let n = x.len().max(y.len());
            (0..n).fold(Tri::T, |acc, i| {
                let p = x.get(i).cloned().unwrap_or_else(|| zero_like(x.first()));
                let q = y.get(i).cloned().unwrap_or_else(|| zero_like(y.first()));
                acc.and(eq(&p, &q))
            })
        }
        _ => match order(a, b) {
            Some(o) => Tri::from_bool(o == Ordering::Equal),
            None => Tri::U,
        },
    }
}

/// Trace of a square matrix.
pub fn trace(m: &V, sem: Sem) -> R<V> {
    match m {
        V::Mat(m) => {
            let mut acc = zero_like(m.e.first());
            for i in 0..m.rows.min(m.cols) {
                acc = arith(Op::Add, &acc, m.get(i, i), sem)?;
            }
            Ok(acc)
        }
        V::Undecided => Ok(V::Undecided),
        other => Err(format!("trace of {other:?}")),
    }
}

/// Trace of a linear map on `n × n` matrices, from its values on the matrix units.
pub fn linear_map_trace(f: &V, n: usize, complex: bool, sem: Sem) -> R<V> {
    let (zero, one) = if complex {
        (V::Cx(num::qi(0), num::qi(0)), V::Cx(num::qi(1), num::qi(0)))
    } else {
        (real_int(0), real_int(1))
    };
    let mut acc = zero.clone();
    for i in 0..n {
        for j in 0..n {
            let mut e = vec![zero.clone(); n * n];
            e[i * n + j] = one.clone();
            let image = apply(f, V::Mat(Rc::new(Mat { rows: n, cols: n, e })))?;
            match image {
                V::Mat(m) => acc = arith(Op::Add, &acc, m.get(i, j), sem)?,
                V::Undecided => return Ok(V::Undecided),
                other => return Err(format!("linear map returned {other:?}")),
            }
        }
    }
    Ok(acc)
}

/// Conjugate transpose (transpose for real matrices).
pub fn conj_transpose(m: &V) -> R<V> {
    match m {
        V::Mat(m) => {
            let mut e = Vec::with_capacity(m.e.len());
            for i in 0..m.cols {
                for j in 0..m.rows {
                    e.push(match m.get(j, i) {
                        V::Cx(r, im) => V::Cx(r.clone(), num::neg(im)),
                        other => other.clone(),
                    });
                }
            }
            Ok(V::Mat(Rc::new(Mat { rows: m.cols, cols: m.rows, e })))
        }
        other => Err(format!("conjugate transpose of {other:?}")),
    }
}

pub fn real_part(v: &V) -> V {
    match v {
        V::Cx(r, _) => real(r.clone()),
        V::Real(_) => v.clone(),
        _ => V::Undecided,
    }
}

impl From<i64> for V {
    fn from(n: i64) -> V {
        real_int(n)
    }
}

pub fn big_abs(n: &BigInt) -> BigInt {
    n.abs()
}
