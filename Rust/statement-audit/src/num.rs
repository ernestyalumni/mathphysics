//! Exact rationals (`qft_symbolics::exact::BigQ`) and rigorous rational intervals.
//!
//! Interval endpoints are exact rationals. Only `exp`, `log` and very long endpoints are rounded,
//! always outward, so every interval contains the true real value. Rounding is relative (a fixed
//! number of significant bits), which keeps tiny values such as `exp δ - 1` for `δ ≈ 2⁻²⁵⁶` narrow.
use num_bigint::{BigInt, Sign};
use num_traits::{One, Signed, Zero};
use qft_symbolics::computation::{Budget, Limits};
use qft_symbolics::exact::BigQ;
use std::cell::RefCell;
use std::cmp::Ordering;

thread_local! {
    static BUDGET: RefCell<Budget> = RefCell::new(Budget::new(Limits {
        max_operations: u64::MAX,
        max_terms: usize::MAX,
        max_degree: usize::MAX,
        max_coefficient_bits: u64::MAX,
        max_input_bytes: usize::MAX,
        max_word_length: usize::MAX,
    }));
}

fn with_budget<T>(f: impl FnOnce(&mut Budget) -> qft_symbolics::computation::Result<T>) -> T {
    BUDGET.with(|b| {
        let mut b = b.borrow_mut();
        // The counter only accumulates; reset it so a long run never reaches the limit.
        if b.used() > u64::MAX / 2 {
            *b = Budget::new(*b.limits());
        }
        f(&mut b).expect("exact arithmetic with unlimited budget")
    })
}

/// Significant bits kept when an endpoint is rounded.
pub const PRECISION: u64 = 1200;
/// Endpoints longer than this are rounded outward to `PRECISION` bits.
const ROUND_THRESHOLD: u64 = 2400;

pub fn q(n: i64, d: i64) -> BigQ {
    with_budget(|b| BigQ::new(BigInt::from(n), BigInt::from(d), b))
}
pub fn qi(n: i64) -> BigQ {
    BigQ::from(n)
}
pub fn qbig(n: BigInt, d: BigInt) -> BigQ {
    with_budget(|b| BigQ::new(n, d, b))
}
pub fn add(a: &BigQ, b: &BigQ) -> BigQ {
    with_budget(|bu| a.checked_add(b, bu))
}
pub fn sub(a: &BigQ, b: &BigQ) -> BigQ {
    with_budget(|bu| a.checked_sub(b, bu))
}
pub fn mul(a: &BigQ, b: &BigQ) -> BigQ {
    with_budget(|bu| a.checked_mul(b, bu))
}
/// Panics on a zero divisor; callers check first.
pub fn div(a: &BigQ, b: &BigQ) -> BigQ {
    with_budget(|bu| a.checked_div(b, bu))
}
pub fn neg(a: &BigQ) -> BigQ {
    with_budget(|bu| a.checked_neg(bu))
}
pub fn pow(a: &BigQ, n: usize) -> BigQ {
    with_budget(|bu| a.checked_pow(n, bu))
}
pub fn sign(a: &BigQ) -> i32 {
    match a.numerator().sign() {
        Sign::Minus => -1,
        Sign::NoSign => 0,
        Sign::Plus => 1,
    }
}
pub fn cmp(a: &BigQ, b: &BigQ) -> Ordering {
    sign(&sub(a, b)).cmp(&0)
}
pub fn abs(a: &BigQ) -> BigQ {
    if sign(a) < 0 {
        neg(a)
    } else {
        a.clone()
    }
}
pub fn min(a: &BigQ, b: &BigQ) -> BigQ {
    if cmp(a, b) == Ordering::Greater { b.clone() } else { a.clone() }
}
pub fn max(a: &BigQ, b: &BigQ) -> BigQ {
    if cmp(a, b) == Ordering::Less { b.clone() } else { a.clone() }
}
/// `floor(log2 |x|)` up to one, from bit lengths; `x` nonzero.
fn log2_estimate(x: &BigQ) -> i64 {
    x.numerator().bits() as i64 - x.denominator().bits() as i64
}
fn floor_div(n: &BigInt, d: &BigInt) -> BigInt {
    // d > 0
    let q = n / d;
    if n.is_negative() && !(n % d).is_zero() {
        q - BigInt::one()
    } else {
        q
    }
}
fn two_pow(k: u64) -> BigInt {
    BigInt::one() << k
}
/// Round `x` down (or up) to `bits` significant bits.
pub fn round_to(x: &BigQ, up: bool, bits: u64) -> BigQ {
    if x.is_zero() {
        return x.clone();
    }
    let e = log2_estimate(x);
    let s = bits as i64 - e; // scale by 2^s
    let (num, den) = if s >= 0 {
        (x.numerator() << (s as u64), x.denominator().clone())
    } else {
        (x.numerator().clone(), x.denominator() << ((-s) as u64))
    };
    let mut m = floor_div(&num, &den);
    if up && !(&m * &den == num) {
        m += BigInt::one();
    }
    if s >= 0 {
        qbig(m, two_pow(s as u64))
    } else {
        qbig(m << ((-s) as u64), BigInt::one())
    }
}
fn tidy(x: BigQ, up: bool) -> BigQ {
    if x.bits() > ROUND_THRESHOLD {
        round_to(&x, up, PRECISION)
    } else {
        x
    }
}

/// A closed interval `[lo, hi]` of reals with rational endpoints.
#[derive(Clone, Debug)]
pub struct Iv {
    pub lo: BigQ,
    pub hi: BigQ,
}

impl Iv {
    pub fn exact(x: BigQ) -> Self {
        Iv { lo: x.clone(), hi: x }
    }
    pub fn int(n: i64) -> Self {
        Self::exact(qi(n))
    }
    pub fn new(lo: BigQ, hi: BigQ) -> Self {
        Iv { lo: tidy(lo, false), hi: tidy(hi, true) }
    }
    pub fn is_exact(&self) -> bool {
        cmp(&self.lo, &self.hi) == Ordering::Equal
    }
    pub fn exact_value(&self) -> Option<&BigQ> {
        if self.is_exact() { Some(&self.lo) } else { None }
    }
    pub fn contains_zero(&self) -> bool {
        sign(&self.lo) <= 0 && sign(&self.hi) >= 0
    }
    pub fn is_exact_zero(&self) -> bool {
        self.lo.is_zero() && self.hi.is_zero()
    }
    pub fn add(&self, o: &Iv) -> Iv {
        Iv::new(add(&self.lo, &o.lo), add(&self.hi, &o.hi))
    }
    pub fn sub(&self, o: &Iv) -> Iv {
        Iv::new(sub(&self.lo, &o.hi), sub(&self.hi, &o.lo))
    }
    pub fn neg(&self) -> Iv {
        Iv { lo: neg(&self.hi), hi: neg(&self.lo) }
    }
    pub fn mul(&self, o: &Iv) -> Iv {
        if self.is_exact() && o.is_exact() {
            return Iv::exact_or_rounded(mul(&self.lo, &o.lo));
        }
        let c = [mul(&self.lo, &o.lo), mul(&self.lo, &o.hi), mul(&self.hi, &o.lo), mul(&self.hi, &o.hi)];
        let mut lo = c[0].clone();
        let mut hi = c[0].clone();
        for x in &c[1..] {
            lo = min(&lo, x);
            hi = max(&hi, x);
        }
        Iv::new(lo, hi)
    }
    /// An exact value, or a rigorous rounded enclosure of it if it is very long.
    pub fn exact_or_rounded(p: BigQ) -> Iv {
        if p.bits() > ROUND_THRESHOLD {
            Iv { lo: round_to(&p, false, PRECISION), hi: round_to(&p, true, PRECISION) }
        } else {
            Iv::exact(p)
        }
    }
    /// `None` if the divisor interval contains zero.
    pub fn div(&self, o: &Iv) -> Option<Iv> {
        if o.contains_zero() {
            return None;
        }
        let inv = if o.is_exact() {
            Iv::exact_or_rounded(div(&qi(1), &o.lo))
        } else {
            Iv::new(div(&qi(1), &o.hi), div(&qi(1), &o.lo))
        };
        Some(self.mul(&inv))
    }
    pub fn powi(&self, n: u64) -> Iv {
        if self.is_exact() {
            return Iv::exact_or_rounded(pow(&self.lo, n as usize));
        }
        let mut acc = Iv::int(1);
        for _ in 0..n {
            acc = acc.mul(self);
        }
        if n % 2 == 0 && self.contains_zero() {
            acc.lo = qi(0);
        }
        acc
    }
    /// Width relative to magnitude is at most 2^-k (exact intervals always qualify).
    pub fn narrow(&self, k: u64) -> bool {
        if self.is_exact() {
            return true;
        }
        let w = sub(&self.hi, &self.lo);
        let m = max(&abs(&self.lo), &abs(&self.hi));
        if m.is_zero() {
            return true;
        }
        cmp(&mul(&w, &qbig(two_pow(k), BigInt::one())), &m) != Ordering::Greater
    }
    pub fn to_f64(&self) -> f64 {
        let f = |x: &BigQ| {
            use num_traits::ToPrimitive;
            let s = x.numerator().bits().max(x.denominator().bits()).saturating_sub(900);
            let n = (x.numerator() >> s).to_f64().unwrap_or(f64::NAN);
            let d = (x.denominator() >> s).to_f64().unwrap_or(f64::NAN);
            n / d
        };
        (f(&self.lo) + f(&self.hi)) / 2.0
    }
}

/// Decision for `a < b` style comparisons: `Some(true)`, `Some(false)`, or `None` (undecided).
/// Intervals that overlap and are both narrow to 600 relative bits are treated as equal.
pub fn compare(a: &Iv, b: &Iv) -> Option<Ordering> {
    if cmp(&a.hi, &b.lo) == Ordering::Less {
        return Some(Ordering::Less);
    }
    if cmp(&a.lo, &b.hi) == Ordering::Greater {
        return Some(Ordering::Greater);
    }
    if a.narrow(600) && b.narrow(600) {
        return Some(Ordering::Equal);
    }
    None
}

// ---------- exp and log ----------

/// Lower and upper bounds of `exp(r)` for rational `0 <= r <= 1/2`.
fn exp_small(r: &BigQ) -> (BigQ, BigQ) {
    if r.is_zero() {
        return (qi(1), qi(1));
    }
    let target = PRECISION + 8;
    let mut lo_sum = qi(1);
    let mut hi_sum = qi(1);
    let mut term_lo = qi(1);
    let mut term_hi = qi(1);
    let rr = round_to(r, false, target);
    let rr_up = round_to(r, true, target);
    let r_small = cmp(&abs(r), &qi(1)) == Ordering::Less;
    let mut i: i64 = 1;
    loop {
        term_lo = round_to(&div(&mul(&term_lo, &rr), &qi(i)), false, target);
        term_hi = round_to(&div(&mul(&term_hi, &rr_up), &qi(i)), true, target);
        lo_sum = add(&lo_sum, &term_lo);
        hi_sum = add(&hi_sum, &term_hi);
        i += 1;
        // stop when the next term is below r * 2^-target (relative to the first-order part)
        let scale = if r_small { r.clone() } else { qi(1) };
        let bound = div(&scale, &qbig(two_pow(target), BigInt::one()));
        if cmp(&term_hi, &bound) == Ordering::Less || i > 4000 {
            break;
        }
    }
    // remainder <= 2 * next term <= 2 * term_hi * r
    let rem = mul(&qi(2), &mul(&term_hi, &rr_up));
    hi_sum = add(&hi_sum, &rem);
    (round_to(&lo_sum, false, PRECISION), round_to(&hi_sum, true, PRECISION))
}

thread_local! {
    static EXP_CACHE: RefCell<std::collections::HashMap<(String, String), (BigQ, BigQ)>> = RefCell::new(Default::default());
    static LOG_CACHE: RefCell<std::collections::HashMap<(String, String), (BigQ, BigQ)>> = RefCell::new(Default::default());
}
fn key(x: &BigQ) -> (String, String) {
    (x.numerator().to_str_radix(32), x.denominator().to_str_radix(32))
}
/// Rigorous bounds of `exp(x)` for rational `x` (memoised: the same arguments recur across instances).
pub fn exp_bounds(x: &BigQ) -> (BigQ, BigQ) {
    let k = key(x);
    if let Some(v) = EXP_CACHE.with(|c| c.borrow().get(&k).cloned()) {
        return v;
    }
    let v = exp_bounds_uncached(x);
    EXP_CACHE.with(|c| {
        let mut c = c.borrow_mut();
        if c.len() > 200_000 {
            c.clear();
        }
        c.insert(k, v.clone());
    });
    v
}
fn exp_bounds_uncached(x: &BigQ) -> (BigQ, BigQ) {
    if x.is_zero() {
        return (qi(1), qi(1));
    }
    if sign(x) < 0 {
        let (l, h) = exp_bounds_uncached(&neg(x));
        return (round_to(&div(&qi(1), &h), false, PRECISION), round_to(&div(&qi(1), &l), true, PRECISION));
    }
    // x > 0: reduce to r = x / 2^s <= 1/2
    let mut s: u64 = 0;
    let e = log2_estimate(x);
    if e >= -1 {
        s = (e + 2) as u64;
    }
    let r = div(x, &qbig(two_pow(s), BigInt::one()));
    let (mut lo, mut hi) = exp_small(&r);
    for _ in 0..s {
        lo = round_to(&mul(&lo, &lo), false, PRECISION);
        hi = round_to(&mul(&hi, &hi), true, PRECISION);
    }
    (lo, hi)
}

/// Bounds of `atanh(t)` for rational `0 <= t <= 1/3`.
fn atanh_small(t: &BigQ) -> (BigQ, BigQ) {
    if t.is_zero() {
        return (qi(0), qi(0));
    }
    let target = PRECISION + 8;
    let t_lo = round_to(t, false, target);
    let t_hi = round_to(t, true, target);
    let t2_lo = round_to(&mul(&t_lo, &t_lo), false, target);
    let t2_hi = round_to(&mul(&t_hi, &t_hi), true, target);
    let mut p_lo = t_lo.clone();
    let mut p_hi = t_hi.clone();
    let mut s_lo = t_lo.clone();
    let mut s_hi = t_hi.clone();
    let mut k: i64 = 1;
    let bound = div(t, &qbig(two_pow(target), BigInt::one()));
    loop {
        p_lo = round_to(&mul(&p_lo, &t2_lo), false, target);
        p_hi = round_to(&mul(&p_hi, &t2_hi), true, target);
        s_lo = add(&s_lo, &round_to(&div(&p_lo, &qi(2 * k + 1)), false, target));
        s_hi = add(&s_hi, &round_to(&div(&p_hi, &qi(2 * k + 1)), true, target));
        k += 1;
        if cmp(&p_hi, &bound) == Ordering::Less || k > 4000 {
            break;
        }
    }
    // tail <= p_hi * t^2 / (1 - t^2) <= (9/8) p_hi t^2
    let tail = mul(&q(9, 8), &mul(&p_hi, &t2_hi));
    (round_to(&s_lo, false, PRECISION), round_to(&add(&s_hi, &tail), true, PRECISION))
}

thread_local! {
    static LOG2: (BigQ, BigQ) = {
        let (l, h) = atanh_small(&q(1, 3));
        (mul(&qi(2), &l), mul(&qi(2), &h))
    };
}
pub fn log2_bounds() -> (BigQ, BigQ) {
    LOG2.with(|v| v.clone())
}

/// Rigorous bounds of `log(x)` for rational `x > 0` (memoised).
pub fn log_bounds(x: &BigQ) -> (BigQ, BigQ) {
    assert!(sign(x) > 0);
    let k = key(x);
    if let Some(v) = LOG_CACHE.with(|c| c.borrow().get(&k).cloned()) {
        return v;
    }
    let v = log_bounds_uncached(x);
    LOG_CACHE.with(|c| {
        let mut c = c.borrow_mut();
        if c.len() > 200_000 {
            c.clear();
        }
        c.insert(k, v.clone());
    });
    v
}
fn log_bounds_uncached(x: &BigQ) -> (BigQ, BigQ) {
    // x = 2^e * y with y in [1/2, 2]
    let mut e = log2_estimate(x);
    let mut y = if e >= 0 {
        div(x, &qbig(two_pow(e as u64), BigInt::one()))
    } else {
        mul(x, &qbig(two_pow((-e) as u64), BigInt::one()))
    };
    while cmp(&y, &qi(2)) == Ordering::Greater {
        y = div(&y, &qi(2));
        e += 1;
    }
    while cmp(&y, &q(1, 2)) == Ordering::Less {
        y = mul(&y, &qi(2));
        e -= 1;
    }
    let y = tidy(y, false);
    // log y = 2 atanh((y-1)/(y+1)); |t| <= 1/3
    let t = div(&sub(&y, &qi(1)), &add(&y, &qi(1)));
    let (a_lo, a_hi) = if sign(&t) >= 0 {
        atanh_small(&t)
    } else {
        let (l, h) = atanh_small(&neg(&t));
        (neg(&h), neg(&l))
    };
    let (l2_lo, l2_hi) = log2_bounds();
    let ee = qi(e);
    let (e_lo, e_hi) = if e >= 0 { (mul(&ee, &l2_lo), mul(&ee, &l2_hi)) } else { (mul(&ee, &l2_hi), mul(&ee, &l2_lo)) };
    (
        round_to(&add(&e_lo, &mul(&qi(2), &a_lo)), false, PRECISION),
        round_to(&add(&e_hi, &mul(&qi(2), &a_hi)), true, PRECISION),
    )
}

pub fn exp_iv(x: &Iv) -> Iv {
    let (lo, _) = exp_bounds(&x.lo);
    let (_, hi) = exp_bounds(&x.hi);
    Iv { lo, hi }
}

/// `log` of an interval that lies strictly on one side of zero, with Lean's `log (-x) = log x`.
/// `None` if the interval contains zero (and is not exactly zero).
pub fn log_iv(x: &Iv, lean: bool) -> Option<Iv> {
    if x.is_exact_zero() {
        return if lean { Some(Iv::int(0)) } else { None };
    }
    if sign(&x.lo) > 0 {
        let (lo, _) = log_bounds(&x.lo);
        let (_, hi) = log_bounds(&x.hi);
        return Some(Iv { lo, hi });
    }
    if lean && sign(&x.hi) < 0 {
        return log_iv(&x.neg(), lean);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exp_log_roundtrip_and_constants() {
        let (l, h) = exp_bounds(&qi(1));
        assert!(Iv { lo: l, hi: h }.to_f64() - std::f64::consts::E < 1e-15);
        let (l, h) = log2_bounds();
        assert!((Iv { lo: l, hi: h }.to_f64() - std::f64::consts::LN_2).abs() < 1e-15);
        let x = Iv::exact(q(7, 3));
        let y = exp_iv(&log_iv(&x, false).unwrap());
        assert_eq!(compare(&x, &y), Some(Ordering::Equal));
        assert!(y.narrow(1000));
        // tiny arguments keep relative precision: exp(d) - 1 for d = 2^-700
        let d = qbig(BigInt::one(), two_pow(700));
        let e1 = exp_iv(&Iv::exact(d.clone())).sub(&Iv::int(1));
        assert!(e1.narrow(300));
        // exp(d) - 1 - d ~ d^2/2 is far below the precision here: undecided or equal, never wrong
        assert_ne!(compare(&e1, &Iv::exact(d.clone())), Some(Ordering::Less));
        assert!(cmp(&e1.lo, &d) != Ordering::Less);
        let d = q(1, 1 << 40);
        let e1 = exp_iv(&Iv::exact(d.clone())).sub(&Iv::int(1));
        assert_eq!(compare(&e1, &Iv::exact(d)), Some(Ordering::Greater));
        // negative arguments and large ones
        let (l, h) = exp_bounds(&qi(-40));
        assert!((Iv { lo: l, hi: h }.to_f64() / (-40f64).exp() - 1.0).abs() < 1e-12);
        let (l, h) = log_bounds(&qi(1_000_000));
        assert!((Iv { lo: l, hi: h }.to_f64() - 1e6f64.ln()).abs() < 1e-12);
    }
}
