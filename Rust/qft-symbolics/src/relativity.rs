//! Time-first Lorentz operations with metric (-+++), positive boost off-diagonals.
use crate::{
    checked::{self, Matrix},
    coefficient::Coefficient,
    computation::{Budget, Error, Resource, Result},
    exact::{BigC, BigQ},
    polynomial::BigPoly,
};
fn domain(message: &'static str) -> Error {
    Error::UnsupportedDomain { message }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(domain("nonfinite input or numerical overflow"))
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rapidity(f64);
impl Rapidity {
    pub fn new(eta: f64) -> Result<Self> {
        Ok(Self(finite(eta)?))
    }
    pub fn from_beta(beta: f64) -> Result<Self> {
        finite(beta)?;
        if beta.abs() >= 1.0 {
            return Err(domain("beta must satisfy |beta|<1"));
        }
        Self::new(beta.atanh())
    }
    pub fn eta(self) -> f64 {
        self.0
    }
    /// Errors if tanh rounds to +/-1; no silent clamp or saturated velocity.
    pub fn beta(self) -> Result<f64> {
        let v = self.0.tanh();
        if v.abs() >= 1.0 {
            Err(domain("beta rounded to light speed"))
        } else {
            finite(v)
        }
    }
    pub fn inverse(self) -> Self {
        Self(-self.0)
    }
    pub fn collinear_add(self, r: Self) -> Result<Self> {
        Self::new(self.0 + r.0)
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}
impl Axis {
    fn index(self) -> usize {
        match self {
            Self::X => 1,
            Self::Y => 2,
            Self::Z => 3,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnitAxis([f64; 3]);
impl UnitAxis {
    /// Absolute tolerance on Euclidean norm; at most 1e-6. Does not normalize.
    pub fn new(n: [f64; 3], tolerance: f64) -> Result<Self> {
        finite(tolerance)?;
        if !(0.0..=1e-6).contains(&tolerance) {
            return Err(domain("axis tolerance outside 0..=1e-6"));
        }
        for x in n {
            finite(x)?;
        }
        let norm = finite(n[0].hypot(n[1]).hypot(n[2]))?;
        if (norm - 1.0).abs() > tolerance {
            return Err(domain("axis is not unit length"));
        }
        Ok(Self(n))
    }
    pub fn axial(axis: Axis) -> Self {
        let mut n = [0.0; 3];
        n[axis.index() - 1] = 1.0;
        Self(n)
    }
    pub fn components(self) -> [f64; 3] {
        self.0
    }
}
pub fn boost_f64(eta: Rapidity, axis: UnitAxis) -> Result<Matrix<f64, 4>> {
    let c = finite(eta.0.cosh())?;
    let s = finite(eta.0.sinh())?;
    let half = finite((eta.0 / 2.0).sinh())?;
    let cm1 = finite(half * half * 2.0)?;
    let mut m = [[0.0; 4]; 4];
    m[0][0] = c;
    for i in 0..3 {
        m[0][i + 1] = finite(axis.0[i] * s)?;
        m[i + 1][0] = m[0][i + 1];
        for j in 0..3 {
            m[i + 1][j + 1] = finite(if i == j { 1.0 } else { 0.0 } + axis.0[i] * axis.0[j] * cm1)?;
        }
    }
    Ok(m)
}
pub fn apply_f64(m: &Matrix<f64, 4>, v: &[f64; 4]) -> Result<[f64; 4]> {
    let mut out = [0.0; 4];
    for i in 0..4 {
        for (j, x) in v.iter().enumerate() {
            finite(*x)?;
            finite(m[i][j])?;
            out[i] = finite(out[i] + finite(m[i][j] * x)?)?;
        }
    }
    Ok(out)
}
pub fn dot_f64(v: &[f64; 4], w: &[f64; 4]) -> Result<f64> {
    let mut s = 0.0;
    for i in 0..4 {
        finite(v[i])?;
        finite(w[i])?;
        s = finite(s + finite(v[i] * w[i])? * if i == 0 { -1.0 } else { 1.0 })?;
    }
    Ok(s)
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    pub absolute: f64,
    pub relative: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Residual {
    pub absolute: f64,
    pub scaled: f64,
    pub scale: f64,
    pub within_tolerance: bool,
}
/// max |L^T g L-g|; scale=max(1,max_ij sum_k |Lki Lkj|).
/// A tolerance pass is a numerical observation, not an exact Lorentz certificate.
pub fn metric_residual_f64(m: &Matrix<f64, 4>, tol: Tolerance) -> Result<Residual> {
    finite(tol.absolute)?;
    finite(tol.relative)?;
    if tol.absolute < 0.0 || tol.relative < 0.0 {
        return Err(domain("negative tolerance"));
    }
    for row in m {
        for x in row {
            finite(*x)?;
        }
    }
    let mut absolute = 0.0_f64;
    let mut scale = 1.0_f64;
    for i in 0..4 {
        for j in 0..4 {
            let mut sum = 0.0;
            let mut magnitude = 0.0;
            for (k, row) in m.iter().enumerate() {
                let p = finite(row[i] * row[j])?;
                sum = finite(sum + if k == 0 { -p } else { p })?;
                magnitude = finite(magnitude + p.abs())?;
            }
            let target = if i == j {
                if i == 0 {
                    -1.0
                } else {
                    1.0
                }
            } else {
                0.0
            };
            absolute = absolute.max(finite(sum - target)?.abs());
            scale = scale.max(magnitude);
        }
    }
    let threshold = finite(tol.absolute + finite(tol.relative * scale)?)?;
    Ok(Residual {
        absolute,
        scaled: absolute / scale,
        scale,
        within_tolerance: absolute <= threshold,
    })
}
fn validate<A: Coefficient>(v: &[A; 4], b: &mut Budget) -> Result<()> {
    b.charge(4)?;
    b.check(Resource::Terms, 4)?;
    for x in v {
        b.check(Resource::CoefficientBits, x.bits())?;
    }
    Ok(())
}
pub fn metric<A: Coefficient>(b: &mut Budget) -> Result<Matrix<A, 4>> {
    let mut m = checked::identity(b)?;
    m[0][0] = A::from_i64(1).checked_neg(b)?;
    Ok(m)
}
pub fn lower<A: Coefficient>(v: &[A; 4], b: &mut Budget) -> Result<[A; 4]> {
    validate(v, b)?;
    let t = v[0].checked_neg(b)?;
    Ok([t, v[1].clone(), v[2].clone(), v[3].clone()])
}
pub fn raise<A: Coefficient>(v: &[A; 4], b: &mut Budget) -> Result<[A; 4]> {
    lower(v, b)
}
pub fn dot<A: Coefficient>(v: &[A; 4], w: &[A; 4], b: &mut Budget) -> Result<A> {
    validate(v, b)?;
    validate(w, b)?;
    let mut s = v[0].checked_mul(&w[0], b)?.checked_neg(b)?;
    for i in 1..4 {
        s = s.checked_add(&v[i].checked_mul(&w[i], b)?, b)?;
    }
    Ok(s)
}
pub fn apply<A: Coefficient>(m: &Matrix<A, 4>, v: &[A; 4], b: &mut Budget) -> Result<[A; 4]> {
    validate(v, b)?;
    for row in m {
        validate(row, b)?;
    }
    let mut out = std::array::from_fn(|_| A::from_i64(0));
    for i in 0..4 {
        for (j, x) in v.iter().enumerate() {
            out[i] = out[i].checked_add(&m[i][j].checked_mul(x, b)?, b)?;
        }
    }
    Ok(out)
}
/// Exact rational c,s parameterization, not a rational representation of eta.
pub fn exact_half_velocity_boost(u: &BigQ, axis: Axis, b: &mut Budget) -> Result<Matrix<BigC, 4>> {
    b.charge(1)?;
    b.check(Resource::CoefficientBits, u.bits())?;
    if u.numerator().magnitude() >= u.denominator().magnitude() {
        return Err(domain("half velocity must satisfy |u|<1"));
    }
    let one = BigQ::from(1);
    let u2 = u.checked_mul(u, b)?;
    let d = one.checked_sub(&u2, b)?;
    let c = one.checked_add(&u2, b)?.checked_div(&d, b)?;
    let s = BigQ::from(2).checked_mul(u, b)?.checked_div(&d, b)?;
    let mut m = checked::identity::<BigC, 4>(b)?;
    let a = axis.index();
    let c = BigC::new(c, 0.into());
    let s = BigC::new(s, 0.into());
    m[0][0] = c.clone();
    m[a][a] = c;
    m[0][a] = s.clone();
    m[a][0] = s;
    Ok(m)
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SymbolicBoost {
    pub matrix: Matrix<BigPoly, 4>,
    /// Required zero polynomial, retained without quotient reduction.
    pub condition: BigPoly,
}
pub fn symbolic_axial_boost(
    c_label: &str,
    s_label: &str,
    axis: Axis,
    b: &mut Budget,
) -> Result<SymbolicBoost> {
    if c_label.is_empty() || s_label.is_empty() || c_label == s_label {
        return Err(domain("boost labels must be nonempty and distinct"));
    }
    b.charge(16)?;
    b.check(Resource::Terms, 16)?;
    let c = BigPoly::try_var(c_label, b)?;
    let s = BigPoly::try_var(s_label, b)?;
    let one = BigPoly::try_constant(1, b)?;
    let condition = c
        .try_mul(&c, b)?
        .try_sub(&s.try_mul(&s, b)?, b)?
        .try_sub(&one, b)?;
    let mut matrix = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i == j {
                one.clone()
            } else {
                BigPoly::default()
            }
        })
    });
    let a = axis.index();
    matrix[0][0] = c.clone();
    matrix[a][a] = c;
    matrix[0][a] = s.clone();
    matrix[a][0] = s;
    Ok(SymbolicBoost { matrix, condition })
}
