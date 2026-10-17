//! Three complementary exact encodings:
//! - Lie elements with Poly coordinates: abstract linear combinations/brackets;
//! - one-angle Trig quotient for closed-form exponentials, not a general CAS;
//! - stereographic unit-quaternion charts with unconstrained REAL parameters.
//!
//! Exact coefficients inherit checked-i128 and polynomial budgets from qft-symbolics.
use crate::{
    scalar::Conjugate, su2, Complex, Error, Matrix, RationalScalar, Scalar, So3Algebra, Su2Basis,
};
use qft_symbolics::{
    polynomial::Poly,
    scalar::{C, Q},
};
use std::{
    collections::BTreeMap,
    fmt,
    ops::{Add, Mul, Neg, Sub},
};
/// Q[s,c]/(s^2+c^2-1), normal form a(s)+c*b(s). ONE common formal angle t.
/// Different independent angles are not representable in this small ring.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Trig {
    a: Poly,
    b: Poly,
}
impl Trig {
    pub fn sin() -> Self {
        Self {
            a: Poly::var("s"),
            b: Poly::constant(0),
        }
    }
    pub fn cos() -> Self {
        Self {
            a: Poly::constant(0),
            b: Poly::constant(1),
        }
    }
    pub fn normal_form(&self) -> (&Poly, &Poly) {
        (&self.a, &self.b)
    }
    pub fn at(&self, s: Q, c: Q) -> Result<Q, Error> {
        if s * s + c * c != Q::ONE {
            return Err(Error::Domain("not a unit circle"));
        }
        let values = BTreeMap::from([("s".into(), C::from(s))]);
        let a = self.a.eval(&values).map_err(Error::Evaluation)?;
        let b = self.b.eval(&values).map_err(Error::Evaluation)?;
        Ok(a.re + c * b.re)
    }
}
impl Scalar for Trig {
    fn from_i64(v: i64) -> Self {
        Self {
            a: Poly::constant(v),
            b: Poly::constant(0),
        }
    }
}
impl RationalScalar for Trig {
    fn from_q(v: Q) -> Self {
        Self {
            a: Poly::constant(v),
            b: Poly::constant(0),
        }
    }
}
impl Conjugate for Trig {
    fn conjugate(&self) -> Self {
        self.clone()
    }
}
impl Add for Trig {
    type Output = Self;
    fn add(self, x: Self) -> Self {
        Self {
            a: self.a + x.a,
            b: self.b + x.b,
        }
    }
}
impl Neg for Trig {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            a: -self.a,
            b: -self.b,
        }
    }
}
impl Sub for Trig {
    type Output = Self;
    fn sub(self, x: Self) -> Self {
        self + (-x)
    }
}
impl Mul for Trig {
    type Output = Self;
    fn mul(self, x: Self) -> Self {
        Self {
            a: self.a.clone() * x.a.clone()
                + (Poly::constant(1) - Poly::var("s").pow(2)) * self.b.clone() * x.b.clone(),
            b: self.a * x.b + self.b * x.a,
        }
    }
}
impl fmt::Display for Trig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}) + c*({})", self.a, self.b)
    }
}
fn unit_axis(axis: [Q; 3]) -> Result<[Q; 3], Error> {
    if axis.iter().fold(Q::ZERO, |s, x| s + *x * *x) != Q::ONE {
        return Err(Error::Domain("symbolic axis must have exact norm one"));
    }
    Ok(axis)
}
/// exp(theta*hat(n)); s=sin(theta), c=cos(theta). Exact rational unit axis.
pub fn so3_exponential(axis: [Q; 3]) -> Result<Matrix<Trig, 3, 3>, Error> {
    let n = unit_axis(axis)?.map(Trig::from_q);
    let k = So3Algebra::new(n).matrix();
    Ok(&(&Matrix::identity() + &k.scaled(Trig::sin()))
        + &(&k * &k).scaled(Trig::one() - Trig::cos()))
}
pub struct Exponential2 {
    pub matrix: Matrix<Complex<Trig>, 2, 2>,
    /// Display/interpret s=sin(argument_scale*theta), c=cos(argument_scale*theta).
    pub argument_scale: Q,
}
/// Full basis: exp(theta*n.F), s=sin(theta). Half basis: exp(theta*n.T),
/// s=sin(theta/2). Matrix expressions must retain this argument metadata.
pub fn su2_exponential(axis: [Q; 3], basis: Su2Basis) -> Result<Exponential2, Error> {
    let n = unit_axis(axis)?.map(Trig::from_q);
    let s = Trig::sin();
    Ok(Exponential2 {
        matrix: su2::quaternion_matrix(&[
            Trig::cos(),
            n[0].clone() * s.clone(),
            n[1].clone() * s.clone(),
            n[2].clone() * s,
        ]),
        argument_scale: match basis {
            Su2Basis::PauliOverI => Q::ONE,
            Su2Basis::HalfPauliOverI => Q::new(1, 2),
        },
    })
}
/// General rational matrix. A nonzero polynomial denominator does NOT exclude
/// all evaluation poles. This type alone makes no group-membership claim.
#[derive(Clone, Debug)]
pub struct RationalMatrix<const N: usize> {
    numerator: Matrix<Poly, N, N>,
    denominator: Poly,
}
impl<const N: usize> RationalMatrix<N> {
    pub fn new(numerator: Matrix<Poly, N, N>, denominator: Poly) -> Result<Self, Error> {
        if denominator.is_zero() {
            return Err(Error::Domain("identically zero denominator"));
        }
        Ok(Self {
            numerator,
            denominator,
        })
    }
    pub fn numerator(&self) -> &Matrix<Poly, N, N> {
        &self.numerator
    }
    pub fn denominator(&self) -> &Poly {
        &self.denominator
    }
    pub fn identity() -> Self {
        Self {
            numerator: Matrix::identity(),
            denominator: Poly::constant(1),
        }
    }
    pub fn multiply(&self, b: &Self) -> Self {
        Self {
            numerator: &self.numerator * &b.numerator,
            denominator: self.denominator.clone() * b.denominator.clone(),
        }
    }
    pub fn transpose(&self) -> Self {
        Self {
            numerator: self.numerator.transpose(),
            denominator: self.denominator.clone(),
        }
    }
    /// All formal symbols are interpreted as real, explicitly.
    pub fn adjoint_real_symbols(&self) -> Self {
        Self {
            numerator: self.numerator.map(Poly::conjugate_coefficients).transpose(),
            denominator: self.denominator.conjugate_coefficients(),
        }
    }
    pub fn residual(&self, b: &Self) -> Matrix<Poly, N, N> {
        &self.numerator.scaled(b.denominator.clone())
            - &b.numerator.scaled(self.denominator.clone())
    }
    /// Identity on the common nonzero-denominator domain, not equality at poles.
    pub fn equals(&self, b: &Self) -> bool {
        self.residual(b).rows().iter().flatten().all(Poly::is_zero)
    }
    pub fn evaluate(&self, values: &BTreeMap<String, C>) -> Result<Matrix<C, N, N>, Error> {
        let denominator = self.denominator.eval(values).map_err(Error::Evaluation)?;
        if denominator.is_zero() {
            return Err(Error::Domain("evaluation pole"));
        }
        let mut result = Matrix::zero();
        for i in 0..N {
            for j in 0..N {
                result[(i, j)] = self.numerator[(i, j)]
                    .eval(values)
                    .map_err(Error::Evaluation)?
                    / denominator;
            }
        }
        Ok(result)
    }
}
#[derive(Clone, Copy, Debug)]
pub enum Chart {
    Plus,
    Minus,
}
/// Unit quaternion rational chart, all symbols REAL. q=(1-|t|^2,2t)/(1+|t|^2).
/// Plus misses -I; Minus misses +I. Together they cover all SU(2).
#[derive(Clone, Debug)]
pub struct SymbolicSu2 {
    q: [Poly; 4],
    d: Poly,
}
impl SymbolicSu2 {
    pub fn stereographic(parameters: [Poly; 3], chart: Chart) -> Result<Self, Error> {
        if parameters.iter().any(|p| !p.has_real_coefficients()) {
            return Err(Error::Domain(
                "chart parameters must be real polynomials in real symbols",
            ));
        }
        let squared = parameters
            .iter()
            .fold(Poly::constant(0), |s, v| s + v.pow(2));
        let sign = Poly::constant(match chart {
            Chart::Plus => 1,
            Chart::Minus => -1,
        });
        Ok(Self {
            q: [
                (Poly::constant(1) - squared.clone()) * sign.clone(),
                parameters[0].clone() * Poly::constant(2) * sign.clone(),
                parameters[1].clone() * Poly::constant(2) * sign.clone(),
                parameters[2].clone() * Poly::constant(2) * sign,
            ],
            d: Poly::constant(1) + squared,
        })
    }
    pub fn identity() -> Self {
        Self {
            q: [
                Poly::constant(1),
                Poly::constant(0),
                Poly::constant(0),
                Poly::constant(0),
            ],
            d: Poly::constant(1),
        }
    }
    pub fn norm_residual(&self) -> Poly {
        self.q.iter().fold(Poly::constant(0), |s, v| s + v.pow(2)) - self.d.pow(2)
    }
    pub fn matrix(&self) -> RationalMatrix<2> {
        let numerator = su2::quaternion_matrix(&self.q)
            .map(|z| z.re.clone() + Poly::constant(C::I) * z.im.clone());
        RationalMatrix {
            numerator,
            denominator: self.d.clone(),
        }
    }
    pub fn so3_matrix(&self) -> RationalMatrix<3> {
        RationalMatrix {
            numerator: su2::quaternion_rotation(&self.q),
            denominator: self.d.pow(2),
        }
    }
    pub fn inverse(&self) -> Self {
        Self {
            q: [
                self.q[0].clone(),
                -self.q[1].clone(),
                -self.q[2].clone(),
                -self.q[3].clone(),
            ],
            d: self.d.clone(),
        }
    }
    pub fn antipode(&self) -> Self {
        Self {
            q: self.q.clone().map(|p| -p),
            d: self.d.clone(),
        }
    }
    pub fn compose(&self, b: &Self) -> Self {
        Self {
            q: su2::hamilton(&self.q, &b.q),
            d: self.d.clone() * b.d.clone(),
        }
    }
    pub fn equals(&self, b: &Self) -> bool {
        (0..4)
            .all(|i| (self.q[i].clone() * b.d.clone() - b.q[i].clone() * self.d.clone()).is_zero())
    }
    pub fn evaluate(&self, values: &BTreeMap<String, C>) -> Result<Matrix<C, 2, 2>, Error> {
        if values.values().any(|c| !c.im.is_zero()) {
            return Err(Error::Domain(
                "chart variables must evaluate to real values",
            ));
        }
        self.matrix().evaluate(values)
    }
}
/// SO(3) represented by a unit-quaternion lift; equality identifies opposite lifts.
/// A single stereographic lift chart covers SO(3), since the omitted -I has
/// the same image as +I. No square-root/angle-branch selection is needed.
#[derive(Clone, Debug)]
pub struct SymbolicSo3 {
    lift: SymbolicSu2,
}
impl SymbolicSo3 {
    pub fn stereographic(parameters: [Poly; 3]) -> Result<Self, Error> {
        Ok(Self::from_lift(SymbolicSu2::stereographic(
            parameters,
            Chart::Plus,
        )?))
    }
    pub fn from_lift(lift: SymbolicSu2) -> Self {
        Self { lift }
    }
    pub fn matrix(&self) -> RationalMatrix<3> {
        self.lift.so3_matrix()
    }
    pub fn inverse(&self) -> Self {
        Self::from_lift(self.lift.inverse())
    }
    pub fn compose(&self, b: &Self) -> Self {
        Self::from_lift(self.lift.compose(&b.lift))
    }
    pub fn equals(&self, b: &Self) -> bool {
        self.matrix().equals(&b.matrix())
    }
}
