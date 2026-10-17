use crate::{
    default_tolerance, tolerance, vector, Complex, Error, Matrix, RationalScalar, Real, Scalar,
    So3, So3Algebra,
};
use qft_symbolics::scalar::Q;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Su2Basis {
    HalfPauliOverI,
    PauliOverI,
}
pub fn pauli<T: Scalar>(index: usize) -> Result<Matrix<Complex<T>, 2, 2>, Error> {
    let z = Complex::zero();
    let o = Complex::one();
    let i = Complex::new(T::zero(), T::one());
    Ok(match index {
        0 => Matrix::from_rows([[z.clone(), o.clone()], [o, z]]),
        1 => Matrix::from_rows([[z.clone(), -i.clone()], [i, z]]),
        2 => Matrix::from_rows([[o.clone(), z.clone()], [z, -o]]),
        _ => {
            return Err(Error::IndexOutOfRange {
                index,
                dimension: 3,
            })
        }
    })
}
pub fn generator<T: RationalScalar>(
    index: usize,
    basis: Su2Basis,
) -> Result<Matrix<Complex<T>, 2, 2>, Error> {
    let denominator = match basis {
        Su2Basis::HalfPauliOverI => 2,
        Su2Basis::PauliOverI => 1,
    };
    Ok(pauli(index)?.scaled(Complex::new(T::zero(), T::from_q(Q::new(-1, denominator)))))
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Su2Algebra<T> {
    coordinates: [T; 3],
}
impl<T: RationalScalar> Su2Algebra<T> {
    /// Coordinates in T_i=-i sigma_i/2, so the bracket is the ordinary cross product.
    pub fn new(coordinates: [T; 3]) -> Self {
        Self { coordinates }
    }
    pub fn basis(index: usize) -> Result<Self, Error> {
        if index >= 3 {
            return Err(Error::IndexOutOfRange {
                index,
                dimension: 3,
            });
        }
        Ok(Self::new(std::array::from_fn(|j| {
            if index == j {
                T::one()
            } else {
                T::zero()
            }
        })))
    }
    pub fn try_from_matrix(m: Matrix<Complex<T>, 2, 2>) -> Result<Self, Error>
    where
        T: crate::scalar::Conjugate,
    {
        if m.adjoint() != m.scaled(Complex::real(-T::one())) || m.trace() != Complex::zero() {
            return Err(Error::Domain(
                "matrix must be exactly traceless and anti-Hermitian",
            ));
        }
        let result = Self::new([
            -(m[(0, 1)].im.clone() + m[(1, 0)].im.clone()),
            m[(1, 0)].re.clone() - m[(0, 1)].re.clone(),
            m[(1, 1)].im.clone() - m[(0, 0)].im.clone(),
        ]);
        if result.matrix() != m {
            return Err(Error::Domain("matrix outside real su2 coefficients"));
        }
        Ok(result)
    }
    pub fn from_full_basis_coordinates(v: [T; 3]) -> Self {
        Self::new(v.map(|x| x * T::from_i64(2)))
    }
    pub fn coordinates(&self) -> &[T; 3] {
        &self.coordinates
    }
    pub fn matrix(&self) -> Matrix<Complex<T>, 2, 2> {
        let mut m = Matrix::zero();
        for i in 0..3 {
            m = &m
                + &generator(i, Su2Basis::HalfPauliOverI)
                    .expect("fixed basis index")
                    .scaled(Complex::real(self.coordinates[i].clone()));
        }
        m
    }
    pub fn bracket(&self, b: &Self) -> Self {
        Self::new(vector::cross(&self.coordinates, &b.coordinates))
    }
    pub fn to_so3(&self) -> So3Algebra<T> {
        So3Algebra::new(self.coordinates.clone())
    }
}
pub fn hamilton<T: Scalar>(a: &[T; 4], b: &[T; 4]) -> [T; 4] {
    let av = [a[1].clone(), a[2].clone(), a[3].clone()];
    let bv = [b[1].clone(), b[2].clone(), b[3].clone()];
    let cross = vector::cross(&av, &bv);
    [
        a[0].clone() * b[0].clone() - vector::dot(&av, &bv),
        a[0].clone() * b[1].clone() + b[0].clone() * a[1].clone() + cross[0].clone(),
        a[0].clone() * b[2].clone() + b[0].clone() * a[2].clone() + cross[1].clone(),
        a[0].clone() * b[3].clone() + b[0].clone() * a[3].clone() + cross[2].clone(),
    ]
}
pub fn quaternion_matrix<T: Scalar>(q: &[T; 4]) -> Matrix<Complex<T>, 2, 2> {
    let [a, b, c, d] = q.clone();
    Matrix::from_rows([
        [
            Complex::new(a.clone(), -d.clone()),
            Complex::new(-c.clone(), -b.clone()),
        ],
        [Complex::new(c, -b), Complex::new(a, d)],
    ])
}
pub fn quaternion_rotation<T: Scalar>(q: &[T; 4]) -> Matrix<T, 3, 3> {
    let vector = [q[1].clone(), q[2].clone(), q[3].clone()];
    let k = So3Algebra::new(vector.clone()).matrix();
    let scalar = q[0].clone() * q[0].clone() - vector::dot(&vector, &vector);
    let mut result = Matrix::identity().scaled(scalar);
    result = &result + &k.scaled(T::from_i64(2) * q[0].clone());
    for i in 0..3 {
        for j in 0..3 {
            result[(i, j)] =
                result[(i, j)].clone() + T::from_i64(2) * vector[i].clone() * vector[j].clone();
        }
    }
    result
}
/// A numerical unit quaternion, not modulo sign: SU(2), not SO(3).
/// ```compile_fail
/// use lie_groups::Su2;
/// let u=Su2{q:[0.0_f64;4]};
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Su2<T> {
    q: [T; 4],
}
impl<T: Real> Su2<T> {
    pub fn identity() -> Self {
        Self {
            q: [T::one(), T::zero(), T::zero(), T::zero()],
        }
    }
    pub fn quaternion(&self) -> &[T; 4] {
        &self.q
    }
    pub fn try_from_unit_quaternion(q: [T; 4], tol: T) -> Result<Self, Error> {
        tolerance(tol)?;
        if q.iter().any(|x| !x.is_finite()) {
            return Err(Error::NonFinite);
        }
        let squared = q.iter().fold(T::zero(), |s, x| s + *x * *x);
        if !squared.is_finite() || (squared - T::one()).abs() > tol {
            return Err(Error::NotSu2);
        }
        Ok(Self { q })
    }
    pub fn matrix(&self) -> Matrix<Complex<T>, 2, 2> {
        quaternion_matrix(&self.q)
    }
    pub fn try_from_matrix(m: Matrix<Complex<T>, 2, 2>, tol: T) -> Result<Self, Error> {
        tolerance(tol)?;
        if m.rows()
            .iter()
            .flatten()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(Error::NonFinite);
        }
        let residual = &(&m.adjoint() * &m) - &Matrix::identity();
        let det = m.determinant() - Complex::one();
        if residual
            .rows()
            .iter()
            .flatten()
            .chain(std::iter::once(&det))
            .any(|z| !z.re.is_finite() || !z.im.is_finite() || z.re.abs() > tol || z.im.abs() > tol)
        {
            return Err(Error::NotSu2);
        }
        let half = T::from_f64(0.5);
        let q = [
            (m[(0, 0)].re + m[(1, 1)].re) * half,
            -(m[(0, 1)].im + m[(1, 0)].im) * half,
            (m[(1, 0)].re - m[(0, 1)].re) * half,
            (m[(1, 1)].im - m[(0, 0)].im) * half,
        ];
        let result = Self::try_from_unit_quaternion(q, tol)?;
        let difference = &m - &result.matrix();
        if difference
            .rows()
            .iter()
            .flatten()
            .any(|z| z.re.abs() > tol || z.im.abs() > tol)
        {
            return Err(Error::NotSu2);
        }
        Ok(result)
    }
    pub fn inverse(&self) -> Self {
        Self {
            q: [self.q[0], -self.q[1], -self.q[2], -self.q[3]],
        }
    }
    pub fn antipode(&self) -> Self {
        Self {
            q: self.q.map(|x| -x),
        }
    }
    pub fn compose(&self, b: &Self) -> Result<Self, Error> {
        Self::try_from_unit_quaternion(hamilton(&self.q, &b.q), default_tolerance())
    }
    pub fn to_so3(&self) -> Result<So3<T>, Error> {
        So3::try_from_matrix(quaternion_rotation(&self.q), default_tolerance())
    }
    pub fn exp_pauli_over_i(v: [T; 3]) -> Result<Self, Error> {
        let r = vector::norm(v)?;
        let mut q = [r.cos(), T::zero(), T::zero(), T::zero()];
        if r != T::zero() {
            for i in 0..3 {
                q[i + 1] = (v[i] / r) * r.sin();
            }
        }
        Self::try_from_unit_quaternion(q, default_tolerance())
    }
}
impl<T: Real> Su2Algebra<T> {
    pub fn exp(&self) -> Result<Su2<T>, Error> {
        let v = self.coordinates;
        let r = vector::norm(v)?;
        let half = r * T::from_f64(0.5);
        let mut q = [half.cos(), T::zero(), T::zero(), T::zero()];
        if r.abs() < T::from_f64(0.125) {
            let h2 = half * half;
            let sinc = T::one()
                + h2 * (T::from_f64(-1. / 6.)
                    + h2 * (T::from_f64(1. / 120.) + h2 * T::from_f64(-1. / 5040.)));
            for i in 0..3 {
                q[i + 1] = v[i] * T::from_f64(0.5) * sinc;
            }
        } else {
            for i in 0..3 {
                q[i + 1] = (v[i] / r) * half.sin();
            }
        }
        Su2::try_from_unit_quaternion(q, default_tolerance())
    }
}
