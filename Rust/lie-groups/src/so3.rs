use crate::{default_tolerance, tolerance, vector, Error, Matrix, Real, Scalar};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct So3Algebra<T> {
    coordinates: [T; 3],
}
impl<T: Scalar> So3Algebra<T> {
    pub fn new(coordinates: [T; 3]) -> Self {
        Self { coordinates }
    }
    pub fn coordinates(&self) -> &[T; 3] {
        &self.coordinates
    }
    pub fn basis(index: usize) -> Result<Self, Error> {
        if index >= 3 {
            return Err(Error::IndexOutOfRange {
                index,
                dimension: 3,
            });
        }
        Ok(Self::new(std::array::from_fn(|j| {
            if j == index {
                T::one()
            } else {
                T::zero()
            }
        })))
    }
    pub fn matrix(&self) -> Matrix<T, 3, 3> {
        let [x, y, z] = self.coordinates.clone();
        Matrix::from_rows([
            [T::zero(), -z.clone(), y.clone()],
            [z, T::zero(), -x.clone()],
            [-y, x, T::zero()],
        ])
    }
    pub fn bracket(&self, b: &Self) -> Self {
        Self::new(vector::cross(&self.coordinates, &b.coordinates))
    }
    pub fn try_from_matrix(m: Matrix<T, 3, 3>) -> Result<Self, Error> {
        if m.transpose() != m.scaled(-T::one()) || (0..3).any(|i| m[(i, i)] != T::zero()) {
            return Err(Error::Domain("matrix must be exactly skew-symmetric"));
        }
        Ok(Self::new([
            m[(2, 1)].clone(),
            m[(0, 2)].clone(),
            m[(1, 0)].clone(),
        ]))
    }
}
/// Approximately SO(3) at a declared floating-point tolerance.
/// Private storage prevents invalid mutation, including determinant-minus-one inputs.
/// ```compile_fail
/// use lie_groups::{So3,Matrix};
/// let r=So3{matrix:Matrix::<f64,3,3>::zero()};
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct So3<T> {
    matrix: Matrix<T, 3, 3>,
}
impl<T: Real> So3<T> {
    pub fn identity() -> Self {
        Self {
            matrix: Matrix::identity(),
        }
    }
    pub fn matrix(&self) -> &Matrix<T, 3, 3> {
        &self.matrix
    }
    pub fn try_from_matrix(matrix: Matrix<T, 3, 3>, tol: T) -> Result<Self, Error> {
        tolerance(tol)?;
        if matrix.rows().iter().flatten().any(|x| !x.is_finite()) {
            return Err(Error::NonFinite);
        }
        let residual = &(&matrix.transpose() * &matrix) - &Matrix::identity();
        let determinant = matrix.determinant();
        if !determinant.is_finite()
            || residual
                .rows()
                .iter()
                .flatten()
                .any(|x| !x.is_finite() || x.abs() > tol)
            || (determinant - T::one()).abs() > tol
        {
            return Err(Error::NotSo3);
        }
        Ok(Self { matrix })
    }
    pub fn from_axis_angle(axis: [T; 3], angle: T) -> Result<Self, Error> {
        if !angle.is_finite() {
            return Err(Error::NonFinite);
        }
        let n = vector::normalized(axis)?;
        So3Algebra::new(n.map(|x| x * angle)).exp()
    }
    pub fn compose(&self, b: &Self) -> Result<Self, Error> {
        Self::try_from_matrix(&self.matrix * &b.matrix, default_tolerance())
    }
    /// Revalidate the transpose: an entrywise M^T M error bound need not also
    /// bound M M^T. Floating group operations can fail rather than leak a bad value.
    pub fn inverse(&self) -> Result<Self, Error> {
        Self::try_from_matrix(self.matrix.transpose(), default_tolerance())
    }
    pub fn apply(&self, v: &[T; 3]) -> Result<[T; 3], Error> {
        if v.iter().any(|x| !x.is_finite()) {
            return Err(Error::NonFinite);
        }
        let result = self.matrix.mul_vector(v);
        if result.iter().any(|x| !x.is_finite()) {
            Err(Error::NonFinite)
        } else {
            Ok(result)
        }
    }
}
impl<T: Real> So3Algebra<T> {
    pub fn exp(&self) -> Result<So3<T>, Error> {
        let v = self.coordinates;
        let r = vector::norm(v)?;
        let matrix = if r.abs() < T::from_f64(0.125) {
            let x = r * r;
            let a = T::one()
                + x * (T::from_f64(-1. / 6.)
                    + x * (T::from_f64(1. / 120.)
                        + x * (T::from_f64(-1. / 5040.) + x * T::from_f64(1. / 362880.))));
            let b = T::from_f64(0.5)
                + x * (T::from_f64(-1. / 24.)
                    + x * (T::from_f64(1. / 720.)
                        + x * (T::from_f64(-1. / 40320.) + x * T::from_f64(1. / 3628800.))));
            let k = self.matrix();
            &(&Matrix::identity() + &k.scaled(a)) + &(&k * &k).scaled(b)
        } else {
            let n = v.map(|x| x / r);
            let k = Self::new(n).matrix();
            &(&Matrix::identity() + &k.scaled(r.sin())) + &(&k * &k).scaled(T::one() - r.cos())
        };
        So3::try_from_matrix(matrix, default_tolerance())
    }
}
