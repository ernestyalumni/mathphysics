//! Numerical and exact symbolic SO(3)/SU(2), with explicit Frankel conventions.
//! Column vectors, active right-handed rotations, epsilon(0,1,2)=+1.
//! Group types contain matrices/coordinates; they do not inherit matrix mutators.
//!
//! ```
//! use lie_groups::{Matrix,So3Algebra,Half,Real};
//! fn quarter_turn<T:Real>() -> lie_groups::So3<T> {
//!     So3Algebra::new([T::from_f64(0.0),T::from_f64(0.0),
//!                     T::from_f64(std::f64::consts::FRAC_PI_2)]).exp().unwrap()
//! }
//! let _: Matrix<f32,3,3> = Matrix::identity();
//! let _: Matrix<f64,3,3> = Matrix::identity();
//! let _: Matrix<Half,3,3> = Matrix::identity();
//! let _ = (quarter_turn::<f32>(),quarter_turn::<f64>(),quarter_turn::<Half>());
//! ```
//!
//! Incompatible matrix dimensions are rejected at compile time:
//! ```compile_fail
//! use lie_groups::Matrix;
//! let a=Matrix::<f64,2,3>::zero();
//! let b=Matrix::<f64,2,2>::identity();
//! let impossible=&a*&b;
//! ```
pub mod goldstein;
pub mod half;
pub mod matrix;
pub mod multitrig;
pub mod rigid_body;
pub mod scalar;
pub mod so3;
pub mod su2;
pub mod symbolic;
pub mod vector;
pub use half::Half;
pub use matrix::Matrix;
pub use scalar::{Complex, RationalScalar, Real, Scalar};
pub use so3::{So3, So3Algebra};
use std::fmt;
pub use su2::{Su2, Su2Algebra, Su2Basis};
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    IndexOutOfRange { index: usize, dimension: usize },
    Domain(&'static str),
    NonFinite,
    NotSo3,
    NotSu2,
    Evaluation(String),
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for Error {}
/// A Copy error keeps alternating-symbol evaluation usable in const contexts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndexError {
    pub index: usize,
    pub dimension: usize,
}
impl fmt::Display for IndexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "index {} outside dimension {}",
            self.index, self.dimension
        )
    }
}
impl std::error::Error for IndexError {}
impl From<IndexError> for Error {
    fn from(e: IndexError) -> Self {
        Self::IndexOutOfRange {
            index: e.index,
            dimension: e.dimension,
        }
    }
}
/// Alternating symbol on zero-based indices. NOT a metric volume form.
/// Evaluates at compile time too. Repetition gives zero; out-of-range is invalid.
pub const fn levi_civita<const N: usize>(indices: [usize; N]) -> Result<i8, IndexError> {
    let mut i = 0;
    while i < N {
        if indices[i] >= N {
            return Err(IndexError {
                index: indices[i],
                dimension: N,
            });
        }
        i += 1;
    }
    let mut sign = 1;
    i = 0;
    while i < N {
        let mut j = i + 1;
        while j < N {
            if indices[i] == indices[j] {
                return Ok(0);
            }
            if indices[i] > indices[j] {
                sign = -sign;
            }
            j += 1;
        }
        i += 1;
    }
    Ok(sign)
}
pub fn default_tolerance<T: Real>() -> T {
    T::from_f64(32.0) * T::epsilon()
}
pub(crate) fn tolerance<T: Real>(value: T) -> Result<(), Error> {
    if !value.is_finite() || value < T::zero() || value >= T::from_f64(0.125) {
        return Err(Error::Domain(
            "tolerance must be finite, nonnegative and <1/8",
        ));
    }
    Ok(())
}
