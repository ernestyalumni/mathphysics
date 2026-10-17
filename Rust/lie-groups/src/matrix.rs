use crate::scalar::{Conjugate, Scalar};
use std::ops::{Add, Index, IndexMut, Mul, Sub};
/// Fixed dimensions and coefficient type are compile-time parameters.
/// General matrices deliberately allow mutation and addition; group types don't.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Matrix<T, const R: usize, const C: usize> {
    data: [[T; C]; R],
}
impl<T, const R: usize, const C: usize> Matrix<T, R, C> {
    pub fn from_rows(data: [[T; C]; R]) -> Self {
        Self { data }
    }
    pub fn rows(&self) -> &[[T; C]; R] {
        &self.data
    }
    pub fn map<U>(&self, mut f: impl FnMut(&T) -> U) -> Matrix<U, R, C> {
        Matrix::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| f(&self.data[i][j]))
        }))
    }
}
impl<T: Scalar, const R: usize, const C: usize> Matrix<T, R, C> {
    pub fn zero() -> Self {
        Self::from_rows(std::array::from_fn(|_| std::array::from_fn(|_| T::zero())))
    }
    pub fn transpose(&self) -> Matrix<T, C, R> {
        Matrix::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| self[(j, i)].clone())
        }))
    }
    pub fn scaled(&self, s: T) -> Self {
        self.map(|v| v.clone() * s.clone())
    }
    pub fn mul_vector(&self, v: &[T; C]) -> [T; R] {
        std::array::from_fn(|i| {
            (0..C).fold(T::zero(), |s, j| s + self[(i, j)].clone() * v[j].clone())
        })
    }
}
impl<T: Conjugate, const R: usize, const C: usize> Matrix<T, R, C> {
    pub fn adjoint(&self) -> Matrix<T, C, R> {
        self.map(Conjugate::conjugate).transpose()
    }
}
impl<T: Scalar, const N: usize> Matrix<T, N, N> {
    pub fn identity() -> Self {
        Self::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| if i == j { T::one() } else { T::zero() })
        }))
    }
    pub fn trace(&self) -> T {
        (0..N).fold(T::zero(), |s, i| s + self[(i, i)].clone())
    }
    pub fn commutator(&self, b: &Self) -> Self {
        &(self * b) - &(b * self)
    }
}
impl<T: Scalar> Matrix<T, 2, 2> {
    pub fn determinant(&self) -> T {
        self[(0, 0)].clone() * self[(1, 1)].clone() - self[(0, 1)].clone() * self[(1, 0)].clone()
    }
}
impl<T: Scalar> Matrix<T, 3, 3> {
    pub fn determinant(&self) -> T {
        self[(0, 0)].clone()
            * (self[(1, 1)].clone() * self[(2, 2)].clone()
                - self[(1, 2)].clone() * self[(2, 1)].clone())
            - self[(0, 1)].clone()
                * (self[(1, 0)].clone() * self[(2, 2)].clone()
                    - self[(1, 2)].clone() * self[(2, 0)].clone())
            + self[(0, 2)].clone()
                * (self[(1, 0)].clone() * self[(2, 1)].clone()
                    - self[(1, 1)].clone() * self[(2, 0)].clone())
    }
}
impl<T, const R: usize, const C: usize> Index<(usize, usize)> for Matrix<T, R, C> {
    type Output = T;
    fn index(&self, (i, j): (usize, usize)) -> &T {
        &self.data[i][j]
    }
}
impl<T, const R: usize, const C: usize> IndexMut<(usize, usize)> for Matrix<T, R, C> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        &mut self.data[i][j]
    }
}
impl<T: Scalar, const R: usize, const C: usize> Add for &Matrix<T, R, C> {
    type Output = Matrix<T, R, C>;
    fn add(self, b: Self) -> Self::Output {
        Matrix::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| self[(i, j)].clone() + b[(i, j)].clone())
        }))
    }
}
impl<T: Scalar, const R: usize, const C: usize> Sub for &Matrix<T, R, C> {
    type Output = Matrix<T, R, C>;
    fn sub(self, b: Self) -> Self::Output {
        Matrix::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| self[(i, j)].clone() - b[(i, j)].clone())
        }))
    }
}
impl<T: Scalar, const R: usize, const C: usize, const K: usize> Mul<&Matrix<T, C, K>>
    for &Matrix<T, R, C>
{
    type Output = Matrix<T, R, K>;
    fn mul(self, b: &Matrix<T, C, K>) -> Self::Output {
        Matrix::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..C).fold(T::zero(), |s, k| {
                    s + self[(i, k)].clone() * b[(k, j)].clone()
                })
            })
        }))
    }
}
