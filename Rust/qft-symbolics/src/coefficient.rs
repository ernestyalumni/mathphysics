//! A fallible coefficient contract shared by bounded and arbitrary-precision polynomials.
use crate::computation::{Budget, Result};
use std::fmt::{Debug, Display};
pub trait Coefficient: Clone + Default + Eq + Debug + Display {
    fn from_i64(value: i64) -> Self;
    fn is_zero(&self) -> bool;
    fn is_real(&self) -> bool;
    fn conjugate(&self) -> Self;
    /// Maximum numerator/denominator component bit length, not total heap usage.
    fn bits(&self) -> u64;
    fn checked_add(&self, rhs: &Self, budget: &mut Budget) -> Result<Self>;
    fn checked_sub(&self, rhs: &Self, budget: &mut Budget) -> Result<Self>;
    fn checked_mul(&self, rhs: &Self, budget: &mut Budget) -> Result<Self>;
    fn checked_div(&self, rhs: &Self, budget: &mut Budget) -> Result<Self>;
    fn checked_neg(&self, budget: &mut Budget) -> Result<Self>;
}
