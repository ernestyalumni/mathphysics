//! Explicit epsilon conventions for commuting two-component spinors.
//!
//! Indices in arrays are (1,2). Raising is epsilon^{ab} v_b; lowering
//! is epsilon_{ab} v^b. The matrices are inverses, not equal. Brackets
//! here are commuting determinants, NOT anticommuting field bilinears (35.25).
//! Conversion holds lower components fixed; upper components and brackets
//! reverse sign between profiles. No existing library default is changed.
use crate::{
    coefficient::Coefficient,
    computation::{Budget, Result},
    polynomial::Polynomial,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EpsilonConvention {
    /// Srednicki (35.1): epsilon^{12}=+1, epsilon_{12}=-1.
    Library,
    /// Mono-notes profile: epsilon^{12}=-1, epsilon_{12}=+1.
    MonoNotes,
}
impl EpsilonConvention {
    pub const fn upper12(self) -> i64 {
        match self {
            Self::Library => 1,
            Self::MonoNotes => -1,
        }
    }
    pub const fn lower12(self) -> i64 {
        -self.upper12()
    }
    fn scale<A: Coefficient>(
        v: &Polynomial<A>,
        sign: i64,
        b: &mut Budget,
    ) -> Result<Polynomial<A>> {
        v.try_mul(&Polynomial::try_constant(A::from_i64(sign), b)?, b)
    }
    /// epsilon^{ab} v_b = upper12 * (v_2,-v_1).
    pub fn raise<A: Coefficient>(
        self,
        v: &[Polynomial<A>; 2],
        b: &mut Budget,
    ) -> Result<[Polynomial<A>; 2]> {
        Ok([
            Self::scale(&v[1], self.upper12(), b)?,
            Self::scale(&v[0], -self.upper12(), b)?,
        ])
    }
    /// epsilon_{ab} v^b = lower12 * (v^2,-v^1).
    pub fn lower<A: Coefficient>(
        self,
        v: &[Polynomial<A>; 2],
        b: &mut Budget,
    ) -> Result<[Polynomial<A>; 2]> {
        Ok([
            Self::scale(&v[1], self.lower12(), b)?,
            Self::scale(&v[0], -self.lower12(), b)?,
        ])
    }
    /// epsilon^{ab} a_a c_b, for commuting components.
    pub fn bracket<A: Coefficient>(
        self,
        a: &[Polynomial<A>; 2],
        c: &[Polynomial<A>; 2],
        b: &mut Budget,
    ) -> Result<Polynomial<A>> {
        let x = a[0].try_mul(&c[1], b)?;
        let y = a[1].try_mul(&c[0], b)?;
        Self::scale(&x.try_sub(&y, b)?, self.upper12(), b)
    }
    /// Keep the lower-index components fixed, validating through the budget.
    pub fn convert_lower<A: Coefficient>(
        self,
        _target: Self,
        v: &[Polynomial<A>; 2],
        b: &mut Budget,
    ) -> Result<[Polynomial<A>; 2]> {
        Ok([Self::scale(&v[0], 1, b)?, Self::scale(&v[1], 1, b)?])
    }
    /// Convert upper components while preserving the lower component array.
    pub fn convert_upper<A: Coefficient>(
        self,
        target: Self,
        v: &[Polynomial<A>; 2],
        b: &mut Budget,
    ) -> Result<[Polynomial<A>; 2]> {
        let sign = self.upper12() * target.upper12();
        Ok([Self::scale(&v[0], sign, b)?, Self::scale(&v[1], sign, b)?])
    }
}
