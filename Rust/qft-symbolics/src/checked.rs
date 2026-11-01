//! Additive fallible capabilities. Budgets bound declared work, not OS resources.
use crate::{
    coefficient::Coefficient,
    computation::{Budget, Error, Resource, Result},
    polynomial::Polynomial,
    scalar::C,
    tensor,
};
use std::collections::BTreeMap;
pub type Matrix<A, const N: usize> = [[A; N]; N];
fn matrix_preflight<A: Coefficient, const N: usize>(
    inputs: &[&Matrix<A, N>],
    b: &mut Budget,
) -> Result<()> {
    let cells = (N as u64).saturating_mul(N as u64);
    b.charge(cells.saturating_mul(inputs.len() as u64).saturating_add(1))?;
    b.check(Resource::Terms, cells)?;
    for m in inputs {
        for row in *m {
            for c in row {
                b.check(Resource::CoefficientBits, c.bits())?;
            }
        }
    }
    Ok(())
}
pub fn identity<A: Coefficient, const N: usize>(b: &mut Budget) -> Result<Matrix<A, N>> {
    matrix_preflight::<A, N>(&[], b)?;
    let one = A::from_i64(1);
    let zero = A::from_i64(0);
    b.check(Resource::CoefficientBits, one.bits().max(zero.bits()))?;
    Ok(std::array::from_fn(|i| {
        std::array::from_fn(|j| if i == j { one.clone() } else { zero.clone() })
    }))
}
pub fn add<A: Coefficient, const N: usize>(
    a: &Matrix<A, N>,
    c: &Matrix<A, N>,
    b: &mut Budget,
) -> Result<Matrix<A, N>> {
    matrix_preflight(&[a, c], b)?;
    let mut out = identity::<A, N>(b)?;
    for i in 0..N {
        for j in 0..N {
            out[i][j] = a[i][j].checked_add(&c[i][j], b)?;
        }
    }
    Ok(out)
}
pub fn sub<A: Coefficient, const N: usize>(
    a: &Matrix<A, N>,
    c: &Matrix<A, N>,
    b: &mut Budget,
) -> Result<Matrix<A, N>> {
    matrix_preflight(&[a, c], b)?;
    let mut out = identity::<A, N>(b)?;
    for i in 0..N {
        for j in 0..N {
            out[i][j] = a[i][j].checked_sub(&c[i][j], b)?;
        }
    }
    Ok(out)
}
pub fn scale<A: Coefficient, const N: usize>(
    a: &Matrix<A, N>,
    c: &A,
    b: &mut Budget,
) -> Result<Matrix<A, N>> {
    matrix_preflight(&[a], b)?;
    b.check(Resource::CoefficientBits, c.bits())?;
    let mut out = identity::<A, N>(b)?;
    for i in 0..N {
        for j in 0..N {
            out[i][j] = a[i][j].checked_mul(c, b)?;
        }
    }
    Ok(out)
}
pub fn mul<A: Coefficient, const N: usize>(
    a: &Matrix<A, N>,
    c: &Matrix<A, N>,
    b: &mut Budget,
) -> Result<Matrix<A, N>> {
    matrix_preflight(&[a, c], b)?;
    b.charge((N as u64).saturating_pow(3))?;
    let mut out = identity::<A, N>(b)?;
    for i in 0..N {
        for j in 0..N {
            let mut s = A::from_i64(0);
            for (k, row) in c.iter().enumerate() {
                s = s.checked_add(&a[i][k].checked_mul(&row[j], b)?, b)?;
            }
            out[i][j] = s;
        }
    }
    Ok(out)
}
pub fn trace<A: Coefficient, const N: usize>(a: &Matrix<A, N>, b: &mut Budget) -> Result<A> {
    matrix_preflight(&[a], b)?;
    let mut s = A::from_i64(0);
    b.check(Resource::CoefficientBits, s.bits())?;
    for (i, row) in a.iter().enumerate() {
        s = s.checked_add(&row[i], b)?;
    }
    Ok(s)
}
pub fn transpose<A: Coefficient, const N: usize>(
    a: &Matrix<A, N>,
    b: &mut Budget,
) -> Result<Matrix<A, N>> {
    matrix_preflight(&[a], b)?;
    Ok(std::array::from_fn(|i| {
        std::array::from_fn(|j| a[j][i].clone())
    }))
}
fn lorentz_index(mu: usize, b: &mut Budget) -> Result<()> {
    b.charge(1)?;
    if mu >= 4 {
        return Err(Error::UnsupportedDomain {
            message: "Lorentz component outside 0..4",
        });
    }
    Ok(())
}
pub fn sigma(mu: usize, b: &mut Budget) -> Result<Matrix<C, 2>> {
    lorentz_index(mu, b)?;
    matrix_preflight::<C, 2>(&[], b)?;
    b.check(Resource::CoefficientBits, 1)?;
    let minus_one = C::ONE.checked_neg(b)?;
    let minus_i = C::I.checked_neg(b)?;
    Ok(match mu {
        0 => [[C::ONE, C::ZERO], [C::ZERO, C::ONE]],
        1 => [[C::ZERO, C::ONE], [C::ONE, C::ZERO]],
        2 => [[C::ZERO, minus_i], [C::I, C::ZERO]],
        _ => [[C::ONE, C::ZERO], [C::ZERO, minus_one]],
    })
}
pub fn barsigma(mu: usize, b: &mut Budget) -> Result<Matrix<C, 2>> {
    let s = sigma(mu, b)?;
    scale(&s, &C::from_i64(if mu == 0 { 1 } else { -1 }), b)
}
pub fn gamma(mu: usize, b: &mut Budget) -> Result<Matrix<C, 4>> {
    lorentz_index(mu, b)?;
    matrix_preflight::<C, 4>(&[], b)?;
    let s = sigma(mu, b)?;
    let t = barsigma(mu, b)?;
    Ok(std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i < 2 && j >= 2 {
                s[i][j - 2]
            } else if i >= 2 && j < 2 {
                t[i - 2][j]
            } else {
                C::ZERO
            }
        })
    }))
}
/// No cancellation: explicit denominators retain the common evaluation domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RationalFunction<A: Coefficient> {
    numerator: Polynomial<A>,
    denominator: Polynomial<A>,
}
impl<A: Coefficient> RationalFunction<A> {
    pub fn new(
        numerator: Polynomial<A>,
        denominator: Polynomial<A>,
        b: &mut Budget,
    ) -> Result<Self> {
        if denominator.is_zero() {
            return Err(Error::ZeroDenominator);
        }
        // Reuse polynomial validation without accessing its representation.
        numerator.try_add(&Polynomial::default(), b)?;
        denominator.try_add(&Polynomial::default(), b)?;
        Ok(Self {
            numerator,
            denominator,
        })
    }
    pub fn numerator(&self) -> &Polynomial<A> {
        &self.numerator
    }
    pub fn denominator(&self) -> &Polynomial<A> {
        &self.denominator
    }
    pub fn try_eval(&self, values: &BTreeMap<String, A>, b: &mut Budget) -> Result<A> {
        let d = self.denominator.try_eval(values, b)?;
        if d.is_zero() {
            return Err(Error::DivisionByZero);
        }
        self.numerator.try_eval(values, b)?.checked_div(&d, b)
    }
    pub fn try_residual(&self, r: &Self, b: &mut Budget) -> Result<Polynomial<A>> {
        self.numerator
            .try_mul(&r.denominator, b)?
            .try_sub(&r.numerator.try_mul(&self.denominator, b)?, b)
    }
    pub fn try_add(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        let n = self
            .numerator
            .try_mul(&r.denominator, b)?
            .try_add(&r.numerator.try_mul(&self.denominator, b)?, b)?;
        Self::new(n, self.denominator.try_mul(&r.denominator, b)?, b)
    }
    pub fn try_neg(&self, b: &mut Budget) -> Result<Self> {
        // Validate denominator before cloning it.
        let d = self.denominator.try_add(&Polynomial::default(), b)?;
        Self::new(self.numerator.try_neg(b)?, d, b)
    }
    pub fn try_sub(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        self.try_add(&r.try_neg(b)?, b)
    }
    pub fn try_mul(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        Self::new(
            self.numerator.try_mul(&r.numerator, b)?,
            self.denominator.try_mul(&r.denominator, b)?,
            b,
        )
    }
    pub fn try_div(&self, r: &Self, b: &mut Budget) -> Result<Self> {
        if r.numerator.is_zero() {
            return Err(Error::DivisionByZero);
        }
        // Extra r.denominator on BOTH sides retains poles of r.
        let n = self
            .numerator
            .try_mul(&r.denominator, b)?
            .try_mul(&r.denominator, b)?;
        let d = self
            .denominator
            .try_mul(&r.numerator, b)?
            .try_mul(&r.denominator, b)?;
        Self::new(n, d, b)
    }
}
pub type Spinor<A> = [Polynomial<A>; 2];
pub fn bracket<A: Coefficient>(
    a: &Spinor<A>,
    c: &Spinor<A>,
    b: &mut Budget,
) -> Result<Polynomial<A>> {
    a[0].try_mul(&c[1], b)?.try_sub(&a[1].try_mul(&c[0], b)?, b)
}
pub fn schouten<A: Coefficient>(
    a: &Spinor<A>,
    c: &Spinor<A>,
    d: &Spinor<A>,
    e: &Spinor<A>,
    b: &mut Budget,
) -> Result<Polynomial<A>> {
    let x = bracket(a, c, b)?.try_mul(&bracket(d, e, b)?, b)?;
    let y = bracket(a, d, b)?.try_mul(&bracket(e, c, b)?, b)?;
    let z = bracket(a, e, b)?.try_mul(&bracket(c, d, b)?, b)?;
    x.try_add(&y, b)?.try_add(&z, b)
}
/// [i,j> shift, staged so every error leaves both inputs unchanged.
pub fn bcfw_shift<A: Coefficient>(
    lambda: &mut [[A; 2]],
    tilde: &mut [[A; 2]],
    i: usize,
    j: usize,
    z: &A,
    b: &mut Budget,
) -> Result<()> {
    b.charge(1)?;
    if lambda.len() != tilde.len() || i == j || i >= lambda.len() || j >= lambda.len() {
        return Err(Error::UnsupportedDomain {
            message: "invalid BCFW legs",
        });
    }
    b.check(Resource::CoefficientBits, z.bits())?;
    for v in lambda[i]
        .iter()
        .chain(&lambda[j])
        .chain(&tilde[i])
        .chain(&tilde[j])
    {
        b.check(Resource::CoefficientBits, v.bits())?;
    }
    b.charge(4)?;
    let l0 = lambda[i][0].checked_add(&z.checked_mul(&lambda[j][0], b)?, b)?;
    let l1 = lambda[i][1].checked_add(&z.checked_mul(&lambda[j][1], b)?, b)?;
    let t0 = tilde[j][0].checked_sub(&z.checked_mul(&tilde[i][0], b)?, b)?;
    let t1 = tilde[j][1].checked_sub(&z.checked_mul(&tilde[i][1], b)?, b)?;
    lambda[i] = [l0, l1];
    tilde[j] = [t0, t1];
    Ok(())
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TensorError {
    Computation(Error),
    Validation(tensor::Error),
}
impl std::fmt::Display for TensorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Computation(e) => e.fmt(f),
            Self::Validation(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for TensorError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Computation(e) => e,
            Self::Validation(e) => e,
        })
    }
}
/// Declaration/metric validator only; no legacy tensor coefficient arithmetic.
pub struct TensorValidator {
    context: tensor::Context,
    orbit_steps: usize,
    declarations: u64,
    group_entries: u64,
    storage_bytes: u64,
}
impl TensorValidator {
    pub fn new(orbit_steps: usize) -> Self {
        Self {
            context: tensor::Context::default().with_budget(orbit_steps),
            orbit_steps,
            declarations: 0,
            group_entries: 0,
            storage_bytes: 0,
        }
    }
    pub fn orbit_steps(&self) -> usize {
        self.orbit_steps
    }
    /// Conservative retained bounds; only successful declarations contribute.
    pub fn retained_bounds(&self) -> (u64, u64, u64) {
        (self.declarations, self.group_entries, self.storage_bytes)
    }
    pub fn context(&self) -> &tensor::Context {
        &self.context
    }
    pub fn require_symmetric_metric(
        &self,
        space: tensor::Space,
        b: &mut Budget,
    ) -> std::result::Result<(), TensorError> {
        b.charge(1).map_err(TensorError::Computation)?;
        self.context
            .require_symmetric_metric(space)
            .map_err(TensorError::Validation)
    }
    pub fn declare(
        &mut self,
        name: &str,
        slots: Vec<tensor::Slot>,
        symmetry: tensor::Symmetry,
        b: &mut Budget,
    ) -> std::result::Result<(), TensorError> {
        let rank = slots.len();
        if rank > 8 {
            return Err(TensorError::Validation(tensor::Error::Unsupported(
                "tensor rank >8".into(),
            )));
        }
        let preflight = |b: &mut Budget| -> Result<(u64, u64, u64)> {
            b.check(Resource::InputBytes, name.len() as u64)?;
            b.check(Resource::WordLength, rank as u64)?;
            b.charge((rank as u64).saturating_add(1))?;
            let factorial = (1..=rank as u64).product::<u64>();
            let (group, generators) = match &symmetry {
                tensor::Symmetry::None => (1, 0),
                tensor::Symmetry::Symmetric | tensor::Symmetry::Antisymmetric => {
                    (factorial, rank.saturating_sub(1) as u64)
                }
                tensor::Symmetry::Generators(g) => {
                    b.check(Resource::Terms, g.len() as u64)?;
                    for (p, _) in g {
                        b.check(Resource::WordLength, p.len() as u64)?;
                        b.charge((p.len() as u64).saturating_add(1))?;
                    }
                    (factorial.saturating_mul(2), g.len() as u64)
                }
            };
            let declarations = self.declarations.saturating_add(1);
            let groups = self.group_entries.saturating_add(group);
            // Retained logical storage bound includes name, slots, map node,
            // permutation vectors, signs and generous per-entry bookkeeping.
            // This is not a promise about allocator implementation overhead.
            let entry_bytes = (rank as u64)
                .saturating_mul(std::mem::size_of::<usize>() as u64)
                .saturating_add(256);
            let added_bytes = (name.len() as u64)
                .saturating_add(
                    (rank as u64).saturating_mul(std::mem::size_of::<tensor::Slot>() as u64),
                )
                .saturating_add(512)
                .saturating_add(group.saturating_mul(entry_bytes));
            let storage = self.storage_bytes.saturating_add(added_bytes);
            b.check(Resource::Terms, declarations.saturating_add(groups))?;
            b.check(Resource::InputBytes, storage)?;
            // Reserve simultaneous queue/set/output and generator clones. The
            // linear group comparison bound avoids relying on BTree internals.
            let transient_entries = group.saturating_mul(4).saturating_add(generators);
            b.check(
                Resource::Terms,
                transient_entries
                    .saturating_add(self.group_entries)
                    .saturating_add(declarations),
            )?;
            b.check(
                Resource::InputBytes,
                self.storage_bytes
                    .saturating_add(added_bytes.saturating_mul(4))
                    .saturating_add(generators.saturating_mul(entry_bytes)),
            )?;
            let work = group
                .saturating_mul(group.saturating_add(1))
                .saturating_mul(generators.saturating_add(1))
                .saturating_mul((rank as u64).saturating_add(1))
                .saturating_mul(8);
            b.charge(work)?;
            Ok((declarations, groups, storage))
        };
        let (declarations, groups, storage) = preflight(b).map_err(TensorError::Computation)?;
        self.context
            .declare(name, slots, symmetry)
            .map_err(TensorError::Validation)?;
        self.declarations = declarations;
        self.group_entries = groups;
        self.storage_bytes = storage;
        Ok(())
    }
}
