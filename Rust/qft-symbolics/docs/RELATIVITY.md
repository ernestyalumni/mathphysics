# Rapidity, Lorentz vectors and boosts

`src/relativity.rs` provides reusable time-first four-vector operations. The
metric is `diag(-1,1,1,1)`; vectors are column components `(t,x,y,z)`.
`boost_f64(eta,n)` uses the Srednicki Chapter 2, (2.37) positive off-diagonal
convention:

\[
L^0{}_0=\cosh\eta,\quad L^0{}_i=L^i{}_0=n_i\sinh\eta,\qquad
L^i{}_j=\delta_{ij}+n_i n_j(\cosh\eta-1).
\]

A positive boost maps `(1,0,0,0)` toward positive spatial components. A frame
transformation using the opposite velocity requires `eta.inverse()` explicitly.
No implicit frame/sign conversion occurs. Composition is left matrix action;
collinear rapidities add. General noncollinear composition remains an ordinary
matrix product and is not represented by `collinear_add`.

## Public API

All fallible functions return `computation::Result`, reusing its existing error
and budget types. Numerical domain/overflow errors use `UnsupportedDomain`.

- `Rapidity::new(f64)`, `from_beta(f64)`; `eta() -> f64`,
  `beta() -> Result<f64>`, `inverse() -> Self`,
  `collinear_add(Self) -> Result<Self>`.
- `Axis::{X,Y,Z}`; `UnitAxis::axial(Axis) -> Self`,
  `UnitAxis::new([f64;3], tolerance) -> Result<Self>`, `components()`.
- `boost_f64(Rapidity, UnitAxis) -> Result<Matrix<f64,4>>`;
  `apply_f64(&Matrix<f64,4>, &[f64;4]) -> Result<[f64;4]>`;
  `dot_f64(&[f64;4], &[f64;4]) -> Result<f64>`.
- `metric_residual_f64(&Matrix<f64,4>, Tolerance) -> Result<Residual>`.
- Generic `metric<A: Coefficient>(&mut Budget) -> Result<Matrix<A,4>>`;
  `lower/raise(&[A;4], &mut Budget) -> Result<[A;4]>`;
  `dot(&[A;4], &[A;4], &mut Budget) -> Result<A>`;
  `apply(&Matrix<A,4>, &[A;4], &mut Budget) -> Result<[A;4]>`.
- `exact_half_velocity_boost(&BigQ, Axis, &mut Budget)` returns
  `Result<Matrix<BigC,4>>`.
- `symbolic_axial_boost(c_label, s_label, Axis, &mut Budget)` returns
  `Result<SymbolicBoost>` with public `matrix: Matrix<BigPoly,4>` and
  `condition: BigPoly`.

`lower` and `raise` both flip the time component, since this metric is its own
inverse. `dot` contracts two contravariant vectors with the metric; it is
bilinear, with no complex conjugation. Applying a matrix is an ordinary
matrix-vector product. Generic operations support coefficient implementations
such as small `C` and arbitrary precision `BigC`; `BigQ` supplies the real
parameter to the exact boost, rather than implementing `Coefficient` itself.

## Numerical policy

Rapidity accepts every finite eta; boost construction can still fail when its
hyperbolic functions overflow. `from_beta` requires finite `|beta|<1`.
`beta()` computes `tanh(eta)` and **returns an error** if rounding produces
`+1` or `-1`. It never clamps to a nearby float or silently reports light speed.
Thus a finite valid rapidity need not have a representable subluminal beta.
Addition rejects overflow; inverse negates the finite eta.

A supplied axis is validated against an absolute tolerance on its Euclidean
norm. The tolerance must be finite and within `[0,1e-6]`. The axis is retained
unchanged: acceptance within tolerance is approximate, not exact normalization.
Use `UnitAxis::axial` for a coordinate axis. The spatial correction uses
`2*sinh(eta/2)^2`, retaining tiny off-diagonal spatial corrections when
`cosh(eta)-1` would round to zero. Underflow and ordinary floating point rounding
can still occur; neither constitutes an exact symbolic calculation.

For `R=L^T g L-g`, diagnostics report `absolute=max_ij |Rij|`,
`scale=max(1,max_ij sum_k |Lki Lkj|)` and `scaled=absolute/scale`.
`Tolerance { absolute, relative }` requires finite nonnegative components.
`within_tolerance` means `absolute <= tolerance.absolute + tolerance.relative*scale`.
Every input, product, sum and tolerance threshold must remain finite; overflow
returns an error, never a NaN tolerance pass. Large boosts suffer cancellation:
a small scaled residual does not establish small absolute error. Select and
report tolerances for the intended numerical regime. The tests use 1e-13 for
moderate boosts, not a universal accuracy guarantee or exact Lorentz proof.

## Exact and conditional symbolic operations

For rational `u=tanh(eta/2)` with `|u|<1`, the exact constructor computes

\[
c=\frac{1+u^2}{1-u^2},\qquad s=\frac{2u}{1-u^2}.
\]

The real matrix entries are embedded in `BigC`. For `u=1/2`, `(c,s)=(5/3,4/3)`;
the inverse uses `-u`. This parameterization produces an exact rational boost,
**not an exact rational rapidity**. It covers axial boosts; arbitrary rational
unit axes are not part of this constructor.

Symbolic construction requires distinct nonempty labels and retains the explicit
zero condition `c^2-s^2-1`. Names are uninterpreted polynomial variables. The
condition is not proved, installed as a quotient relation or silently reduced.
The matrix is Lorentz only after imposing that condition. No transcendental
parser, exact evaluation of hyperbolic functions, or polynomial algorithm copy
is involved: all polynomial operations use existing checked `BigPoly` methods.

Exact operations accept a shared mutable `Budget`, validate coefficient sizes
before cloning and use existing checked coefficient/matrix primitives. Failures
return no partial output; consumed budget is not rolled back. Budgets bound
logical algebraic work and declared sizes, not elapsed time or allocator usage.

## Verification and ownership

`tests/rapidity.rs` independently checks explicit component values, metric
raising/lowering, inverse, intervals, composition, sign controls, the rational
3-4-5 boost, symbolic condition evaluation and invalid/budget/overflow cases.
Run from the repository root:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-ch2-rapidity cargo test --offline --locked \
  --manifest-path Rust/qft-symbolics/Cargo.toml --test rapidity
```

Specification: `prompts/relativity_rust.prompt`, written before implementation.
This is manual prompt-first implementation; no provider generation, PDD apply,
whole-project regeneration or formal proof is claimed. Module exports and PDD
mapping integration belong to the parent task. Legacy modules and mono-notes
are outside this change's ownership.
