# Checked reusable capabilities

`src/checked.rs` is additive. Parent integration exports `pub mod checked`;
legacy types and callers remain unchanged. No new dependencies are required.
The prompt is `prompts/checked_rust.prompt`; independent acceptance tests are
`tests/checked_research.rs`. This is manual prompt-first implementation, not
provider-backed generation or verified regeneration.

## Interfaces

- `Matrix<A, N> = [[A; N]; N]`, with borrowed-input `identity`, `add`, `sub`,
  `mul`, `scale`, `trace`, `transpose`, all returning computation `Result`.
  Matrix entries use `Coefficient`, supporting bounded `C` and exact `BigC`.
- `sigma`, `barsigma`, `gamma` return checked bounded `C` matrices and reject
  indices outside 0..4. Convention: eta=(-,+,+,+), {gamma_mu,gamma_nu}=-2eta_mu_nu.
  These constructors use checked negation and no legacy operator path.
- `RationalFunction<A>::new`, `numerator`, `denominator`, `try_add`, `try_sub`,
  `try_neg`, `try_mul`, `try_div`, `try_residual`, `try_eval` use the existing
  generic polynomial engine. Fields are private; the denominator is nonzero.
  No factor cancellation occurs. A zero numerator does not remove its poles.
  For (n/d)/(m/e), use n*e*e/(d*m*e): the extra e on both sides retains the
  original divisor's domain restriction. Reject an identically zero divisor.
  A zero cross-product residual establishes equality only on the common
  nonzero domain, not equality of domains. Evaluation checks denominator first.
- `Spinor<A>`, `bracket`, `schouten` operate on commuting polynomial components;
  bracket(a,b)=a0*b1-a1*b0. No Grassmann or conjugacy assumptions are imposed.
- `bcfw_shift<A>` acts on coefficient-valued spinors with [i,j> convention:
  lambda_i += z lambda_j, tilde_j -= z tilde_i. Validate indices and compute
  all four updates before mutation; every returned error leaves inputs intact.
- `TensorValidator::new(orbit_steps)`, `orbit_steps`, `context`, `declare`,
  `require_symmetric_metric` delegate existing Context validation.
  `TensorError::Validation(tensor::Error)` retains structured native errors;
  `TensorError::Computation(computation::Error)` distinguishes shared limits.
  The wrapped error is available through `std::error::Error::source`.

## Limits and boundaries

Every arithmetic call receives the caller's shared mutable Budget. Matrix
entry count is bounded by Terms, coefficient bits are checked before cloning,
and multiplication reserves N^3 logical work before allocating output.
Polynomial term/degree/name-byte/bit/work limits remain owned by Polynomial.
BCFW validates only the selected legs it reads and modifies, not untouched legs.
Tensor declaration preflights closure work and transient/retained storage through
Budget before invoking Context. Rank is at most eight; group bounds are 1 for
None, n! for built-in symmetric/antisymmetric groups and 2*n! for arbitrary
signed generators. Negative identity generators can produce both signs.

`retained_bounds()` reports successful declarations, conservative group entries,
and conservative logical storage bytes. Combined declarations + group entries
consume Terms; cumulative names, slots and group storage consume InputBytes.
These private counters survive renewed Budget instances. Transient closure
queue/set/output entries and generator clones receive additional Terms/byte
checks. Work reserves 8*g*(g+1)*(k+1)*(n+1), with group bound g, generator count
k and rank n, before any closure expansion; this deliberately conservative
comparison bound can reject large groups even when the native orbit cap permits
them. Each entry reserves n*sizeof(usize)+256 logical storage bytes plus
512 per declaration and retained name/slot bytes; transient storage reserves
four declaration footprints plus generator entries. These are declared storage
units, not measurements or guarantees of allocator overhead.

Only successful Context declarations commit counters; any validation, shared
limit or orbit failure leaves registry/counters unchanged. Budget work already
spent does not roll back. Its separate `orbit_steps` remains an additional
per-declaration cap. Caller-owned vectors already exist on entry.

Limits bound declared algebraic work and component sizes, not total heap,
allocator failure, stack size, elapsed time or OS resources. Coefficient trait
implementations must obey their arithmetic/default/zero contracts. Generic
matrix/spinor algebra assumes commuting scalar coefficients; rational-function
interpretation assumes the usual integral coefficient domains C or BigC.
Conservative intermediate limits can reject expressions that would cancel.

Tensor expression construction, contraction and arithmetic are deliberately
not wrapped: their existing implementation invokes legacy polynomial operators.
This wrapper provides typed declaration/metric validation, not a checked tensor
CAS migration. Sigma/gamma constructors are bounded C; general matrix arithmetic
supports BigC. No chapter solvers, claim certificates or formal proofs are added.

## Verification

Run after parent export:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-next-arithmetic cargo test --offline --locked --test checked_research
```

Nine independent tests exercise literal chiral matrices and the Clifford form,
a 2^128 matrix entry with a separately specified exact square, rational poles
and quotient domains, symbolic Schouten plus an independent determinant value,
integer momentum conservation and atomic BCFW failure, typed tensor invalid /
unsupported / orbit-budget cases, and bit/degree/byte/term/operation exhaustion.
Before parent export, the same source/tests were compiled through a temporary
external harness reexporting the existing crate primitives. That harness checks
implementation behavior but is not evidence that the parent export is integrated.

Independent review R1/R2 correction: regression tests reproduce the exact
rank-eight symmetric closure with operations=9, terms=1, input_bytes=1;
isolate operations exhaustion with relaxed storage; check signed identity
closure's two-sign group bound; and verify renewed budgets cannot bypass
aggregate registry Terms/bytes. Native duplicate/orbit errors also preserve
registry and counters. No legacy tensor implementation was changed.
