# Exact polynomial differential operators

`qft_symbolics::differential` provides reusable representations, with no
textbook problem numbers in production APIs. All constructors and algebra use
`computation::Result` and a final mutable `Budget`. Coefficients are `BigPoly`
over Gaussian rationals; floats, infinity and NaN have no representation.

`differentiate(&BigPoly, variable: &str, budget)` computes an exact first
partial. Every other polynomial symbol is a constant parameter. Multiplicity
is exact, including repeated coordinate occurrences and arbitrary precision
rational/complex coefficients.

`DifferentialOperator<const N: usize>` stores private ordered coordinates
`[String; N]` and a private `BTreeMap<[u16; N], BigPoly>`. A key `alpha`
represents ordinary lower partials, with all derivatives to the right:

```text
A = sum_alpha a_alpha(x) D^alpha.
```

The borrowed `coords()` and `terms()` accessors expose the normalized result.
Empty/duplicate coordinates, zero dimensions, dimensions above 16 and invalid
axes fail. Coordinate order and spelling must match exactly for binary
operations; no implicit renaming or permutation occurs. Additional coefficient
symbols are allowed and differentiate as constant parameters. Zero coefficients
are removed, including explicit zeros passed to `from_terms`.

The frozen API is:

```text
from_terms(coords, BTreeMap<[u16; N], BigPoly>, budget) -> Result<Self>
partial(coords, axis: usize, budget) -> Result<Self>
multiplication(coords, BigPoly, budget) -> Result<Self>
try_add / try_sub / compose / commutator(&self, &Self, budget) -> Result<Self>
try_scale(&self, &BigPoly, budget) -> Result<Self>
apply(&self, &BigPoly, budget) -> Result<BigPoly>
coords(&self) -> &[String; N]
terms(&self) -> &BTreeMap<[u16; N], BigPoly>
```

Scaling is polynomial multiplication on the left. `A.compose(&B, budget)`
means apply B first, then A. It uses the complete multivariate Leibniz rule:

```text
(a_alpha D^alpha)(b_beta D^beta)
  = sum_gamma<=alpha binomial(alpha,gamma)
      a_alpha (D^gamma b_beta) D^(alpha-gamma+beta).
```

The multiindex binomial is the product of the coordinate binomials.
`commutator` computes AB−BA in this same normal form, retaining all intermediate
orders until exact cancellation. It does not merely compare polynomial samples.

## Minkowski constructors

For `DifferentialOperator<4>`, coordinates are upper `x^mu`, indexed 0 through 3.
With `eta = diag(-1,1,1,1)`, raised derivatives are
`partial^mu = eta[mu,mu] D_mu`. The constructors take owned `hbar: BigQ`:

```text
momentum(coords, mu, hbar, budget):
  P^mu = -i hbar partial^mu
orbital_lorentz(coords, mu, nu, hbar, budget):
  L^{mu nu} = (hbar/i)(x^mu partial^nu - x^nu partial^mu).
```

In particular, P^0=+i hbar D_0 and P^j=−i hbar D_j. The signs give:

```text
[L^{mu nu}, P^rho]
  = i hbar (eta^{mu rho} P^nu - eta^{nu rho} P^mu)
[L^{mu nu}, L^{rho sigma}]
  = i hbar (eta^{mu rho} L^{nu sigma} - eta^{nu rho} L^{mu sigma}
            + eta^{nu sigma} L^{mu rho} - eta^{mu sigma} L^{nu rho})
[P^mu, P^nu] = 0.
```

These are identities in the formal algebra of commuting partials and polynomial
coefficients, hence apply to arbitrary smooth scalar fields with those
commuting partials. Polynomial witnesses additionally check executable action.
This does not establish general quantum operator domains, self-adjointness,
Hermitian/unitary group representations or absence of central extensions in
other representations.

## Bounds and failures

All input operands are validated even for zero results. The hard dimension cap
is 16; the hard total derivative order cap is 32. Both `Degree` and `WordLength`
also limit derivative order. Thus the default WordLength=12 limits operators to
order 12 unless explicitly raised. Composition conservatively rejects a pair
whose sum of orders exceeds a limit, even if its eventual commutator cancels.
The representation and enumeration are flat and iterative: no expression-tree
recursion or unbounded recursive depth is introduced.

`Degree` also bounds coefficient monomial lengths. `Terms` bounds each input,
intermediate and output polynomial; for each operator it additionally bounds
operator entries plus all coefficient monomials. `InputBytes` bounds aggregate
UTF-8 coordinate bytes plus all coefficient variable occurrences in each
operator, and the existing polynomial aggregate separately for polynomial
inputs/intermediates. Explicit zero entries in an input map are validated before
normalization. No limits count combined live heap usage or cumulative allocation.

`Operations` charges validation, iterative enumeration, cloning and exact
coefficient work; exact arithmetic reuses bit-weighted checked coefficient
operations. `CoefficientBits` checks inputs and conservative arithmetic growth.
Known resource failures precede owned key construction/insertion; local
polynomial intermediates obey their own existing bounds before operator insertion.
Cancellation may fail under conservative budgets. Failed calls may consume work
budget but never mutate borrowed operators or polynomials. Allocation OOM and
wall-clock time are outside the deterministic logical budget contract.

## Verification

`tests/differential_operators.rs` covers product and mixed derivatives, explicit
multivariate binomial coefficients, large rational coefficients, composition
versus successive action, cancellation of second partials, all 256 Lorentz,
64 Lorentz–momentum and 16 momentum brackets, polynomial witnesses, metric signs,
context mismatches and resource/invalid-input cases. Run:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-ch2-differential cargo test --offline --test differential_operators
```

This unit was implemented manually from its prompt. Passing tests and prompt
registration do not claim provider regeneration or whole-project PDD sync.
