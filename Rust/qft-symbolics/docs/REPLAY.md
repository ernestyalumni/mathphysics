# Independent bounded algebra replay

`qft_symbolics::replay` checks exact identities over commuting scalar symbols
with Gaussian-rational coefficients. It independently expands the scalar AST
into a private sparse `BTreeMap` of sorted monomials and `BigC` coefficients.
It deliberately does **not** call the production polynomial or rational-function
normalizers. Trusted code still includes the shared `BigC` arithmetic, Rust,
and the expression validator/canonical serializer. This is a second algebra
implementation sharing a coefficient core, not an independently verified proof
kernel or a Lean/formal QFT proof.

## Public API and runnable integration

```rust
use qft_symbolics::{
    computation::Budget, expression::Expr,
    replay::{replay, Certificate, ReplayStep, Rule,
             CERTIFICATE_VERSION, SCALAR_DOMAIN},
};
let x = Expr::Symbol("x".into());
let before = Expr::Quotient(Box::new(x.clone()), Box::new(x.clone()));
let after = Expr::Constant(1.into());
let certificate = Certificate {
    version: CERTIFICATE_VERSION,
    domain: SCALAR_DOMAIN.into(),
    start: before.clone(),
    end: after.clone(),
    steps: vec![ReplayStep { before, after, rule: Rule::RationalIdentity }],
    nonzero_assumptions: vec![x],
};
let report = replay(&certificate, &mut Budget::default()).unwrap();
assert_eq!(report.steps_checked, 1);
assert_eq!(report.conditional_domain.len(), 1);
assert!(!report.formal_proof);
```

`replay(&Certificate, &mut Budget) -> computation::Result<ReplayReport>` uses
one caller-owned budget throughout. `ReplayReport` has public fields
`steps_checked: usize`, `conditional_domain: Vec<String>`,
`statement: &'static str`, and `formal_proof: bool` (always false on success).
Failure is a typed computation error, never a partial-success report.
`UnsupportedDomain` also describes certificate/chain/rule rejection; arithmetic
and resource failures retain their shared error variants.

Version is `1`; the exact domain string is
`gaussian-rational-commuting-scalars/explicit-nonzero-poles/v1`. Changing either
rejects the certificate. These constants fix coefficient, commutation and pole
conventions; there is no automatic conversion to another convention.

## Rules, chain binding and assumptions

`PolynomialIdentity` accepts `Constant`, `Symbol`, `Sum`, and `Product` and
compares independently expanded sparse polynomials. Quotients are rejected
under this rule. Empty sums/products mean zero/one. `RationalIdentity` also
accepts `Quotient`, constructs uncancelled numerator/denominator pairs and
compares exact cross products. Every expression is validated, including unused
assumptions, endpoints and every intermediate side. Ordered products, operators,
tensors and contractions/binders are rejected, even when their resulting
signature could be scalar. No implicit operator commutation is permitted.

Each original quotient denominator must appear in `nonzero_assumptions`, matched
by `Expr::canonical` syntax. Even a constant denominator must be explicitly
supplied. Algebraically equivalent but structurally different assumptions need
separate entries; no premise implication engine is provided. All assumptions
are validated and independently expanded; normalized-zero assumptions and
normalized-zero denominators reject. Nested rational assumptions need their
own denominator assumptions. A nonzero assumption is a supplied hypothesis,
not a fact proved by replay. All supplied assumptions, including unused ones,
are sorted/deduplicated into the report's conditional domain. Joint consistency
of nonzero hypotheses is not a separately certified result.

Poles survive cancellation and zero masking. `x/x = 1` requires `x != 0`;
`0 * (1/x) = 0` still requires it. `1/(1/x) = x` requires both `x != 0` and
`(1/x) != 0`. A denominator `x + (-1)*x` rejects even if supplied as nonzero.

At least one step is required; an identity certificate needs an explicit step.
The canonical syntax of `start` must match the first `before`, every `after`
must match the next `before`, and the final `after` must match `end`.
Canonical syntax sorts commuting children but does not simplify/flatten them.
These syntax comparisons only bind the chain. Algebraic validity comes from
independent expansion of the actual before/after expressions; there is no
comparison against an author-supplied normalized expected answer.

## What tampering tests establish

`tests/replay_checker.rs` checks `(x+2)^2 = x*x+4*x+4`, `i*i = -1`, and
`1/2+1/3 = 5/6` using independently specified exact expected values. Replacing
`4*x` by `3*x` fails exact coefficient comparison. Forging `end` or inserting
a syntactically disconnected intermediate fails chain binding, even if that
replacement is algebraically equal. Removing a denominator assumption fails
domain checking; changing domain/version fails convention binding; swapping
ordered operators is unsupported. These are exact symbolic checks, with no
floating-point sampling or numerical witness.

Run the standalone integration suite:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-next-replay cargo test --offline --test replay_checker
CARGO_TARGET_DIR=/tmp/mathphysics-next-replay cargo test --offline --release --test replay_checker
```

## Resource boundary and evidence limits

A pre-scan bounds aggregate AST nodes and conservative input bytes over the
entire certificate, then validates every input. AST depth is capped at 128 even
if user limits are huge. `max_degree` also bounds traversal depth, so conservative
rejection can occur for shallow algebra with verbose structure. Sparse operations
preflight terms, monomial degree, conservative bytes, dense pair count and
coefficient work before cross-product growth. Checked `BigC` operations bound
intermediate coefficient bits/work. Syntax, sparse traversal and comparison
work consume the same budget; repeated replay consumes its remaining allowance.
Cancellation may not rescue a certificate that exceeds a preflight bound.

These limits constrain declared algebra work and conservative representations,
not OS allocation/RSS or elapsed time. The caller constructs the AST before
replay; this module supplies no parser or untrusted-input deserializer and does
not bound construction/destruction outside the checker. Huge caller budgets
can permit expensive wide computations despite the hard recursion cap.

A successful report means independently replayed bounded algebra conditional
on the supplied nonzero domain. It does not establish tensor identities,
operator relations, analytic continuation through poles, source fidelity,
physical assumptions, or a formal proof. The prompt is authored before code;
provider-backed regeneration and whole-project PDD synchronization are not
claimed. Parent integration owns module export and external-consumer checks.
