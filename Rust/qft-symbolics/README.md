# qft-symbolics — exact Rust starting point

A small Rust library and JSON CLI for exact QFT calculations and independent checks, started 2026-10-16. This is **not a full Rust translation of
Cadabra2**, not API-compatible with it, and not a replacement for the Python
package yet. Existing Python, Cadabra source, manuscripts and book sources are
unchanged.

## Walk through the checks

From the mathphysics root:

```sh
export CARGO_TARGET_DIR=/tmp/mathphysics-rust-target
cargo run --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --example checks_walkthrough
```

This is a runnable explanation, with real assertions. It checks two-gamma traces
against explicit matrices, translates epsilon conventions, rejects evaluation at
a rational pole, preserves operator order and dummy-index scope, and independently
replays scalar identities with explicit nonzero assumptions. Deliberately wrong
results and missing assumptions must be rejected. It never labels the result a
formal QFT proof.

For the actual Srednicki acceptance checks:

```sh
cargo test --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --test textbook_research -- --nocapture
cargo test --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --test replay_checker
cargo test --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --test typed_expression
cargo test --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --test checked_research
```

Read [SREDNICKI-CHECKS.md](docs/SREDNICKI-CHECKS.md) alongside those tests: it
explains each claim, public API, independent oracle, negative control and meaning
of success. The ten new test functions include 5,461 exact component trace words;
that count is not 5,461 textbook equations or completed book items.

The new public capability modules are [checked arithmetic](docs/CHECKED-ARITHMETIC.md),
[typed structural expressions](docs/EXPRESSION.md), and [independent algebra replay](docs/REPLAY.md).
The replay checker independently expands a scalar ring rather than calling the
production polynomial normalizer, but shares trusted BigC arithmetic. Its conditional
algebra evidence is distinct from the numerical Python comparisons and from formal proof.
The separately locked research consumer exercises all modules without a book,
frontend, server or AI harness.

Mono notes remain theory-only. Implementation mappings and test explanations live
here near the code and use stable mathematical references. No general expression
parser, tensor rewrite replay, gamma5 prescription, OS resource isolation, full CAS
or all-book completion is claimed. Existing tensor expression arithmetic and legacy
transport protocols retain their separately documented limits.

## Lorentz invariance and rapidity

```sh
# From the mathphysics root, keeping build products outside the repository:
export CARGO_TARGET_DIR=/tmp/mathphysics-rust-target
cargo run --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --example lorentz_walkthrough
cargo test --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --test chapter_two -- --nocapture
```

[The Chapter 2 guide](docs/CHAPTER-TWO.md) explains the per-problem premises,
calculations, independent expectations and limits. It distinguishes algebraic
coefficient derivations, representation checks, finite numerical witnesses and
open upstream statements. In particular, the unitary group-law proof upstream
of 2.2 is not automatically reconstructed, and 2.8(f)'s absent-central-charge
statement is not inferred from a field-action calculation.

Reusable modules:

- [relativity](docs/RELATIVITY.md): rapidity/velocity conversion and collinear
  composition, explicit positive-off-diagonal Srednicki boosts, metric/vector
  operations, exact rational boosts and conditional symbolic boost matrices.
- [ordered algebra](docs/ORDERED.md): noncommuting words with exact commuting
  polynomial coefficients, commutators, Jacobi expansion, substitution and
  selected-parameter degree truncation/coefficient extraction.
- [differential operators](docs/DIFFERENTIAL.md): exact differentiation and
  full multivariate Leibniz composition, with Lorentz/momentum representation
  constructors and strict coordinate contexts.

The metric is diag(-1,1,1,1); numerical boost sign follows (2.37).
Rapidity is eta=atanh(beta). Numeric results retain declared tolerances and reject
nonfinite/overflow/saturated velocity results. The exact boost parameter
u=tanh(eta/2) can be rational even when eta is not; u=1/3 gives cosh(eta)=5/4
and sinh(eta)=3/4. Symbolic c/s entries retain c^2-s^2-1=0 as a required condition,
not a fact established from their names. No general transcendental parser,
quantum-domain theorem, unitarity certificate or central-extension proof is claimed.

## Run

From this directory:

```bash
export CARGO_TARGET_DIR="${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-rust/target"
cargo test --offline
cargo run --offline -- demo
cargo run --offline -- schouten
cargo run --offline -- tensor-cases
cargo run --offline -- schouten --wrong-sign  # intentionally FAIL, exit 1
cargo run --offline -- gamma-trace a b c d
cargo run --offline -- smga-fixture 6 0
```

## Checked and arbitrary-precision arithmetic

The optional `BigQ`/`BigC` coefficient types use the same generic polynomial and Clifford engine as the existing small types. Checked APIs return structured errors under explicit work limits. Legacy matrix/tensor arithmetic remains bounded.

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-rust-target cargo run --offline --locked \
  --manifest-path Rust/qft-symbolics/Cargo.toml --example arithmetic
```

See [the foundation audit and evidence](../../docs/plans/rust-physics-platform/ARITHMETIC-FOUNDATION.md) for domains, budgets, independent-client qualification and remaining migrations.

## Reusable library example from the repository

The source is [examples/reuse.rs](examples/reuse.rs). From the mathphysics root:

```bash
CARGO_TARGET_DIR=/tmp/mathphysics-rust-target cargo run --offline --locked \
  --manifest-path Rust/qft-symbolics/Cargo.toml --example reuse
```

This compiles the repository’s library and example directly. Source remains in
mathphysics; Cargo’s executables and intermediate build products stay outside it.
The example calls shared Clifford/spinor operations with fresh inputs and does not
require textbook files, a server or the frontend. The earlier outputs-directory
client is an archived demonstration; this is the maintained example location.

**Complete validation**, including release-mode tests, Clippy, real Cadabra
transformations, and differential parity with the old Python package:

```bash
bash validation/run.sh
```

It uses the existing Python environment at
`../../amplitudes/MHVamplitudes/.venv/bin/python` and existing
`cadabra2-ubuntu: 24.04` image. Override `PYTHON`, `CADABRA_IMAGE`,
`CARGO_TARGET_DIR` or `MATHPHYSICS_RUN_DIR` if needed. No package installation,
downloads, model calls or network access are required; Docker runs with
`--network none`. Build products and validation logs stay in the external cache.

## What is implemented

| Module | Capability | Domain/evidence |
| --- | --- | --- |
| `scalar` | Reduced rational and Gaussian-rational arithmetic | Exact **within checked i128 intermediate range**, not arbitrary precision |
| `polynomial` | Commuting sparse polynomials; rational-function cross-multiplication | Q(i); explicit nonzero-denominator domain |
| `clifford` | Ordered gamma-word normalisation, trace pairing, explicit formal-d sandwich contractions | Symbolic Clifford form h=-eta; labels alone never imply sums |
| `tensor` | Typed spaces/variance, free signatures, hygienic dummies, signed commuting-factor orbits, metric/delta contractions | Bounded independent implementation; see [contract](docs/TENSOR.md) |
| `lorentz` | Evaluated chiral four-slash traces and epsilon contraction | Four dimensions, arbitrary commuting vectors, explicit orientation |
| `matrices` | Pauli, epsilon, sigma/bar-sigma, chiral gamma, Lorentz generators | Exact finite matrices, eta=(-,+,+,+) |
| `twistor` | Srednicki brackets, wavefunctions, momenta, slash, polarization numerators | Polynomial identities in independent phi and phi-star components |
| `spinor` | Component determinant algebra, Parke–Taylor ratio, BCFW shift | Bare color-ordered formula; no coupling/phase normalization inferred |
| `smga` | Closed R1 formula and interval-DP preamplitude/BG recurrence | **Half-collinear (2, 2)**, conserved rational vectors; not ordinary Lorentzian BG |

Start with `corpus.json` for named textbook checks, source anchors, direct
chapter prerequisites, and explicit partial coverage. The current suite has
**54 integration tests plus two compile-fail doctests**, in debug and release.
Thirteen named textbook entries now include formal-d 47.2, angular 50.2a and
fully evaluated 50.4. Live Cadabra agrees on 12 positive/negative tensor cases
and both formal-d gamma identities; a wrong-sign control stays nonzero.
See [extended derivations](docs/TEXTBOOK-EXTENDED.md).
See [derivations](docs/DERIVATIONS.md), [migration audit](docs/MIGRATION.md)
and [verification record](docs/VERIFICATION.md).

## Conventions are part of the API

- `matrices` and `twistor`: eta=diag(-1,+1,+1,+1),
  {gamma^mu,gamma^nu}=-2 eta^mu nu,
  epsilon^12=+1, epsilon_12=-1,
  p_(a,dota)=-phi_a phi-star_dota.
- `clifford`: {gamma(a),gamma(b)}=2 h(a,b). **Set h=-eta** for Srednicki.
  Labels are opaque free slots/component placeholders, not bound indices.
- Generic `spinor::bracket(a,b)` is det(a,b)=a0 b1-a1 b0.
  `twistor::square` is **minus** that determinant of lower phi components;
  `twistor::angle` is the determinant of phi-star components.
  `smga::bracket` is also minus the determinant, explicitly matching Python.
- `twistor::polarization` returns **sqrt(2) times epsilon**. This avoids
  falsely approximating sqrt(2) by a rational. All denominator-domain
  restrictions remain in the returned ratios. q=k is rejected.
- `spinor::parke_taylor` returns <ij>^4 / product <a,a+1>, without i, g,
  color factors or sqrt(2) vertex normalizations. It tests algebraic structure,
  not derivation from Yang–Mills Feynman rules.
- Independent spinors in polynomial/BCFW algebra do not impose a Lorentzian
  reality condition. Real rational SMGA vectors belong to the separate (2, 2)
  calculation. Do not transfer a sign-function formula to real Lorentzian QFT.
- Cyclicity is **not** a gauge Ward identity. Polarization transversality and
  cyclicity are separate checks.

## Agent-facing verdicts

CLI stdout is JSON. Exit 0: PASS (or COMPUTED for a trace, **not a proof**).
Exit 1: a checked residual disagreed (FAIL). Exit 2: unsupported input,
singular domain, exhausted budget, arithmetic overflow or internal panic
(INCONCLUSIVE). Stderr may contain the panic diagnostic.

Library arithmetic operators deliberately panic on unrepresentable
intermediates; callers must catch that boundary or let their process fail.
**Never promote a process failure to verified.** Polynomials have 8, 192-term
and degree-64 limits; gamma operations accept at most 12 factors with a
100, 000-rewrite budget for normal ordering. SMGA supports 3..12 legs (CLI
fixtures 3..10). Complexity still grows rapidly.

An identically zero residual proves an identity only in the stated algebra
and assumptions. Pointwise agreement with Python, or finite multiplicity
recurrence tests, is not an all-n theorem or a novelty claim.

## Not implemented / next gates

1. Arbitrary-precision coefficients and uniformly structured arithmetic errors.
   Tensor validation/work-budget errors already return Result.
2. Scalable/general tensor canonicalisation beyond the implemented bounded
   commuting-factor layer. Plain Clifford labels still never imply sums.
3. General Grassmann algebra, general Fierz/Majorana identities, derivative
   barriers, multiterm and dimension-dependent tensor relations.
4. General gamma/tensor contractions and epsilon/gamma5 regularization
   schemes. Specialized formal-d 47.2 is implemented; gamma5 remains 4D.
5. Full Python-module parity; chapter 50.2b, 60.2b–60.4; nonabelian color algebra,
   three-/four-gluon vertices, ordinary Yang–Mills BG, then independent
   4-/5-point MHV verification and BCFW recursion.
6. Expression parser / arbitrary candidate ingestion, LLM harness, held-out
   evals and an unattended research loop. The CLI is a tool seed, not those.

No original physics is claimed. All shipped identities are known results.

## Provenance and licensing

Repository MIT license applies to newly authored Rust and the port of this
repository's Python calculation. The algebra is independently written from
mathematical definitions; **no Cadabra C++ source, xPerm code or test text has
been translated or vendored**. Cadabra runs only as an external verification
backend. Its source inspection informed scope and missing features.

Local Cadabra checkout inspected:
`67cdebce798dbb2e230ae3246221ca3302a03129`.
Cadabra core is GPL-3.0-or-later (see its source notices); a future translation
of that code needs a separate license decision and retained notices. Changing
implementation language does not remove GPL obligations. This crate does not
assert a legal exemption for arbitrary future FFI or process integrations.
