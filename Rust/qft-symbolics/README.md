# qft-symbolics — exact Rust starting point

A small, dependency-free Rust library and JSON CLI for **verified QFT
calculations**, started 2026-10-16. This is **not a full Rust translation of
Cadabra2**, not API-compatible with it, and not a replacement for the Python
package yet. Existing Python, Cadabra source, manuscripts and book sources are
unchanged.

## Run

From this directory:

```bash
export CARGO_TARGET_DIR="${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-rust/target"
cargo test --offline
cargo run --offline -- demo
cargo run --offline -- schouten
cargo run --offline -- schouten --wrong-sign  # intentionally FAIL, exit 1
cargo run --offline -- gamma-trace a b c d
cargo run --offline -- smga-fixture 6 0
```

**Complete validation**, including release-mode tests, Clippy, real Cadabra
transformations, and differential parity with the old Python package:

```bash
bash validation/run.sh
```

It uses the existing Python environment at
`../../amplitudes/MHVamplitudes/.venv/bin/python` and existing
`cadabra2-ubuntu:24.04` image. Override `PYTHON`, `CADABRA_IMAGE`,
`CARGO_TARGET_DIR` or `MATHPHYSICS_RUN_DIR` if needed. No package installation,
downloads, model calls or network access are required; Docker runs with
`--network none`. Build products and validation logs stay in the external cache.

## What is implemented

| Module | Capability | Domain/evidence |
| --- | --- | --- |
| `scalar` | Reduced rational and Gaussian-rational arithmetic | Exact **within checked i128 intermediate range**, not arbitrary precision |
| `polynomial` | Commuting sparse polynomials; rational-function cross-multiplication | Q(i); explicit nonzero-denominator domain |
| `clifford` | Ordered gamma-word normalisation and trace pairing recurrence | Symbolic symmetric Clifford form h; **no implicit Einstein sums** |
| `matrices` | Pauli, epsilon, sigma/bar-sigma, chiral gamma, Lorentz generators | Exact finite matrices, eta=(-,+,+,+) |
| `twistor` | Srednicki brackets, wavefunctions, momenta, slash, polarization numerators | Polynomial identities in independent phi and phi-star components |
| `spinor` | Component determinant algebra, Parke–Taylor ratio, BCFW shift | Bare color-ordered formula; no coupling/phase normalization inferred |
| `smga` | Closed R1 formula and interval-DP preamplitude/BG recurrence | **Half-collinear (2,2)**, conserved rational vectors; not ordinary Lorentzian BG |

Start with `corpus.json` for named textbook checks, source anchors, direct
chapter prerequisites, and explicit partial coverage.
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
  reality condition. Real rational SMGA vectors belong to the separate (2,2)
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
**Never promote a process failure to verified.** Polynomials have 8,192-term
and degree-64 limits; gamma operations accept at most 12 factors with a
100,000-rewrite budget for normal ordering. SMGA supports 3..12 legs (CLI
fixtures 3..10). Complexity still grows rapidly.

An identically zero residual proves an identity only in the stated algebra
and assumptions. Pointwise agreement with Python, or finite multiplicity
recurrence tests, is not an all-n theorem or a novelty claim.

## Not implemented / next gates

1. Arbitrary-precision coefficients and structured resource-error types.
2. Typed Lorentz/dotted/undotted indices, free-index signature validation,
   dummy renaming, symmetric/antisymmetric tensor canonicalisation, and
   flat metric/delta elimination. Do **not** pass a gamma label named "mu"
   and assume contraction occurs.
3. General Grassmann algebra, general Fierz/Majorana identities, derivative
   barriers, dimension-dependent tensor relations.
4. d-dimensional gamma contractions, epsilon/gamma5 regularization schemes;
   the current gamma5 checks are strictly four-dimensional.
5. Full Python-module parity; chapter 50.2, 60.2b–60.4; nonabelian color algebra,
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
