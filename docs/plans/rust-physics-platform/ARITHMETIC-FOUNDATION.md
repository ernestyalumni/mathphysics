# Reusable arithmetic foundation — 2026-10-31

This delivery implements a bounded slice of LIB01/LIB02/F06/F07. It does not complete the Srednicki corpus or a general CAS. Textbook chapter numbering is not part of the production API.

## Capability ownership and conventions

| Capability | Reusable owner | Adapter or independent evidence |
| --- | --- | --- |
| Small rational/Gaussian arithmetic | qft-symbolics::scalar::{Q,C} | Existing matrix, twistor, spinor and tensor modules remain bounded callers |
| Arbitrary-precision rational/Gaussian coefficients | qft-symbolics::exact::{BigQ,BigC} | Explicit opt-in; no silent replacement of Copy types |
| Shared polynomial algorithms | Polynomial<A: Coefficient>, aliases Poly/BigPoly | Same algorithm for both coefficient domains |
| Formal Clifford algebra | qft-symbolics::clifford checked generic functions | Formal anticommutator 2h, trace(I)=4; no implicit sums; existing h=-eta bridge retained |
| Pure SU(2) observations | lie-groups::su2_observation::observe | physics-gateway parses/serializes; independent Pauli matrix and HTTP tests |
| Tool transport | physics-tools / physics-gateway | Legacy protocol unchanged; older kernel_panic errors retained for compatibility |
| Source fixtures | qft-symbolics::tensor_cases, corpus and tests | Verification surfaces; not a separate production symbolic solver |

SU(2) observation explicitly converts the observatory plus-i Pauli coordinates to the existing quaternion convention. Numerical unitary/determinant residuals are not formal proofs. Existing group element constructors retain their own documented tolerances.

## New public surface

`computation::{Error,Resource,Limits,Budget}` provides typed errors and explicit deterministic work limits. Scalar checked operations, polynomial try_* operations and checked Clifford operations return Result. `try_trace_exact` and `BigPoly` use arbitrary-precision coefficients through the same generic algorithms. The CLI adds `rational add|sub|mul|div n d n d` and `gamma-trace-checked labels...`, with structured failures and formal_proof=false.

The small backend uses checked i128 intermediates and can reject intermediate overflow even when cancellation would make the final mathematical value small. Legacy infallible operations remain bounded and may panic. Ratio, matrix, tensor and downstream transport APIs have not all migrated to arbitrary precision or one unified error type. Those follow-up migrations remain explicit work.

Budgets bound declared work and per-object algebraic size, not wall time or total live heap. InputBytes is aggregate symbol bytes per polynomial (including repetitions), separately for input/intermediate/output objects and evaluation names; it is not a request-wide allocation counter. BigQ parsing bounds combined input text. CoefficientBits measures the largest component numerator/denominator. Conservative intermediate bounds may reject cancellable expressions. Default work limits need adjustment for larger exact calculations. OS allocation failure recovery and arbitrary third-party coefficient behavior are not guaranteed.

## Independent reuse qualification

`Rust/qft-symbolics/tests/consumers/research-client` is a separately locked Rust package depending on both public libraries by path. It uses fresh non-textbook symbols, hard-coded decimal rational expectations, independently expanded Pauli entries and failure controls, without HTTP, browser, Python, AI or chapter dependencies. It is run by the external_consumer integration test. The book and other papers can call these same APIs; this fixture is not a claim to validate an arXiv paper.

## Run from the repository

```sh
export CARGO_TARGET_DIR=/tmp/mathphysics-rust-target
cargo run --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --example arithmetic
cargo run --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml -- rational add 340282366920938463463374607431768211456 3 2 3
cargo run --offline --locked --manifest-path Rust/lie-groups/Cargo.toml --example observation
cargo test --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml
```

Sources stay in mathphysics; generated executables stay in the external target directory.

## Dependencies and adoption

The opt-in exact backend pins pure-Rust num-bigint 0.4.6, num-rational 0.4.2 and num-traits 0.2.19 (MIT/Apache-2.0). No GMP/C dependency is introduced. serde_json is a development-only CLI assertion dependency in qft-symbolics. Existing consumers' lockfiles include the new transitive arithmetic dependencies. This is a deliberate pinned baseline, not a claim of newest versions or GPU acceleration.

The reviewed request and read-only intent plan are retained in [foundation-intent](foundation-intent/); the pre-change characterization copy/log is `/tmp/mathphysics-foundation-baseline*`. Seven bounded qft units have versioned prompts and an acyclic ownership graph. Implementation was manual; no provider-backed intent apply, successful regeneration or full PDD sync is claimed. Independent tests/examples and unmapped modules remain outside generated ownership.

Final verification evidence is recorded alongside this document after integration. Independent resource review: `/tmp/mathphysics-foundation-review.txt`; its allocation sequencing findings were assigned for correction before acceptance.

## Final verification receipt

- qft-symbolics: 84 unit/integration tests plus two doctests in both debug and release, including the external-consumer and allocation-observer regressions. Formatting and all-target Clippy pass.
- lie-groups: 57 unit/integration tests plus four doctests; physics-contracts: 19; physics-tools: 11; physics-gateway: eight. All pass debug/release and all-target Clippy. Across five crates: 185 per profile including doctests.
- Existing browser suite: 30 pass against the refactored gateway, covering standalone Explore and the separate notebook/PDF interface. Current repository tests ran from an external staging directory with the existing pinned Playwright installation; the unchanged frontend build was reused.
- Live independent Python parity: 200 cases pass. Cadabra Clifford/Schouten, tensor cases and formal-d sandwich identities pass, including deliberately wrong signs that remain nonzero.
- Seven prompt preprocesses, ownership paths and acyclic dependency checks pass without model invocation. Both resource-review allocation findings are corrected, with a 2 MB allocation-observing regression and positive control.
- Repository-owned arithmetic and SU(2) examples run successfully. Sources stay in the repo and targets/logs stay outside it.

Local evidence: `/tmp/mathphysics-foundation-validation/`, `/tmp/mathphysics-foundation-validation.log`, `/tmp/mathphysics-foundation-regression.log`, `/tmp/mathphysics-foundation-browser.log`, `/tmp/mathphysics-foundation-pdd.log`. The validation directory includes source hashes/provenance. Initial formatting and browser-runner setup failures were corrected before these final passing runs; no solver expectation was changed to manufacture a pass.

Follow-up work: unify remaining typed matrix/tensor/transport errors as separately reviewed migrations; adopt exact coefficients where needed; specify typed expressions/domain semantics and replay checking before expanding reusable research capabilities. Existing 13 textbook checks and chapter inventories remain partial. No unattended workers or all-book completion are implied.

Final dependency cleanup removes the scalar-to-exact MIN-normalization dependency using native unsigned GCD reduction. QFT debug/release each still pass 86 tests including doctests, plus formatting and all-target Clippy: `/tmp/mathphysics-foundation-final-native-min.log`. The exact prompt clarifies that From<C> belongs to BigC. Independent final allocation probes confirm both review findings resolved: `/tmp/mathphysics-foundation-su2/review-probes-final.log`.
