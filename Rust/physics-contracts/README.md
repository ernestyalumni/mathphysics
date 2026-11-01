# Mathematical item and evidence contracts

A Rust library and JSON CLI for the first shared contract layer of the mathematical physics platform. It validates source/item identities, declared assumptions and conventions, provenance metadata, evidence links and consistent evidence classifications.

**`record_valid` means metadata is consistent. It does not mean a mathematical result is correct, a source was certified, a digest was authenticated, or an artifact/proof was checked.** Every CLI report carries this boundary. The implementation never reads referenced artifacts or source files.

## Run

From the repository root, place build products outside the source tree:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-contracts-target cargo test --offline --locked --manifest-path Rust/physics-contracts/Cargo.toml
CARGO_TARGET_DIR=/tmp/mathphysics-contracts-target cargo run --offline --locked --manifest-path Rust/physics-contracts/Cargo.toml -- validate Rust/physics-contracts/examples/valid-bundle.json
```

The CLI also accepts `validate -` or `validate` to read stdin. It emits JSON and exits 0 for consistent records, 2 for malformed/invalid input or I/O failure. Input is limited to 1 MiB before decoding. The library supports typed Rust callers through `validate(&Bundle)` or bounded JSON through `validate_json(&[u8])`.

Dependencies are pinned serde and serde_json. No Python, C++, external CAS, provider API, GPU or network call is needed to execute validation. Offline builds require the recorded Cargo packages to exist locally.

## Record model

- IDs distinguish `equation`, `exercise`, `prose`, `diagram`, `example`, `definition` and `theorem`. Srednicki chapter numbers are 001–097; explicit `frontmatter` and `backmatter` locations support other source obligations.
- Sources have declared SHA-256 digests and one-based PDF/text locators. `candidate`, `text_checked` and `visually_checked` are declared review states; this validator does not certify them.
- Items identify an explicit statement/domain, a complete assumption set, a convention, an implementation state and evidence links.
- Evidence records bind to exactly one item with matching statement, domain, assumptions and convention. Numerical evidence requires a precision/error contract. Exact, formal, replay, external and source evidence have distinct constraints.
- Claimed formal/replay evidence needs checker/version/digest metadata linked to a declared artifact. This is structural consistency, not execution of that checker.
- Duplicate IDs, unknown references, conflicting evidence links and duplicate artifact paths are rejected. Arithmetic exhaustion or unsupported computations belong to an inconclusive result with a reason; metadata validation cannot make those computations pass.

JSON schemas live in `../../docs/schemas`. Schemas validate shapes and some local constraints; the Rust validator additionally checks registry membership, cross-record equality/backlinks, aggregate budgets and UTF-8 byte limits. Prefer the CLI for complete metadata checks. The examples are synthetic fixtures, not evidence for textbook equations.

Resource fields in evidence describe intended execution budgets. They do not implement process or GPU isolation. The typed caller owns preexisting allocations; the validator rejects excessive record counts before descending through oversized collections.

## PDD and verification

`architecture.json` maps two versioned prompts to `src/lib.rs` and `src/main.rs`. `docs/accepted-request.txt` and `docs/intent-plan.json` preserve the approved F01 scope. Initial code was implemented by the current coding agent; a full automated PDD regeneration is not claimed. Independent tests, JSON schemas and source corpora are outside generated outputs.

Tests cover real CLI invocation, bad sources, duplicate/unknown references, inconsistent evidence links, numerical overclaims, fake formal metadata, input/resource limits and nonfinite tolerances through the Rust API. An independent review prompted fixes for approximate external-comparison metadata, duplicate artifact identities, front/back-matter IDs and early resource rejection.
