# Mathematical physics platform execution

The user approved starting the plan and spawning subagents on 2026-10-31. The first implementation wave has delivered the foundation below; final review and evidence are recorded here. This is execution of the approved program, not a claim that the full CAS or all textbook mathematics is implemented.

## Delivered foundation

- `Rust/physics-contracts`: typed Rust metadata library and bounded JSON CLI. It checks source/statement/domain/reference/evidence consistency while explicitly refusing to equate metadata validity with mathematical truth, source certification or replayed proof.
- `docs/schemas`: four portable JSON Schema 2020-12 documents corresponding to the Rust contract. Rust adds semantic cross-record checks.
- `docs/capabilities/registry.json`: 46 initial capability families, current partial coverage, acceptance examples, dependencies and optional reference backends. These are specifications, not 46 implemented engines.
- `source-review/page-map.json`: 97 chapter starts checked against PDF text and printed folios; +24 PDF/printed offset independently observed at each start. Only chapter12's heading required a visual extraction correction. This is not mathematical certification of 97 chapters.
- `source-review/quarantine.json`: all 72 suspicious OCR records accounted for and retained in quarantine; 766 historical resolutions cross-reconciled, without retroactively certifying their mathematics or historical execution.
- `items/chapter-003.jsonl`: initial chapter3 inventory, 102 items including 38 numbered equations, 12 intermediate steps, 35 prose/definition/example items and 17 exercise records. Every item is unassessed for implementation and has no evidence claim. Independent inventory review is tracked separately.

Additional first-wave delivery:

- `Rust/physics-tools`: native Rust JSON discovery and caller-supplied gamma-word trace bridge to the existing bounded exact kernel. Eleven tests pass in debug/release, with formatting and Clippy. An independent Codex subagent made 162 actual CLI probes with fresh labels, independent matrices, parser/budget negatives and a real expansion failure; all expected outcomes passed. This qualifies one Codex exec route, not all harnesses or MCP.
- `docs/capabilities/conventions.json`: 23 versioned convention records and 13 assumptions, with source anchors and incompatible-interpretation warnings. An explicit alias map supplies versioned IDs accepted by F01. Six new cross-module convention tests pass in debug/release and Clippy.
- `items/chapter-035.jsonl`: 61 initial source items (29 numbered equations, nine intermediate steps, 19 prose/definition/example records and four exercises). The author text-read and visually inspected all five pages; independent inventory review remains pending. Existing exercise35.1/35.2 checks remain a separate partial crosswalk.

Chapter3's independent source review confirmed all38 numbered equations and all five exercises/subparts, including3.35. Review refinements add explicit nonrelativistic mass restrictions, clarify the staged field-reality condition, and preserve a fermionic constant claim as unresolved rather than endorsing cancellation. Both chapter bundles pass metadata validation, with every implementation state unassessed and zero mathematical evidence claims.

## Verification

Fresh existing-code baseline: qft-symbolics passes 55 integration tests plus two compile-fail doctests in debug/release; lie-groups passes three unit, 48 integration and four doctests in debug/release. The unmodified qft validation runner passes live Cadabra transformations, negative controls and 200 Python parity cases. The baseline source snapshot is unchanged by those runs.

physics-contracts passes 19 integration/CLI tests in debug/release and Clippy with warnings denied. Independent review reproduced and resolved four findings: approximate external evidence without numerical budgets, conflicting artifact identities, front/back-matter scope, and delayed typed-record budget rejection. Thirty independent JSON cases and six typed Rust checks passed their expected outcomes after fixes. Schema tests pass 14 test groups including 432 cross-field combinations.

Scoped PDD preprocessing passes. Corrected path checks resolve the actual Rust code and test files with no unsafe path failures; dry-run state is UNBASELINED. Full automated PDD intent application, a generated story contract, certification and fresh regeneration are not claimed. A request to implement the approved plan does not require repeating approval for paths/commands. Each bounded unit preserves its reviewed request and prompt graph.

## Dispatch and continuation

`tasks.json` is the central task list. The user approval remains effective for the program; agents inspect actual dependencies, claim exclusive paths, characterize existing behavior and maintain bounded prompts without making the user choose mechanical details. Parallel workers produce isolated deliverables; the integrator reviews and copies them to the repository. Runtime logs, page images and build products stay outside the source tree.

The first checked arithmetic and opt-in arbitrary-precision slices (F06/F07) are delivered; remaining migrations are recorded in [ARITHMETIC-FOUNDATION.md](ARITHMETIC-FOUNDATION.md). Next core dependencies include typed expressions and domain semantics (F08/F09), replay checking (F10), and qualified agent transport (F12–F15). Continue chapter inventories against the shared contracts, selecting both foundational chapters and existing Rust coverage. An initial bounded tool bridge does not mean the general expression/checker prerequisites are complete.

No new harness accounts/provider credentials, package installs, GPU jobs, external posts, commits or pushes are part of this wave. Subagent identities and evidence locations are in the first-wave verification receipt. Do not infer unattended work persists merely because tasks are queued.

## Use the first callable Rust tool

From the repository root:

```sh
printf '%s\n' '{"command":"gamma_trace","labels":["mu","nu"],"domain":"formal_symmetric_h","budget":1}' | CARGO_TARGET_DIR=/tmp/mathphysics-tools-target cargo run --offline --locked --quiet --manifest-path Rust/physics-tools/Cargo.toml
```

This returns a computed exact expression `4 h(mu,nu)` in the declared formal algebra, with `formal_proof=false`. Use `{"command":"discover"}` to inspect the actual supported surface. Word length, pairing budget, arithmetic and expansion bounds remain explicit. General expression parsing, transport adoption of the opt-in arbitrary-precision backend, certificate replay, MCP and remaining harness integrations remain queued.

## Unified UX story and parallel Srednicki refinement

The user requested the third [unified notebook story](../../../user_stories/story__unified-notebook-workspace.md). Browser is primary; Linux/macOS dedicated apps, shared notebook/TypeScript frontend/Rust services, closely Mathematica-style and readable SNES pixel themes, fitting sound controls and every named harness are requirements. UX01–UX10 add concrete dependencies without replacing the Rust/Srednicki work.

Reference source acquisition downloaded pinned SNES.css, selected CSSwitch files, a pinned pixel font with its notice and official notebook documentation. Provenance and101 hashes are in docs/design/notebook-references; downloaded code was not executed and no frontend was published. CSSwitch source repository/license was not established, so those files remain visual/reference material rather than a product dependency.

In parallel an independent worker checked all five original chapter35 pages, accepted the61-item initial inventory and confirmed29 numbered equations/four exercises. Refinements make generator premises explicit and give all28 unnumbered items batch ownership. These remain source/specification work, not completed mathematical proofs or chapter implementation.

## Reusable arithmetic delivery — 2026-10-31

[Foundation audit, public APIs and reproducible commands](ARITHMETIC-FOUNDATION.md): shared checked/generic polynomial and Clifford arithmetic, opt-in arbitrary-precision coefficients, typed SU(2) library extraction and an independent Rust client. Five crates pass 185 tests per debug/release profile (including doctests), all-target Clippy passes, 30 browser tests pass, and the live Cadabra/200-case Python comparison passes with negative controls. This advances F06/F07 without claiming all legacy modules or transports migrated. Seven bounded PDD prompts preprocess without model calls; full regeneration remains unverified.

## Explainable QFT research batches — 2026-10-31

All five approved workstreams now have a bounded implementation: checked generic
matrices/rational functions/spinors and declaration validation; explicit epsilon
profiles; ten additional individually assessed Chapter 35/47/50 research checks;
a typed structural-expression/domain slice; and independent scalar algebra replay.
The source test IDs remain fixtures rather than production solver entry points.

Start with the repository-owned `checks_walkthrough` example described in
[the library README](../../../Rust/qft-symbolics/README.md#walk-through-the-checks).
It performs actual assertions while explaining two-gamma traces, convention
conversion, rational poles, ordered operators/binder hygiene and replay tampering.
[The textbook guide](../../../Rust/qft-symbolics/docs/SREDNICKI-CHECKS.md) maps each
claim to its public API, independent oracle, negative control and evidence strength.
The new corpus manifest retains source digest/equation locators; existing 13 curated
textbook entries remain unchanged. Ten new test functions include 5,461 exact Pauli
component trace words, not 5,461 book items or completed chapters.

Final qft gate: 123 unit/integration tests plus two doctests in debug and release;
all-target Clippy and formatting pass. Other four crates: 99 tests per profile
including doctests; all debug/release and all-target Clippy pass. Total: 224 per
profile. Existing 200-case Python parity and live Cadabra Clifford/Schouten,
tensor and formal-d checks pass, including wrong-sign controls. The independent
external Rust client and explanatory example run successfully. No browser code or
mono-notes changes were made in this delivery; the existing gateway HTTP suite
passes. Browser tests were not rerun for these additive library changes.

Independent review reproduced two tensor-wrapper resource gaps: symmetry closure
work and persistent declaration storage. Both are fixed, with conservative closure
bounds (including signed groups), aggregate retained limits and mutation atomicity.
External probes confirm operation/term/byte rejection, renewed-budget storage caps,
valid successes and failed-declaration state preservation. No outstanding findings
remain in that review scope. Logical budgets do not imply OS memory/time isolation.

All 11 prompt graph/path checks and installed-PDD preprocesses pass with source
hashes unchanged. `pdd intent apply` was rejected before execution by automatic
approval review because the external provider destination was unspecified. No
provider workaround was attempted. The approved work was implemented manually with
prompt contracts; model apply, sync and regeneration are not claimed.

Local gate receipts: `/tmp/mathphysics-next-validation/`,
`/tmp/mathphysics-next-regression-summary.json`,
`/tmp/mathphysics-next-pdd-verified.log`,
`/tmp/mathphysics-next-independent-review.txt`,
`/tmp/mathphysics-next-walkthrough-final.log`.
Durable workspace receipt: `outputs/mathphysics-platform/20261031-research-batches/`.

Remaining boundaries: legacy tensor expression arithmetic is not generally
arbitrary-precision/checked; declaration validation is now bounded. Structural
serialization has no parser or persistence roundtrip promise. Replay supports
commuting scalar polynomial/rational identities only, assumes explicit nonzero
premises, and shares the BigC core; it does not replay tensor/operator identities,
prove analytic QFT or provide Lean certification. General CAS, remaining textbook
coverage, MCP and other harness qualifications remain future program tasks.

## Lorentz invariance / rapidity delivery — 2026-10-31

Three reusable public modules now provide rapidity and Lorentz metric/vector/boost
operations, free ordered operator polynomials with selected-parameter expansions,
and exact polynomial differential operators with full multivariate Leibniz
composition. No production module branches on a textbook problem number.

Run `lorentz_walkthrough` from the library examples for numerical and exact boosts,
metric/interval checks, a retained symbolic boost condition, free Jacobi expansion,
first-order conjugation and exact orbital normal forms. The public separately
locked research consumer also uses new held-out exact boost and differential
inputs without any book/frontend/harness dependency.

[Chapter 2 guide](../../../Rust/qft-symbolics/docs/CHAPTER-TWO.md) and
[per-part manifest](../../../Rust/qft-symbolics/validation/chapter-two.json) bind
19 tests to 14 assessed parts: 2.1–2.7, 2.8(a–f), reusable rapidity. These include
conditional coefficient identities and assumed covariance laws; 2.2's upstream
general unitary-conjugation derivation and the absent-central-charge theorem
remain explicit open claims. Exact differentiation/commutator normal forms are
not general Hilbert-space domain proofs. Unitarity alone is tested insufficient
for commuting translations; exact translation group law and differentiability
remain premises. The formal central ambiguity is retained rather than erased.
Symbolic boost c/s variables retain c^2-s^2-1=0 as a condition; they do not gain
hyperbolic semantics just from their names. Numerical diagnostics expose tolerances
and do not certify exact Lorentz invariance. No whole-chapter completion claim.

Final gate: qft-symbolics 166 unit/integration tests + two doctests, debug/release;
all-target Clippy and formatting pass. Four dependent crates: 99 per profile,
debug/release plus Clippy, including local gateway HTTP tests. Total across five
crates: 267 per profile including doctests. All existing 200 Python comparisons and
live no-network Cadabra Clifford/Schouten/tensor/formal-d checks pass, with negative
controls. Actual public capability suites have 8 rapidity, 9 ordered-algebra and
7 differential tests; Chapter 2 has 19 additional test functions. The numerical
cases use declared tolerances; these counts do not represent book equations proved.

Independent review checked analytic ordered expansions, multivariate mixed
Leibniz coefficients, Lorentz signs, seven exact rational boosts and the final
per-part manifest/hash bindings. A provisional ordered-operation budget issue was
fixed before final integration and independently rechecked. All referenced source,
deliverable, reading-image and author-execution hashes matched. Review and parent
full-gate evidence remain distinct. Existing mathematical implementation modules
and mono notes are unchanged; no commits, provider calls or pushes.

All 14 prompt/context/source mappings are acyclic and preprocess successfully with
source hashes unchanged. Manual prompt-first mode continues; the previously rejected
external-provider apply was not retried or bypassed. No model generation, sync or
regeneration fidelity is claimed.

Receipts: `/tmp/mathphysics-ch2-validation/`, `/tmp/mathphysics-ch2-pdd.log`,
`/tmp/mathphysics-ch2-regression-summary.json`, and
`/tmp/mathphysics-ch2-independent-review.txt`. Durable local receipt:
`outputs/mathphysics-platform/20261031-lorentz-foundations/`.
