# Rust platform first wave verification

User authorization: approved program execution and subagents on 2026-10-31. Existing branch: `feat/pdd-rotations-adoption`; earlier rotation setup, planning work and unrelated generated build report are preserved. No commits or pushes.

Runtime evidence is retained outside the source repository:
`/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-wave1/`.
Subdirectories contain commands, logs, source/artifact hashes, source-review images and independent review reproductions. Temporary build targets remain outside Git; source manifests and curated data are versioned candidates in the repository.

## Existing mathematical kernels

The unmodified `Rust/qft-symbolics/validation/run.sh` passed offline debug/release, formatting, Clippy with denied warnings, corpus mappings, exact positive/negative controls, live no-network Cadabra transformations and 200 NumPy parity fixtures. Fresh count: 55 integration tests and two compile-fail doctests per mode. This is the kernel's stated bounded algebra, not general theorem proving.

`Rust/lie-groups` debug/release passed three unit tests, 48 integration tests and four doctests per mode. Repeated modes are not additional unique tests. Baseline has source hashes and Git snapshot, not just an old receipt.

## Metadata contracts

`Rust/physics-contracts` debug/release: 19 tests passed, including two real CLI tests. Clippy passed with warnings denied. Runtime is Rust with pinned serde/serde_json, using the local Cargo cache; no CAS, Python or provider is required. Unknown/duplicate fields, missing sources, bad hashes, invalid locators/IDs, unknown references, statement/domain/assumption mismatches, broken backlinks, duplicate artifact paths, numerical overclaims, fake formal metadata and budget limits have negative tests.

Independent review identified four defects and re-reviewed fixes. Thirty JSON reproductions and six direct Rust cases met expected outcomes. No proof-assurance bypass was identified in that reviewed scope. This is not an exhaustive fuzzing audit.

The four schemas pass Draft 2020-12 metaschema checks and 14 test groups, including 432 evidence-kind/arithmetic/outcome/budget combinations. An existing Python JSON Schema validator was used only as an independent development oracle; reproduction scripts are retained in the evidence directory. Schema validity is weaker than the Rust cross-record semantic checks, and both are weaker than mathematical verification.

## Sources and inventory

97/97 chapter starts checked against PDF text/printed folios, with a separately observed +24 offset at each. One heading image (chapter12) resolved a broken text-layer hbar; the other heading starts are text checks. All ten source hashes match the planning inputs. All 72 unscoped Nougat-only tags remain quarantined. Historical 766 conflict resolutions reconcile as 435 automatic and 331 recorded vision decisions; no fresh mathematical certification is inferred from those records.

Chapter3's author reports reading/visually inspecting all nine pages and inventories 102 items: 38 numbered equations, 12 intermediate equations, 35 prose/definition/example items and 17 problem/subpart/unlettered exercise-goal records. Equation3.35 was absent from the original OCR candidate ledger and is explicitly inventoried. Cutoff approximation, delta support, dimensional restoration and a translation-index issue are recorded for review. All items remain unassessed, with no evidence claims. The Rust bundle validator returns record_valid only for metadata; independent source review is tracked separately.

## PDD boundary

Read-only intent plan for the new contract scope: `rust-mathematical-item-and-evidence-contracts-9a1b05f4`. User approval of the parent plan carries to this bounded unit without a new mechanical checkpoint. Initial implementation is agent-authored under versioned prompts. Rust preprocessing succeeds; path dry runs resolve correct code/tests and report UNBASELINED. No full `intent apply`, generated story contract, fresh regeneration or baseline certification is asserted.

## Delegation

The actual session subagent launcher was used, inheriting the parent model. Assignments used isolated `/tmp` output directories with exclusive writes; only the integrator modified the repository. Initial agents: Descartes (fresh baseline and independent review), Bernoulli (source map and chapter inventories), Locke (capability audit and conventions), Newton (schema implementation and tests), Halley (bounded Rust tool bridge). Completed workers are closed when no longer needed. Task status must not imply a running background agent after this turn ends.

## Integrated convention gate and callable tool

The six new qft-symbolics convention tests pass in debug and release after integration, with all-target Clippy warnings denied. They cross-check Lorentz/Clifford sign, epsilon/gamma5 orientation, spinor bracket identities, scaled polarizations/reference-domain rejection, coefficient conjugation versus the physical real slice, and bare-amplitude/SMGA interpretation. Existing source implementations are unchanged; this adds characterization tests, not new physics functionality.

physics-tools passes 11 tests in debug and release: two unit, five library integration and four actual CLI tests. Formatting and Clippy pass. The integration uses a portable `../qft-symbolics` dependency. Independent qualification performed 162 actual CLI invocations: fresh-label formulas, directly constructed 4x4 Dirac matrices, repeated/odd words, strict JSON/UTF-8/field handling, byte/pairing/label limits and a real 8192-term kernel failure. That failure returns inconclusive, null evidence_kind, formal_proof=false and exit5. No false-computed or false-proof case was reproduced in this bounded review.

The tool returns computed exact-algebra expressions within the existing checked-i128 kernel bounds. It is not the general expression engine, a proof checker, an MCP server or a complete F13 release. Only this Codex-subagent/exec route was exercised; no other harness is silently marked qualified. Corrected PDD dry runs resolve both code/test pairs with zero path failures, still UNBASELINED.

Chapter35 has61 metadata-valid, unassessed items and a separate partial exercise crosswalk. Chapter3 review refinements preserve102 IDs/locators while correcting assumptions and retaining unresolved source claims. Initial source inventory is separate from executable chapter implementation.
