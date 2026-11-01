# Assigning work across AI agents

The plan is at `README.md`; concrete dependencies and acceptance criteria are in `tasks.json`. The file contains proposed tasks, not an executing scheduler. Use Hermes Agent, Pi, Claude Code, Codex or Grok Build according to the capability actually installed. No background agent, paid call, worktree or CI job has been launched by this planning package.

## First assignments

Start these planning/verification assignments independently after the user accepts the program's meaning. Source inspection already authorized by the planning request need not wait for product implementation approval.

| Assignment | Owner | Work that can proceed | Output |
|---|---|---|---|
| F01 | Integrator/schema specialist | Define namespaces, source anchors, assumptions, evidence and state transitions | Reviewable schemas and tiny examples |
| F02 | Verification agent | Rerun existing Rust suites and available external oracles | Fresh baseline, versions/hashes and negative-control outcomes |
| F17 | CAS feature auditor, after F01 | Inventory useful feature families across reference systems | Capability/parity registry, original examples and evidence requirements |
| F03 | Source specialist, after F01 | Review chapter/page map and OCR provenance, including front/back matter | Corrected source map and quarantine disposition |
| F04 | Physics reviewer, after F01/F02 | Freeze the shared convention registry | Explicit mathematical conventions and incompatible-mode tests |
| CH-001 through CH-097 | Chapter specialists, after F01/F03/F04 | Inventory complete assigned chapters | Item ledgers and proposed implementation batches |

For a small first chapter-assignment set, choose CH-003, CH-035, CH-047, CH-050, CH-060 and CH-069. These exercise scalar foundations, existing Rust coverage and the nonabelian extension. All remaining chapters still have explicit tasks; these six are scheduling examples, not a reduction of scope.

Next, assign F05 for scoped Rust PDD adoption, F08 for the expression design, F12 for the tool contract, F11 for formal-checker evaluation and G01 for bounded GPU feasibility. Respect their dependencies in `tasks.json`. Keep implementation of shared coefficient, expression and tensor types serialized by ownership even when research runs in parallel.

Do not start every capability task at once. Two to four authors plus one independent reviewer is an initial scheduling proposal, adjusted to measured CPU/RAM/VRAM and available provider limits. The RTX 3060's memory must be scheduled with any local LLM use. Only explicitly authorized providers/harnesses may execute; record actual usage rather than guessing it.

## Mandatory library reuse gate

Before dispatching any chapter/paper implementation, read REUSABLE-LIBRARY.md and the reusable library story. Each assignment must name existing APIs it reuses, the actual missing capability (if any), its single owning module, generic input/domain/error contract, source-specific fixture scope and independent acceptance. A chapter assignment does not authorize a chapter-specific production solver. Inspect/characterize existing clients and convention adapters before extracting or renaming code. Reference code is inspiration, not an instruction to port C++ structures literally.

Required return fields: `existing_api_reused`, `capability_gap`, `owning_module`, `thin_fixture_paths`, `independent_oracle`, `external_consumer_check`, `compatibility_or_migration`, `remaining_limits`. A source inventory alone cannot satisfy these implementation gates.

## Assignment template

Copy the following into the chosen harness with a concrete task ID and an exclusive file scope:

```text
Work in the mathphysics repository. Read AGENTS.md and PDD.md, then
  docs/plans/rust-physics-platform/README.md
  docs/plans/rust-physics-platform/tasks.json
  docs/plans/rust-physics-platform/source-inventory.json
and the existing implementation/verification documents relevant to TASK_ID.

Assignment: TASK_ID.
Role: its owner_role.
Exclusive write scope: ASSIGNED_PATHS.
Dependency evidence: DEPENDENCY_ARTIFACTS_AND_HASHES.
Mode: inventory/design, implementation, or independent review as explicitly assigned.

Reinspect Git state and PDD ownership. Preserve concurrent work. Check that task
dependencies have actual reviewed evidence, not merely a checked status box.
Do not widen the existing Python rotation prompt's ownership. For a new Rust
production unit, prepare its bounded intent plan and review card; apply only
meaning already approved by the user. Do not introduce provider/auth changes.

Use the original PDF for source ambiguities. Treat OCR and matching CAS output
as evidence requiring assessment. Never edit source manuscripts to match code.
Keep equations, exercises, diagrams and prose claims in distinct namespaces.
Prefer Rust; justify foreign-language code and mark external runtime calls.
Never turn overflow, unsupported input, timeout or backend failure into PASS.
Do not alter independent acceptance tests to accommodate a solver's answer.

Return: changed paths; assumptions and mathematical domains; reusable APIs;
source/item IDs; derivation; tests and negative controls; reproducible commands;
evidence hashes/locations; unresolved questions; and next dependency-ready tasks.
Do not commit, push, merge, post externally or spawn additional agents unless
your assignment explicitly authorizes those actions.
```

## Chapter inventory item

Use a stable namespaced ID such as `srednicki:ch035:exercise:35.1` or `srednicki:ch035:equation:35.1`. These are different items even when the printed numbers match. Give unnumbered displays, prose claims and diagram obligations stable IDs tied to source page/block anchors. Do not rely solely on a mutable OCR line number or automatically generated expression hash.

Each item needs:

- Primary PDF locator with explicit page convention, source digest, book chapter and equation/exercise/figure label when present; secondary solution/DGDT references kept separate.
- Reviewed statement and assumption/domain/convention references; ambiguities and unresolved source readings.
- Mathematical prerequisites and computational dependencies, including items in other chapters.
- Reusable Rust operations, exact versus approximate behavior, required calculus/distribution/operator/regularization semantics.
- Implementation state and evidence state as separate fields; native versus external backend and conditional/research-gap flags.
- Positive, boundary and wrong-sign/type/domain controls; independent checker/oracle and its independence limits.
- Proposed implementation batch with owned paths, acceptance gate, resource budget and scope suitable for one reviewed PDD unit.

Equation candidates and the 13 existing exercise checks are initial cross-references, not a complete chapter inventory. A chapter task reads every section and every exercise subpart, including text between displays. Do not generate one opaque monolithic “implement chapter” prompt.

## Integration and independent review

A reviewer should reproduce the derivation or choose an independent route, try to falsify the result, verify source readings, inspect side conditions and replay certificates with a separately implemented checker where feasible. Reading another agent's explanation and agreeing is insufficient.

The integrator reconciles shared APIs and conventions, runs cross-component tests and a held-out regression slice, then updates manifests only after evidence exists. Unknown/unsupported/research-gap outcomes remain visible. Runtime logs/build artifacts stay outside Git; retain portable scripts/manifests and approved sources. Never mark all 97 chapters implemented just because all 97 inventory tasks are complete.

Use local claims with an owner, task ID, base commit/dirty hashes and explicit paths to prevent two active agents taking the same write scope. A claim/lease mechanism is a planned engineering detail, not something this package pretends to implement. Do not use a shared mutable worktree for simultaneous source edits.

## PDD and remaining initial rotation work

The previous Python rotation setup remains separate and uncommitted on `feat/pdd-rotations-adoption`. Its evidence record lists unfinished baseline/regeneration/contract work. This program does not mark those steps complete, move their files, or reuse their intent ID. Plan each new Rust unit against its actual directory and the applicable installed CLI; root discovery of a nested configuration does not establish product ownership.

`proposed_write_paths` provides concrete candidate ownership for each task. Existing files retain their current owners until an assignment explicitly transfers scope. Proposed new crate names are reviewable boundaries, not a requirement to create every directory now. Foundational changes can affect downstream modules; enumerate those consumers and integration tests before dispatch. Budget values are set for the actual operation and hardware, using each task's `required_budget_fields`; a task cannot treat an omitted budget as unlimited.
