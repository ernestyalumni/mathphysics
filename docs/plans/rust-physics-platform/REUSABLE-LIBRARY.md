# Capability-first reusable Rust library

This is the engineering interpretation of the new [reusable research library story](../../../user_stories/story__reusable-rust-research-library.md). It strengthens the existing platform scope: Srednicki drives requirements/validation; a chapter outline is not a production module tree. Product meaning lives in PRODUCT_INTENT and the story; implementation details remain in bounded module prompts when adopted.

## Boundaries and dependency direction

| Layer | Ownership | Permitted source-specific content |
|---|---|---|
| Reusable mathematical kernel | Existing qft-symbolics / lie-groups and justified extensions | Definitions/algorithms and explicit convention profiles; no book filesystem or item-number dispatch |
| Domain capabilities | Typed expressions, indices/operators, kinematics, diagrams/amplitudes as their gates are implemented | Mathematical theory/context parameters; no one-off chapter solvers |
| Evidence and provenance | physics-contracts metadata now; independent replay/checkers later | Book/paper identifiers, version/hash, claim/assumption/evidence bindings |
| Adapters | physics-tools CLI, physics-gateway, future MCP; UI callers | Decode/encode, authorize/limit, map explicit conventions, call public library operations |
| Validation corpus | Source inventories, fixtures, examples and independent tests | Chapters/equations/exercises/papers, expected answers and independently justified checks |

Dependency direction is consumer → adapter/domain capability → library primitive; source data references the APIs it needs. The library must not depend on a chapter fixture, UI, server, particular AI harness or external CAS. Independent oracle implementations remain separate on purpose. Existing internal module boundaries are the starting point; these layers do not authorize creating a crate for every noun or a single enormous generic engine.

## Per-item implementation procedure

1. Record the source statement and explicit domain/conventions; find its reusable operations in existing libraries and the capability registry.
2. If supported, add only a thin source adapter/fixture and independent expected evidence. Do not generate a function for its subsection.
3. If unsupported, record a capability gap with a concrete second input/consumer or clear mathematical generality; design the smallest extension at the existing owning module before proposing a new crate. A capability need not have arbitrary generality to be reusable.
4. Freeze its typed public interface, failure/resource behavior, convention conversion and compatibility obligations. Run the bounded PDD plan; update its prompt before implementation. Preserve characterization tests and unrelated sources.
5. Validate using fresh labels/parameters and independent representations/negative cases, then exercise the API from a separate Rust consumer and a real adapter. Add workload benchmarks only where meaningful.
6. Link the source item to that capability/version and evidence. An inventory row or valid metadata is not chapter completion or mathematical truth.

Chapter batches may group source review and acceptance work; their CH/section IDs are scheduling/provenance labels, never instructions to create chapter-specific production solvers.

## Current audit and next extraction candidates

- `Rust/qft-symbolics/src/lib.rs` exports capability modules: scalar, polynomial, clifford, tensor, matrices, lorentz, twistor, spinor and smga. Preserve this useful starting structure. Current bounded exact arithmetic and specialized domains are not yet general production research infrastructure.
- `clifford::trace(labels)` and `normal_order(labels)` take mathematical inputs. `sandwich_srednicki(labels, dimension)` documents a convention and accepts general words/dimension; its name alone is not a subsection-specific defect. Future convention abstraction must preserve its behavior/adapter compatibility.
- `tests/textbook_extended.rs` imports shared lorentz/matrices/polynomial/scalar/twistor operations. Its `ch50_...` tests preserve traceability. Independently constructed matrix/determinant expectations should not be replaced by calls to the same production operation.
- `physics-tools` already wraps the qft-symbolics kernel rather than duplicating gamma traces. This is a good library/client direction; it is still a bounded seed interface.
- `lie-groups::su2_observation` now owns pure typed observation mathematics; `physics-gateway/src/su2.rs` parses and serializes its result. The plus-i Pauli lab and minus-i quaternion convention have an explicit tested conversion. Preserve that boundary.
- `qft-symbolics::tensor_cases` is an exported verification-case surface. Audit fixture/verification exports separately from the production API and plan compatibility before moving them into examples/tests. Do not mass-rename or delete working consumers just to achieve a tidy directory tree.

## Research-grade gates still open

The first F06/F07 checked coefficient/matrix/spinor/rational slices, F08/F09 structural expression/domain model and F10 independent scalar replay now exist. Remaining legacy tensor/transport migration, wider domain semantics and independently replayable QFT transformations remain prerequisites for broader research use. Existing i128/panic boundaries and specialized assumptions stay explicit until replaced. General color/Feynman/amplitude engines and formal checker integration are not delivered by this design change.

A paper extension records a versioned arXiv identifier (including legacy forms where applicable) or content digest, source locator, precise claim, assumptions/conventions and evidence classification. The runtime should execute the mathematical operations without knowing where that input originated. No automatic PDF reading-to-proof pipeline or exhaustive paper support is claimed.

## Idiomatic adaptation of C++ practice

Borrow small cohesive components, explicit public contracts, reusable numerical kernels, separate clients/tests, reference comparisons and measured performance. In Rust use normal modules/crates, ownership/borrowing, composition, suitable traits/generics and explicit errors. Dynamic dispatch, unsafe/FFI, CUDA specialization and additional dependencies require actual workloads and bounded ownership; C++ header splitting or inheritance hierarchies are not a template to copy. Prefer recent maintained reference code, but audit tests and caveats rather than treating commit dates as correctness evidence.

Concrete dated repository examples follow after the read-only audit. No source from those repositories is changed or automatically copied.

## Dated source examples actually inspected

These are architectural inferences from read-only source review, not reproduced test results or copied implementations. The referenced repositories stayed unchanged. Dates below describe Git history, not filesystem mtimes. Prefer the October2026 Rust examples over older C++ patterns when they address the same boundary.

| Reference | Git evidence | Useful analogue for mathphysics |
|---|---|---|
| [Cosmos Rust workspace](../../../../Propulsion/Cosmos/Rust/Cargo.toml), [dynamics dependencies](../../../../Propulsion/Cosmos/Rust/cosmos_dynamics/Cargo.toml) | 2026-10-25, `bd66277c` / `96500056` | Separate coherent capability/domain consumers with explicit dependency direction; not a library per book section. |
| [Embedded ODE interfaces](../../../../Propulsion/Cosmos/Rust/cosmos_numerical/src/ode/runge_kutta/embedded.rs), lines55–128 | 2026-10-25, `bd66277c` | Shared driver, traits/associated types and interchangeable policies. Borrow composition boundaries, not ODE-specific abstraction names for algebra. |
| [Order conditions](../../../../Propulsion/Cosmos/Rust/cosmos_numerical/src/ode/runge_kutta/order_conditions.rs), [reference tests](../../../../Propulsion/Cosmos/Rust/cosmos_numerical/src/ode/runge_kutta/hairer_reference_tests.rs), lines213–281 | 2026-10-25, `96500056` | Mathematical residuals plus independent reference fixtures, explicit outcomes/budgets/counters and tolerance cases. Reference provenance was read, not regenerated in this audit. |
| [Quaternion convention library](../../../../Propulsion/Cosmos/QuaternionConventionLab/rust/src/lib.rs), lines1–39 and143–246 | 2026-08-17, `999ca8c0` | Conventions on reusable operations/types; meaningful quarter-turn/antipode/layout tests. SO(3) antipodal equivalence must not erase distinct SU(2) elements. |
| [UnionFind interface](../../../../HrdwCCppCUDA/Voltron/Source/DataStructures/UnionFind.h), [implementation](../../../../HrdwCCppCUDA/Voltron/Source/DataStructures/UnionFind.cpp), [library target](../../../../HrdwCCppCUDA/Voltron/Source/DataStructures/CMakeLists.txt) | 2026-06-01, `05e8725c` | A reusable capability owns its algorithm; clients supply inputs. Rust should protect invariants rather than copy the exposed mutable fields/unchecked indices. |
| [Capability tests](../../../../HrdwCCppCUDA/Voltron/Source/UnitTests/QuestionsDEntrevue/Streams/FireAlertParser_tests.cpp), [test-target organization](../../../../HrdwCCppCUDA/Voltron/Source/UnitTests/CMakeLists.txt) | 2026-06-01, `9b467b83` | Tests align with capabilities and import public interfaces; source identifiers need not determine production APIs. Existing TODO/behavior tests are not correctness certification. |
| [Error-handling contract](../../../../HrdwCCppCUDA/Constructicon/Source/Utilities/ErrorHandling/HandleReturnValue.h), lines25–37 | 2023-12-03, `26c16032`; implementation/tests last changed2023-11-30, `ce74fd6f` | Separate underlying operations from error interpretation; adapt to Rust Result/domain errors, not virtual functors, global errno or mutable last-error state. This example is older, not presented as recent source. |

The HrdwCCppCUDA dependency-maintenance commit `ea0880c4` has author date2026-05-31 and committer date2026-06-01; dependency maintenance does not make all underlying algorithms recent. Propulsion HEAD was `95249d8` (2026-10-27 notes), while the relevant reusable Rust implementation changes were October25.

Cautions: an empty success-only construction test is not behavioral evidence; platform-specific bitwise reference assertions need explicit platform/tolerance contracts; conventions and partial functions must become tested, documented public boundaries. Do not import broad dependency/test coupling or blindly migrate old indexing/ownership patterns. Full audit notes with exact source lines and hashes are in [Propulsion evidence](reuse-intent/propulsion-reference-audit.txt) and [HrdwCCppCUDA evidence](reuse-intent/hrdw-reference-audit.txt).

## Concrete acceptance example

One shared `clifford::trace(labels)` operation should serve textbook gamma-trace checks, fresh non-book words and the existing JSON bridge. A source fixture supplies labels, a declared Clifford form/convention and an independently derived expected identity; it does not implement a second trace algorithm or branch on `47.2`. Similarly, a group library owns group operations, while an Explore adapter translates the documented coordinate convention and transports inputs/results. This example points to existing APIs; the wider expression/replay/research gates remain queued.

## Verification/status of this refinement

The exact new request and successful read-only PDD plan are in reuse-intent/. This turn changes the human-authoritative story, product/architecture/dispatch guidance and task metadata only. No Rust implementation was refactored, no external repository was modified, and no production/arXiv-proof readiness was claimed. No runtime tests were necessary for these documentation-only changes; JSON/reference/dependency checks below validate the updated plan, not mathematical claims. Existing fresh Rust baseline evidence remains in the previous status receipt.
