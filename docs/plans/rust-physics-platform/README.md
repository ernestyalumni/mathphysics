# Rust mathematical physics platform and Srednicki implementation plan

Status: approved for execution on 2026-10-31. The user explicitly authorized starting the plan and spawning subagents. Current implementation and evidence are tracked in [EXECUTION.md](EXECUTION.md) and tasks.json; the design below remains the approved program scope.

Build an extensible Rust mathematical physics environment that Ernest can use through his choice of AI agent, with independently checked calculations, explicit assumptions, traceable derivations and book-wide coverage. Use Srednicki's entire Quantum Field Theory as the first comprehensive capability driver. Expand toward the useful functionality of Cadabra, SymPy, SageMath, Maple, Mathematica and other systems through individually specified capabilities. The ambition includes capabilities outside Srednicki; those remain in an extensible feature backlog rather than being silently excluded.

Start with the existing `Rust/qft-symbolics` and `Rust/lie-groups` implementations. Develop capabilities from mathematical definitions and compare them with independent references. Replacing every feature of several mature CAS products is an open-ended program; delivering a tested, book-complete capability set is a concrete milestone within that program.

## Original platform stories and current refinements

### Story A Agent accessible mathematical physics tools

As a mathematical physics researcher using mathphysics, I want my chosen AI agents to call a shared, primarily Rust collection of symbolic, numerical and verification tools, so that they can derive results, test propositions, find counterexamples and produce auditable evidence with capabilities comparable to and eventually extending the tools available around Lean, Cadabra, SymPy, SageMath, Maple and Mathematica.

Acceptance:

1. A Rust library and stable versioned tool interface expose discovered capabilities, typed requests, assumptions, resource limits, outputs and evidence. A user need not write Python or C++ to use supported calculations.
2. A Rust CLI is usable from Hermes Agent, Pi, Claude Code, Codex and Grok Build through their supported tool mechanism. An MCP server is supplied for compatible harnesses. Actual integration is tested per installed harness; unsupported integrations are reported accurately.
3. Results distinguish derived expressions, exact identities in a declared algebra, numerically supported results, independently replayed certificates, formal proofs, counterexamples and inconclusive attempts. An LLM judgment or matching CAS output cannot become a formal proof.
4. Rust is the first implementation choice; CUDA C++ or other C++ is a documented bounded fallback, followed by Python only where needed. Existing external CAS installations can be optional reference backends. The core Rust-only test lane must not import them or quietly delegate to them.
5. Any frontend is written in TypeScript. Source references, derivation steps, assumptions and verification status remain visible to the user.
6. Parallel agents can work on independent tasks without overwriting shared mathematics, counting peer agreement as proof, or exhausting shared GPU/provider resources.
7. GPU implementations are adopted where measured workloads benefit, with a reproducible CPU reference and an explicit exactness or numerical-error contract.
8. The feature registry includes broadly useful algebra, analysis, geometry, proof and numerical capabilities beyond the first textbook. “Better” is demonstrated by held-out correctness, reproducibility, usable agent interfaces and measured performance on named tasks.

### Story B Executable Srednicki with traceable derivations

As a researcher studying Srednicki's Quantum Field Theory, I want every chapter's equations, intermediate derivations, worked examples, exercises and substantive physics claims mapped to reusable tools and reproducible checks, so that my chosen AI agent can explain, reproduce and investigate the book's calculations using mathphysics while showing the assumptions and the strength of the evidence.

Acceptance:

1. Every chapter 1–97 has a source-reviewed inventory. Include unnumbered formulas, prose claims, diagrams and exercises; an equation-tag index alone is insufficient. Front matter, references and index are reviewed for additional computational material and cross-references.
2. Each item has a stable identifier, source location, source-review status, statement, explicit conventions/assumptions, prerequisites, reusable operations, implementation location, tests, evidence and remaining gaps.
3. Preserve Ernest's solution sources and notation. Compare his derivations with the book; differences become reviewable records rather than automatic edits to his manuscripts. DGDT material supplies relevant geometry prerequisites and a separate extension backlog.
4. A derivation can be replayed from inputs and justified transformations; storing the book's answer as an output is not implementation. Boundary cases and intentionally incorrect variants expose sign, index, normalization, domain and regularization mistakes.
5. Distinguish ordinary functions, distributions, operator identities, perturbative expansions, regulated expressions, approximations and physical models. Each requires the appropriate mathematical contract.
6. Never claim a physical assertion or open problem has been proved because a finite simulation, formal manipulation or special example agrees. Retain difficult items as explicit research/formalization gaps. An accounted-for item is not automatically an implemented item.
7. Book completion requires per-item review and a fresh independent chapter audit. Report coverage separately for source inventory, implementations, exact checks, numerical checks, formal proofs and unresolved items.

## Current evidence and source boundaries

Inspected 2026-10-31:

- `Rust/qft-symbolics/corpus.json` records 13 selected textbook entries. Its historical verification receipt records 54 integration tests and two compile-fail doctests; these counts are documented history, not tests rerun for this planning task. The current engine is bounded: checked i128 coefficients, a limited commuting tensor layer and selected Clifford/spinor/amplitude operations.
- `Rust/lie-groups` supplies numerical and exact SO(3)/SU(2) and related algebra. Neither crate is currently PDD-managed merely because Python/rotations is managed.
- The existing reconciler has 3,011 equation records: 1,366 agree, 766 conflict, 651 marker-only and 228 nougat-only. All 766 conflict records have a resolution field in `equations_resolved.json`. That does not establish a complete visual source audit. Of all records, 2,939 have syntactically plausible chapter/equation tags; 72 require tag triage. The ledger preserves all records, including suspicious ones.
- The full parsed contents lists all 97 chapters and their printed starting pages and stated prerequisites. `chapters.csv` transcribes those entries and links existing partial checks; its titles and prerequisites still need chapter-owner source review. No constant PDF-page offset is assumed.
- The original textbook PDF is authoritative for textbook transcription. Ernest's `LaTeX_and_pdfs/Srednicki_QFT.tex` is a distinct solution/reference source; `DGDT_dump.tex` and its PDF are separate mathematical context. Existing OCR derivatives stay untouched.
- Docker inspection confirms `cadabra2-ubuntu:24.04`, image `sha256:48e6f90b4716e1e107deebe52da8b83e572b6ef6df36a1eabcf8fc72ba78ef50`. Availability does not mean its transformations were executed in this planning session.
- Sage checkout exists at the supplied TQFT path, version file `11.0.beta5`, commit `36e49d91b05`. This establishes source availability, not a working Sage runtime.
- GPU inspection reports GTX 980 Ti 6,144 MiB and RTX 3060 12,288 MiB; CUDA compiler 13.0.48 and driver 580.76.05. Use the RTX 3060 for the initial CUDA 13 lane. CUDA 13 removed offline compilation/library support for Maxwell; using the older card needs a separately tested compatible toolchain. No GPU benchmark was run.
- Lean clone location is pending Ernest's follow-up; planning and Rust algebra work do not depend on it. Lean's documentation describes its trusted kernel and links an independent Rust checker. Evaluate that option before proposing a new general proof foundation.

This planning pass enumerates the full table of contents and existing equation records. It is not a claim to have mathematically audited every derivation or visually certified every page. The 97 chapter inventory tasks explicitly perform that exhaustive work.

## Requirements that govern all work

- The product MUST prioritize Rust, then justified C++/CUDA C++, then Python. Native dependencies and foreign-runtime calls MUST be reported separately from the language of our own code.
- Frontend code MUST be TypeScript. Rust services and WebAssembly may supply computation. Browser graphics use browser graphics APIs; CUDA is a native/server compute option, not a browser graphics API.
- Every transformation MUST preserve or explicitly update its domain and side conditions. Division, square roots, logarithms, limits, integration by parts, Wick rotation and regularization MUST NOT silently discard restrictions.
- Results MUST identify the exact convention set: metric signature, Clifford sign, epsilon orientation, index spaces, conjugation/reality rules, units, incoming/outgoing momenta, factors of i, Fourier transforms and field/vertex normalization.
- Four-dimensional gamma5 identities MUST NOT be extrapolated to symbolic dimension without an explicit scheme. Lorentzian, complexified and (2,2) spinor calculations MUST remain distinguishable.
- Resource exhaustion, overflow, unsupported expressions, parser ambiguity and backend failure MUST produce an inconclusive/unsupported result, never PASS. Exact arithmetic MUST NOT silently fall back to floating point.
- Symbolic exactness, numerical evidence and formal proof MUST remain distinct. A proof certificate MUST identify its checker and assumptions; checking a custom algebra certificate does not automatically establish Lean-equivalent foundational assurance.
- Existing mathematical sources and independent tests MUST be preserved. Features MUST have fresh evidence before being advertised. Native Rust capability and external-backend capability MUST have different labels.
- Original mathematical implementation and optional oracle comparison are the default. Inspect dependency/source licenses and provenance before reusing upstream code; changing language is not a provenance strategy. Proprietary packages are not required to run the core product.

## Architecture to develop incrementally

```text
Hermes / Pi / Claude Code / Codex / Grok Build
                   |
       Rust CLI and optional Rust MCP server
                   |
         typed operation and capability registry
          /                |                \
 exact Rust algebra   numerical Rust     derivation/certificate checker
          \                |                /
     book item registry and reproducible evidence store
                   |
 optional external reference adapters and CUDA acceleration
                   |
          optional TypeScript research interface
```

Proposed bounded units, not existing crates or implemented APIs:

| Unit | Responsibility | Primary independent gate |
|---|---|---|
| Coefficients and domains | Big integers/rationals, Gaussian rationals, polynomial/rational domains, algebraic numbers, assumptions | Large exact witnesses, overflow/cancellation/domain controls |
| Expression and binders | Typed immutable expression graph, scoped indices/variables, canonical serialization, substitutions and source maps | Alpha equivalence, capture rejection, round trips and malformed input |
| Tensor and statistics | Mono-/multi-term symmetries, contractions, differential operators, Grassmann/noncommuting products | Cadabra comparisons plus independently derived identities and negative controls |
| Scalar analysis | Polynomial ideals, series, special functions, symbolic differentiation/integration, branch-aware limits | Independent analytic identities, differentiated primitives, certified numerics |
| Fields and perturbation | Wick contractions, Gaussian/Grassmann functionals, diagrams, symmetry factors, LSZ, Feynman rules | Explicit low-order derivations and crossing/normalization checks |
| Loops and schemes | Integral families, parameterization, Laurent expansions, renormalization and RG | Known exact residues and independently derived one-loop examples |
| Groups and geometry | Lie algebra/representations, color, forms, covariant derivatives, bundles and topology | Jacobi/Bianchi/Cartan tests and source-reviewed geometry examples |
| Proof evidence | Small replay checker for declared rewrite rules and certificates; optional formal-proof adapters | Forged-step/certificate rejection, assumption audit and independent replay |
| Agent tools | Rust library/CLI, capability discovery, schemas, MCP, asynchronous jobs | Real harness calls and unsupported-input contracts |
| Evidence and book registry | Stable source identities, derivation dependencies, test manifests and completion reports | Missing-evidence rejection and independent chapter audit |

Do not split crates solely to match this table. First identify ownership boundaries in the existing crates; migrate only with characterization tests and scoped PDD adoption. Exact arithmetic, domain handling and the typed expression contract are shared dependencies, not separate reinventions by chapter agents.

## Book coverage waves

These are delivery groupings. `chapters.csv` carries each chapter's actual prerequisite edges, including prerequisites in other groups. Each chapter retains an explicit task even when its implementation shares an engine with another group.

| Chapters | Capabilities and calculations to deliver | Evidence focus |
|---|---|---|
| 1–5 | Relativistic kinematics, Lorentz tensors, free scalar fields, commutators, distributions, spin-statistics assumptions, scalar LSZ | On-shell/domain checks, explicit commutators, distributional identities |
| 6–12 | Gaussian/oscillator path integrals, generating functionals, Wick expansion, scalar rules, cross sections, phase space and units | Finite Gaussian witnesses, symmetry factors, exact low-order diagrams |
| 13–21 | Spectral representations, scalar loops, 1PI functions, counterterms, all-orders structure and effective action | Regulator/scheme declarations, poles, discontinuities, perturbative order |
| 22–29 | Noether/discrete/internal symmetries, resonances, IR behavior, schemes, RG and EFT | Conservation assumptions, decoupling/matching and scale consistency |
| 30–32 | Scalar symmetry breaking, vacuum shifts, loop corrections and Goldstone structure | Hessians, stationarity, mode counting, perturbative qualifications |
| 33–40 | Lorentz representations, Weyl/Dirac fields, spinor indices, quantization, discrete symmetries | Existing checks extended to typed index/Grassmann contracts |
| 41–50 | Fermion LSZ/propagators/functionals, Dirac/Majorana rules, spin sums, Clifford traces and helicity | Massive/massless limits, signs, completeness, independent amplitudes |
| 51–53 | Yukawa loops, beta functions and functional determinants | Matched schemes, fermion-loop signs and explicit expansions |
| 54–61 | Maxwell/Coulomb gauge, photon LSZ/functionals, spinor/scalar QED and helicity | Gauge dependence, transversality, tree amplitudes and Ward identities |
| 62–68 | QED loops, vertices, magnetic moment, beta functions and Ward identities | Scheme-aware one-loop coefficients and physical limits |
| 69–74 | Yang–Mills, representations/color, ghosts, gauge-fixed rules, beta function and BRST | Jacobi, Slavnov–Taylor/BRST domains, ghost signs, nilpotency |
| 75–81 | Chiral/global anomalies, background/Gervais–Neveu gauges, matrix fields and QCD scattering | Regulator-dependent identities, anomaly coefficients, independent gluon amplitudes |
| 82–86 | Wilson loops/lattice, confinement discussion, chiral breaking and Higgs mechanisms | Finite lattice/controlled approximations explicitly separated from continuum claims |
| 87–91 | Standard Model gauge/Higgs/lepton/quark sectors, hadrons and neutrinos | Representation/charge/anomaly checks, mass mixing and phenomenological assumptions |
| 92–94 | Solitons, monopoles, instantons, topology, theta vacua and quark effects | Boundary/topological charge checks and controlled analytic/numerical solutions |
| 95–97 | Supersymmetry, MSSM and grand unification | Graded algebra, component expansions, representation decompositions and matching |

The spectrum includes nonperturbative and conceptual arguments that are not finite symbolic equalities. Those still receive implementation/research tasks; the evidence type and any open mathematical premises remain explicit.

## Execution sequence and first useful release

1. **Inventory and conventions:** adopt a versioned item/evidence schema, review sources, rerun the current Rust baseline, freeze shared conventions and enumerate every chapter item. Multiple chapter inventories can proceed concurrently; the shared schema has one owner.
2. **Trust and foundation:** structured arithmetic errors, arbitrary precision, expression/binder/domain semantics, trace replay and a Rust tool interface. Preserve the existing textbook checks during each migration.
3. **First end-to-end release:** a chosen harness calls a real Rust operation for a newly supplied expression, obtains a traceable exact result and independent replay, then submits a wrong-sign or invalid-domain variant and sees rejection. Cover existing Schouten/gamma examples and one new scalar/field calculation. Fixed demo selection alone is insufficient.
4. **Parallel capability development:** tensor/statistics, scalar analysis, fields/diagrams, groups/geometry and formal-checking investigations. Integration follows dependency contracts and chapter needs.
5. **Rolling book releases:** release source-reviewed chapter items continuously across the entire chapter map. Avoid waiting until a universal CAS is finished. Complete foundations in a prerequisite-aware order while selecting independent later-chapter slices to exercise new engines.
6. **Agent and interface qualification:** qualify each installed harness against the same tool/evidence contract; add TypeScript browsing and derivation inspection once the underlying contract is stable.
7. **GPU optimization and broad CAS growth:** profile representative workloads, accelerate proven bottlenecks, and expand the capability registry beyond the book. A limited parallel GPU feasibility task can start earlier without changing core semantics.
8. **Independent completion audit:** certify every chapter's inventory and coverage, rerun held-out tasks, enumerate all unresolved items and assess remaining CAS feature gaps.

This is a staged research and software program. Estimate concrete implementation batches after their chapter inventories and domain contracts are reviewed; do not assign a fictitious deadline to “all of Mathematica” or to unresolved mathematical results.

## Evidence contract and definition of done

A tool response should carry a semantic status plus evidence kind, not a single green check. Proposed semantic statuses: `established`, `counterexample`, `computed`, `inconclusive`, `unsupported`, `invalid_input`. Evidence kinds are separate: `exact_algebra`, `numerical`, `formal_kernel`, `certificate_replay`, `external_comparison`, `source_derivation`. These names are design proposals; preserve existing CLI compatibility with adapters until reviewed.

Every evidence bundle records statement/item ID, source hash/location, assumptions, mathematical domain, operation/input digest, result/residual, trace/certificate, checker/backend version and digest, exact versus approximate arithmetic, tolerance/error budget, resource limits, code commit/dirty hashes, seed, commands, test results and counterexamples. Large runtime outputs live outside Git; versioned manifests identify them. No credentials enter task or evidence files.

A chapter item is implemented only when its reviewed statement has an executable construction or check, side conditions and reusable operations, positive/boundary/negative tests, replayable evidence and an independent reviewer. A “source inventory complete” label is a different state. Formal proof requires a stated foundation and checked proof artifact. An all-chapters audit cannot silently drop unsupported cases or average them into a passing score.

Compare backends critically: Cadabra and other CAS can share dependencies or algorithms. Agreement is corroboration, not automatically independent proof. Use independently constructed small examples, algebraic certificates and held-out mutations as additional checks. Numerical comparison must state conditioning, precision and tolerance.

## GPU strategy

Start on the RTX 3060, with reproducible Rust CPU implementations. Candidate work: batched finite-field polynomial evaluation, sparse coefficient arithmetic and modular reconstruction, batched tensor-component contractions/linear algebra, independent amplitude phase-space points, Monte Carlo/lattice experiments and graphics generation. Finite-field sampling alone is a probabilistic identity test; exact claims require a certified reconstruction/bound or independent exact checker. Floating-point GPU agreement never becomes exact proof.

Benchmark sorting/canonicalization and graph-heavy symbolic rewriting before moving them to CUDA; irregular control flow and expression-size variation may favor CPU execution. A Rust-owned backend interface may call narrowly scoped CUDA C++ kernels. Include transfer, allocation, warmup, host orchestration and checking in timings. Report crossover batch sizes, peak VRAM, correctness, speed and end-to-end workloads. Gate acceleration with the CPU reference and leave a functioning CPU path.

Native CUDA does not render directly in a browser frontend. TypeScript uses browser graphics or a Rust/WASM component; heavy GPU jobs execute in the local service. GPU scheduling must coordinate with local models and external CAS jobs rather than assuming the full 12 GB is always free.

## Agent dispatch and shared ownership

`tasks.json` contains foundational/capability tasks and 97 chapter inventory tasks. `HANDOFF.md` contains the dispatch procedure and a reusable assignment template. No harness-specific subagent is launched by writing these files.

- Allocate one owner to shared schema, conventions, coefficient types, registry and integration; chapter authors do not independently change them.
- Give parallel implementers isolated branches/worktrees and exclusive path ownership. A separate reviewer owns acceptance evidence. Only the integrator merges shared changes after dependency checks; do not auto-merge to master.
- Use existing installed harnesses and authorized providers. Capability-probe MCP support; a CLI fallback is valid only after testing. Log actual harness/model/tool versions, timeouts and usage. Do not assume all named tools support the same hooks, subagents or authentication.
- Explicit budgets govern concurrency, API usage, jobs, cancellation and GPU allocation. Never use multiple agents voting for the same answer as a verification method.
- New production units require bounded PDD intent plans and meaning approval. This plan does not widen ownership of the existing rotation prompt or change its conventions.

## PDD planning record

The original request is preserved in `request.txt`. A read-only planner run is preserved in `intent-plan.json`, ID `rust-mathematical-physics-tools-and-complete-sre-5f45bd76`. Its root-level discovery sees the existing nested Python rotation PDD files and recommends an existing-PDD change. That is discovery, not authorization to put this program under the rotation prompt. Before implementation, re-plan bounded Rust units against their actual scopes and preserve this request as the parent intent.

The user approved the two stories and architecture and authorized execution. New scoped Rust contracts and tools are being implemented with versioned prompt mappings; their actual verification is recorded in EXECUTION.md. No automatic `intent apply`, generated story contract or full regeneration is claimed without a successful recorded run. The planner's extracted summary is not authoritative over omitted user constraints.

## References

- Local primary inputs and hashes: `source-inventory.json`; all-chapter prerequisites and ownership: `chapters.csv`; candidate equation accounting: `equation-candidates.csv`.
- Lean kernel and Rust checker reference: https://lean-lang.org/doc/reference/latest/Elaboration-and-Compilation/#the-kernel and https://github.com/ammkrn/nanoda_lib . Rust checker applicability, version compatibility, axiom reporting and supported export formats require a dedicated evaluation; it is not assumed to formalize QFT automatically.
- Cadabra feature scope: https://cadabra.science/ . Symbolic scalar and wider algebra/geometry scope: https://www.sympy.org/en/index.html and https://doc.sagemath.org/html/en/reference/index.html . These inform feature families, not a claim of feature parity.
- Official Rust MCP SDK: https://github.com/modelcontextprotocol/rust-sdk . Pin an inspected version when implementing, then test actual harness compatibility.
- CUDA 13 architecture support: https://docs.nvidia.com/cuda/archive/13.0.0/cuda-toolkit-release-notes/index.html#deprecated-architectures .

## Third story Unified notebook UX

The user added [the unified notebook story](../../../user_stories/story__unified-notebook-workspace.md) and explicitly asked it be worked alongside the Srednicki tools. [UX-DELIVERY.md](UX-DELIVERY.md) defines a shared TypeScript browser/desktop workspace, candidate Rust/Tauri native shell, two themes, audio preferences and live harness integration. UX01–UX10 extend the central task graph; U01 tracks the complete delivery. This is an additive stream, preserving the approved mathematical-tooling and all-book scope.

## Game-like exploration (new story)

See [game exploration plan](game-exploration/README.md) and [user story](../../../user_stories/story__mathematical-physics-exploration.md). GX01–GX07 propose a distinct interaction experience, beginning with SU(2), alongside Notebook/Pixel and the ongoing Srednicki program. Story/design only; runtime implementation is pending.

## Capability-first architecture refinement

The [reusable research-library story](../../../user_stories/story__reusable-rust-research-library.md) and [engineering policy](REUSABLE-LIBRARY.md) govern all chapter work. Chapter/subsection batches organize source review and acceptance; production code belongs to reusable capabilities. New items first reuse an existing public API, or extend its owning library once. Textbook/paper-specific expected outputs stay in fixtures. Existing independent validators may intentionally use different mathematics; do not deduplicate away their independence.
