# YM-T09 — Matrix integration: tutorial, measures and reusable exact/numerical tools

Status: ready; no implementation worker started. Tier T2–T3; independent mathematical and software review required. This activates the existing Haar-integration roadmap ID, with a beginner-facing first milestone. Later general SU(N) Weingarten/character integration remains future work, not promised by this packet.

## Intent and first assignment

Ernest wants to understand what it means to integrate a matrix, which measure to use, and how explicit hand calculations become reusable symbolic and numerical Rust tools, followed by CUDA C++ where useful. He is reading the Wilson-measure lesson and Gattringer–Lang §3.1, with Srednicki as continuum background. Preserve that learning sequence; explain each new term and show intermediate steps.

First deliverable: phase A below and one bounded phase-B exact integration unit. Submit these for independent review before widening implementation. Phases C and D are concrete follow-on milestones; do not label the entire task complete after A/B. General group-word integration may use YM-T15 once its interface is reviewed, but the U(1)/quaternion-coordinate first unit does not depend on it.

The original request is in [request-2026-11-16-matrix-integration.txt](../request-2026-11-16-matrix-integration.txt); its read-only planner output is alongside it. The planner's suggested trace/tracker targets are lexical matches, not ownership decisions. This is a manual task specification, not a generated PDD contract or an applied intent.

## A. Mathematical write-up and checked calculations

Read `documents/mono-notes/STYLE.md`, the existing lattice topic and the Wilson reading guide first. Add a coherent prerequisite subsection/topic, cross-reference it from the Wilson measure, and integrate through `master.tex`; `screen.tex` already inherits the master. Keep software/agent discussion in a companion note. Build both editions, check citations/references, and inspect representative pages. Coordinate shared master/topic edits with the tracker.

Distinguish these operations explicitly:

1. **Integrate a matrix-valued function over another variable.** For finite matrices, integrate each integrable entry. Work through
   `A(t) = [[t,t²],[1,2t]]`, `0≤t≤1`, obtaining `[[1/2,1/3],[1,1]]`.
   Explain when trace commutes with integration. Show that determinant generally does not: `A(t)=t I₂` gives `det(∫A)=1/4`, whereas `∫det A=1/3`.
2. **Integrate a function over a space of matrices.** Specify the domain, its independent real coordinates and the measure before writing an integral. For real 2×2 matrices, `dM=∏ᵢⱼdMᵢⱼ`; on the unit four-cube `∫tr(M)dM=1`. On the whole space,
   `∫exp[-tr(MᵀM)/2]dM=(2π)²`, by four scalar Gaussian integrals. Distinguish Lebesgue measure from its normalized Gaussian probability density. Briefly compare complex matrices (2n² real coordinates), Hermitian matrices (n²), and constrained group matrices. There is no universal matrix measure.
3. **Integrate over a compact matrix group with normalized Haar measure.** Define left/right invariance and total mass one; distinguish this measure from flat independent entries and from the subsequent Wilson density. Show the U(1) Fourier identity `∫exp(ikθ)dθ/(2π)=δₖ₀` for integer k.
4. **Calculate SU(2) examples in full.** Use the explicitly stated Pauli convention
   `U=q₀I+i(q₁σ₁+q₂σ₂+q₃σ₃)`, `q∈S³`.
   Explain why normalized Haar measure is `dΩ₃/(2π²)`. With `q₀=cosχ` and `q⃗=sinχ n(θ,φ)`, derive
   `dh=(1/(2π²)) sin²χ sinθ dχ dθ dφ`, for `χ,θ∈[0,π]`, `φ∈[0,2π)`; explain endpoint coordinate degeneracy and check normalization. Work out `∫U dh=0`, `∫qαqβ dh=δαβ/4`, `∫q₀⁴ dh=1/8`, `∫q₀²q₁² dh=1/24`, `∫|tr U|²dh=1`, and `∫|U₁₁|²dh=1/2`. Also derive `∫U₁₁U₂₂ dh=1/2` and `∫U₁₂U₂₁ dh=-1/2`: determinant-one constraints mean U(N) phase-balance selection rules cannot simply be copied to SU(N).
5. **Return to the Wilson measure.** Explain product Haar over links, scalar/matrix-valued observables, the Boltzmann density, Z and normalized expectations. Under normalized product Haar, `Z(β=0)=1`. Existing Z₂ density-of-states histograms are raw counts, so divide the count partition function by `2^link_count` to compare absolute Z; ratios may cancel that factor. Finite-lattice integration is not a continuum construction or a mass-gap proof.

## B. Exact symbolic Rust first

Inspect existing crate boundaries and PDD mappings before selecting the owning API. Read the reusable Rust research-library story. Reuse `Rust/qft-symbolics/src/sphere.rs` (exact rational sphere monomial averages), its exact coefficient/budget infrastructure and `Rust/lie-groups/src/su2.rs` where appropriate. Do not duplicate the sphere kernel or build a chapter-specific calculator.

First supported domains: finite U(1) Laurent/Fourier polynomials and SU(2) polynomials in quaternion coordinates with exact coefficients. Decide whether these are two prompt units or one coherent bounded unit after inspection; land one reviewed unit first. Make domain, normalized measure and integrand type explicit. Integrate matrix-valued inputs entrywise through a generic adapter where needed. Exact output must stay exact; symbolic expressions outside the supported grammar, such as general exponentials, must return an explicit unsupported/formal result rather than a guessed simplification.

SU(2) coordinate convention is load-bearing: the existing scalar-first `quaternion_matrix` maps to `q₀I-i q⃗·σ`, while the lesson above uses plus-i. Verify the sign bridge and multiplication order with explicit Pauli matrices before sharing Hamilton-product helpers. Do not silently change the existing API convention.

Specify degree, term, coefficient-growth and work limits; reject invalid exponents/domains and expose budget exhaustion. Reuse existing resource errors where possible. Tests must include fresh input polynomials/matrix entries and independent hand-derived or scalar-quadrature oracles, not just calls to the same sphere helper. Include odd/cross moments, the nonconjugate SU(2) controls, U(1) nonzero modes, linearity, normalization, trace/entry agreement, unsupported input, invalid dimensions and budget failures. A fixture or textbook equation ID must never select a production answer.

Use the actual owning crate's prompt-first PDD workflow: inspect installed CLI help, preserve request/acceptance, run read-only intent planning there, update the bounded `.prompt`/mapping before code, characterize affected behavior and run scoped checks. Do not run the root provider-backed architecture generator, alter authentication or incur implicit extra-model usage. Do not claim regeneration/synchronization succeeded without evidence.

## C. Numerical Rust with honest error estimates

Implement one reusable CPU estimator/quadrature interface over the same domains and integrands. Keep sample generation, weights, integrand evaluation and reduction distinct. Record sample count, seed/algorithm, precision, error estimate and assumptions; handle zero/invalid samples, nonfinite values and work limits explicitly. For independent random sampling, uncertainty must reflect the actual variance, not an arbitrary tolerance.

An SU(2) option is four independent standard normal coordinates divided by their Euclidean norm: justify the rotationally invariant S³ distribution and handle the zero-norm event. Deterministic quadrature must use the Haar Jacobian. Uniformly sampling χ is a deliberately wrong sampler: it predicts `E[q₀²]=1/2`, instead of Haar's `1/4`. Use it as a negative control. Verify analytic moments independently and report convergence with reproducible experiments; do not present numerical agreement as a proof. Keep gauge-ensemble heat-bath work under YM-N01 rather than expanding this packet silently.

## D. CUDA C++ only after a measured workload

Once phase C yields a representative batch/profile, specify a reusable C++ CUDA evaluation/reduction backend with explicit data layout, ownership, precision and resource errors, retaining a stable domain/integrand contract with Rust. Benchmark end-to-end against the CPU including initialization, transfers and reductions; record hardware/toolchain, workload sizes and crossover. Use identical caller-supplied samples for CPU/GPU parity, plus analytic mean checks and justified tolerances. Do not claim a speedup from a successful kernel launch. The Mac Mini may coordinate or develop CPU work; execute CUDA on a verified NVIDIA CUDA host. If no useful workload or hardware is available, record D as deferred with the evidence and remaining steps.

## Sources, deliverables and completion

- Main reading remains Gattringer–Lang, *Quantum Chromodynamics on the Lattice*, §3.1: [publisher](https://link.springer.com/book/10.1007/978-3-642-01850-3). Only publisher contents/previews were checked in the preceding lesson; verify exact chapter claims against a legitimate copy before attributing equations.
- Optional agent reference, not another required book: Elizabeth Meckes, *The Random Matrix Theory of the Classical Compact Groups*, chs.1–2, [author-hosted manuscript](https://case.edu/artsci/math/mwmeckes/elizabeth/Haar_book.pdf). Check the precise group before applying a moment identity.
- Preserve mathematical sources and CAS comparisons. Lean can provide statement/proof context; Rust symbolic algebra and numerical estimates have separate evidence and do not replace a proof kernel.

Deliver: the manuscript prerequisite with derivations; `notes/YM-T09-matrix-integration.md` with API ownership, mathematical-to-code map, commands/results and unsupported cases; bounded mapped source prompts, reusable Rust APIs, maintained examples and independent tests; phase-C evidence and a measured phase-D decision as separate milestones. A/B may be submitted independently; mark the full task done only with each milestone's disposition and independent review. Report exact branch/commit, checks and remaining work on the shared tracker for Fable. No public push or LAN identifiers in source.
