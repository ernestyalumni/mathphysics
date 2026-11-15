# Tooling roadmap (Rust first)

This extends the approved Rust mathematical-physics platform ([PRODUCT_INTENT.md](../../PRODUCT_INTENT.md), [REUSABLE-LIBRARY.md](../rust-physics-platform/REUSABLE-LIBRARY.md)). Same rules: capabilities are implemented once behind reusable public APIs; this program supplies **fixtures and requirements**, not a parallel engine. Language order: Rust, then justified C++/CUDA, then Python. Python (SymPy, the existing `.venv`) stays as an independent oracle.

**PDD status: approved by Ernest on 2026-11-15, with amendments** (review card at the end). Intents: `yang-mills-construction-audit-rust-lattice-gauge-d3c6da77` (original request, `request.txt`) and `yang-mills-program-amendment-local-network-agent-fb139368` (amendment, `request-2026-11-15-amendment.txt`). What `pdd intent apply` did and did not do is recorded in `HANDOFF.md`; never mark a unit as generated or synchronized unless that file says so. Each unit gets its own prompt, `.pddrc` context and `architecture.json` entry in its crate, as in `Rust/qft-symbolics`.

## Why these tools

Each one is tied to a place where the proof or the research plan needs a computation, and to a known answer that validates it.

| Unit | Crate | Capability | Validated against | Serves |
|---|---|---|---|---|
| T01 | `Rust/qft-symbolics` (extend) | Exact angular averages of monomials on S^{d−1}; radial shell integrals; rational × π^k results | ⟨n₁²⟩ = 1/4, ⟨n₁²n₀²⟩ = 1/24 on S³; shell integral 2π² log L; κ_G = 11C/(24π²). Oracle: `ledger/evidence/check_oneloop.py` | Audit #7; any one-loop computation |
| T02 | `Rust/lie-groups` (extend) | Compact simple Lie algebras: structure constants, Killing form, Casimir constant in a stated normalization, for su(N), so(N), sp(N); exceptional algebras from Cartan data | C_𝔤 = 2N² for su(N) with ⟨X,Y⟩ = −(1/N) Re tr(XY); dual Coxeter numbers | Audit #7; every group-dependent constant |
| T03 | `Rust/qft-symbolics` (extend) | Free-field Wick engine in position space for curvature composites of dim 𝔤 Maxwell fields | ⟨F²(x) F²(0)⟩ = 3 dim 𝔤 / (π⁴|x|⁸) in the proof's normalization (take the normalization from `YM/docs/examples` LeadingTrace worked example) | Audit #7; trace anomaly checks |
| T04 | new `Rust/lean-graph` | Read Armstrong's `data/*.tsv` and (later) a Lean export; answer "what does X depend on", "which unit owns X", "shortest path from the main theorem to X" | Reachable-set counts published by Armstrong and by OpenAI's `verification/declaration-usage.tsv` | Every reading and audit task; lets weak agents navigate a million lines |
| T05 | new `Rust/lattice-gauge` | Lattice geometry in d dimensions: sites, links, plaquettes, cubes, tori, twisted boundary conditions, co-boundary operators with exact integer arithmetic | d² = 0; counts 4V links, 6V plaquettes; `3V+1` free links in a tree gauge (the proof's determinant budget) | Everything below |
| T06 | `Rust/lattice-gauge` | SU(2) (unit quaternions) and SU(3) link variables; Wilson action; heat-bath and over-relaxation updates; reproducible RNG streams | Published SU(2) plaquette averages at standard β; strong-coupling series for ⟨P⟩ | YM-N01 |
| T07 | `Rust/lattice-gauge` | Observables: Wilson loops, Creutz ratios, plaquette–plaquette correlators, Polyakov loops; jackknife errors | Published string-tension estimates in the scaling window | YM-N01, N04 |
| T08 | `Rust/lattice-gauge` | Centre twists and vortex containers; ratio Z_tw/Z by a chain of intermediate ensembles | 2D Z₂ and U(1) cases with closed forms; published SU(2) twist free energies | YM-N02 (the key experiment) |
| T09 | new `Rust/haar` or in `lie-groups` | Exact Haar integrals of polynomial link functions: SU(2) closed forms, Weingarten calculus for SU(N); character expansion coefficients of the Wilson weight | ∫ tr U tr U† = 1; SU(2) character coefficients I_{2j+1}(β)/I_1(β); low orders of the strong-coupling expansion | R3 baseline; exact small-lattice work |
| T10 | `Rust/lattice-gauge` | Exact evaluation on tiny lattices: enumerate Z₂ and Z_N gauge theories fully; SU(2) by character expansion; exact rational or interval output | 2D lattice gauge theory closed forms; duality of 3D Z₂ gauge theory with Ising | YM-G03, YM-N03 |
| T11 | `Rust/qft-symbolics` (extend) | Noncommutative differentiation of products of group elements along Lie-algebra directions: nth variation of tr(U₁ e^{B} U₂ …) as an ordered-algebra expression | N1 eqs. (2.10)–(2.12) reproduced for O(3) from the paper's LaTeX | YM-G03 (the key symbolic computation) |
| T12 | later | Interval arithmetic and certified bounds for polymer-expansion convergence criteria (Kotecký–Preiss sums) | The criterion of summary §"Cluster expansion": 4(2B+1)⁴γ/δ ≤ q < 1 | R4; any computer-assisted step |
| T13 | later, CUDA | GPU Monte Carlo for T06–T08 on the RTX 3060 | Bitwise or statistical agreement with the CPU path | Scale for YM-N02 |

Order: T01, T02, T05 first (small, exact, unblock the rest). T04 in parallel (independent). T11 and T10 next (they decide route R1). T06–T08 after T05. T03, T09 when needed. T12, T13 deferred; do not start unasked.

## What "replacement for Sage, SymPy, Cadabra2, Mathematica" gains from this

- T01, T03, T11 are symbolic-algebra capabilities the platform does not have: exact integration over spheres, position-space Wick calculus, and calculus on noncommuting group-valued expressions. T11 in particular is something Cadabra and SymPy do badly.
- T02 and T09 are the start of a representation-theory layer (LiE, Sage's `WeylCharacterRing`).
- T04 makes a formal library queryable by agents, which is the missing link between the CAS and Lean.

## Engineering rules for workers

- Read `PDD.md` at the repository root and the crate's own `PDD.md` and `README.md` first.
- `export CARGO_TARGET_DIR=/tmp/mathphysics-rust-target`. No build products in the repository.
- No new crates.io dependency without a line in the packet saying it is approved. Offline, locked builds are the norm here (`--offline --locked`).
- Every public function has: a positive test against an independent oracle, a negative control that must fail, and a doc comment stating normalization conventions.
- Exact arithmetic wherever the mathematics is exact. Floating point only in Monte Carlo and clearly marked.
- snake_case, spelled-out names. No abbreviations like `numel`.
- Monte Carlo runs record seed, parameters, code revision and wall time in their output file. Outputs are not committed; a summary receipt is.

## Rust first, CUDA C++ where it pays

Ernest asked (2026-11-15) whether Rust-first tools with CUDA C++ for the parallel parts will work. Yes, with the split below. The rule: every capability has a Rust CPU reference path behind a public API; a CUDA path is added behind the same API only for a workload that is data-parallel **and measured** to be the bottleneck.

| Workload | Shape | Rust CPU | CUDA C++ |
|---|---|---|---|
| Symbolic algebra: sphere averages, Wick contractions, calculus on words of group elements (T01, T03, T11) | Irregular trees, big rationals, hashing | Yes. Parallelize across independent terms with `std::thread::scope` | **No.** Branchy, allocation-heavy code does not vectorize |
| Lie-algebra constants, Haar integrals, character coefficients (T02, T09) | Small exact linear algebra | Yes | No |
| Lattice geometry, cochains (T05) | Integer index arithmetic | Yes | Shared index conventions only |
| Monte Carlo: heat bath and over-relaxation (T06), observables (T07) | Checkerboard: all links of one direction and one parity are conditionally independent | Reference path, small lattices | **Yes.** The classic GPU workload. An SU(2) lattice of side 16 is 4 MB of link data |
| Twist ratios by a chain of intermediate ensembles (T08, YM-N02) | Many replicas (β values, twist strengths, seeds), each a small lattice | Slow | **Yes, the best case:** batch the replicas, one kernel launch updates all of them |
| Exact small lattices by transfer matrix (T10) | Matrix-free power iteration over 2^k slice states, bit operations | First version, threaded | Candidate later, if the state space is the limit |
| Interval arithmetic, certified bounds (T12) | Needs controlled rounding | Yes | No. Rounding control on GPUs is not dependable |
| Declaration-graph queries (T04) | Graph traversal on 77k nodes | Yes | No |

Engineering constraints:

- The same counter-based random number generator on both paths, keyed by (seed, replica, site, direction, sweep), so a GPU run can be compared with the CPU run configuration by configuration on a small lattice. Acceptance of any CUDA unit: that comparison, plus a measured speed-up on the YM-N02 workload written into the receipt.
- Two GPUs are present: an RTX 3060 (12 GB) and a GTX 980 Ti (6 GB). **The local Qwen server also uses the RTX 3060.** Monte Carlo and local-model workers contend for it; do not run both without checking `nvidia-smi`. Whether the installed CUDA toolkit still supports the 980 Ti must be checked before planning on it.
- CUDA C++ lives in its own directory with a C ABI; Rust calls it through a thin `extern "C"` layer. No Rust CUDA crates, no new crates.io dependencies, without Ernest's approval.
- T13 stays deferred until T06–T08 exist on the CPU and a profile shows where the time goes.

## Tools as the discovery loop

Ernest's direction (2026-11-15): do not rely completely on Lean; win by having new tools. The division of labour between the tools:

| Stage | Tool | Output |
|---|---|---|
| Conjecture | Rust symbolic and exact engines (T10, T11) | An identity or inequality tested exactly on small lattices, with a certificate file (exact rationals) |
| Falsification | Rust and CUDA Monte Carlo (T06–T08) | Numbers with error bars, on sizes the exact engine cannot reach |
| Proof | A strong model plus a human-readable argument in the mono-notes | LaTeX, reviewed by a second model family |
| Certification | Lean, in Docker, on top of `OpenAI-YM-20261028` | A kernel-checked statement |

A result is not called proved before the last row. It can be called **verified on lattices up to size n** after the first.

## Review card (approved by Ernest 2026-11-15, with amendments)

- **Scope.** Repository `mathphysics`. Documentation under `docs/plans/yang-mills-mass-gap/`. Mathematical account in `documents/mono-notes/topics/yang-mills-lattice.tex`. New crates `Rust/lattice-gauge`, `Rust/lean-graph` (and possibly `Rust/haar`). Extensions to `Rust/qft-symbolics` and `Rust/lie-groups`. The four Yang–Mills repositories and `openai/math` are read-only inputs (`SOURCES.md`).
- **Durable requirements.**
  - The program MUST keep audit, tooling and mass-gap research as separate workstreams with separate evidence.
  - Every tool MUST be validated against an independent oracle and a named number or identity from the proof or the literature.
  - Mathematical capabilities MUST live in the reusable libraries; Yang–Mills material supplies fixtures only.
  - No result MUST be called proved unless it is kernel-checked or accepted by a named human expert.
  - The tools MUST be Rust first. CUDA C++ MUST be used only where a measured data-parallel workload benefits.
  - The program MUST NOT rely completely on Lean: the Rust symbolic and exact-computation tools are the discovery and falsification loop; Lean is the final certification.
  - Lean MUST run in a Docker image, with the toolchain and compiled build stored on the Expansion drive.
  - **Agents MUST be allowed to talk to other local agents, including agents running on other machines on the same local network** (amendment; see `06-AGENT-PROTOCOL.md` §10).
  - Agents MUST NOT post, publish, upload or message **outside the local network** without Ernest's explicit approval. One standing exception: an agent that Ernest has given GitHub access MAY contribute to Hairer's public repository `hairer/humanise-yang-mills` on its own branch, following that repository's coordination rules (`06-AGENT-PROTOCOL.md` §11).
  - Agents MUST NOT edit the pinned source repositories.
  - Work MUST be handed off through files in this directory so that any agent can resume it.
- **For example**, a Qwen3.5 worker on a second machine on the local network claims packet YM-T05 through the shared tracker, pushes its branch to the canonical mathphysics checkout over the local network, and a Codex reviewer on the desktop accepts or rejects it.
- **Worker harnesses.** Pi, Hermes Agent, Codex, Claude Code, Grok Build; local Qwen3.5 models (llama.cpp server profiles in `claw-work/Deployments/Scripts/llama-cpp-server`) when Ernest is away. Tier mapping in `06-AGENT-PROTOCOL.md` §12.
- **Proposed prompt units.** T01–T11 above, one prompt each, in the crate that owns the capability.
- **Story.** One user story, "Yang–Mills audit and mass-gap research workspace", with acceptance: (a) a newcomer agent can pick a packet and finish it using only the packet and the protocol; (b) every numerical constant quoted in `01-WHAT-OPENAI-PROVED.md` §1 is recomputed by a Rust test; (c) the twist-ratio experiment runs end to end and writes a receipt.
- **Decisions taken by Ernest on 2026-11-15.** (1) Meaning approved; run `pdd intent apply`. (2) Contributions to Hairer's public repository are allowed for an agent with GitHub access. (3) Harnesses as listed above. (4) Lean in Docker; install on the Expansion drive.
