# Board

Status values: `ready` (can be claimed now), `blocked:<reason>`, `claimed`, `in-review`, `done`. The tracker thread is authoritative for who holds a claim; this file is the summary. Updated by the coordinator or strategist only.

Last updated: 2026-11-16 (strategist, session 3). Local GPU model serving on the desktop is blocked: `notes/local-models.md` §3a.

## Wave 1 — ready now (reading and analysis; no installs, no code in the product crates)

Tracker items: YM-R01 = local:1, R02 = 2, R03 = 3, R04 = 4, R05 = 5, R08 = 6, A03 = 7, A04 = 8.

| ID | Title | Tier | Status | Depends on | Parallel-safe with |
|---|---|---|---|---|---|
| [YM-R01](tasks/YM-R01.md) | Digest N1 (2D O(n) exponential decay) | T3 | ready | — | all |
| [YM-R02](tasks/YM-R02.md) | Digest N2 §2: defect ⇒ gap | T3 | ready | — | all |
| [YM-R03](tasks/YM-R03.md) | Digest N2 §3, §9: mixing and transport | T3 | ready | — | all |
| [YM-R04](tasks/YM-R04.md) | Digest N3 §8, §§12–15: gap and uniqueness | T3 | ready | — | all |
| [YM-R05](tasks/YM-R05.md) | Define the unit-scale class 𝒞(β̄) from the Lean | T3 | ready | — | all |
| [YM-R08](tasks/YM-R08.md) | Infrared end of the Lean: IV.1, V.1, V.3 | T2–T3 | ready | — | all |
| [YM-A03](tasks/YM-A03.md) | Red-team the formal statement (run twice, independently) | T3 | ready | — | all |
| [YM-A04](tasks/YM-A04.md) | Independent transcription of the public theorem | T3 | ready | — | not the same agent as A03 |

## Wave 1b — ready (reading public sources is allowed by the workspace rules; send nothing private in a request)

| ID | Title | Tier | Status |
|---|---|---|---|
| [YM-R06](tasks/YM-R06.md) | Tomboulis 2007 and Ito–Seiler | T2 | ready |
| [YM-R07](tasks/YM-R07.md) | Literature baseline table | T2 | ready |
| [YM-W01](tasks/YM-W01.md) | Upstream watch (recurring) | T1 | ready |
| [YM-W02](tasks/YM-W02.md) | Jaffe–Witten text and Clay rules, quoted | T1–T2 | ready |

## Wave 2 — after wave 1 deliverables

| ID | Title | Tier | Status | Depends on |
|---|---|---|---|---|
| [YM-G02](tasks/YM-G02.md) | Defect ⇒ gap lemma for lattice gauge theory | T3 | ready (see above) | — |
| [YM-G03](tasks/YM-G03.md) | Gauge fourth variation: is there a positive square? (run twice) | T3 | blocked:YM-R01 | R01 |
| [YM-G04](tasks/YM-G04.md) | Embedded Z₂ gauge theory and frustration | T3 | blocked:YM-R01 | R01 |
| [YM-G05](tasks/YM-G05.md) | Warm-up lemma: clustering collapses the phases | T3+T2 | blocked:YM-R08 | R08 |
| [YM-G06](tasks/YM-G06.md) | Feasibility of a computer-assisted finite-size criterion | T3 | blocked:YM-R02,R07 | R02, R07 |
| [YM-A05](tasks/YM-A05.md) | Spot-check outlines III.9, III.12, IV.1, V.1, V.3 | T2 | blocked:YM-R05,R08 | R05, R08 |
| [YM-G01](tasks/YM-G01.md) | Problem sheet (milestone G0) | T4+T3 | blocked:YM-R02,R03,R05,R08 | R02, R03, R05, R08 |

## Submitted by the strategist on 2026-11-15 — need review (reviewer: not Claude if possible)

Tracker items: YM-O01 = local:9, T01 = 10, T02 = 11, T05 = 12, G02 = 13.

| ID | Title | Tier | Status | What to review |
|---|---|---|---|---|
| [YM-O01](tasks/YM-O01.md) | Lean in Docker + compiled build | T1–T2 | in-review; packet step 4 open | `notes/YM-O01-lean-setup.md`. Rerun the axiom print; run the three verification scripts. |
| [YM-T01](tasks/YM-T01.md) | Exact sphere averages, one-loop slope | T2 | in-review | `Rust/qft-symbolics/src/sphere.rs`, `tests/sphere_averages.rs` |
| [YM-T02](tasks/YM-T02.md) | Lie algebra constants | T2 | in-review; su(N) only, so(N)/sp(N)/exceptional open | `Rust/lie-groups/src/special_unitary.rs`, `tests/special_unitary_casimir.rs` |
| [YM-T05](tasks/YM-T05.md) | Lattice geometry, exact | T2 | in-review; twisted boundary conditions and integer (not mod 2) coboundaries open | `Rust/lattice-gauge` (`torus`, `cochain`) |
| [YM-G02](tasks/YM-G02.md) | Defect ⇒ gap lemma | T3 | **ready** (no longer blocked: proofs written) | Mono-notes Props. `ym:prop:purity`, `ym:prop:defect`, `ym:prop:conditional-gap`: check every line, compare with N2 §2, then attack the recurrence. Red-team pass required. |

## Added 2026-11-16 (strategist, session 6)

| ID | Title | Tier | Status | Tracker |
|---|---|---|---|---|
| [YM-W04](tasks/YM-W04.md) | Read the Faizal–Shabir mass-gap claim (arXiv:2606.19362), locate its crossover step | T3 | ready | local:19 |
| YM-G02 addendum | Review the transfer-positivity sketch (`notes/spectral-theory-study-route.md` §3, claim C23) together with the G02 proofs | T3, non-Claude | ready | local:13 |

## Ready now, no approvals needed (all decisions were taken on 2026-11-15)

| ID | Title | Tier | Status |
|---|---|---|---|
| [YM-A01](tasks/YM-A01.md) | Clean rebuild on Linux, in Docker | T2 | ready (after O01 review); many hours of CPU |
| [YM-A02](tasks/YM-A02.md) | Comparator with real isolation; second kernel | T2–T3 | blocked:YM-A01 |
| [YM-T04](tasks/YM-T04.md) | `lean-graph` query tool | T2 | ready |
| [YM-T10](tasks/YM-T10.md) | Exact small-lattice evaluation | T2–T3 | ready; first slice exists (`Rust/lattice-gauge/src/z2_gauge.rs`): extend to Z_N, transfer-matrix method, SU(2) characters |
| [YM-T11](tasks/YM-T11.md) | Calculus on words of group elements | T3 | blocked:YM-R01 |
| [YM-T14](tasks/YM-T14.md) | Micro-task runner for local 9B models (`referee` and `index` types first; see `notes/local-models.md`) | T2 | merged (64c2c2b) |
| [YM-T16](tasks/YM-T16.md) | Weight edit of the local 9B referee (claw-work `Rust/weight-space`, phase 2) | T2 | done: no difference established (`notes/YM-T16-weight-space-referee.md`) |
| [YM-T17](tasks/YM-T17.md) | Proof steps as numeric tests written by a local model | T2 | merged (64c2c2b); verdict: local models are not reliable judges (best: 3-of-3 agreement, 0/54 false alarms, 15% recall; YM-O04 §5) |
| [YM-T18](tasks/YM-T18.md) | Verified samplers (`ymlib`) and agreement voting for numeric tests | T2 | merged (fb9dcaa) |
| [YM-N01](tasks/YM-N01.md) | SU(2) Monte Carlo baseline | T2 | merged (a18a4dc): Rust CPU + CUDA C++; 2D exact, C_V maximum, Creutz ratios, glueball and chi_t as cross-checks against published values (YM-O04 §3); also U(1), monopoles, vortices, topology |
| YM-F01 | Lean formalisation of the conditional-gap chain (casimir-sun, purity, defect telescope, conditional gap, rectangle recurrence); statements read against the tex by hand | T2–T3 | claimed: claude-opus55-q7nd, `ym/YM-F01-claude-opus55-lean`; recurrence proved (`Rectangle.lean`, 2ed351b); next proposed: certified seed Δ ≤ 192βn⁴ (`notes/YM-strategy-2026-11-21.md`) |
| YM-T19 | Algebraic topology library (homology, cup/cup-1, Pontryagin square, persistence) and vortex/monopole/instanton topology | T2–T3 | merged (4fc88e0); corrections on `ym/YM-N03b-claude-opus55-recursion` |
| [YM-N02](tasks/YM-N02.md) | Twist ratios on vortex containers | T3/T2 | blocked:G04; low priority (estimator ready: `examples/twist_z.rs`) |
| [YM-N03](tasks/YM-N03.md) | Defect criterion in solvable models | T2 | SU(2) n = 2, 4 and U(1) n = 2, 4 measured; recurrence proved (YM-O03 §1); see `notes/YM-N03-defect-measurement.md` |
| [YM-N04](tasks/YM-N04.md) | Z₂ monopole density at weak coupling | T2 | blocked:G04 (monopole tools in `Rust/lattice-gauge/src/monopole.rs`) |
| YM-R06, R07, W01, W02 | Literature and watch packets (wave 1b above) | T1–T2 | ready |


## Suggested first launch

Eight agents in parallel, one per wave-1 row, T3 where marked. Put two **different** model families on YM-A03. Give YM-R01 to the strongest model available: G03 and G04, which decide the primary route, both wait on it. One coordinator pass after the first SUBMIT, then reviewers from a different family for each.

## Dependency sketch

    R01 ─┬─> G03 ──> (T11) ──> G3 verdict on route R1
         └─> G04 ──> N02, N04 ──> G4
    R02 ───> G02 ──> N03
    R02,R03,R05,R08 ──> G01 (problem sheet, G0) ──> conditional theorem G5
    R08 ───> G05 (first Lean extension)
    O01 ──> A01 ──> A02 (G1 with A03, A04)

## Continuation requested 2026-11-15 — Mac first, then source reading

| ID | Title | Tracker | Status |
|---|---|---|---|
| [YM-O02](tasks/YM-O02.md) | First Mac Mini connection and Fable continuity | local:14 | done: technical SSH/tracker acceptance independently checked; original Pi model acknowledgement not observed |
| [YM-R09](tasks/YM-R09.md) | Srednicki-to-lattice source and mono-notes reading map | local:15 | ready; no worker claimed |

Codex report for Fable: [CODEX-HANDOFF-2026-11-15.md](CODEX-HANDOFF-2026-11-15.md). These additions preserve the existing review-first research priorities. Private machine values are outside the public repository.

Daily source watch: [YM-W03](tasks/YM-W03.md), local:16. Collector active at/after09:00 Los Angeles; AI review queued, not automatic. See [watch/README.md](watch/README.md) and [review baseline](notes/YM-W03-baseline-2026-11-16.md). Initial source coverage8/10 with explicit403/429 failures.

Wilson-measure study handoff: [YM-T15](tasks/YM-T15.md), local:17, ready for scoped implementation planning; no worker claimed. Start with [the reading guide](notes/wilson-measure-reading-guide.md) and [the bounded symbolic plan](notes/wilson-symbolics-plan.md). Implementation contract is draft, not a generated/mapped unit.


Matrix-integration learning/tooling handoff: [YM-T09](tasks/YM-T09.md), local:18, ready; no worker claimed or implementation started. [Pasteable agent prompt](prompts/matrix-integration-worker.md). First submission: entrywise versus matrix-domain integration tutorial and one bounded exact Rust unit; numerical Rust and profiled CUDA C++ are separate follow-on milestones. Reuses existing exact sphere moments; general word integration may follow T15 without blocking the first coordinate-polynomial unit.


## Jaffe–Witten extensions — 21 November 2026 UTC

These extend the research agenda, not the Millennium prize conditions. [Roadmap](08-EXTENSIONS-ROADMAP.md), [shared worker prompt](prompts/yang-mills-extensions-worker.md). Core G02 reviews retain priority. No new proof worker claimed.

| Packet | First bounded milestone | Tracker | Status |
|---|---|---|---|
| [YM-E01](tasks/YM-E01.md) | Isolated particle mass shell | local:22 | ready: definition/source milestone |
| [YM-E02](tasks/YM-E02.md) | Specified confinement criterion | local:23 | ready: definition/source milestone |
| [YM-E03](tasks/YM-E03.md) | Gauge–matter existence / ultraviolet assumptions | local:24 | ready: definition/source milestone |
| [YM-E04](tasks/YM-E04.md) | Matter infrared phases / chiral symmetry | local:25 | ready: definition/source milestone |
| [YM-E05](tasks/YM-E05.md) | Curved-manifold existence specification | local:26 | ready: definition/source milestone |
