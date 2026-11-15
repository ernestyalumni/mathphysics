# Board

Status values: `ready` (can be claimed now), `blocked:<reason>`, `claimed`, `in-review`, `done`. The tracker thread is authoritative for who holds a claim; this file is the summary. Updated by the coordinator or strategist only.

Last updated: 2026-11-14.

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

## Wave 1b — ready once Ernest allows network reads for public papers

| ID | Title | Tier | Status |
|---|---|---|---|
| [YM-R06](tasks/YM-R06.md) | Tomboulis 2007 and Ito–Seiler | T2 | blocked:network-ok |
| [YM-R07](tasks/YM-R07.md) | Literature baseline table | T2 | blocked:network-ok |
| [YM-W01](tasks/YM-W01.md) | Upstream watch (recurring) | T1 | blocked:network-ok |
| [YM-W02](tasks/YM-W02.md) | Jaffe–Witten text and Clay rules, quoted | T1–T2 | blocked:network-ok |

## Wave 2 — after wave 1 deliverables

| ID | Title | Tier | Status | Depends on |
|---|---|---|---|---|
| [YM-G02](tasks/YM-G02.md) | Defect ⇒ gap lemma for lattice gauge theory | T3 | blocked:YM-R02 | R02 |
| [YM-G03](tasks/YM-G03.md) | Gauge fourth variation: is there a positive square? (run twice) | T3 | blocked:YM-R01 | R01 |
| [YM-G04](tasks/YM-G04.md) | Embedded Z₂ gauge theory and frustration | T3 | blocked:YM-R01 | R01 |
| [YM-G05](tasks/YM-G05.md) | Warm-up lemma: clustering collapses the phases | T3+T2 | blocked:YM-R08 | R08 |
| [YM-G06](tasks/YM-G06.md) | Feasibility of a computer-assisted finite-size criterion | T3 | blocked:YM-R02,R07 | R02, R07 |
| [YM-A05](tasks/YM-A05.md) | Spot-check outlines III.9, III.12, IV.1, V.1, V.3 | T2 | blocked:YM-R05,R08 | R05, R08 |
| [YM-G01](tasks/YM-G01.md) | Problem sheet (milestone G0) | T4+T3 | blocked:YM-R02,R03,R05,R08 | R02, R03, R05, R08 |

## Needs Ernest's decision before it can start

| ID | Title | Tier | Decision needed |
|---|---|---|---|
| [YM-O01](tasks/YM-O01.md) | Lean toolchain + compiled build | T1–T2 | Install on Expansion drive (~15 GB, drive is 97% full) |
| [YM-A01](tasks/YM-A01.md) | Clean rebuild on Linux | T2 | after O01; a day of CPU |
| [YM-A02](tasks/YM-A02.md) | Comparator with real isolation; second kernel | T2–T3 | after A01 |
| [YM-T01](tasks/YM-T01.md) | Exact sphere averages, one-loop coefficient | T2 | PDD intent approval |
| [YM-T02](tasks/YM-T02.md) | Lie algebra constants | T2 | PDD intent approval |
| [YM-T04](tasks/YM-T04.md) | `lean-graph` query tool | T2 | PDD intent approval |
| [YM-T05](tasks/YM-T05.md) | Lattice geometry, exact | T2 | PDD intent approval |
| [YM-T10](tasks/YM-T10.md) | Exact small-lattice evaluation | T2–T3 | after T05 |
| [YM-T11](tasks/YM-T11.md) | Calculus on words of group elements | T3 | PDD intent approval; after R01 |
| [YM-N01](tasks/YM-N01.md) | SU(2) Monte Carlo baseline | T2 | after T05 |
| [YM-N02](tasks/YM-N02.md) | Twist ratios on vortex containers | T3/T2 | after N01, G04 |
| [YM-N03](tasks/YM-N03.md) | Defect criterion in solvable models | T2 | after R02, T10 |
| [YM-N04](tasks/YM-N04.md) | Z₂ monopole density at weak coupling | T2 | after N01, G04 |

Not yet written as packets (see `05-TOOLING-ROADMAP.md`): T03 free-field Wick engine, T06–T08 (covered by N01, N02 for now), T09 Haar integrals, T12 interval arithmetic, T13 CUDA.

## Suggested first launch

Eight agents in parallel, one per wave-1 row, T3 where marked. Put two **different** model families on YM-A03. Give YM-R01 to the strongest model available: G03 and G04, which decide the primary route, both wait on it. One coordinator pass after the first SUBMIT, then reviewers from a different family for each.

## Dependency sketch

    R01 ─┬─> G03 ──> (T11) ──> G3 verdict on route R1
         └─> G04 ──> N02, N04 ──> G4
    R02 ───> G02 ──> N03
    R02,R03,R05,R08 ──> G01 (problem sheet, G0) ──> conditional theorem G5
    R08 ───> G05 (first Lean extension)
    O01 ──> A01 ──> A02 (G1 with A03, A04)
