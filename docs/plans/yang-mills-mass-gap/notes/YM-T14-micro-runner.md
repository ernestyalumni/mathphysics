# YM-T14 — micro-task runner for local models, and the first unattended ten-hour run

Author: strategist (Claude Opus 5.5 at the time of the launch; this session started on Fable 5.2), 2026-11-17 UTC. Branch `ym/YM-T14-claude-fable52-q7nd`. Requested by Ernest: give the local models ten hours of Yang–Mills work (research, `.tex` writing), on the desktop and on the Mac Mini.

## 1. Result

- `bin/ym-micro` (standard library only) runs a queue of micro-tasks against one local OpenAI-compatible endpoint, resumes after interruption (one result file per task and pass), waits out a dead endpoint, and stops on a `STOP` file in its output directory. `bin/ym_micro_plan.py` builds the units, briefs and queues. `tests/test_ym_micro.py`: 14 offline tests, all pass.
- Task types, each with a check a script can run (protocol §12):

| Type | What the model does | What the runner checks | What counts |
|---|---|---|---|
| `referee` | Judges a proof split into numbered steps, several seeds, the probe's `REVIEW_QUESTION` | Verdict parsed | INCORRECT flags only |
| (calibration) | Same, on *mutants*: one token changed at a known step (`\le`↔`\ge`, an exponent) | Was it flagged; was the right step named | Detection rate per model |
| `quote` | Answers a reading-packet question with ITEM/QUOTE pairs from LaTeX source | The quote occurs verbatim in the source, whitespace ignored; the line number is computed by the runner | Only verbatim quotes |
| `tex-draft` | Writes one subsection from a brief, under the manuscript's style rules | No forbidden commands; balanced environments; compiles with `pdflatex -no-shell-escape`, paranoid file access, in a temporary directory; hidden oracles (expected formulas) are matched by regex, never shown to the model | Raw material for a human or T3 rewrite |
| `draft-referee` | Referees the drafts produced by either machine | as `referee` | flags |
| `lean-reading`, `lean-compare` | Reads one Lean declaration (the packet's original task); compares it with Armstrong's one-line summary | Identifiers occur in the snippet; verdict parsed | Cards; INCONSISTENT flags |
| `pi` | Headless Pi with read-only tools builds a source-locator table (YM-R09) | Every cited `path:line` exists | Locator table |
| `probe` | The admission test (`probe.py`) | keyword grading | Admission evidence |

- Admission test of the Mac model (first run, 2026-11-16, one run per row): `qwen38-27b-turbo-iq2m` (a community 27B Qwen3.8 merge at 2.7 bits) passed **rows b, c, d, e, f, g, h, i, j, k**; row a was cut off at the 4,000-token limit. The 9B models fail f, g, h and j (`notes/local-models.md` §4a). By the rule of §6 this makes it a *candidate* for T2 review work. One run is not a measurement; two more admission runs are in its queue. Speed: about 8.7 generated tokens/s and 60–80 prompt tokens/s, i.e. five times slower than the 9B on the RTX 3060.

## 2. Labels

All model output of the run: `unreviewed (local model output; not evidence)`. The runner, the checks and the tests: done by the author, not reviewed by another family.

## 3. What is queued (ten hours each, priority order; the tail will not be reached)

**Desktop, Qwythos 9B Q5 (`127.0.0.1:8080`, Ernest's `launch.sh`):** referee of the ten proofs of the mono-notes topic plus the transfer-positivity sketch C23 (three seeds, repeatable passes) → 14 LaTeX briefs → 24 proof mutants → YM-R01 and YM-R02 quote sheets (16 questions) → the ten subsections of the topic, refereed in windows of ten steps → drafts refereed → 15 proofs of N1 and N2 §2 as a control (source text, presumed correct) → Lean cards and summary comparison for the 150 main results → two Pi locator maps for YM-R09 → section mutants → YM-R03, R04 quote sheets → the remaining 839 Lean rows.

**Mac, Qwen3.8-27B IQ2_M (Mac loopback `:8080`, reached through an SSH forward to desktop `:8082`):** the same ten proofs and C23 (two seeds) → admission run 2 → briefs 1–7 → eight mutants → the ten subsections → briefs 8–14 → N2 §2 and N1 Prop. 2.2 → YM-R02 quote sheet → six questions over **the whole of N2** in one context (about 100k tokens; only the Mac has the context) → YM-R04 quote sheet → admission run 3 → referee of the desktop's drafts.

Briefs (`bin/ym_micro_plan.py`, `BRIEFS`): compact self-adjoint integral operators and Mercer; Perron–Frobenius and Jentzsch; Peter–Weyl and convolution by a class function; SU(2) character coefficients of the Wilson weight (Bessel functions); unbounded spectral theorem and the definition of the gap; reflection positivity of the Wilson measure; SU(2) Haar moments (the YM-T09 values); three meanings of matrix integration (YM-T09); the Laplacian on a compact group and the heat-kernel weight; masses from correlators; exact 2D U(1); exercises on the Wilson measure; exact 2D Z₂ on a torus; the 2D SU(2) Wilson loop. Each has two to four hidden oracle patterns.

## 4. Evidence

- Smoke runs of every task type, 2026-11-17 (outside git, `~/.local/state/mathphysics/ym-overnight/2026-11-16/smoke*`): referee 171 s for three seeds of Prop. `ym:prop:purity` (no flags); a mutant of Prop. `ym:prop:defect` was flagged once in three runs, at the wrong step; a quote sheet for YM-R01 returned three items, all verbatim; a Pi locator table returned five Srednicki locators, four of which the strategist opened and found correct (§§ "The beta function in nonabelian gauge theory", "Wilson loops, lattice theory, and confinement", "Spontaneous breaking of gauge symmetries", and the asymptotic-freedom paragraph); a Z₂ draft derived Z = cosh^N β + sinh^N β correctly but contains a false sentence ("deg_S(ℓ) … can never be 1") and failed to compile (`Z_2` in a section title outside math mode).
- Run state, logs, the private launcher (it names the Mac) and results: `~/.local/state/mathphysics/ym-overnight/2026-11-16/` on the desktop. `REPORT.md` there is regenerated every 30 minutes by `ym-micro report`.

## 5. Not done / not checked

- The packet's 20-row `lean-reading` report and a reviewer's usefulness judgment: the Lean rows run late in the desktop queue; judge them from `REPORT.md` after the run.
- Splitting LaTeX into steps is heuristic (sentence boundaries, displays kept inside their sentence); long subsections are refereed in windows of ten steps, so an error that needs two windows is invisible.
- Mutants are presumed false, not proved false (a flipped `\le` inside a hypothesis can leave a true statement).
- The Mac server is single-slot and is shared with Ernest's Pi sessions there; so is the desktop server with his interactive Pi. Both are slower while the runs go.
- No review of the runner by another model family.

## Resume here

1. Read `~/.local/state/mathphysics/ym-overnight/2026-11-16/REPORT.md` (or run `bin/ym-micro report --state <that dir>`).
2. Calibration first: if a model's detection rate on mutants is low, its silence on our proofs means nothing (it means nothing anyway) and its flags need care.
3. Check every flag on our proofs and on C23 by hand; post real findings on local:13 (G02). Flags on N1/N2 (presumed-correct source) measure false alarms.
4. Drafts: read those whose oracles hit; none go into the manuscript without a T3 rewrite and review.
5. Admission runs 2 and 3 of the Mac model: record in `notes/local-models.md` §4a; if f, g, h, j hold, give it T2 referee packets.
