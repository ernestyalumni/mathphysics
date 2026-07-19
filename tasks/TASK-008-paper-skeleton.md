# TASK-008 — arXiv paper skeleton wired to computed results

Status: TODO
Priority: P1
Depends on: TASK-002, TASK-003 (results to cite)
Estimated size: one session (3–4 h)
Branch suggestion: `feat/task-008-paper-skeleton`

## Goal

Turn `amplitudes/99-master.tex` + chapter files into an arXiv-shaped paper skeleton where
every numerical/symbolic claim is pulled from, or explicitly cites, a script and result
file in this repo. Target: a submittable "methods + reproduction + extensions" paper even
before any new discovery lands; discovery results (TASK-006/007) slot in when ready.

## Structure (target)

1. Introduction — the verified-discovery-loop thesis (crib from `../TMP-MISSION.md`,
   tone: sober, no AI hype).
2. Conventions & spinor-helicity (condense existing `02-helicity.tex`, `03-spinor-helicity.tex`).
3. MHV & Parke-Taylor with verification (existing `04-mhv.tex` + TASK-001 test results).
4. Single-minus amplitudes in half-collinear regimes: SMGA reproduction (TASK-002) and
   extended verification of Eq.(16) to n = 10 (TASK-003).
5. Consistency-check battery as a methodology section (TASK-004).
6. Open directions / new results (TASK-006, TASK-007 as available).
7. Reproducibility appendix: exact commands, seeds, environment.

## Steps

1. Reorganize `99-master.tex` includes to the structure above; don't delete old chapter
   files — re-include or mark unused.
2. Create `amplitudes/results/latex/` where scripts export LaTeX tables (add an export
   flag to the TASK-002/003 scripts or a small converter from their .md/.json output).
3. Compile cleanly: `latexmk -pdf 99-master.tex` with zero errors (warnings logged).
4. Add a `make paper` or `scripts/build_paper.sh` one-shot build.

## Definition of done

PDF builds from a clean checkout of the branch with one documented command; every table in
sections 3–5 is generated, not hand-typed; a TODO list of missing pieces sits at the top of
the introduction as a LaTeX comment.
