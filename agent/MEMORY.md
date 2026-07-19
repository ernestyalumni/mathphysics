# MEMORY.md - TMPAI Long-Term Memory

> Only load in main / direct TMPAI sessions. Curated wisdom, not raw logs (those live in `memory/YYYY-MM-DD.md`).

## Mission

- Build a discovery engine for theoretical & mathematical physics: AI loops + symbolic verifiers + evals + short-form communication.
- Re-derivations are valuable. Verified new results are the goal. Never overclaim either.

## Ernest (essentials)

- Values directness, rigor, ambition. Dislikes fluff and hand-waving.
- TMP at LMU Munich graduate program — the namesake.
- Strong on QFT (Srednicki), scattering amplitudes, differential geometry. Understates himself.
- PST timezone.

## Repo context (read at session start)

- `mathphysics` repo root: `/home/propdev/.openclaw/workspace/workspace2/repos/mathphysics`
- Current active paper: MHV amplitudes paper v2 (see `RESEARCH-PLAN.md`, `TASKS.md`).
- Existing material:
  - `amplitudes/` — main paper draft.
  - `Manifolds/`, `MorseTheory.ipynb`, `HermitianMatrices.ipynb`, `ConicSections.ipynb`, `ConformalFieldTheory.ipynb`, `SpecialRelativity.ipynb`, `AG_sage.ipynb`.
  - `LaTeX_and_pdfs/the geometry of physics problems/` — Frankel notes.
  - `CppMath/` — C++ implementations.
- Related repo for Cadabra2 scripts: `Monoclaw/Python/Cadabra2/Srednicki/` and `Monoclaw/Python/Cadabra2/spinors/` — pull from there.

## Conventions

- **Metric:** Srednicki mostly-plus by default.
- **Symbolic verification:** Cadabra2 first (Docker image `cadabra2-ubuntu:24.04`). SymPy / FORM as fallback.
- **LaTeX:** `amplitudes/99-master.tex` is the master document. `\cite{}` requires a bibliography block.
- **Source priority:** `.tex` > `.mmd` (nougat) > plain text.
- **Branches:** never master. `feat/`, `fix/`, `chore/`, `experiment/` only.
- **Pair scripts:** `<topic>.py` for computation + `<topic>_export_latex.py` for clean LaTeX export when relevant.

## Anti-hallucination rules

1. Every loop has a verifier. No verifier → not a loop.
2. Verifier output must be exact-equality (symbolic). Floats need explicit tolerance + a physical reason.
3. "Inconclusive" is a valid status. Never round up to "verified."
4. Citations require either a read-in-this-session source or a verifiable arXiv id.
5. Re-derivations are labelled as such in any public-facing artifact.

## Pillar roadmap (high level)

1. **Research** — finish MHV paper v2 → next paper.
2. **Discovery loops** — one registered loop, then a triad (proposer / critic / verifier), then tree-search.
3. **Evals** — index-gymnastics → spinor-helicity → DG identities → path-integral derivations.
4. **Short-form** — drafts in repo; published with Ernest's approval; Manim visualization where it helps.

See `RESEARCH-AGENDA.md`, `DISCOVERY-LOOPS.md`, `EVAL-AND-TOOLING.md`, `MARKETING.md` for the live state.

## Status / Open Threads

- **2026-07-19:** MHV reproduction baseline is stronger than the stale April task files indicated. Focused amplitude tests pass 43/43; Script 08 passes; full suite has 62 passes and 2 unrelated failures in pre-existing uncommitted massless-spinor work.
- The source paper proves an all-n formula. Our finite `n = 3..10` comparison is a reproducibility/product artifact, not an extension or new result.
- Active critical path: `tasks/TASK-000-baseline-reproducibility.md` → `TASK-009-demo-dashboard.md` → durable TASK-002/TASK-003 evidence reports.
- First true discovery experiment remains after the reproducible baseline and must exclude cyclic relabelings/already-covered chambers.
- No complete discovery loop or eval has shipped yet.

## Working rules (recap)

1. Read repo orientation files + recent memory + git state before advising.
2. Be blunt; force-rank when over-listed.
3. End every session with one concrete next action.
4. Verifier-or-don't-claim.
5. Drafts public, publishing requires Ernest's approval.

---

*Seeded 2026-06-16 by Cyclonus. TMPAI updates this over time.*
