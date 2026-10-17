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
- **Symbolic verification:** Cadabra2 first (Docker image `cadabra2-ubuntu: 24.04`). SymPy / FORM as fallback.
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

- **2026-10-17:** New Frankel tooling on feat/lie-groups-so3-su2:
  Rust/lie-groups and CppMath/LieGroups implement numerical/symbolic SO3/so3,
  SU2/su2. Compile-time scalars/dimensions; validated private group state;
  f32/f64/software binary16 Rust; float/double/long double C++; exact polynomial,
  one-angle quotient and stereographic group charts. Normalize explicitly:
  T=-i sigma/2 maps to E; F=-i sigma maps to 2E. Opposite SU2 lifts remain distinct.
  Joint debug/release, C++ ASan/UBSan, live CAS and negative controls passed.
  Important review fixes: revalidate floating SO3 inverse; reset moved-from C++
  symbolic denominators/group values to valid states.
  Sage 10.8 generic SO(3,QQ) admitted a reflection; independently require det+1.
  Source audit/derivations and continuation are local to Rust/lie-groups.
  No original Propulsion code, book text, mono-notes or CAS source was changed.

- **2026-10-16:** Ernest directly authorized Rust QFT migration and TheoPhysEval
  lessons. Current implementation: Rust/qft-symbolics on feat/rust-qft-verifier.
  Two bounded slices: checked exact arithmetic, spinors/Clifford, strict typed
  tensor/index kernel, formal-d 47.2, evaluated 50.4, angular 50.2a, specialized
  SMGA recurrence.54 integration tests +2 compile-fail doctests; 13 textbook
  entries (partial chapters); 200 exact SMGA fixtures; live Cadabra differential
  checks. No arbitrary precision, ordinary Yang–Mills amplitudes, full BCFW
  recursion, held-out LLM eval, or background discovery loop.
  Canonical next steps: Rust/qft-symbolics/CONTINUE.md. This supersedes the July
  task order below for this directly requested tool-building project.

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
