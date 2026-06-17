# RESEARCH-AGENDA.md — TMPAI Multi-Quarter Agenda

Last updated: 2026-06-16
Owner: Ernest Yeung
Curator: TMPAI

This is a living document. Edit when threads close or new ones open. Force-rank when there are more than 4 active threads.

## Ranking criteria

Use these in this order:

1. **Verifiability** — can the central claim of this thread be checked with Cadabra2 / SymPy / FORM, or by a clean limit comparison?
2. **Leverage** — does it unlock multiple downstream results (e.g. a recursion that yields many amplitudes)?
3. **Communicability** — can a verified result be turned into a 60–90s short with a clean visual beat?
4. **Existing footprint** — do we already have `.tex` / Cadabra2 scaffolding? Don't restart from zero if a 60% sketch already exists.
5. **Ernest's appetite** — TMP is a long-horizon project; threads Ernest is *actually excited about this quarter* beat strategically-correct-but-boring ones.

## Active threads (Q3 2026)

### 1. MHV amplitudes paper v2 (`feat/mhv-paper-v2` style work)

- **Status:** primary thread. See `../RESEARCH-PLAN.md` and `../TASKS.md`.
- **Verifier:** Cadabra2 scripts for Parke-Taylor, BCFW recursion relations.
- **Surface:** `amplitudes/Ch.05/bcfw.tex`, `Cadabra2/MHV/`, `99-master.tex`.
- **Reading list (priority):** Parke-Taylor 1986; Britto-Cachazo-Feng-Witten 2005; Elvang-Huang review; Dixon MHV review.
- **Next concrete actions:** see `../TASKS.md`.
- **Shorts potential:** "Why MHV is the simplest amplitude in QCD" — strong visual beat with helicity diagrams.

### 2. Cadabra2 verification scaffolding

- **Status:** infrastructure thread, blocks pillar 2 (discovery loops) until at least 3 verifiers are stable.
- **Verifier itself:** unit tests over canonical known-answer pairs.
- **Surface:** `Cadabra2/` (to be created in this repo; templates exist in `Monoclaw/Python/Cadabra2/`).
- **Next concrete actions:**
  - Port `mhv_parke_taylor.py` + `mhv_export_latex.py` template from Monoclaw to `Cadabra2/MHV/`.
  - Add a `Cadabra2/_verifier.py` wrapper that returns `{status: pass|fail|inconclusive, residual, notes}`.
  - Add a `Cadabra2/_test/` with 3 known-answer tests.
- **Shorts potential:** "Computer-verified physics in 90 seconds" — process video.

### 3. Differential geometry / Frankel program (`LaTeX_and_pdfs/the geometry of physics problems/`)

- **Status:** evergreen thread. Slow burn.
- **Verifier:** SymPy for explicit coordinate-chart computations; Cadabra2 for index gymnastics; manual for proofs.
- **Surface:** existing `.tex`, `Manifolds/`, `MorseTheory.ipynb`.
- **Next concrete actions:** decide which Frankel chapter to ship next (Lie derivatives? Curvature? Cartan structure equations?). Pick one and produce a derivation + verifier pair.
- **Shorts potential:** very high — "What is curvature really?" "Why differential forms?" — visual-first, narrative-driven.

### 4. Conformal field theory (`ConformalFieldTheory.ipynb`)

- **Status:** dormant. Re-evaluate Q4.
- **Verifier:** Cadabra2 for OPE manipulations; SymPy for explicit Virasoro algebra checks.
- **Next concrete actions (if reactivated):** scope a small ch. on 2d CFT primary operators + OPE; verifier checks.

## Backlog (not active this quarter)

- Path-integral derivations: Srednicki ch. 7–10 redo with Cadabra2 verifier where applicable.
- BCJ / double-copy: Bern-Carrasco-Johansson relations — heavy compute, big leverage if shipped.
- Holography / AdS-CFT primer: only after CFT thread is alive.
- C++ HPC for numerical amplitude eval: `CppMath/` exists; cross over with `Cosmos` if relevant.
- Special relativity refresh: `SpecialRelativity.ipynb` is a stub.

## Threads explicitly NOT in scope for TMPAI

- Career, interview, recruitment, networking → SpaceXAI persona.
- General SWE / DevOps / infrastructure → Cyclonus main session.
- Generic marketing (non-physics) → marketer persona.

## Quarterly review checklist

At end of each quarter, TMPAI should:

1. Re-rank active threads per the criteria above.
2. Move closed / shipped threads to a `closed-YYYYQN.md` archive.
3. Promote one backlog item if a slot opened.
4. Count: how many verified-discovery files in `discoveries/` this quarter? Trend?
5. Count: how many evals built? How many runs? Best / worst model?
6. Count: how many short-form drafts? How many published (after Ernest's approval)?

Ship the review to `MEMORY.md` as a 5-line summary.
