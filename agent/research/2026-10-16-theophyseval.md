# TheoPhysEval: lessons for our physics research engine

Read 2026-10-16. Supplied artifact: arXiv-2610.16733v1,
**TheoPhysEval: Assessing AI agents at autonomous theoretical physics problem
solving**. Main LaTeX and substantive included appendices/results were read,
including a delegated full-source review. No benchmark experiments were
reproduced, external dataset downloaded, or authenticity check performed.
The source, rather than new PDF OCR, was the reading backbone.

**Concrete outcome:** [Rust QFT kernel](../../Rust/qft-symbolics/README.md),
with exact symbolic residuals, negative controls, live Cadabra and Python
differential checks. This is research infrastructure and reproduction of
known results, not an autonomous discovery.

## Source provenance

Supplied paper directory:
`workspace2/Data/Public/papers/AI/arXiv-2610.16733v1/`.

SHA-256:
- PDF `2610.16733v1.pdf`:
  `720423de39bcab14b31fc3a746ba0215c6ae6feaf9dda8708584560ad55995cb`
- `paper.tex`:
  `c2addb9bf6e5c798609e2c37951e960fdbb8903db26a50868bfd5a91797a5b2f`
- `appendix.tex`:
  `afcbc6aeb5dd71c07c4ec15b58e6132f24c541216ce06928f38a834f52cd6922`
- `appsolutions.tex`:
  `4b13d7cd7d98496b52531746edf715789b4e210b849eb0403eb011bb0ffe4c81`

`00README.json` is compilation metadata, not authentication. The local
source lacks the external benchmark's graders, full solution bundle, task
sandboxes, auxiliary data and raw agent trajectories
(`appendix.tex:6`, `appsolutions.tex:1–4`). Reproducibility would require
obtaining a pinned external task/dataset revision, separately authorized.

## What the paper actually measures

The authors select 50 verifiable, mostly "turn-the-crank" research
calculations: adapt a known method, generalize a parameter/dimension, and
return a checkable answer. Significant open research problems were excluded.
They report 49/50 for their best realistic-track agent, and 46/50 for the
best isolated-track agents. These are the manuscript's **reported** results,
not results independently established here.

Important restrictions:

- One full run per model/harness/track, not repeated reliability estimates
  (`paper.tex:299`).
- Numerical **final answers** are graded, even when prompts ask for reasoning;
  symbolic answers are converted into numerical tasks
  (`paper.tex:271–281`). This does not validate every derivation step.
- No controlled head-to-head human study establishes the broader "surpassed
  experts" interpretation.
- Some tasks are linked: an upstream sign error propagated through an
  electronic-structure chain (`paper.tex:353`).
- The authors acknowledge AI-assisted problem preparation, possible leakage,
  non-adversarial selection, and uneven subject coverage.
- Small search/no-search changes do not by themselves prove memorization;
  isolated runs sometimes received local references and only one run was made.
- Choosing theories, important questions, observables and tractable
  formulations is distinct from solving a fully specified calculation
  (`paper.tex:68–72,385–393`). Observations about agents lacking research
  "taste" are explicitly anecdotal.
- The introduction asserts extraordinary external mathematical breakthroughs
  and names model versions. These claims were **not authenticated or
  independently verified**; do not repeat them as established facts.

## Most useful failure modes to encode as safeguards

1. **Reference formulas override the task.** Reported failures include copying
   reference initial conditions instead of the changed prompt conditions,
   using 5D formulas in 4D, and calculating geodesic lengths instead of
   minimal areas (`paper.tex:347–351`).
   **Action:** before deriving, write a reference-to-target delta: dimension,
   signature, normalization, initial/boundary data, observable, approximation
   order and branch. The task wins over retrieved examples.

2. **Numerical convergence is not correctness.** The atom prompt requires a
   hydrogenic known-spectrum check before self-consistency and warns about
   ill-conditioned logarithmic grids (`appendix.tex:1390–1423`).
   **Action:** validate analytic limits and residuals before production
   calculations, then perform precision/grid/parameter stability studies.

3. **The answer key can be wrong.** Authors describe correcting expert
   references after apparent model errors (`paper.tex:341`).
   **Action:** preserve model/reference disagreement; independently derive
   and version any corrected reference. Never treat supplied solutions or
   existing scripts as infallible.

4. **A correct scalar is insufficient.** Appendix checks include entire-cell
   band-node searches, Chern mesh and gap checks, SU(2) spin-sector agreement,
   pseudopotential norm/smoothness/nodelessness, and ghost-spectrum tests
   (`appendix.tex:958–998,1131–1197,1306–1357,1509–1539,1646–1676`).
   **Action for QFT:** exact bracket identities, on-shell/momentum constraints,
   free-index signatures, little-group weights, Ward identities,
   permutation/color relations, residues and soft/collinear limits.
   Keep cyclicity distinct from a Ward identity.

5. **Correlated mistakes survive multi-agent consensus.**
   The pass matrix offers no evidence that simply combining the listed agents
   would beat the best realistic run: every agent failed that run's sole
   unsolved case (`detailed_scores.tex:21–72`).
   **Action:** use heterogeneous verifiers, not a vote over explanations:
   component matrices vs Clifford rewriting; recurrence vs closed form;
   symbolic residual vs fresh numerical witnesses.

6. **Output/resource exhaustion is a result, not an invitation to loop forever.**
   Opus reportedly repeatedly hit an output limit; costs/tokens varied widely.
   **Action:** bound expansion size, runtime and attempts; checkpoint;
   record INCONCLUSIVE rather than silently weakening acceptance.

## A caution about the inverse-problem example

The only universally failed example asks for an autonomous differential
equation with low-order derivatives from a curve. Its proposed answer uses
Riemann–Liouville fractional derivatives and a delayed term, with a prescribed
history. The authors explicitly acknowledge wordplay and visual-fit grading
(`appendix.tex:21–29`, `appsolutions.tex:6–37`).
This is not a clean, uniquely posed ordinary-ODE reconstruction task.

**Lesson:** state the admissible model class, observations, fit criterion and
identifiability assumptions. Do not "solve" a problem by quietly relaxing
its mathematical meaning.

## Operational loop to adopt

**Human scientific objective → calculation contract → literature adaptation
audit → proposer → independent critic → executable verifier → result ledger
→ human direction/novelty review.**

A calculation contract should contain:

- Question/observable; why it matters; what changes relative to known work.
- Exact input/source revision and problem identifier.
- Signature, field normalization, units, index spaces/statistics, regulator,
  expansion order, branch, boundary/initial conditions.
- Allowed algebra and excluded singular domains.
- At least one known-positive and known-negative oracle test.
- Acceptance checks, tolerances with justification, resource/attempt budget.
- Required output artifacts and reproducibility metadata.

A result record should contain:

- Candidate and assumptions, cited starting equations and derivation.
- Backend/version, command, seed/fixture, git revision **plus source hashes
  if dirty**, branch and input dependencies.
- Symbolic residual, numerical residual/tolerance, failed and skipped checks.
- Status: reproduced / exact-within-domain / finite-witness-only /
  refuted / inconclusive.
- Whether the result is known; literature search status; next falsifying test.

Invalidate descendants when an upstream convention or result changes.
Do not let a downstream agent consume an artifact merely because a file exists.

## How this changes our present work

The Rust implementation now distinguishes:

- **Exact symbolic:** Schouten; slash products; polarization dot products;
  momentum-conservation factorization.
- **Exact finite representation:** sigma/gamma identities and all components
  of selected Lorentz/gamma relations.
- **Finite exact witnesses:** 200 SMGA formula/recurrence comparisons, n=3..10.
- **Numerical differential parity:** same fixtures through pre-existing NumPy,
  tolerance 1e-10.
- **Independent external symbolic backend:** actual Cadabra distribution,
  canonicalisation/collection and gamma joining; known-false sign rejected.

This closes some verification gaps; it does **not** establish an all-n
SMGA proof, general Cadabra replacement, full Python parity, or autonomous
scientific problem selection.

Next high-leverage step: typed tensor/index validation and ordinary Yang–Mills
vertices, followed by an independent 4-/5-point MHV recurrence.
The current half-collinear single-minus recurrence is a separate reproduction,
not a substitute for that Lorentzian/complexified Yang–Mills calculation.
