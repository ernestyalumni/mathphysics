# Continue the Rust QFT corpus

Branch: `feat/rust-qft-verifier`. Start by reading README.md, corpus.json,
docs/MIGRATION.md, docs/DERIVATIONS.md, docs/VERIFICATION.md and the paper
lessons at ../../agent/research/2026-10-16-theophyseval.md.

## Current stopping boundary

The first bounded implementation slice is complete and tested, not the
user's whole long-term request. It provides actual exact algebra, selected
Srednicki checks and an exact port of the existing specialized SMGA routines.
No background autonomous job was installed or left running.

## Next implementation units, in order

1. **Typed-index kernel:** separate Lorentz/dotted/undotted spaces, variance,
   validated free/dummy signatures, constant metric/delta contractions.
   Reject triple occurrences and free-index mismatch. Start with bounded,
   independently implemented symmetry-orbit enumeration; do not pretend
   alphabetical names implement Cadabra/xPerm.
   Positive/negative Cadabra differential checks before extending the corpus.
2. **Coefficient/resource layer:** arbitrary-precision rationals and structured
   budget errors. Obtain any needed new dependencies under local policy.
   Current checked i128 arithmetic is deliberately finite-range.
3. **Finish the selected textbook base:** formal d in 47.2 only after index
   semantics; 50.2, evaluated 50.4; 60.2b–60.4 with explicit conservation.
   Preserve original problem IDs and statement/hint text where legible.
   Consult source PDF for malformed transcriptions before adopting formulas.
4. **Nonabelian/QCD branch:** 69/70/71/72, then 79/80/81, following the
   prerequisite subgraph. Color traces and antisymmetric structure constants,
   derive and verify 3-/4-gluon vertices.
5. **Independent amplitudes:** ordinary complexified/Lorentzian Yang–Mills
   BG or explicit diagrams at four/five points versus the Parke–Taylor ratio;
   then full BCFW amplitude recursion. The existing BCFW function is only a
   kinematic shift, and the SMGA BG port is a different recurrence.
6. Only then propose bounded generalizations. Human still owns scientific
   relevance; exact verifier, fresh-point tests and literature audit precede
   any originality claim.

## Guardrails

- Read root AGENTS.md; inspect PDD ownership again if adding a managed project.
- Never stage existing mono-notes PDF/build-output changes.
- All build/output logs go to external cache; never commit target/, Python
  caches, LaTeX intermediates or generated per-run JSON.
- Keep negative controls, poles, invalid inputs and resource exhaustion
  first-class. Library process failure is not verified.
- Do not silently turn wall/domain errors into skipped-success examples.
- MIT crate code was independently authored. A future translation of Cadabra's
  GPL C++ requires a separate licensing decision, not merely Rust syntax.
