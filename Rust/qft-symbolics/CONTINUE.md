# Continue the Rust QFT corpus

Branch: feat/rust-qft-verifier. Read root AGENTS.md, this crate's README.md,
corpus.json, docs/{MIGRATION,DERIVATIONS,TENSOR,TEXTBOOK-EXTENDED,VERIFICATION}.md,
and ../../agent/research/2026-10-16-theophyseval.md.

## Current stopping boundary

Two bounded implementation slices are complete, not the user's entire project.
Current suite: 54 integration tests plus two compile-fail doctests; 13 named
textbook entries (partial chapter coverage); 200 exact SMGA fixtures.
No autonomous background job was installed or left running.

Second slice implements:
- Strict Lorentz/dotted/undotted indices and variance; validated free signatures
  and typed zeros; hygienic local dummy binders and multiplication.
- Independent bounded exhaustive commuting-factor/signed-symmetry orbits,
  including Lorentz dummy variance orientation and sign-stabilizer zeros.
- Constant symmetric Lorentz metric and typed delta contractions; formal d.
- Twelve positive/negative tensor residuals agreeing with live Cadabra.
- Formal-d47.2 via a derived explicit Clifford sandwich operation, independently
  checked by Cadabra and 498 finite matrix component assignments.
- Fully evaluated 50.4 chiral trace, scalar invariants and epsilon orientation.
- Validated rational angular 50.2a constructor, general derivation and 48 exact
  witnesses; invalid-angle and direct spinor-phase controls.

Read-only review corrections are recorded in docs/VERIFICATION.md.

## Next implementation units

1. **Coefficient/resource layer:** arbitrary-precision rational coefficients
   and uniformly structured arithmetic errors. No dependency downloads or
   system installation have been approved/performed. Current checked i128
   coefficients and legacy Clifford operations can panic on exhaustion;
   tensor validation/orbit errors already return Result.
2. **Remaining selected textbook base:**50.2b's explicit boosted-spinor
   massless limit; 60.2b–60.4 with explicit conservation. Preserve problem IDs,
   statements and hints; consult source PDF for malformed transcriptions.
3. **Nonabelian/QCD branch:**69/70/71/72 then 79/80/81, following the corpus
   prerequisite subgraph. Derive color traces/Jacobi and 3-/4-gluon vertices.
   Extend tensor facilities only as required; current factors all commute.
4. **Independent amplitudes:** ordinary complexified/Lorentzian Yang–Mills BG
   or explicit diagrams at4/5 points against Parke–Taylor, including Ward
   identities and all normalization conventions; only then BCFW recursion.
   Existing BCFW is a kinematic shift; SMGA BG is a different recurrence.
5. **Tensor followups:** scalable signed group algorithms; explicit spinor
   epsilon raising/lowering; multi-term identities; noncommuting statistics
   and derivative/dependency semantics. Do not imply these already exist.
6. Only then propose bounded generalizations. Human owns scientific relevance;
   fresh verification, literature audit and held-out tests precede novelty.

## Guardrails

- Reinspect PDD ownership when introducing a managed project.
- Never stage existing mono-notes PDFs/build outputs.
- Keep all build products, generated logs/JSON and Python caches outside git.
- Preserve negative controls and all designated SMGA fixtures; domain/budget
  errors are not skipped-success points. A failed process is not verified.
- Clifford label strings still do not imply tensor summation. Formal-d gamma
  sandwich is an explicit operation, not generic gamma/tensor contraction.
- MIT crate algebra is independently authored. Translating GPL Cadabra C++
  would require a separate licensing decision, not merely Rust syntax.
