# TASK: Dixon–Srednicki Bridge Document

**Priority:** High  
**Status:** In Progress (sub-agent `tidy-sage` working on initial draft)  
**Output:** `amplitudes/09-dixon-srednicki-bridge.tex` + Cadabra2/Python scripts  
**Branch:** `refactor/amplitudes-structure`

## Goal

Create a detailed LaTeX document that bridges Dixon (2013) arXiv:1310.5353 with Srednicki QFT, with **full step-by-step derivations** (no hand-waving) and **Cadabra2/Python scripts** for every computation.

## Source Materials

| Resource | Path |
|----------|------|
| Dixon paper (LaTeX) | `workspace2/Data/Public/papers/physics/arXiv-1310.5353v1/ModAmpIntro.tex` |
| Srednicki Ch.50 | `amplitudes/Srednicki/individual_chapters/ch50_*.tex` |
| Srednicki Ch.60 | `amplitudes/Srednicki/individual_chapters/ch60_*.tex` |
| Srednicki Ch.34-36 | `amplitudes/Srednicki/individual_chapters/ch34_*.tex` etc. |
| Existing preamble | `amplitudes/00-preamble.tex` |
| SMGA paper | `workspace2/Data/Public/papers/physics/arXiv-2602.12176v2/SMGA.tex` |

## Sections Required

### Section 1: Color Decomposition (Dixon §2 ↔ Srednicki Ch.79-80)
- SU(N) generators with Tr(TᵃTᵇ) = δᵃᵇ (Dixon normalization)
- **Derive**: [Tᵃ, Tᵇ] = i√2 f^{abc} Tᶜ — explain the √2 comes from doubled trace norm vs standard Tr = ½δ
- **Derive**: Dixon (2.2): i√2 f^{abc} = Tr(TᵃTᵇTᶜ) - Tr(TᵃTᶜTᵇ). Show ALL steps from commutator.
- **Derive**: Srednicki 80.17 — trace back to what leads to it
- SU(Nₑ) Fierz identity Dixon (2.3) — full derivation
- **Prove**: tracelessness statement. Contract (2.3) with δ^{ī₁}_{i₁}, show every index contraction explicitly
- Compare Srednicki 81.2 with Dixon (2.4) — are they the same? (Yes, single-trace tree decomposition)
- Double-trace terms: only at loop level (Dixon discusses, Srednicki doesn't)
- **Cadabra2 script**: Verify Fierz identity numerically for SU(2), SU(3)

### Section 2: Spinor Variables (Dixon §3.1 ↔ Srednicki Ch.48, 50)
- **Notation comparison table**: Dixon vs Srednicki conventions side by side
  - Dixon: |i⁺⟩ ≡ λᵢᵅ (right-handed), |i⁻⟩ ≡ λ̃ᵢ^{α̇} (left-handed)
  - Srednicki: different labeling — document precisely
- Dixon (3.1-3.2) spinor products ↔ Srednicki angle/square brackets
- **Derive**: Dixon (3.5) positive energy projector in 2-component form ↔ Srednicki equivalent
- **Show explicitly**: Dixon (3.6) ≡ Srednicki 50.15 — write both, match term by term
- Complex square-root property (Dixon after 3.12) — explain why important for collinear singularities and connection to SMGA half-collinear regime
- **Cadabra2 script**: Construct explicit spinors for given momenta, verify p² = 0

### Section 3: Spinor Product Identities (Dixon §3 ↔ Srednicki Ch.50)
- **Derive each** with full steps:
  - Anti-symmetry: Dixon (3.13) ↔ Srednicki
  - Squaring: Dixon (3.14) ↔ Srednicki  
  - Momentum conservation: Dixon (3.15) ↔ Srednicki
  - Schouten identity: Dixon (3.16) ↔ Srednicki Eq. 50.36 / Problem 50.3
- For Schouten: Ernest has solved Problem 50.3 — check his solutions file at `workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/Srednicki_QFT.tex`
- **Cadabra2 script**: Verify all identities symbolically. If Cadabra2 can't handle it, use SymPy or another package and document why.
- **Python script**: Numerical verification with random momenta

### Section 4: Four-Point Example (Dixon §3.2 ↔ Srednicki Ch.60)
- Dixon: e⁺e⁻ → massless fermion pair. Srednicki Ch.60: spinor QED with photons.
- Explain connection and differences
- **Derive**: Dixon Eq. 3.24 using Srednicki notation and identities. Every step.
- **Derive**: Dixon Eq. 3.25 using Srednicki notation and identities. Every step.
- **Cadabra2/Python script**: Compute the amplitude numerically for sample kinematics

### Section 5: Helicity Formalism for Massless Vectors (Dixon §3.3 ↔ Srednicki Ch.81)
- How Srednicki Ch.81 relates to Dixon §3.3
- Important results with full derivations and rationale
- Polarization vectors in spinor-helicity form
- **Cadabra2 script**: Verify polarization vector properties (transversality, completeness)

### Section 6: Soft Gluon Limit (Dixon §4.1)
- k₄ → 0 means all components of p⁴_μ → 0 (soft limit)
- **Derive**: Dixon Eq. 4.1 step by step in Srednicki notation. Show every algebraic step.
- Connect to SMGA paper: their Weinberg soft theorem (Eq. 19) is the same physics
- **Python script**: Numerical check of soft limit for 4-point and 5-point amplitudes

## Style Requirements

- **Be explicit and thorough** — no "it can be shown that" or "one can verify"
- Show EVERY algebraic step in derivations
- Write expressions in BOTH Dixon and Srednicki notation
- Every equation that can be verified computationally should have a corresponding script
- Scripts should be standalone Python files (using cadabra2 or sympy) that can be run directly
- Save scripts in `amplitudes/scripts/` directory (create if needed)
- Each script should print clear output showing what it verifies

## After Completion

1. Update `amplitudes/99-master.tex` to include the new chapter
2. Test compile: `cd amplitudes && pdflatex 99-master.tex`
3. Commit all new files on `refactor/amplitudes-structure` branch

## For the Agent

- Read the actual source `.tex` files, don't guess at equations
- Cross-reference equation numbers carefully (e.g. "Srednicki 50.15" means eq 50.15 in the textbook)
- If you can't find an equation in Srednicki, say so explicitly
- The nougat parse of Srednicki is at: `workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/[Mark_Srednicki]_Quantum_Field_Theory(BookSee.org).mmd`

---

Created: 2026-04-05
