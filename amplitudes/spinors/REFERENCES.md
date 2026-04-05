# MHV Amplitudes and Spinor-Helicity Formalism — References

Key papers and resources for the physics covered in this directory.
Ordered from most accessible to most advanced.

---

## Most Accessible for Physics Grad Students (Srednicki background)

### [1] Mangano & Parke (1991) — Classic Review
**Title:** Multiparton processes in gauge theories
**Authors:** M.L. Mangano, S.J. Parke
**Journal:** Phys. Rept. 200, 301 (1991)
**arXiv:** hep-th/0509223 (later preprint version available)

**Summary:** The definitive pedagogical review of spinor-helicity methods,
color ordering of amplitudes, and MHV amplitudes in gauge theories. Covers
the Parke-Taylor formula, Berends-Giele recursion, and explicit examples.
**Recommended starting point** for anyone with Srednicki QFT background.

---

### [2] Parke & Taylor (1986) — Original MHV Paper
**Title:** Amplitude for n-Gluon Scattering
**Authors:** S.J. Parke, T.R. Taylor
**Journal:** Phys. Rev. Lett. 56, 2459 (1986)
**DOI:** 10.1103/PhysRevLett.56.2459
**INSPIRE:** https://inspirehep.net/literature/227338

**Summary:** The original two-page letter presenting the Parke-Taylor formula
for n-gluon MHV amplitudes as a conjecture (later proved by Berends-Giele).
The formula A_n = ⟨ij⟩^4 / (⟨12⟩⟨23⟩...⟨n1⟩) compresses what would be
thousands of Feynman diagrams into a single expression. Historically landmark.

---

### [3] Elvang & Huang (2015) — Modern Textbook
**Title:** Scattering Amplitudes in Gauge Theory and Gravity
**Authors:** Henriette Elvang, Yu-tin Huang
**Publisher:** Cambridge University Press (2015)
**arXiv:** 1308.1697 [hep-th] (2013 preprint)
**Cambridge:** https://www.cambridge.org/core/books/scattering-amplitudes-in-gauge-theory-and-gravity/34C045B8331E6FF229E2496F6D8321C5

**Summary:** The best modern pedagogical textbook on amplitudes. Starts from
spinor-helicity formalism (Chapter 2), covers MHV amplitudes and Parke-Taylor
(Chapter 2-3), BCFW recursion (Chapter 3), superamplitudes (Chapter 4), and
modern topics (Grassmannians, amplituhedron). Freely available on arXiv.
**Best book** for a systematic modern treatment after Srednicki.

---

## BCFW Recursion

### [4] Britto, Cachazo, Feng (2005) — Original Recursion
**Title:** New Recursion Relations for Tree Amplitudes of Gluons
**Authors:** R. Britto, F. Cachazo, B. Feng
**Journal:** Nucl. Phys. B 715, 499 (2005)
**arXiv:** hep-th/0412308

**Summary:** Discovered on-shell recursion relations for tree-level gluon
amplitudes by deforming two external momenta by a complex parameter z and
using complex analysis (residue theorem). The [j,i⟩ shift: p̂_i → p̂_i + z q,
p̂_j → p̂_j − z q factorizes amplitudes on physical poles.

---

### [5] Britto, Cachazo, Feng, Witten (2005) — Proof of BCFW
**Title:** Direct Proof of Tree-Level Recursion Relation in Yang-Mills Theory
**Authors:** R. Britto, F. Cachazo, B. Feng, E. Witten
**Journal:** Phys. Rev. Lett. 94, 181602 (2005)
**arXiv:** hep-th/0501052

**Summary:** Witten joined and provided the proof that the BCFW recursion
follows from basic properties of tree-level Yang-Mills amplitudes (their
rational structure and behavior at large complex momentum). With Parke-Taylor
as seed (3-point amplitudes), BCFW generates all tree-level gauge amplitudes
without Feynman diagrams. One of the most impactful papers in modern amplitudes.

---

## Twistor Theory Connection

### [6] Witten (2004) — Twistor String Theory
**Title:** Perturbative Gauge Theory As A String Theory In Twistor Space
**Authors:** E. Witten
**Journal:** Commun. Math. Phys. 252, 189 (2004)
**arXiv:** hep-th/0312171

**Summary:** Showed that n-gluon MHV amplitudes in Yang-Mills theory, when
Fourier transformed to twistor space (Penrose 1967), are supported on degree-1
algebraic curves (lines). Higher degree curves correspond to non-MHV (NMHV etc.)
amplitudes. Proposed that N=4 SYM is equivalent to a topological B-model string
on CP^{3|4}. Initiated the modern era of amplitude methods and the CSW rules.

---

### [7] Cachazo, Svrcek, Witten (2004) — CSW/MHV Vertices
**Title:** MHV Vertices And Tree Amplitudes In Gauge Theory
**Authors:** F. Cachazo, P. Svrcek, E. Witten
**arXiv:** hep-th/0403047

**Summary:** Constructed Feynman-like rules where MHV amplitudes themselves
serve as vertices (off-shell continued), reproducing all tree-level gauge
theory amplitudes. The "CSW rules" or "MHV vertex expansion" is the direct
application of Witten's twistor string insight to practical computation.

---

## Two-Component Spinor Techniques

### [8] Dreiner, Haber, Martin (2010) — 2-Component Spinor Bible
**Title:** Two-component spinor techniques and Feynman rules for quantum field theory and supersymmetry
**Authors:** H.K. Dreiner, H.E. Haber, S.P. Martin
**Journal:** Phys. Rept. 494, 1 (2010)
**arXiv:** 0812.1594 [hep-ph]

**Summary:** Comprehensive 230-page reference for 2-component (Weyl) spinor
calculations. Covers van der Waerden notation, Fierz identities, Feynman rules
in 2-component notation, SUSY, and explicit translations between 2- and
4-component formalisms. The definitive reference for Fierz identities and
sigma matrix algebra beyond what Srednicki covers.

---

## Primary Physics Text Reference

### [9] Srednicki — Quantum Field Theory (Chapters 34-38)
**Title:** Quantum Field Theory
**Author:** M. Srednicki
**Publisher:** Cambridge University Press (2007)
**Free draft:** https://web.physics.ucsb.edu/~mark/qft.html

**Summary:**
- **Ch. 34:** Two-component spinors, van der Waerden notation, ε tensors, SL(2,C) reps
- **Ch. 35:** Dirac equation from 2-component spinors; σ^μ and σ̄^μ matrices
- **Ch. 36:** σ^μ completeness, momentum spinors, helicity
- **Ch. 37:** Spinor-helicity formalism, polarization vectors as spinor bilinears
- **Ch. 38:** Helicity Feynman rules for Yang-Mills, derivation of 4-gluon MHV amplitude
- **App. B:** Sigma matrix identities, Fierz rearrangement
