# Amplitudes Study Guide — What to Read and Why

This guide maps out the minimum reading needed to make discoveries in MHV scattering amplitudes at the level of arXiv:2602.12176v2 (SMGA / OpenAI).

## Srednicki QFT — Required Chapters

### Tier 1: CRITICAL (read carefully, do computations)

| Ch. | Title | Why |
|-----|-------|-----|
| **34** | Left- and Right-Handed Spinor Fields | Foundation: Weyl spinors, ε-tensor, σ-matrices |
| **35** | Manipulating Spinor Indices | σ-algebra, index gymnastics you'll use constantly |
| **36** | Lagrangians for Spinor Fields | Weyl/Majorana/Dirac connection |
| **38** | Spinor Technology for Spin-1/2 | Spin sums, trace technology |
| **48** | Spin-1/2 Particles and Spinor-Helicity | **THE chapter** — λ, λ̃, angle/square brackets, polarization vectors |
| **50** | Massless Particles and Spinor Helicity | Massless limit, completeness relation -k̸ = |k]⟨k| + |k⟩[k| |
| **60** | Spinor Helicity for Spinor QED | MHV amplitudes, Parke-Taylor formula |

### Tier 2: IMPORTANT (read, understand key results)

| Ch. | Title | Why |
|-----|-------|-----|
| **10** | Scattering Amplitudes and Feynman Rules | Basic S-matrix, how amplitudes work |
| **27** | Other Renormalization Schemes | MS-bar, needed for loop-level later |
| **37** | Canonical Quantization of Spinor Fields | Mode expansion, anticommutation |
| **69** | Nonabelian Gauge Theory | SU(N) gauge fields, gluon vertices |
| **70** | Group Representations | Color algebra, Tr[TᵃTᵇ] = ½δᵃᵇ |
| **72** | Feynman Rules for Nonabelian Gauge Theory | Gluon propagator, 3/4-gluon vertices |
| **80** | The Feynman Rules for N × N Matrix Fields | Color-ordered amplitudes |

### Tier 3: SKIM (useful background, read as needed)

| Ch. | Title | Why |
|-----|-------|-----|
| **1–5** | Attempts at Relativistic QM, Lorentz invariance | Basics (you probably know this) |
| **39–43** | More spinor tech, Coulomb gauge, Majorana rules | Deepens Tier 1; Ch.43 connects to SUSY |
| **81–84** | Loop amplitudes in gauge theory | Needed if we go to loop level |

### What You Can Skip

Most of Part I scalar field theory (Ch.3–9, 12–26) and Part II renormalization details unless you need them for loop calculations later.

---

## Papers — Required Reading

### Must Read (before attempting discovery)

| # | Paper | arXiv | Why |
|---|-------|-------|-----|
| 1 | Dixon, "A brief introduction to modern amplitude methods" (2013) | [1310.5353](https://arxiv.org/abs/1310.5353) | Pedagogical review: spinor-helicity, color decomposition, BCFW, unitarity |
| 2 | Guevara et al., "Single-minus gluon tree amplitudes are nonzero" (2025) | [2602.12176](https://arxiv.org/abs/2602.12176) | **Our benchmark** — the discovery we want to match/exceed |
| 3 | Elvang & Huang, "Scattering Amplitudes in Gauge Theory and Gravity" (2013) | [1308.1697](https://arxiv.org/abs/1308.1697) | THE modern amplitudes textbook. Covers everything Dixon does but deeper, plus BCFW, superamplitudes |
| 4 | Witten, "Perturbative Gauge Theory As A String Theory In Twistor Space" (2003) | [hep-th/0312171](https://arxiv.org/abs/hep-th/0312171) | Started the modern amplitudes revolution. SMGA explicitly references it |

### Should Read (to find new directions)

| # | Paper | arXiv | Why |
|---|-------|-------|-----|
| 5 | Britto, Cachazo, Feng, Witten, "Direct Proof Of Tree-Level Recursion Relation..." (2005) | [hep-th/0501052](https://arxiv.org/abs/hep-th/0501052) | The original BCFW paper. Short, elegant |
| 6 | Parke & Taylor, "Amplitude for n-Gluon Scattering" (1986) | Phys. Rev. Lett. 56, 2459 | The original MHV formula |
| 7 | Berends & Giele, "Recursive Calculations for Processes with n Gluons" (1988) | Nucl. Phys. B306, 759 | The recursion SMGA uses as their proof backbone |
| 8 | Arkani-Hamed et al., "Scattering Amplitudes and the Positive Grassmannian" (2012) | [1212.5605](https://arxiv.org/abs/1212.5605) | Amplituhedron, geometric approach |

### Nice to Have (deeper exploration)

| # | Paper | arXiv | Why |
|---|-------|-------|-----|
| 9 | Bern, Dixon, Kosower, "One-Loop Amplitudes..." review (2011) | [1103.1559](https://arxiv.org/abs/1103.1559) | One-loop methods |
| 10 | Strominger, "Lectures on the Infrared Structure of Gravity and Gauge Theory" (2017) | [1703.05448](https://arxiv.org/abs/1703.05448) | Celestial holography (SMGA mentions as extension) |

---

## Recommended Reading Order

1. **Finish Dixon** (1310.5353) — you're in it now
2. **Srednicki Ch.48 + 50 + 60** carefully — the spinor-helicity core
3. **Read SMGA** (2602.12176) thoroughly — understand every equation
4. **Get Elvang & Huang** (1308.1697) — replaces needing many individual papers
5. **Skim Witten 2003 and BCFW 2005** for the key ideas
6. Then start reproducing SMGA results computationally

---

## Local File Locations

| Resource | Path |
|----------|------|
| Srednicki textbook (LaTeX) | `workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/Srednicki_QFT.tex` |
| Srednicki textbook (nougat) | `workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/[Mark_Srednicki]_Quantum_Field_Theory(BookSee.org).mmd` |
| SMGA paper (LaTeX) | `workspace2/Data/Public/papers/physics/arXiv-2602.12176v2/SMGA.tex` |
| Dixon review (LaTeX) | `workspace2/Data/Public/papers/physics/arXiv-1310.5353v1/ModAmpIntro.tex` |
| Our notes | `amplitudes/08-literature-review.tex` |
| Discovery plan | `amplitudes/DISCOVERY-PLAN.md` |
| Reading list | `amplitudes/READING-LIST.md` |
| Cadabra2 computations | `amplitudes/Srednicki/individual_chapters/` |

---

Last updated: 2026-04-05
