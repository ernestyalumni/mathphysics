# Research Plan: MHV Amplitudes & Modern Methods in Scattering Theory

**Repository:** `mathphysics`
**Goal:** Produce an arXiv-quality review + original work paper on MHV amplitudes, spinor-helicity formalism, and BCFW recursion, with heavy symbolic computation using Cadabra2.

We are building on the style of the recent MHV paper that leveraged LLMs (ChatGPT), but we will make ours more rigorous by grounding everything in Cadabra2 symbolic computations and Srednicki.

## Core Philosophy
- **Cadabra2 First:** Every major result in the paper should have a corresponding Cadabra2 script that computes or verifies it.
- **Modular:** Keep LaTeX chapters clean and pull in computed expressions from Cadabra2 via `Ex._latex_()` and custom export scripts.
- **Reproducible:** All computations must be runnable in the `cadabra2-ubuntu:24.04` Docker environment.

## Reading List (Priority Order)

### Must-Read Papers
1. **Parke-Taylor (1986)** — Original MHV paper. Short and foundational.
2. **BCFW Recursion (Britto, Cachazo, Feng, Witten 2005)** — Core of Ch.05.
3. **The recent MHV paper that used ChatGPT** — (provide arXiv number when available).
4. **Elvang-Huang Review** (if available) — Excellent modern reference on spinor-helicity.
5. **Dixon's MHV Review** — Classic reference.

### Srednicki Chapters (Target Completion Order)
- Ch. 48–52: Spinor-Helicity basics
- Ch. 60+: Scattering amplitudes and advanced topics
- Revisit Ch. 36, 37, 38 (as we did previously with Cadabra2)

## Paper Structure (Target)

- **Part I:** Foundations (Spinor-Helicity)
- **Part II:** MHV Amplitudes (Parke-Taylor, single-minus, why MHV?)
- **Part III:** BCFW On-Shell Recursion
- **Part IV:** Computational Tools (Cadabra2, FORM, SymPy integration)
- **Part V:** New Results / Open Directions

## Cadabra2 Integration Strategy

We will create two scripts per major topic:
- `chNN_topic.py` — performs the symbolic computation
- `chNN_export_latex.py` — converts results to clean LaTeX using `Ex._latex_()` and `mat2pmatrix()`

All scripts will live in `Cadabra2/` and be referenced from the LaTeX document.

## Immediate Next Steps (Next 24h)

1. Create `feat/mhv-paper-v2` branch
2. Finish reading Parke-Taylor + BCFW papers
3. Expand Ch.05 with BCFW recursion (using Cadabra2 to verify recursion relations)
4. Integrate existing Cadabra2 MHV results from Monoclaw into the LaTeX document
5. Update `99-master.tex` with better structure and cross-references

## Milestones

- **Week 1:** Complete Ch.05 draft + Cadabra2 verification of BCFW
- **Week 2:** Full draft of Parts I–III with all symbolic results embedded
- **Week 3:** Tools chapter + clean compilation
- **Week 4:** Polish, figures, and submission readiness

---

**Last Updated:** 2026-04-04  
**Owner:** Cyclonus + Ernest  
**Status:** Active Planning Phase