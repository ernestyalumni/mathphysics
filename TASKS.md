# Current Research Tasks — MHV Amplitudes Paper

**Primary Goal:** Produce a high-quality arXiv paper on MHV amplitudes, spinor-helicity formalism, and BCFW recursion, with heavy use of Cadabra2 for symbolic verification.

**Philosophy:** Every major claim in the paper must be backed by a runnable Cadabra2 computation.

## Immediate Priorities

### 1. Reading & Foundation (Next 6–8 hours)
- Read Parke-Taylor (1986) original MHV paper
- Read Britto-Cachazo-Feng-Witten (BCFW) 2005 paper
- Review Srednicki Ch. 48–52 (spinor-helicity)

### 2. Cadabra2 Integration (Active)
- Pull existing MHV/Parke-Taylor computations from `Monoclaw/Python/Cadabra2/Srednicki/`
- Create dedicated `Cadabra2/MHV/` directory in this repo
- Write `mhv_parke_taylor.py` + `mhv_export_latex.py` pair

### 3. Document Work (This Week)
- Expand `amplitudes/Ch.05/bcfw.tex` with BCFW recursion
- Add Cadabra2 verification of the recursion relations
- Improve cross-referencing in `99-master.tex`

### 4. Long-term Paper Goals
- Clean, publication-quality LaTeX (using the existing Preamble macros)
- Include computed expressions directly from Cadabra2
- Add figures (Feynman diagrams, collinear limits, etc.)
- Discuss computational tools (Cadabra2, FORM, SymPy)

## Branching Rule
- Work on `feat/mhv-paper-v2` (or similar)
- Never commit directly to master

---

**Last Updated:** 2026-04-04  
**Status:** Planning & Reading Phase