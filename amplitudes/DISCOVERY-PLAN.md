# Amplitudes Research — Agent Task Board

## Mission
Make a physics discovery in MHV scattering amplitudes, similar to arXiv:2602.12176v2 (Guevara, Lupsasca, Skinner, Strominger, Weil / OpenAI).

That paper showed single-minus gluon tree amplitudes are nonzero in the half-collinear regime and found a closed-form formula conjectured by GPT-5.2 Pro. We want to do something analogous.

## Background Reading (in priority order)

### Currently Reading
1. **Dixon 2013** — arXiv:1310.5353 — "A brief introduction to modern amplitude methods"
   - LaTeX source: `workspace2/Data/Public/papers/physics/arXiv-1310.5353v1/ModAmpIntro.tex`
   - Our notes: `amplitudes/08-literature-review.tex`
   - Covers: color decomposition, spinor-helicity, factorization, Parke-Taylor, BCFW, unitarity

2. **SMGA 2025** — arXiv:2602.12176v2 — "Single-minus gluon tree amplitudes are nonzero"
   - LaTeX source: `workspace2/Data/Public/papers/physics/arXiv-2602.12176v2/SMGA.tex`
   - Our notes: `amplitudes/08-literature-review.tex` (to be added)
   - THE BENCHMARK: this is what we're trying to match/exceed

### Background (Srednicki QFT)
- Chapters 34–43: spinor formalism, Weyl/Dirac/Majorana, gamma technology
- Chapters 48, 50, 60: spinor-helicity, MHV amplitudes
- Our Cadabra2 computations: `amplitudes/Srednicki/individual_chapters/`

## Discovery Strategy

### Phase 1: Build Foundation (Current)
- [ ] Finish reading Dixon review (understand color decomposition, BCFW deeply)
- [ ] Read SMGA paper carefully (understand half-collinear regime, recursion, sign functions)
- [ ] Complete Srednicki chapters 48–60 computations in Cadabra2
- [ ] Implement Berends-Giele recursion numerically (Python)

### Phase 2: Reproduce & Verify
- [ ] Reproduce SMGA's stripped amplitudes A_{123} through A_{123456} numerically
- [ ] Verify their formula Eq.(16) for n=3..10
- [ ] Verify consistency conditions (cyclicity, KK relations, soft theorem)
- [ ] Implement in Cadabra2 where possible

### Phase 3: Explore & Discover
- [ ] Extend to other kinematic regions beyond R_1
- [ ] Look for patterns in multi-minus amplitudes (NMHV, N²MHV)
- [ ] Explore graviton amplitudes (paper mentions direct generalization)
- [ ] Explore supersymmetric extensions
- [ ] Look at loop-level structure
- [ ] Try to find simpler expressions via different variable choices
- [ ] Connect to celestial holography / Mellin transforms (paper mentions Lauricella functions)

### Phase 4: Write It Up
- [ ] Document findings in `amplitudes/` LaTeX notes
- [ ] Prepare for arXiv submission if discovery is significant

## Key Equations to Internalize

### Parke-Taylor MHV (Dixon/Srednicki)
A_n(1⁻,2⁻;3⁺...n⁺) = ⟨12⟩⁴ / (⟨12⟩⟨23⟩...⟨n1⟩)

### SMGA Key Formula (Eq. 16)
A_{1...n}|_{R₁} = (1/2^{n-2}) ∏_{m=2}^{n-1} (sg_{m,m+1} + sg_{1,2...m})

### BCFW Recursion
A_n = Σ_P A_L(z_P) · (1/P²) · A_R(z_P)

### Berends-Giele Recursion (SMGA Eq. A7)
F_{1...m} = (1/(p²_{1...m}+iε)) Σ_j [λ̃_{1...j} λ̃_{j+1...m}] F_{1...j} F_{j+1...m}

## Repo Locations
- Paper notes: `amplitudes/08-literature-review.tex`
- Reading list: `amplitudes/READING-LIST.md`
- Cadabra2 scripts: `amplitudes/Srednicki/individual_chapters/`
- Dashboard: `amplitudes/dashboard/qft-dashboard/`
- Spinor computations: `amplitudes/spinors/`

## For Other Agents
If you're an OpenClaw agent picking up this work:
1. Read this file first
2. Read `amplitudes/READING-LIST.md` for paper status
3. Read the SMGA paper LaTeX source (it's short and self-contained)
4. Check `amplitudes/Srednicki/` for existing computations
5. Never push to master — use feature branches

---

Last updated: 2026-04-05
