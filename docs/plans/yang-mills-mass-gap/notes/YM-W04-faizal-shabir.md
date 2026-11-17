# YM-W04 — Faizal–Shabir, arXiv:2606.19362v1: where the crossover step should be

Author: strategist, Claude Opus 5.5 (claude-opus55-strat), 2026-11-17. Status: **unreviewed; needs review by a different model family** (packet YM-W04, tracker local:19).

Source: M. Faizal, A. Shabir, "Reflection-Positive Construction of a Four-Dimensional SU(N) Yang-Mills Theory with Mass Gap and Confinement", arXiv:2606.19362v1 (submitted 9 June 2026; the only version on 2026-11-17). The PDF (593 pages: the main paper followed by four companion papers) and the LaTeX source (one file, `main.tex`, 21,093 lines) are in `Data/Public/papers/physics/` in the workspace (not in this repo):

    2606.19362v1.pdf            sha256 0648f30f67719a6c83976f1da53131fc6f9a9ddc9ba1c406345f5ed959b32939
    arXiv-2606.19362v1.tar.gz   sha256 9e1d40ae4e65880e6b90af41de4f718993cd9f741a1fff8b9cca28e0c379f620
    arXiv-2606.19362v1/main.tex (extracted)

Citations below are `main.tex:line` and PDF page (p.).

## 1. Result

1. **The paper starts at strong coupling and coarse-grains.** The input is the strong-coupling gap at a fixed lattice spacing a₀ and a fixed β < β⋆(N) (main.tex:4411; companion paper II, main.tex:13897–13930). The renormalization-group step makes the lattice coarser: a_{k+1} = b·a_k with b ≥ 2 (main.tex:2267, p.54–55; main.tex:4461). The interlacing identity T_{k+1} = V_k* T_k^b V_k − D_k + E_k (main.tex:2249) compares the coarse transfer operator with b fine steps, so it only makes sense in that direction. Coarse-graining one fixed strong-coupling lattice theory moves towards long distances at fixed bare β. It never reaches the continuum limit (a → 0, β → ∞), which is where the strong-coupling expansion stops converging. The gap of the fixed-β lattice theory was already known (Osterwalder–Seiler 1978, cited by the paper as `OS-gauge`; from memory, unverified). Nothing in this mechanism goes from strong coupling to weak coupling.
2. **The paper contradicts itself on the direction of the scales.** Elsewhere it takes a_k = b^{−k} a₀ → 0 (main.tex:2671, p.68; main.tex:3532, p.96) and builds the continuum limit on a_k ↓ 0 (main.tex:3111, main.tex:4551). The gap-transport lemma closes only by asserting that "the sequence (a_k) is geometric with common ratio b, so a_{k+1} is uniformly comparable to 1 … a_{k+1}^{−1} ≤ c₁ and e^{a_{k+1}} ≤ c₂" (main.tex:2277, p.55). A geometric sequence with ratio b ≠ 1 has no such two-sided bound. If the lattice gets finer (a_k → 0), c₁ does not exist. If it gets coarser (a_k → ∞), c₂ does not exist and there is no continuum limit. In the main theorem the same step appears as "c_k := e^{a_kΔ_k}/a_{k+1} uniformly bounded along admissible steps" (main.tex:4435).
3. **The weak-coupling section reverses asymptotic freedom.** It claims that after finitely many *coarse-graining* steps the flow enters a neighbourhood of the Gaussian fixed point, with g_{k+1} = g_k − β₀g_k³ + R_k, β₀ > 0 of the one-loop sign, and a coupling that decreases monotonically along the trajectory (main.tex:3857–3866, p.109–110). For an asymptotically free theory the coupling *grows* under coarse-graining: to one loop, g_{k+1}^{−2} = g_k^{−2} − 2b₀ log b. A coupling that decreases towards the infrared is infrared freedom, the behaviour of φ⁴₄, not of Yang–Mills. So the claimed meeting of the strong-coupling and weak-coupling halves rests on a flow with the wrong sign.
4. **The spectral-gap theorem is conditional.** As stated (main.tex:4397–4408, p.130), it assumes (ii): a one-step interlacing inequality with summable defects along the whole flow to the continuum. That assumption is the hard part of the problem.

Conclusion: the paper does not cross from strong to weak coupling. Results 1–3 are each enough to stop its continuum argument. This is not a judgement on every lemma in the 593 pages.

## 2. Labels

| Statement | Label |
|---|---|
| Paper's claims (gap, area law, universality, AF matching) | source-claims |
| 1: base input is fixed-β strong coupling; RG step coarsens | checked (read main.tex at cited lines) |
| 1: fixed-β coarse-graining cannot produce the continuum limit | argued (standard; no computation needed) |
| 1: Osterwalder–Seiler already gives the fixed-β lattice gap | from memory, unverified |
| 2: internal contradiction a_{k+1}=b a_k vs a_k=b^{−k}a₀; false uniform bound on a_{k+1}^{±1} | checked |
| 3: sign of the coarse-graining flow is opposite to one-loop asymptotic freedom | checked against the paper's own display; the AF sign is standard (Gross–Wilczek, Politzer; from memory) |
| 4: theorem assumes (ii) | checked |

FINDINGs for the ledger: results 2 and 3 (exact lines, mechanical to verify). Result 1 is the structural reason why fixing 2 would not save the argument.

## 3. Evidence

    curl -sSL -o 2606.19362v1.pdf https://arxiv.org/pdf/2606.19362v1
    curl -sSL -o arXiv-2606.19362v1.tar.gz https://arxiv.org/e-print/2606.19362v1
    tar xzf arXiv-2606.19362v1.tar.gz -C arXiv-2606.19362v1
    grep -n -E '\\(section|subsection)\{' main.tex          # outline
    pdftotext -f N -l N -layout 2606.19362v1.pdf -          # page location of quoted passages

Passages read in full: abstract and introduction (main.tex:758–765); §2 opening and OS positivity (772–870); the gap-transport lemma and its proof (2247–2284); §10 opening and Hypothesis "Stated for convenience" (3853–3880); §11 "Theorems and Full Explicit Proofs", all three proofs (4391–4602); conclusion (4604–4629); companion paper II §"Transport from strong coupling to the scaling window", opening and strong-coupling lemma (13889–13960).

The paper does not mention OpenAI or the OpenAI Lean construction (`grep -ci openai main.tex` → 0). It cites Bałaban's 1984 papers only as an analogy (main.tex:2284).

## 4. Not done / not checked

- About 95% of the text was not read: the FRD construction (§4), the defect-summability proofs (§5.2), the continuum OS verification (§§7–8), companion papers I, III, IV, and the appendices.
- Not checked whether the "horizon projector" P_σ changes the gauge-invariant measure. If it does, the theory constructed is not Wilson's Yang–Mills even at fixed a. The paper says it acts on gauge-fixed slices (main.tex:4607); this needs one careful read.
- Not checked: whether a different choice of a_k direction makes some version of the interlacing argument consistent. My view is that none does (see result 1), but this is argued, not proved.
- The Osterwalder–Seiler statement in result 1 is from memory.
- Style note, not evidence: the abstract switches to the first person singular ("in my view") in a two-author paper.

## 5. Questions for the strategist

- Should results 2–3 go to `ledger/findings/` once a different family has reviewed them? An upstream report would be an external action and needs Ernest's go-ahead; I recommend none: Clay still lists the problem as open, and the errors are visible to any expert reader.
- Suggest closing YM-W04 after review, with claim C25 relabelled "refuted as a proof (crossover step absent; see YM-W04)".
