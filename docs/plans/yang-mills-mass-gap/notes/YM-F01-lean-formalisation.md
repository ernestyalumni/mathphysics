# YM-F01 — Lean formalisation of the conditional-gap chain, and a local-prover test

Author: Claude Opus 5.5, 2026-11-20 UTC. Branch `ym/YM-F01-claude-opus55-lean`. Request: Ernest, 2026-11-20, "Ok let's give 1,2,3 a try please" (`request-2026-11-20-lean-formalisation.txt`).

## Result

**Steps chosen and proved.** All 18 theorems compile in Lean 4.35 with Mathlib, in Docker with no network (`formal/check.sh --strict`). There is no `sorry`, and the main theorems use only the axioms propext, Classical.choice and Quot.sound. The theorems formalise these parts of `documents/mono-notes/topics/yang-mills-lattice.tex`:

| Source | Lean (`formal/YMFormal/`) | What is proved |
|---|---|---|
| `ym:prop:purity`, Eq. `ym:purity-bounds` | `Purity.purity_bounds` | `𝖯(n) = (1+q_n)/(1+s_n)²`, `𝖯 ≤ 1`, `s_n ≤ 𝖯⁻¹-1 ≤ 2s_n+s_n²`, `s_{2n} ≤ s_n²` for any summable spectrum `0 ≤ λ_j ≤ λ₁`, `λ₁ > 0` |
| Eq. `ym:gap-from-purity` | `Purity.exp_neg_gap_le`, `Purity.gap_from_purity` | `e^{-n𝗆} ≤ 𝖯⁻¹-1` and `𝗆 ≥ (1/n) log(1/(𝖯⁻¹-1))` |
| `ym:prop:defect`, Eq. `ym:defect-telescope` | `Defect.defect_telescope`, `Defect.defect_nonneg_and_last` | `e^{-Δ} = 𝖯₁⁸𝖯₂⁴𝖯₃²𝖯₄`, `Δ ≥ 0`, `𝖯₄ ≥ e^{-Δ}` |
| Eq. `ym:gap-from-defect` | `Defect.gap_from_defect` | `𝗆 ≥ (1/n) log(1/(e^Δ-1))` |
| `ym:prop:conditional-gap` | `ConditionalGap.conditional_gap` | under the recurrence, `𝗆(2^{k+1}n) ≥ log 2 / n` for all `k` |
| `ym:prop:casimir-sun` | `Casimir.trace_sandwich`, `Casimir.trace_adjoint_comp`, `Casimir.casimir_of_killing` | `tr(Z ↦ AZB) = tr A tr B`; `tr(ad X ad Y) = 2N tr(XY) - 2 tr X tr Y` on `M_N`; `C_𝔤 = 2N²` from the Killing relation |
| values after Eq. `ym:slope` | `Casimir.slope_values` | `b_* = 11/(3π²)` for SU(2), `33/(4π²)` for SU(3) |

**Statements read against the source.** I wrote the statements myself and checked them against the tex. No other reviewer has checked them.
- **Differences from the source:**
  - **`λ₂ > 0` is a hypothesis of the gap bounds.** It is unstated in the source, which needs it for `𝗆 = log(λ₁/λ₂)` to be finite.
  - **Positive defects are a hypothesis of `conditional_gap`.** Zero would mean a purity of 1 and an infinite gap.
  - **The source's side conditions `𝖯(n) > 1/2` and `Δ < log 2` are not needed.** They only make the bounds positive.
- **Hypotheses, not proved here:**
  - **Spectrum.** That the spectrum is that of the transfer operator, and `Z(n) = Tr 𝖳ⁿ` (`ym:trace-formula`).
  - **Purities.** That each `𝖯_k` is a purity, which uses the symmetry of `Z`.
  - **Recurrence.** The defect recurrence, Conjecture `ym:conj:recurrence`, which is open.
  - **Casimir.** For `casimir-sun`, that the real trace over `su(N)` equals the complex trace over `sl_N(ℂ)`. The strategist checked `K(X,Y) = 2N tr(XY)` on `su(N)` numerically for N = 2 to 5 (200 instances); that is evidence, not a proof.
- **What this establishes.** The finite and algebraic steps that carry the conditional gap are correct as stated. The source has no error in them.
- **The `ym-s04` step 9 values.** The local 9B flagged them as wrong (YM-T16). Given `C_𝔤 = 2N²`, which is proved here up to the complexification step, the values are right.

## Step 3: local prover (pre-registered before any generation)

**Question.** Can a local Lean prover on the RTX 3060 fill in any of these proofs?

**Model.**
- Goedel-Prover-V2-8B (Goedel-LM, Apache-2.0), GGUF Q6_K from `mradermacher/Goedel-Prover-V2-8B-GGUF`, SHA-256 `43a23c09…71aff3a`.
- Served by the pinned llama.cpp image used for YM-T16, with context 24,576, flash attention and one slot.

**Benchmark.** The 18 statements of commit `a8a31e3`, frozen before any proof was written. `formal/prover.py tasks` builds the tasks file, SHA-256 `304173d9…db914a3`.
- **What the model sees for each theorem:**
  - The module up to that theorem, with my proofs of the earlier theorems.
  - Imported modules with their proofs omitted.
  - The theorem with `sorry`, inside the model card's prompt, verbatim.
- **What is kept from its answer:** only the proof of that theorem, and any new `theorem`/`lemma` helpers. They go under my statement text, so the model cannot change what is proved.

**Sampling.** 4 attempts per theorem, temperature 0.6, top-p 0.95, top-k 20, min-p 0, at most 16,384 new tokens. Seeds are `20261120 + 1000·attempt + index`.

**Judge.** Lean only. A candidate passes if it compiles in Docker with no network, and `#print axioms` for the theorem lists nothing beyond propext, Classical.choice and Quot.sound. That makes `sorry`, new axioms and `native_decide` failures.
- Controls, run before the model:
  - My own proofs pass.
  - A bare `sorry`, a helper with `sorry`, and an injected `axiom` all fail (3 theorems × 4 controls).

**Primary outcome.** Theorems with at least one passing attempt, out of 18 (pass@4). Secondary: passes per theorem, and why attempts fail (no code block, cut off at the token limit, Lean error).

**Expected difficulty.** The model was trained on an older Mathlib, and some lemma names have changed since (`div_le_iff` → `div_le_iff₀`, `tsum_le_tsum` → `Summable.tsum_le_tsum`). Name errors count as failures; that cost is part of what is measured.

## Labels
- The Lean proofs are machine-checked. The match between the Lean statements and the source is my own reading, listed above.
- Model output is labelled `unreviewed (local model output; not evidence)`. The Lean verdict on it is a check, not a review.

## Resume here
Rerun the proofs with `formal/check.sh --strict`. The prover run: `formal/prover.py generate`, then `judge`, then `report`, writing to a directory outside git. Results go in Step 3 above.
