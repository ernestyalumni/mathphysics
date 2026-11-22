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
  - **Recurrence.** The defect recurrence, Conjecture `ym:conj:recurrence`, which is open. The 4D SU(2) measurements in `notes/YM-N03-defect-measurement.md` (branch `ym/YM-N01-claude-opus55-su2`, be8e89a) find Δ small at n = 2 only at strong coupling, where the 2⁴ torus is confined. So the hypothesis "Δ small at some scale" is roughly finite-size confinement at that scale, and the theorem's content is the step n → 2n, which is this recurrence.
  - **Casimir.** For `casimir-sun`, that the real trace over `su(N)` equals the complex trace over `sl_N(ℂ)`. The strategist checked `K(X,Y) = 2N tr(XY)` on `su(N)` numerically for N = 2 to 5 (200 instances); that is evidence, not a proof.
- **What this establishes.** The finite and algebraic steps that carry the conditional gap are correct as stated. The source has no error in them.
- **The `ym-s04` step 9 values.** The local 9B flagged them as wrong (YM-T16). Given `C_𝔤 = 2N²`, which is proved here up to the complexification step, the values are right.

## Comparison with Codex's review (YM-O03 §1, fdaccd7 on `ym/YM-O03-codex-review`)

Requested by Ernest, 2026-11-21: "compare §1 with F01, especially which hypotheses Lean assumes and how to cover the rank-one case at β=0".

**The review's reading of the four modules agrees with the code.**

**Hypotheses the Lean theorems assume:**
- **Spectrum.** The Lean spectrum `lam` is only `Spectrum lam`: `0 < lam 0`, `0 ≤ lam j ≤ lam 0`, summable. The tail is not ordered. To read `lam 1` as the physical λ₂, enumerate the spectrum in decreasing order with multiplicity. That is a bridge hypothesis, added to the list above.
- **Log forms.** These add `0 < lam 1`, `0 < Δ` or `0 < δ k`.
- **Trace identifications.** Two are inputs, not theorems: `hP : exp (-Δ) ≤ purity lam n` and the `𝖯_k ≤ 1` bounds. No transfer operator is constructed.
- **Conditional gap.** `conditional_gap` takes abstract sequences, `from_defect` at every k, and the recurrence.
- **Casimir.** The Lean proves the matrix identities over any commutative ring. Restriction to sl and complexification of su are manuscript steps. The statement checker's self-check confirms `K(X,Y) = 2N tr(XY)` exactly on 1,066 random su(N) elements, N = 2–4. That is evidence, not a proof.

**Rank one (β = 0) is now covered.** The log forms exclude it by hypothesis; they are not wrong there. `formal/YMFormal/RankOne.lean` (in `check.sh --strict` and `Axioms.lean`) restates the chain for `r = λ₂/λ₁` with no positivity assumptions:

| Theorem | Statement |
|---|---|
| `ratio_pow_le_impurity` | `(λ₂/λ₁)ⁿ ≤ 𝖯⁻¹ - 1` |
| `ratio_pow_le_defect` | `e^{-Δ} ≤ 𝖯(n) → (λ₂/λ₁)ⁿ ≤ e^Δ - 1`, for any Δ |
| `rank_one_of_zero_defect` | Δ = 0 forces λ₂ = 0 |
| `conditional_gap_ratio` | `0 ≤ δ_k`, recurrence, start, and `r_kᵉ ≤ e^{δ_k} - 1` with `e = 2ᵏ n` imply `r_kⁿ ≤ 1/2` for all k, with no `n ≥ 1` |

The gap is `-log r` when `r > 0` and infinite when `r = 0`, so `log 2 / n ≤ 𝗆` follows. A zero defect propagates zero under the recurrence, as the review says.

**Codex's rectangle recurrence, read against F01.** I checked it by hand. It is not formalised.
- **Estimate 1 comes from `purity_bounds`.** The three bounds needed are `s ≤ 𝖯⁻¹ - 1`, `s₂ₙ ≤ sₙ²` and `𝖯⁻¹ - 1 ≤ 2s + s²`, the last giving `log 𝖯⁻¹ ≤ 2s`. Then `2s² ≤ 2p²e^{2p} ≤ 4p²` needs `e^{2p} ≤ 2`, which holds for `p ≤ 1/4`.
- **Estimate 2 is an identity of F.** Both sides expand to `2F(Dⱼx) - F(DᵢDⱼx)`.
- **The induction holds.**
  - The path weights 8, 4, 2, 1 sum the edge costs to Δ(n), so each edge cost is at most δ.
  - In the step, dropping `-2pₕ(y) ≤ 0` gives `pᵢ(x) ≤ 2pᵢ(y) + pₕ(Dᵢy)`.
  - The second term is a doubling of a rectangle `z` with `zₕ = n` and j − 1 long sides, so Estimate 1 bounds it by `4A_{j-1}²`.
  - That gives `A_j ≤ 3A_{j-1}`.
  - The constant is `4(8 + 36 + 162 + 729) = 3740`.
  - `δ ≤ 1/7480` implies `27δ ≤ 1/4`.
- **The form matches F01's `recurrence` exactly.** With `c₁ = 3740` and threshold `δ ≤ 1/7480 = 1/(2c₁)`, it is `δ k ≤ 1/(2c₁) → δ (k+1) ≤ c₁ δ k²`. If the derivation stands, the open hypothesis is discharged.
- **The cost: more spectral input than F01 now assumes.** F01 uses one spectral purity per scale. The rectangle argument needs every `pᵢ(x)` to be the log impurity of a positive transfer operator in direction i, for every rectangle with sides in {n, 2n, 4n}. That is reflection positivity in all four directions (C23).
- **Formalising it.** With `F : (Fin 4 → ℕ) → ℝ` abstract and Estimate 1 and `pᵢ ≥ 0` as hypotheses, the induction is finite and could be checked in Lean. Not started.
- **Unchanged.** A small defect at bounded physical size (`start`) remains open.

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

### Step 3 result (run 2026-11-20 13:04–18:08Z; unreviewed local model output, checked by Lean)

**Primary outcome: 4 of 18 theorems solved (pass@4).** 13 of 72 attempts pass.

| Solved | Passes | What the theorem is |
|---|---|---|
| `Casimir.casimir_of_killing` | 4/4 | solve `2Nt = C t/N` for `C` |
| `Casimir.slope_values` | 4/4 | evaluate `11·C/(48π²)` for N = 2, 3 |
| `Defect.defect_telescope` | 3/4 | `e^{-Δ} = z₅/z₁¹⁶` written as the product of the `𝖯_k` |
| `Purity.excited_double` | 2/4 | `(x^{2n}) = (x^n)²` inside a `tsum` |

**Unsolved, 14 theorems.**
- **Every theorem that needs an inequality over a `tsum`, a chain of real inequalities, or the trace of a linear map.** That includes all four steps the conditional gap rests on: `purity_bounds`, `gap_from_purity`, `gap_from_defect` and `conditional_gap`.
- **Why the 59 failing attempts fail:**

  | Reason | Attempts | Detail |
  |---|---|---|
  | Lean error | 38 | 3 of them were cut off at 16,384 tokens |
  | No code block | 15 | all cut off at 16,384 tokens |
  | Target theorem missing | 6 | the answer proves a different statement, or stops before the target |

- **Cut-offs:** 21 of 72 answers reached the token limit. The median answer used 6,366 tokens.

**Name drift (exploratory, defined after seeing the data).**
- **In 12 of the 38 Lean-error attempts, every reported error is an unknown identifier.** These are renamed or removed Mathlib lemmas: `le_div_iff` (9), `lt_div_iff` (8), `pow_le_one` (5), `div_le_div_iff` (4) and others.
- **In 18 of the 38, the first error is an unknown name.**
- **This does not mean those 12 proofs would compile under the old names.** Lean stops checking a branch at the first unknown name, so later errors stay hidden.
- **Experiment 1 tests it directly:** it feeds those errors back to the model.

**Reading.**
- **The model proves the routine algebraic steps.** I could also do those by hand in a few lines.
- **It proves none of the analytic steps of the chain.**
- **Lean caught every failure.** No `sorry`, new axiom or false statement got through, because of the judge's design. That is the case for a prover over an informal referee: a failed proof is visibly a failure, not a false alarm.

## Follow-ups, pre-registered before running (Ernest, 2026-11-20: "pre register and run those 3 that order")

All three start after the step-3 run above. Experiment 2 runs on the laptop, so it may overlap experiment 1. If any of these models proves effective, Ernest wants it deployed with our own Rust and CUDA C++ inference, not llama.cpp; nothing here depends on that.

**Experiment 1: Goedel-Prover-V2-8B with Lean error feedback.** This is the model card's self-correction mode.
- **Setup:**
  - Same model file, server and 18 tasks (tasks file `304173d9…`). Context 32,768, at most 20,480 new tokens per answer.
  - 4 attempts per theorem. After each failing answer, Lean's first three errors, each with the source line it points at, go back to the model, up to 2 revision rounds (`prover.py refine`, message `REVISE`). An answer with no usable code gets the `MISSING` message instead.
  - Round 0 reuses the step-3 seeds, so it is paired with step 3 apart from the larger token limit.
- **Judge:** unchanged.
- **Primary outcome:** theorems solved by any attempt and round, out of 18, against step 3's pass@4.
- **Secondary outcomes:** generations used, the round of each pass, and why attempts fail.
- **Time:** at most 216 generations, likely 6–10 h.

**Experiment 2: Goedel-Formalizer-V2-8B**, Q5_K_M, SHA-256 `123e5ad7…`, on the laptop's RTX 3070 with context 16,384.
- **Input:** the 7 natural-language statements in `formal/informal-statements.json`. Each matches Lean statements in `YMFormal` that I checked by hand.
- **Prompt and sampling:** the model card's prompt verbatim, temperature 0.9, top-k 20, top-p 0.95, at most 14,336 new tokens, 4 attempts each.
- **Scoring:**
  - **Compiles:** the last ```lean4 block compiles with `import Mathlib` in our Docker image (`formalizer.py judge`). `sorry` as the proof is expected.
  - **Faithful:** graded by me against the informal statement and the reference, as equivalent, stronger hypotheses, weaker hypotheses, or wrong. A statement made vacuous by contradictory hypotheses counts as wrong. The grading is not blind; grades go in `grades.json` with a reason for each.
- **Primary outcome:** statements with at least one attempt that compiles and is equivalent, out of 7.

#### Experiment 1 result (run 2026-11-20 18:08 – 2026-11-21 11:08Z; unreviewed local model output, checked by Lean)

**Primary outcome: 5 of 18 theorems solved, against 4 of 18 without feedback.** 16 passing answers from 189 generations; at most 216 were allowed, and attempts stop once they pass.

**What feedback added.** `ConditionalGap.defect_double_exponential`, solved by one attempt in its second revision. It is the induction that carries the conditional gap: under the recurrence, `c₁ δ_k ≤ (1/2)^(2^k)`. It is the first step of the chain that is more than routine algebra to be proved by a local model.

**Where passes came from.**

| Round | Passing answers |
|---|---|
| 0, no feedback, 20,480 tokens instead of 16,384 | 12 |
| 1 | 1 (`excited_double`) |
| 2 | 2 (`defect_telescope`, `defect_double_exponential`) |

**Unsolved: the 13 others.** These are every tsum inequality, the gap bounds, `conditional_gap` itself, `exp_sub_one_le_two_mul` and the trace lemmas. 56 of 189 generations reached the 20,480-token limit.

**Reading.** Error feedback helps a little, at the margin. Nearly all the passes are first answers. One hard step fell after two rounds of Lean errors, out of 216 allowed generations. That is not a prover that can carry this chain.

#### Experiment 2 result (run 2026-11-21 04:04–04:47Z on the RTX 3070; unreviewed local model output, graded by me)

**Primary outcome: 2 of 7 statements** have an attempt that compiles and is equivalent: `conditional_gap` (4 of 4 attempts) and `defect_telescope` (1 of 4). 9 of 28 attempts compile; answers were short (median about 950 tokens, none cut off).

| Statement | Compiles | Equivalent | Equivalent, ignoring compile failure |
|---|---|---|---|
| `casimir_su` | 0/4 | 0 | 0 |
| `trace_ad_ad` | 1/4 | 0 | 1 |
| `purity_bounds` | 0/4 | 0 | 2 |
| `gap_from_purity` | 0/4 | 0 | 1 |
| `defect_telescope` | 4/4 | 1 | 1 |
| `gap_from_defect` | 0/4 | 0 | 1 (+1 with weaker hypotheses) |
| `conditional_gap` | 4/4 | 4 | 4 |

**How the equivalent ones were checked.** For `conditional_gap` attempts 0 and 3 and `defect_telescope` attempt 0, I derived the model's statement from my proved theorems in Lean (and `conditional_gap` attempt 0 also in the other direction). The rest are graded by reading.

**How they go wrong.**
- **`λ` as a variable name, 12 attempts.** It is reserved in Lean 4, so every statement about the spectrum fails to parse. Read as if renamed, 5 of these 12 are faithful.
- **Optional parameters used as definitions, 7 attempts.** `(P₁ : ℝ := z₂ / z₁^2)` in a theorem header is a universally quantified variable with a default, not a definition. Three such `defect_telescope` statements compile and are false. Lean derives `False` from attempt 1 with P₁…P₄ = 0 and Δ = 1. A compile check alone would have passed them.
- **A compiling statement that is false.** `trace_ad_ad` attempt 2 states that the trace of a commutator equals `2N tr(XY) − 2 tr X tr Y`.
- **The Killing form of su(N) not stated in any attempt.** Every `casimir_su` attempt defines K as something other than a trace over su(N).
- **Wrong index ranges and a precedence slip.** In the index ranges, `Z` summed over j = 0 with λ₀ unconstrained, or dropped λ₁. The precedence slip: `log(1/(P n)⁻¹ - 1)` is `log(P − 1)`.

**Reading.** The formalizer states the simple real-number statements correctly, `conditional_gap` every time. For anything involving sequences, sums or linear maps, it needs a human to check each statement, and compiling is no evidence of faithfulness: 4 of the 9 compiling attempts are false statements. Grades and reasons: `grades.json` in the run directory (outside git).

**Experiment 3: Phi-4-reasoning-plus (14B)** as an informal-referee control.
- **Setup:**
  - Q4_K_M from `bartowski/microsoft_Phi-4-reasoning-plus-GGUF`, SHA-256 `fcb2cdd1…`.
  - The YM-T16 held-out cases (`referee-test.json`, `b90342fb…`), 2 samples each, run by claw-wsp `referee` with the YM-T16 protocol unchanged: temperature 0.6, top-k 20, top-p 0.95, cap 12,288, context 24,576, the same seeds.
  - It is compared with the YM-T16 control arm, the base Qwen3.8-9B Q5, by `referee-compare` and per case.
  - The model's recommended settings (temperature 0.8, top-k 50, up to 32K tokens) are not used, so that only the model changes. Cap hits are reported as the cost of that choice.
  - First, a one-case smoke test checks that its chat template works with the referee: a single `</think>` token, and thinking on.
  - **Deviation, 2026-11-21, before any generation.** The first launch failed to load: the f16 KV cache at context 24,576 needs 4.8 GB on top of 8.6 GB of weights, more than the 3060's 11.9 GB free. My estimate of about 10.5 GB was wrong. The run uses a q8_0 KV cache instead (`claw-wsp referee --cache-type q8_0`, claw-work 2026-11-21); everything else is unchanged. The control arm used an f16 cache, so the comparison changes the cache precision as well as the model.
- **Outcomes:**
  - **Primary:** errors located plus originals passed.
  - **Secondary:** mutants located, false alarms, cap hits, and two-sample agreement.

#### Experiment 3 result (run 2026-11-21 16:00–20:40Z; unreviewed local model output)

**Two harness fixes were needed before the result could be read.** Both are in claw-work `feat/weight-space-discovery` (f01d4f6, 0a83cc0).
- **Turn tags.** The referee looked only for Qwen's `<|im_start|>assistant`; Phi-4 renders `<|user|>…<|end|><|assistant|>`. The smoke test caught it.
- **Answer split.** Phi-4 ends its thinking with a special `</think>` token (id 100351) that the server leaves out of the returned text. Under the text split, every answer was unparsed, although it closed its thinking in 79 of 86 generations.
- **The fix.** The referee now takes the answer as the detokenized tokens after that token. `claw-wsp referee-regrade` applied it to the finished run, changing 79 rows and generating nothing new; Qwen runs, whose text contains `</think>`, parse as before.
- **The comparison.** `referee-compare --cross-model` compared the runs. It requires the same cases, seeds, sampling, cap, context and runtime image, but allows a different model, `</think>` token and tool build.

| 43 held-out cases × 2 samples | Qwen3.8-9B Q5 (YM-T16 control) | Phi-4-reasoning-plus Q4 |
|---|---|---|
| Correct (primary) | 39 / 86 | 53 / 86 |
| Originals passed | 10 / 28 | 25 / 28 |
| False alarms on originals | 11 | 2 |
| Planted errors located | 29 / 58 | 28 / 58 |
| Planted errors flagged | 46 / 58 | 30 / 58 |
| Cap hits | 12 | 7 |

**Primary difference: +16 points.** The 95 % bootstrap interval is 3–30 points and McNemar's p is 0.03. Both treat the two samples of a case as independent, which overstates certainty. By case:
- Phi-4 is right on both samples where Qwen is wrong on both in 6 cases.
- The reverse holds in 3.

**Reading.**
- **Where the gain comes from.** Phi-4 is a far less trigger-happy referee: 2 false alarms against 11.
- **Errors found.** It locates no more planted errors than the 9B, 28 against 29 of 58, and flags only half of them.
- **Caveats.** The q8_0 KV cache, and the chat template's built-in system prompt, differ from the control. A referee that misses half the planted errors still cannot certify a proof.

## Statement checker (pre-registered 2026-11-21, before the checker existed or ran on any graded statement)

Ernest, 2026-11-21: "we've gotta automate this, really can't rely upon me", then "ok go ahead"; and "banish the use of sympy/numpy" (`request-2026-11-21-statement-checker.txt`, PDD intent `automatic-lean-statement-checker-in-rust-3c8e1031`).

**Question.** Can a program decide, without a human reader, whether a Lean statement says what an informal statement says?

**Method: two translations compared on instances.**
- **Lean side, mechanical.** A Lean command (`formal/Audit/Export.lean`) exports the elaborated statement as JSON: binders, hypotheses, conclusion, with `YMFormal` definitions unfolded from their Lean values. The Rust evaluator applies Lean's own semantics: `x/0 = 0`, `log 0 = 0`, `log (-x) = log x`, truncated `ℕ` subtraction, and a `tsum` that is only decided for finitely supported terms. No text is parsed, so notation, precedence, coercions and `optParam` binders mean what Lean says they mean.
- **Informal side, written by hand.** Executable predicates in Rust for the 7 statements of `formal/informal-statements.json`: their hypotheses, and their conclusion split into conjuncts. I wrote these, and I also wrote the reference Lean statements, so the two sides are not independent readers. That is a limit of this test.
- **Arithmetic: exact or rigorous, never floating point.** All values are exact rationals (`qft-symbolics` `BigQ`/`BigC`), and `exp`/`log` are rigorous rational intervals. Each comparison is true, false or undecided. Crate `Rust/statement-audit`, CLI `ym-audit`; no SymPy, NumPy or Python.
- **Instances.** Fixed-seed random and edge cases:
  - finitely supported spectra with ties, zeros or a single nonzero term;
  - doubly exponential sequences;
  - small integer and Gaussian-rational matrices, N = 1–4;
  - special reals (0, 1/4, log 2 bounds).

  Universal statements over `ℕ` inside a statement are evaluated for the first 8 values.
- **Variable matching.** Lean variables are matched to informal ones by name with fixed aliases (`λ`, `lam`: the spectrum; `m`, `gap`; `δ`; …). For sequences, index offsets 0 and 1 are both tried, and the mapping with the fewest disagreements is kept.
  - Variables in the Lean statement only are searched for values that satisfy the Lean hypotheses.
  - A variable fixed by a defining hypothesis (`v = e`, `∀ x, v x = e`) is set from it.
- **Verdicts for each Lean statement.**
  - **FALSE:** an instance satisfies the Lean hypotheses but not the Lean conclusion.
  - **VACUOUS:** no instance found satisfying the Lean hypotheses.
  - **STRONGER:** an instance satisfies the informal hypotheses, and no value of the Lean-only variables satisfies the Lean ones.
  - **MISMATCH:** a Lean conclusion conjunct agrees with no informal conjunct on every instance where both are decided, including instances outside the hypotheses.
  - **MISSING:** an informal conjunct has no Lean conjunct matching it (formalizer statements only; for the reference theorems, coverage is reported across the set).
  - **UNMATCHED:** an informal variable has no Lean counterpart.
  - **UNDECIDED:** fewer than 50 instances decided.
  - **PASS:** none of the above. A note, **WEAKER**, records Lean hypotheses that hold where the informal ones fail.
- **Hypothesis scan.** Each hypothesis of each informal statement is dropped in turn, and a counterexample is searched for. The result is compared with the hand list above: λ₂ > 0 is needed; `P(n) > 1/2` and `Δ < log 2` are not.
- **Self-check of the informal statements.** Any instance satisfying an informal statement's hypotheses but not its conclusion is reported as an error in the statement (or in my predicate). For `casimir_su`, this tests `K(X,Y) = 2N tr(XY)` exactly on su(N), for N = 2–4, which none of the Lean theorems states.

**Benchmark, all fixed now.**
- **(a) The 28 Goedel-Formalizer statements, with my grades** (`grades.json`). `λ` is first renamed to `lam` automatically, and statements that still do not elaborate are excluded. Graded equivalent should PASS; graded weaker hypotheses should PASS with WEAKER; graded wrong should be flagged.
- **(b) The 9 reference theorems with informal counterparts:** `purity_bounds`, `exp_neg_gap_le`, `gap_from_purity`, `defect_telescope`, `defect_nonneg_and_last`, `gap_from_defect`, `conditional_gap`, `trace_sandwich`, `trace_adjoint_comp`. All should PASS, and together they should cover every conjunct of their informal statements except `casimir_su`'s.
- **(c) Mutants of (b)**, generated mechanically by single text edits to the statement:
  - `≤`↔`<`;
  - a numeral changed;
  - one hypothesis dropped;
  - an index shifted (`j + 1` → `j`, `lam 0` → `lam 1`, `k + 1` → `k`);
  - a binary `-` swapped with `+`;
  - the sign inside `exp (-…)` removed.

  Mutants that fail to elaborate are discarded. A mutant counts as detected if its report differs from the unmutated statement's report in any flag or WEAKER note.

**Primary outcome.**
- On (a): misses (wrong but PASS) and false alarms (equivalent but flagged).
- On (c): the detection rate.

**Useful if:** (a) has no misses and at most one false alarm, (b) all pass, and (c) detects at least 90 %.

**Allowed after this point:** fixes for constructs the evaluator does not support, or crashes, logged with the reason. **Not allowed:** changes to the verdict rules, the matching rule, the thresholds or the instance generators after any benchmark result is seen. The checker is developed on toy statements in `formal/Audit/Toy.lean`, not on the benchmark.

### Statement checker result (benchmark 2026-11-21 21:45–22:15Z, 900 instances per mapping, seed 20261121)

**Against the pre-registered bar: not useful as it stands.** It fails two of the three criteria. **It passed no wrong statement.**

| Part | Result | Bar |
|---|---|---|
| (a) 20 elaborating formalizer statements (8 do not elaborate even after `λ → lam`) | **0 misses**, 4 false alarms; 10 of 10 wrong statements flagged | 0 misses, ≤ 1 false alarm |
| (b) 9 reference theorems | 8 pass; `gap_from_purity` flagged MISMATCH | all pass |
| (c) 167 mutants (all elaborate) | 152 detected, 91 % | ≥ 90 % |

**Coverage (b).** Together the references cover every informal conjunct except `gap_from_purity`'s second (`P > 1/2 → …`) and both of `casimir_su`'s. That matches the Lean.

**Informal self-check and hypothesis scan.**
- **Self-check.** No failures. `K(X,Y) = 2N tr(XY)` holds exactly on all 1,066 su(N) instances, N = 2–4.
- **First run.** It reported 206 failures for `trace_ad_ad`. All came from ill-typed instances, where perturbation changed N but not the matrices; see fix 3 below.
- **Hypothesis scan, compared with the hand list.**
  - `Δ < log 2` is not needed.
  - `0 ≤ λ_j ≤ λ₁`, `P(n) ≥ e^{-Δ}`, the start and the recurrence are needed.
  - `λ₂ > 0` gives no counterexample, because without it the statement is undefined rather than false.
  - `c₁ ≥ 2`, `n ≥ 1` and `δ_k > 0` were never violated alone by the instances, so they are untested.

**Why the false alarms happen (post hoc, read off the report).** All three mechanisms are conservative: they make the checker refuse to vouch, never vouch for a false statement.
1. **Defining hypotheses are re-checked: 3 false alarms.** `gap_from_purity#3` and `gap_from_defect#1/#2` define `Z` by `∀ m, Z m = ∑' j, …`. The checker sets `Z` from that equation, then also evaluates the equation for `m = 0…7`. At `m = 0` the sum `∑ λ_j⁰` diverges, so it is undecided, the Lean hypotheses are never true, and the statement is reported VACUOUS.
2. **Random search for Lean-only values: 1 false alarm.** `purity_bounds#0` bounds an unused `λ₀` by `λ₁`. With `λ₁ ≈ 6 · 10⁻⁸`, 20 random tries did not find a value in `[0, λ₁]`, so it is reported STRONGER.
3. **Conclusions must agree exactly: 1 reference flagged.** `gap_from_purity` proves the bound without the premise `P > 1/2`, a documented and legitimate strengthening. The rule compares truth values, so a stronger Lean conclusion is a mismatch.

**Undetected mutants.** 15 of 167. Most make the statement weaker or strengthen a hypothesis only near a boundary:
- `0 < Δ` → `0 ≤ Δ`;
- dropping `hn`;
- the recurrence threshold `1/(2c₁)` → `2/(2c₁)`.

Either the instances rarely reach that boundary, or the parent already carried the WEAKER note, so the signature did not change.

**What would fix it.** None of this is applied: the benchmark is spent, and a fixed checker needs a fresh pre-registered test set.
- Treat a defining hypothesis as satisfied by construction.
- Search Lean-only variables over boundary values (0, the endpoints, the informal values) before random ones.
- Accept a Lean conclusion that implies the informal one under the informal hypotheses, reported as STRONGER-CONCLUSION.

**Changes after pre-registration.** All are fixes for crashes, unsupported constructs or ill-typed instances; none changed a verdict rule, the matching rule, a threshold or the instance generators.
1. **Mutant file imports.** The Defect mutant file lacked `import Audit.Export`.
2. **Mutant line ranges.** These are used to attribute Lean errors; they were off by one.
3. **Informal side, ill-typed matrices.** Instances whose matrices are not N × N are now undecided. Found by the self-check.
4. **Unknown constants.** Instance arguments are now opaque placeholders instead of errors. Found after seeing 3 ERROR results, for `trace_adjoint_comp` mutants, in an aborted first benchmark run that produced 4 results.
5. **`LinearMap.trace` module argument.** It is now located by type, not position.
6. **Lean side, ill-typed matrices.** These instances are now skipped. Found after seeing reference `trace_sandwich` wrongly flagged FALSE on an N = 0 instance with 1 × 1 matrices, during development.
7. **Memoised `exp`/`log`.** A performance fix with no effect on values.

The 600-bit tie rule was set during development on toy statements, before the benchmark.

## Rectangle recurrence in Lean (2026-11-21, Ernest: "ok sounds good go ahead please")

**What is proved.** `formal/YMFormal/Rectangle.lean` proves Codex's derivation (YM-O03 §1), and `check.sh --strict` passes with axioms propext, Classical.choice and Quot.sound only. Boxes are indexed by exponents: `G e` is `log Z` of the box with sides `n·2^{e_k}`.
- **Hypotheses.** `nonneg`: every doubling cost `pᵢ = 2G(e) - G(Dᵢe)` is at least 0. `estimate`, the review's Estimate 1: `pᵢ ≤ 1/4 → pᵢ(Dᵢe) ≤ 4pᵢ²`.
- **`rectangle_recurrence`.** `Δ(n) ≤ 1/7480 → Δ(2n) ≤ 3740 Δ(n)²`. The steps:
  - Estimate 2, `cost_swap`, is proved as an identity.
  - The `{n, 2n}` base case uses a path potential, `φ(e) = 16G(e)/2^{Σe}`.
  - The induction gives `pᵢ ≤ 3^j Δ(n)` for j sides of length 4n (`cost_le_pow_three`).
  - The final path gives the constant `4(8 + 36 + 162 + 729) = 3740`.
- **`box_recurrence`.** The same statement for side lengths, with the hypotheses required only for boxes with positive sides.
- **`conditional_gap_of_boxes`.** The conditional gap with `ym:conj:recurrence` discharged. The inputs are `nonneg`, `estimate`, the ratio-form spectral bound at every scale, and the seed `Δ(n) ≤ 1/7480`. The output is `(λ₂/λ₁)ⁿ ≤ 1/2`, that is `𝗆 ≥ log 2 / n`, at every dyadic scale.

**What is not proved.**
- **The two hypotheses.** Both rest on reflection positivity in all four directions, with every box's `Z` the trace of a positive transfer operator (C23). From there Estimate 1 follows by `purity_bounds`, as noted above, but that identification is not formalised.
- **The seed.** A small defect at bounded physical size is not proved either.

## Labels
- The Lean proofs are machine-checked. The match between the Lean statements and the source is my own reading, listed above.
- Model output is labelled `unreviewed (local model output; not evidence)`. The Lean verdict on it is a check, not a review.

## Resume here
Rerun the proofs with `formal/check.sh --strict`. The prover run: `formal/prover.py generate`, then `judge`, then `report`, writing to a directory outside git. Results go in Step 3 above.
