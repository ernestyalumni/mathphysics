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

## Labels
- The Lean proofs are machine-checked. The match between the Lean statements and the source is my own reading, listed above.
- Model output is labelled `unreviewed (local model output; not evidence)`. The Lean verdict on it is a check, not a review.

## Resume here
Rerun the proofs with `formal/check.sh --strict`. The prover run: `formal/prover.py generate`, then `judge`, then `report`, writing to a directory outside git. Results go in Step 3 above.
