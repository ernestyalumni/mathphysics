# EVAL-AND-TOOLING.md — Registered Evals & Harnesses

Last updated: 2026-06-16

This file registers TMPAI's evaluation benchmarks and harnesses. Evals measure whether an LLM-driven loop is actually doing physics, not just producing plausible physics text.

## Why evals matter for TMP

LLMs are excellent at producing physics-shaped output. They are not yet reliably correct on:

- Index placement (raised vs lowered, anti-symmetrization).
- Metric-sign conventions (mostly-plus vs mostly-minus — Srednicki vs Peskin).
- Spinor-helicity bookkeeping (angle vs square brackets, contractions, weights).
- Gauge-invariance edge cases (Ward identities under cyclic momentum shifts).
- Small / soft / collinear limits.
- Differential-geometry identities (sign of curvature tensor, Bianchi structure).

A purpose-built eval that catches these failure modes is the foundation of the discovery engine. Without it, the loops in `DISCOVERY-LOOPS.md` produce noise.

## Universal rules

1. **Symbolic ground truth.** Where possible, the answer is a canonicalized symbolic expression Cadabra2 / SymPy can compare exactly. Free-form text answers are last resort.
2. **Pinned models.** Each eval run records model id + provider + prompt version + harness version.
3. **Hallucination detector.** Every eval scores both *accuracy* (right answer) and *honesty* (admitted-unknown vs asserted-wrong). An LLM that says "I don't know" beats one that confidently states the wrong sign.
4. **Reproducibility.** Pinned seeds where the provider supports it. Deterministic decoder settings recorded.
5. **No leaking the test set into the loop training/prompt pipeline.** Hold out at least 20% per eval as a private test set; never expose it to a proposer model.

## Eval directory layout

```
evals/
  <name>/
    README.md           — what this eval covers, scoring rubric, baseline scores
    problems.jsonl      — one problem per line, with structured ground truth
    verifier.py         — symbolic equality check
    harness.py          — runs a pinned model, scores responses, writes results
    results/
      YYYY-MM-DD-<provider>-<model>.json
    _test/              — known-answer sanity tests for the verifier itself
```

## Scoring rubric (default)

Each problem yields one row:

```json
{
  "problem_id": "string",
  "model_answer": "string",
  "canonical_form": "string after symbolic canonicalization, or null",
  "verifier_status": "PASS | FAIL | INCONCLUSIVE",
  "admitted_unknown": false,
  "latency_ms": 1234,
  "tokens_in": 567,
  "tokens_out": 89
}
```

Aggregate scores per eval run:

- `accuracy` = PASS / total.
- `confident_wrong` = FAIL where the model did NOT admit uncertainty. **This is the headline metric.**
- `honest_unknown` = `admitted_unknown` count. (Higher is fine; better than `confident_wrong`.)
- `inconclusive_rate` = INCONCLUSIVE / total. Investigate if > 10%.

## Registered evals

### EVAL-001 — `index-gymnastics` (planned)

- **Skill:** raise / lower indices, contract, anti-symmetrize, work with metric of given signature.
- **Ground truth:** canonical tensor expression after `Cadabra2.canonicalise_indices`.
- **Failure modes targeted:** Srednicki-vs-Peskin sign flip, summation-convention mistakes, dummy-index collisions.
- **Problem count target:** 30 (20 public, 10 held out).
- **Baseline to publish:** at minimum claude-sonnet-4-6, claude-opus-4-7, gpt-5.5, grok-4.3, kimi-k2.6.

### EVAL-002 — `spinor-helicity-basics` (planned)

- **Skill:** Schouten identity, momentum conservation in `⟨ij⟩[ji]`, Fierz rearrangements, weight tracking.
- **Ground truth:** Cadabra2 simplification of LHS − RHS down to zero.
- **Failure modes targeted:** confusing `⟨ij⟩` vs `[ij]`, dropping prefactors, sign on antiholomorphic side.
- **Problem count target:** 25.

### EVAL-003 — `mhv-amplitudes-textbook` (planned)

- **Skill:** state the Parke-Taylor amplitude, derive single-minus vanishing, compute 5-pt MHV from 4-pt via BCFW shift.
- **Ground truth:** explicit symbolic expressions, normalised.
- **Problem count target:** 15.

### EVAL-004 — `differential-geometry-identities` (planned)

- **Skill:** Bianchi (1st & 2nd), Lie-derivative formula on forms, Cartan structure equations, exterior derivative chain.
- **Ground truth:** SymPy-checkable equalities in given coordinate charts, plus Cadabra2 abstract checks.
- **Problem count target:** 25.

### EVAL-005 — `srednicki-redo` (planned)

- **Skill:** given a Srednicki problem statement, produce the derivation and final symbolic answer.
- **Ground truth:** exists in the Srednicki solutions (or already-verified `.tex` in this repo).
- **Problem count target:** 20, drawn from chapters currently covered in `amplitudes/`.

## Harness

`harnesses/` (to be created) holds:

- `model_adapter.py` — uniform `call(prompt, model_id, provider) -> response` over anthropic, openai, xai, openrouter, local llama-server.
- `runner.py` — given an eval directory and a model id, runs all problems, writes `results/`.
- `judge.py` — symbolic comparator; calls `evals/<name>/verifier.py` per problem.
- `report.py` — produces a markdown summary of a run (PASS/FAIL counts, confident-wrong examples).

Use the same harness for eval scoring AND for verifier-in-the-loop pipelines in `DISCOVERY-LOOPS.md`. Keep them in sync.

## Publication intent

When at least 3 of the planned evals are stable with public baselines:

1. Write a short paper / blog post: "Symbolic Evals for LLM-Driven Theoretical Physics: First Results."
2. Release `evals/` under a permissive license.
3. Publish baselines on HuggingFace eval-style leaderboard.
4. Submit a short to YouTube explaining one failure mode visually ("Here's why GPT-5 gets MHV signs wrong half the time").

All publication requires Ernest's explicit approval. Drafts live in `marketing/drafts/` per `MARKETING.md`.

## Anti-pattern catalog

- **Eval grades model output instead of verifier output.** Always canonicalize before scoring.
- **Eval grows around models' strengths.** Curate problems against known failure modes, not problems where everyone scores 90%+.
- **"Inconclusive" hidden as "PASS".** Surface it.
- **Held-out set leaks into prompts.** Strict separation. Anyone touching prompts cannot see the held-out problems.

## Cross-references

- `../TMP-MISSION.md` — pillar 3.
- `DISCOVERY-LOOPS.md` — verifiers used here are the same as loop verifiers; reuse code.
- `RESEARCH-AGENDA.md` — evals serve as regression tests for the active research thread.
