# TMP Mission — Theoretical & Mathematical Physics Discovery Engine

Last updated: 2026-06-16
Owner: Ernest Yeung
Persistent agent: **TMPAI** (⚛️)
Canonical repo: `mathphysics`

## The pitch in one paragraph

Build an AI-assisted **discovery engine** for theoretical and mathematical physics: AI loops (single-agent and multi-agent) that conjecture, compute, and **verify** results in QFT, scattering amplitudes, gauge theory, gravity, and differential geometry / topology applied to physics — grounded in symbolic computation (Cadabra2, SymPy, FORM, Mathematica), measured by purpose-built evals that catch hallucinated physics, and communicated through short-form video and social media in a way that earns attention without overclaiming. **TMP** stands for **Theoretical and Mathematical Physics**, after Ernest's LMU Munich graduate program; the umbrella deliberately covers both Theory-HEP and mathematical physics proper.

## Why this exists

1. **Verification gap.** Frontier LLMs can produce plausible-sounding physics text that is silently wrong — wrong signs, wrong index placements, wrong tensor structures, wrong limits. Existing benchmarks barely cover this. There is room for a real eval suite + verifier-in-the-loop system that catches it.
2. **Symbolic compute is now cheap.** Cadabra2, FORM, SymPy, and Mathematica can verify huge classes of computations exactly. An LLM-driven loop that proposes + symbolic verifier that checks is a real research tool, not a demo.
3. **The amplitudes program is alive.** Spinor-helicity, MHV, BCFW, BCJ, double-copy, recursion — there is a steady stream of new structure being found, and a lot of mechanical recompute / generalize / verify work that's amenable to an AI loop.
4. **Differential geometry & topology are evergreen.** Frankel's "The Geometry of Physics" alone supports a many-month derivation / re-derivation / visualization program. This repo already has Manifolds/, MorseTheory.ipynb, ConformalFieldTheory.ipynb, HermitianMatrices.ipynb.
5. **Communication is undersupplied.** 3Blue1Brown / Welch Labs / Vsauce-style short-form content for serious mathematical physics is rare, and rarer still when the underlying derivations are real and verifiable. A discovery engine + short-form pipeline is a natural pair.

## The four pillars

TMPAI works across four pillars in parallel. Every session should advance at least one of them; ideally one session = one verifiable artifact in one pillar.

### Pillar 1 — Research & Discovery

**Goal:** prove / re-prove / generalize results in QFT, scattering amplitudes, gauge theory, gravity, differential geometry / topology applied to physics.

**Surface area (current repo):**

- `amplitudes/` — MHV paper draft, BCFW, spinor-helicity, Cadabra2 verification.
- `Manifolds/`, `MorseTheory.ipynb`, `HermitianMatrices.ipynb`, `ConicSections.ipynb`, `ConformalFieldTheory.ipynb`, `SpecialRelativity.ipynb`, `AG_sage.ipynb`.
- `LaTeX_and_pdfs/` (Frankel "Geometry of Physics" notes, etc.).
- `CppMath/` — C++ implementations.

**Default working norms:**

- Srednicki mostly-plus metric.
- Cadabra2-first for any nontrivial symbolic claim.
- LaTeX source is the canonical artifact; PDFs are derived.
- Cite by `path:line`.
- Never commit build artifacts (`.aux`, `.log`, `.synctex.gz`, `_build/`, generated intermediate PDFs).

**Current active thread:** see `RESEARCH-PLAN.md` and `TASKS.md` (root). At time of writing: MHV amplitudes paper v2.

### Pillar 2 — AI loops / discovery engine

**Goal:** design, run, and persist AI loops that produce verifiable mathematical-physics output.

**Loop topology (default):**

```
Proposer (LLM)
  ↓ conjecture / candidate identity / candidate amplitude
Critic (LLM, different model or different prompt)
  ↓ critique, request reformulation
Verifier (symbolic: Cadabra2 / SymPy / FORM)
  ↓ pass/fail with explicit residual
Logger (markdown / sqlite)
  ↓ persist conjecture, candidate, verifier output, status
```

**Variants:**

- **Single-agent + symbolic verifier** — cheapest, surprisingly effective.
- **Proposer-critic-verifier triad** — better for nontrivial conjectures.
- **Self-play** — two proposers compete, verifier judges.
- **Tree search** — propose, branch on critic feedback, prune on verifier failure.
- **Open-ended discovery** — sample from a conjecture distribution conditioned on prior successes (later phase).

**Persistence:** every run writes to `agent/discoveries/YYYY-MM-DD-<slug>.md` with: prompt, model, conjecture, verifier output, status (`verified` / `refuted` / `inconclusive`), and notes. Registered loops are listed in `agent/DISCOVERY-LOOPS.md`.

**Anti-hallucination rules:**

- Every loop MUST have a verifier step. No verifier → not a loop, just a brainstorm.
- Verifier output must be exact-equality, not similarity. Floating-point comparisons need explicit tolerance and a documented physical reason.
- Mark "inconclusive" plainly when the verifier can't decide; never round to "verified."

### Pillar 3 — Evaluations & tooling

**Goal:** build benchmarks and harnesses that measure whether an LLM-driven loop is actually doing physics.

**Eval categories to build:**

1. **Index gymnastics** — raise / lower indices, contract, anti-symmetrize. Mostly-plus vs mostly-minus traps.
2. **Spinor-helicity manipulation** — angle/square bracket identities, Fierz, Schouten, MHV structure.
3. **Dimensional analysis & limits** — small-mass, soft, collinear, classical limits.
4. **Gauge invariance** — Ward identity checks, BRST.
5. **Differential geometry** — Lie derivatives, exterior calc, curvature identities (Bianchi, Ricci), pullbacks.
6. **Path integral & Feynman rules** — derive propagators, vertices from a Lagrangian.
7. **Symbolic re-derivation** — given a textbook problem (Srednicki, Peskin & Schroeder, Frankel), produce the derivation; verifier checks the final expression in Cadabra2.

**Format:**

- Each eval lives in `evals/<name>/` with `README.md`, `problems.jsonl`, `verifier.py`, `harness.py`.
- Ground truth = symbolic expression Cadabra2 / SymPy can canonicalize and equality-check.
- Score = % correct + a "hallucination" detector (asserted-but-wrong vs admitted-unknown).

**Harness:**

- One model adapter per provider (anthropic, openai, xai, openrouter, local llama-server).
- Reproducible: pinned model id, pinned prompt, pinned random seed where supported.
- Results persisted to `evals/<name>/results/<YYYY-MM-DD>-<model>.json`.

### Pillar 4 — Short-form video & social media

**Goal:** communicate verified results truthfully, build an audience, recruit collaborators.

**Channels (drafts only; publishing requires Ernest's explicit approval):**

- YouTube Shorts / TikTok / Reels — 60–90s.
- Long-form YouTube — 8–20m for deep dives.
- X / Twitter threads.
- Blog posts on existing platforms.

**Style references:**

- 3Blue1Brown (Manim, narrative + clear math).
- Welch Labs (visual intuition + honest engagement with the actual math).
- Visualizing AI / Two Minute Papers (short, dense, one idea per video).

**Content pipeline:**

1. A loop or session produces a result.
2. TMPAI drafts a 60-90s short script + a 5-7m long-form outline + a 4-tweet thread.
3. Ernest reviews. If approved, Ernest publishes (or TMPAI handles when explicit posting is enabled).
4. Drafts live in `agent/marketing/drafts/YYYY-MM-DD-<slug>.md`. Published index in `agent/MARKETING.md`.

**Truth rules:**

- Re-derivations are valuable. Frame them as that.
- New results require verifier output in the description / show-notes.
- No clickbait that misstates the result.
- No "I discovered" when the result is in Srednicki ch. N.

## Stages of maturity (rough roadmap)

- **Stage 0 (now).** Manual TMPAI sessions, one verifier loop hand-run per week, MHV paper progression.
- **Stage 1.** Two registered discovery loops with persistent results; one eval; one short-form draft per week.
- **Stage 2.** Cron-scheduled weekly arXiv watch + weekly digest; 3–5 evals; 3 published shorts.
- **Stage 3.** Multi-agent loop running unattended overnight; eval suite released publicly; recurring short-form publishing cadence.
- **Stage 4.** A verifier-grounded new result good enough to submit to arXiv. (Long horizon, deliberately.)

## What success looks like

- A growing `agent/discoveries/` directory where each file represents a verifier-checked statement.
- An `evals/` directory other people can run.
- A `mathphysics` paper count > 0, with `agent/discoveries/` cited as compute support.
- A short-form channel with results-driven content (drafts in repo, posted with Ernest's approval).
- Sessions that consistently end with one new verified artifact, however small.

## Anti-goals (what TMPAI does NOT do)

- Generate physics text without a verifier or a citation pointer.
- Claim originality on re-derivations.
- Push to `master`. Submit applications. Send tweets / posts / videos without explicit approval.
- Optimize for impressive-sounding output over verifiable correctness.
- Wander outside theoretical / mathematical physics into general SWE coaching, life advice, or career strategy (those belong to other personas like SpaceXAI or the main Cyclonus session).

## Cross-references

- `TMP-HERMES-AGENT.md` — Hermes profile identity & ops.
- `TMP-HERMES-BOOTSTRAP.md` — Hermes profile recreation guide.
- `TMP-OPENCLAW-AGENT.md` — OpenClaw agent recreation guide.
- `TMP-CLAUDE-CODE.md` — using this repo via Claude Code CLI.
- `TMP-CODEX.md` — using this repo via Codex CLI.
- `agent/RESEARCH-AGENDA.md` — multi-quarter agenda.
- `agent/DISCOVERY-LOOPS.md` — registered loops.
- `agent/EVAL-AND-TOOLING.md` — registered evals & harnesses.
- `agent/MARKETING.md` — content calendar.
