# DISCOVERY-LOOPS.md — Registered AI Discovery Loops

Last updated: 2026-06-16

This file is the registry of TMPAI's AI loops — the workhorses of the discovery engine. Each loop has: a goal, a topology, prompts, a verifier, persistence, and a run log.

## Universal rules

1. **Every loop has a verifier.** No verifier step → not a loop. Brainstorms go elsewhere (e.g. session notes).
2. **Symbolic equality only.** Verifier output is exact-equality after canonicalization, or a documented numerical limit with explicit tolerance + a physical reason.
3. **Persist every run** to `discoveries/YYYY-MM-DD-<slug>.md` regardless of pass / fail. Failure logs are data.
4. **No publishing from a loop.** Loops produce drafts; Ernest approves anything that leaves the repo.
5. **Pin models.** Loops record model id + provider + prompt version + verifier version + git sha at run time.

## Standard topology

```
Proposer (LLM)
  ↓ conjecture / candidate identity / candidate amplitude / candidate derivation
Critic (LLM, different model or different prompt)
  ↓ critique, request reformulation
Verifier (symbolic: Cadabra2 / SymPy / FORM)
  ↓ pass / fail / inconclusive with explicit residual
Logger
  ↓ writes discoveries/YYYY-MM-DD-<slug>.md
```

Single-agent variant: proposer + verifier only. Use when critic-cost outweighs the conjecture cost.

## Run-log file format (`discoveries/YYYY-MM-DD-<slug>.md`)

```markdown
# <slug>

**Loop:** <loop-id from this registry>
**Date:** 2026-06-16
**Models:** proposer=<provider/model>, critic=<provider/model>
**Prompt version:** <hash or version>
**Verifier version:** <Cadabra2 docker tag + verifier git sha>
**Git sha:** <repo HEAD at run time>
**Branch:** <branch name>

## Conjecture
<the conjecture in math notation>

## Critic notes
<critic output summary>

## Verifier output
- Status: PASS | FAIL | INCONCLUSIVE
- Residual: <expression after canonicalization>
- Notes: <anything the verifier could not decide>

## Verdict
<one sentence: is this a known result (cite source), a re-derivation, or genuinely new?>

## Follow-up
- [ ] <next action, if any>
```

---

## Registered loops

### LOOP-001 — `mhv-ward-identity-checker` (planned)

- **Status:** not yet implemented. Skeleton only.
- **Goal:** given an n-point MHV amplitude expression, propose Ward-identity checks (gauge invariance under cyclic shifts of momenta) and verify them in Cadabra2.
- **Topology:** single-agent + verifier (proposer + Cadabra2).
- **Proposer prompt skeleton:** `prompts/loop-001-proposer.md` (to be created).
- **Verifier:** `Cadabra2/MHV/ward_check.py` (to be created).
- **Stopping criterion:** verifier passes OR 3 reformulations exhausted.
- **Persistence:** `discoveries/YYYY-MM-DD-mhv-ward-<n>.md`.
- **First milestone:** verify n=4 Parke-Taylor under cyclic shifts.

### LOOP-002 — `spinor-identity-conjecturer` (planned)

- **Status:** not yet implemented.
- **Goal:** generate candidate spinor-helicity identities (angle / square brackets, Schouten, Fierz) and verify them.
- **Topology:** proposer + critic + verifier.
  - Proposer: small / cheap model (high-throughput conjecturing).
  - Critic: stronger model that rejects obviously-trivial / known-equivalent restatements.
  - Verifier: Cadabra2 + a curated identity database.
- **Open question:** how to encode "obviously trivial" without hard-coding the database the loop is supposed to extend?
- **First milestone:** rediscover Schouten identity from a blank prompt; log it as INCONCLUSIVE-NEW until verified against the database.

### LOOP-003 — `bcfw-recursion-runner` (planned)

- **Status:** not yet implemented.
- **Goal:** given a tree-level amplitude target, drive BCFW recursion symbolically and verify the result matches Parke-Taylor for the MHV case.
- **Topology:** single-agent + verifier. Mostly mechanical, but useful as a regression test for Cadabra2 infrastructure.
- **First milestone:** automate the n=5 MHV case end-to-end.

### LOOP-004 — `dg-identity-checker` (planned, longer horizon)

- **Status:** not yet implemented.
- **Goal:** verify differential-geometry identities (Bianchi, Ricci, structure equations, pullback / pushforward) symbolically.
- **Topology:** proposer + verifier (SymPy / Cadabra2 depending on identity).
- **First milestone:** Bianchi identity in a fixed connection.

---

## Anti-pattern catalog

These mistakes are easy to make in a hurry. Avoid them.

- **"Verified" without canonicalization.** Two expressions that look different may be equal — and two that look equal may differ by a sign in mostly-plus vs mostly-minus. Always canonicalize before comparing.
- **"New result" without literature check.** Almost everything has been derived before. Default verdict: re-derivation. Only escalate after a real literature search.
- **Loop without persistence.** If a run isn't logged to `discoveries/`, it didn't happen.
- **Prompt drift.** Update prompt version when changed. Don't silently mutate a registered loop.
- **Verifier silently downgraded.** If you weaken the verifier (looser tolerance, skipped step), version-bump it and note why.

## How to add a new loop

1. Add a new section `### LOOP-NNN — <slug>` to this file.
2. Specify topology, prompts location, verifier location, stopping criterion, persistence path.
3. Implement prompts under `prompts/loop-NNN-*.md`.
4. Implement verifier under `Cadabra2/<topic>/` (or `evals/<name>/verifier.py` if it's an eval).
5. Add a `_test/` with at least one known-answer case before the loop is allowed to run on novel input.
6. First run goes to `discoveries/YYYY-MM-DD-<slug>-first.md` regardless of outcome.

## Cross-references

- `../TMP-MISSION.md` — pillar 2 (AI loops / discovery engine).
- `EVAL-AND-TOOLING.md` — verifier infrastructure overlaps with eval harnesses.
- `RESEARCH-AGENDA.md` — which threads each loop serves.
