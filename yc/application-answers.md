# YC Application — Public-Safe Draft Language

**Status:** DRAFT
**Note:** The exact authenticated form, personal/company facts, and character limits are maintained outside this public repo. No text here should be pasted without Ernest's review.

## Company name

TMP *(working name)*

## Short description

AI agents that verify theoretical-physics calculations.

## What will the company make?

TMP is an AI research system for theoretical physics. It generates candidate calculations, routes them to independent symbolic or numerical checks, and records an explicit pass, failure residual, or inconclusive result.

We are starting with scattering amplitudes, where plausible output can still be wrong by a sign, convention, symmetry, or physical limit. The long-term product is a research workspace in which agents can explore aggressively because verification and provenance are built into the loop.

## Why this problem?

Language models can produce convincing physics that is subtly wrong. Computer algebra systems can reliably execute specified calculations, but an expert still has to decide what to try, maintain conventions, run independent checks, and preserve negative results.

TMP connects those two modes: a creative proposer and a deliberately boring verifier with a persistent evidence ledger.

## Founder fit

Ernest has a B.S. in Physics from Caltech, an M.Sc. in Physics from LMU Munich's Theoretical and Mathematical Physics program, and production aerospace software experience across embedded C++, Rust, CUDA, and Python. That combination supports both the domain work and the reliability discipline required to build the product.

## Key insight

Proposal and verification should be separate systems with different incentives. The agent can be creative and cheap to rerun; the verifier should be deterministic wherever possible and must be allowed to reject the claim or return “inconclusive.”

Scattering amplitudes are a useful first wedge because exact identities, recurrence relations, symmetries, and soft/collinear limits make correctness measurable.

## Current progress

- A Python implementation of the closed-form expression from arXiv:2602.12176v2.
- A separate implementation of the specialized recurrence described in the same paper.
- Focused unit/integration suite: 43 tests passed on 2026-07-19.
- Code-path comparison through `n = 10`.
- Standalone checks of explicit low-multiplicity expressions, discrete structure, cyclicity, and a soft theorem.

This is a verified reproduction of a known result. The paper itself gives an all-`n` formula and proof, so our finite-range test is an engineering/reproducibility milestone—not evidence that we extended the physics.

## Competitors and substitutes

- General LLM assistants propose derivations but do not provide dependable ground truth.
- Mathematica, FORM, Cadabra2, SymPy, and theorem provers verify bounded tasks but generally require expert-directed workflows.
- Scientific-agent systems automate broader research tasks; TMP's starting point is physics-specific verifier contracts and evidence persistence.

These systems are inputs and complements, not merely competitors.

## Users and go to market

The first potential users are theoretical-physics researchers and scientific-AI teams with calculations that can be rejected by exact identities or physical constraints. The working plan is to release a reproducible benchmark, recruit design partners directly, and integrate one painful existing calculation per team.

This is a hypothesis until real user conversations support it. Do not claim users or demand in this public draft.

## Business model hypothesis

Paid private team workspaces plus managed model/compute usage, with public local verifiers and benchmarks supporting adoption. Collaboration, private research context, audit logs, and domain integrations are candidate paid features.

No pricing or revenue claim is made here.

## Why now?

Models can now write code and use tools over multi-step workflows, while mature symbolic systems provide deterministic checks for bounded classes of claims. The new opportunity is to connect them into a closed research loop and measure it with domain-specific evals.

## Near-term milestones

1. Ship a repeatable proposer/verifier/evidence-ledger demo.
2. Recruit serious design partners and observe real workflows.
3. Release a small spinor-helicity eval with measured baselines.
4. Explore a bounded research question outside the reproduced baseline, logging failures as well as successes.
5. Convert a repeated user workflow into a paid pilot.

## Risks

- The initial research wedge may not become a product people pay for.
- Verifiers cover bounded claim classes and must expose inconclusive cases.
- Original results are rare; verification and reproducibility must create value before novelty.
- Founder commitment, legal structure, users, and revenue must be answered from private facts, not inferred here.
