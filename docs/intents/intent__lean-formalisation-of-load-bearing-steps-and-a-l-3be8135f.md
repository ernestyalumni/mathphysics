# Intent: Lean formalisation of load-bearing steps and a local-prover test

<!-- pdd-intent-id: lean-formalisation-of-load-bearing-steps-and-a-l-3be8135f -->
<!-- pdd-intent-sha256: 3be8135f77a18262b8297122194e553e91c5383eeb962ae960a5f2919c96e6b4 -->

## Record

- Intent ID: `lean-formalisation-of-load-bearing-steps-and-a-l-3be8135f`
- Kind: `add`
- Supersedes: none
- Approval ID: `lean-formalisation-of-load-bearing-steps-and-a-l-3be8135f`
- Source kind: `file`
- Source reference: `docs/plans/yang-mills-mass-gap/request-2026-11-20-lean-formalisation.txt`
- Request SHA-256: `3be8135f77a18262b8297122194e553e91c5383eeb962ae960a5f2919c96e6b4`
- Project scope: `repository`
- Adoption scenario: `existing_pdd_change`

- Technology: `docker`

## Original Request

> Ernest, 2026-11-20: "Ok let's give 1,2,3 a try please", approving these three steps proposed by the agent:
>
> 1. Pick 3-5 steps the Yang-Mills mass-gap argument actually depends on (in documents/mono-notes/topics/yang-mills-lattice.tex).
> 2. Formalise them in Lean (Lean 4 + Mathlib, run in Docker, never on the host).
> 3. Test whether a local Lean prover model on the RTX 3060 can fill in any of the proofs, with Lean as the only judge.
>
> Context he accepted: a compiled Lean proof checks the Lean statement, not the paper, so each statement must be read against the source by hand; expect a few verified steps, not broad coverage; the analytic core of the mass-gap argument is out of reach for now.

## Must Stay Unchanged

- Formalise them in Lean (Lean 4 + Mathlib, run in Docker, never on the host).
- Context he accepted: a compiled Lean proof checks the Lean statement, not the paper, so each statement must be read against the source by hand; expect a few verified steps, not broad coverage; the analytic core of the mass-gap argument is out of reach for now.

## Examples

- None stated.

## Candidate Product Areas

- Collector
- Wilson Group Trace Contract
- Ym Tracker Bash
- Exact density of states of Z2 lattice gauge theory on small tori, with twists
- Exact integer check of the purity bounds for a nonnegative transfer operator
