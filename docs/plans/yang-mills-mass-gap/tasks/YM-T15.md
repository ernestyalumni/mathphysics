# YM-T15 — Signed group words and trace semantics for Wilson observables

TierT2–T3 · Symbolic tools · First bounded packet from Ernest's Wilson-measure study request · Review: independent T2–T3.

Read notes/wilson-symbolics-plan.md, the reusable Rust library story and PDD instructions, Rust/qft-symbolics's ordered/coefficient/budget units and their prompts, then this packet's tracker thread. Claim there and work in your own ym/* branch/worktree. The source request and draft contract are context, not evidence of successful generation.

## Do only this first unit

Design and implement a bounded reusable signed group-word and trace kernel at the owning qft-symbolics boundary. Inputs are caller-supplied typed generator IDs, orientations, representation dimension/domain assumptions and declared trace kind. Preserve multiplication order, inverse-of-product order reversal, safe adjacent cancellations and cyclic trace equivalence. Separate complex trace from normalized real trace; retain a conservative equality/unknown distinction. Reuse existing exact scalar/budget infrastructure without duplicating ordered algebra. A lattice-specific fixture may supply the plaquette word and gauge factors; production lattice adapters, Haar integrals, Monte Carlo and CUDA are out of this packet.

Before mutation run read-only intent planning against the actual owning scope and read installed help. Use the established bounded prompt-first workflow; do not launch the root19-step provider workflow, change authentication, or claim a generated/synchronized unit without evidence. Promote the reviewed contract to a properly mapped crate prompt and preserve characterization of existing units. Ask only for unresolved mathematical meaning, not filenames or commands.

## Independent acceptance

- Words with fresh generator names reduce adjacent inverse pairs; word inversion reverses order and orientation; original inputs remain unchanged.
- Cyclic trace equality works, including justified cyclic inverse cancellation, without sorting noncommuting factors.
- Exact Pauli3-factor oracle distinguishes Tr(ABC)=2 from Tr(ACB)=-2. Two-factor traces intentionally remain equal under cyclic rotation.
- Exact SU3 diagonal oracle distinguishes trace from trace of inverse; normalized real traces agree under inversion when unitarity is declared.
- Normalized real trace of identity is1; raw complex trace of identity is positive representation dimension m. Invalid dimensions/domains, excessive words/coefficient growth and budget exhaustion fail explicitly.
- Gauge-transformed four-link fixture reduces to endpoint conjugation; closed-loop trace is invariant. No generic open-path invariant claim.
- Existing affected algebra/CLI tests still pass; no hardcoded equation IDs/answers in production.

## Deliverable

One source prompt/API implementation with independent tests, a maintained small example, and notes/YM-T15-group-trace-kernel.md recording owner/mapping, commands, limits and what remains unsupported. A new dependency or semantic widening requires explicit justification. No mass-gap or general-CAS-parity claim. Submit for independent review before done.
