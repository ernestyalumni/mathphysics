# physics-tools: bounded exact Rust seed bridge

Early F12/F13 CLI over the existing qft-symbolics Clifford kernel. No general
F09/F10 expression or proof engine, MCP server, arbitrary execution or formal
proof is implemented. This is not complete F13 acceptance.

Build/test entirely offline with cached dependencies:

```sh
export CARGO_TARGET_DIR=/tmp/mathphysics-wave1-tools-target
cargo test --offline
cargo test --offline --release
cargo clippy --offline --all-targets -- -D warnings
printf '%s\n' '{"command":"discover"}' | cargo run --offline --quiet
printf '%s\n' '{"command":"gamma_trace","labels":["new_mu","new_nu"],"domain":"formal_symmetric_h","budget":1}' | cargo run --offline --quiet
```

One JSON stdin object (<=1 MiB), one JSON stdout response. No arguments. Discovery
publishes the schema, native runtime, bounds, conventions and unsupported scope.
Gamma words have 0..12 labels, each an ASCII identifier of 1..64 bytes. All
unknown/duplicate fields, including null-first duplicates, and trailing JSON
are rejected. Domain must be exactly `formal_symmetric_h`; budget is a positive
integer <=10395 covering the oriented pairing count `(n-1)!!` (empty/odd=1).
Budget counts pairing leaves, not time, memory or instructions. Byte-bounded input
does not impose a wall-clock timeout on a caller that keeps stdin open.

`computed` has evidence_kind `exact_algebra`, formal_proof false, input labels,
domain, assumptions, kernel identity and exact polynomial terms. Coefficients
are decimal integer strings, avoiding JSON consumer integer precision loss.
Factors are opaque length-prefixed symmetric metric variables such as
`h(6:new_mu;6:new_nu)`. Term order follows the kernel canonical order; zero is
an empty term array. Display is also preserved. Serialization relies only on
the current trace's integral real display grammar and checks that invariant;
an unexpected format becomes inconclusive. It does not parse arbitrary expressions.

Conventions: `{gamma(a),gamma(b)}=2h(a,b)`; Srednicki `h=-eta`, eta=(-,+,+,+).
Four-component Dirac trace, trace(identity)=4; no gamma5. Every label is free,
never implicitly summed. Repetition identifies the same slot, leaving `h(a,a)`
formal: trace(gamma(a)gamma(a))=4h(a,a), with no normalization or dimension sum.
The domain is a symbolic declaration, not a certificate checking physical premises.

Exit codes: 0 discovered/computed, 2 invalid_input, 3 unsupported, 4 over_budget,
5 inconclusive. Errors have distinct `error.code` values. Kernel panic (including
checked-i128 overflow and polynomial budget failure) is caught around computation
and serialization; the default hook may emit stderr diagnostics. There is no
floating-point fallback or PASS. OS termination, allocation abort, panic=abort
build overrides and broken stdout cannot be guaranteed a structured response.

The unchanged kernel has polynomial bounds 8192 terms and degree64. Twelve
independent distinct labels generate 10395 pairing monomials and therefore yield
`inconclusive` at its term bound, despite an adequate declared pairing budget.
Twelve repeated labels can compute. Both behaviors are tested. No general
unbounded arithmetic, tensor canonicalizer, implicit contraction or proof is implied.

Implementation provenance and handoff: [docs/IMPLEMENTATION.md](docs/IMPLEMENTATION.md).
Approved intent: [docs/PRODUCT_INTENT.md](docs/PRODUCT_INTENT.md). The versioned
prompt/config graph owns src/lib.rs and src/main.rs only. Tests are independent
and retained across future generation. This code was manually implemented by
an agent; automated PDD apply/sync/regeneration/baseline success is not claimed.
