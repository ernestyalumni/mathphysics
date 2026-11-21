# Intent: Automatic Lean statement checker in Rust

<!-- pdd-intent-id: automatic-lean-statement-checker-in-rust-3c8e1031 -->
<!-- pdd-intent-sha256: 3c8e1031cd1769e1cfc7b20c52516209828acc0fded863b8344b5c74da2139cc -->

## Record

- Intent ID: `automatic-lean-statement-checker-in-rust-3c8e1031`
- Kind: `add`
- Supersedes: none
- Approval ID: `automatic-lean-statement-checker-in-rust-3c8e1031`
- Source kind: `inline`
- Source reference: `not applicable`
- Request SHA-256: `3c8e1031cd1769e1cfc7b20c52516209828acc0fded863b8344b5c74da2139cc`
- Project scope: `repository`
- Adoption scenario: `existing_pdd_change`

- Technology: `Rust`

## Original Request

> Ernest, 2026-11-21, after the Goedel-Formalizer result (4 of 9 compiling statements were false):
> "what about writing enough symbolic computation tools to help out? tbh we've gotta automate this, really can't rely upon me"
> "ok go ahead" (approving the agent's plan: an automatic checker that a Lean statement says what the paper says, by two independent translations - paper claim to an executable predicate, Lean statement mechanically from Lean's elaborated term - compared on random and edge-case instances; plus vacuity, counterexample and hypothesis-necessity checks; validated on the 28 graded formalizer statements and mutants of the 18 reference statements, pre-registered.)
> "also anytime we use sympy / numpy we need to start writing a port of those tools into Rust / CUDA C++."
> "i really want to banish the use of sympy/numpy because we have a full Rust port already"

## Must Stay Unchanged

- None stated.

## Examples

- None stated.

## Candidate Product Areas

- Bounded reusable arithmetic foundation; manual implementation
- Loopback read-only file/static gateway and existing physics-tools bridge
- Approved reusable research batches; bounded manual implementation
- Approved reusable research batches; bounded manual implementation
- Approved reusable research batches; bounded manual implementation
