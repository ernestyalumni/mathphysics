# Typed tensor kernel: implemented boundary

Implemented 2026-10-16 as an independent mathematical implementation.
It is not a translation of Cadabra/xPerm and does not claim their scalability.

## Mathematical contract

- Three distinct index spaces: Lorentz, undotted, dotted.
- Lorentz slots represent all symmetric-metric-related variance views of the
  **same tensor**. A slot permutation moves the complete variance-bearing
  index. Spinor slots have fixed space and variance.
- Each index in a raw monomial occurs once (free), or exactly twice with
  opposite variance (dummy). Same-variance repetitions and triple/fourfold
  occurrences are errors, even when the coefficient is zero.
- Free signatures include space, name and variance. Addition requires identical
  signatures; free names are not alpha-renamed.
- A zero tensor retains its free signature. A covector zero is not silently
  accepted as a scalar zero.
- All tensor factors in this layer **commute**. Grassmann fields and gamma
  chains must not be represented as ordinary tensor factors here.
- Lorentz metric and inverse are constant and symmetric. Mixed metric is the
  identity delta. No derivative/dependency machinery is implied.
- Deltas are typed identities, including spinor identities. Spinor epsilon may
  be declared as a generic antisymmetric tensor, but automatic epsilon
  raising/lowering is explicitly unsupported.

## What canonicalisation does

1. Validate declarations, slots, occurrences and every summand's free signature.
2. Eliminate constant metric/delta contractions, preserving that signature.
   Closed deltas give the appropriate dimension: Lorentz 4 by default (or
   explicit formal d), undotted/dotted 2.
3. Enumerate the joint finite orbit of:
   - declared signed slot symmetries;
   - permutations of commuting factors;
   - renamings of dummy binders;
   - simultaneous up/down flips at the two ends of each **Lorentz** dummy.
4. Alpha-normalise dummy names in each candidate traversal. Internal bound IDs
   have a different Rust type from user free labels, preventing name collision.
5. Choose the least representative. Reaching it with both signs means exact
   zero, e.g. antisymmetric K_ab S^ab or K_ab V^a V^b.
6. Collect coefficient polynomials over Q(i).

This is mono-term symmetry reduction, not general tensor equivalence.
Multi-term Bianchi/Schouten identities and dimension-specific antisymmetry
identities are not inferred.

Canonical representatives remain a sum of tensor monomials; scalar coefficient
variables do not acquire tensor semantics merely because they are named "mu".

## Hygienic multiplication

Separately built contractions have local binders:

    (A^mu B_mu) * (C^mu D_mu)
       = A^mu B_mu C^nu D_nu.

The RHS operand's internal bound IDs are freshened before combination.
Free indices may contract across operands when they have the same name/space
and opposite variance. A raw flat A^mu B_mu C^mu D_mu, on the other hand,
is rejected; the kernel never guesses a pairing for a fourfold occurrence.

Contexts are immutable snapshots inside validated expressions. Combining
different snapshots is rejected rather than silently changing the algebra.

## Resource and error semantics

Tensor operations use Result with:
- Invalid: malformed declarations/indices.
- SignatureMismatch or ContextMismatch.
- Unsupported: algebra outside the implemented domain.
- Budget: no complete representative established.

Defaults: at most eight commuting factors, twenty slots per monomial,
rank eight per declared tensor, 100, 000 work steps. Group closure,
variant preparation and complete orbit search consume explicit budgets.
A high-symmetry orbit is charged **before** cloning all variants.
These are work bounds, not an OS sandbox or a byte-exact memory guarantee.

Coefficient operations retain the crate's checked-i128 panic boundary;
arbitrary precision and uniformly fallible arithmetic are not implemented.
A failure/partial orbit must never be used as an equality proof.

## Use

See src/tensor_cases.rs for complete small Rust constructions.
From the crate directory:

```bash
cargo run --offline -- tensor-cases
cargo test --offline --test tensor
```

The JSON command computes twelve named residual verdicts, including three
known-nonzero controls. validation/tensor_oracle.py applies actual Cadabra
transformations to the same equations and compares the zero/nonzero outcomes.

Cadabra adapter details are explicit: fixed Lorentz index declarations,
symmetry declarations for applicable variance views, inverse-metric identity
g_a^b=delta^b_a, product sorting, dummy renaming, and collection.
A mixed metric does not disappear merely because it was displayed.

There are twenty Rust tests covering these semantics, including invalid-zero
inputs, capture avoidance, signed stabilizers, formal dimension, contexts,
epsilon limitations, canonical fixed points and exhausted high-symmetry work.
