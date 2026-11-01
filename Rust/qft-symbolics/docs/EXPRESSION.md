# Typed structural expressions

`expression::Expr` exposes `Constant(BigC)`, `Symbol(String)`, `Sum(Vec<Expr>)`,
`Product(Vec<Expr>)`, `OrderedProduct(Vec<Expr>)`, and
`Quotient(Box<Expr>, Box<Expr>)`. Additional variants are `Operator(String)`,
`Tensor(TensorAtom)` and `Contract { binder, space, body }`.
`validate(&mut Budget) -> Result<Signature>` and `canonical(&mut Budget) ->
Result<String>` always validate raw public fields. `Signature` contains `domain`
and sorted `free_indices`. Replay can use `nonzero_obligations(&mut Budget)`
to obtain sorted, deduplicated canonical `nonzero(...)` predicates.

Constants are exact Gaussian rationals. Scalar symbols commute. Product accepts
only scalars, including fully contracted tensors. OrderedProduct supports tensor
outer products or operators with scalar factors, preserving factor order; mixing
operators and tensors is rejected. Sum requires identical domains and free-index
signatures. Empty Sum and Product represent structural scalar zero and one;
empty OrderedProduct is a scalar identity. Structural nodes are retained: there
is no distribution, arithmetic folding, associativity flattening or cancellation.

Tensor atoms declare an ordered list of typed slots and indices; these must
match exactly. Space and Variance are reused from the existing tensor module.
No tensor symmetry, metric contraction or existing tensor canonicalizer is
reimplemented. Use the existing tensor engine for those algebraic operations.
Repeated free names are rejected, even across spaces. Contraction is explicit:
`IndexName::Bound(name)` must resolve to a lexical binder of matching space,
with two opposite-variance occurrences in every additive branch. Free names
are a different variant and cannot be captured. Empty/shadowed/dangling binders
and branch-signature mismatches fail. Independent contractions may reuse names.

Quotient denominators must be scalar. A literal zero denominator fails; symbolic
nonzero conditions are obligations, not proven facts. Every quotient is retained,
including zero numerators; all descendants are validated before serialization.
No rewrite may infer validity solely because a product contains zero.

Canonical strings start with `qft-expr.v1:` and use structural tags C/S/A/M/P/Q/O/T/B, length-prefixed UTF-8
names, fixed Space/Variance tags, and de Bruijn distances for bound references.
Sums and commuting products sort serialized children; ordered products retain
sequence. Alpha-equivalent binder renamings serialize equally. This is a
structural equality key, not mathematical equivalence and not a persistence
format with a parser: no parser or roundtrip claim is made.

Preflight traverses borrowed input before cloning signatures or allocating output.
AST depth uses max_degree (root depth one) with an additional hard cap of 256;
aggregate node count, child count, and atom rank use max_terms. Aggregate
conservative output bytes use max_input_bytes, including coefficient bit-based
bounds and per-node overhead. Coefficient bits use max_coefficient_bits.
Operations charge traversal, serialization bytes and conservative sorting work.
These are deterministic structural limits, not OS memory/time limits. Caller
construction, Clone/Debug/Drop on arbitrarily deep unvalidated trees are outside
these APIs. Node structure is preserved and output may conservatively reject
small final strings. Full CAS behavior, tensor symmetry equivalence, arbitrary
operator representations and automatic nonzero proofs are outside this slice.

The preflight also reserves conservative work for signature copies and sorting
before allocating them. Obligation extraction checks an aggregate upper bound
(node count times conservative tree bytes) before constructing its output.
Obligation payloads include the same `qft-expr.v1:` version marker.
