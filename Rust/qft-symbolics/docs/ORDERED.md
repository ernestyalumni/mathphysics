# Ordered polynomials

`ordered::OrderedPolynomial` is the free associative algebra over commuting
`polynomial::BigPoly` coefficients. Its private `BTreeMap<Vec<String>, BigPoly>`
has one entry per operator word. Empty words are scalars, an empty map is zero,
and zero coefficients disappear. Words concatenate without sorting: `AB != BA`.
There are no implicit metric, Lorentz, adjoint, unitarity, or commutator rules.
Generator and commuting coefficient names inhabit separate namespaces.

All operations take a caller-owned `computation::Budget` and return
`computation::Result<Self>`, except the borrowed `terms()` accessor:

```text
scalar(BigPoly, &mut Budget)
generator(&str, &mut Budget)
terms(&self) -> &BTreeMap<Vec<String>, BigPoly>
try_add/try_sub/try_mul(&self, &Self, &mut Budget)
try_scale(&self, &BigPoly, &mut Budget)
try_neg(&self, &mut Budget)
commutator(&self, &Self, &mut Budget)
truncate(&self, parameters: &[&str], max_degree: usize, &mut Budget)
coefficient(&self, parameter: &str, power: usize, &mut Budget)
substitute(&self, &BTreeMap<String, Self>, &mut Budget)
```

`commutator(X,Y)` means `XY-YX`. Truncation retains monomials whose total degree
in the selected commuting parameters is at most the bound; duplicate selections
count only once. Other parameters and operator words are unchanged.
`coefficient("epsilon",1,...)` selects exactly one power of epsilon and removes
it. This is coefficient extraction, without factorials or derivatives.
Substitution replaces operator generators simultaneously, distributing each
replacement in original word order. Unmapped generators remain. Replacement
expressions are not recursively substituted and scalar parameters are untouched.

For `U=1+epsilon J` and a separately justified inverse through first order,
`U^-1=1-epsilon J+O(epsilon^2)`, the library expands and truncates
`U^-1 A U = A+epsilon [A,J]+O(epsilon^2)`. This is a conditional algebraic
identity. Identifying the generators with a physical representation requires
separate assumptions and evidence. These checks do not exclude central terms.

## Bounds and validation

Names must be nonempty and have no Unicode whitespace or control characters.
Names are validated in every supplied polynomial, selected parameter list,
and replacement map, including unused replacements and zero-input operations.
All coefficient arithmetic uses checked BigPoly APIs; no legacy panic operators
are called. Inputs remain unchanged on failure, but consumed budget is retained.

Each input, intermediate and result is bounded independently; limits do not
bound the sum of simultaneously live heap allocations or allocator overhead.
`Terms` counts aggregate coefficient monomials across all operator words.
`InputBytes` counts every operator-name occurrence and every commuting-name
occurrence across coefficients, including repeated names. Replacement-map keys
and all replacement payloads also receive one combined byte check.
`CoefficientBits` bounds the sum of the maximum rational-component bit lengths
of all stored Gaussian coefficients. `WordLength` limits each operator word;
aggregate word storage is bounded by `max_terms * max_word_length`, charged
validation work, and aggregate symbol bytes. Commuting degree is bounded per
monomial. Both word length and commuting degree have an additional hard cap of
256 even when caller limits are larger. Algorithms are iterative with no
recursive substitution or expression trees.

Admission is deliberately conservative: a pending insertion counts existing
and incoming terms/bytes/bits before merging, and multiplication admits the
Cartesian coefficient term expansion before allocating its concatenated word.
A conservative four-times component-bit estimate admits coefficient products;
checked BigPoly arithmetic also enforces its own intermediate limits. A small
budget can therefore reject an expression whose fully cancelled result would
fit. This behavior is intentional; increase an explicit limit for that work.
Resource counters bound declared algebraic work, not wall time or physical RAM.
Validation and candidate key work are charged before growth. Failed operations
return typed errors rather than partially modifying their inputs.

## Verification and provenance

`tests/ordered_algebra.rs` checks full nested-commutator word expansion and
Jacobi cancellation, noncommutativity, commutator signs, first-order conjugation,
selected total degree, coefficient extraction, simultaneous substitution,
invalid names, zero-input validation, aggregate resource limits, hard caps,
and exact coefficients larger than i128.

Implementation is manual and prompt-first, as authorized after the earlier
external-provider approval rejection. The prompt is
`prompts/ordered_rust.prompt`. No provider generation or automatic PDD sync
success is claimed. Parent integration owns module export and graph mapping.
