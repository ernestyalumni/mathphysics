# Srednicki research checks: what passing means

Ernest, the book supplies claims to investigate. The reusable library supplies
arithmetic, sparse polynomials, Clifford words and explicit convention adapters.
Independent component calculations supply evidence. Neither your notes nor an
OCR transcription is an infallible answer key. These ten tests assess individual
claims and boundaries; they do not complete chapters 35, 47 or 50.

## Convention API and signs

`conventions::EpsilonConvention::{Library, MonoNotes}` has no default. Both use
ordered arrays `(1,2)` and contract the second epsilon index. For lower components
`v`, `raise(v) = upper12 * (v[1], -v[0])`; for upper components `u`,
`lower(u) = lower12 * (u[1], -u[0])`. Thus lower-after-raise is identity, whereas
two raises give minus identity. Library has `upper12=+1, lower12=-1`; MonoNotes
has `upper12=-1, lower12=+1`.

`bracket(a,c) = upper12 * (a[0]*c[1] - a[1]*c[0])` uses commuting components.
It is not the anticommuting field bilinear in (35.25). Reversing the profile
reverses one bracket; a product of two brackets is unchanged. The existing
`spinor`, `twistor`, and matrix APIs keep their existing signs. In particular,
existing twistor angle uses the star-component determinant, while square uses
minus the phi-component determinant. The adapter does not silently rename those.

Conversion preserves the underlying lower component array:
`source.convert_lower(target,v,b)` returns those same components;
`source.convert_upper(target,u,b)` multiplies by
`source.upper12()*target.upper12()`. Lowering converted upper components with the
target recovers the same lower components as lowering the input with the source.
Holding both upper and lower arrays fixed while reversing epsilon is impossible.

All five polynomial methods (`raise`, `lower`, `bracket`, `convert_lower`,
`convert_upper`) are generic over `Coefficient`, borrow their input arrays, take
`&mut Budget`, and return `computation::Result`. Even identity conversion validates
through checked arithmetic. `upper12` and `lower12` are const methods returning
`i64`. There is no dependency on expression, checked-expression, or replay APIs.

## Claim → API → independent oracle → negative control → meaning

All named tests below live in `tests/textbook_research.rs` and import the public
`qft_symbolics::conventions` export. No production chapter solver is introduced.

| Test | Claim and library API | Independent oracle | Negative control | What pass means |
|---|---|---|---|---|
| `ch35_profiles_preserve_components_and_expose_signs` | (35.1), second-index contraction prose: all convention methods on `Polynomial<BigC>` | Explicit `(2,3) -> (3,-2)` or `(-3,2)`, determinant of `(2,3),(5,7)` is -1 | Mixing profiles without conversion changes recovered components; zero-work identity conversion fails | Exact symbolic round trips for both profiles and conversion directions; chosen signs remain explicit |
| `ch35_sigma_contraction_all_components_independent_epsilon` | (35.4), (35.19)-(35.20): `BigPoly` arithmetic | Locally specified Pauli entries and explicit lower epsilon, all 16 contraction components and 16 raised-sigma entries per profile | Nonzero epsilon contraction distinguishes +2 from -2 | Complete finite component checks of these displayed tensors, including transpose placement; not all ch35 representation theory |
| `ch47_symbolic_four_trace_pairings_and_rename_control` | (47.13): `try_trace_exact`, `try_metric_with<BigC>` | Independently written three pairings with signs +,-,+ | Flip middle sign; distinguish punctuation-bearing slot names | Exact polynomial identity with formal metric variables for two label sets; no gamma5 assertion |
| `ch47_component_trace_words_through_six_pauli_oracle` | (47.1), (47.4), (47.6), (47.8), (47.13) and recurrence extension: `try_trace_exact`, `try_eval` | Locally built 4x4 Pauli blocks and ordinary matrix multiplication; explicit h=-eta substitution | Timelike two-gamma trace is +4, not -4 | 5,461 evaluated component words, all lengths 0 through 6; finite witnesses beyond the symbolic four-word check |
| `ch47_formal_dimension_and_component_sandwiches` | (47.18)-(47.21): `try_sandwich_srednicki_with<BigC>` | Formal one-word coefficient d-2; independently contracted Pauli matrices for every component word of lengths 0 through 3 (85 words) at d=4 | Formal empty-word d=5 gives -5 rather than matrix -4; trace word of length 14 rejected | Symbolic one-word relation and finite d=4 sandwich witnesses, not a d-dimensional matrix representation or formal proof of all of (47.20)-(47.21) |
| `ch50_schouten_symbolic_both_profiles_and_renamed_components` | Problem 50.3, (50.36): `bracket`, `BigPoly` | Three determinant products in the book's stated order, with independent component variables | Flip middle sign; rename every component using punctuation and Unicode; fixed determinant witnesses | Exact commuting polynomial identity, both profiles; wrong-sign polynomial remains nonzero |
| `ch50_exact_rational_boundary_and_invalid_domains` | Arithmetic/domain support for using twistor products: `BigQ::parse`, `bracket`, `try_eval`; legacy `Ratio` boundary | `(N/3,0)` bracket `(0,3)` equals N, N=340282366920938463463374607431768211457 | Narrowing fails, zero rational denominator fails, missing variable fails; x/x has zero cross residual against 1 but cannot evaluate at x=0 | Arbitrary-precision exact witness and explicit invalid-domain behavior; the Ratio control uses the existing bounded coefficient backend |
| `ch50_massless_factorization_symbolic_and_rational_pauli_witness` | (50.12), (50.15): `BigPoly`, `try_eval` | Invert the displayed 2x2 block; separately reconstruct it using local Pauli entries and a complex rational spinor | Nonreal off-diagonal distinguishes a transpose/wrong sigma2 sign | Exact symbolic null residual with independent phi/star variables, plus one Lorentzian conjugate rational factorization witness; does not verify the angular formula (50.9) or all of problem 50.2 |
| `ch50_four_bracket_chiral_trace_rational_component_witness` | Problem 50.4, (50.37): four convention brackets and `BigPoly::try_eval` | Four independently inverted momenta, local Pauli-block slash product, upper chiral block trace | Opposite chiral block gives a different answer; trace is nonzero | One complex-component witness of the four-bracket/chiral-trace relation; not an all-variable symbolic proof |
| `ch50_nonnull_determinant_sign_control` | Boundary of (50.15): constant `BigPoly` | Elementary determinant of the displayed matrix at p=(2,1,0,0) | +3 differs from eta(p,p)=-3 | Guards the sign that a massless-only check cannot observe; does not certify the adjacent OCR prose |

The matrix oracle never calls production `matrices`, `spinor`, `twistor`, or
`lorentz` helpers. It uses the existing small exact Gaussian rationals only for
these modest component witnesses. Arbitrary-precision polynomials and Clifford
functions are the systems under test. Symbolic pairing expectations share the
public metric constructor, so the independently encoded component substitution
and Pauli calculations provide a separate check of metric interpretation.
The null-momentum inversion is test-only mathematics, not a production solver.

## Source locators and OCR limits

Source root, relative to the workspace:
`workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/ocr-compare/reconciled/`.
Read the canonical Monoclaw `Python/ocr-compare/README.md` and `AGENTS.md` before
using these existing outputs. No OCR, provider, or vision-judge run was started.
The existing `resolved.md` and `equations_resolved.json` were consulted directly.

Verified `resolved.md` locators (one-based lines in the recorded source):

- (35.1): 6644-6648, including the second-index contraction sentence.
  The JSON labels (35.1) Nougat-only; the complete sign statement was read in the
  resolved text, not inferred from that status. (35.4): 6664; (35.19): 6741;
  (35.20): 6746, each marked agreement.
- (47.1)-(47.4): 9516-9528; odd traces (47.6): 9542-9543;
  four trace (47.13): 9587-9588, custom high-confidence correction in JSON.
  The existing image `pages_small/page_319_eq47.13.png` was also directly inspected:
  four distinct slashed vectors and pairings `(ad)(bc)-(ac)(bd)+(ab)(cd)`.
- (47.18)-(47.21): 9618-9648. (47.20)-(47.21) still have garbled slash glyphs
  in resolved text; the tests do not claim to recover all vector ordering from
  that text. Formal d-2 and the independently contracted d=4 components are the
  bounded evidence. No gamma5 dimensional prescription is asserted.
- (50.12): 10024 (agreement); (50.15): 10042-10046 (both-correct resolution,
  high confidence). Existing `pages_small/page_334_eq50.15.png` was directly
  inspected for the displayed matrix signs. This crop cuts off its lower row;
  the full matrix is provided by resolved text/JSON. The following OCR prose
  assigns determinant `-p0^2+|p|^2`; the displayed matrix instead gives
  `p0^2-|p|^2`. We record that discrepancy without changing the source or notes
  or claiming a fresh full-page inspection of the printed prose.
- Problem 50.3 / (50.36): 10185-10190; problem 50.4 / (50.37): 10192-10197.
  Both are Marker-only in the reconciliation JSON. Their statements were located
  in the actual resolved book; no new source-image verification is claimed.

The original reconciliation summary contains 1,366 agreeing, 766 conflicting,
651 Marker-only and 228 Nougat-only entries, with many Nougat dropped ranges and
repetitions, including dropped (47.2)-(47.3) and (50.29)-(50.30). These are pipeline
counts, not verified book coverage. No fresh page-coverage audit was performed;
we do not treat the local corpus as complete or the original conflict count as
remaining unresolved conflicts. Strip filenames carry pipeline page identifiers;
we do not assert an unverified printed-folio conversion.

Source SHA-256:

```text
resolved.md
98c3ae6f4f278701d15a164290c7c6dff071f38ef7922aa5808a6447625791f5
equations.json
1a6dfcc3fd2f61fd854cd9dcab16191d5f509813b99ea1c994221db54f7f90c9
equations_resolved.json
8438cff68f74fc52535a11222de46ef926f44053ef68916e0280668f50ae35da
```

## Reproduction and limits

From the mathphysics repository root:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-next-corpus cargo test --offline --locked --manifest-path Rust/qft-symbolics/Cargo.toml --test textbook_research
CARGO_TARGET_DIR=/tmp/mathphysics-next-corpus cargo test --offline --locked --release --manifest-path Rust/qft-symbolics/Cargo.toml --test textbook_research
```

Default checked limits are 1,000,000 work units, 8,192 terms, degree 64,
65,536 coefficient bits, 1 MiB input bytes, and word length 12. Each of the 5,461
trace witnesses gets a fresh default-sized work budget; the large-rational check
explicitly allows 10,000,000 work units. Exhaustion is an error, never an identity
result. The generic polynomial operators used in test expectations retain legacy
panic-on-budget-error semantics; this does not make the new production adapter
infallible. Clifford hard caps remain 12 for traces/normal ordering and 8 inner
slots for sandwiches. The finite matrix oracle is four-dimensional.

`validation/research-conventions.json` inventories these tests individually for
parent review; it is not the corpus ledger, a chapter-completion metric, or a
replacement for running the tests. `corpus.json`, notes, Cargo, library exports,
and PDD maps are parent-owned. Manual source/prompt implementation only; no
provider-backed apply, regeneration, or synchronization fidelity is claimed.
Read-only intent plan: `approved-mathphysics-reusable-qft-work-own-only--a2a32e2a`.

Verification: all ten new public-library tests passed in debug and release.
The existing convention, textbook and extended-textbook suites also passed in
release. The two owned Rust files pass rustfmt check. Build products and run logs
are outside the repository under `/tmp`; no commits or provider calls were made.
