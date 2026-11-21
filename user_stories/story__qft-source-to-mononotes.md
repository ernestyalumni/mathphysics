# Relearn QFT through traceable mono-notes

## Story

As Ernest reacquainting myself with quantum field theory, I want a Srednicki-based reading path and carefully checked derivations connecting it to Wilson lattice Yang–Mills and the mass-gap question, so that the mono-notes become a coherent account I can study and extend.

Source: `docs/plans/yang-mills-mass-gap/request-2026-11-15-continuation.txt`, outcome e. Manually drafted acceptance source; no generated contract or semantic model check.

Refined by Ernest's reading corrections in
`docs/plans/yang-mills-mass-gap/request-2026-11-21-explicit-exposition.txt`.
He must be able to follow the derivation within this handbook, with local
definitions and links, without retrieving a cited paper or guessing an
unstated substitution. This is a manually maintained independent acceptance
source; the LaTeX manuscript remains outside the generated prompt graph.

## Acceptance

1. Inventory existing source and note coverage before requesting new OCR. Cite the actual source's chapter/section/equation/page and extraction version; do not guess chapter numbers.
2. Prefer original TeX, then existing audited extracts. If OCR is required, follow Monoclaw's ocr-compare instructions and retain coverage/conflict information. Keep original books and papers unchanged; bulk book extracts remain outside version control.
3. Explain the free-field/path-integral prerequisites, gauge symmetry and curvature, nonabelian interactions, renormalization/asymptotic freedom, Euclidean reconstruction and spectral-gap question in an explicit dependency order. Srednicki does not by itself replace the constructive/lattice literature.
4. Derivations display the algebra and conventions, separate known results from hypotheses, and identify any unresolved source disagreement. A second reviewer checks mathematical content before incorporation.
5. Integrate reviewed content in `documents/mono-notes` using its STYLE.md, keeping code maps/agent instructions in companion docs. No duplicated production solver per subsection; use the canonical reusable-library story for tool work.
6. Every externally cited formula actually used in a derivation MUST be displayed in the handbook and have a stable local equation label. Uses MUST link to that local display; the source's original equation/page citation remains for provenance. A sentence containing only "Jaffe–Witten Eq.(2)" does not meet acceptance. A distant local reference alone is insufficient when a short repeated formula is needed to follow the current step.
7. Calculation verbs such as "substitute," "specialize," "contract" and "rearrange" MUST accompany the actual displayed work. Identify the starting expression, define introduced quantities, display replacement assignments and intermediate algebra, and retain the explicit terms that vanish. Prose may explain the steps; it MUST NOT replace them.
8. For example, from the connected two-point function, define the local observable O and its arguments, display A(x)=O(x) and B(y)=O(y), replace both operators and both one-point expectations, show why translation invariance makes the mean independent of position, and evaluate the zero-mean subtraction as (0)(0)=0. Show the resulting correlation and link directly to the reproduced clustering inequality in the handbook.
9. When declaring a specialization, the added display MUST preserve the required assumptions, sign/normalization conventions and source meaning. Do not silently claim the whole cited theorem was proved merely because its equation was reproduced.

## Planned evidence

YM-R09 source inventory and short reading map, followed by one checked derivation and a successful notes build when manuscript edits are actually made. Reading-map creation alone is not completion of the manuscript.

For the exposition refinement, inspect the displayed clustering definition,
substitution and linked inequality, then follow those links in both built
editions. Spot-check another words-only substitution (for example the Wilson
loop expectation or Hamilton equations recovering Klein–Gordon), and check
that every introduced symbol has a nearby meaning. This evidence must come
from the rendered mathematics and link targets, not only a text search for
`\eqref`. The current repair targets the shared QFT and Yang–Mills topics;
it does not assert an exhaustive audit of all imported legacy chapters.
