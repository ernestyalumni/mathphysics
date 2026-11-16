# Symbolic Wilson tooling: bounded implementation plan

Source: request-2026-11-16-wilson-reading.txt. This turn plans future code; no Rust/CUDA implementation or PDD regeneration is claimed. The existing Rust-first programme and reusable-library policy remain authoritative. Original wording and read-only PDD intake are preserved; no provider-backed apply. A generic code-generation prompt is a draft below, not a registered production unit.

## Scope and ownership

Start with symbolic group words, inverses and trace semantics. Reuse Rust/qft-symbolics/src/ordered.rs for the existing free associative polynomial capability and its exact coefficients/resource budget where appropriate. That module does not currently implement group inverses or trace quotients; do not encode inverse semantics by inventing special string names. Keep the generic algebra in qft-symbolics and lattice geometry/adapters in Rust/lattice-gauge. No production solver per textbook subsection and no new crate/dependency without inspection/justification.

## Stages

| Stage | Bounded deliverable | Acceptance |
|---|---|---|
|1: YM-T15| Typed signed group-word and trace kernel in the owning algebra library | Adjacent inverse reduction, word inverse, cyclic trace rules, raw complex vs normalized real trace, resource errors and exact negative controls |
|2: extend T05/T07 boundary| Oriented-lattice adapter consuming that kernel | Plaquette signs +,+,-,-; shifted reverse edge; endpoint gauge cancellation; no geometry clone |
|3: existing T09| Narrow exact Haar integration capability | Specify group/representation before using orthogonality; normalized integral of1=1; independent analytic oracle |
|4: existing T10| Symbolic finite-group partition/observable evaluation | Rational polynomial or other explicitly declared exact representation; normalize raw histograms; compare independent enumeration and solvable case |
|5: existing T11| Group-word directional derivatives | Reuse accepted word/trace kernel; preserve ordering; verify source-located identities independently |
|6: gated T13| CUDA C++ evaluation backend only after a CPU benchmark | Same mathematical expression/input schema, declared numeric error, CPU/GPU parity, timings incl transfer; no GPU speedup or exact-proof claim inferred |

Stage1 is ready for a coding agent to scope through the existing PDD workflow. Later stages are dependencies, not an instruction to implement them all at once. A Mac worker can implement/test Rust; CUDA validation needs a CUDA-capable host. Small symbolic simplifications may have no GPU payoff. Do not turn CPU parallelism or CUDA agreement into proof of an infinite-volume/continuum theorem.

## Acceptance traps to preserve

1. Noncommutativity: Tr(AB)=Tr(BA) is cyclicity, so it is a bad negative control. A=i sigma1,B=i sigma2,C=i sigma3 gives Tr(ABC)=2 and Tr(ACB)=-2, using exact2x2 matrix arithmetic.
2. Path reversal: reverse order **and invert each factor**. ReTr(W^-1)=ReTr(W) for unitary W; mere reversal is false. For W=diag(i,i,-1) in SU3, Tr(W)=-1+2i and Tr(W^-1)=-1-2i; complex traces must not be identified.
3. Real trace is only real-linear. Do not pull arbitrary complex coefficients through Re. First scope can restrict its scalar coefficients to exact real rationals.
4. Existing Torus::plaquette_links returns ordered indices, not orientation signs. Its constructor rejects side lengths1; preserve that explicit software domain although the theory permits degenerate one-site tori.
5. Existing z2_gauge histograms are raw configuration counts. Z_count(beta)=sum_j N_j exp[-2beta(P-j)], so Z_count(0)=2^E. Normalized-Haar Z=Z_count/2^E. The2x2 histogram[32,0,192,0,32] sums256; its normalized Z(0)=1. Twist ratios cancel a common normalization but absolute partition functions/defects need explicit conventions.
6. Different normalized words are not a decision procedure for all finite-dimensional trace identities. Return equality under declared rewrite rules or unknown; never infer non-equality from differing normal forms alone.
7. No silent group/representation/normalization change, commuting of links, open-path gauge-invariance claim, numerical fallback or unchecked resource exhaustion.

The independent symbolic/matrix oracles may intentionally repeat mathematics in a distinct representation. They must not call the same production normalizer to generate expected answers. Existing T11's broad swapped-factor negative-control sentence needs the cyclicity qualification above; do not inherit it literally.
