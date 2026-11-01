# Chapter 2: independent acceptance and source audit

Ernest, these fixtures check the algebraic steps in problems **2.1–2.7 and 2.8(a–f)** with public reusable operations. They do not call a Chapter 2 solver. A passing Rust test is exact algebra evidence within the declared algebra, or a numerical witness where stated; it is not a theorem about all quantum operators, their domains, or all representations.

## Source and conventions

The reading source is `workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/ocr-compare/reconciled/resolved.md`, lines 596–852, SHA-256 `98c3ae6f4f278701d15a164290c7c6dff071f38ef7922aa5808a6447625791f5`. Its OCR reconciliation is input evidence, not mathematical certification. The Monoclaw `Python/ocr-compare/{README,AGENTS}.md` policy was read: Marker supplies the backbone; Nougat is an untrusted equation cross-check; marker-only equations remain visible and are not silently treated as agreement. No OCR/provider job was run and no source was changed.

This worker independently viewed the parent's renders of PDF pages 43–45, printed folios 19–21: `/tmp/mathphysics-ch2-source-page19.png`, `page20.png`, and `page21.png` (each with the full `mathphysics-ch2-source-` prefix). These confirm every requested exercise statement, the scalar-covariance premise (2.26), the ordering in 2.8(b–d), `hbar/i` in (2.30), the qualification in 2.8(f), and the positive off-diagonal entries in (2.37). This is a bounded source check; this worker has not independently viewed earlier chapter pages containing (2.3)–(2.20). Those upstream equations were read in the reconciled source. Corpus inventories and their acceptance are parent/librarian responsibilities.

The metric is time-first `eta=diag(-1,+1,+1,+1)`. Upper partials are `partial^mu=eta[mu] partial_mu`. We use `[A,B]=AB-BA`, `K_i=M^{i0}`, `J=(M^{23},M^{31},M^{12})`, `H=P^0`, and

`L^{mu nu}=-i hbar (x^mu partial^nu-x^nu partial^mu)`.

No silent `hbar=1` convention is applied. The coefficient derivations retain a commuting symbolic `hbar`; the conjugation fixture uses `hbar=3`; differential relation witnesses use 1, the full orbital bracket fixture uses 2, and the independent constructor/sign maps use 3. Coordinates `tau,east,north,up` are fresh caller labels for four spacetime coordinates; `tau` stands for `x^0=ct`, without an additional dynamics/unit-system assertion.

## What each exercise assessment means

| Part | Source anchors in resolved.md | Executed evidence | Boundary |
|---|---|---|---|
| 2.1 | (2.3), (2.7–8), lines 608–644; exercise 779 | Exact general first-order coefficient calculation with all 16 mixed entries; six antisymmetric lowered parameters | First-order metric preservation is imposed; no finite group classification |
| 2.2 | (2.11–14), lines 668–686; exercise 780 | Conditional coefficient identity for all six independent antisymmetric parameters and transformed tensor antisymmetry | Exact group law, inverse/metric invariance and differentiable expansion are premises; general operator conjugation proof is open in executable coverage |
| 2.3 | (2.14),(2.16), lines 685–697; exercise 781 | Exact symbolic tensor linearization for 16 components and six perturbations, matched to expanded (2.16) | Assumes (2.14); separate free-algebra conjugation test checks first-order order/sign |
| 2.4 | (2.16–17), lines 697–708; exercise 783 | General symbolic component contractions for all nine i,j pairs, plus differential representation witness | Assumes Lorentz bracket law; representation equality is not universal operator proof |
| 2.5 | (2.15),(2.18), lines 692–715; exercise 784 | Exact symbolic vector linearization, four components and six perturbations | Four-vector covariance (2.15) is an independent physical premise, not derived from a label |
| 2.6 | (2.18–19), lines 715–720; exercise 785 | General symbolic contractions of all J/H, J/P, K/H, K/P relations and differential representation witness | Assumes (2.18); K=M^{i0} is essential |
| 2.7 | (2.20),(2.23),(2.25), lines 723–753; exercise 786 | Exact mixed expansion coefficient identity; all 16 differential momentum brackets; unitary and central negative controls | Assumes exact `T(a)T(b)=T(a+b)` and differentiability; unitarity/projective composition alone do not suffice |
| 2.8(a) | (2.26),(2.29–30), lines 758,787–795 | Conditional scalar Taylor coefficients for all six pairs, plus exact independent L30 map/application | Assumes scalar covariance and formal differentiability |
| 2.8(b) | line 797 | Exact composition/application and independently written second-order normal form for L01 L12 | Assumes coordinate coefficients and differentiation pass through commutation with Hilbert-space M |
| 2.8(c) | line 798 | General associative free-algebra Jacobi expansion with distinct arbitrary generators | Associative products must exist on a common algebra/domain |
| 2.8(d) | (2.31), lines 799–802 | Exact free-algebra nested-commutator identity with order-reversal negative control | Identification with L actions is conditional on (a,b) |
| 2.8(e) | line 804 | All 256 exact differential bracket normal forms with independent sparse input maps | Formal polynomial differential algebra; no analytic operator-domain theorem |
| 2.8(f) | line 805 | Algebraic conclusion modulo the commutant; nonzero 7I remains invisible to field/derivative multiplication and orbital/momentum commutator actions | Source's parenthetical absence of Lorentz central charges remains open here |
| Reusable rapidity | (2.37–38), lines 839–852; printed 21 | Exact rational finite boost witness, conditional symbolic metric polynomial, numerical collinear composition/inverse/interval witnesses | Rapidity tests do not establish (2.38) as a general operator exponential identity or solve all of 2.9 |

## Derivations, independent oracles and negative controls

**2.1.** Expanding `Lambda^T eta Lambda` with `Lambda=I+epsilon W` gives coefficient `eta_rr W^r_s+eta_ss W^s_r`. The test obtains this from the full four-index contraction, then compares to the independently lowered expression. Requiring the coefficient to vanish forces antisymmetry, including zero diagonal, over characteristic zero. An unrestricted W does not vanish. A boost has equal mixed time-space entries but opposite lowered entries; the wrong mixed-antisymmetry choice fails the control. Second-order terms are outside this infinitesimal claim.

**2.2.** Starting from the exact group law, let `Lambda'=I+omega'` and set `tilde omega=Lambda^{-1} omega' Lambda`. Metric invariance gives `Lambda^{-1}=eta Lambda^T eta`, hence `tilde omega_rs=omega'_mu nu Lambda^mu_r Lambda^nu_s`. Substituting the infinitesimal U expansion gives (2.13). With antisymmetric coefficients, contraction is `2 sum_{mu<nu} omega_mu nu A^{mu nu}`; equality for every independent parameter determines each antisymmetric coefficient, yielding (2.14). The executable fixture verifies the latter extraction and the transformed tensor's antisymmetry with entirely symbolic Lambda entries. It does **not** implement the whole general invariant-conjugation derivation or establish a representation U, unitarity, or Hermiticity. Those upstream steps are the explicitly open executable gap. Adding a symmetric coefficient is invisible to the probes: the negative control prevents extracting arbitrary non-antisymmetric components. The tensor-antisymmetry calculation itself does not need Lorentzness; that assumption belongs in the upstream conjugation step.

**2.3 and 2.5.** Generic ordered multiplication verifies `(I-epsilon B) A (I+epsilon B)=A+epsilon[A,B]+O(epsilon^2)`; a reversed bracket gives a different answer. With (2.12), a single independent lowered parameter has coefficient `(i/hbar)[A,M^{rs}]`; the factor 1/2 disappears because both antisymmetric index orders contribute. On the covariance side, use `delta omega^mu_r=eta[mu] omega_mu r`: each tensor index contributes separately, and a vector contributes once. Multiplying each coefficient by `hbar/i=-i hbar` yields the explicit four-term Lorentz formula and two-term vector formula. The tests construct the transformation-side coefficients independently from the expanded expected bracket expressions, retaining symbolic hbar. A nonzero boost-bracket sign reversal is rejected. These are exact conditional coefficient derivations, not numerical matrix tests of a general quantum field theory.

**2.4 and 2.6.** Antisymmetry reduces `J_i=epsilon_ijk M^{jk}/2` to the three oriented pairs above. Substituting these and `K_i=M^{i0}` in the expanded source brackets gives, for every i,j,

`[J_i,J_j]=i hbar epsilon_ijk J_k`, `[J_i,K_j]=i hbar epsilon_ijk K_k`, `[K_i,K_j]=-i hbar epsilon_ijk J_k`,

`[J_i,H]=0`, `[J_i,P_j]=i hbar epsilon_ijk P_k`, `[K_i,H]=i hbar P_i`, `[K_i,P_j]=i hbar delta_ij H`.

The symbolic contraction fixture verifies these consequences generally under the input laws, using free linear generator labels and an independently written 3D epsilon rule. A second fixture computes their actual brackets in the polynomial differential representation with `P^mu=-i hbar partial^mu`. It also computes every `[P^mu,P^nu]=0`. Negative controls reject `[K_x,K_y]=+i J_z` and the sign from replacing `M^{i0}` by `M^{0i}`. These differential operators need no positive-definite Hilbert inner product; no physical Hermiticity claim is inferred.

**2.7.** The needed extra property is an **exact representation of the additive translation group**, `T(a)T(b)=T(a+b)`, not merely unitarity. With `A=-i a_mu P^mu/hbar` and `B=-i b_mu P^mu/hbar`, the coefficient bilinear in independent translation parameters is `AB` on the product side and `(AB+BA)/2` on the exponential side. Equality forces `[A,B]=0`, hence every momentum bracket by arbitrary a,b. The common nonzero scalar prefactor is omitted in the fixture and stated here. The free-algebra calculation does not silently impose the group law: the mismatch is actually nonzero before the assumption. Exact Pauli X,Z are unitary but anticommute, demonstrating why unitarity alone is insufficient. They are a discrete phase witness, not a constructed continuous translation representation. Separately `[d_x,x]=I` gives a nonzero central bracket in the differential algebra. These controls do not exclude projective central extensions. No general operator exponential/convergence or strong-continuity theorem is claimed.

**2.8(a).** Taylor expansion gives `delta phi=-omega^r_s x^s partial_r phi`. Its six independent coefficients, multiplied by `hbar/i`, yield the L expression with raised partials. A second test independently writes the exact sparse map `L30=3i(up*d_tau+tau*d_up)` and applies it to `tau*up`, obtaining `3i(up^2+tau^2)`. The opposite sign fails. All 16 constructor components are also compared with independent sparse maps at hbar=3.

**2.8(b–d).** When coordinates are external labels and derivatives commute with Hilbert-space commutation, `[L_a phi,M_b]=L_a[phi,M_b]=L_a L_b phi`; the differential composition acts on the right operand first. The b fixture compares the five nonzero slots in L01 L12 to a hand-expanded product-rule map, including the first-order slot `-tau*d_north`. Dropping coefficient differentiation or reversing composition cannot give that map. The c fixture expands the twelve signed ordered words in cyclic Jacobi; they cancel in the free associative algebra. It also checks associativity and that distinct words do not commute. The d fixture independently verifies `[phi,[A,B]]=[[phi,A],B]-[[phi,B],A]`. Substitution gives exactly the action ordering in (2.31), not its negative. These formal identities apply to associative algebras; undefined unbounded-operator products are not addressed.

**2.8(e,f).** The product rule cancels all second-derivative terms in the bracket, leaving

`[L^{mu nu},L^{rho sigma}]=i hbar(eta^{mu rho}L^{nu sigma}-eta^{nu rho}L^{mu sigma}+eta^{nu sigma}L^{mu rho}-eta^{mu sigma}L^{nu rho})`.

The e fixture compares exact normal forms for all 256 ordered quadruples at hbar=2, including repeated/diagonal indices, using direct sparse maps for the input operators rather than a production Lie-bracket oracle. This is an identity on arbitrary polynomial arguments in this declared differential algebra, stronger than a few sample evaluations. Combining this with (a–d) shows that the difference between `[M,M]` and the source RHS commutes with phi and its derivatives. It does not show that difference is zero. The f fixture constructs nonzero 7I and explicitly shows its invisible commutator action on a sample scalar field, all its first derivative multiplication operators, and all orbital/momentum operators. It also compares R and R+7I through momentum commutator actions. A commutant of the field algebra need not equal the center of every larger operator algebra. We do not claim these shifts produce an arbitrary consistent Lorentz central extension; proving the source's absence statement would require additional representation/cohomology reasoning and is open.

**Rapidity.** For `u=tanh(eta/2)=1/2`, independent rational identities give `cosh eta=5/3`, `sinh eta=4/3`. The exact boost is compared entry-by-entry to a literal matrix, its image of the unit time vector to a literal vector, and its interval to -1. Symbolic c,s retain the nonzero condition `c^2-s^2-1`; there is no unconditional metric-preservation PASS for unconstrained c,s. Numerically beta=0.6 gives gamma=1.25 and positive gamma*beta=0.75 on each axis, matching source (2.37). Collinear rapidities compose by addition and yield `(beta_a+beta_b)/(1+beta_a beta_b)`; inverse/application and a nontrivial interval are checked with 1e-12 absolute/relative metric tolerance and 1e-10 interval tolerance. A deliberately damaged boost fails the metric check. Invalid/lightlike/nonfinite beta and a zero axis are rejected. Finite numerical witnesses are not an exact theorem for all boosts. The source exponential sign (2.38) is visually confirmed, not proved by these fixtures.

## Reproduction and limits

Run from this crate, keeping build products external:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-ch2-fixtures/build cargo test --offline --test chapter_two
CARGO_TARGET_DIR=/tmp/mathphysics-ch2-fixtures/build cargo test --offline --release --test chapter_two
CARGO_TARGET_DIR=/tmp/mathphysics-ch2-fixtures/build cargo clippy --offline --test chapter_two -- -D warnings
```

This worker ran these checks in an isolated staged crate containing copies of the public implementation modules; the manifest records hashes and the actual executed commands/results. The parent owns copying the three deliverables into the repository and its final whole-crate gate. The manifest does not substitute an unexecuted parent integration result. Debug/release counts are repeated profiles, not distinct mathematical results.

Fixtures use checked BigPoly/BigC/ordered/differential operations, with up to 100 million operation units **per helper operation**, resetting the budget between operations; this is not an end-to-end workload/time bound. Current low degrees/orders remain within default terms, degree, coefficient-bit and word-length limits. Exactness means rational/Gaussian-rational polynomial and ordered-word equality in this kernel. There is no formal proof assistant certificate, general noncommuting exponential engine, operator-domain analysis, source-wide completion claim, or proof of absent central charges. Existing book manuscripts and production modules are not owned by this acceptance assignment.
