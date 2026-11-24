# YM-A04: independent transcription of the public theorem

Label: unreviewed (claude-opus55-strat). Source: `yang-mills-group/OpenAI-YM-20261028` at 3456e5e: `PublicTheory.lean`, `PublicWilson.lean`, `ChallengeSpec.lean`, `PublicView.lean`, `Challenge.lean`.

**Protocol.** This transcription was written from the Lean files alone. `ARM/statement.md` had not been opened by this agent when this section was committed. The comparison with `statement.md` is a separate, later commit.

Conventions:
- Mathlib's inner product ⟨·,·⟩ is conjugate-linear in the first slot.
- "Eventually as r → 0⁺" means for all sufficiently small r > 0.
- Path:line references are to the files above.

## 1. Gauge group (`PublicTheory.lean:25–45`)

A gauge group is given by data and one condition:
- an integer n ≥ 0;
- a **closed** subgroup G of U(n+1);
- a real subspace 𝔤 ⊂ M_{n+1}(ℂ) with X ∈ 𝔤 ⟺ exp(tX) ∈ G for all t ∈ ℝ.

So 𝔤 is the Lie algebra of the closed subgroup G, and G is compact.

**Simple** means two things:
- 𝔤 is non-abelian: some X, Y ∈ 𝔤 have [X, Y] ≠ 0;
- every real subspace I ⊂ 𝔤 with [𝔤, I] ⊂ I is 0 or 𝔤.

G is **not** required to be connected.

## 2. Spacetime, tests, Wightman data (`PublicTheory.lean:51–321`)

**Spacetime and tests.**
- Spacetime is ℝ⁴ with the Euclidean norm, and time is the last coordinate.
- The closed forward cone is V̄₊ = {x : |x_spatial| ≤ x_time}. The Minkowski square is x² = x_time² − |x_spatial|², and spacelike means x² < 0.
- Test functions are complex Schwartz functions on ℝ⁴, and f̄ is the pointwise conjugate.

**Field system.**
- A complex Hilbert space H and a complex vector space Λ of labels, of any dimension.
- A complex subspace D ⊂ H, a vector Ω ∈ D, and a map φ: Λ × 𝒮 → End(D), linear in each argument.
- Write φ(λ₁, f₁)⋯φ(λ_k, f_k)Ω for nested application, the first factor outermost.
- The vacuum expectation is ⟨Ω, φ(λ₁, f₁)⋯φ(λ_k, f_k)Ω⟩.

**Domain axioms.**
- ‖Ω‖ = 1 and D is dense.
- D **equals** the linear span of all vectors φ(λ₁, f₁)⋯φ(λ_k, f_k)Ω with k ≥ 0, so Ω is cyclic and D is exactly the polynomial domain.
- Temperedness: for each λ and ψ ∈ D, the map f ↦ φ(λ, f)ψ is continuous from 𝒮 to H in norm.

**Adjoint axiom.** There is an antilinear J: Λ → Λ with ⟨φ(λ, f)ψ, χ⟩ = ⟨ψ, φ(Jλ, f̄)χ⟩ for all ψ, χ ∈ D. J is not required to be an involution.

**Locality.** If every point of supp f minus every point of supp g is spacelike, then φ(λ, f)φ(λ′, g)ψ = φ(λ′, g)φ(λ, f)ψ on D. This is commutation only; there is no fermionic grading.

**Faithfulness** (used later). For every λ ≠ 0 there is some f with φ(λ, f) ≠ 0.

**Translations.** Bounded operators T(a), a ∈ ℝ⁴, with:
- T(0) = 1 and T(a)T(b) = T(a + b);
- norm-preserving, hence unitary;
- strongly continuous.

**Energy–momentum measure.** A map P from subsets of ℝ⁴ to bounded operators:
- P(∅) = 0 and P(ℝ⁴) = 1;
- multiplicative and self-adjoint on measurable sets;
- strongly countably additive;
- with finite measures μ_ψ satisfying μ_ψ(A) = ⟨ψ, P(A)ψ⟩ for measurable A.

**Proper orthochronous Lorentz group.** Invertible linear maps of ℝ⁴ with determinant 1 that preserve the Minkowski square and satisfy g x ∈ V̄₊ ⟺ x ∈ V̄₊. Poincaré coordinates are pairs (Λ, a).

**Poincaré action.** A type S with a unit, a multiplication and a topology. There is no separate group axiom; the coordinates below make S a copy of the Poincaré group. The data:
- unitaries U(g): U(1) = 1, U(gh) = U(g)U(h), strongly continuous;
- linear automorphisms g ↦ L(g) of Λ, with no continuity and no finite-dimensionality;
- a homeomorphism c: S → SO⁺(1,3) × ℝ⁴ with c(1) = (1, 0) and c(gh) = (Λ_g Λ_h, Λ_g a_h + a_g);
- an embedding of translations a ↦ τ(a) with c(τ(a)) = (1, a), additive, τ(0) = 1;
- test functions transform by (g·f)(x) = f(Λ_g⁻¹(x − a_g)).

A **Wightman theory** consists of the field system, J, T, P and the Poincaré action, with U(τ(a)) = T(a).

**W0.**
- (a) **Spectral representation.** For all ψ and a ∈ ℝ⁴: ∫ e^{i a·p} dμ_ψ(p) = ⟨ψ, T(a)ψ⟩, with the **Euclidean** pairing a·p = Σ_k a_k p_k.
- (b) **Positive energy.** P(ℝ⁴ ∖ V̄₊) = 0.
- (c) **Unique vacuum.** Every vector fixed by all T(a) is a multiple of Ω.

**W1.** The domain axioms and the adjoint axiom.

**W2.**
- L is a representation.
- D is invariant under every U(g), and U(g)Ω = Ω.
- U(g)φ(λ, f)ψ = φ(L(g)λ, g·f)U(g)ψ for ψ ∈ D.

**W3.** Locality.

## 3. The Wilson lattice and its limit (`PublicWilson.lean`)

**Lattice and Wilson law.**
- For an integer ℓ ≥ 1 the lattice is (ℤ/ℓ)⁴, and a configuration assigns an element of G to each pair (site, axis).
- The plaquette is U_p = U(x,i) U(x+e_i, j) U(x+e_j, i)⁻¹ U(x, j)⁻¹ for i < j. The energy is E(U) = Σ_{x, i<j} (1 − Re tr U_p/(n+1)).
- The Wilson law at β, with a measure h on G, is the product measure h^{⊗links} with density e^{−βE}/∫e^{−βE}.

**A Wilson family** is chosen data:
- strictly increasing integer sequences L(p) ≥ 1 (the periods) and m(N) ≥ 1 (the denominators);
- a **two-index** coupling β(N, p);
- a probability Haar measure h on G;
- a type of tokens 𝒯;
- for each N, p, a token v and a site x, an **arbitrary** complex function insertion_{N,p}(v, x, U) on configurations of the lattice of side ℓ = L(p) m(N).

From this data:
- The spacing is a_N = 1/m(N), and site y sits at the point a_N · valMinAbs(y) ∈ (−L/2, L/2]⁴. So the physical torus has the integer period L(p).
- The finite moment is the Wilson-law expectation of ∏ᵢ insertion_{N,p}(vᵢ, xᵢ, U) at β(N, p).
- **Tests.** Write Sep_n ⊂ (ℝ⁴)ⁿ for the open set of n-tuples of pairwise distinct points. A test is a C^∞ function with compact support in Sep_n. It is **suitable** for p if its support lies in the open box (−L(p)/2, L(p)/2)^{4n}.
- **Smearing.** smear(N, p, v, t) = a_N^{4n} Σ_{x ∈ sitesⁿ} t(places of x) · moment(x), a Riemann sum.
- **Convergence** along an order k ↦ p_k, to limit functions S_v on (ℝ⁴)ⁿ: for all n, tokens v and tests t,
  - (i) t is eventually suitable for p_k;
  - (ii) there is a row R(p) with smear(N, p, v, t) → R(p) as N → ∞ for every suitable p;
  - (iii) R(p_k) → ∫ t S_v as k → ∞.

  First the mesh limit at fixed physical period, then the volume limit along the order.

## 4. Regularity of the family (`ChallengeSpec.lean:12–22`)

- **Running.** There is a slope b > 0 with β(N, p)/log(1/a_N) → b as N → ∞, uniformly in p. The slope is **existential**: any positive number.
- **Insertions.** Each is continuous in U and invariant under lattice gauge transformations U_e ↦ t(x) U_e t(x + e_μ)⁻¹.
- **Locality.** For each token a finite stencil of integer offsets, independent of N and p, such that insertion(v, x, U) depends only on links with both ends in x + stencil. In physical units the support therefore shrinks like a_N.

## 5. The hierarchy (`ChallengeSpec.lean:24–80, 117–143`)

A hierarchy over a family C has these parts.

**The Wightman part.**
- A Hilbert space, labels and a Wightman theory.
- For each d ∈ ℕ a finite-dimensional complex normed space F_d, the fields of scaling degree d, with linear maps ι_d: F_d → Λ.
- A map π: 𝒯 → ⊔_d F_d from tokens to tagged fields.

**Phases.**
- A compact metric (Borel) space X of "characters", a measure ν on X, and a selected point s₀ ∈ X.
- Functions corr(s; a; x) for s ∈ X ∪ {none}, tagged fields a₁…aₙ and points x ∈ (ℝ⁴)ⁿ.
- A continuation cont(a; z; 𝐱) for complex times z ∈ ℂⁿ and spatial points 𝐱 ∈ (ℝ³)ⁿ.

**The OPE part.**
- For each j a finite set of channels; the actor of a channel is a field in F_j.
- Coefficients c_{a,j,i}(x, o).
- A field Θ ∈ F₄, the centred trace.

**Reconstructed** requires:
- **Axioms.** W0–W3 hold, and the fields are faithful.
- **Labels.** Each ι_d is injective. Some φ(ι(a), f) is nonzero. Every label is some ι_d(a); π is surjective.
- **Phase measure.** ν is a probability measure with full support (every open set has positive measure).
- **Correlations.**
  - corr is symmetric under permuting slots, and C^∞ on Sep_n for each s.
  - s ↦ corr(s; a; x) is continuous.
  - ∫ corr(s; a; x) dν(s) = corr(none; a; x): the **mixture**.
- **Analyticity.**
  - For fixed spatial points, cont is holomorphic in z on the tube {Im z₁ < ⋯ < Im zₙ}.
  - **Root.** For x with strictly increasing real times, cont(a; i t; 𝐱) = corr(s₀; a; x).
  - **Boundary values.** For every strictly increasing h and Schwartz tests f_i, as r → 0⁺,

    ∫ cont(a; t + i r h; 𝐱) ∏ f_i(x_i) dx → ⟨Ω, φ(ι a₁, f₁)⋯φ(ι aₙ, fₙ)Ω⟩.

**Operator product expansion.** Write w(a) = Σ degrees.
- **Within(r, M, η).** All x_i lie within Mr of o, pairwise at least ηr apart.
- **Clear(c).** The y_j lie at least c from o and from each other.
- **Coefficients.** For all M ≥ 0, η, ε > 0 and every channel i of degree j: eventually as r → 0⁺, |c_{a,j,i}(x, o)| ≤ r^{j − w(a) − ε} for all x, o in Within.
- **Remainder.** For all M ≥ 0, η, c, ε > 0 and D ∈ ℕ: eventually as r → 0⁺, for all x, o, y with Within and Clear and every s ∈ X ∪ {none},

  |corr(s; a@x, b@y) − Σ_{j ≤ D} Σ_i c_{a,j,i}(x, o) corr(s; actor_{j,i}@o, b@y)| ≤ r^{D + 1 − w(a) − ε}.

  The slots of (a, b) are enumerated by a fixed bijection, which the permutation symmetry makes harmless.

**Trace test.**
- **Centred.** ⟨Ω, φ(ιΘ, f)Ω⟩ = 0 for all f.
- **Amplitude.** 0 < 3 dim_ℝ 𝔤/π⁴.
- **Short distance.** r⁸ (log 1/r)² · corr(s₀; Θ, Θ; (0, −r/2), (0, r/2)) → 3 dim_ℝ 𝔤/π⁴ as r → 0⁺.
- **Excludes.** That two-point function does **not** agree, on any interval (0, r₀), with a finite Wick profile. A Wick profile is a finite sum Σ c ∏ᵢ ∫ (2E)⁻¹ E^{k_i} N_i(p) e^{−rE} dp over modes, where:
  - E = √(m_i² + p²) with m_i ≥ 0;
  - N_i is continuous and positively homogeneous of integer degree.

**Properties along an order k ↦ p_k.**
- The family converges along p_k to S_v = corr(none; π∘v).
- The hierarchy is reconstructed.
- The OPE holds.
- The trace test holds.

## 6. The theorem (`ChallengeSpec.lean:186–188`, `Challenge.lean`)

**JointExistence.** For every compact simple gauge group G there is a Wilson family C satisfying the regularity of §4 such that:
- for every strictly increasing p: ℕ → ℕ there are a strictly increasing q and a hierarchy T;
- T has the properties of §5 along k ↦ p(q(k)).

`Challenge.lean` states it with `sorry`; the solution supplies the proof.

**Quantifier order.** ∀G ∃C ∀p ∃q ∃T.
- The family is fixed for the group.
- The subsequence and the whole hierarchy, including the Hilbert space, the phases and the selected phase, may depend on p.

`PublicView.lean` is not part of the statement. It builds the construction's Wightman theory from `YM.Boundary.Wightman`:
- Hilbert space, domain, Ω, fields and J are taken verbatim;
- labels come from the construction's field table;
- the Poincaré type is `YM.Mink.Poincare`, with identity coordinates;
- labels transform by the construction's spin representation of Λ.

It proves W0–W3 for that view from a `Realizes` hypothesis (`view_axioms`), and identifies its vacuum expectations with the construction's `W` (`publicMoment`).

## 7. What I notice from the Lean alone (before reading `statement.md`)

These are observations, not findings. Each may already be in `statement.md` or the conceded non-claims.

1. **Amplitude is vacuous.** Given simplicity, dim 𝔤 ≥ 3, so 0 < 3 dim 𝔤/π⁴ always holds. The real content is in `shortDistance`.
2. **The running slope is free.** Any b > 0 is allowed. The statement does not require the asymptotic-freedom coefficient of G (in Wilson normalization, β ≈ (2N)·2b₀ log(1/a), with b₀ fixed by G). The trace-test amplitude 3 dim 𝔤/π⁴ does not depend on the coupling, so it does not fix b either.
3. **β depends on the volume index p**, controlled only by the uniform slope limit. The defect of the conditional-gap chain needs two volumes at equal bare coupling (YM-O05 §3).
4. **The lattice limit is the mixture, not the selected phase.**
   - Lattice convergence is stated for corr(none), the ν-average over phases.
   - The Wightman theory is tied to corr(s₀) through root and boundary values.
   - That the selected phase's Schwinger functions are lattice limits is therefore not stated. Only their ν-average is, and continuity plus full support relate s₀ to the average only weakly.
5. **Insertions are arbitrary within the regularity class**, and may depend on N and on p.
   - This holds renormalization, including p-dependent renormalization.
   - Π must be surjective onto every tagged field, so every field of the hierarchy must come from some lattice token.
6. **Temperedness and the domain are stronger than the classical axioms.**
   - Norm continuity of f ↦ φ(f)ψ is stronger than tempered matrix elements.
   - D equals the polynomial domain, rather than containing it.
   - J is not required to satisfy J² = 1.
   - Label transformations need no continuity and no finite-dimensionality. That is weaker than the usual finite-dimensional representation of SL(2, ℂ).
7. **The spectral pairing is Euclidean** (W0(a)). This is equivalent to the Minkowski form, because V̄₊ is invariant under p_spatial ↦ −p_spatial. Positive energy corresponds to the convention T(a) = e^{iHa₀}⋯. This matches the boundary-value direction Im z₁ < ⋯ < Im zₙ and the pullback (g·f)(x) = f(Λ⁻¹(x − a)).
8. **cont is determined by corr(s₀).** The root condition fixes it on {i t : t strictly increasing}, an open piece of a totally real n-dimensional subspace inside the connected tube. By the identity theorem, cont is determined there, so it carries no independent data.
9. **No condition at coincident points.** Tests have compact support in Sep_n, so nothing is said there; this is standard.
10. **Conjunct (i) of convergence is automatic**, since L(p_k) → ∞.
11. **Scaling degrees are natural numbers.** The OPE uses integer degrees; anomalous dimensions must be absorbed by the ε in the exponents.
12. **The OPE constants are 1.** Both OPE bounds are stated with constant 1, and "eventually in r" makes this equivalent to a constant ≤ r^{−ε′}.

The candidates for YM-A03, the clause red-team, are 2, 3, 4 and 5. Of these, 4 is listed in `01-WHAT-OPENAI-PROVED.md` §2 as "no phase is shown to be a limit of lattice states". I have not checked whether 2 is listed.
