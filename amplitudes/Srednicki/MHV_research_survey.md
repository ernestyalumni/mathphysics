# MHV Amplitudes: Comprehensive Research Survey
### For Ernest Yeung — February 2026
### Context: Post-Srednicki Ch. 34, building toward publishable MHV work with Cadabra2

---

## Table of Contents

1. [MHV Amplitudes: State of the Art (2024–2026)](#1-mhv-amplitudes-state-of-the-art-2024-2026)
2. [The Parke-Taylor Formula and Extensions](#2-the-parke-taylor-formula-and-extensions)
3. [Spinor-Helicity Formalism Tools: Cadabra2 + Python](#3-spinor-helicity-formalism-tools-cadabra2--python)
4. [Low-Hanging Fruit: Specific Research Opportunities](#4-low-hanging-fruit-specific-research-opportunities)
5. [Recommended Reading Order](#5-recommended-reading-order)
6. [Concrete Next Steps](#6-concrete-next-steps)

---

## 1. MHV Amplitudes: State of the Art (2024–2026)

### 1.1 What Are MHV Amplitudes?

A **Maximally Helicity Violating (MHV)** amplitude is an *n*-gluon scattering amplitude
where exactly 2 external legs carry negative helicity (all others positive). The name
comes from the fact that helicity is maximally violated relative to the all-same-sign case
(which vanishes).

The key result is the **Parke-Taylor formula** (1986):

$$\mathcal{A}_n(1^+,\ldots,j^-,\ldots,k^-,\ldots,n^+) = g^{n-2} \frac{\langle jk\rangle^4}{\langle 12\rangle\langle 23\rangle\cdots\langle n1\rangle}$$

where `⟨ij⟩ = ε_{αβ} λ_i^α λ_j^β` are holomorphic spinor brackets. This is exact at
tree level for all *n*, a stunning compression of what would otherwise be exponentially
complex Feynman diagram calculations.

---

### 1.2 MHV in N=4 SYM vs QCD vs Gravity

#### N=4 Super-Yang-Mills (still the laboratory)

N=4 SYM remains the primary arena for MHV amplitude research because:
- **Dual conformal symmetry** (Drummond–Henn–Korchemsky–Sokatchev, ~2008) constrains amplitudes tightly
- **Yangian symmetry** — infinite-dimensional symmetry of tree-level amplitudes
- **Wilson loop / amplitude duality** (Alday–Maldacena 2007): MHV amplitudes = null polygon Wilson loops
- **BDS ansatz** and remainder functions for multi-loop MHV

**2024 status:**
- 2-loop MHV amplitudes for n=6,7 are known analytically; n≥8 partially via cluster algebras
- The **momentum amplituhedron** (Damgaard–Ferro–Lukowski–Moze 2019) reformulates the
  **loop integrand** of MHV amplitudes geometrically in dual-momentum space

Recent paper: **arXiv:2407.12906** (Glew–Lukowski–Morales, July 2024)
"The Two-loop MHV Momentum Amplituhedron from Fibrations of Fibrations"
— Computes 2-loop n-point MHV integrands for N=4 SYM via a fibration-of-fibrations
construction over tree-level kinematic data.

Recent paper: **arXiv:2411.14989** (Glew et al., November 2024)
"Positive and Negative Ladders in Loop Space"
— Finds canonical forms for L-loop ladder contributions to MHV_n amplitudes;
results factorize into chiral pentagons reminiscent of the 1-loop/2-loop momentum
amplituhedron.

#### QCD

MHV amplitudes in QCD differ from N=4 SYM in two important ways:
1. **No supersymmetric cancellations** — loop corrections are harder
2. **Beyond-leading-color** contributions (subleading in 1/N_c) are computationally nontrivial

Recent work: **arXiv:PRD 2025** (Maybee et al.)
"Serendipitous syzygies of scattering amplitudes"
— Shows that MHV tree-level amplitudes and one-loop all-plus amplitudes satisfy unexpected
algebraic relations (syzygies), found using rational function representations.

For phenomenology, MHV amplitudes feed into **NLO and NNLO collider predictions** (W+jets,
Z+jets, multi-parton production). The framework is mature but new simplifications continue
to appear.

#### Gravity

Gravity MHV amplitudes have the **Hodges formula** (2012):

$$\mathcal{M}_n = \sum_{\text{spanning trees}} \prod_{(ij)\in T} \frac{[ij]^2}{\langle ij\rangle^2} \cdot (\text{Parke-Taylor factors})$$

(in the CHY / RSV representation)

Key open issue: **The double copy structure is not manifest in the Hodges formula.**
An explicit double copy construction of the Cachazo-Skinner graviton amplitude
representation has remained elusive.

Recent paper: **arXiv:2406.04539 / Comm.Math.Phys. 2025** 
"The KLT Kernel in Twistor Space"
— Shows that KLT double copy (Yang-Mills)² → gravity is present but obscured in the
graviton MHV sector. Some recent progress on color-kinematic numerators directly from
the Hodges formula.

Recent paper: **arXiv:2412.08713** (December 2024, Le Matematiche 2025)
"Uniqueness of MHV Gravity Amplitudes"
— Investigates MHV gravity amplitudes on the spinor-helicity variety. Shows they do NOT
have logarithmic singularities (unlike gluons) and do NOT admit an amplituhedron-like
construction. Instead, numerators have interesting zeroes. Conjectures uniqueness of
the numerator — new direction.

---

### 1.3 BCFW Recursion: 2024–2025 Developments

BCFW recursion (Britto-Cachazo-Feng-Witten, 2005) gives a recursion:

$$\mathcal{A}_n = \sum_{i,h} \mathcal{A}_L(\hat{1},\ldots,\hat{P}_i^h) \frac{1}{P_i^2} \mathcal{A}_R(\hat{-P}_i^{-h},\ldots,\hat{n})$$

where the shift is `|1⟩ → |1⟩ + z|n⟩`, `|n] → |n] - z|1]`.

The Parke-Taylor formula itself is a fixed point of BCFW: each channel contributes an
MHV×(3-pt anti-MHV) or (3-pt MHV)×MHV, which assembles back to the same form.

**Recent 2024–2025 applications:**

- **"2-Split of Form Factors via BCFW Recursion"** arXiv:2509.12564 (Dec 2025)
  Uses BCFW to prove 2-split factorization for form factors of composite operators
  like Tr(F²) and Tr(φ³ + (∂φ)²). Shows BCFW is still generating new results for
  non-amplitude objects.

- **"Hidden Zeros and 2-split via BCFW"** arXiv:2504.14215 (Apr 2025)
  Proves a new "hidden zero" structure in multi-particle amplitudes via BCFW induction.
  These zeros (when two adjacent momenta go collinear in a specific way) are related to
  the "2-split" factorization discovered in 2024 by Arkani-Hamed et al.

- **Boundary contributions to BCFW** remain an active topic for theories where the
  standard shift gives nonzero boundary terms (e.g., gravity, higher-derivative theories).

---

### 1.4 Amplituhedron (2024 Frontiers of Science Award)

The **amplituhedron** (Arkani-Hamed & Trnka, arXiv:1312.2007) encodes the MHV
integrand of N=4 SYM as the **volume of a geometric object** in Grassmannian space,
with no reference to locality or unitarity.

In 2024, the original paper won the **Frontiers of Science Award in Theoretical Physics**.

**Current state of the amplituhedron:**
- Proven for tree-level MHV amplitudes (Arkani-Hamed et al.)
- Momentum amplituhedron extends to dual-momentum space (works better at loop level)
- 1-loop and 2-loop: canonical forms computed (the fibrations papers above)
- **Open:** Proof of equivalence to all-loop integrand remains conjectural; extension
  to NMHV beyond n=6 is technically hard
- **Open:** Non-planar amplituhedron — there is no known geometric structure for
  non-planar (subleading color) contributions

The amplituhedron is arguably the deepest reformulation of MHV physics but requires
significant mathematical machinery (positive Grassmannians, cluster algebras, tropical
geometry).

---

### 1.5 Celestial Amplitudes and Holographic MHV

Celestial holography rewrites S-matrix elements in terms of conformal primary operators
on the 2D celestial sphere (the sky). The Mellin-transformed amplitudes become correlation
functions in a 2D CFT.

**Key 2024 results:**

- **arXiv:2403.18896** (Melton–Strominger–Wang, Mar 2024)
  "A Celestial Dual for MHV Amplitudes"
  Shows a 2D CFT = (Liouville theory c) + (level-1 Kac-Moody rank-N) + (weight -3/2
  free fermion) holographically generates 4D MHV tree-level scattering amplitudes.
  The celestial amplitudes arise in a large-N, large-c limit as leaf amplitudes.

- **arXiv:2408.10944** (Mol, Aug 2024)
  "A Holographic Construction of MHV Graviton Amplitudes in Celestial CFT"
  Extends the Melton et al. construction to generate graviton MHV amplitudes, using
  AdS₃/CFT₂ duality analytically continued to Klein space.

- **arXiv:2411.14311** (Nov 2024)
  "An AdS₃ Dual for Supersymmetric MHV Celestial Amplitudes"
  Holographic reproduction of tree-level MHV celestial amplitudes for gravitons using
  AdS₃ string theory.

**Open problems in celestial MHV:**
- The 2D CFT dual is not fully understood beyond MHV
- Extension to loop-level celestial amplitudes is largely open
- The role of the w_{1+∞} symmetry (BMS + area-preserving diffeomorphisms) at loop level

---

### 1.6 Soft Theorems and MHV

Soft theorems (Weinberg 1965, and Cachazo-Strominger subleading extensions) encode the
behavior of amplitudes when one external leg goes soft (p → 0):

$$\mathcal{A}_{n+1}(p_1,\ldots,p_n, q) \xrightarrow{q\to 0} \left[S^{(0)} + S^{(1)} + S^{(2)} + \ldots\right] \mathcal{A}_n(p_1,\ldots,p_n)$$

MHV amplitudes factorize cleanly under soft limits, and this factorization is
strongly constrained by BMS symmetry.

**Key recent result:**

- **arXiv:2505.16639 / JHEP 2025** 
  "Soft Factor Structure of MHV Amplitudes for Massless Charged Particles"
  — Shows MHV amplitudes in massless spinor and scalar QED are **fully determined by
  their soft photon behavior** and admit a factorized form. A simple derivation of
  MHV amplitudes from first principles using only soft factorization.
  **This is pedagogically important** — it shows soft data determines MHV completely.

- **arXiv:2411.13652** (Nov 2024)
  "Soft limits of gluon amplitudes in holography and cosmology"
  — Studies soft gluon limits in AdS/holographic settings, analyzes intrinsic curved-
  space contributions and their BMS-like structure.

The connection soft theorems → MHV via BMS Ward identities is well-established at tree
level. At loop level, the soft theorems receive IR-log corrections which modify the
simple factorization form.

---

### 1.7 Color-Kinematics Duality (BCJ) and MHV

The **BCJ duality** (Bern-Carrasco-Johansson, 2008) states that gauge theory amplitudes
can be written so that kinematic numerators satisfy the same Jacobi identities as color
factors. The double copy then gives gravity amplitudes.

For MHV specifically:
- At tree level, BCJ numerators for MHV are known: they take a particularly simple form
  related to the Parke-Taylor factors
- The double copy (MHV gluon)² = (MHV graviton) gives the Hodges formula

**2024–2026 developments:**

- **arXiv:2601.09815** (Jan 2026) — "The Entire Four-Graviton EFT from Color-Kinematics"
  Uses BCJ to classify all four-point graviton EFT operators from gauge-theory input.

- **Self-dual cosmology and BCJ in curved backgrounds** (JHEP Oct 2024)
  BCJ extended to de Sitter/cosmological backgrounds with a kinematic algebra.

- **Key open problem:** BCJ numerators for multi-loop MHV are not unique — there is a
  "generalized gauge freedom." Finding canonical BCJ representations for loop MHV is
  an active research problem.

---

### 1.8 Loop-Level MHV: 1-Loop and 2-Loop

**1-loop MHV in pure Yang-Mills** was first computed by Bern, Dixon, Dunbar, and Kosower
(BDDK, 1994) using **generalized unitarity cuts** — a major early success of modern
methods. The result involves rational functions of spinor brackets and box/triangle
integral functions.

The **BDS ansatz** (Bern-Dixon-Smirnov, 2005) for planar N=4 SYM proposes that
all-loop MHV amplitudes exponentiate:

$$\ln \mathcal{M}_n^{MHV} = \text{BDS}_n + R_n$$

where R_n is the **remainder function** (zero at 1-loop, nontrivial at 2-loop and beyond).

**2024 state of the art:**
- 2-loop remainder function for n=6 (hexagon): known in terms of harmonic polylogarithms
  and classical polylogarithms (Goncharov–Spradlin–Vergu–Volovich 2010)
- 2-loop n=7 remainder: known analytically via cluster coordinates
- 2-loop n≥8: only known numerically or in special kinematic limits
- The **symbol** of the remainder function is a tool for tracking its complexity
  (Goncharov–Spradlin–Vergu–Volovich 2010)

---

## 2. The Parke-Taylor Formula and Extensions

### 2.1 Original 1986 Result

**Paper:** S.J. Parke and T.R. Taylor, "Amplitude for n-Gluon Scattering"
*Phys. Rev. Lett.* 56, 2459 (1986)

The formula was initially **guessed** from low-n numerical results and then conjectured
for all n. At the time, no proof existed — it was confirmed numerically for n=5,6,7.

The formula (in modern notation):

$$\mathcal{A}_n(1^+, 2^+, \ldots, j^-, \ldots, k^-, \ldots, n^+) = i\,g^{n-2}\,(2\pi)^4\delta^4\!\left(\sum p_i\right) \frac{\langle j\,k\rangle^4}{\langle 1\,2\rangle\langle 2\,3\rangle\cdots\langle n\,1\rangle}$$

The power of the formula: it's one term, regardless of n. The corresponding Feynman
diagram expansion has (n-2)! terms.

### 2.2 First Proof: Berends-Giele Recursion (1988)

Berends and Giele proved the formula using **off-shell recursion relations** for
color-ordered partial amplitudes. The key ingredient was a recursion on the current:

$$J^\mu(1,\ldots,n) = \frac{1}{(\sum p_i)^2} \left[J^\mu_{\text{3-pt}} + J^\mu_{\text{4-pt}}\right]$$

Proof by induction on n.

### 2.3 Proof via BCFW Recursion (2005)

After BCFW was established, a clean proof emerged:

**BCFW shift:** `|1⟩ → |1̂⟩ = |1⟩ + z|n⟩`, `|n] → |n̂] = |n] - z|1]`

The shifted Parke-Taylor amplitude has only one pole in z (from channel ⟨1̂|···|n̂] = 0),
and the factorization at that pole gives:

$$A_n^{MHV}(\hat{1},2,\ldots,\hat{n}) = A_3^{\overline{MHV}}(\hat{1},2,\hat{P}) \frac{1}{P_{12}^2} A_{n-1}^{MHV}(-\hat{P},3,\ldots,\hat{n})$$

By induction (base case n=3), the Parke-Taylor formula follows. This is the cleanest
known proof.

### 2.4 Gravity MHV via KLT and Double Copy

**KLT relations** (Kawai-Lewellen-Tye, 1986) from string theory:

$$\mathcal{M}_n^{\text{grav}} = \sum_{\sigma,\tau \in S_{n-3}} A_n^{\text{YM}}(1,\sigma,n-1,n) \cdot S[\sigma|\tau] \cdot A_n^{\text{YM}}(1,\tau,n,n-1)$$

where S[σ|τ] is the **KLT kernel** (momentum-kernel matrix).

At MHV level, this simplifies considerably. For n=4 and n=5, explicit formulas are
tractable by hand. For general n, the Hodges formula (2012) gives:

$$\mathcal{M}_n^{MHV}(1^{--},2^{--},3^{++},\ldots,n^{++}) = \langle 12 \rangle^8 \cdot \det'(H)$$

where H is a (n-2)×(n-2) matrix with entries `H_{ij} = [ij]/⟨ij⟩` for i≠j,
and `H_{ii}` is determined by row-sum = 0.

**Recent:** The KLT kernel in twistor space (arXiv:2406.04539, 2024) explores
why the double copy structure is hidden in the Hodges formula — a problem that
remained open for a decade.

### 2.5 Extension to Massive Particles

For massive particles, the simple `⟨ij⟩` notation breaks down because the little group
is SU(2) rather than U(1). The **massive spinor-helicity formalism** was developed by
Arkani-Hamed, Huang, Huang (2019).

In this formalism, a massive momentum is decomposed as:
$$p^{\alpha\dot\alpha} = \lambda^\alpha \tilde\lambda^{\dot\alpha} + \mu^\alpha\tilde\mu^{\dot\alpha} = p^{\flat,\alpha\dot\alpha} + \frac{m^2}{2p\cdot q} q^{\alpha\dot\alpha}$$

There are **massive MHV amplitudes** but the Parke-Taylor simplicity is partially lost.
New results continue to appear:

- **arXiv:2501.09062** (Jan 2025) — "Massive Helicity-Chirality Spinor Formalism from
  Massless Amplitudes with On-shell Mass Insertion"
  New helicity-chirality formalism extending massive spinors.

- **arXiv:2401.09781** (Jan 2024) — "Higher-Point Gauge-Theory Couplings of Massive
  Spin-2 States"
  "MHV-like" helicity configurations for massive spin-2 states in string theory.

### 2.6 Extension to Loop Level

At loop level, the Parke-Taylor formula is modified by:
1. **UV divergences** (in non-SUSY theories): regulated by dimensional regularization
2. **IR divergences**: universal soft/collinear divergences, captured by Catani's formula
3. **Finite remainder**: theory-dependent

The **BDS ansatz** (for N=4 SYM) claims the loop-corrected log of the MHV amplitude
takes the form:
$$\ln A_n^{MHV,\ell-\text{loop}} = \text{(1-loop)} \cdot f(\epsilon) + C + O(\epsilon)$$

This is an infinite resummation. It fails at n=6 at 2-loops — the correction (remainder
function) is non-trivial and highly constrained by integrability.

---

## 3. Spinor-Helicity Formalism Tools: Cadabra2 + Python

### 3.1 What Cadabra2 Can Do

Cadabra2 is a field-theory-motivated CAS with first-class support for:
- **Index-carrying objects** with specific symmetry properties
- **Spinor algebra**: γ-matrices, σ-matrices, Fierz identities
- **Tensor canonicalization**: index symmetrization, antisymmetrization
- **Pattern matching** with implicit index contractions
- **Python integration**: can call numpy/scipy for numerical checks

Unlike Mathematica, Cadabra2 has **native spinor index handling**, making it
particularly suited for van der Waerden (2-component) spinor notation.

### 3.2 Cadabra2 for MHV: Specific Capabilities

#### Symbolic Verification of Amplitude Identities

```python
# Example: Verifying the Schouten identity in Cadabra2
# ⟨ij⟩⟨kl⟩ + ⟨ik⟩⟨lj⟩ + ⟨il⟩⟨jk⟩ = 0
# Can be encoded as:
#   eps_{ab} lambda_i^a lambda_j^b * eps_{cd} lambda_k^c lambda_l^d + ... = 0
# Cadabra2 handles the epsilon tensor contractions automatically
```

Key tasks:
- **Schouten identity checking**: εαβ εγδ has only 3 independent components in 2D →
  any product of 4 spinors satisfies the Schouten identity. Cadabra2 can verify these
  symbolically for specific or general configurations.

- **Momentum conservation**: `Σᵢ λᵢ^α λ̃ᵢ^α̇ = 0`. Cadabra2 can impose this as a
  constraint and reduce expressions.

- **Spinor completeness relations**: `λ^α λ̃^α̇ = p^{αα̇}`. Cadabra2 can substitute
  and simplify.

#### Automatic Fierz Rearrangements

A Fierz identity for Weyl spinors:
$$\chi_\alpha (\psi \bar\sigma^\mu \bar\xi) = -(\chi\sigma^\mu\bar\xi)\psi_\alpha - (\chi\sigma^\mu\bar\psi)\bar\xi_{\dot\alpha}...$$

(Martin–Moroi two-component spinor paper is the reference here.)

Cadabra2 can automate Fierz by:
1. Expressing the LHS in terms of `eps * (spinor products)`
2. Using the completeness of the Pauli matrices basis
3. Reading off coefficients

This is currently under-automated in the community — most researchers do Fierz by hand.

#### BCFW Recursion for Low n

For n=4 to n=7, BCFW can be implemented as a **recursive Python function** calling
Cadabra2 symbolic objects for the sub-amplitudes. The output can be compared against
the Parke-Taylor formula to verify correctness at each step.

This would be an excellent pedagogical paper: "Automated BCFW Recursion with Cadabra2"

### 3.3 Gaps in Available Tools

| Tool | Strength | Weakness |
|------|----------|----------|
| Cadabra2 | Spinor algebra, index handling | No built-in spinor-helicity bracket system |
| S@M (Mathematica) | Spinor brackets, numerical checks | Proprietary, Mathematica-only |
| Caravel (C++) | Multi-loop unitarity cuts | No symbolic manipulation |
| FeynCalc | Feynman diagrams, Dirac traces | No spinor-helicity notation |
| FORM | Fast polynomial reduction | Low-level, no spinor automation |

**Cadabra2's unique niche:** It's the only open-source CAS with field-theory-native
spinor index handling. Building a **spinor-helicity module for Cadabra2** would be
a genuine contribution to the tools ecosystem.

---

## 4. Low-Hanging Fruit: Specific Research Opportunities

---

### 4.1 Opportunity A: "Cadabra2 Implementation of Spinor-Helicity Formalism"

**Title suggestion:**
"A Cadabra2 Package for Spinor-Helicity Computations and Automated Verification of MHV Amplitudes"

**What it involves:**
1. Define the spinor-helicity bracket notation in Cadabra2:
   - Angle bracket ⟨ij⟩ = ε_{αβ} λᵢ^α λⱼ^β as an antisymmetric bilinear
   - Square bracket [ij] = ε_{α̇β̇} λ̃ᵢ^α̇ λ̃ⱼ^β̇
   - Momentum p_{αα̇} = λ_α λ̃_{α̇}
2. Implement: Schouten identity simplification, momentum conservation reduction
3. Implement: BCFW recursion for n=4,5,6,7 with output verified against Parke-Taylor
4. Verify spinor completeness, Fierz identities, etc.
5. Release as open-source package

**Computation involved:**
- Mostly Python/Cadabra2 implementation work
- Some physics: understanding the 2-component van der Waerden formalism (Srednicki Ch. 34-38)
- Verification against known results (Parke-Taylor, Dixon's MHV examples)

**Estimated difficulty:** ★★★ (3/5 — substantial implementation work, but the physics is clear)

**Tools/background needed:**
- Cadabra2 (Peeters, JOSS 2018)
- Srednicki Ch. 34-38 for 2-component spinors
- Elvang-Huang textbook for spinor-helicity notation
- S@M paper (arXiv:0709.1377) for comparison: this is the Mathematica version

**Similar existing papers:**
- S@M: "A Mathematica implementation of the spinor-helicity formalism" (Maitre-Mastrolia, 2007)
  — old Mathematica implementation, no Cadabra2 equivalent exists
- SpinorHelicity6D (Mathematica package) — different formalism
- caravel (C++): numerical, not symbolic

**Gap:** No Cadabra2 spinor-helicity package exists. There is no open-source symbolic
implementation in Python-based CAS with proper index handling.

**Potential venues:**
- *Computer Physics Communications* (software papers)
- *Journal of Open Source Software* (JOSS)
- *SciPost Physics Codebases*

---

### 4.2 Opportunity B: "Automated BCFW Proof of Parke-Taylor via Cadabra2"

**Title suggestion:**
"Symbolic Proof of the Parke-Taylor Formula for All n via BCFW Recursion in Cadabra2"

**What it involves:**
1. Set up BCFW recursion in Cadabra2 symbolically
2. Carry out the inductive step: assume Parke-Taylor for n-1, prove for n
3. Show the residue at the pole is exactly the product of sub-amplitudes
4. Verify the large-z behavior (boundary terms vanish) symbolically

**Computation:**
- The key algebraic step is verifying that under the shift
  `|1⟩ → |1⟩ + z|n⟩`, `|n] → |n] - z|1]`,
  the n-point Parke-Taylor has a single pole at `z = -⟨1n⟩/⟨2n⟩` (for negative
  helicities on legs j=1 and k=2), and the residue factors correctly.
- This is a symbolic manipulation in the spinor bracket algebra — exactly what Cadabra2 does.

**Estimated difficulty:** ★★ (2/5 — the physics is textbook; the implementation is the work)

**Tools needed:** Cadabra2, symbolic spinor bracket package (opportunity A)

**Similar papers:**
- The BCFW proof is in many textbooks (Elvang-Huang, Srednicki)
- No automated/computer-verified symbolic proof in Cadabra2 exists
- Most pedagogical proofs are done "by hand" for specific cases

**Novelty:** Not new physics, but:
(a) First machine-verified proof of Parke-Taylor via BCFW in an open-source CAS
(b) Template for automated proofs of other amplitude identities
(c) Highly citable as a tools paper

**Potential venues:**
- *European Physical Journal C* (Springer, tools/methods papers accepted)
- *Computers & Mathematics with Applications*
- *SciPost Physics Codebases*

---

### 4.3 Opportunity C: "MHV Amplitudes in QED via Soft Factorization — A Pedagogical Study with Symbolic Verification"

**Title suggestion:**
"MHV Amplitudes in Scalar and Spinor QED from Soft Photon Factorization: Symbolic Derivation and Verification"

**Context:**
The recent paper arXiv:2505.16639 (JHEP Oct 2025) showed MHV amplitudes in massless
QED are *fully determined* by soft photon behavior. This is a clean, elegant result
that admits a **symbolic verification** in Cadabra2.

**What it involves:**
1. Implement the QED soft factor operator in Cadabra2
2. Derive MHV amplitudes for (qq̄→Nγ) processes using only soft factorization
3. Compare to direct computation from Feynman rules
4. Extend: fermion pair production with multiple photons — can you get a QED Parke-Taylor?

**The "QED Parke-Taylor":**
For massless QED, the MHV amplitude for `f f̄ + Nγ` (fermion pair + N photons) is:

$$\mathcal{A}(p^-, q^+; k_1^+, \ldots, k_N^+) \propto \frac{\langle pq\rangle^3}{\langle p k_1\rangle \langle k_1 k_2\rangle \cdots \langle k_N q\rangle}$$

This is known but not widely discussed in the pedagogical literature.

**Estimated difficulty:** ★★★ (3/5 — need to implement the QED soft factor; result is known but derivation is fresh)

**Tools needed:** Cadabra2, basic QED knowledge, spinor-helicity package

**Similar papers:** arXiv:2505.16639 (the source); Dixon's QED amplitude paper (old)

**Potential venues:** *Physical Review D*, *JHEP*

---

### 4.4 Opportunity D: "Schouten Identity Checker and Amplitude Simplification via Gröbner Bases"

**Title suggestion:**
"Automated Simplification of Spinor Bracket Expressions Using Schouten Identities and Momentum Conservation: A Cadabra2 Approach"

**What it involves:**
The algebra of spinor brackets is:
- Antisymmetry: ⟨ij⟩ = -⟨ji⟩
- Schouten identity: ⟨ij⟩⟨kl⟩ + ⟨ik⟩⟨lj⟩ + ⟨il⟩⟨jk⟩ = 0
- Momentum conservation: Σ pᵢ = 0 → Σ λᵢ^α λ̃ᵢ^{α̇} = 0
- Spinor products: ⟨ij⟩[ji] = 2pᵢ·pⱼ = sᵢⱼ

These form a polynomial ring with relations → **Gröbner basis** methods apply.

The idea:
1. Build a Gröbner basis for the ideal generated by Schouten + momentum conservation
2. Use it to find canonical forms of spinor bracket expressions
3. Automatically simplify complex multi-term expressions to minimal form

This is related to the recent paper on "Serendipitous Syzygies" (PRD 2025) which found
unexpected relations, but that paper used a different approach (rational function
representation). A Gröbner approach might find different/new syzygies.

**Estimated difficulty:** ★★★★ (4/5 — requires computational algebraic geometry + physics)

**Tools needed:** Cadabra2 (for spinor algebra) + SageMath or Singular (for Gröbner bases)

**Similar papers:**
- "Serendipitous syzygies of scattering amplitudes" PRD 2025 (uses rational functions)
- Mathematica-based Schouten identity tools in various amplitude codes
- No Gröbner basis approach to spinor bracket simplification is published

**Potential venues:** *Physical Review D*, *Journal of Symbolic Computation*

---

### 4.5 Opportunity E: "BCFW Recursion for MHV Form Factors via Cadabra2"

**Title suggestion:**
"BCFW Recursion Relations for MHV Form Factors of Gauge-Invariant Operators: Automated Computation in Cadabra2"

**Context:**
**Form factors** are matrix elements ⟨p₁,...,pₙ|O(0)|0⟩ where O is a local operator.
They generalize amplitudes and have their own MHV sector. The recent papers on
"2-Split of Form Factors via BCFW" (arXiv:2509.12564, 2025) show this is an active area.

**What it involves:**
1. Implement BCFW recursion for MHV form factors in Cadabra2
2. Compute form factors for Tr(F²), Tr(φ²), Tr(φ³) for n=3,4,5 external gluons
3. Verify the 2-split property (factorization when adjacent momenta become collinear)
4. Compare to the analytical results of Yang and collaborators

**Estimated difficulty:** ★★★★ (4/5 — requires learning form factor technology on top of BCFW)

**Tools needed:** Cadabra2, knowledge of N=4 SYM operator content

**Similar papers:** arXiv:2509.12564, earlier Yang-type form factor papers (2012-2016)

**Potential venues:** *JHEP*, *Physical Review D*

---

### 4.6 Opportunity F: "MHV Amplitudes in de Sitter Space: Symbolic Verification"

**Title suggestion:**
"Symbolic Verification of MHV Amplitudes in the de Sitter Static Patch via Cadabra2"

**Context:**
arXiv:2105.07572 (updated Mar 2024) by Albrychiewicz and Neiman derives MHV amplitudes
and BCFW recursion for Yang-Mills in the de Sitter static patch. The computations are
done analytically but are intricate.

**What it involves:**
1. Implement the de Sitter static patch momentum variables in Cadabra2
2. Verify the MHV amplitude formulas from arXiv:2105.07572 symbolically
3. Extend to N-1 MHV from the self-dual YM solution (Section III of that paper)

**Estimated difficulty:** ★★★★ (4/5 — requires curved space background knowledge)

**Tools needed:** Cadabra2, de Sitter geometry, arXiv:2105.07572

**Similar papers:** arXiv:2105.07572 only (very new, very few citations)

**Potential venues:** *Physical Review D*, *Class. Quantum Grav.*

---

### 4.7 Opportunity G: "Quantum Algorithm for MHV Amplitudes — Symbolic Complexity Analysis"

**Title suggestion:**
"Classical vs. Quantum Complexity of n-Gluon MHV Amplitude Evaluation: A Symbolic Analysis"

**Context:**
arXiv:2507.14252 (July 2025) presents a quantum algorithm for the n-gluon MHV
scattering amplitude. The paper focuses on the quantum circuit, but **the classical
symbolic analysis of the algorithm's complexity** using Cadabra2 would be a natural companion.

**What it involves:**
1. Analyze the complexity of exact symbolic evaluation of Parke-Taylor for large n
2. Study the polynomial structure and how it scales
3. Compare: Feynman diagram complexity O((n-2)!) vs Parke-Taylor O(n) vs
   quantum algorithm O(polylog n)

**Estimated difficulty:** ★★★ (3/5 — interdisciplinary; requires quantum computation knowledge)

**Potential venues:** *Quantum*, *PRD*, *European Physical Journal C*

---

## 5. Recommended Reading Order

Below is a tiered reading list, prioritized for Ernest's current position
(post-Srednicki Ch. 34) and goal of publishable MHV work.

### Tier 0: Must Read Before Anything (Days 1–3)

| Priority | Paper | What It Covers | Link |
|----------|-------|----------------|------|
| 🔴 ESSENTIAL | Srednicki Ch. 35–38 | Weyl spinors, spinor-helicity in QFT notation | textbook |
| 🔴 ESSENTIAL | Elvang & Huang Ch. 1–3 | Spinor-helicity formalism, Parke-Taylor, BCFW | arXiv:1308.1697 |
| 🔴 ESSENTIAL | Dixon CERN Lecture 1995 | Introduction to modern amplitude methods | cds.cern.ch/record/1613349 |

### Tier 1: Core Formalism Papers (Week 1–2)

| Priority | Paper | Key Result | arXiv/DOI |
|----------|-------|-----------|-----------|
| 🔴 HIGH | Parke-Taylor 1986 | Original MHV formula | Phys.Rev.Lett. 56, 2459 |
| 🔴 HIGH | Britto-Cachazo-Feng 2004 | BCFW factorization | hep-th/0412308 |
| 🔴 HIGH | Britto-Cachazo-Feng-Witten 2005 | BCFW recursion proof | hep-th/0501052 |
| 🟡 MED | Mangano-Parke review 1991 | Multi-parton review, spinor tech | Phys.Rep. 200, 301 |
| 🟡 MED | SAGEX Review (Travaglini et al.) | 2022 comprehensive review | arXiv:2203.13011 |
| 🟡 MED | Elvang-Huang textbook | Best modern introduction | arXiv:1308.1697 |

### Tier 2: State of the Art (Week 2–4)

| Priority | Paper | Key Result | arXiv |
|----------|-------|-----------|-------|
| 🟡 MED | Arkani-Hamed–Trnka (2013) | Amplituhedron | 1312.2007 |
| 🟡 MED | Glew–Lukowski–Morales (2024) | 2-loop MHV momentum amplituhedron | 2407.12906 |
| 🟡 MED | Melton–Strominger–Wang (2024) | Celestial dual for MHV | 2403.18896 |
| 🟡 MED | 2412.08713 (Dec 2024) | Uniqueness of MHV gravity amplitudes | 2412.08713 |
| 🟢 LOW | 2505.16639 (May 2025) | MHV in QED from soft factorization | 2505.16639 |
| 🟢 LOW | KLT in twistor space (2024) | Double copy for gravity MHV | 2406.04539 |

### Tier 3: Tools & Implementation References

| Paper | What It Covers | Reference |
|-------|----------------|-----------|
| S@M Mathematica package | Spinor-helicity implementation | arXiv:0709.1377 |
| Cadabra2 JOSS paper | Cadabra2 capabilities | JOSS 2018, doi:10.21105/joss.01118 |
| Peeters Cadabra paper (2007) | Original Cadabra, spinor examples | hep-th/0701238 |
| Martin-Moroi-Wells 2-component spinor review | Van der Waerden notation | arXiv:hep-ph/9506380 |
| Learning Simplicity (2024) | ML for amplitude simplification | arXiv:2408.04720 |

### Tier 4: Advanced Topics (Month 2+)

| Paper | Topic |
|-------|-------|
| Bern-Dixon-Smirnov 2005 | BDS ansatz, loop MHV | hep-th/0505205 |
| Drummond-Henn-Korchemsky-Sokatchev | Dual conformal symmetry | arXiv:0709.2368 |
| Goncharov-Spradlin-Vergu-Volovich | Symbols and the hexagon remainder | arXiv:1006.5703 |
| Hodges 2012 | Gravity MHV formula | arXiv:1204.1930 |
| BCJ original paper | Color-kinematics duality | arXiv:0805.3993 |

---

## 6. Concrete Next Steps

Given Ernest's current position (done with Srednicki Ch. 34, working on Ch. 35–36),
here is a **concrete 6-week plan** building toward publishable MHV work in Cadabra2.

### Week 1: Spinor-Helicity Foundations in Cadabra2

**Goal:** Get comfortable with 2-component spinors in Cadabra2 notation.

**Cadabra2 exercises:**
```python
# 1. Define 2-component spinors
{lambda_{alpha}::Spinor.}
{lambda_tilde^{alphadot}::Spinor.}

# 2. Define the epsilon tensors
eps_{alpha beta}::AntiSymmetric.
eps^{alpha beta}::AntiSymmetric.

# 3. Verify: eps_{alpha gamma} eps^{gamma beta} = delta_alpha^beta
# Write this in Cadabra2 and simplify

# 4. Define angle bracket: <ij> = eps_{alpha beta} lambda_i^alpha lambda_j^beta
# Verify antisymmetry: <ij> = -<ji>
```

**Reference:** Srednicki Ch. 34–35 (just finished + working on)

**Deliverable:** A `.cnb` Cadabra2 notebook demonstrating 2-component spinor algebra
and epsilon tensor identities.

---

### Week 2: Spinor Products and the Schouten Identity

**Goal:** Implement spinor products and verify the Schouten identity symbolically.

**Key identity to verify:**
$$\langle ij\rangle\langle kl\rangle + \langle ik\rangle\langle lj\rangle + \langle il\rangle\langle jk\rangle = 0$$

**Cadabra2 approach:**
```python
# The Schouten identity follows from:
# eps_{alpha beta} eps_{gamma delta} + eps_{alpha gamma} eps_{delta beta}
#     + eps_{alpha delta} eps_{beta gamma} = 0
# (which is a 2D identity: a 2x2 determinant with 3 antisymmetric slots)

# Write this explicitly and verify using Cadabra2's epsilon contraction
```

Also implement:
- Momentum conservation constraint: `sum_i lambda_i^alpha lambda_tilde_i^alphadot = 0`
- Mandelstam variables: `s_{ij} = <ij>[ji]`

**Deliverable:** Cadabra2 notebook with Schouten identity proof and momentum conservation.

---

### Week 3: The 4-Point MHV Amplitude

**Goal:** Compute the 4-point MHV gluon amplitude A(1⁻,2⁻,3⁺,4⁺) two ways.

**Method 1: Direct from Parke-Taylor**
$$\mathcal{A}_4(1^-,2^-,3^+,4^+) = \frac{\langle 12\rangle^4}{\langle 12\rangle\langle 23\rangle\langle 34\rangle\langle 41\rangle} = \frac{\langle 12\rangle^3}{\langle 23\rangle\langle 34\rangle\langle 41\rangle}$$

**Method 2: From BCFW recursion**
Shift |3⟩ → |3⟩ + z|4⟩, |4] → |4] - z|3].

The only pole is at z₀ = -⟨34⟩/⟨14⟩ (from the channel (1+2) channel).

At the pole: A_4 = A_3(1⁻, 2⁻, P̂⁺) × (1/s₁₂) × A_3(-P̂⁻, 3⁺, 4⁺)

Verify both give the same Parke-Taylor result.

**Cadabra2 implementation:**
- Implement the BCFW shift symbolically
- Compute the factorization
- Verify algebraic equality with Parke-Taylor output

**Deliverable:** Cadabra2 notebook verifying BCFW for n=4.

---

### Week 4: 5-Point MHV and the General BCFW Proof

**Goal:** Implement the 5-point case and understand the structure of the general proof.

**The 5-point case** A(1⁻,2⁻,3⁺,4⁺,5⁺):
$$\mathcal{A}_5 = \frac{\langle 12\rangle^4}{\langle 12\rangle\langle 23\rangle\langle 34\rangle\langle 45\rangle\langle 51\rangle}$$

BCFW shift |5⟩ → |5⟩ + z|1⟩, |1] → |1] - z|5]:

Two channels contribute: (23) and (34). Verify that each channel gives a term,
and the sum equals the Parke-Taylor result using the Schouten identity.

This is the key step — the Schouten identity is essential for combining the BCFW channels.

**Deliverable:** Cadabra2 notebook verifying BCFW for n=5 with explicit Schouten identity usage.

---

### Week 5: Begin the Cadabra2 Spinor-Helicity Package

**Goal:** Start developing a reusable Cadabra2 module for spinor-helicity computations.

**Module design:**
```python
# spinor_helicity.py (Cadabra2 module)

class SpinorBracket:
    """Symbolic angle and square brackets"""
    def angle(self, i, j): ...  # <ij>
    def square(self, i, j): ...  # [ij]
    def schouten_reduce(self, expr): ...  # Apply Schouten id
    def momentum_conserve(self, expr, excluded): ...  # Apply mom cons
    def mandelstam(self, i, j): ...  # s_ij = <ij>[ji]

class BCFWRecursor:
    """Symbolic BCFW recursion"""
    def shift(self, expr, i, j): ...  # Apply [i,j> shift
    def find_poles(self, expr, z): ...  # Find poles in z
    def residue(self, expr, z, pole): ...  # Compute residue
    def recurse(self, n, helicities): ...  # Full recursion
```

**Deliverable:** First working version of `spinor_helicity.py` for Cadabra2.

---

### Week 6: Draft Paper for Opportunity A

**Goal:** Write up the Cadabra2 spinor-helicity package as a paper.

**Paper outline:**
1. Introduction: why a Cadabra2 implementation is useful
2. Cadabra2 spinor-helicity notation: mapping from standard formalism
3. Implemented algorithms: Schouten reduction, momentum conservation, BCFW
4. Verification: Parke-Taylor for n=4,5,6,7
5. Applications: symbolic checks of amplitude identities
6. Conclusions and future directions (extension to superamplitudes, NMHV)

**Target:** ~15–20 pages, submit to *Computer Physics Communications* or *SciPost Physics Codebases*

---

### Longer-Term Path to High-Impact Research

After the tools paper, the natural progression is:

```
Cadabra2 spinor-helicity package (tools paper)
    ↓
Automated symbolic proof of Parke-Taylor (pedagogy paper)
    ↓
MHV in QED from soft factorization (physics paper)
    ↓
Schouten/Gröbner syzygy approach to amplitude simplification (math-physics paper)
    ↓
One of: celestial MHV, de Sitter MHV, massive MHV, form factors (frontier paper)
```

---

## Summary: Key Papers from 2024–2026

| Paper | Topic | arXiv | Year |
|-------|-------|-------|------|
| Two-loop MHV momentum amplituhedron | MHV loop integrands | 2407.12906 | 2024 |
| Celestial dual for MHV amplitudes | Celestial holography | 2403.18896 | 2024 |
| Holographic MHV graviton amplitudes | AdS₃/MHV | 2408.10944 | 2024 |
| AdS₃ dual for SUSY MHV celestial | Celestial SUSY | 2411.14311 | 2024 |
| Uniqueness of MHV gravity amplitudes | Gravity spinor-helicity variety | 2412.08713 | 2024 |
| KLT kernel in twistor space | Double copy for gravity MHV | 2406.04539 | 2024 |
| Positive/negative ladders in loop space | Multi-loop MHV N=4 SYM | 2411.14989 | 2024 |
| Soft factor structure of MHV (QED) | MHV from soft factorization | 2505.16639 | 2025 |
| 2-Split of form factors via BCFW | BCFW for form factors | 2509.12564 | 2025 |
| Learning simplicity of amplitudes (ML) | Neural network re-derives PT formula | 2408.04720 | 2024/25 |
| Serendipitous syzygies of amplitudes | Unexpected amplitude relations | PRD 2025 | 2025 |
| Four-graviton EFT from color-kinematics | BCJ for gravity EFT | 2601.09815 | 2026 |
| Quantum algorithm for n-gluon MHV | QC applied to Parke-Taylor | 2507.14252 | 2025 |

---

## Honest Assessment of the Research Landscape

**What's saturated:**
- The Parke-Taylor formula is 40 years old and proven many ways. A new proof alone has
  no novelty value unless it reveals new structure.
- Basic MHV in N=4 SYM at tree level is extremely well-understood.
- Simple BCFW recursion is textbook material.

**What's genuinely open (as of Feb 2026):**
- The explicit double copy structure of the Hodges (gravity MHV) formula
- Uniqueness of MHV gravity amplitude numerators (arXiv:2412.08713 — very recent)
- Celestial MHV at loop level
- Non-planar amplituhedron
- MHV in curved spacetime beyond tree level
- **Automated symbolic tools** — there is no open-source, Python-based, proper CAS
  implementation of spinor-helicity computations. This is a **genuine gap**.

**Best strategy for Ernest:**
Start with the tools paper (Opportunity A). It is:
1. Definitely publishable (tools papers are highly cited)
2. Builds all needed skills for subsequent physics papers
3. Creates infrastructure that makes all other opportunities easier
4. Fits perfectly with his current position in Srednicki Ch. 34–36

The physics papers (Opportunities B–G) will be significantly more tractable once the
Cadabra2 spinor-helicity infrastructure is in place.

---

*Survey compiled by Cyclonus (OpenClaw AI) — February 26, 2026.*
*Sources: 15+ arxiv papers from 2024–2026, web searches across Brave search API.*
*All arxiv links are in the format `arxiv.org/abs/XXXXXXXXX` for easy lookup.*
