# Spectral theory for the mass gap: a study route for a differential geometer

Author: strategist (Claude Opus 5.5), 2026-11-16, for Ernest. Question: since spec(H) is the heart of the mass gap, should we study functional analysis or related pure mathematics, coming from differential geometry?

## 1. Answer

Yes, but a targeted slice, not a course. Five ideas carry almost all the weight, and two of them are differential geometry in disguise (Laplacians on compact Lie groups, Peter–Weyl). The Hamiltonian is **not a matrix**. On a finite lattice the transfer operator is an integral operator with a continuous kernel on a compact group manifold: a "matrix with continuous indices", compact, with discrete spectrum. In infinite volume and in the continuum, H is an unbounded self-adjoint operator whose spectrum is mostly continuous. The mass gap is the statement that this continuous spectrum stays a fixed distance Δ above the isolated eigenvalue 0. Finite-dimensional matrix theory is the warm-up; the real objects need the operator versions below.

Skip for now: general Banach-space theory, distributions beyond what the Lean statement uses, operator algebras, Kato's perturbation theory (later, if we need continuity of gaps).

## 2. The route

Each stage names what to learn, why this program needs it, and where it already appears in our work. Book chapter numbers are given only where checked (see §5).

| Stage | Learn | Why we need it | Where it bites in our work |
|---|---|---|---|
| 0 (known) | Riemannian geometry, Laplace–Beltrami, compact Lie groups, bi-invariant metrics (Hall GTM 222, Frankel, Lee; already in the mono-notes) | The link variables live on a compact Lie group | `su2-so3-covering.tex`, `quaternions-spin.tex` |
| 1 | Hilbert spaces L²(X, μ); bounded operators and adjoints; compact operators; Hilbert–Schmidt integral operators; spectral theorem for compact self-adjoint operators; trace class; **Mercer's theorem** | The finite-lattice transfer operator T on L²(G^{spatial links}, Haar) is of exactly this type, and Z = Tr Tⁿ needs trace class | Mono-notes Eqs. `ym:transfer-kernel`, `ym:trace-formula`; Prop. `ym:prop:purity` |
| 2 | Positivity: Perron–Frobenius for matrices, then **Jentzsch's theorem** (Krein–Rutman) for positive kernels; positive-definite functions on compact groups; **Peter–Weyl**; convolution by a class function is diagonal in characters (Schur orthogonality) | Gives T ≥ 0, a unique vacuum and a gap **on every finite lattice**; shows precisely what the Clay problem adds (uniformity) | §3 below; the open hypothesis of Prop. `ym:prop:defect` |
| 3 | Unbounded self-adjoint operators; spectral theorem with projection-valued measures; Stone's theorem; semigroups e^{−tH}; point vs continuous spectrum; isolated eigenvalues | This is the language of the theorem: the Lean's `EnergyMomentumMeasure` is a projection-valued measure on ℝ⁴, and (W0) is a spectral condition | Lean `PublicTheory.lean:180`; mono-notes Prop. `ym:prop:gap-decay`, Lemma `ym:lem:layer`, Eq. `ym:massless-rate` |
| 4 | Euclidean side: reflection positivity, Osterwalder–Schrader reconstruction, transfer matrices of lattice gauge theory | How a Euclidean lattice measure produces a Hilbert space and a positive H at all | Lean units V.1–V.3 (Armstrong's map); packet YM-R08 |
| 5 (later) | Cluster and polymer expansions (Kotecký–Preiss); multiscale analysis | The analytic engine of constructive QFT; needed to read OpenAI's ultraviolet part or any strong-coupling argument | Roadmap T12; route R3 |

## 3. Worked example: where stages 1–2 pay off immediately

Mono-notes Prop. `ym:prop:defect` assumes the transfer operator is nonnegative in every direction. Gattringer–Lang (p. 21) cites Lüscher, *Commun. Math. Phys.* 54, 283 (1977), "Construction of a selfadjoint, strictly positive transfer matrix for Euclidean lattice gauge theories". The following sketch shows how the hypothesis follows for the Wilson action. **Label: strategist's argument, unreviewed; Lüscher's paper not read; the four classical theorems are quoted from memory.** It is posted for the G02 reviewer (tracker local:13).

Notation of the notes: spatial configurations U, U′ on one time slice; temporal links V; the kernel 𝖳(U,U′) = e^{−βE_s(U)/2} K(U,U′) e^{−βE_s(U′)/2} with K(U,U′) = ∫ dV ∏_{x,j} k(U(x,j) V(x+e_j) U′(x,j)⁻¹ V(x)⁻¹) and k(g) = exp(−β(1 − (1/m) Re tr g)).

1. **Gauge average.** With (U′)^V(x,j) = V(x) U′(x,j) V(x+e_j)⁻¹, the argument of k is U(x,j) ((U′)^V(x,j))⁻¹. So K = K₀P, where K₀(U,W) = ∏_ℓ k(U_ℓ W_ℓ⁻¹) and (Pf)(W) = ∫ dV f(W^V). Haar invariance makes P an orthogonal projection (onto gauge-invariant functions). Since k is a class function, K₀(U^γ, W^γ) = K₀(U,W), so K₀ commutes with P. Hence K = PK₀P.
2. **One link.** Convolution by a class function f = Σ_λ a_λ χ_λ acts on the Peter–Weyl block V_λ ⊗ V_λ* of L²(G) by the scalar a_λ/d_λ (Schur orthogonality). So it is a positive operator iff all a_λ ≥ 0.
3. **The Wilson weight has a_λ ≥ 0.** e^{(β/m) Re tr g} = Σ_n (β/m)ⁿ (Re tr g)ⁿ/n!, and Re tr g = (χ_F + χ_F̄)/2. Products of characters are characters of tensor products, which decompose with nonnegative integer multiplicities. So every a_λ ≥ 0 for β ≥ 0. Checked numerically for SU(2): a_j = 2(2j+1) I_{2j+1}(β)/β > 0 by quadrature and closed form, β ∈ {0.5, 2, 8}, j ≤ 3 (`ledger/evidence/transfer-positivity/check_su2_character_coefficients.py`).
4. **Assemble.** K₀ is a tensor product of positive operators, so K₀ ≥ 0; then K = PK₀P ≥ 0; then 𝖳 = MKM ≥ 0 with M the positive multiplication operator e^{−βE_s/2}. The kernel is continuous on a compact space, so by **Mercer** T is trace class and Tr Tⁿ = ∫ 𝖳ⁿ(U,U) dU = Z, as used in `ym:trace-formula`. By the symmetry of the hypercubic torus, the same holds in every direction.
5. **Finite-volume gap.** 𝖳 is strictly positive pointwise (an integral of a positive function). **Jentzsch:** the largest eigenvalue λ₀ is simple, its eigenfunction is positive, and every other eigenvalue has |λ| < λ₀. So on *every* finite torus there is a unique vacuum and a gap log(λ₀/λ₁) > 0.

The lesson of step 5 is the most important point for orientation: **a gap on each finite lattice is automatic.** The whole difficulty of the mass gap is that log(λ₀/λ₁), measured in physical units, must stay bounded below as the volume grows and the lattice spacing shrinks. Spectral theory tells you exactly where the problem is; it does not solve it.

## 4. The differential-geometry bridge

In the Hamiltonian form of lattice gauge theory (Kogut–Susskind 1975, from memory, unverified), the Hilbert space for one link is L²(G), and the "electric" energy of a link is a multiple of minus the Laplace–Beltrami operator of the bi-invariant metric on G. Its spectrum is given by Peter–Weyl: on V_λ ⊗ V_λ* the Laplacian acts by the Casimir eigenvalue of λ. For G = SU(2) ≅ S³ with the round metric of radius 1, the eigenvalues of −Δ are ℓ(ℓ+2) with multiplicity (ℓ+1)², and ℓ = 2j matches V_j ⊗ V_j* of dimension (2j+1)².

At strong coupling the electric term dominates: the ground state is the constant function on every link, and the cheapest gauge-invariant excitation is one plaquette loop in the fundamental representation. Its energy is a positive multiple of the fundamental Casimir eigenvalue, so the gap is positive and is spectral geometry of a compact Lie group. At weak coupling the magnetic (plaquette) term dominates and no such picture survives. The open problem lives between these regimes. Step 3 above (positivity of a_λ) is the Euclidean twin of the same mathematics: for the heat-kernel action, where k is the heat kernel of G, a_λ = d_λ e^{−t·Casimir(λ)} and positivity is immediate.

## 5. Books, in the order to use them

| Book | Use | Verified |
|---|---|---|
| B. C. Hall, *Quantum Theory for Mathematicians*, GTM 267 (Springer 2013) | Stages 1 and 3, written for mathematicians by the author of the Lie-groups book you already use. Spectral theory is chs. 6–10: bounded spectral theorem, then unbounded self-adjoint operators (ch. 9) and their spectral theorem (ch. 10). Ends with Hamiltonian mechanics on manifolds and geometric quantization | Chapter titles from library and publisher listings (web search 2026-11-16); not on disk |
| M. Reed, B. Simon, *Methods of Modern Mathematical Physics I: Functional Analysis* | Reference for stage 1–3 facts (compact operators, trace class, spectral theorem). Cited by Gattringer–Lang for the continuous spectrum | Citation verified in G–L ref. [6] of ch. 1; contents not checked; not on disk |
| J. Glimm, A. Jaffe, *Quantum Physics: A Functional Integral Point of View*, 2nd ed. (1987) | Stage 4: Osterwalder–Schrader reconstruction, reflection positivity, transfer matrices, cluster expansions | Cited by G–L as the textbook for OS reconstruction; contents not checked; not on disk |
| E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, Lect. Notes Phys. 159 (1982) | The closest rigorous text to our problem: lattice gauge theory, positivity, strong-coupling results | Citation verified in G–L ref. [14]; contents not checked; not on disk |
| M. Lüscher, *Commun. Math. Phys.* 54, 283–292 (1977) | The positivity of the transfer matrix (§3 above) | Bibliographic data confirmed by web search; full text not read |
| K. Osterwalder, E. Seiler, *Ann. Phys.* 110, 440 (1978) | Reflection positivity of the Wilson action; rigorous strong-coupling confinement | Citation as given in G–L refs.; volume and page not independently confirmed |
| Gattringer–Lang (on disk, being parsed) | Physics context only: §1.1 for the transfer matrix in quantum mechanics, §6.3 for measuring masses from correlators | Parsed text layer; see `notes/gattringer-lang-evaluation.md` |

## 6. Exercises that move the program

Each is small, checkable, and feeds a packet.

1. Peter–Weyl on SU(2): derive the spectrum of −Δ on S³ from V_j ⊗ V_j* (feeds YM-T09's SU(2) examples).
2. Prove step 2 of §3 (convolution by a class function is a_λ/d_λ on each block).
3. Write steps 1–5 of §3 as a proposition with proof for the mono-notes, replacing the hypothesis in Prop. `ym:prop:defect`, **after** the G02 reviewer has checked §3.
4. Truncate the SU(2) character basis on a 2×2×2 spatial torus and compute the top two eigenvalues of T numerically as β varies (feeds YM-N03 and T10): watch the finite-volume gap shrink as β grows.
5. Identify which facts of stage 3 the proof of Prop. `ym:prop:gap-decay` uses (answer: spectral theorem for H ≥ 0, the spectral measure of a vector, Tonelli).
