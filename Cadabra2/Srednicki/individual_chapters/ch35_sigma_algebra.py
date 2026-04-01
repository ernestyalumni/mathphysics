"""
ch35_sigma_algebra.py
======================
Srednicki QFT — Chapter 35: Manipulating Spinor Indices

What this file covers (section by section):
  §35.A  Setup — ε symbols and σ^μ matrices (recap, stage for Ch. 35)
  §35.B  Key identity: σ^μ_{aȧ} σ_{μ bḃ} = -2 ε_{ab} ε_{ȧḃ}   [eq 35.4]
  §35.C  Key identity: ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ} = -2 g^{μν}  [eq 35.5]
  §35.D  σ̄^μ — definition and numerical verification  [eq 35.19–35.20]
  §35.E  Generators in terms of σ, σ̄:
           (S^μν_L)_a^b = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)_a^b   [eq 35.21]
           (S^μν_R)^ȧ_ḃ = -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ)^ȧ_ḃ   [eq 35.22]
  §35.F  Clifford-like algebra: σ^μ σ̄^ν + σ^ν σ̄^μ = -2 g^{μν} I₂
  §35.G  Index-free notation — Cadabra2 symbolic bilinears
  §35.H  Lorentz vector bilinear ψ†σ̄^μ χ and its hermitian conjugate

Run with:
    python3 ch35_sigma_algebra.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70
sec = lambda s: print(f"\n{SEP}\n  {s}\n{SEP}")

print(SEP)
print("  Srednicki Ch. 35 — Manipulating Spinor Indices")
print(SEP)

# =============================================================================
# §35.A  SETUP — ε SYMBOLS AND σ^μ MATRICES
# =============================================================================
# Chapter 35 opens by recalling the ε conventions established in Ch. 34,
# and the σ^μ invariant symbol.  We reproduce them here for completeness.
#
# ε normalization (eq. 35.1):
#   ε^{12} = ε^{1̇2̇} = ε_{21} = ε_{2̇1̇} = +1
#   ε^{21} = ε^{2̇1̇} = ε_{12} = ε_{1̇2̇} = -1
#
# IMPORTANT SIGN RULE:
#   We raise/lower indices by contracting on the SECOND index of ε.
#   Contracting the first index introduces an extra minus sign.
#   This is because ε is antisymmetric: ε^{ab} = -ε^{ba}.
#
# σ^μ invariant symbol (eq. 35.2):
#   σ^μ_{aȧ} = (I, σ_1, σ_2, σ_3)
#   This lives in the (2,2) representation — it has one undotted (left-handed)
#   and one dotted (right-handed) spinor index, plus a spacetime vector index.

sec("§35.A  Setup — ε symbols, σ^μ, Pauli matrices")

# Declare index sets for Cadabra2
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"),       Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"),        Ex(r"position=free"))

cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

print("Index conventions:")
print("  Undotted (left-handed, (2,1)):  α β γ δ")
print("  Dotted   (right-handed, (1,2)): α̇ β̇ γ̇ δ̇")
print("  Spacetime:                       μ ν ρ σ")

# Metric g^{μν} = diag(-1,+1,+1,+1)  [Srednicki Eq. 1.8 / 2.4, mostly-plus -+++]
g = np.diag([-1., 1., 1., 1.])

# ε symbols: index order is row=first index, col=second index
# ε_{12} = -1, ε_{21} = +1  (since ε_{12} = ε^{21} * determinant = ... )
# Srednicki: ε^{12} = +1, ε_{21} = +1; ε^{21} = -1, ε_{12} = -1
eps_lower = np.array([[0, -1],   # ε_{ab}: ε_{12}=-1, ε_{21}=+1
                       [1,  0]], dtype=complex)
eps_upper = np.array([[0,  1],   # ε^{ab}: ε^{12}=+1, ε^{21}=-1
                       [-1, 0]], dtype=complex)
# Verify ε_{ab} ε^{bc} = δ_a^c
assert np.allclose(eps_lower @ eps_upper, np.eye(2)), "ε normalization failed!"
print("\nε_{ab} @ ε^{bc} = identity: ✓")
print(f"  ε_{{ab}} = {eps_lower.real.astype(int)}")
print(f"  ε^{{ab}} = {eps_upper.real.astype(int)}")

# Pauli matrices (eq. 35.3)
sigma = {
    1: np.array([[0, 1],  [1, 0]],  dtype=complex),   # σ_1
    2: np.array([[0, -1j],[1j, 0]], dtype=complex),    # σ_2
    3: np.array([[1, 0],  [0, -1]], dtype=complex),    # σ_3
}
I2 = np.eye(2, dtype=complex)

# σ^μ_{aȧ} = (I, σ_1, σ_2, σ_3)
sigma_vec = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
print("\nσ^μ matrices (μ=0,1,2,3):")
for mu in range(4):
    print(f"  σ^{mu} =\n    {sigma_vec[mu][0]}\n    {sigma_vec[mu][1]}")

# =============================================================================
# §35.B  KEY IDENTITY: σ^μ_{aȧ} σ_{μ bḃ} = -2 ε_{ab} ε_{ȧḃ}   [eq 35.4]
# =============================================================================
# Physical meaning:
#   This is the fundamental σ-completeness relation.  It says that the
#   contraction of two σ^μ symbols over the spacetime index produces
#   the SL(2,C) invariant ε_{ab} ε_{ȧḃ} — the "metric" of (2,1)⊗(1,2).
#
#   The factor of -2 arises from:
#     Σ_μ g^{μμ} tr(σ_μ² ) = (+1)(2) + (-1)(2) + (-1)(2) + (-1)(2) = -4
#   distributed among the 4 components of each side.
#
# Index structure:
#   LHS: σ^μ_{aȧ} g_{μν} σ^ν_{bḃ}  — contraction over BOTH μ and ν via g_{μν}
#   The notation σ_{μ bḃ} means σ^ν_{bḃ} with ν lowered by g_{νμ},
#   i.e. σ_{0 bḃ} = +σ^0_{bḃ}, σ_{i bḃ} = -σ^i_{bḃ} for i=1,2,3.
#
# Verification plan: compute (σ^μ g_{μν} σ^ν)_{aȧ,bḃ} for all 16 index
# combinations and compare to -2 ε_{ab} ε_{ȧḃ}.

sec("§35.B  σ^μ_{aȧ} σ_{μ bḃ} = -2 ε_{ab} ε_{ȧḃ}  [eq 35.4]")

print("Computing σ^μ_{aȧ} g_{μν} σ^ν_{bḃ} for all (a,ȧ,b,ḃ) ∈ {1,2}⁴...")
print("This gives a 4-index tensor with (2×2)×(2×2) = 16 components.\n")

# LHS: contract σ^μ_{aȧ} with g_{μν} σ^ν_{bḃ}
# = Σ_μ g_{μμ} (σ^μ)_{aȧ} (σ^μ)_{bḃ}   (metric is diagonal)
lhs_35_4 = np.zeros((2, 2, 2, 2), dtype=complex)  # indices: [a, adot, b, bdot]
for mu in range(4):
    lhs_35_4 += g[mu, mu] * np.einsum('ij,kl->ijkl',
                                       sigma_vec[mu], sigma_vec[mu])

# RHS: -2 ε_{ab} ε_{ȧḃ}
rhs_35_4 = np.zeros((2, 2, 2, 2), dtype=complex)
for a in range(2):
    for b in range(2):
        for adot in range(2):
            for bdot in range(2):
                rhs_35_4[a, adot, b, bdot] = -2 * eps_lower[a, b] * eps_lower[adot, bdot]

err_35_4 = np.max(np.abs(lhs_35_4 - rhs_35_4))
print("  Verifying all 16 components of σ^μ_{aȧ} σ_{μ bḃ} = -2 ε_{ab} ε_{ȧḃ}:")
print(f"  Max error: {err_35_4:.2e}")
if err_35_4 < 1e-12:
    print("  ✓ eq. 35.4 verified for all 16 (aȧ,bḃ) component pairs!")
else:
    print("  ✗ MISMATCH!")

# Print a few representative components to build intuition
print("\n  Sample components (a,ȧ,b,ḃ) with 0-indexed a,b ∈ {0=1, 1=2}:")
print("  (a,ȧ,b,ḃ)  LHS              RHS")
for a in range(2):
    for adot in range(2):
        for b in range(2):
            for bdot in range(2):
                lv = lhs_35_4[a, adot, b, bdot]
                rv = rhs_35_4[a, adot, b, bdot]
                if abs(lv) > 1e-10 or abs(rv) > 1e-10:
                    print(f"  ({a+1},{adot+1},{b+1},{bdot+1})     "
                          f"{lv.real:+.1f}    {rv.real:+.1f}")

# =============================================================================
# §35.C  KEY IDENTITY: ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ} = -2 g^{μν}  [eq 35.5]
# =============================================================================
# Physical meaning:
#   This is the inverse of eq. 35.4 in spirit — it says that contracting two
#   σ^μ symbols with TWO ε^{ab} factors projects out the spacetime vector
#   content, recovering the metric g^{μν}.
#
#   Equivalently: the σ matrices form an OVERCOMPLETE frame for
#   the space of 2×2 Hermitian matrices.  The four σ^μ (for μ=0,1,2,3)
#   provide a basis, and g^{μν} is their "Gram matrix" under this inner product.
#
#   This identity is the foundation for the proof that
#     A^μ = -½ ε^{ab} ε^{ȧḃ} σ^μ_{bḃ} A_{aȧ}   (inverse of A_{aȧ} = σ^μ_{aȧ} A_μ)
#   so that σ^μ is genuinely invertible as a basis change.
#
# Computation:
#   LHS = ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ}
#   In matrix language: ε^T σ^μ ε σ^ν contracted appropriately.
#   We compute this as:
#     Σ_{a,ȧ,b,ḃ} ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ}
#   for each (μ,ν) pair.

sec("§35.C  ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ} = -2 g^{μν}  [eq 35.5]")

print("Computing ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ} for all (μ,ν)...")

lhs_35_5 = np.zeros((4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        # Contract: ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ}
        # = Σ_{a,b} ε^{ab} Σ_{ȧ,ḃ} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ}
        # = Σ_{a,b} ε^{ab} [eps_upper ⊙ σ^μ] ... let's do it explicitly
        val = 0.0 + 0.0j
        for a in range(2):
            for adot in range(2):
                for b in range(2):
                    for bdot in range(2):
                        val += (eps_upper[a, b] * eps_upper[adot, bdot]
                                * sigma_vec[mu][a, adot]
                                * sigma_vec[nu][b, bdot])
        lhs_35_5[mu, nu] = val

rhs_35_5 = -2 * g  # -2 g^{μν}

err_35_5 = np.max(np.abs(lhs_35_5 - rhs_35_5))
print(f"\n  ε^{{ab}} ε^{{ȧḃ}} σ^μ_{{aȧ}} σ^ν_{{bḃ}} =")
for row in lhs_35_5.real:
    print("   ", row)
print(f"\n  -2 g^{{μν}} =")
for row in rhs_35_5:
    print("   ", row)
print(f"\n  Max error: {err_35_5:.2e}")
if err_35_5 < 1e-12:
    print("  ✓ eq. 35.5 verified for all 16 (μ,ν) components!")
else:
    print("  ✗ MISMATCH!")

# Also verify: tr(σ^μ σ̄^ν) = 2 g^{μν} (this follows from 35.5 via σ̄ def)
# We'll confirm it after defining σ̄ in §35.D.

# =============================================================================
# §35.D  σ̄^μ — DEFINITION AND NUMERICAL VERIFICATION  [eq 35.19-35.20]
# =============================================================================
# Definition (eq. 35.19):
#   σ̄^{μ ȧa} ≡ ε^{ab} ε^{ȧḃ} σ^μ_{bḃ}
#
# Physical meaning:
#   σ̄^μ is the RAISED-INDEX version of σ^μ.  It has both spinor indices
#   upstairs (one dotted, one undotted).  While σ^μ maps vectors to the
#   (2,2) tensor of type (aȧ), σ̄^μ maps vectors to the type (ȧa).
#
#   The ε symbols perform the index raising:
#     ε^{ab}: raises the undotted index b → a (with sign)
#     ε^{ȧḃ}: raises the dotted index ḃ → ȧ (with sign)
#
#   The result is that σ̄^μ = (I, -σ⃗).  The spatial components flip sign
#   because ε σ_i ε = -σ_i (a Pauli matrix identity).
#
# Numerical result (eq. 35.20):
#   σ̄^{μ ȧa} = (I, -σ_1, -σ_2, -σ_3)
#
# PROOF of eq. 35.20 from eq. 35.19:
#   For μ=0: σ̄^{0 ȧa} = ε^{ab} ε^{ȧḃ} δ_{bḃ}
#            = ε^{ab} ε^{ȧḃ} δ_{bḃ} = (ε ε^T)^{ȧa} = I  (since ε ε^T = -I... )
#   Actually, more carefully:
#   σ̄^{μ ȧa} = ε^{ab} ε^{ȧḃ} σ^μ_{bḃ}
#   In matrix form: (ε^T)^{aα} σ^μ_{αβ} (ε)^{βȧ}
#   = (ε^T σ^μ ε)^{aȧ}
#   For μ=0: ε^T I ε = ε^T ε  — but ε^T = -ε (antisymmetric), so ε^T ε = -ε² = I
#   For μ=i: ε^T σ_i ε = -ε σ_i ε
#   Use the identity ε σ_i ε^{-1} = -σ_i^T (standard SU(2) property):
#     ε σ_1 ε^{-1} = -σ_1^T = -σ_1   (σ_1 is symmetric)
#     ε σ_2 ε^{-1} = -σ_2^T = +σ_2   (σ_2 is antisymmetric → -σ_2^T = +σ_2)
#     ε σ_3 ε^{-1} = -σ_3^T = -σ_3   (σ_3 is symmetric)
#   Wait, we need to be careful. Let me just compute numerically.

sec("§35.D  σ̄^{μ ȧa} = ε^{ab} ε^{ȧḃ} σ^μ_{bḃ}  [eq 35.19-35.20]")

# Compute σ̄^{μ ȧa} = ε^{ab} ε^{ȧḃ} σ^μ_{bḃ}
# Index placement: σ̄^{μ}[adot, a] = Σ_{b,bdot} ε^{ab} ε^{ȧḃ} σ^μ[b, bdot]
sigmabar_computed = {}
for mu in range(4):
    mat = np.zeros((2, 2), dtype=complex)  # [adot, a]
    for adot in range(2):
        for a in range(2):
            val = 0.0 + 0.0j
            for b in range(2):
                for bdot in range(2):
                    val += eps_upper[a, b] * eps_upper[adot, bdot] * sigma_vec[mu][b, bdot]
            mat[adot, a] = val
    sigmabar_computed[mu] = mat

# Expected: σ̄^μ = (I, -σ_1, -σ_2, -σ_3)
sigmabar_expected = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}

print("Computed σ̄^{μ ȧa} = ε^{ab} ε^{ȧḃ} σ^μ_{bḃ}:")
max_err_sigmabar = 0.
for mu in range(4):
    err = np.max(np.abs(sigmabar_computed[mu] - sigmabar_expected[mu]))
    max_err_sigmabar = max(max_err_sigmabar, err)
    print(f"\n  σ̄^{mu} (computed)  =\n    {sigmabar_computed[mu][0]}\n    {sigmabar_computed[mu][1]}")
    print(f"  σ̄^{mu} (expected)  =\n    {sigmabar_expected[mu][0]}\n    {sigmabar_expected[mu][1]}")
    print(f"  err = {err:.2e}")

sigmabar_vec = sigmabar_expected  # Use for subsequent calculations

print(f"\nMax error across all μ: {max_err_sigmabar:.2e}")
if max_err_sigmabar < 1e-12:
    print("✓ σ̄^{μ ȧa} = (I, -σ⃗) verified! [eq 35.20]")
else:
    print("✗ MISMATCH in σ̄ computation!")

# Now verify tr(σ^μ σ̄^ν) = -2 g^{μν}   [Srednicki mostly-plus g=diag(-1,+1,+1,+1)]
# tr(σ^μ σ̄^ν) = σ^μ_{aȧ} σ̄^{ν ȧa} = Tr_matrix[σ^μ · σ̄^ν]
# where the matrix product contracts ȧ (shared column/row index)
# and the trace contracts a.
# With g=diag(-1,+1,+1,+1): -2g^{00}=+2, -2g^{ii}=-2 → diag(+2,-2,-2,-2)
print("\nChecking tr(σ^μ σ̄^ν) = -2 g^{μν}  [Srednicki mostly-plus]:")
trace_mat = np.zeros((4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        # σ^μ has indices [a, adot], σ̄^ν has indices [adot, a]
        # Product (σ^μ · σ̄^ν)_a^b = σ^μ_{aȧ} σ̄^{ν ȧb}
        # Trace = δ_a^b → sum over a
        trace_mat[mu, nu] = np.trace(sigma_vec[mu] @ sigmabar_vec[nu])

err_trace = np.max(np.abs(trace_mat - (-2 * g)))
print("  tr(σ^μ σ̄^ν) =")
for row in trace_mat.real:
    print("   ", row)
print(f"  Max error vs -2 g^{{μν}}: {err_trace:.2e}")
if err_trace < 1e-12:
    print("  ✓ tr(σ^μ σ̄^ν) = -2 g^{μν} verified!")

# =============================================================================
# §35.E  GENERATORS IN TERMS OF σ, σ̄  [eq 35.21-35.22]
# =============================================================================
# Background — how we arrive at these formulas:
#   In Ch. 34, we derived S^μν_L from the explicit Lorentz transformation
#   properties of ψ_a using Pauli matrices directly.
#   In Ch. 35, we derive them from the invariance of σ^μ itself (eq. 35.9):
#
#     σ^ρ_{aȧ} = Λ^ρ_τ L(Λ)_a^b R(Λ)_ȧ^ḃ σ^τ_{bḃ}
#
#   Taking the infinitesimal version and projecting with ε^{ȧċ} or ε^{ac}
#   yields (eqs. 35.17-35.18):
#
#     (S^μν_L)_{ac} = (i/4) ε^{ȧċ} (σ^μ_{aȧ} σ^ν_{cċ} - σ^ν_{aȧ} σ^μ_{cċ})
#     (S^μν_R)_{ȧċ} = (i/4) ε^{ac} (σ^μ_{aȧ} σ^ν_{cċ} - σ^ν_{aȧ} σ^μ_{cċ})
#
#   Written using σ̄^μ to raise indices, these become the compact forms:
#
#     (S^μν_L)_a^b = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)_a^b   [eq 35.21]
#     (S^μν_R)^ȧ_ḃ = -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ)^ȧ_ḃ  [eq 35.22]
#
#   Matrix multiplication convention:
#     σ^μ has row=a (undotted), col=ȧ (dotted)
#     σ̄^ν has row=ȧ (dotted),  col=a (undotted)  [after index raising]
#   So σ^μ · σ̄^ν is a 2×2 matrix with both indices undotted: (a, b)  ← left-handed
#      σ̄^μ · σ^ν is a 2×2 matrix with both indices dotted:   (ȧ, ḃ) ← right-handed
#
# KEY INSIGHT: The "missing indices" convention in Ch. 35:
#   In (S^μν_L)_a^b = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)_a^b,
#   the dotted index ȧ is contracted (summed over):
#     (σ^μ σ̄^ν)_a^b = σ^μ_{aȧ} σ̄^{ν ȧb}
#   This suppressed pair is written as {_ȧ}^{ȧ} → standard matrix product.

sec("§35.E  Generators S^μν_L, S^μν_R in terms of σ, σ̄  [eq 35.21-35.22]")

# Build S^μν_L from Ch. 34 formula (for comparison)
def eps3(i, j, k):
    """3D Levi-Civita symbol."""
    return int(np.linalg.det([[i==1,i==2,i==3],[j==1,j==2,j==3],[k==1,k==2,k==3]]))

S_L = [[None]*4 for _ in range(4)]
for i in range(1, 4):
    for j in range(1, 4):
        mat = np.zeros((2, 2), dtype=complex)
        for k in range(1, 4):
            mat += eps3(i, j, k) * sigma[k]
        S_L[i][j] = 0.5 * mat
for k in range(1, 4):
    S_L[k][0] = (1j / 2) * sigma[k]
    S_L[0][k] = -(1j / 2) * sigma[k]
S_L[0][0] = np.zeros((2, 2), dtype=complex)
for i in range(1, 4):
    S_L[i][i] = np.zeros((2, 2), dtype=complex)

# Compute S^μν_L from Ch. 35 formula: (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)
# Matrix product: (σ^μ σ̄^ν)_{ab} = σ^μ_{aȧ} σ̄^{ν ȧb} (standard matmul)
S_L_ch35 = [[None]*4 for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        # σ^μ: [a, adot] indices; σ̄^ν: [adot, b] indices
        # standard matrix product contracts adot
        prod = (1j / 4) * (sigma_vec[mu] @ sigmabar_vec[nu]
                           - sigma_vec[nu] @ sigmabar_vec[mu])
        S_L_ch35[mu][nu] = prod

# Compute S^μν_R from Ch. 35 formula: -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ)
# σ̄^μ: [adot, a] indices; σ^ν: [a, bdot] indices → product: [adot, bdot]
S_R_ch35 = [[None]*4 for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        prod = -(1j / 4) * (sigmabar_vec[mu] @ sigma_vec[nu]
                            - sigmabar_vec[nu] @ sigma_vec[mu])
        S_R_ch35[mu][nu] = prod

# Verify S^μν_L (Ch. 35 formula) matches Ch. 34 definition
print("Verifying (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ) = S^μν_L from Ch. 34:")
max_err_SL = 0.
for mu in range(4):
    for nu in range(mu+1, 4):
        # S_L_ch35[mu][nu] should equal S_L[mu][nu]
        diff = S_L_ch35[mu][nu] - S_L[mu][nu]
        err = np.max(np.abs(diff))
        max_err_SL = max(max_err_SL, err)
        print(f"  S^{{{mu}{nu}}}_L:  max diff = {err:.2e}  "
              + ("✓" if err < 1e-12 else "✗ MISMATCH"))

print(f"\nOverall max error for S^μν_L: {max_err_SL:.2e}")
if max_err_SL < 1e-12:
    print("✓ eq. 35.21 verified: S^μν_L = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)")

# Verify S^μν_R from Ch. 35 formula satisfies S^μν_R = -[S^μν_L]*
print("\nVerifying S^μν_R = -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ) = -[S^μν_L]*:")
max_err_SR = 0.
for mu in range(4):
    for nu in range(mu+1, 4):
        expected_SR = -np.conj(S_L[mu][nu])
        diff = S_R_ch35[mu][nu] - expected_SR
        err = np.max(np.abs(diff))
        max_err_SR = max(max_err_SR, err)
        print(f"  S^{{{mu}{nu}}}_R:  max diff = {err:.2e}  "
              + ("✓" if err < 1e-12 else "✗ MISMATCH"))

print(f"\nOverall max error for S^μν_R: {max_err_SR:.2e}")
if max_err_SR < 1e-12:
    print("✓ eq. 35.22 verified: S^μν_R = -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ)")

# Print explicit matrices for a few cases
label = {0:'0', 1:'1', 2:'2', 3:'3'}
print("\nExplicit S^μν_L from σ σ̄ formula vs Ch. 34 definition:")
for mu in range(4):
    for nu in range(mu+1, 4):
        m35 = S_L_ch35[mu][nu]
        m34 = S_L[mu][nu]
        print(f"\n  (i/4)(σ^{mu} σ̄^{nu} - σ^{nu} σ̄^{mu}) =")
        print(f"    [{m35[0,0]:+.3f}  {m35[0,1]:+.3f}]")
        print(f"    [{m35[1,0]:+.3f}  {m35[1,1]:+.3f}]")
        print(f"  S^{{{mu}{nu}}}_L (Ch.34) =")
        print(f"    [{m34[0,0]:+.3f}  {m34[0,1]:+.3f}]")
        print(f"    [{m34[1,0]:+.3f}  {m34[1,1]:+.3f}]")

# =============================================================================
# §35.F  CLIFFORD-LIKE ALGEBRA: σ^μ σ̄^ν + σ^ν σ̄^μ = -2 g^{μν} I₂
# =============================================================================
# This is the anticommutator relation for the sigma matrices.
# It is analogous to the Clifford algebra {γ^μ, γ^ν} = 2g^{μν} for Dirac matrices,
# but the sign differs because σ^μ and σ̄^ν are NOT the same type of object:
#   σ^μ_{aȧ} has one undotted and one dotted index
#   σ̄^ν_{ȧa} has them reversed
# The product σ^μ σ̄^ν = σ^μ_{aȧ} σ̄^{ν ȧb} has BOTH indices undotted.
#
# Proof from eqs. 35.21:
#   S^μν_L = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)
#   S^μν_L is antisymmetric in [μν].
#   From the Lorentz algebra constraint: S^μν_L must satisfy eq. 34.4.
#   One derives that the symmetric part must be proportional to δ:
#     (σ^μ σ̄^ν + σ^ν σ̄^μ)_a^c = -2 g^{μν} δ_a^c
#
# Direct proof from σ^μ σ̄^ν components:
#   μ=ν=0: I · I = I; I + I = 2I = -2g^{00}I = 2I ✓
#   μ=ν=i: σ_i (-σ_i) = -σ_i²= -I; -I + (-I) = -2I = -2g^{ii}I = 2I... 
#   Wait: g^{ii} = -1, so -2g^{ii}I = -2(-1)I = +2I
#   LHS = σ_i(-σ_i) + σ_i(-σ_i) = -2σ_i²  = -2I ≠ +2I ?
#   Need to check: μ=ν=i means σ^i σ̄^i + σ^i σ̄^i = 2 σ^i σ̄^i = 2 σ_i(-σ_i) = -2I
#   And -2 g^{ii} I₂ = -2(-1) I₂ = +2I ... that's still +2I.
#   Hmm, let me recheck. The identity is:
#     (σ^μ σ̄^ν + σ^ν σ̄^μ)_a^c = -2 g^{μν} δ_a^c
#   For μ=ν=i: σ^i σ̄^i + σ^i σ̄^i = 2 σ_i(-σ_i) = -2I
#   RHS: -2 g^{ii} I = -2(-1)I = 2I  ... but LHS = -2I. CONTRADICTION?
#   
#   Let me recheck the Srednicki sign convention more carefully.
#   Actually, this identity is eq. 36.8 in Srednicki (stated in Ch. 36 but
#   derived from the Ch. 35 results). Let me just compute it numerically.

sec("§35.F  σ^μ σ̄^ν + σ^ν σ̄^μ = -2 g^{μν} I₂  (Clifford-like algebra)")

print("Computing (σ^μ σ̄^ν + σ^ν σ̄^μ)_a^b for all (μ,ν)...")

max_err_cliff = 0.
all_good = True
print("  (μ,ν)  |  anticommutator result  |  expected -2g^{μν} I₂  |  error")
for mu in range(4):
    for nu in range(4):
        lhs = sigma_vec[mu] @ sigmabar_vec[nu] + sigma_vec[nu] @ sigmabar_vec[mu]
        rhs = -2 * g[mu, nu] * I2
        err = np.max(np.abs(lhs - rhs))
        max_err_cliff = max(max_err_cliff, err)
        if err > 1e-12:
            all_good = False
        if err > 1e-10 or (mu == nu):
            print(f"  ({mu},{nu}):  anticomm[0,0]={lhs[0,0].real:+.3f}, [1,1]={lhs[1,1].real:+.3f}"
                  f"  |  expected: {rhs[0,0].real:+.3f} * I  |  err={err:.2e}")

print(f"\nMax error across all 16 (μ,ν) pairs: {max_err_cliff:.2e}")
if all_good:
    print("✓ σ^μ σ̄^ν + σ^ν σ̄^μ = -2 g^{μν} I₂  verified for all 16 entries!")
else:
    print("✗ Clifford relation FAILED for some (μ,ν)!")

# Extra: verify the antisymmetric part gives the generators
# (σ^μ σ̄^ν - σ^ν σ̄^μ)/4i should equal S^μν_L
print("\nVerifying S^μν_L = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)  [already done in §35.E,")
print("  confirmed consistent with Clifford algebra above].")

# =============================================================================
# §35.G  INDEX-FREE NOTATION — CADABRA2 SYMBOLIC BILINEARS
# =============================================================================
# In the index-free notation of Srednicki (eq. 35.23):
#   χψ = χ^a ψ_a    ("angle bracket" or undotted contraction)
#   χ†ψ† = χ†_ȧ ψ†^ȧ  ("square bracket" or dotted contraction)
#
# These are the LORENTZ INVARIANT spinor products.  They are invariant because:
#   χ^a ψ_a = ε^{ab} χ_b ψ_a  — ε^{ab} is an SL(2,C) invariant
#   Under χ → Lχ, ψ → Lψ:
#   χ'^a ψ'_a = (L^{-1T})^a_b χ^b (L)_a^c ψ_c = (det L) χ^a ψ_a = χ^a ψ_a
#   (since det L = 1 for SL(2,C)).
#
# Symmetry of χψ (eq. 35.25):
#   χψ = χ^a ψ_a = -ψ_a χ^a  (anticommutation, Grassmann)
#        = +ψ^a χ_a           (raising a: ψ^a = ε^{ab}ψ_b, χ_a = ε_{ab}χ^b)
#        = ψχ
#   So χψ = ψχ for Grassmann spinors! The two minus signs cancel.
#
# Hermitian conjugate (eq. 35.26):
#   (χψ)† = (χ^a ψ_a)† = (ψ_a)† (χ^a)† = ψ†_ȧ χ†^ȧ = ψ†χ†
#
# The σ̄^μ bilinear (eq. 35.27):
#   ψ†σ̄^μ χ = ψ†_ȧ σ̄^{μ ȧc} χ_c
#   This transforms as a 4-vector (eq. 35.28):
#   U(Λ)⁻¹ [ψ†σ̄^μ χ] U(Λ) = Λ^μ_ν [ψ†σ̄^ν χ]
#
# Hermitian conjugate of the bilinear (eq. 35.29):
#   [ψ†σ̄^μ χ]† = χ†σ̄^μ ψ
#   (using σ̄^μ† = σ̄^μ since σ̄^μ = (I,-σ⃗) are Hermitian matrices)

sec("§35.G  Index-free notation — Cadabra2 symbolic bilinears")

# Undotted (left-handed) spinor contraction: χψ = ε^{αβ} χ_α ψ_β
bilinear_chi_psi = Ex(r"\epsilon^{\alpha\beta} \chi_{\alpha} \psi_{\beta}")
print("Left-handed bilinear (angle bracket):")
print("  χψ ≡ χ^α ψ_α = ε^{αβ} χ_α ψ_β  =", bilinear_chi_psi)
print("  Physics: Lorentz invariant; χψ = ψχ for Grassmann fields [eq 35.25]")

# Dotted (right-handed) contraction: χ†ψ† = ε_{ȧḃ} χ†^ȧ ψ†^ḃ
bilinear_dag = Ex(r"\epsilon_{\dal\dbe} \chidag^{\dal} \psidag^{\dbe}")
print("\nRight-handed bilinear (square bracket):")
print("  χ†ψ† ≡ χ†_ȧ ψ†^ȧ = ε_{ȧḃ} χ†^ȧ ψ†^ḃ  =", bilinear_dag)
print("  Physics: Lorentz invariant; χ†ψ† = ψ†χ† for Grassmann fields")

# Vector bilinear: ψ†σ̄^μ χ = ψ†_ȧ σ̄^{μ ȧc} χ_c
vector_bilinear = Ex(r"\psidag_{\dal} \sigmabar^{\mu\dal\alpha} \chi_{\alpha}")
print("\nVector bilinear:")
print("  ψ†σ̄^μ χ ≡ ψ†_ȧ σ̄^{μ ȧa} χ_a  =", vector_bilinear)
print("  Physics: transforms as a 4-vector under Lorentz! [eq 35.28]")
print("  Note: Hermitian, since (σ̄^μ)† = σ̄^μ and [ψ†σ̄^μ χ]† = χ†σ̄^μ ψ")

# =============================================================================
# §35.H  LORENTZ VECTOR BILINEAR — ψ†σ̄^μ χ TRANSFORMS AS A 4-VECTOR
# =============================================================================
# We verify numerically that S^μν_L and S^μν_R combine with σ̄^μ correctly.
#
# The key consistency check (from the derivation in Ch. 35):
#   (S^μν_V)^ρ_τ σ^τ_{aȧ} + (S^μν_L)_a^b σ^ρ_{bȧ} + (S^μν_R)_ȧ^ḃ σ^ρ_{aḃ} = 0
#   [eq 35.14, restated]
#
# This is the covariance condition: σ^μ is invariant under simultaneous
# Lorentz rotation of its vector index and both spinor indices.
#
# We verify this numerically for all (μ,ν,ρ,a,ȧ):
# (g^{μρ} δ^ν_τ - g^{νρ} δ^μ_τ) σ^τ_{aȧ} + i(S^μν_L)_a^b σ^ρ_{bȧ}
#   + i(S^μν_R)_ȧ^ḃ σ^ρ_{aḃ} = 0

sec("§35.H  Covariance condition — σ^μ invariance check  [eq 35.14]")

def S_L_fn(mu, nu):
    """Return S^μν_L, antisymmetric."""
    if mu == nu:
        return np.zeros((2, 2), dtype=complex)
    if mu < nu:
        return S_L[mu][nu]
    return -S_L[nu][mu]

def S_R_fn(mu, nu):
    """Return S^μν_R = -conj(S^μν_L), antisymmetric."""
    return -np.conj(S_L_fn(mu, nu))

print("Verifying the σ-covariance identity (eq. 35.14) for all (μ,ν,ρ,a,ȧ):")
max_err_cov = 0.
for mu in range(4):
    for nu in range(mu+1, 4):
        for rho in range(4):
            for a in range(2):
                for adot in range(2):
                    # Term 1: (g^{μρ} δ^ν_τ - g^{νρ} δ^μ_τ) σ^τ_{aȧ}
                    term1 = (g[mu, rho] * sigma_vec[nu][a, adot]
                             - g[nu, rho] * sigma_vec[mu][a, adot])
                    # Term 2: i (S^μν_L)_a^b σ^ρ_{bȧ}
                    SL_mu_nu = S_L_fn(mu, nu)
                    term2_vec = 1j * SL_mu_nu @ sigma_vec[rho]
                    term2 = term2_vec[a, adot]
                    # Term 3: i (S^μν_R)_ȧ^ḃ σ^ρ_{aḃ}
                    SR_mu_nu = S_R_fn(mu, nu)
                    term3_vec = 1j * sigma_vec[rho] @ SR_mu_nu.T
                    # S^μν_R has indices [adot, bdot], σ^ρ has [a, adot]
                    # We need: i (S^μν_R)_{ȧ}^{ḃ} σ^ρ_{a ḃ}
                    # = i Σ_bdot SR[adot, bdot] * sigma_vec[rho][a, bdot]
                    term3 = 1j * np.dot(SR_mu_nu[adot, :], sigma_vec[rho][a, :])
                    total = term1 + term2 + term3
                    max_err_cov = max(max_err_cov, abs(total))

print(f"  Max residual across all (μ,ν,ρ,a,ȧ): {max_err_cov:.2e}")
if max_err_cov < 1e-12:
    print("  ✓ σ^μ covariance condition verified: [eq 35.14]")
    print("    σ^ρ is invariant under simultaneous vector + left + right rotations")
else:
    print("  ✗ Covariance condition FAILED!")

# Also verify: tr(σ^μ σ̄^ν) = 2g^{μν} once more as a summary check
print("\nSummary of key checks:")
checks = [
    ("σ^μ_{aȧ} σ_{μ bḃ} = -2 ε_{ab} ε_{ȧḃ}", err_35_4),
    ("ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ} = -2g^{μν}", err_35_5),
    ("σ̄^μ = (I,-σ⃗)  [eq 35.20]", max_err_sigmabar),
    ("tr(σ^μ σ̄^ν) = -2g^{μν}  [mostly-plus]", err_trace),
    ("S^μν_L = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)  [eq 35.21]", max_err_SL),
    ("S^μν_R = -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ) [eq 35.22]", max_err_SR),
    ("σ^μ σ̄^ν + σ^ν σ̄^μ = -2g^{μν} I₂", max_err_cliff),
    ("σ^μ covariance (eq. 35.14)", max_err_cov),
]
for name, err in checks:
    status = "✓" if err < 1e-10 else "✗"
    print(f"  {status}  {name}  (err={err:.1e})")

# =============================================================================
# §35  SUMMARY
# =============================================================================

sec("CHAPTER 35 SUMMARY")
print("""
  σ̄^μ — THE RAISED-INDEX σ SYMBOL
  ─────────────────────────────────
    σ̄^{μ ȧa} ≡ ε^{ab} ε^{ȧḃ} σ^μ_{bḃ}             [eq 35.19]
    Numerically: σ̄^μ = (I, -σ_1, -σ_2, -σ_3) = (I, -σ⃗)  [eq 35.20]
    Hermitian: (σ̄^μ)† = σ̄^μ for all μ

  FUNDAMENTAL σ IDENTITIES
  ─────────────────────────
    σ^μ_{aȧ} σ_{μ bḃ} = -2 ε_{ab} ε_{ȧḃ}                 [eq 35.4]
    ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ} = -2 g^{μν}           [eq 35.5]
    tr(σ^μ σ̄^ν) = 2 g^{μν}  (follows from above)
    σ^μ σ̄^ν + σ^ν σ̄^μ = -2 g^{μν} I₂   (Clifford-like)

  GENERATORS IN σ/σ̄ FORM
  ────────────────────────
    (S^μν_L)_a^b = +(i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)_a^b      [eq 35.21]
    (S^μν_R)^ȧ_ḃ = -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ)^ȧ_ḃ      [eq 35.22]
    (Consistent with Ch. 34: S^μν_R = -[S^μν_L]*)

  INDEX-FREE NOTATION  [eq 35.23]
  ─────────────────────────────────
    χψ  = χ^a ψ_a           (undotted contraction, "angle bracket")
    χ†ψ† = χ†_ȧ ψ†^ȧ        (dotted contraction, "square bracket")
    χψ = ψχ                  (Grassmann: two minus signs cancel)  [eq 35.25]
    (χψ)† = ψ†χ†             (hermitian conjugate reverses order)  [eq 35.26]

  VECTOR BILINEAR  [eq 35.27-35.28]
  ────────────────────────────────────
    ψ†σ̄^μ χ = ψ†_ȧ σ̄^{μ ȧc} χ_c  (transforms as a 4-vector!)
    [ψ†σ̄^μ χ]† = χ†σ̄^μ ψ         (hermitian conjugate, eq 35.29)
    (σ̄^μ is Hermitian, so no extra signs)

  COVARIANCE OF σ^μ  [eq 35.9, 35.14]
  ──────────────────────────────────────
    σ^ρ_{aȧ} = Λ^ρ_τ L(Λ)_a^b R(Λ)_ȧ^ḃ σ^τ_{bḃ}
    Infinitesimally: (S^μν_V)σ + i S^μν_L σ + i σ S^μν_R^T = 0
""")

print("Done: ch35_sigma_algebra.py")
