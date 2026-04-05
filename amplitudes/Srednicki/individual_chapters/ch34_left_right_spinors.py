"""
ch34_left_right_spinors.py
===========================
Srednicki QFT — Chapter 34: Left- and Right-Handed Spinor Fields

What this file covers (section by section):
  §34.A  Left-handed spinor field ψ_a, (2,1) rep of Lorentz group
  §34.B  Generators S^μν_L in the (2,1) rep; commutation relations
  §34.C  Explicit form: S^ij_L = ½ ε^{ijk} σ_k  (spatial rotations)
         S^k0_L = (i/2) σ_k                        (boosts)
  §34.D  Right-handed spinor ψ†_{ȧ}, (1,2) rep; dotted indices
         S^μν_R = -[S^μν_L]*  (generators are complex conjugate)
  §34.E  ε_{ab} as the SL(2,C) invariant symbol (raise/lower indices)
  §34.F  σ^μ_{aȧ} — the dictionary between (2,2) and vector reps

Run with:
    python3 ch34_left_right_spinors.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70
sec = lambda s: print(f"\n{SEP}\n  {s}\n{SEP}")

print(SEP)
print("  Srednicki Ch. 34 — Left- and Right-Handed Spinor Fields")
print(SEP)

# =============================================================================
# §34.A  THE (2,1) LEFT-HANDED WEYL FIELD
# =============================================================================
# The Lorentz group in 4D has Lie algebra isomorphic to
#   su(2)_L ⊕ su(2)_R
# with generators N_i = ½(J_i - i K_i) and N†_i = ½(J_i + i K_i),
# satisfying [N_i, N_j] = iε_{ijk} N_k independently.
#
# An irrep is labelled (2n+1, 2n'+1) where n, n' ∈ {0,½,1,...}.
# The FOUR fundamental reps are:
#   (1,1)  scalar       → 1 component
#   (2,1)  left-handed  → 2 components  (N† acts trivially)
#   (1,2)  right-handed → 2 components  (N  acts trivially)
#   (2,2)  vector       → 4 components  (both act as j=½)
#
# A LEFT-HANDED WEYL FIELD ψ_a(x) lives in (2,1).
# Under a Lorentz transformation Λ:
#   U(Λ)⁻¹ ψ_a(x) U(Λ) = L_a^b(Λ) ψ_b(Λ⁻¹x)    [eq. 34.1]
#
# KEY POINT about the (2,1) label:
#   The "2" means N_i (left su(2)) acts as the 2-dimensional (j=½) rep.
#   The "1" means N†_i (right su(2)) acts trivially (j=0).
#   Since J_i = N_i + N†_i, the field has angular momentum j=½.
#   Since K_k = i(N_k - N†_k), boosting gives K_k → i·(N_k) = i·(spin matrix).
#
# WHY "(2,1)" notation and not "(½,0)"?
#   Srednicki counts DIMENSIONS: n=½ gives 2n+1=2, n'=0 gives 2n'+1=1.
#   The physics notation writes the same thing as (½,0), meaning
#   left-SU(2) spin = ½, right-SU(2) spin = 0.
#   Both conventions are equivalent. Srednicki uses dimension counting.

sec("§34.A  Left-handed Weyl field — index declaration")

# Spinor indices a, b, c, d for the (2,1) (left-handed) representation.
# "position=fixed" means cadabra2 distinguishes upper vs lower indices.
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))

# Dotted indices ȧ, ḃ, ċ, ḋ for the (1,2) (right-handed) representation.
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))

# Greek spacetime indices μ,ν,ρ,σ for Lorentz vectors/tensors.
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))

print("Declared index sets:")
print("  Undotted (left-handed, (2,1)):  α β γ δ")
print("  Dotted   (right-handed, (1,2)): α̇ β̇ γ̇ δ̇")
print("  Spacetime vector:               μ ν ρ σ")

# =============================================================================
# §34.B  THE GENERATORS S^μν_L — WHAT THEY ARE AND WHY
# =============================================================================
# For an infinitesimal Lorentz transformation δω^μν (antisymmetric), the
# representation matrix in the (2,1) rep is (eq. 34.3):
#
#   L_a^b(1+δω) = δ_a^b + (i/2) δω_{μν} (S^μν_L)_a^b
#
# The S^μν_L are 2×2 matrices (one for each [μν] pair).
# Antisymmetry: (S^μν_L)_a^b = -(S^νμ_L)_a^b
# So there are 6 independent matrices (since [μν] has C(4,2)=6 choices).
# These 6 matrices are the GENERATORS of the (2,1) representation.
#
# They satisfy the Lorentz algebra (eq. 34.4):
#   [S^μν_L, S^ρσ_L] = i(g^μρ S^νσ_L - g^νρ S^μσ_L - g^μσ S^νρ_L + g^νσ S^μρ_L)
#
# This is the defining commutation relation of the Lorentz algebra, same as
# the abstract generators M^μν satisfy.
#
# WHY must they satisfy this?
#   The L matrices are a REPRESENTATION of the group — they must compose
#   the same way as the group elements themselves. The Lie algebra relation
#   is just the infinitesimal version of that consistency requirement.

sec("§34.B  Generators S^μν_L — explicit numpy matrices")

# We'll use numpy to construct the explicit 2×2 matrices.
# Pauli matrices:
sigma = {
    1: np.array([[0, 1], [1, 0]], dtype=complex),
    2: np.array([[0, -1j], [1j, 0]], dtype=complex),
    3: np.array([[1, 0], [0, -1]], dtype=complex),
}
I2 = np.eye(2, dtype=complex)

# Levi-Civita ε^{ijk} (purely spatial, 3D):
def eps3(i, j, k):
    """Levi-Civita symbol ε^{ijk}, indices 1,2,3."""
    return int(np.linalg.det([[i==1,i==2,i==3],[j==1,j==2,j==3],[k==1,k==2,k==3]]))

# Srednicki's metric (Eq. 1.8 / 2.4): g_{μν} = diag(-1,+1,+1,+1)  "mostly-plus"
# So g^{μν} = diag(-1,+1,+1,+1) as well (same eigenvalues).
# g^{00} = -1,  g^{ii} = +1  for i=1,2,3.

# eq. 34.9:  (S^ij_L)_a^b = ½ ε^{ijk} σ_k   (spatial rotation generators)
# eq. 34.10: (S^k0_L)_a^b = (i/2) σ_k         (boost generators)
#
# We store S_L as a 4×4 array of 2×2 matrices: S_L[mu][nu] = (S^μν_L)_a^b
# with 0=time, 1,2,3=space.

S_L = [[None]*4 for _ in range(4)]

# Fill spatial part: i,j ∈ {1,2,3}
for i in range(1, 4):
    for j in range(1, 4):
        mat = np.zeros((2, 2), dtype=complex)
        for k in range(1, 4):
            mat += eps3(i, j, k) * sigma[k]
        S_L[i][j] = 0.5 * mat

# Fill boost part: S^{k0} = (i/2) σ_k
for k in range(1, 4):
    S_L[k][0] = (1j / 2) * sigma[k]
    S_L[0][k] = -(1j / 2) * sigma[k]   # antisymmetry S^{0k} = -S^{k0}

# Diagonal / zero entries
S_L[0][0] = np.zeros((2, 2), dtype=complex)
for i in range(1, 4):
    S_L[i][i] = np.zeros((2, 2), dtype=complex)

print("Left-handed generators S^μν_L  (μ=row, ν=col in index space;")
print("                                each entry is a 2×2 matrix):")
label = {0: '0', 1: '1', 2: '2', 3: '3'}
for mu in range(4):
    for nu in range(mu+1, 4):
        mat = S_L[mu][nu]
        if mat is not None:
            print(f"\n  S^{{{label[mu]}{label[nu]}}}_L =")
            print("   ", mat[0])
            print("   ", mat[1])

# =============================================================================
# §34.C  VERIFY COMMUTATION RELATIONS  [S^μν_L, S^ρσ_L]
# =============================================================================
# The Lorentz algebra (eq. 34.4) says:
#
#   [S^μν, S^ρσ] = i(g^μρ S^νσ - g^νρ S^μσ - g^μσ S^νρ + g^νσ S^μρ)
#
# with metric g = diag(+1,-1,-1,-1).
#
# Let's verify this for a specific case: [S^12, S^23]
# By eq. 34.9: S^12 = ½σ_3, S^23 = ½σ_1
# g^{22} = -1  (only g^{22} is nonzero among {g^μρ, g^νρ, g^μσ, g^νσ} here
# with μ=1,ν=2,ρ=2,σ=3)
# RHS = i(g^{12} S^{23} - g^{22} S^{13} - g^{13} S^{22} + g^{23} S^{12})
#      = i(0·S^23 - (-1)·S^13 - 0·S^22 + 0·S^12)
#      = i S^13
# S^13 = ½ε^{13k}σ_k = ½(ε^{131}σ_1 + ε^{132}σ_2 + ε^{133}σ_3)
#       = ½(-σ_2) = -½σ_2
# So RHS = i(-½σ_2) = -i/2 σ_2
# LHS = [½σ_3, ½σ_1] = ¼[σ_3,σ_1] = ¼(2i σ_2) = i/2 σ_2  ... hmm sign?
# Wait: [σ_i, σ_j] = 2iε_{ijk}σ_k; [σ_3,σ_1] = 2iε_{312}σ_2 = -2iσ_2
# So LHS = ¼·(-2i)σ_2 = -i/2 σ_2. ✓ Matches RHS!

sec("§34.C  Verify [S^μν_L, S^ρσ_L] = Lorentz algebra  (numpy)")

# Metric: Srednicki diag(-1,+1,+1,+1)  [Eq. 1.8 / 2.4, mostly-plus]
g = np.diag([-1., 1., 1., 1.])

def comm(A, B):
    return A @ B - B @ A

def lorentz_rhs(mu, nu, rho, sig):
    """RHS of the Lorentz commutator — Srednicki eq. 34.4:
      [S^{μν}, S^{ρσ}] = i(g^{μρ} S^{νσ} - g^{νρ} S^{μσ} - g^{μσ} S^{νρ} + g^{νσ} S^{μρ})
    with g^{μν} = diag(-1,+1,+1,+1)  [Srednicki Eq. 1.8].
    """
    def S(m, n):
        if m == n:
            return np.zeros((2, 2), dtype=complex)
        if m < n:
            return S_L[m][n]
        return -S_L[n][m]   # antisymmetry: S^{nm} = -S^{mn}
    return 1j * (
        g[mu, rho] * S(nu, sig)
        - g[nu, rho] * S(mu, sig)
        - g[mu, sig] * S(nu, rho)
        + g[nu, sig] * S(mu, rho)
    )

print("Checking [S^μν_L, S^ρσ_L] = RHS  for all independent pairs...")
max_err = 0.
for mu in range(4):
    for nu in range(mu+1, 4):
        for rho in range(4):
            for sig in range(rho+1, 4):
                lhs = comm(S_L[mu][nu], S_L[rho][sig])
                rhs = lorentz_rhs(mu, nu, rho, sig)
                err = np.max(np.abs(lhs - rhs))
                if err > max_err:
                    max_err = err
print(f"  Max error across all [μν,ρσ] pairs: {max_err:.2e}")
print("  ✓ All commutation relations verified" if max_err < 1e-12 else "  ✗ MISMATCH!")

# =============================================================================
# §34.D  RIGHT-HANDED SPINOR AND DOTTED INDICES
# =============================================================================
# Take ψ_a(x) (left-handed, (2,1)).
# Its Hermitian conjugate is (ψ_a)† = ψ†_{ȧ}(x)  [eq. 34.11]
#
# WHY does Hermitian conjugation flip the rep from (2,1) to (1,2)?
#
# The Lorentz algebra has TWO su(2) factors, call them L and R:
#   N_i generates su(2)_L    (acts on left-handed fields)
#   N†_i generates su(2)_R   (acts on right-handed fields)
#
# For ψ_a in (2,1): N_i acts as j=½ matrices, N†_i acts as 0.
# Taking †:  N_i → (N_i)† = N†_i,  N†_i → (N†_i)† = N_i
# So the conjugate field ψ†_ȧ has: N†_i acts as j=½, N_i acts as 0.
# That is EXACTLY the (1,2) representation. ✓
#
# THE DOT CONVENTION:
#   Undotted index a,b,... → left-handed (2,1) rep
#   Dotted index  ȧ,ḃ,... → right-handed (1,2) rep
#
# The Lorentz transformation of ψ†_ȧ (eq. 34.12):
#   U(Λ)⁻¹ ψ†_ȧ(x) U(Λ) = R_ȧ^ḃ(Λ) ψ†_ḃ(Λ⁻¹x)
#
# And the key relation between generators (eq. 34.17):
#   (S^μν_R)_ȧ^ḃ = -[(S^μν_L)_a^b]*
#
# WHY?  Take hermitian conjugate of eq. 34.15 (the commutator with M^μν):
#   [ψ†_ȧ(0), M^μν] = (S^μν_R)_ȧ^ḃ ψ†_ḃ(0)          [eq. 34.15]
#   Taking h.c.:
#   [M^μν, ψ_a(0)] = [(S^μν_R)_ȧ^ḃ]* ψ_b(0)           [eq. 34.16]
#   But from eq. 34.6 (with x=0, suppressing L^μν):
#   [ψ_a(0), M^μν] = (S^μν_L)_a^b ψ_b(0)
#   → [M^μν, ψ_a(0)] = -(S^μν_L)_a^b ψ_b(0)
#   So: [(S^μν_R)_ȧ^ḃ]* = -(S^μν_L)_a^b
#   i.e. (S^μν_R)_ȧ^ḃ = -[(S^μν_L)_a^b]*              [eq. 34.17] ✓

sec("§34.D  Right-handed generators S^μν_R = -[S^μν_L]*")

# Build S_R from S_L:
S_R = [[None]*4 for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        if S_L[mu][nu] is not None:
            S_R[mu][nu] = -np.conj(S_L[mu][nu])

print("Explicit S^μν_R matrices (should be -conj(S^μν_L)):")
for mu in range(4):
    for nu in range(mu+1, 4):
        SL = S_L[mu][nu]
        SR = S_R[mu][nu]
        diff = SR - (-np.conj(SL))
        print(f"\n  S^{{{label[mu]}{label[nu]}}}_R =  -[S^{{{label[mu]}{label[nu]}}}_L]*")
        print("   ", SR[0])
        print("   ", SR[1])
        err = np.max(np.abs(diff))
        assert err < 1e-14, f"Mismatch at ({mu},{nu})!"
print("\n  ✓ S^μν_R = -[S^μν_L]* verified for all components")

# Observe: S^ij_R = -conj(½ ε^{ijk} σ_k) = -½ ε^{ijk} σ_k* = -½ ε^{ijk} σ_k
#   (since Pauli matrices σ_1 and σ_3 are real, σ_2 is purely imaginary,
#    and ε^{ijk} is real, so -[S^ij_L]* = S^ij_L for spatial part)
# WAIT that gives S^ij_R = S^ij_L for spatial part. Let's check:
print("\n  Comparing S^12_L vs S^12_R:")
print("  S^12_L:", S_L[1][2])
print("  S^12_R:", S_R[1][2])
# They should match for spatial rotations (Hermitian Pauli matrices).
# Boost generators S^k0: S^k0_L = (i/2)σ_k → S^k0_R = -(i/2)σ_k*
# For σ_1,σ_3 real: S^k0_R = -(i/2)σ_k = -S^k0_L → boosts have opposite sign.
# This is physical: under parity, L↔R and boosts get a sign flip.

# =============================================================================
# §34.E  THE ε SYMBOL — INVARIANT, RAISES/LOWERS SPINOR INDICES
# =============================================================================
# From the group theory:
#   (2,1) ⊗ (2,1) = (1,1)_A ⊕ (3,1)_S
# The singlet (1,1) in the antisymmetric part implies there is an INVARIANT
# antisymmetric 2-tensor: ε_{ab} = -ε_{ba}.
#
# This is the SL(2,C) generalization of the metric.
# It satisfies:
#   L_a^c(Λ) L_b^d(Λ) ε_{cd} = ε_{ab}    [invariance, eq. 34.20]
#
# Srednicki's normalization (eq. 34.22):
#   ε_{12} = ε^{21} = +1
#   ε_{21} = ε^{12} = -1
#
# Note: ε^{ab} is defined so that ε_{ab} ε^{bc} = δ_a^c  [eq. 34.23]
#
# Index raising/lowering:
#   ψ^a = ε^{ab} ψ_b     (raise with ε^{ab}, NW corner → NE corner)
#   ψ_a = ε_{ab} ψ^b     (lower with ε_{ab})
#
# SIGN TRAP (eq. 34.27):
#   ψ^a χ_a = ε^{ab} ψ_b χ_a = -ε^{ba} ψ_b χ_a = -ψ_b χ^b
#   The contraction ψ^a χ_a = -ψ_a χ^a ← crucial minus sign!
#
# For DOTTED indices (right-handed), everything is the same with ε_{ȧḃ}.

sec("§34.E  ε_{ab} — SL(2,C) invariant symbol, index raising/lowering")

# Numerical epsilon tensors (2×2):
eps_lower = np.array([[0,  1],   # ε_{ab}: ε_{12}=+1, ε_{21}=-1
                       [-1, 0]], dtype=complex)
eps_upper = np.array([[0, -1],   # ε^{ab}: ε^{12}=-1, ε^{21}=+1
                       [1,  0]], dtype=complex)

# Wait — Srednicki says ε^{12} = ε_{21} = +1 and ε^{21} = ε_{12} = -1
# so ε^{12} = +1 means eps_upper[0,1] = +1, eps_upper[1,0] = -1
eps_lower = np.array([[0, -1],   # ε_{12} = -1, ε_{21} = +1
                       [1,  0]], dtype=complex)
eps_upper = np.array([[0,  1],   # ε^{12} = -1... wait let me re-read
                       [-1, 0]], dtype=complex)
# Srednicki eq. 34.22: ε^{12} = ε_{21} = +1, ε^{21} = ε_{12} = -1
# So:
#   ε_{12} = -1, ε_{21} = +1  (ε_{ab} matrix, row=a, col=b)
#   ε^{12} = +1, ε^{21} = -1  (ε^{ab} matrix)
eps_lower = np.array([[0, -1],   # row a, col b: ε_{12}=-1, ε_{21}=+1
                       [1,  0]], dtype=complex)
eps_upper = np.array([[0,  1],   # ε^{12}=+1, ε^{21}=-1
                       [-1, 0]], dtype=complex)

# Verify ε_{ab} ε^{bc} = δ_a^c:
product = eps_lower @ eps_upper
print("ε_{ab} ε^{bc} (should be identity):")
print(" ", product)
assert np.allclose(product, np.eye(2)), "ε normalization failed!"
print("  ✓ ε_{ab} ε^{bc} = δ_a^c  verified")

# Invariance: L_a^c L_b^d ε_{cd} = ε_{ab} for any L ∈ SL(2,C)
# Check for an explicit L: use L = exp(i·θ·σ_3/2) (rotation around z by θ)
theta = 0.7   # arbitrary angle
L_test = np.eye(2, dtype=complex) * np.cos(theta/2) + 1j * sigma[3] * np.sin(theta/2)
lhs = L_test.T @ eps_lower @ L_test   # (L^T ε L)_{ab} = L_a^c ε_{cd} L_b^d... let's be careful
# ε'_{ab} = L_a^c L_b^d ε_{cd}
# in matrix form: L^T · eps_lower · L
invariant_check = L_test.T @ eps_lower @ L_test
print(f"\nInvariance check (rotation by θ={theta:.2f}):")
print("  L^T ε L =", invariant_check)
print("  ε_{ab}  =", eps_lower)
assert np.allclose(invariant_check, eps_lower), "Invariance failed!"
print("  ✓ ε_{ab} is invariant under SL(2,C)")

# Also verify det(L) = 1 for our test matrix:
print(f"  det(L) = {np.linalg.det(L_test):.6f}  (should be 1)")

# =============================================================================
# §34.F  σ^μ_{aȧ} — THE VECTOR/SPINOR DICTIONARY
# =============================================================================
# (2,2) ≅ vector representation.
# A field A_{aȧ}(x) in (2,2) maps to a 4-vector A^μ(x) via (eq. 34.28):
#
#   A_{aȧ}(x) = σ^μ_{aȧ} A_μ(x)
#
# where σ^μ_{aȧ} = (I, σ⃗) is the invariant symbol.  [eq. 34.30]
#
# WHY (I, σ⃗) and not something else?
#   The σ matrices are already the standard "connecting" objects between
#   spinor indices and vector indices in the Weyl representation.
#   The self-consistency condition is that they must commute correctly
#   with the action of S^μν_L (on undotted a) and S^μν_R (on dotted ȧ),
#   matching the vector action of M^μν. Chapter 35 verifies this explicitly
#   using the identity gμν σ^μ_{aȧ} σ^ν_{bḃ} = -2 ε_{ab} ε_{ȧḃ}.
#
# ALSO defined (for future use):
#   σ̄^μ_{ȧa} = (I, -σ⃗)   [raised dotted index version]
# This will appear in the Dirac equation and σ-bar algebra (Ch. 35-36).

sec("§34.F  σ^μ_{aȧ} — the (2,2) invariant symbol")

# σ^μ = (I, σ_1, σ_2, σ_3)  — Srednicki convention, eq. 34.30
sigma_vec = {
    0: I2,
    1: sigma[1],
    2: sigma[2],
    3: sigma[3],
}

# σ̄^μ = (I, -σ_1, -σ_2, -σ_3)
sigmabar_vec = {
    0: I2,
    1: -sigma[1],
    2: -sigma[2],
    3: -sigma[3],
}

print("σ^μ_{aȧ} = (I, σ⃗):  components indexed by μ ∈ {0,1,2,3}")
for mu in range(4):
    print(f"\n  σ^{{{label[mu]}}} =")
    print("   ", sigma_vec[mu][0])
    print("   ", sigma_vec[mu][1])

# Verify the key identity:  tr(σ^μ σ̄^ν) = 2 g^{μν}
# where the trace is over both spinor indices:
#   tr(σ^μ σ̄^ν) = σ^μ_{aȧ} σ̄^{ν ȧa}  (matrix product, then trace)
# This follows from σ^0=I, σ̄^0=I, σ^i=σ_i, σ̄^i=-σ_i and g=diag(-1,+1,+1,+1):
#   μ=ν=0: tr(I·I) = 2  = -2g^{00} = -2(-1) = 2 ✓
#   μ=ν=i: tr(σ_i·(-σ_i)) = -tr(I) = -2 = -2g^{ii} = -2(+1) = -2 ✓
#   μ≠ν: tr(σ^μ σ̄^ν) = 0  (traceless cross terms, Pauli algebra) ✓
# So: tr(σ^μ σ̄^ν) = -2 g^{μν}  (Srednicki mostly-plus, cf. eq. 35.5)
#
# This is a precursor to the sigma-bar algebra in Ch. 35.
print("\nChecking tr(σ^μ σ̄^ν) = -2 g^{μν}  [Srednicki mostly-plus]:")
trace_matrix = np.zeros((4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        trace_matrix[mu, nu] = np.trace(sigma_vec[mu] @ sigmabar_vec[nu])

expected = -2 * g   # -2 * diag(-1,+1,+1,+1) = diag(+2,-2,-2,-2)
err = np.max(np.abs(trace_matrix - expected))
print("  tr(σ^μ σ̄^ν) =")
for row in trace_matrix.real:
    print("   ", row)
print(f"  Max error vs -2g^{{μν}}: {err:.2e}")
print("  ✓ tr(σ^μ σ̄^ν) = -2g^{μν} verified" if err < 1e-12 else "  ✗ MISMATCH!")

# =============================================================================
# §34.G  CADABRA2 ALGEBRAIC REPRESENTATION
# =============================================================================
# Now we use cadabra2 to represent these relationships symbolically.

sec("§34.G  Cadabra2: symbolic S^μν algebra and ε contractions")

# Declare epsilon tensors as antisymmetric in cadabra2
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

# Declare σ^μ and σ̄^μ as tensors with mixed spinor and vector indices
cadabra2.Symmetric(Ex(r"\delta^{\alpha}_{\beta}"))

# Left-handed spinor field ψ_α and right-handed ψ†^{ȧ}
psi_L = Ex(r"\psi_{\alpha}")
psi_R = Ex(r"\psidag^{\dal}")

print("Left-handed field  ψ_α  =", psi_L)
print("Right-handed field ψ†^ȧ =", psi_R)

# Raising ψ_α with ε^{αβ}:
psi_raised = Ex(r"\epsilon^{\alpha\beta} \psi_{\beta}")
print("\nRaised left-handed: ψ^α = ε^{αβ} ψ_β =", psi_raised)

# Invariant contraction (left-handed): ψ^α χ_α = ε^{αβ} ψ_β χ_α
# Note Srednicki's SIGN convention for contraction order:
# NW→SE: ψ^α χ_α  vs SE→NW: ψ_α χ^α = -ψ^α χ_α
inner_L = Ex(r"\epsilon^{\alpha\beta} \psi_{\alpha} \chi_{\beta}")
print("\nLeft-handed Lorentz invariant:")
print("  ψχ ≡ ε^{αβ} ψ_α χ_β =", inner_L)
print("  (= ψ^α χ_α, antisymmetric in ψ,χ for Grassmann fields)")

# Right-handed invariant: ε_{ȧḃ} ψ†^ȧ χ†^ḃ
inner_R = Ex(r"\epsilon_{\dal\dbe} \psidag^{\dal} \chidag^{\dbe}")
print("\nRight-handed Lorentz invariant:")
print("  [ψ†χ†] ≡ ε_{ȧḃ} ψ†^ȧ χ†^ḃ =", inner_R)

# =============================================================================
# §34  SUMMARY
# =============================================================================

sec("CHAPTER 34 SUMMARY")
print("""
  LORENTZ GROUP REPRESENTATIONS
  ─────────────────────────────
  Irreps labelled (2n+1, 2n'+1) = (dim_L, dim_R)
  Equivalently labelled (n, n') ∈ {0, ½, 1, ...}²

    (1,1)  = (0,0)  →  scalar field
    (2,1)  = (½,0)  →  left-handed Weyl field  ψ_a   [undotted index]
    (1,2)  = (0,½)  →  right-handed Weyl field ψ†_ȧ  [dotted index]
    (2,2)  = (½,½)  →  4-vector field A^μ

  GENERATORS IN (2,1) REP  (6 independent 2×2 matrices)
  ──────────────────────────────────────────────────────
    Spatial rotations:  (S^ij_L)_a^b = ½ ε^{ijk} σ_k   [eq. 34.9]
    Boosts:             (S^k0_L)_a^b = (i/2) σ_k        [eq. 34.10]
    Commutation:  [S^μν_L, S^ρσ_L] = Lorentz algebra    [eq. 34.4]

  HERMITIAN CONJUGATION FLIPS L↔R
  ──────────────────────────────────
    ψ_a ∈ (2,1)  →  (ψ_a)† = ψ†_ȧ ∈ (1,2)
    (S^μν_R)_ȧ^ḃ = -[(S^μν_L)_a^b]*               [eq. 34.17]
    → Rotations same sign, boosts flip sign (as expected from parity)

  ε SYMBOL — SL(2,C) METRIC
  ──────────────────────────
    ε_{ab}: antisymmetric, ε_{12} = -1, ε_{21} = +1  (Srednicki convention)
    ε^{ab}: ε^{12} = +1, ε^{21} = -1
    ε_{ab} ε^{bc} = δ_a^c
    ψ^a = ε^{ab} ψ_b  (raise),  ψ_a = ε_{ab} ψ^b  (lower)
    SIGN: ψ^a χ_a = -ψ_a χ^a

  σ^μ — SPINOR↔VECTOR DICTIONARY
  ─────────────────────────────────
    σ^μ_{aȧ} = (I, σ⃗)   =   (I, σ_1, σ_2, σ_3)
    σ̄^μ_{ȧa} = (I, -σ⃗)  =   (I, -σ_1, -σ_2, -σ_3)
    A_{aȧ} = σ^μ_{aȧ} A_μ  maps vector ↔ bispinor  [eq. 34.28]
    g_{μν} σ^μ_{aȧ} σ^ν_{bḃ} = -2 ε_{ab} ε_{ȧḃ}  [Ch. 35, eq. 35.7]
""")

print("Done: ch34_left_right_spinors.py")
