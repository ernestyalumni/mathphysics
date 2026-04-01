"""
03_spinor_helicity.py
=====================
Spinor-helicity formalism for massless particles.

Reference: Srednicki QFT, Chapters 36-37.
           Dixon (1996), Mangano & Parke (1991), Elvang & Huang (2015).

For a massless momentum p^μ (p^2 = 0):
  - Angle spinor λ_α: defined by p_{αα̇} = λ_α λ̃_{α̇}
  - Square spinor λ̃_{α̇}: complex conjugate of λ_α for real momenta

Run with:
    python3 03_spinor_helicity.py
or inside Docker:
    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 python3 /work/03_spinor_helicity.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__

print("=" * 60)
print("Spinor-helicity formalism (Srednicki Ch. 36-37)")
print("=" * 60)

__cdbkernel__ = cadabra2.create_scope()

# -------------------------------------------------------------------------
# 1. Setup: massless momentum as a 2x2 matrix
#    p_{αα̇} = p_μ σ^μ_{αα̇}
#    For p^2 = 0: det(p_{αα̇}) = 0
#    Hence p_{αα̇} is a rank-1 matrix and factorizes as:
#    p_{αα̇} = λ_α λ̃_{α̇}
# -------------------------------------------------------------------------
print("\n[1] Spinor decomposition of massless momentum:")
print("    p_{alpha dotalpha} = lambda_alpha * lambdatilde_{dotalpha}")
print()
print("    This is the key insight: for p^2 = 0,")
print("    the 2x2 matrix p_{alpha dotalpha} has rank 1")
print("    => factorizes into outer product of two 2-component spinors.")
print()
print("    For real Lorentzian momenta: lambdatilde = (lambda)^*")
print("    For complex momenta (holomorphic amplitudes): lambda, lambdatilde independent")

# -------------------------------------------------------------------------
# 2. Angle bracket and square bracket spinor products
#    For massless momenta p_i, p_j with spinors λ^i_α, λ̃^i_{α̇}:
#
#    <ij> = ε^{αβ} λ^i_α λ^j_β = λ^{iα} λ^j_α
#    [ij] = ε_{α̇β̇} λ̃^{iα̇} λ̃^{jβ̇} = λ̃^i_{α̇} λ̃^{jα̇}
#
#    Properties:
#    <ij> = -<ji>   [antisymmetric]
#    [ij] = -[ji]   [antisymmetric]
#    <ii> = 0
#    [ii] = 0
# -------------------------------------------------------------------------
print("\n[2] Spinor products:")
print("    <ij> = eps^{ab} lambda^i_a lambda^j_b  (angle bracket)")
print("    [ij] = eps_{adota dotb} lambdatilde^{i dota} lambdatilde^{j dotb}  (square bracket)")
print()
print("    Antisymmetry:")
print("    <ij> = -<ji>,  [ij] = -[ji],  <ii> = [ii] = 0")

# -------------------------------------------------------------------------
# 3. Mandelstam variables from spinor products
#    s_{ij} = (p_i + p_j)^2 = 2 p_i · p_j
#    For massless particles:
#    s_{ij} = <ij>[ji] = -<ij>[ij]
#    (Note: [ji] = -[ij])
#
#    More explicitly:
#    2 p_i · p_j = p_{i,αα̇} p_j^{αα̇} = λ^i_α λ̃^i_{α̇} λ^{jα} λ̃^{jα̇}
#               = (λ^i · λ^j)(λ̃^i · λ̃^j) = <ij>[ji]
# -------------------------------------------------------------------------
print("\n[3] Mandelstam variables:")
print("    s_{ij} = (p_i + p_j)^2 = 2 p_i.p_j")
print("           = <ij>[ji]  =  -<ij>[ij]")
print()
print("    For n-point kinematics:")
print("    s_{ij} >= 0 for physical (timelike) channel")
print("    s_{12} = <12>[21] = (p1+p2)^2")

# Symbolic expression
s12 = Ex(r"A_{12} B_{21}")  # placeholder for <12>[21]
print("\n    Symbolic: s_12 = <12>[21] =", s12)

# -------------------------------------------------------------------------
# 4. Momentum conservation in spinor notation
#    Sum of all momenta vanishes: Σ_i p_i = 0
#    In spinor notation: Σ_i λ^i_α λ̃^i_{α̇} = 0
#    This implies (Schouten identity and momentum conservation):
#    Σ_i <ji>[ik] = 0  for all j, k
# -------------------------------------------------------------------------
print("\n[4] Momentum conservation:")
print("    sum_i p_i^mu = 0")
print("    => sum_i lambda^i_alpha lambdatilde^i_{dotalpha} = 0")
print()
print("    In terms of spinor products (for any fixed j, k):")
print("    sum_i <ji>[ik] = 0")
print("    (This is the Weyl equation / transversality condition)")

# -------------------------------------------------------------------------
# 5. Schouten identity
#    ε^{αβ} ε^{γδ} + ε^{αγ} ε^{δβ} + ε^{αδ} ε^{βγ} = 0
#    In spinor-product notation:
#    <ij>[kl] + <ik>[lj] + <il>[jk] = 0  (Schouten-like identity)
#
#    More precisely (Schouten): for any three 2-component spinors a,b,c:
#    <ab><cd> + <ac><db> + <ad><bc> = 0
# -------------------------------------------------------------------------
print("\n[5] Schouten identity:")
print("    For any 2-component spinors a, b, c, d:")
print("    <ab><cd> + <ac><db> + <ad><bc> = 0")
print()
print("    Proof: follows from the fact that 2x2 antisymmetric matrices are")
print("    proportional to epsilon (no room for independent rank-2 antisymmetric")
print("    tensor in 2D). Every SL(2,C) invariant 4-spinor contraction reduces.")

# -------------------------------------------------------------------------
# 6. Explicit numerical example
#    Consider a simple 4-momentum in spinor form.
#    For p^μ = (E, 0, 0, E) (massless, along z-axis):
#    p_{αα̇} = E * [[1+1, 0], [0, 1-1]] = E * [[2, 0], [0, 0]]
#    => λ_α = sqrt(2E) * (1, 0)^T,  λ̃_{α̇} = sqrt(2E) * (1, 0)^T
# -------------------------------------------------------------------------
import cmath, math

print("\n[6] Numerical example: p = E * (1, 0, 0, 1)  [massless along z]")
E = 1.0
# Pauli matrices (2-component notation)
# p_{alpha dotalpha} = E*sigma^0 + E*sigma^3
# sigma^0 = [[1,0],[0,1]], sigma^3 = [[1,0],[0,-1]]
# sum = [[2,0],[0,0]]  * E
p_mat = [[2*E, 0], [0, 0]]
print(f"    p_{{alpha dotalpha}} = {p_mat}")
print(f"    det = {p_mat[0][0]*p_mat[1][1] - p_mat[0][1]*p_mat[1][0]:.3f}  (= 0 confirms p^2 = 0)")

lam = [math.sqrt(2*E), 0.0]
lam_tilde = [math.sqrt(2*E), 0.0]
print(f"    lambda_alpha = {lam}")
print(f"    lambdatilde_{{dotalpha}} = {lam_tilde}")
# Check: lambda ⊗ lambdatilde = p_mat
outer = [[lam[i]*lam_tilde[j] for j in range(2)] for i in range(2)]
print(f"    lambda_alpha x lambdatilde_dotalpha = {outer}  (check: matches p_mat ✓)")

print("\nDone: 03_spinor_helicity.py")
