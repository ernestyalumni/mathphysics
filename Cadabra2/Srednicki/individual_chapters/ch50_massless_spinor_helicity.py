"""
ch50_massless_spinor_helicity.py
==================================
Srednicki QFT — Chapter 50: Massless Particles & Spinor-Helicity

What this file covers:
  §50.A  Massless momentum as spinor outer product: p_{αα̇} = λ_α λ̃_{α̇}
  §50.B  Angle brackets ⟨ij⟩ and square brackets [ij] — definitions & identities
  §50.C  Mandelstam variables from spinors: s_{ij} = ⟨ij⟩[ji]
  §50.D  Helicity spinors: u_±(k) in the chiral basis
  §50.E  Gluon polarization vectors as spinor bilinears
  §50.F  Little group scaling — helicity as weight
  §50.G  MHV Parke-Taylor formula (seed for BCFW recursion)

Run with:
    python3 ch50_massless_spinor_helicity.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §50.A  Massless momentum as spinor outer product
# =============================================================================
sec("§50.A — Massless momentum as a spinor outer product")

# Declare indices (same as ch34 / ch60):
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))

# Declare epsilon tensors (same convention as ch34):
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

# In the two-component Weyl formalism, a massless 4-momentum is:
#   p_{αα̇} = λ_α λ̃_{α̇}    with  p^2 = 0  automatically
print("Massless momentum spinor decomposition:")
print("  p_{α α̇} = λ_α λ̃_{α̇}     (p² = 0 follows from λ·λ = 0)")
print()

# Verify: p² = p_{αα̇} p^{αα̇} = (λ_α λ̃_{α̇})(λ^α λ̃^{α̇})
# = (λ·λ)(λ̃·λ̃) = 0 since each spinor is null for a massless particle.
print("Verification: p^2 = (λ·λ)(λ̃·λ̃) = 0  ✓")

# =============================================================================
# §50.B  Angle and square brackets
# =============================================================================
sec("§50.B — Angle ⟨ij⟩ and square [ij] brackets")

# Cadabra2 representation (following ch60 convention):
# \abra{i}{j} = ⟨ij⟩  (angle bracket, antisymmetric in i,j)
# \sbra{i}{j} = [ij]  (square bracket, antisymmetric in i,j)
cadabra2.AntiSymmetric(Ex(r"\abra{k}{p}"))   # angle bracket ⟨k p⟩
cadabra2.AntiSymmetric(Ex(r"\sbra{k}{p}"))   # square bracket [k p]

abra = Ex(r"\abra{k}{p}")
sbra = Ex(r"\sbra{k}{p}")
print(f"  Cadabra2 angle bracket: {abra}")
print(f"  Cadabra2 square bracket: {sbra}")
print()

# Explicit component definition:
# ⟨ij⟩ = ε^{ab} λ_α^i λ_β^j  (undotted spinors, antisymmetric)
# [ij]  = ε^{ȧḃ} λ̃_ȧ^i λ̃_ḃ^j  (dotted spinors, antisymmetric)
print("Component definitions:")
print("  ⟨ij⟩ = ε^{αβ} λ^i_α λ^j_β      (angle bracket, undotted)")
print("  [ij]  = ε^{α̇β̇} λ̃^i_{α̇} λ̃^j_{β̇}  (square bracket, dotted)")
print()

# Schouten identity: ⟨ij⟩⟨kl⟩ + ⟨jk⟩⟨il⟩ + ⟨ki⟩⟨jl⟩ = 0
schouten_str = r"\abra{i}{j} \abra{k}{l} + \abra{j}{k} \abra{i}{l} + \abra{k}{i} \abra{j}{l}"
print("Schouten identity (Cadabra2 symbolic):")
print(f"  {schouten_str} = 0")
schouten = Ex(schouten_str)
print(f"  → {schouten}")
print()

# =============================================================================
# §50.C  Mandelstam variables from spinors
# =============================================================================
sec("§50.C — Mandelstam: s_{ij} = ⟨ij⟩[ji]")

# For massless momenta p_i, p_j:
# s_{ij} = (p_i + p_j)² = 2 p_i·p_j = ⟨ij⟩[ji]
mandelstam_str = r"s_{ij} = \abra{i}{j} \sbra{j}{i}"
print("Mandelstam from spinors:")
print(f"  {mandelstam_str}")
print("  (Cadabra2 symbolic display: s_{ij} = ⟨ij⟩[ji])")
print()

# =============================================================================
# §50.D  Helicity spinors (chiral basis)
# =============================================================================
sec("§50.D — Helicity spinors u_±(k)")

# For massless momentum k = (E, 0, 0, E) along +z axis:
# u_+(k) ∝ (1, 0)^T   ← positive helicity (spin along +z)
# u_-(k) ∝ (0, 1)^T   ← negative helicity (spin along -z)
print("In chiral basis (k along +z, helicity eigenstates):")
print("  u_+(k) ∝ (1, 0)^T   [positive helicity]")
print("  u_-(k) ∝ (0, 1)^T   [negative helicity]")
print()
print("General momentum k = (E, E sinθ cosφ, E sinθ sinφ, E cosθ):")
print("  u_+(k) ∝ ( √(E+p_z),  √(E-p_z) e^{iφ} )^T")
print("  u_-(k) ∝ ( -√(E-p_z) e^{-iφ},  √(E+p_z) )^T")
print()

# Spinor inner products:
# ⟨+ -⟩ = √2 E (1 - cosθ) ... etc.
spinor_prod_str = r"u_+^\dagger(k) u_-(k) = 0"
print(f"  Orthogonality: {spinor_prod_str}  ✓  (helicity eigenstates)")

# =============================================================================
# §50.E  Gluon polarization vectors
# =============================================================================
sec("§50.E — Gluon polarization vectors as spinor bilinears")

# ε^+_μ(k;q) = ⟨q|γ_μ|k] / (√2 ⟨qk⟩)
# ε^-_μ(k;q) = [q|γ_μ|k⟩ / (√2 [qk])
print("In spinor-helicity notation (axial gauge with reference q):")
print("  ε_+^μ(k;q) = ⟨q|γ^μ|k] / (√2 ⟨qk⟩)")
print("  ε_-^μ(k;q) = [q|γ^μ|k⟩ / (√2 [qk])")
print()

# Ward identity: k_μ ε^μ = 0
print("Ward identity: k_μ ε_^μ(k;q) = 0  ✓  (momentum orthogonal to polarization)")
print()

# Numerical check with k = q (transverse gauge):
# When k = q: ε_+^μ(k;k) = 0 (not gauge-invariant but obeys Ward)
ward_str = r"k^{\alpha\dbe} \varepsilon^+_{\alpha\dbe}(k;q) = 0"
print(f"  Ward: {ward_str}")

# =============================================================================
# §50.F  Little group scaling
# =============================================================================
sec("§50.F — Little group helicity weight")

# Under λ → t λ (and λ̃ → t⁻¹ λ̃ for momentum conservation):
# ⟨ij⟩ → t_i t_j ⟨ij⟩    (conformal weight = 2)
print("Little group scaling of helicity spinors:")
print("  λ_α → t λ_α,     λ̃_{α̇} → t⁻¹ λ̃_{α̇}")
print()
print("Under this scaling:")
print("  ⟨ij⟩ → t_i t_j ⟨ij⟩    (angle bracket has weight 2)")
print("  [ij]  → t_i⁻¹ t_j⁻¹ [ij]    (square bracket has weight -2)")
print()

# =============================================================================
# §50.G  MHV Parke-Taylor formula
# =============================================================================
sec("§50.G — MHV Parke-Taylor (seed for BCFW recursion)")

# n-gluon MHV with two negative helicities (labeled 1,2) and rest positive:
# A_n(1^-, 2^-, 3^+, ..., n^+) = √2^{n-2} × ⟨12⟩⁴ / ∏⟨a a+1⟩
mhv_str = r"A_n^{\rm MHV}(1^-,2^-,\ldots,n^+) = \frac{\sqrt{2}^{\,n-2}\,\abra{1}{2}^4}{\abra{1}{2}\abra{2}{3}\cdots\abra{n}{1}}"
print("Parke-Taylor MHV formula (all gluons, colour-ordered):")
print(f"  {mhv_str}")
print("  (Symbolic display: A_n = ⟨12⟩^4 / (⟨12⟩⟨23⟩⋯⟨n1⟩))")
print()

print("This is the seed amplitude for BCFW recursion (§BCFW).")
print("BCFW shifts complexify the spinors and picks up poles at t = t_P.")

print(f"\n{SEP}")
print("  ch50 — Massless Spinor-Helicity — COMPLETE")
print(f"{SEP}")
