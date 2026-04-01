"""
ch50_spinor_helicity_core.py
========================
Srednicki QFT — Chapter 50: Massless Particles & Spinor-Helicity

THIS IS THE HEART OF MHV RESEARCH.

What this file covers:
  §50.A  Massless momentum as spinor outer product p_{αα̇} = λ_α λ̃_{α̇}
  §50.B  Little group and helicity
  §50.C  Angle brackets ⟨ij⟩ and square brackets [ij]
  §50.D  Spinor products and dot products: p_i · p_j = ½ ⟨ij⟩[ji]
  §50.E  Helicity spinors u_±(k), v_±(k)
  §50.F  Polarization vectors ε^μ_±(k;q)
  §50.G  Gauge invariance and little group scaling

Run with:
    python3 ch50_spinor_helicity_core.py
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

cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))

cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

# Declare angle/square brackets as in ch60:
cadabra2.AntiSymmetric(Ex(r"\abra{k}{p}"))   # angle bracket ⟨kp⟩
cadabra2.AntiSymmetric(Ex(r"\sbra{k}{p}"))   # square bracket [kp]

print("Key decomposition (Srednicki eq. 50.1):")
print("  p_{α α̇} = λ_α λ̃_{α̇}    with  p² = 0  automatically")
print()
print("This is the foundational identity of the spinor-helicity formalism.")
print("Every massless 4-momentum is a rank-1 outer product of two Weyl spinors.")
print()

# Cadabra2 symbolic:
print("  Symbolic: p_{alpha alphȧ} = λ_α λ̃_{α̇} = angle bracket form")

# =============================================================================
# §50.B  Little group and helicity
# =============================================================================
sec("§50.B — Little group and helicity")

print("Little group for massless particles: U(1) (for real momenta)")
print()
print("Under λ → t λ (t ∈ ℂ*), the momentum is invariant:")
print("  p = λ ⊗ λ̃  →  (t λ) ⊗ (t⁻¹ λ̃) = λ ⊗ λ̃ = p")
print()
print("Helicity h: eigenvalue of Ŝ_z / |p|:")
print("  h = +½  for u_+(p) ∝ (1,0)  (right-handed, positive helicity)")
print("  h = -½  for u_-(p) ∝ (0,1)  (left-handed, negative helicity)")
print()
print("Under λ → t λ, λ̃ → t⁻¹ λ̃:")
print("  Angle bracket: ⟨ij⟩ → t_i t_j ⟨ij⟩     (weight +1 per λ)")
print("  Square bracket: [ij]  → t_i⁻¹ t_j⁻¹ [ij]   (weight -1 per λ̃)")

# =============================================================================
# §50.C  Angle and square brackets
# =============================================================================
sec("§50.C — Angle ⟨ij⟩ and square [ij] brackets")

abra = Ex(r"\abra{k}{p}")
sbra = Ex(r"\sbra{k}{p}")
print(f"  Cadabra2 angle bracket: {abra}")
print(f"  Cadabra2 square bracket: {sbra}")
print()

print("Definitions (Srednicki eq. 50.16):")
print("  ⟨ij⟩ = ε^{αβ} λ_α^i λ_β^j     (undotted spinors)")
print("  [ij]  = ε^{α̇β̇} λ̃_α̇^i λ̃_β̇^j  (dotted spinors)")
print()

print("Antisymmetry:")
print("  ⟨ij⟩ = -⟨ji⟩,   [ij] = -[ji]")
print("  ⟨ii⟩ = [ii] = 0  (null spinor)")
print()
print("Schouten identity:")
print("  ⟨ij⟩⟨kl⟩ + ⟨jk⟩⟨il⟩ + ⟨ki⟩⟨jl⟩ = 0")
schouten = Ex(r"\abra{i}{j} \abra{k}{l} + \abra{j}{k} \abra{i}{l} + \abra{k}{i} \abra{j}{l}")
print(f"  → {schouten}")

# =============================================================================
# §50.D  Dot products from spinors
# =============================================================================
sec("§50.D — Dot products: 2 p_i·p_j = ⟨ij⟩[ji]")

print("Key identity (Srednicki eq. 50.17):")
print("  2 p_i · p_j = ⟨ij⟩ [ji]    (massless momenta)")
print()
print("Special cases:")
print("  p_i² = 0  ✓  (since ⟨ii⟩ = [ii] = 0)")
print("  p_i · p_i = 0")
print()

print("Mandelstam variables:")
print("  s_{ij} = (p_i + p_j)² = 2 p_i · p_j = ⟨ij⟩[ji]")
print("  t_{ij} = (p_i - p_j)² = -⟨i j̃⟩[j ĩ]  etc.")
print("  u_{ij} = -⟨i j⟩[i j] - ⟨i j̃⟩[i j̃]  (momentum conservation)")

# =============================================================================
# §50.E  Helicity spinors
# =============================================================================
sec("§50.E — Helicity spinors u_±(k), v_±(k)")

print("Massless helicity spinors (positive-energy solutions):")
print("  u_+(p)  ∝ ( λ_α , 0 )^T   ← right-handed, h = +½")
print("  u_-(p)  ∝ ( 0 ,  λ̃_α̇ )^T  ← left-handed, h = -½")
print()
print("For momentum along +z (p = (E,0,0,E)):")
print("  u_+(p) = √(2E) (1, 0)^T")
print("  u_-(p) = √(2E) (0, 1)^T")
print()
print("Inner products (orthonormality):")
print("  ū_+ u_+ = 2E,   ū_- u_- = 2E,   ū_+ u_- = 0")
print("  [±|±] = 2E,   [±|∓] = 0   (square bracket inner products)")

# =============================================================================
# §50.F  Polarization vectors
# =============================================================================
sec("§50.F — Gluon polarization vectors ε^μ_±(k;q)")

print("Gauge field polarization vectors (axial gauge, reference q):")
print()
print("  ε^μ_+(k;q) = ⟨q| γ^μ |k] / (√2 ⟨qk⟩)")
print("  ε^μ_-(k;q) = [q| γ^μ |k⟩ / (√2 [qk])")
print()
print("Properties:")
print("  k_μ ε^μ_±(k;q) = 0           (transverse)")
print("  q_μ ε^μ_±(k;q) = 0           (reference gauge)")
print("  ε^μ_± ε_μ^± = -1")
print()
print("Ward identity: Σ_± ε^μ_± ε^ν_± = -η^{μν} + k^μ q^ν/(k·q) + k^ν q^μ/(k·q)")

# =============================================================================
# §50.G  Gauge invariance and little group
# =============================================================================
sec("§50.G — Gauge invariance and little group scaling")

print("Under gauge transformation ε_μ → ε_μ + α k_μ:")
print("  The spinor form is manifestly gauge-invariant because k/|k⟩ projects out the gauge term.")
print()
print("Little group scaling of amplitude:")
print("  λ_i → t_i λ_i,  λ̃_i → t_i⁻¹ λ̃_i")
print("  A_n → t_1^{-2h₁} ... t_n^{-2h_n} A_n")
print()
print("For a gluon with helicity h = ±1:")
print("  Numerator: ε_± · J ∝ ⟨q| γ·J |k] or [q| γ·J |k⟩")
print("  Under λ_k → t λ_k: ⟨qk⟩ → t ⟨qk⟩, [qk] → t⁻¹ [qk]")
print("  Net: h = +1: ⟨qk⟩⁻¹ → t⁻¹ ⟨qk⟩⁻¹ (weight -1)")
print("       h = -1: [qk]⁻¹ → t [qk]⁻¹ (weight +1)")

print(f"\n{SEP}")
print("  ch50 — Spinor-Helicity (CORE) — COMPLETE")
print(f"{SEP}")
