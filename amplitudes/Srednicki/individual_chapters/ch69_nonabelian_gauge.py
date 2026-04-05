"""
ch69_nonabelian_gauge.py
==========================
Srednicki QFT — Chapter 69: Nonabelian Gauge Theory

What this file covers:
  §69.A  SU(N) gauge fields and generators T^a
  §69.B  Covariant derivative D_μ = ∂_μ + ig A_μ
  §69.C  Field strength F^a_{μν}
  §69.D  Nonabelian Bianchi identity
  §69.E  Gauge-invariant Lagrangian L = -¼ F^a_{μν} F^{a μν}
  §69.F  QCD Lagrangian and quark-gluon vertices

Run with:
    python3 ch69_nonabelian_gauge.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §69.A  SU(N) gauge fields and generators
# =============================================================================
sec("§69.A — SU(N) gauge fields and generators")

print("Gauge field A_μ = A^a_μ T^a  (a = 1,..., N²-1 for SU(N))")
print()
print("Generator properties:")
print("  [T^a, T^b] = i f^{abc} T^c    (structure constants)")
print("  Tr[T^a T^b] = ½ δ^{ab}")
print("  (T^a)_ij (T^a)_kl = ½(δ_i^l δ_j^k - (1/N) δ_i^k δ_j^l)  [Fierz for SU(N)]")
print()
print("SU(3) generators (Gell-Mann matrices):")
print("  λ¹ = [[0,1,0],[1,0,0],[0,0,0]],  λ² = [[0,-i,0],[i,0,0],[0,0,0]],  λ³ = [[1,0,0],[0,-1,0],[0,0,0]]")
print("  ... (8 Gell-Mann matrices for SU(3))")

# =============================================================================
# §69.B  Covariant derivative
# =============================================================================
sec("§69.B — Covariant derivative D_μ = ∂_μ + ig A_μ")

print("Gauge covariant derivative:")
print("  D_μ = ∂_μ + ig A_μ = ∂_μ + ig A^a_μ T^a")
print()
print("Quark covariant derivative:")
print("  D_μ q = ∂_μ q + ig A^a_μ T^a q")
print()
print("Transformation law:")
print("  A_μ → U A_μ U⁻¹ + (i/g) U ∂_μ U⁻¹")
print("  D_μ → U D_μ U⁻¹   (transforms homogeneously)")

# =============================================================================
# §69.C  Nonabelian field strength
# =============================================================================
sec("§69.C — Nonabelian field strength F^a_{μν}")

print("Nonabelian field strength (curvature):")
print("  F_μν = ∂_μ A_ν - ∂_ν A_μ + ig[A_μ, A_ν]")
print("  F^a_μν = ∂_μ A^a_ν - ∂_ν A^a_μ - g f^{abc} A^b_μ A^c_ν")
print()
print("Unlike U(1): F_μν ≠ ∂_μ A_ν - ∂_ν A_μ  (has self-interaction!)")
print()
print("The cubic and quartic gluon terms arise from the commutator [A_μ, A_ν].")
print("This is what makes QCD nonabelian and gives gluon self-interactions.")

# =============================================================================
# §69.D  Nonabelian Bianchi identity
# =============================================================================
sec("§69.D — Nonabelian Bianchi identity")

print("Bianchi identity:")
print("  D_{[μ} F^a_{νρ]} = 0   (analogous to ∂_{[μ} F_{νρ]} = 0 in U(1))")
print()
print("The D_ covariant derivative ensures gauge covariance:")
print("  D_ρ F^a_{μν} + D_μ F^a_{νρ} + D_ν F^a_{ρμ} = 0")
print()
print("This is a key identity that constrains the gluon vertex tensor.")

# =============================================================================
# §69.E  Gauge-invariant Lagrangian
# =============================================================================
sec("§69.E — Gauge-invariant Lagrangian L = -¼ F^a_{μν} F^{a μν}")

print("Yang-Mills Lagrangian:")
print("  L_YM = -¼ F^a_{μν} F^{a μν}")
print("        = -½ Tr[F_{μν} F^{μν}]")
print()
print("Variation → Yang-Mills equations:")
print("  D^μ F^a_{μν} = g J^a_ν   (nonabelian Gauss's law)")
print()
print("Expanding the Lagrangian:")
print("  L = -½ (∂_μ A^a_ν ∂^μ A^a_ν - ∂_μ A^a_ν ∂^ν A^a_μ)  [kinetic]")
print("      + g f^{abc} (∂_μ A^a_ν) A^b_μ A^c_ν               [3-gluon]")
print("      - ¼ g² f^{abc} f^{ade} A^b_μ A^c_ν A^d_μ A^e_ν  [4-gluon]")

# =============================================================================
# §69.F  QCD Lagrangian and quark-gluon vertices
# =============================================================================
sec("§69.F — QCD Lagrangian and quark-gluon vertices")

print("QCD Lagrangian (quarks + gluons):")
print("  L_QCD = ψ̄(i γ^μ D_μ - m)ψ - ¼ F^a_{μν} F^{a μν}")
print()
print("Quark-gluon vertex (Feynman rule):")
print("  -ig γ_μ T^a    (color matrix T^a acting on quark color index)")
print()
print("Ghost-gluon vertex:  g f^{abc} p̄_μ   (from Fadeev-Popov)")
print()
print("Color flow in Feynman diagrams:")
print("  Gluon lines carry adjoint color index a = 1,...,8")
print("  Quark lines carry fundamental color i = 1,2,3 (SU(3))")
print("  Color factor: Tr[T^a T^b] = ½ δ^{ab},  f^{abc} f^{abc} = N_c δ^{ab}")

print(f"\n{SEP}")
print("  ch69 — Nonabelian Gauge Theory — COMPLETE")
print(f"{SEP}")
