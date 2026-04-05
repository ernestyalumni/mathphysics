"""
ch70_group_representations.py
==============================
Srednicki QFT — Chapter 70: Group Representations (SU(N))

What this file covers:
  §70.A  SU(N) fundamental and adjoint representations
  §70.B  Casimir operators and Dynkin indices
  §70.C  Color factor algebra for QCD amplitudes
  §70.D  Fierz identities for SU(N) generators
  §70.E  Decomposition of products of representations

Run with:
    python3 ch70_group_representations.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cddkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §70.A  SU(N) fundamental and adjoint representations
# =============================================================================
sec("§70.A — SU(N) fundamental and adjoint representations")

print("SU(N) group: N² - 1 generators T^a")
print()
print("Fundamental representation (quarks, N = 3 for QCD):")
print("  Dimension: N")
print("  Generators: (T^a)_ij  (N×N Hermitian matrices)")
print("  Normalization: Tr[T^a T^b] = ½ δ^{ab}")
print()
print("Adjoint representation (gluons, for QCD N_c = 3):")
print("  Dimension: N² - 1")
print("  Generators: (T^a)_bc = -i f^{abc}  (structure constants)")
print("  Fierz identity: (T^a)_ij (T^a)_kl = ½(δ_i^l δ_j^k - (1/N) δ_i^k δ_j^l)")

# =============================================================================
# §70.B  Casimir operators
# =============================================================================
sec("§70.B — Casimir operators and Dynkin indices")

print("Quadratic Casimir:")
print("  C_2(R) T^a T^a = C_2(R) I_R   (proportional to identity in rep R)")
print()
print("Values:")
print("  C_2(F) = (N²-1)/(2N)    [fundamental]")
print("  C_2(A) = N_c               [adjoint, for N_c = 3: C_2(A) = 3]")
print()
print("Trace anomaly coefficient:")
print("  T(R) δ^{ab} = Tr[T^a T^b]  (normalization: T(F) = ½ for fundamental)")
print("  T(A) = N_c")

# =============================================================================
# §70.C  Color factor algebra for QCD amplitudes
# =============================================================================
sec("§70.C — Color factor algebra for QCD amplitudes")

print("Color structures in QCD:")
print()
print("  Tr[T^a T^b]       = ½ δ^{ab}    (U(1) factor, Abelian part)")
print("  d^{abc} = 2 Tr[{T^a,T^b}T^c]  (symmetric d-symbol)")
print()
print("Important identities:")
print("  f^{abc} f^{abd} = N_c δ^{cd}")
print("  d^{abc} d^{abd} = (N_c² - 4)/N_c δ^{cd}")
print("  f^{abc} d^{abd} = 0")
print()
print("For n-gluon tree amplitudes: color decomposition into:")
print("  - Single-trace (primitive) terms: Tr[T^{a1}...T^{an}]")
print("  - Multi-trace terms (suppressed at large N_c)")

# =============================================================================
# §70.D  Fierz identities for SU(N)
# =============================================================================
sec("§70.D — Fierz identities for SU(N) generators")

print("SU(N) Fierz rearrangement:")
print()
print("  (T^a)_ij (T^a)_kl = ½(δ_i^l δ_j^k - (1/N) δ_i^k δ_j^l)")
print()
print("Applications:")
print("  - qq̄ → gg color factor: C_F = (N_c²-1)/(2N_c) = 4/3 for N_c=3")
print("  - gg → gg color factors: C_A = N_c = 3, C_F = 4/3")
print()
print("qq̄ annihilation color structure:")
print("  |M_qq̄→gg|² ∝ C_F C_A  (interference between s-channel and t/u-channel)")

# =============================================================================
# §70.E  Representation products
# =============================================================================
sec("§70.E — Decomposition of products of representations")

print("SU(N) representation theory:")
print()
print("  3 ⊗ 3̄ = 1 ⊕ 8        (quark-antiquark = singlet + octet)")
print("  3 ⊗ 3 = 3̄ ⊕ 6         (two quarks = antisextet + triplet)")
print("  8 ⊗ 8 = 1 ⊕ 8 ⊕ 8 ⊕ 10 ⊕ 10̄ ⊕ 27")
print()
print("For QCD amplitudes:")
print("  - Color flows along quark lines")
print("  - Gluon exchange between quark lines gives color factors C_F, C_A")
print("  - Four-gluon vertex gives C_A² terms")

print(f"\n{SEP}")
print("  ch70 — Group Representations — COMPLETE")
print(f"{SEP}")
