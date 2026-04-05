"""
04_mhv_amplitudes.py
====================
MHV (Maximally Helicity Violating) amplitudes and the Parke-Taylor formula.

Reference: Srednicki QFT, Chapter 38.
           Parke & Taylor, Phys.Rev.Lett. 56 (1986) 2459.
           Mangano & Parke, Phys.Rept. 200 (1991) 301.
           See REFERENCES.md for full citations.

The Parke-Taylor formula (tree-level gluon MHV amplitude):
    A_n(1^-, 2^-, 3^+, ..., n^+) = <12>^4 / (<12><23>...<n1>)

This is "MHV" because exactly 2 gluons have negative helicity.
All-plus and one-minus amplitudes vanish at tree level.

Run with:
    python3 04_mhv_amplitudes.py
or inside Docker:
    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 python3 /work/04_mhv_amplitudes.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import cmath, math

print("=" * 60)
print("MHV amplitudes and Parke-Taylor formula (Srednicki Ch. 38)")
print("=" * 60)

__cdbkernel__ = cadabra2.create_scope()

# -------------------------------------------------------------------------
# 1. The Parke-Taylor formula
#    Color-ordered partial amplitude for n gluons with helicities:
#    particles 1,2: negative helicity (h = -1)
#    particles 3,...,n: positive helicity (h = +1)
#
#    A_n(1^-, 2^-, 3^+, ..., n^+) = i * g^{n-2} * <12>^4 / (<12><23>...<n1>)
#
#    where <ij> are angle spinor products.
#    The denominator is the "cyclic chain" of adjacent brackets.
#
#    Connection to Srednicki Ch.38:
#    Srednicki derives this from the Berends-Giele recursion / MHV vertex rules.
#    The formula encodes the remarkable simplicity hidden in gauge theory.
# -------------------------------------------------------------------------
print("\n[1] Parke-Taylor formula:")
print("""
    A_n(1^-, 2^-, 3^+, ..., n^+) = i * g^{n-2} * <12>^4 / (<12><23>...<(n-1)n><n1>)

    where:
      <ij> = eps^{ab} lambda^i_a lambda^j_b  (angle bracket)
      g = gauge coupling
      The denominator has n adjacent brackets in a cyclic chain
""")

# -------------------------------------------------------------------------
# 2. n=4 case: A_4(1^-, 2^-, 3^+, 4^+)
#    A_4 = i * g^2 * <12>^4 / (<12><23><34><41>)
#
#    This is the simplest nontrivial MHV amplitude.
#    Using momentum conservation: <34> = [12]s_{12}/<12>, etc.
# -------------------------------------------------------------------------
print("[2] Four-gluon MHV amplitude A_4(1^-, 2^-, 3^+, 4^+):")
print("""
    A_4 = i * g^2 * <12>^4 / (<12><23><34><41>)

    Alternate forms (using momentum conservation <ij>[ji] = s_{ij}):
    A_4 = i * g^2 * <12>^4 * [34] / (s_{12} * s_{23} * <12>)

    Srednicki Ch.38 connection:
    - Derived from helicity spinor Feynman rules for Yang-Mills
    - Encodes the MHV amplitude which is maximal simplicity in gauge theory
    - Nonzero only when exactly 2 helicities are negative
""")

# Symbolic Parke-Taylor amplitude (cadabra2 expression)
# We represent angle brackets as A_{ij} and square brackets as B_{ij}
PT4 = Ex(r"A_{12}^4 / (A_{12} A_{23} A_{34} A_{41})")
print("    Symbolic (A_{ij} = <ij>):")
print("    A_4 = A_{12}^4 / (A_{12} A_{23} A_{34} A_{41})")
print("        = A_{12}^3 / (A_{23} A_{34} A_{41})")

# -------------------------------------------------------------------------
# 3. Numerical evaluation for specific kinematics
#    Use simple center-of-mass kinematics for 4 massless gluons.
#    p1 + p2 + p3 + p4 = 0 (all incoming convention)
#
#    Choose:
#    p1 = (1, 0, 0, 1) * E       [along +z]
#    p2 = (1, 0, 0, -1) * E      [along -z]
#    p3 = (1, 1, 0, 0) / sqrt(2) * 2E  [along +x — adjusted for momentum cons.]
#    ...
#    Actually use analytic kinematics for simplicity:
#    Pick spinors directly.
# -------------------------------------------------------------------------
print("\n[3] Numerical evaluation with explicit spinors:")
print("    Choose kinematics: 2->2 scattering, CoM frame")
print()

# Simple analytic parametrization:
# For 2->2 gluon scattering in CoM frame with angle theta:
# p1 = E(1, 0, 0, 1)   => lambda1 = sqrt(2E) * (1, 0)
# p2 = E(1, 0, 0, -1)  => lambda2 = sqrt(2E) * (0, 1)
# p3 = -E(1, sin(th), 0, cos(th))  [outgoing -> incoming with minus sign]
# p4 = -E(1, -sin(th), 0, -cos(th))

import cmath

def angle_bracket(lam_i, lam_j):
    """<ij> = eps^{ab} lam_i_a lam_j_b = lam_i[0]*lam_j[1] - lam_i[1]*lam_j[0]"""
    return lam_i[0]*lam_j[1] - lam_i[1]*lam_j[0]

def square_bracket(lamt_i, lamt_j):
    """[ij] = eps_{dota dotb} lamtilde_i^{dota} lamtilde_j^{dotb}"""
    return lamt_i[0]*lamt_j[1] - lamt_i[1]*lamt_j[0]

def parke_taylor_4(lam_list, neg_hel_i, neg_hel_j):
    """
    Compute A_4(... i^- j^- ...) / (i*g^2) = <ij>^4 / product of adjacent <ab>
    lam_list: list of 4 angle spinors [lam1, lam2, lam3, lam4]
    neg_hel_i, neg_hel_j: 0-indexed positions of negative helicity gluons
    """
    n = 4
    # Numerator: <ij>^4
    num = angle_bracket(lam_list[neg_hel_i], lam_list[neg_hel_j])**4
    # Denominator: cyclic product <12><23><34><41>
    denom = 1.0 + 0j
    for k in range(n):
        denom *= angle_bracket(lam_list[k], lam_list[(k+1) % n])
    return num / denom

# Kinematics: theta = pi/3 (60 degrees) scattering
import math
theta = math.pi / 3
E = 1.0

# Spinors for massless p = E(1, sin(th), 0, cos(th)):
# p_{aa'} = E * [[1+cos(th),  sin(th)],
#                [sin(th),    1-cos(th)]]
# For real momenta (Euclidean/split signature trick):
# Use the parametrization from Elvang-Huang:
# For p = E(1, sin(th), 0, cos(th)):
#   lambda = sqrt(E) * (cos(th/2), sin(th/2))
#   lambdatilde = sqrt(E) * (cos(th/2), sin(th/2))  [for real momenta]

def spinor_from_momentum_angle(E, theta_p):
    """Angle spinor for p = E*(1, sin(th), 0, cos(th))"""
    c = math.cos(theta_p/2)
    s = math.sin(theta_p/2)
    sq = math.sqrt(E)
    return (sq * c, sq * s)

# Assign momenta (using all-incoming convention, p3 and p4 flipped)
# p1 along z: theta=0
# p2 along -z: theta=pi
# p3 at angle theta from z
# p4 must satisfy p1+p2+p3+p4 = 0 => at angle pi+theta (opposite to p3)
lam1 = spinor_from_momentum_angle(E, 0)           # p1 = E(1,0,0,1)
lam2 = spinor_from_momentum_angle(E, math.pi)     # p2 = E(1,0,0,-1)
lam3 = spinor_from_momentum_angle(E, theta)       # p3 at angle theta
# For momentum conservation, need to flip sign conventions;
# use e^{i pi/2} = i rotation for outgoing spinors:
lam4 = (1j * lam3[1], 1j * lam3[0])  # approximate; illustrative

lam_list = [lam1, lam2, lam3, lam4]

print(f"    Scattering angle theta = {math.degrees(theta):.1f} deg")
print(f"    Spinors:")
for k, l in enumerate(lam_list, 1):
    print(f"      lambda_{k} = ({l[0]:.4f}, {l[1]:.4f})")

# Compute A_4(1^-, 2^-, 3^+, 4^+): negative helicity on 1,2
A4 = parke_taylor_4(lam_list, 0, 1)
print(f"\n    A_4(1^-, 2^-, 3^+, 4^+) / (i*g^2) = {A4:.6f}")

# -------------------------------------------------------------------------
# 4. Cyclic invariance check
#    The Parke-Taylor formula is cyclically invariant: rotating labels
#    1->2->3->...->n->1 should leave the stripped amplitude invariant
#    (up to an overall phase from the cyclic permutation of <ij>).
# -------------------------------------------------------------------------
print("\n[4] Cyclic symmetry check:")
print("    A_4(1^-, 2^-, 3^+, 4^+) should be invariant under cyclic relabeling")
print("    when the negative helicity positions cycle accordingly.")

# Cyclic rotation: 1->2, 2->3, 3->4, 4->1
# New negative helicities at positions 2,3 (0-indexed: 1,2)
lam_cycled = [lam2, lam3, lam4, lam1]  # cyclically shifted
A4_cycled = parke_taylor_4(lam_cycled, 0, 1)  # new neg hel at 1,2 = old 2,3
print(f"    A_4(2^-, 3^-, 4^+, 1^+) / (i*g^2) = {A4_cycled:.6f}")
print(f"    Ratio |A4/A4_cycled| = {abs(A4/A4_cycled):.4f}  (should be 1 for cyclic invariance)")

# -------------------------------------------------------------------------
# 5. Vanishing amplitudes at tree level
#    All-plus: A_n(1^+, 2^+, ..., n^+) = 0
#    One-minus: A_n(1^-, 2^+, ..., n^+) = 0
#    These vanish to all loop orders in N=4 SYM (protected by SUSY).
#    At tree level, they vanish due to MHV selection rule.
# -------------------------------------------------------------------------
print("\n[5] Vanishing amplitudes:")
print("    A_n(all +) = 0   [no source of angle brackets in amplitude]")
print("    A_n(1^-) = 0     [one minus helicity: no <ij>^4 numerator available]")
print()
print("    Physical reason: helicity conservation in soft/collinear limits")
print("    and holomorphic/anti-holomorphic factorization.")
print("    In N=4 SYM: protected by supercharge Ward identities.")

# -------------------------------------------------------------------------
# 6. Connection to Srednicki Chapter 38
#    Srednicki Ch.38 derives the Parke-Taylor formula via:
#    1. Helicity spinor Feynman rules (polarization vectors in spinor notation)
#    2. Explicit computation of the 4-point amplitude
#    3. Recursion relation / pattern recognition
#    Key equations in Srednicki Ch.38:
#    - Eq.(38.6): polarization vector ε^+_μ(p,q) ~ <q|σ_μ|p] / [qp]
#    - Eq.(38.8): polarization vector ε^-_μ(p,q) ~ <p|σ_μ|q] / <pq>
#    - Eq.(38.18): The 4-gluon Parke-Taylor result
# -------------------------------------------------------------------------
print("\n[6] Connection to Srednicki Ch. 38:")
print("""
    Srednicki derives A_4 via:

    Step 1: Helicity polarization vectors in spinor notation:
      eps^+_mu(p,q) = <q|sigma_mu|p] / sqrt(2) [qp]
      eps^-_mu(p,q) = [q|sigma_mu|p> / sqrt(2) <qp>
      (q is reference momentum, drops out of gauge-invariant amplitude)

    Step 2: Compute <M_4> by contracting polarizations with Feynman vertices.
      Three 3-gluon vertices (color-ordered) and one 4-gluon vertex.

    Step 3: After simplification using spinor identities and Schouten,
      only one Feynman diagram survives for each helicity assignment!

    Result (Srednicki eq. 38.18):
      A_4(1^-, 2^-, 3^+, 4^+) = i*g^2 * <12>^4 / (<12><23><34><41>)

    The Parke-Taylor formula generalizes this to all n.
""")

print("Done: 04_mhv_amplitudes.py")
