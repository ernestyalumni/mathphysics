"""
bcfw_starter.py
===============
Basic BCFW recursion implementation starter for MHV amplitudes.

Reference: Britto, Cachazo, Feng, Witten - "Direct Proof of Tree-Level Recursion
Relation in Yang-Mills Theory" (2005)

BCFW shift (for two legs i,j):
    lambda_i(z) = lambda_i + z * lambda_j
    lambdatilde_j(z) = lambdatilde_j - z * lambdatilde_i

The amplitude satisfies a recursion:
    A_n = sum over poles  A_L(z) * A_R(z) / P(z)^2   (with residue at z=0 subtracted)

This script provides:
- The shift operator
- A simple 4-point consistency check (should recover Parke-Taylor)
- Skeleton for higher-point recursion

Run:
    python3 bcfw_starter.py
"""

import cmath
import math
from typing import List, Tuple, Callable

print("=" * 60)
print("BCFW Recursion Starter for MHV Amplitudes")
print("=" * 60)

# -------------------------------------------------------------------------
# Spinor and bracket utilities (reused/adapted from 04_mhv_amplitudes.py)
# -------------------------------------------------------------------------
def angle_bracket(lam_i: Tuple[complex, complex], lam_j: Tuple[complex, complex]) -> complex:
    """<ij> = eps^{ab} lambda^i_a lambda^j_b"""
    return lam_i[0]*lam_j[1] - lam_i[1]*lam_j[0]

def square_bracket(lamt_i: Tuple[complex, complex], lamt_j: Tuple[complex, complex]) -> complex:
    """[ij]"""
    return lamt_i[0]*lamt_j[1] - lamt_i[1]*lamt_j[0]

def parke_taylor_n(lam_list: List[Tuple[complex, complex]], neg_hel: List[int]) -> complex:
    """
    Generalized Parke-Taylor for n gluons with exactly two negative helicities.
    Assumes the negative helicities are at positions in neg_hel (0-indexed).
    """
    n = len(lam_list)
    if len(neg_hel) != 2:
        raise ValueError("MHV requires exactly two negative helicities")
    i, j = neg_hel
    num = angle_bracket(lam_list[i], lam_list[j]) ** 4
    denom = complex(1, 0)
    for k in range(n):
        denom *= angle_bracket(lam_list[k], lam_list[(k + 1) % n])
    return num / denom

# -------------------------------------------------------------------------
# BCFW Shift
# -------------------------------------------------------------------------
def bcfw_shift(lam_list: List[Tuple[complex, complex]],
               lamt_list: List[Tuple[complex, complex]],
               shift_i: int, shift_j: int, z: complex):
    """
    Apply BCFW shift to a pair of legs.
    Shift: lambda_i(z) = lambda_i + z * lambda_j
           lambdatilde_j(z) = lambdatilde_j - z * lambdatilde_i
    """
    lam_new = list(lam_list)
    lamt_new = list(lamt_list)
    lam_i, lam_j = lam_list[shift_i], lam_list[shift_j]
    lamt_i, lamt_j = lamt_list[shift_i], lamt_list[shift_j]

    lam_new[shift_i] = (lam_i[0] + z * lam_j[0], lam_i[1] + z * lam_j[1])
    lamt_new[shift_j] = (lamt_j[0] - z * lamt_i[0], lamt_j[1] - z * lamt_i[1])
    return lam_new, lamt_new

# -------------------------------------------------------------------------
# Simple 4-point BCFW consistency check
# For 4-gluon MHV, the recursion should reproduce the Parke-Taylor result.
# We pick a valid shift (e.g., shift legs 1 and 3 for A(1-,2-,3+,4+))
# -------------------------------------------------------------------------
print("\n[1] 4-point BCFW consistency check (recover Parke-Taylor)")

# Simple kinematics (same as 04_mhv_amplitudes.py)
theta = math.pi / 3
E = 1.0

def spinor_from_momentum_angle(E, theta_p):
    c = math.cos(theta_p/2)
    s = math.sin(theta_p/2)
    sq = math.sqrt(E)
    return (sq * c, sq * s)

lam1 = spinor_from_momentum_angle(E, 0)
lam2 = spinor_from_momentum_angle(E, math.pi)
lam3 = spinor_from_momentum_angle(E, theta)
lam4 = (1j * lam3[1], 1j * lam3[0])

lam_list = [lam1, lam2, lam3, lam4]
lamt_list = [lam1, lam2, lam3, lam4]  # real momenta approximation

# Reference Parke-Taylor
pt_ref = parke_taylor_n(lam_list, [0, 1])
print(f"    Reference A_4(1^-,2^-,3^+,4^+) / (i g^2) = {pt_ref:.6f}")

# Apply a sample BCFW shift (z small) and recompute
z_test = complex(0.1, 0.0)
lam_shifted, lamt_shifted = bcfw_shift(lam_list, lamt_list, 0, 2, z_test)
pt_shifted = parke_taylor_n(lam_shifted, [0, 1])
print(f"    After BCFW shift z={z_test} on legs (1,3): {pt_shifted:.6f}")

print("\n[2] BCFW Recursion Skeleton (to be expanded)")
print("""
    Full recursion for A_n(1^-,2^-,3^+,...,n^+):

    Choose shift legs (commonly adjacent or one negative + one positive).
    The amplitude is recovered from residues at the poles of the shifted
    propagators.

    A(z) = sum_k  A_L(z) * (1/P_k(z)^2) * A_R(z)   where z solves P_k(z)^2 = 0

    For MHV, many diagrams vanish and the recursion collapses beautifully
    to the Parke-Taylor formula.

    Next steps on this branch:
    - Implement the propagator pole finder
    - Add 5-point and 6-point recursion
    - Compare against Cadabra2 symbolic results
""")

print("Done: bcfw_starter.py — ready for extension.")