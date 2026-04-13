"""
kinematics.spinor_helicity — spinor-helicity bracket operations.

Implements the basic spinor contractions in the SMGA frame (arXiv:2602.12176).

Frame conventions (SMGA eq. 2):
    |i> = lambda_i = (1, z_i)
    |i] = tilde_lambda_i = omega_i * (1, tilde_z_i)

Bracket definitions (SMGA eq. 4):
    <ij> = z_{ij} = z_i - z_j
    [ij] = omega_i * omega_j * tilde_z_{ij}

where z_{ij} = z_i - z_j and tilde_z_{ij} = tilde_z_i - tilde_z_j.

The sign function sg(x) = +1 if x > 0, -1 if x < 0, 0 if x = 0.
"""

import numpy as np


def sg(x: float) -> int:
    """Sign function: sg(x) = +1 if x > 0, -1 if x < 0, 0 if x = 0."""
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def square_bracket(a: np.ndarray, b: np.ndarray) -> float:
    """
    Compute the square bracket [a, b] for two 2-component tilde-lambda vectors.

    SMGA convention: [a, b] = a[1]*b[0] - a[0]*b[1]

    This gives [ij] = omega_i * omega_j * (tilde_z_i - tilde_z_j) when
    a = omega_i * (1, tilde_z_i) and b = omega_j * (1, tilde_z_j).
    """
    return a[1] * b[0] - a[0] * b[1]


def angle_bracket(lam_i: np.ndarray, lam_j: np.ndarray) -> float:
    """
    Compute the angle bracket <ij> for two 2-component lambda vectors.

    SMGA convention: <ij> = lam_i[1]*lam_j[0] - lam_i[0]*lam_j[1] = z_i - z_j
    when lam_i = (1, z_i) and lam_j = (1, z_j).

    Same contraction convention as square_bracket (epsilon_{12} = -1).
    """
    return lam_i[1] * lam_j[0] - lam_i[0] * lam_j[1]


def construct_spinors(p: np.ndarray):
    """
    Construct spinor-helicity variables from a 4-momentum p = (p0, px, py, pz).

    Uses mostly-plus metric g = diag(-1,+1,+1,+1).
    Returns (lambda, tilde_lambda) as 2-component complex arrays.

    For generic (non-half-collinear) kinematics. Uses the standard
    decomposition via p0 + pz.
    """
    p0, px, py, pz = p
    r = complex(p0 + pz)
    if abs(r) < 1e-10:
        r = r + 1e-8
    sqrt_r = np.sqrt(r)
    lam = np.array([sqrt_r, (px + 1j * py) / sqrt_r], dtype=complex)
    lam_tilde = np.array([sqrt_r, (px - 1j * py) / sqrt_r], dtype=complex)
    return lam, lam_tilde
