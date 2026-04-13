"""
amplitudes.smga — SMGA closed-form formula for single-minus amplitudes.

Implements the closed-form expression from arXiv:2602.12176 (Guevara,
Lupsasca, Skinner, Strominger, Weil) for stripped single-minus gluon
tree amplitudes in the kinematic region R_1.

Main formula (SMGA eq. 16 / eq. 39):

    A_{1...n}|_{R_1} = (1/2^{n-2}) * prod_{m=2}^{n-1} (sg_{m,m+1} + sg_{1,2...m})

where:
    sg_{m,m+1} = sg([tilde_lambda_m, tilde_lambda_{m+1}])
    sg_{1,2...m} = sg([tilde_lambda_1, tilde_lambda_2 + ... + tilde_lambda_m])

Particles are 1-indexed in the physics convention but 0-indexed in code:
    paper particle 1 -> code index 0, etc.
"""

from ..kinematics.spinor_helicity import sg


def sg_ij(omega_i, tilde_z_i, omega_j, tilde_z_j):
    """
    Sign of the square bracket [ij].

    sg_{ij} = sg(omega_i * omega_j * (tilde_z_i - tilde_z_j))
    """
    return sg(omega_i * omega_j * (tilde_z_i - tilde_z_j))


def sg_i_set(omega_i, tilde_z_i, omegas_S, tilde_zs_S):
    """
    Sign of [tilde_lambda_i, sum_{j in S} tilde_lambda_j].

    In the SMGA frame:
        [lam_i, sum_S lam_j] = omega_i * (Omega_S * tilde_z_i - Z_S)

    where Omega_S = sum omega_j and Z_S = sum omega_j * tilde_z_j.
    """
    Omega_S = sum(omegas_S)
    Z_S = sum(o * z for o, z in zip(omegas_S, tilde_zs_S))
    val = omega_i * (Omega_S * tilde_z_i - Z_S)
    return sg(val)


def smga_formula_r1(n, omegas, tilde_zs):
    """
    Compute A_{1...n}|_{R_1} using the SMGA closed-form formula (eq. 16).

    A_{1...n}|_{R_1} = (1/2^{n-2}) * prod_{m=2}^{n-1} (sg_{m,m+1} + sg_{1,2...m})

    Parameters
    ----------
    n : int
        Number of particles (>= 3).
    omegas : list of float
        Energy parameters. 0-indexed: omegas[0] < 0 (particle 1 in paper).
    tilde_zs : list of float
        Holomorphic coordinates.

    Returns
    -------
    float
        The stripped amplitude value (piecewise constant: 0, +/-1).
    """
    assert n >= 3
    assert len(omegas) == n and len(tilde_zs) == n

    result = 1.0
    for m in range(2, n):
        # Paper m runs 2..n-1 (1-indexed). Code: i1=m-1, i2=m (0-indexed).
        i1 = m - 1
        i2 = m
        sg_adjacent = sg_ij(omegas[i1], tilde_zs[i1],
                            omegas[i2], tilde_zs[i2])

        # sg_{1,2...m}: particle 1 (index 0) vs sum of particles 2..m
        # (indices 1..m-1)
        sg_1_partial = sg_i_set(
            omegas[0], tilde_zs[0],
            omegas[1:m], tilde_zs[1:m]
        )

        result *= (sg_adjacent + sg_1_partial)

    return result / (2 ** (n - 2))
