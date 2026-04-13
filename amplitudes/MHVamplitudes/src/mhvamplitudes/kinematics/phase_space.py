"""
kinematics.phase_space — momentum configuration generators.

Generates kinematic configurations for scattering amplitude calculations,
particularly in the half-collinear regime of SMGA (arXiv:2602.12176).

In the SMGA frame:
    tilde_lambda_i = omega_i * (1, tilde_z_i)

Region R_1 (SMGA eq. 25):
    omega_1 < 0,  omega_a > 0  for a >= 2

Momentum conservation (half-collinear):
    sum_i omega_i = 0
    sum_i omega_i * tilde_z_i = 0
"""

import numpy as np


def make_tilde_lambdas(omegas, tilde_zs):
    """
    Construct tilde-lambda 2-vectors from omega and tilde_z parameters.

    Parameters
    ----------
    omegas : list of float
        Energy-like parameters for each particle.
    tilde_zs : list of float
        Holomorphic coordinates for each particle.

    Returns
    -------
    list of np.ndarray
        Each element is a 2-component array omega_i * (1, tilde_z_i).
    """
    return [np.array([w, w * z]) for w, z in zip(omegas, tilde_zs)]


def make_r1_config(n, rng, conserve_momentum=True):
    """
    Generate a random kinematic configuration in region R_1.

    R_1: omega_0 < 0, omega_a > 0 for a >= 1 (0-indexed).
    tilde_z_i: random distinct real values.

    Parameters
    ----------
    n : int
        Number of particles (>= 3).
    rng : numpy.random.Generator
        Random number generator.
    conserve_momentum : bool
        If True, enforce sum omega_i = 0 and sum omega_i * tilde_z_i = 0.

    Returns
    -------
    omegas : list of float
        Energy parameters. omegas[0] < 0, rest > 0.
    tilde_zs : list of float
        Holomorphic coordinates. All distinct.
    """
    assert n >= 3

    omegas_rest = [rng.uniform(0.5, 3.0) for _ in range(n - 1)]
    tilde_zs_rest = list(rng.uniform(-5.0, 5.0, size=n - 1))

    # Ensure tilde_z values among particles 2..n are distinct
    for i in range(1, n - 1):
        while any(abs(tilde_zs_rest[i] - tilde_zs_rest[j]) < 1e-6
                  for j in range(i)):
            tilde_zs_rest[i] += rng.uniform(-0.1, 0.1)

    if conserve_momentum:
        omega_0 = -sum(omegas_rest)
        sum_wz = sum(o * z for o, z in zip(omegas_rest, tilde_zs_rest))
        tilde_z_0 = sum_wz / (-omega_0)
    else:
        omega_0 = -rng.uniform(0.5, 3.0)
        tilde_z_0 = rng.uniform(-5.0, 5.0)

    omegas = [omega_0] + omegas_rest
    tilde_zs = [tilde_z_0] + tilde_zs_rest

    # Ensure tilde_z_0 is distinct from the rest
    for i in range(1, n):
        while abs(tilde_zs[0] - tilde_zs[i]) < 1e-6:
            # Perturb one of the positive-omega particles and recompute
            tilde_zs_rest[i - 1] += rng.uniform(-0.1, 0.1)
            tilde_zs[i] = tilde_zs_rest[i - 1]
            if conserve_momentum:
                sum_wz = sum(o * z for o, z
                             in zip(omegas_rest, tilde_zs_rest))
                tilde_zs[0] = sum_wz / (-omega_0)

    return omegas, tilde_zs
