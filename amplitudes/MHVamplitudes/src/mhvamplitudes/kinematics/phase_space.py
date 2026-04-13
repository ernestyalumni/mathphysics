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

Convention note: SMGA half-collinear frame vs Srednicki
=======================================================

This module generates phase-space points in the SMGA half-collinear frame,
which is specific to (2,2) Klein signature. This frame has NO direct analogue
in Srednicki's Lorentzian (1,3) framework.

In the half-collinear limit (SMGA sec. 2):
    - All holomorphic spinors become proportional: |i> ~ |1> (z_i -> z for all i)
    - The angle brackets <ij> = z_i - z_j -> 0
    - Momenta are parametrized entirely by (omega_i, tilde_z_i)
    - The 4-momentum is p_i^{mu} ~ omega_i * q^{mu} + O(epsilon) where
      q is the common collinear direction

In Srednicki's Lorentzian framework:
    - Collinear limits exist but single-minus amplitudes A(1^-,2^+,...,n^+) = 0
      identically (by helicity selection rules / holomorphy of MHV amplitudes)
    - Phase-space generators would use explicit 4-momenta satisfying
      p^2 = 0 and momentum conservation sum p_i = 0
    - Spinors are complex: lambda_i and tilde_lambda_i are related by
      complex conjugation for real momenta

The omega_i parameters here are NOT the same as energies in Srednicki. They
are components of the tilde-lambda spinors along a reference direction. The
sign of omega_i determines the "region" (R_1, R_2, ...) and hence which
helicity configurations give nonzero amplitudes in (2,2) signature. In
Lorentzian signature, all omega_i would be positive for physical momenta.
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
