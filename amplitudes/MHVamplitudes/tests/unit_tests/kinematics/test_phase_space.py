"""Tests for kinematics.phase_space module."""

import numpy as np
import pytest
from mhvamplitudes.kinematics.phase_space import make_r1_config, make_tilde_lambdas


class TestMakeTildeLambdas:
    def test_shape(self):
        tls = make_tilde_lambdas([1.0, 2.0, 3.0], [0.5, 1.0, 1.5])
        assert len(tls) == 3
        for tl in tls:
            assert tl.shape == (2,)

    def test_values(self):
        """tilde_lambda_i = omega_i * (1, tilde_z_i)."""
        tls = make_tilde_lambdas([2.0], [3.0])
        np.testing.assert_allclose(tls[0], [2.0, 6.0])


class TestMakeR1Config:
    def test_r1_signs(self):
        """omega_0 < 0, rest > 0."""
        rng = np.random.default_rng(42)
        for n in [3, 4, 5, 6, 7]:
            omegas, tilde_zs = make_r1_config(n, rng)
            assert omegas[0] < 0, f"omega_0 should be < 0 for n={n}"
            for a in range(1, n):
                assert omegas[a] > 0, f"omega_{a} should be > 0 for n={n}"

    def test_momentum_conservation(self):
        """sum omega_i = 0 and sum omega_i * tilde_z_i = 0."""
        rng = np.random.default_rng(123)
        for n in [3, 5, 8]:
            omegas, tilde_zs = make_r1_config(n, rng)
            assert abs(sum(omegas)) < 1e-10, "sum omega_i != 0"
            sum_wz = sum(o * z for o, z in zip(omegas, tilde_zs))
            assert abs(sum_wz) < 1e-10, "sum omega_i*tilde_z_i != 0"

    def test_distinct_tilde_z(self):
        """All tilde_z values should be distinct."""
        rng = np.random.default_rng(999)
        omegas, tilde_zs = make_r1_config(6, rng)
        for i in range(6):
            for j in range(i + 1, 6):
                assert abs(tilde_zs[i] - tilde_zs[j]) > 1e-8

    def test_n3_minimum(self):
        rng = np.random.default_rng(0)
        omegas, tilde_zs = make_r1_config(3, rng)
        assert len(omegas) == 3
        assert len(tilde_zs) == 3
