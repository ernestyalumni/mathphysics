"""Tests for kinematics.spinor_helicity module."""

import numpy as np
import pytest
from mhvamplitudes.kinematics.spinor_helicity import (
    sg, square_bracket, angle_bracket,
)


class TestSignFunction:
    def test_positive(self):
        assert sg(3.14) == 1

    def test_negative(self):
        assert sg(-2.7) == -1

    def test_zero(self):
        assert sg(0.0) == 0

    def test_small_positive(self):
        assert sg(1e-20) == 1

    def test_small_negative(self):
        assert sg(-1e-20) == -1


class TestSquareBracket:
    def test_smga_convention(self):
        """[ij] = omega_i * omega_j * (tilde_z_i - tilde_z_j)."""
        omega_i, omega_j = 2.0, 3.0
        tz_i, tz_j = 1.0, 4.0
        a = np.array([omega_i, omega_i * tz_i])
        b = np.array([omega_j, omega_j * tz_j])
        expected = omega_i * omega_j * (tz_i - tz_j)
        assert abs(square_bracket(a, b) - expected) < 1e-12

    def test_antisymmetry(self):
        a = np.array([1.0, 2.0])
        b = np.array([3.0, 5.0])
        assert abs(square_bracket(a, b) + square_bracket(b, a)) < 1e-12

    def test_self_bracket_zero(self):
        a = np.array([2.0, 7.0])
        assert abs(square_bracket(a, a)) < 1e-12

    def test_linearity(self):
        """[a, b+c] = [a,b] + [a,c]."""
        a = np.array([1.0, 3.0])
        b = np.array([2.0, 5.0])
        c = np.array([4.0, 1.0])
        lhs = square_bracket(a, b + c)
        rhs = square_bracket(a, b) + square_bracket(a, c)
        assert abs(lhs - rhs) < 1e-12


class TestAngleBracket:
    def test_smga_convention(self):
        """<ij> = z_i - z_j when lam = (1, z)."""
        z_i, z_j = 2.0, 5.0
        lam_i = np.array([1.0, z_i])
        lam_j = np.array([1.0, z_j])
        expected = z_i - z_j
        assert abs(angle_bracket(lam_i, lam_j) - expected) < 1e-12

    def test_antisymmetry(self):
        a = np.array([1.0, 2.0])
        b = np.array([1.0, 5.0])
        assert abs(angle_bracket(a, b) + angle_bracket(b, a)) < 1e-12
