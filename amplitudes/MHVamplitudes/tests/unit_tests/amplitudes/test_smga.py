"""Tests for amplitudes.smga module."""

import numpy as np
import pytest
from mhvamplitudes.amplitudes.smga import sg_ij, sg_i_set, smga_formula_r1
from mhvamplitudes.kinematics.phase_space import make_r1_config


class TestSgIj:
    def test_positive_omegas(self):
        """For positive omegas, sg_{ij} = sg(tilde_z_i - tilde_z_j)."""
        assert sg_ij(1.0, 3.0, 2.0, 1.0) == 1   # 1*2*(3-1)=4 > 0
        assert sg_ij(1.0, 1.0, 2.0, 3.0) == -1   # 1*2*(1-3)=-4 < 0

    def test_mixed_omega_signs(self):
        """In R_1, omega_1 < 0, which flips the sign."""
        # omega_1=-1, omega_2=2, tilde_z_1=3, tilde_z_2=1
        # sg(-1 * 2 * (3-1)) = sg(-4) = -1
        assert sg_ij(-1.0, 3.0, 2.0, 1.0) == -1

    def test_antisymmetry_of_sign(self):
        """sg_{ij} = -sg_{ji} (when omegas have same sign)."""
        assert sg_ij(1.0, 2.0, 3.0, 5.0) == -sg_ij(3.0, 5.0, 1.0, 2.0)


class TestSgISet:
    def test_single_element_set(self):
        """sg_{i,{j}} should equal sg_{ij}."""
        w_i, z_i = -1.0, 0.0
        w_j, z_j = 2.0, 1.0
        assert sg_i_set(w_i, z_i, [w_j], [z_j]) == sg_ij(w_i, z_i, w_j, z_j)

    def test_two_element_set(self):
        """sg_{1,{2,3}} with known values."""
        # omega_1=-3, tz_1=0, omega_2=1, tz_2=1, omega_3=2, tz_3=-0.5
        # Omega_S = 1+2 = 3, Z_S = 1*1 + 2*(-0.5) = 0
        # val = -3 * (3*0 - 0) = 0
        assert sg_i_set(-3.0, 0.0, [1.0, 2.0], [1.0, -0.5]) == 0


class TestSmgaFormulaR1:
    def test_n3_specific(self):
        """A_{123}|_{R1} with specific config. From script 08."""
        # omegas=[-1,1,2], tilde_zs=[3,1,2]
        # sg_{23} = sg(1*2*(1-2)) = sg(-2) = -1
        # sg_{1,{2}} = sg(-1*(1*3 - 1)) = sg(-2) = -1
        # A = (1/2)(-1 + (-1)) = -1
        A = smga_formula_r1(3, [-1.0, 1.0, 2.0], [3.0, 1.0, 2.0])
        assert abs(A - (-1.0)) < 1e-12

    def test_n3_zero_case(self):
        """A_{123}|_{R1} can be zero."""
        # omegas=[-1,1,2], tilde_zs=[0,1,2]
        # sg_{23} = sg(1*2*(1-2)) = -1
        # sg_{1,{2}} = sg(-1*(1*0 - 1)) = sg(1) = +1
        # A = (1/2)(-1 + 1) = 0
        A = smga_formula_r1(3, [-1.0, 1.0, 2.0], [0.0, 1.0, 2.0])
        assert abs(A) < 1e-12

    def test_values_are_discrete(self):
        """A_{1...n}|_{R1} * 2^{n-2} should be an integer."""
        rng = np.random.default_rng(42)
        for n in [3, 4, 5, 6]:
            for _ in range(20):
                omegas, tilde_zs = make_r1_config(n, rng)
                A = smga_formula_r1(n, omegas, tilde_zs)
                val = A * 2 ** (n - 2)
                assert abs(val - round(val)) < 1e-10, \
                    f"A*2^(n-2) = {val} not integer for n={n}"

    def test_soft_theorem(self):
        """
        Weinberg soft theorem (SMGA eq. 8):
        A_{1...n} = (1/2)(sg_{n-1,n} + sg_{n,1}) * A_{1...n-1}

        Since the amplitude is piecewise constant, this is exact
        (not a limit) within each chamber.
        """
        rng = np.random.default_rng(200)
        for n in [4, 5, 6]:
            for _ in range(5):
                omegas, tilde_zs = make_r1_config(n, rng)
                A_n = smga_formula_r1(n, omegas, tilde_zs)
                A_nm1 = smga_formula_r1(n - 1, omegas[:n - 1],
                                        tilde_zs[:n - 1])

                last = n - 1  # 0-indexed
                sg_prev_last = sg_ij(omegas[last - 1], tilde_zs[last - 1],
                                     omegas[last], tilde_zs[last])
                sg_last_first = sg_ij(omegas[last], tilde_zs[last],
                                      omegas[0], tilde_zs[0])
                soft_factor = 0.5 * (sg_prev_last + sg_last_first)

                predicted = soft_factor * A_nm1
                assert abs(A_n - predicted) < 1e-10, \
                    f"Soft theorem failed for n={n}: A_n={A_n}, predicted={predicted}"
