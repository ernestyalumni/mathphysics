"""
Integration test: Berends-Giele recursion vs SMGA closed-form formula.

This is the core Phase 2 verification: compute the stripped single-minus
amplitude A_{1...n}|_{R_1} two independent ways:

1. Berends-Giele recursion (from first principles, SMGA eq. 8 / eq. A17)
2. SMGA closed-form formula (SMGA eq. 16 / eq. 39)

and verify they agree for n = 3, 4, 5, 6, 7 with random kinematics in R_1.

For n >= 8 the recursion becomes slow (exponential in n) but should still
work. Mark those as slow tests.

Reference: arXiv:2602.12176 (Guevara, Lupsasca, Skinner, Strominger, Weil)
"""

import numpy as np
import pytest

from mhvamplitudes.amplitudes.smga import smga_formula_r1
from mhvamplitudes.amplitudes.berends_giele import compute_stripped_amplitude
from mhvamplitudes.kinematics.phase_space import make_r1_config


TOL = 1e-10


class TestBGvsSMGA:
    """Verify Berends-Giele recursion reproduces the SMGA closed-form formula."""

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_bg_matches_smga(self, n):
        """BG recursion == SMGA formula for multiple random configs."""
        rng = np.random.default_rng(42 + n)
        num_trials = 10

        for trial in range(num_trials):
            omegas, tilde_zs = make_r1_config(n, rng)

            A_smga = smga_formula_r1(n, omegas, tilde_zs)
            A_bg = compute_stripped_amplitude(
                n, omegas, tilde_zs, minus_particle=0
            )

            assert abs(A_bg - A_smga) < TOL, (
                f"n={n}, trial={trial}: BG={A_bg}, SMGA={A_smga}, "
                f"diff={abs(A_bg - A_smga)}"
            )

    @pytest.mark.parametrize("n", [7])
    def test_bg_matches_smga_n7(self, n):
        """n=7: fewer trials since recursion is slower."""
        rng = np.random.default_rng(77)
        num_trials = 5

        for trial in range(num_trials):
            omegas, tilde_zs = make_r1_config(n, rng)

            A_smga = smga_formula_r1(n, omegas, tilde_zs)
            A_bg = compute_stripped_amplitude(
                n, omegas, tilde_zs, minus_particle=0
            )

            assert abs(A_bg - A_smga) < TOL, (
                f"n={n}, trial={trial}: BG={A_bg}, SMGA={A_smga}"
            )

    @pytest.mark.slow
    @pytest.mark.parametrize("n", [8, 9, 10])
    def test_bg_matches_smga_large_n(self, n):
        """n=8,9,10: slow due to exponential growth of partitions."""
        rng = np.random.default_rng(100 + n)
        num_trials = 3

        for trial in range(num_trials):
            omegas, tilde_zs = make_r1_config(n, rng)

            A_smga = smga_formula_r1(n, omegas, tilde_zs)
            A_bg = compute_stripped_amplitude(
                n, omegas, tilde_zs, minus_particle=0
            )

            assert abs(A_bg - A_smga) < TOL, (
                f"n={n}, trial={trial}: BG={A_bg}, SMGA={A_smga}"
            )


class TestBGConsistencyChecks:
    """Verify properties of BG-computed amplitudes independently of SMGA."""

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_amplitude_is_discrete(self, n):
        """A * 2^{n-2} should be an integer."""
        rng = np.random.default_rng(500 + n)
        for _ in range(10):
            omegas, tilde_zs = make_r1_config(n, rng)
            A = compute_stripped_amplitude(
                n, omegas, tilde_zs, minus_particle=0
            )
            val = A * 2 ** (n - 2)
            assert abs(val - round(val)) < TOL, \
                f"A*2^(n-2)={val} not integer for n={n}"

    @pytest.mark.parametrize("n", [4, 5, 6])
    def test_soft_theorem_bg(self, n):
        """
        Soft theorem verified: BG-computed A_n matches soft prediction.

        A_n = (1/2)(sg_{n-1,n} + sg_{n,1}) * A_{n-1}

        Note: A_{n-1} is computed via SMGA formula (not BG) because
        the sub-config omegas[:n-1] doesn't satisfy momentum conservation,
        and the BG recursion's collapse to the SMGA formula requires it.
        The BG == SMGA match is verified separately in TestBGvsSMGA.
        """
        from mhvamplitudes.amplitudes.smga import sg_ij
        rng = np.random.default_rng(600 + n)
        for _ in range(5):
            omegas, tilde_zs = make_r1_config(n, rng)

            A_n = compute_stripped_amplitude(
                n, omegas, tilde_zs, minus_particle=0
            )
            A_nm1 = smga_formula_r1(
                n - 1, omegas[:n - 1], tilde_zs[:n - 1]
            )

            last = n - 1
            sg_prev = sg_ij(omegas[last - 1], tilde_zs[last - 1],
                            omegas[last], tilde_zs[last])
            sg_wrap = sg_ij(omegas[last], tilde_zs[last],
                            omegas[0], tilde_zs[0])
            soft_factor = 0.5 * (sg_prev + sg_wrap)
            predicted = soft_factor * A_nm1

            assert abs(A_n - predicted) < TOL, \
                f"Soft theorem failed: A_n={A_n}, predicted={predicted}"

    def test_n3_equals_sg12(self):
        """A_{123} = sg_{12} (SMGA eq. 9)."""
        from mhvamplitudes.kinematics.spinor_helicity import sg, square_bracket
        from mhvamplitudes.kinematics.phase_space import make_tilde_lambdas

        rng = np.random.default_rng(999)
        for _ in range(10):
            omegas, tilde_zs = make_r1_config(3, rng)
            A = compute_stripped_amplitude(
                3, omegas, tilde_zs, minus_particle=0
            )
            tl = make_tilde_lambdas(omegas, tilde_zs)
            expected = sg(square_bracket(tl[0], tl[1]))
            assert abs(A - expected) < TOL, \
                f"A={A}, sg_12={expected}"
