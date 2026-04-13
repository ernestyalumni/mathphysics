"""Tests for amplitudes.berends_giele module."""

import numpy as np
import pytest
from mhvamplitudes.amplitudes.berends_giele import (
    ordered_partitions,
    BerendsGieleRecursion,
    compute_stripped_amplitude,
)
from mhvamplitudes.kinematics.spinor_helicity import sg, square_bracket
from mhvamplitudes.kinematics.phase_space import make_tilde_lambdas


class TestOrderedPartitions:
    def test_single_element(self):
        parts = list(ordered_partitions((1,)))
        assert parts == [((1,),)]

    def test_two_elements_all(self):
        parts = list(ordered_partitions((1, 2)))
        # k=1: ((1,2),)
        # k=2: ((1,),(2,))
        assert ((1, 2),) in parts
        assert ((1,), (2,)) in parts
        assert len(parts) == 2

    def test_three_elements_all(self):
        parts = list(ordered_partitions((1, 2, 3)))
        # k=1: ((1,2,3),)
        # k=2: ((1,),(2,3)), ((1,2),(3,))
        # k=3: ((1,),(2,),(3,))
        assert len(parts) == 4

    def test_count_formula(self):
        """Number of ordered partitions of n elements = 2^{n-1}."""
        for n in range(1, 8):
            seq = tuple(range(n))
            count = len(list(ordered_partitions(seq)))
            assert count == 2 ** (n - 1), f"n={n}: got {count}"

    def test_min_parts_filter(self):
        parts = list(ordered_partitions((1, 2, 3), min_parts=3))
        assert len(parts) == 1
        assert parts[0] == ((1,), (2,), (3,))

    def test_min_parts_2(self):
        parts = list(ordered_partitions((1, 2, 3), min_parts=2))
        # k=2 and k=3 only
        assert len(parts) == 3

    def test_max_parts(self):
        parts = list(ordered_partitions((1, 2, 3, 4), max_parts=2))
        # k=1 and k=2 only: 1 + C(3,1) = 4
        assert len(parts) == 4

    def test_partition_covers_sequence(self):
        """Each partition should reconstruct the original sequence."""
        seq = (10, 20, 30, 40)
        for partition in ordered_partitions(seq):
            reconstructed = sum(partition, ())
            assert reconstructed == seq


class TestVertex:
    def _make_bg(self, omegas, tilde_zs):
        tl = make_tilde_lambdas(omegas, tilde_zs)
        tl_dict = {i: tl[i] for i in range(len(tl))}
        return BerendsGieleRecursion(tl_dict), tl

    def test_V_single_vector(self):
        """V of a single vector = 1."""
        bg, tl = self._make_bg([1.0], [2.0])
        assert bg.vertex_V([tl[0]]) == 1.0

    def test_V_two_vectors_zero(self):
        """V of two vectors = 0 (Theta(-1) = 0)."""
        bg, tl = self._make_bg([1.0, 2.0], [1.0, 3.0])
        assert bg.vertex_V([tl[0], tl[1]]) == 0.0

    def test_Vbar_two_vectors(self):
        """V-bar of two vectors = sg_{12} (Theta(+1) = 1)."""
        bg, tl = self._make_bg([1.0, 2.0], [1.0, 3.0])
        vbar = bg.vertex_Vbar([tl[0], tl[1]])
        expected_sg = sg(square_bracket(tl[0], tl[1]))
        assert vbar == expected_sg

    def test_pt_hat_two_vectors(self):
        """PT-hat of two vectors = -sg_{12}."""
        bg, tl = self._make_bg([1.0, 2.0], [1.0, 3.0])
        pt = bg.pt_hat([tl[0], tl[1]])
        expected_sg = sg(square_bracket(tl[0], tl[1]))
        assert pt == -expected_sg

    def test_V_vanishes_in_r1(self):
        """In R_1, V_{lam_2...lam_n} = 0 (SMGA eq. 33)."""
        from mhvamplitudes.kinematics.phase_space import make_r1_config
        rng = np.random.default_rng(42)
        for n in [4, 5, 6]:
            omegas, tilde_zs = make_r1_config(n, rng)
            bg, _ = self._make_bg(omegas, tilde_zs)
            tl = make_tilde_lambdas(omegas, tilde_zs)
            # V applied to particles 2..n (indices 1..n-1), all positive omega
            plus_vectors = tl[1:]
            assert bg.vertex_V(plus_vectors) == 0.0, \
                f"V should vanish in R_1 for n={n}"


class TestPreamplitude:
    def _make_bg(self, omegas, tilde_zs):
        tl = make_tilde_lambdas(omegas, tilde_zs)
        tl_dict = {i: tl[i] for i in range(len(tl))}
        return BerendsGieleRecursion(tl_dict)

    def test_singleton(self):
        bg = self._make_bg([1.0], [0.0])
        assert bg.preamplitude((0,)) == 1.0

    def test_pair(self):
        bg = self._make_bg([1.0, 2.0], [0.0, 1.0])
        assert bg.preamplitude((0, 1)) == 0.0

    def test_triple_in_r1(self):
        """In R_1, A-bar_S = 0 for |S| >= 2 subsets of {2,...,n}."""
        from mhvamplitudes.kinematics.phase_space import make_r1_config
        rng = np.random.default_rng(42)
        omegas, tilde_zs = make_r1_config(5, rng)
        bg = self._make_bg(omegas, tilde_zs)
        # Preamplitudes for subsets of the positive-omega particles
        assert bg.preamplitude((1, 2, 3)) == 0.0
        assert bg.preamplitude((2, 3, 4)) == 0.0


class TestStrippedAmplitude:
    def _make_bg(self, omegas, tilde_zs):
        tl = make_tilde_lambdas(omegas, tilde_zs)
        tl_dict = {i: tl[i] for i in range(len(tl))}
        return BerendsGieleRecursion(tl_dict)

    def test_n3_specific(self):
        """A_{123} = sg_{12} for general kinematics (SMGA eq. 9)."""
        # Use momentum-conserving config
        omegas = [-3.0, 1.0, 2.0]
        tz_1 = 0.0
        tz_2 = 1.0
        tz_3 = -(omegas[0] * tz_1 + omegas[1] * tz_2) / omegas[2]
        tilde_zs = [tz_1, tz_2, tz_3]

        bg = self._make_bg(omegas, tilde_zs)

        # A_{1,2,3} with particle 0 minus => A_{1,2,0} via cyclicity
        A = bg.stripped_amplitude((1, 2, 0))
        tl = make_tilde_lambdas(omegas, tilde_zs)
        expected = sg(square_bracket(tl[0], tl[1]))
        assert abs(A - expected) < 1e-10, f"A={A}, expected sg_12={expected}"

    def test_n3_via_convenience(self):
        """Test compute_stripped_amplitude convenience function."""
        omegas = [-3.0, 1.0, 2.0]
        tz_1 = 0.0
        tz_2 = 1.0
        tz_3 = -(omegas[0] * tz_1 + omegas[1] * tz_2) / omegas[2]
        tilde_zs = [tz_1, tz_2, tz_3]

        A = compute_stripped_amplitude(3, omegas, tilde_zs, minus_particle=0)
        tl = make_tilde_lambdas(omegas, tilde_zs)
        expected = sg(square_bracket(tl[0], tl[1]))
        assert abs(A - expected) < 1e-10
