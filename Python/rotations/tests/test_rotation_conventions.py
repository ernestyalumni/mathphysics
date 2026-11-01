"""Independent acceptance checks for the approved rotation conventions."""
import math
from pathlib import Path
import subprocess
import sys
import unittest

import rotations as r


class RotationConventions(unittest.TestCase):
    def test_zyx_product_order_against_explicit_formula(self):
        h, p, b = 0.4, 0.7, -0.3
        ch, sh, cp, sp, cb, sb = math.cos(h), math.sin(h), math.cos(p), math.sin(p), math.cos(b), math.sin(b)
        expected = ((ch*cp, ch*sp*sb-sh*cb, ch*sp*cb+sh*sb),
                    (sh*cp, sh*sp*sb+ch*cb, sh*sp*cb-ch*sb),
                    (-sp, cp*sb, cp*cb))
        matrix = r.from_euler('zyx', h, p, b)
        for row, target in zip(matrix, expected):
            for value, reference in zip(row, target):
                self.assertAlmostEqual(value, reference, places=14)
        angles = r.to_euler(expected, 'zyx')
        self.assertFalse(angles.singular)
        for value, reference in zip(angles.values, (h, p, b)):
            self.assertAlmostEqual(value, reference, places=14)

    def test_quaternion_sign_and_independent_quarter_turn(self):
        q = r.to_quaternion(((0., -1., 0.), (1., 0., 0.), (0., 0., 1.)))
        for value, reference in zip(q, (math.sqrt(0.5), 0., 0., math.sqrt(0.5))):
            self.assertAlmostEqual(value, reference, places=14)
        # Negative scalar representatives and all largest-diagonal branches.
        for axis in ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.), (2/3, -1/3, 2/3)):
            for angle in (-5., -3.5, -math.pi, 0., math.pi, 3.5, 5.):
                self.assertGreaterEqual(r.to_quaternion(r.from_axis_angle(axis, angle))[0], 0.)

    def test_identity_axis_and_canonical_angle_ranges(self):
        self.assertEqual(r.to_axis_angle(r.I3), ((1., 0., 0.), 0.))
        for convention in ('zxz', 'zyz', 'zyx'):
            for middle in (-8., -math.pi, -1.8, 0., 1.8, math.pi, 8.):
                matrix = r.from_euler(convention, .4, middle, -.3)
                angles = r.to_euler(matrix, convention)
                lo, hi = (-math.pi/2, math.pi/2) if convention == 'zyx' else (0., math.pi)
                self.assertTrue(lo <= angles.middle <= hi)
                self.assertTrue(0. <= r.to_axis_angle(matrix)[1] <= math.pi)

    def test_ck_components_and_goldstein_transpose(self):
        q = (.5, -.5, .5, -.5)
        ck = (.5-.5j, .5-.5j, -.5-.5j, .5+.5j)
        self.assertEqual(r.quaternion_to_ck(q), ck)
        self.assertEqual(r.ck_to_quaternion(ck), q)
        # Hamilton q maps ex -> -ey, ey -> -ez, ez -> ex; Goldstein uses the transpose.
        self.assertEqual(r.ck_to_goldstein(ck), ((0., -1., 0.), (0., 0., -1.), (1., 0., 0.)))
        for value in (complex(float('nan'), 0), complex(0, float('inf'))):
            with self.assertRaises(ValueError):
                r.ck_to_quaternion((value, 0j, 0j, value.conjugate()))

    def test_no_silent_projection_or_input_mutation(self):
        accepted = [[1.+1e-12, 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        original = [row[:] for row in accepted]
        self.assertEqual(r.validate_matrix(accepted), tuple(map(tuple, original)))
        self.assertEqual(accepted, original)
        bad = ((1.01, 0., 0.), (0., 1., 0.), (0., 0., 1.))
        for convert in (r.validate_matrix, r.to_quaternion, r.to_axis_angle,
                        lambda a: r.to_euler(a, 'zyx')):
            with self.assertRaises(ValueError):
                convert(bad)
        with self.assertRaises(ValueError):
            r.from_axis_angle((2., 0., 0.), .4)
        with self.assertRaises(ValueError):
            r.from_quaternion((2., 0., 0., 0.))

    def test_runtime_works_without_site_packages(self):
        root = str(Path(r.__file__).resolve().parent)
        script = """
import sys
sys.path.insert(0, sys.argv[1])
import rotations as r
m = r.from_euler('zyx', .4, .7, -.3)
r.to_euler(m, 'zyx')
r.to_axis_angle(m)
r.ck_to_goldstein(r.quaternion_to_ck(r.to_quaternion(m)))
assert not any('site-packages' in getattr(mod, '__file__', '') for mod in tuple(sys.modules.values()) if getattr(mod, '__file__', None))
"""
        subprocess.run([sys.executable, '-I', '-S', '-c', script, root], check=True, capture_output=True, text=True)
