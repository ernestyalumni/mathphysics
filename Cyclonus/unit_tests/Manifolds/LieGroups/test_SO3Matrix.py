from Cyclonus.Manifolds.LieGroups.SO3Matrix import SO3Matrix

import numpy as np
import pytest

@pytest.fixture
def setup_tests_fixture():
    setup_values = {
        "unit vector x": np.array([1., 0., 0.]),
        "unit vector y": np.array([0., 1., 0.]),
        "unit vector z": np.array([0., 0., 1.])
    }    

    return setup_values


def test_rotation_about_x_by_45degrees(setup_tests_fixture):

    Ox45deg = SO3Matrix.rotation_about_x(np.deg2rad(45))

    x = 1. / np.sqrt(2)

    expected = np.array([[1., 0., 0.],
        [0., x, -x],
        [0., x, x]])

    assert np.allclose(Ox45deg.mat_, expected)

    unit_vector_x = setup_tests_fixture["unit vector x"]
    unit_vector_y = setup_tests_fixture["unit vector y"]
    unit_vector_z = setup_tests_fixture["unit vector z"]

    result = Ox45deg.dot(unit_vector_x)

    assert np.allclose(result, unit_vector_x)

    result = Ox45deg.dot(unit_vector_y)

    assert np.allclose(result, np.array([0., x, x]))

    result = Ox45deg.dot(unit_vector_z)

    assert np.allclose(result, np.array([0., -x, x]))
