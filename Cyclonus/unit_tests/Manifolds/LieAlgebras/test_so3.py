from Cyclonus.Manifolds.LieAlgebras.so3 import so3

import numpy as np
import pytest


@pytest.fixture
def setup_tests_fixture():
    setup_values = {
        "input omega": [1, 2, 3],
        "input omegas": np.array([[1, 2, 3], [4, 5, 6]]),
        "expected omega hat 1": np.array([
            [0., -3., 2.],
            [3., 0., -1.],
            [-2., 1., 0.]]),
        "expected omega hat 2": np.array([
            [0., -6., 5.],
            [6., 0., -4.],
            [-5., 4., 0.]])
    }    

    return setup_values


def test_hat_map(setup_tests_fixture):
    omega = setup_tests_fixture["input omega"]
    omegas = setup_tests_fixture["input omegas"]
    expected1 = setup_tests_fixture["expected omega hat 1"]
    expected2 = setup_tests_fixture["expected omega hat 2"]

    omega_hat = so3.hat_map(omega)

    omega_hats = so3.hat_map(omegas)

    assert np.array_equal(omega_hat, expected1)
    assert np.array_equal(omega_hats[0], expected1)
    assert np.array_equal(omega_hats[1], expected2)


def test_inverse_hat_map(setup_tests_fixture):
    omega = setup_tests_fixture["input omega"]
    omegas = setup_tests_fixture["input omegas"]
    omega_hat_1 = setup_tests_fixture["expected omega hat 1"]
    omega_hat_2 = setup_tests_fixture["expected omega hat 2"]

    omega_1 = so3.inverse_hat_map(omega_hat_1)

    resulting_omegas = so3.inverse_hat_map(np.array([omega_hat_1, omega_hat_2]))

    assert np.array_equal(omega_1, omega)
    assert np.array_equal(resulting_omegas[0], omegas[0])
    assert np.array_equal(resulting_omegas[1], omegas[1])
