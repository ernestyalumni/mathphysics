import numpy as np

pauli_matrices = {
    1: np.array([[0, 1], [1, 0]], dtype=complex),
    2: np.array([[0, -1j], [1j, 0]], dtype=complex),
    3: np.array([[1, 0], [0, -1]], dtype=complex),
}

I2 = np.eye(2, dtype=complex)