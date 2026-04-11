from mhvamplitudes.spinors.invariants import (
    EpsilonTensor,
)

def test_verify():
    eps = EpsilonTensor()
    assert eps.verify()
