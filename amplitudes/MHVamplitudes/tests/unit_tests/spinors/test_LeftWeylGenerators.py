from mhvamplitudes.spinors.LeftWeylGenerators import LeftWeylGenerators

def test_check_lorentz_algebra():
    S = LeftWeylGenerators()
    S.print_generators()

    assert S.check_lorentz_algebra()
    print("✓ All commutation relations satisfied")
