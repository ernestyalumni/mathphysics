"""
Srednicki.spinors.representations — High-level interface to Weyl
representations.
"""

from .LeftWeylGenerators import LeftWeylGenerators


class WeylRepresentation:
    """Convenience wrapper around the numerical generators."""

    def __init__(self):
        self.generators = LeftWeylGenerators()

    def verify(self):
        print("Verifying Lorentz algebra in (2,1) representation...")
        if self.generators.check_lorentz_algebra():
            print("✓ Lorentz algebra satisfied")
        else:
            print("✗ Algebra verification failed")
        self.generators.print_generators()
        return True
