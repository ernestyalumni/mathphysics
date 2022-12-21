#!/usr/bin/env sage

"""
@url https://doc.sagemath.org/html/en/tutorial/programming.html
@details See Creating Compiled Code and Standalone Python/Sage Scripts for Sage
Math documentation.

cf. https://stackoverflow.com/questions/70622228/how-to-remove-this-error-in-python-no-module-sage-all
To run this, you'll have to use Sage's Python

$ ~/TQFT/sage/sage --python SageManifoldTutorial.py 45
"""
from pathlib import Path

import sys

# Not needed. Instead we are running the sage binary.
#sys.path.append("/home/topolo/TQFT/sage")

# This works but we'll try to import specifically.
# from sage.all import *
from sage.all import factor, sage_eval

# To look up where modules might possibly be, go I look at the source code on
# github:
# https://github.com/sagemath/sage/blob/develop/src/sage/manifolds/differentiable/manifold.py

from sage.manifolds.manifold import *

if __name__ == "__main__":

    print(sys.path)

    M = Manifold(3, 'M', latex_name=r'\mathcal{M}', start_index=1)

    if len(sys.argv) != 2:
        print("Usage: %s <n>" % sys.argv[0])
        print("Outputs the prime factorization of n.")

        print(factor(42))
        sys.exit(1)

    print(sys.argv[1])
    print(factor(sage_eval(sys.argv[1])))


