"""
ch42_export_latex.py
=====================
Export script for ch42_coulomb_gauge.tex
Srednicki QFT — Chapter 42: Electrodynamics in Coulomb Gauge

This file contains the LaTeX source exported from the Cadabra2/Python analysis.
Run pdflatex on the output: pdflatex ch42_coulomb_gauge.tex

Usage:
    python3 ch42_export_latex.py > ch42_coulomb_gauge_exported.tex
    # or just use ch42_coulomb_gauge.tex directly
"""

import subprocess
import os
import sys

TEX_FILE = os.path.join(os.path.dirname(__file__), "ch42_coulomb_gauge.tex")


def compile_latex():
    """Compile the tex file with pdflatex."""
    if not os.path.exists(TEX_FILE):
        print(f"ERROR: {TEX_FILE} not found", file=sys.stderr)
        return False

    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", os.path.basename(TEX_FILE)],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(TEX_FILE),
    )

    if result.returncode == 0:
        print("pdflatex: ch42_coulomb_gauge.tex compiled successfully")
        return True
    else:
        print("pdflatex ERROR:")
        print(result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
        return False


def print_tex_content():
    """Print the tex file content."""
    with open(TEX_FILE, "r") as f:
        print(f.read())


if __name__ == "__main__":
    if "--compile" in sys.argv:
        success = compile_latex()
        sys.exit(0 if success else 1)
    else:
        print_tex_content()
