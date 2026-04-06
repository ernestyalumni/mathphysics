#!/usr/bin/env python3
"""
Strip preamble from each chapter .tex file and save as _clean.tex.
Content between \begin{document} and \end{document} is extracted.
Preamble lines (\maketitle, \tableofcontents, \newpage immediately after
\begin{document}) are also stripped.
"""
import os
import re

CHAPTERS_DIR = os.path.join(os.path.dirname(__file__), "individual_chapters")

# Files to process (in order)
FILES = [
    "ch34_left_right_spinors.tex",
    "ch35_sigma_algebra.tex",
    "ch36_weyl_lagrangian.tex",
    "ch37_canonical_quantization.tex",
    "ch38_spinor_technology.tex",
    "ch39_canonical_quantization_II.tex",
    "ch40_spin_sums.tex",
    "ch41_gamma_technology.tex",
    "ch41_feynman_rules_dirac.tex",
    "ch42_coulomb_gauge.tex",
    "ch43_feynman_rules_majorana.tex",
    "ch50_spinor_helicity.tex",
    "ch50_massless_spinor_helicity.tex",
    "ch60_spinor_helicity.tex",
]

# Lines to strip from the body (immediately after \begin{document})
STRIP_BODY_LINES = {
    r'\maketitle',
    r'\tableofcontents',
    r'\newpage',
    r'\bibliographystyle{plain}',
}

def strip_preamble(content):
    """Extract content between \\begin{document} and \\end{document}."""
    # Find \begin{document}
    begin_match = re.search(r'\\begin\{document\}', content)
    if not begin_match:
        raise ValueError("No \\begin{document} found")

    # Find \end{document}
    end_match = re.search(r'\\end\{document\}', content)
    if not end_match:
        raise ValueError("No \\end{document} found")

    body = content[begin_match.end():end_match.start()]

    # Strip leading whitespace/newlines
    body = body.lstrip('\n')

    # Strip top-level preamble-like lines (maketitle, tableofcontents, newpage)
    lines = body.split('\n')
    # Remove leading preamble-like lines
    # We'll remove them wherever they appear as standalone lines at the top
    clean_lines = []
    stripping_top = True
    for line in lines:
        stripped = line.strip()
        if stripping_top and stripped in STRIP_BODY_LINES:
            continue  # skip these top-level lines
        elif stripping_top and stripped == '':
            continue  # skip blank lines at top
        else:
            stripping_top = False
            clean_lines.append(line)

    # Also strip \bibliographystyle from end
    while clean_lines and clean_lines[-1].strip() in STRIP_BODY_LINES | {''}:
        clean_lines.pop()

    return '\n'.join(clean_lines) + '\n'


for filename in FILES:
    src = os.path.join(CHAPTERS_DIR, filename)
    dst = os.path.join(CHAPTERS_DIR, filename.replace('.tex', '_clean.tex'))

    with open(src, 'r', encoding='utf-8') as f:
        content = f.read()

    try:
        clean = strip_preamble(content)
        with open(dst, 'w', encoding='utf-8') as f:
            f.write(clean)
        print(f"OK: {filename} -> {os.path.basename(dst)} ({len(clean.splitlines())} lines)")
    except ValueError as e:
        print(f"ERROR: {filename}: {e}")

print("Done.")
