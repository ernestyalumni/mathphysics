# Srednicki QFT — Chapter-by-Chapter Cadabra2 Computations

**Goal:** Build up spinor technology chapter by chapter, verifying every key
equation numerically and exporting clean LaTeX/PDF notes — from Weyl spinors
(Ch. 34) through the full spinor-helicity and MHV amplitude machinery.

---

## Directory Structure

```
Srednicki/
├── README.md                        ← this file
│
├── ch34_left_right_spinors.py       ← computation: Ch. 34
├── ch34_export_latex.py             ← LaTeX exporter: Ch. 34
├── ch34_left_right_spinors.tex      ← generated .tex
├── ch34_left_right_spinors.pdf      ← compiled PDF
│
├── ch35_sigma_algebra.py            ← computation: Ch. 35
├── ch35_export_latex.py
├── ch35_sigma_algebra.tex / .pdf
│
├── ch36_weyl_lagrangian.py          ← Ch. 36 …
│   …
└── chNN_*.py                        ← successive chapters
```

Each chapter pair:
- `chNN_<topic>.py` — all computation, verification, cadabra2 symbolic algebra
- `chNN_export_latex.py` — pulls results from the above, generates `.tex`, compiles to `.pdf`

---

## Prerequisites

### 1. Cadabra2 Docker image

All cadabra2 code runs inside a Docker container.
The image is pre-built locally as `cadabra2-ubuntu:24.04`.

Check it exists:
```bash
docker image ls | grep cadabra
# cadabra2-ubuntu    24.04    ...
```

If missing, rebuild (from the Cadabra2 build directory):
```bash
docker build -t cadabra2-ubuntu:24.04 .
```

### 2. Python packages inside the image

The container already has `cadabra2`, `numpy`, `sympy`.
No extra installs needed for the Srednicki scripts.

### 3. LaTeX (on the host, for PDF compilation)

```bash
which pdflatex     # /usr/bin/pdflatex
which latexmk      # /usr/bin/latexmk
```

Packages used: `amsmath`, `amssymb`, `mathtools`, `booktabs`, `geometry`,
`hyperref`, `fancyhdr`, `xcolor`.  All standard in a full TeXLive install.

---

## Running a Chapter

### Step 1 — run the computation script

```bash
cd Python/Cadabra2/Srednicki

docker run --rm \
    -v "$(pwd)":/work \
    cadabra2-ubuntu:24.04 \
    python3 /work/ch34_left_right_spinors.py
```

This prints everything to stdout: numerical verifications, matrix values,
cadabra2 symbolic expressions.

### Step 2 — generate LaTeX and compile PDF

```bash
# Generate .tex inside Docker:
docker run --rm \
    -v "$(pwd)":/work \
    cadabra2-ubuntu:24.04 \
    python3 /work/ch34_export_latex.py

# Compile on host:
pdflatex ch34_left_right_spinors.tex
pdflatex ch34_left_right_spinors.tex   # second pass for ToC/refs

# Or use latexmk for fully automatic multi-pass:
latexmk -pdf ch34_left_right_spinors.tex
```

### One-liner (generate + compile, any chapter NN):

```bash
NN=34
docker run --rm -v "$(pwd)":/work cadabra2-ubuntu:24.04 \
    python3 /work/ch${NN}_export_latex.py && \
pdflatex -interaction=nonstopmode ch${NN}_*.tex && \
pdflatex -interaction=nonstopmode ch${NN}_*.tex
```

---

## Cadabra2 Cheat Sheet (for these scripts)

```python
import cadabra2
from cadabra2 import Ex, __cdbkernel__

__cdbkernel__ = cadabra2.create_scope()

# Declare index types
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe}"),            Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho}"),        Ex(r"position=free"))

# Declare tensor symmetries
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.Symmetric(Ex(r"g_{\mu\nu}"))

# Build an expression
e = Ex(r"\epsilon^{\alpha\beta} \psi_{\alpha} \chi_{\beta}")

# Get LaTeX string (key for the exporter):
print(e._latex_())   # → \epsilon^{\alpha \beta} \psi_{\alpha} \chi_{\beta}

# Apply cadabra2 algorithms:
cadabra2.sort_product(e)
cadabra2.canonicalise(e)
cadabra2.substitute(e, Ex(r"\psi_{\alpha} -> \chi_{\alpha}"))
cadabra2.collect_terms(e)
```

Key things cadabra2 does that sympy doesn't:
- **Index position tracking** (upper vs lower, undotted vs dotted)
- **Young tableau canonicalization** of tensor symmetries
- **Automatic dummy index renaming**
- **Spinor algebra with Fierz rearrangements** (Ch. 35+)
- Exports every expression via `._latex_()` → drop straight into `.tex`

---

## Chapter Roadmap: Ch. 34 → MHV

| Chapter | Topic | Key outputs |
|---------|-------|-------------|
| **34** | Left/right-handed Weyl fields, dotted/undotted indices, ε symbol, σ^μ | ✅ done |
| **35** | σ-bar algebra, index-free dot/bar notation, Weyl Lagrangian, equations of motion | |
| **36** | Quantization of Weyl fields; Majorana fermions | |
| **37** | Dirac fields, Dirac Lagrangian, γ matrices | |
| **38** | Spinor summation, trace technology | |
| **39** | Parity, time reversal, charge conjugation for spinors | |
| **48** | Spinor helicity formalism — angle/square bracket variables λ, λ̃ | |
| **50** | Spinor-helicity amplitudes for massless particles | |
| **MHV** | Parke-Taylor formula, BCFW recursion, MHV amplitudes | |

Chapters 40–47 (QED, renormalization of spinors) can be skipped or done in
parallel — they're not needed to reach MHV.

---

## LaTeX Export Pattern

Every `chNN_export_latex.py` follows this pattern:

```python
# 1. Run the same numpy/cadabra2 setup as chNN_<topic>.py
# 2. Collect results: matrices as pmatrix, cadabra exprs as _latex_() strings
# 3. Build a big doc string with f-string interpolation
# 4. Write to /work/chNN_<topic>.tex
# 5. Host compiles: pdflatex chNN_<topic>.tex
```

The helper functions used in every exporter:

```python
def cdb(expr_str):
    """Return LaTeX string for a cadabra2 expression."""
    return Ex(expr_str)._latex_()

def fmt_complex(z, tol=1e-10):
    """Format complex number as clean LaTeX fraction."""
    # → "0", r"\tfrac{1}{2}", r"\tfrac{i}{2}", etc.
    ...

def mat2pmatrix(mat):
    """2×2 numpy array → LaTeX pmatrix."""
    ...
```

---

## Research Direction: Toward MHV

The chapters above build toward the **spinor-helicity formalism** in which
amplitudes become rational functions of angle/square brackets:

```
⟨ij⟩ = ε^{αβ} λ_iα λ_jβ      (left-handed helicity spinors)
[ij]  = ε_{α̇β̇} λ̃_i^α̇ λ̃_j^β̇  (right-handed)
```

The **MHV (Maximally Helicity Violating)** n-gluon tree amplitude is:

```
A_n(1⁻ 2⁻ 3⁺ … n⁺) = i ⟨12⟩⁴ / (⟨12⟩⟨23⟩⟨34⟩…⟨n1⟩)
```
(Parke-Taylor formula, 1986)

Key tools in scope for Cadabra2 code:
- **BCFW recursion** — on-shell recursion building MHV from 3-point amplitudes
- **Grassmannian / momentum twistors** — geometric interpretation of MHV
- **Soft/collinear limits** — universal factorization properties verifiable symbolically
- **Gravity MHV** — KLT relations connecting gravity and gauge amplitudes

---

## Git Workflow

All Srednicki work lives on `feat/srednicki-ch34` (will become
`feat/srednicki-spinors` or similar as chapters accumulate).

```bash
# Ernest merges to master when satisfied
git push origin feat/srednicki-ch34

# For new chapters, branch off current:
git checkout -b feat/srednicki-ch35
```

Commit convention:
```
feat(Cadabra2): Srednicki Ch.35 — sigma-bar algebra and Weyl Lagrangian
```

---

## Notes on the MMD Source

The Srednicki book is parsed via nougat OCR to `.mmd` at:
```
workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/
  [Mark_Srednicki]_Quantum_Field_Theory(BookSee.org).mmd
```

The LaTeX source (`.tex`) is preferred for reading if available.
The MMD parse is mostly clean but occasionally has ambiguous index positions
in tensor formulas — always verify equations numerically as a sanity check.

---

*Last updated: 2026-02-26*
