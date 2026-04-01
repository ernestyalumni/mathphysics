# Cadabra2 Spinor Technology Examples

Cadabra2 Python API examples for spinor/tensor calculations in quantum field
theory, following Srednicki's QFT notation (Chapters 34-38).

## Files

| File | Topic |
|------|-------|
| `01_weyl_spinors.py` | 2-component Weyl spinors, epsilon tensors, index raising/lowering, spinor inner products |
| `02_sigma_matrices.py` | σ^μ and σ̄^μ matrices, completeness, momentum matrices, Clifford algebra |
| `03_spinor_helicity.py` | Spinor-helicity formalism for massless particles, angle/square brackets, Mandelstam variables |
| `04_mhv_amplitudes.py` | MHV amplitudes, Parke-Taylor formula, 4-gluon amplitude, cyclic invariance |
| `05_fierz_identities.py` | Fierz rearrangement, sigma completeness, SUSY spinor identities, Schouten identity |
| `REFERENCES.md` | Key papers on MHV amplitudes and spinor-helicity formalism |

## Background

**2-component (Weyl) spinors** are the natural language for massless particle
physics and supersymmetry. The van der Waerden notation uses:
- Undotted indices α, β, γ, δ for the left-handed (1/2, 0) representation of SL(2,C)
- Dotted indices α̇, β̇, γ̇, δ̇ for the right-handed (0, 1/2) representation
- Epsilon tensors ε_{αβ}, ε^{αβ} to raise and lower spinor indices

**Spinor-helicity formalism** exploits the fact that massless 4-momenta
factorize as p_{αα̇} = λ_α λ̃_{α̇}. The angle brackets ⟨ij⟩ and square
brackets [ij] provide Lorentz-invariant building blocks that make gauge
theory amplitudes dramatically simpler than Feynman diagrams suggest.

**MHV amplitudes** (Maximally Helicity Violating) are gluon scattering
amplitudes with exactly 2 negative-helicity gluons. The Parke-Taylor formula
gives a closed-form expression:

```
A_n(1^-, 2^-, 3^+, ..., n^+) = i g^{n-2} ⟨12⟩^4 / (⟨12⟩⟨23⟩...⟨n1⟩)
```

This is one formula where n-gluon tree amplitudes that would otherwise require
thousands of Feynman diagrams reduce to a single expression.

## Srednicki References

- **Chapter 34**: Two-component (Weyl) spinors, van der Waerden notation
- **Chapter 35**: Dirac spinors from Weyl spinors; σ^μ matrices
- **Chapter 36**: Helicity spinors, completeness, momentum spinors
- **Chapter 37**: Spinor-helicity formalism, polarization vectors
- **Chapter 38**: MHV amplitudes, Parke-Taylor formula, helicity Feynman rules
- **Appendix B**: Fierz identities, sigma matrix products

## How to Run

### With cadabra2 installed locally:
```bash
python3 01_weyl_spinors.py
python3 04_mhv_amplitudes.py
```

### With Docker:
```bash
# Build the image first (from Deployments/DockerBuilds/Physics/Cadabra2/):
docker build -t cadabra2-ubuntu:24.04 .

# Run a Python example:
docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 python3 /work/04_mhv_amplitudes.py

# Interactive session:
docker run --rm -it -v $(pwd):/work cadabra2-ubuntu:24.04 python3

# Cadabra2 GUI notebook (requires X11):
../../Deployments/DockerBuilds/Physics/Cadabra2/run_gui.sh $(pwd)
```

## Physics Notes

### Why 2-component spinors?

4-component Dirac spinors contain redundancy for massless particles.
The left- and right-handed Weyl spinors are irreducible representations of
the Lorentz group and are the natural building blocks. For massive particles,
the Dirac spinor is a pair (ψ_α, χ̄^{α̇}) joined by a mass term.

### The spinor-helicity miracle

For n gluons at tree level, Feynman diagrams give O(n!) terms that combine
into the Parke-Taylor formula — a single fraction. This simplicity is not
accidental: it reflects the twistor space structure of Yang-Mills theory
(Witten 2004) and the existence of on-shell recursion relations (BCFW 2005).

### BCFW recursion

The Britto-Cachazo-Feng-Witten recursion builds n-point amplitudes from
products of lower-point amplitudes. Together with Parke-Taylor as the
seed (3-point MHV and MHV-bar), BCFW generates all tree-level gauge
theory amplitudes from first principles without Feynman diagrams.

See `REFERENCES.md` for the key papers.
