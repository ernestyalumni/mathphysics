# mhvamplitudes

Python implementation of single-minus gluon amplitude (SMGA) calculations from
[arXiv:2602.12176](https://arxiv.org/abs/2602.12176) (Guevara, Lupsasca, Skinner, Strominger,
Weil), used as a verification testbed: two independently implemented code paths — a
closed-form formula and the paper's specialized Berends-Giele-style recurrence — are compared
against each other on random kinematics.

## Truth label

This package reproduces a known result from a published paper. Agreement between the two code
paths through `n = 10` is a reproducibility/engineering milestone, not a new physics result and
not evidence beyond the paper's own all-`n` formula and proof.

## Package layout

- `src/mhvamplitudes/kinematics/` — spinor-helicity kinematics and phase-space point generation.
- `src/mhvamplitudes/amplitudes/` — the two independent code paths (`smga.py` closed form,
  `berends_giele.py` recurrence).
- `src/mhvamplitudes/spinors/` — spinor/gamma-matrix machinery, including a `cadabra2/`
  subpackage used only by the (separate, Docker-based) Cadabra2-backed verification path. Not
  required for the tests below.

## Quick start

From this directory, in a fresh checkout:

```bash
uv run pytest -q tests/unit_tests/amplitudes tests/integration_tests/test_smga_reproduction.py
```

No Docker, Cadabra2, or manual `--with` flags required — `numpy` and `pytest` are declared as
project/dev dependencies and `uv` resolves them automatically.

See [`TESTING.md`](TESTING.md) for the full-suite command and the separate Docker workflow
that will be needed when Cadabra2-backed verifier tests are added.
