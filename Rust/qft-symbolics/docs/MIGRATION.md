# Python → Rust audit and remaining port map

2026-10-16. Static source audit plus the explicitly recorded runtime check.
Old Python files have not been "fixed" or removed.

## What was wrong with the previous verification story

Paths below are relative to
`amplitudes/Srednicki/individual_chapters/`.

### Display is not symbolic verification

- `ch50_massless_spinor_helicity.py: 84–88`: constructs and displays Schouten
  in Cadabra, but never transforms a residual to zero. There are no numerical
  assertions in this file.
- Same file, lines 139–143: prints Ward success and comments that choosing
  q=k gives zero polarization. The displayed formula divides by <qk>; q=k
  is singular, not a valid zero-polarization check.
- The script's "COMPLETE" ending cannot establish chapter completion.
- `ch35_sigma_algebra.py` contains real numerical residual diagnostics,
  but most failures only print a mismatch. Only epsilon normalization is
  asserted. Its line 646 trace summary has the wrong sign compared with
  its own implemented -2g check at lines 296–316.
- `ch41_gamma_technology.py: 46, 63, 90`: defines mostly-plus eta and gamma
  matrices satisfying {gamma,gamma}=-2eta, then compares against +2eta.
  Lines 134–147 use the corresponding wrong two-gamma trace sign. These are
  nonfatal diagnostics, not enforced acceptance checks.
- That file's full Dirac Fierz ansatz at 293–303 is not a valid general
  identity; a diagonal-component counterexample gives 0 on the gamma side
  and -2 on the proposed delta side. At 315–320, PL gamma PL is zero, but
  is compared with a nonzero chiral tensor. Do not port those expected values.

### Chapter 60 really does contain assertions—but does not finish

`ch60_spinor_helicity.py` has substantial NumPy checks: Dirac equations,
rank-one factorization, bracket/trace relations, transversality, Fierz and
Schouten assertions. It must not be dismissed as display-only.

However, its general constructor at 585–587 puts the same two-component
spinor into opposite chiral blocks. It disagrees with the earlier specialized
constructor. At line 603 it asserts that the resulting angle spinor solves
the massless Dirac equation.

Observed command (read-only source mount; Docker network disabled):

```bash
docker run --rm --network none \
  -v "$PWD/amplitudes/Srednicki/individual_chapters:/chapters:ro" \
  cadabra2-ubuntu: 24.04 python3 /chapters/ch60_spinor_helicity.py
```

**Observed exit 1:** `AssertionError: q/ |q⟩ ≠ 0!`.
Earlier bracket/trace checks passed; later assertions were not reached.
Log kept outside the repository at `/tmp/mathphysics-original-ch60.log`.

Additional static issues: the line-648 normalization success text contradicts
the cross-product assertions immediately above it; the polarization
completeness expression is displayed without an expected-tensor assertion;
line 928 says "Verified algebraically" without performing that verification.
Do not inherit the final "ALL VERIFIED" summary.

## Port map

| Existing Python | Rust starting point | Status |
| --- | --- | --- |
| ch35 epsilon/sigma and left generators | matrices + algebra tests | Exact Q(i) components, not generic index canonicalisation |
| ch41 gamma technology | clifford + matrices + textbook tests | Rebuilt from Clifford definitions; old sign/Fierz expectations rejected |
| ch50 display-only overview | spinor + twistor | Actual polynomial Schouten, completeness and slash checks |
| ch60 genuine but blocked calculations | twistor + textbook tests | Polarization identities 60.1a/b exact; no singular q=k |
| MHVamplitudes/kinematics/spinor_helicity.py brackets | smga::bracket | Explicit preserved sign, no copied explanatory claims |
| MHVamplitudes/amplitudes/smga.py | smga::formula_r1 | Rational/conserved R1 domain checked |
| MHVamplitudes/amplitudes/berends_giele.py | smga::recursion | Interval dynamic programming, not line-by-line translation |
| Angular massless spinor constructor | twistor::from_rational_angles | Validated exact rational parameterization, 48 witnesses; not full Lorentz-package parity |
| Typed index/property core | tensor | Independently implemented bounded symmetry orbit, metric/delta contractions and hygienic dummies |
| Ordinary QCD vertices / BG / BCFW amplitude recursion | — | Not ported |

The old package's SMGA docstrings incorrectly cite chapter 60 as containing
nonabelian color decomposition/Parke–Taylor gluon formulas, and make overly
broad three-point single-minus vanishing claims. Preserve the tested
algorithms, **not those explanations**. Chapter 60 is QED. Three-point
complexified branches must be distinguished from generic real Lorentzian
higher-point tree statements.

## Why this isn't a full Cadabra port

Cadabra separates expression storage, contextual properties, index
classification, signed permutation-group canonicalisation, metric/delta
elimination, distribution, gamma joins and term collection.
Relevant source paths inspected: `core/Storage.hh`,
`core/IndexClassifier.cc`, `core/algorithms/canonicalise.cc`,
`eliminate_metric.cc`, `eliminate_kronecker.cc`, `join_gamma.cc`,
`sort_spinors.cc`, `fierz.cc` and `core/Exchange.cc`.

The hard parts are interacting dummy renamings, slot symmetries, statistics,
variance, index spaces and derivative dependencies—not parsing LaTeX.
Sorting symbols alphabetically would not reproduce that machinery.

This crate starts with the narrower algebra needed for spinor-helicity and
small gamma calculations. The second slice adds the independent bounded typed
tensor kernel and specialized formal-d gamma sandwich identities. See
[TENSOR.md](TENSOR.md) and [TEXTBOOK-EXTENDED.md](TEXTBOOK-EXTENDED.md).
No general tensor equivalence, arbitrary Fierz, Majorana flip rule or full
dimensional-regularization scheme is claimed; no Cadabra source was translated.
