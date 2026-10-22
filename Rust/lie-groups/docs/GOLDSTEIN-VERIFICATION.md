# Goldstein rotation verification — 2026-10-21

The Rust, C++ and Python extensions implement the convention dictionaries in
Goldstein §§4.2–4.7 and the shared mono-notes. They extend the existing toolkit;
Propulsion/Cosmos source and the book corpus were not changed.

## Reproduce

From the mathphysics repository root:

```sh
bash Python/rotations/validation/run.sh
bash Rust/lie-groups/validation/run.sh
bash documents/mono-notes/build.sh
```

The Python runner uses the existing uv-created root `.venv`; its README gives
an isolated `uv sync` alternative. SymPy 1.14.0 ran in that environment with
Python 3.13.3. Cadabra ran separately in the existing, network-disabled Docker
image: its CPython 3.12 native binding cannot import into the root uv interpreter.
No packages were installed or environments changed during this work.

Completed independent root runs:

- `/tmp/mathphysics-goldstein-python-root-20261021`: new integrated checks;
  `provenance.json` binds 44 tested source/binary paths. All hashes were rechecked.
- `/tmp/mathphysics-goldstein-full-20261021`: full existing regression plus new
  native tests; all 58 recorded input hashes were rechecked.
- Durable text/JSON receipts are under
  `~/.cache/mathphysics-lie-groups/runs/20261021-goldstein/{integrated,full}/`.
  Build trees/binaries remain external to the repository; the recorded binary
  hashes identify the executions even after temporary build trees are removed.

## Passed checks

| Layer | Executed evidence |
|---|---|
| Rust | 3 unit tests, 40 integration tests, 4 doctests in debug/release; formatting and Clippy with warnings denied |
| C++ | All 3 CTests in Debug, Release and ASan+UBSan with LeakSanitizer enabled; 10 new Goldstein groups include float/double/long double |
| Python | 10 groups, including 3021 Euler reconstructions and 240 axis-angle/vector cases, screw actions, shape/domain rejection |
| New cross-language fixtures | Rust 24 and C++ 18 records compared independently with Python; largest component error `3.9429924734139306e-16`, bound `5e-11` |
| New SymPy oracle | 34 exact identity gates plus 15 comparisons against Rust-owned multi-angle exports; three negative controls |
| New Cadabra oracle | Three exact zero indexed reductions; deliberately wrong triple-product sign remains nonzero |
| Existing independent checks | 25 Taylor-exponential fixtures, largest error `5.551115123125783e-16`; SymPy 80 checks; Sage 86 checks; existing Cadabra and full QFT regression pass |
| Shared notes | Both editions rebuilt: 492 portrait / 267 screen pages, zero LaTeX warnings; pre-existing overflow counts unchanged at 96 / 52 |

The source/PDF hashes match the notes' final build report. Screen pp.146–147
were visually reviewed. An independent mathematical review checked Euler's
fixed-axis proof, normal-form orientation, Goldstein's clockwise sign,
identity/half-turn axis recovery, the principal half-angle branch, and the
Chasles screw-point formula.

The tests distinguish passive body-to-inertial coordinate conversion from
active vector rotation, enforce product-order Euler arguments, and compare
reconstructed matrices rather than nonunique angle triples or quaternion signs.
They exercise identity, half-turns, near-zero/near-pi angles, all Euler poles,
large finite angles, reflections, nonfinite data, malformed shapes and CK sign
errors. Tests found and fixed a ZYZ extraction sign at theta=pi and Python
conversion tolerance-closure failures. Matrix helpers reject silent zip truncation.

## Mathematical domains

- Numerical values approximate SO(3)/SU(2) within explicit tolerances. Rust
  inverse extraction is f64; forward maps support f32/f64/software Half. C++
  supports float/double/long double; Python uses standard floating scalars.
- Inverse Euler maps select principal branches and mark chart singularities.
  Matrix-to-quaternion extraction cannot recover a lost lift sign. Direct Euler
  quaternion products and CK imports preserve distinct q/-q lifts.
- Rust `MultiTrig<N>` is a bounded exact real circle quotient for up to eight
  independent angles; coefficients inherit checked-i128 limits. C++ retains
  arbitrary-precision polynomial and one-common-angle quotient support.
- SymPy checks use explicitly stated circle/unit-quaternion/unit-axis ideals.
  Cadabra checks indexed Euclidean epsilon algebra. These are exact identities
  within those domains, not all-input floating-point proofs or a general CAS.

The original Sage SO/QQ reflection diagnostic remains independently gated by
determinant+1; class names are not accepted as group-membership evidence.
