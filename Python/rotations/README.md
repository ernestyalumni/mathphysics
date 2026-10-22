# Rotations: Goldstein §§4.2–4.7

Reusable numerical Python functions and independent CAS/native-language checks.
The rotation library uses only the standard library; exact verification uses
SymPy. The ideal [6-DOF rigid-body simulator](DYNAMICS.md) uses NumPy; its
independent continuous reference uses SciPy DOP853. It complements `Rust/lie-groups` and `CppMath/LieGroups`.

```python
from rotations import from_euler, to_euler, to_axis_angle, to_quaternion
R = from_euler('zyx', 0.4, 0.7, -0.3)  # heading, pitch, roll, radians
angles = to_euler(R, 'zyx')
axis, phi = to_axis_angle(R)
q = to_quaternion(R)                  # scalar-first Hamilton
```

Conventions:

- Right-handed column vectors, `x_I = R x_B`; the same numerical matrix acts
  actively on a vector in a fixed frame. Goldstein's reverse coordinate matrix
  is `A_G = R.T`.
- Intrinsic product order: `ZXZ(phi,theta,psi)`, `ZYZ(phi,theta,psi)`,
  `ZYX(heading,pitch,roll)`. The last is the usual engineering yaw/pitch/roll
  product `Rz(heading) Ry(pitch) Rx(roll)`; signs still depend on frame axes.
- Proper Euler extraction: middle angle `[0,pi]`. ZYX extraction: middle angle
  `[-pi/2,pi/2]`. At lock choose last angle zero and return `singular=True`.
  Near-lock reconstruction error is `O(singular_tolerance)`; angles themselves
  are ill-conditioned there. Compare reconstructed matrices.
- Axis-angle extraction: angle `[0,pi]`; at identity choose `(1,0,0)` as arbitrary
  axis. The axis sign at pi is arbitrary. Largest-diagonal quaternion extraction
  avoids dividing by `sin(phi)` near pi; quaternion sign is geometrically free.
- CK `(alpha,beta,gamma,delta)` obeys `gamma=-conj(beta)`,
  `delta=conj(alpha)`, `|alpha|²+|beta|²=1`;
  `alpha=w+iz`, `beta=y+ix`, `K=U(q).dagger`, `Phi(K)=R(q).T`.
- Imports reject nonfinite values, reflections and invalid constraints. Matrix
  conversion revalidates output and can reject accumulated tolerance drift.
  Accepted data is not silently projected onto SO(3), except the documented
  normalization in matrix-to-quaternion extraction.

All three languages use **product-order angle arrays**; the ZYX tuple is
`[yaw,pitch,roll]`, not `[roll,pitch,yaw]`.

## Run

From the repository root, reuse its existing **uv-created** `.venv` (Python
3.13.3 and SymPy 1.14.0 in the tested environment):

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=Python/rotations .venv/bin/python -m unittest discover -s Python/rotations/tests -v
bash Python/rotations/validation/run.sh
```

No packages were installed or environment contents changed during this work.
For another machine, an isolated uv environment is supported by the local
`pyproject.toml`:

```bash
uv sync --project Python/rotations --group verification
ROTATIONS_PYTHON="$PWD/Python/rotations/.venv/bin/python" bash Python/rotations/validation/run.sh
```

The full runner builds native producer binaries offline under a fresh `/tmp`
receipt directory, runs cross-language comparison and exact SymPy checks, then
uses the existing **local Docker Cadabra2 image**, without network. Cadabra2 is
not installed in the root uv environment: no host Python binding/wheel was found.
The image binding is CPython 3.12-specific (`cadabra2.cpython-312`), whereas
the root uv interpreter is CPython 3.13; it cannot be imported into that venv.
The Docker backend is a separate external oracle, not a uv-installed Cadabra
binding. `WITH_CADABRA=0` explicitly records SKIPPED if that image is unavailable.
`ROTATIONS_RUN_DIR` selects another receipt path; build files stay outside the
source tree. Full pre-existing Rust/C++ regression is separate:
`bash Rust/lie-groups/validation/run.sh`.

## Evidence and bounds

Ten Python test groups include 3021 Euler branch reconstructions, 240
axis-angle/physical-vector cases, screw actions, conjugation, invalid imports,
shape rejection and wrong-sign/order controls. Rust and C++ producers are
compared independently against Python formulas; native extraction outputs are
reconstructed rather than compared as ambiguous angle tuples.

The SymPy oracle explicitly constructs matrices and Hamilton quaternion products.
Three independent circle ideals verify all Euler products. Separate unit-
quaternion/unit-axis polynomial ideals verify orthogonality, determinant,
Goldstein CK matrix entries, hat equivariance, the SO(3) characteristic polynomial,
actual Cayley–Hamilton, Rodrigues and screw identities. It also checks unrestricted
3×3 Cayley–Hamilton and the Rust-owned `MultiTrig` exports. These are exact
residual checks in the stated polynomial domains; finite precision stress tests
are not proofs for all floating inputs.

The Cadabra oracle reduces indexed Euclidean epsilon contractions for
`n × (n × x)`, `hat(n)^3 = -|n|²hat(n)`, and the fixed axis, and retains a wrong
sign as nonzero. It does not infer group membership from a symbol's name.

`validation/provenance.py` hashes tested source inputs and native producer
binaries, records interpreter/SymPy versions and Git HEAD. The receipt's
logs and hashes identify actual executed artifacts; no generated receipt is
committed. No general-purpose CAS, globally continuous Euler chart, dynamics
integrator or complete textbook implementation is claimed.
