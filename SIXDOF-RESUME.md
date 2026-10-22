# Resume: rotations, mono-notes, and the six-DOF lab

Saved 2026-10-22 at Ernest's request before closing the session.

## Start here

Repository: `/home/propdev/.openclaw/workspace/workspace2/repos/mathphysics`
Branch: `feat/lie-groups-so3-su2`.
Ernest requested a local feature-branch commit before closing. The relevant
source, LaTeX notes, tests, documentation, and this handoff are included in the
commit containing this file; identify it with `git log -1 --oneline`.
No push or merge was requested. Generated PDFs and the local build report are
left on disk outside the commit; rebuild them with
`bash documents/mono-notes/build.sh`. Inspect `git status` before editing or
changing branches, and preserve any local changes.

```sh
cd /home/propdev/.openclaw/workspace/workspace2/repos/mathphysics
./run_rigid_body_lab.sh
```

Open <http://127.0.0.1:8765>. Check the URL before starting another server;
`--port 8766` selects a different port. Normal launch rebuilds the frontend
into `/tmp`; `--no-build` requires that previous build still to exist.
The launcher uses the existing root uv environment when available, otherwise
`uv run` with the Python project. npm installs pinned dependencies if missing.
The simulation is one shared pilot state across browser tabs.

Confirmed reachable at handoff: {"status": "ok", "model": "ideal rigid body", "fixed_dt": 0.004166666666666667}. Left running; closing the chat does not request a shutdown.

Server metadata, if still applicable: `/tmp/mathphysics-rigid-body-lab.pid`
and `/tmp/mathphysics-rigid-body-lab.log`. A saved PID can be stale; verify its
command before stopping a process. No tests or implementation tasks were left
running. The interactive server was intentionally left available.

## Conventions and writing preferences

- `a_I = R a_B`: R converts body components to inertial components. This is a
  passive coordinate change; the same numeric matrix has an active vector
  interpretation in a fixed frame. Goldstein's passive matrix is `A_G = R^T`.
- Hamilton scalar-first `q=(w,x,y,z)`, body to inertial;
  `qdot = (1/2) q ⊗ (0,omega_B)`. Preserve continuous signed lifts; `q` and
  `-q` describe the same SO(3) attitude. Three.js receives `(x,y,z,w)`.
- Euler arrays are product order: ZXZ/ZYZ `[phi,theta,psi]`;
  standard aerospace 3-2-1 ZYX `[yaw psi,pitch theta,roll phi]`.
  "Our notation" does not mean a custom engineering convention.
- Cayley–Klein uses `gamma=-beta*`, `delta=alpha*`, with normalization;
  distinguish the signed SU(2) lift from its SO(3) image.
- Prefer phi/theta first for angles; use psi for a third angle. Notes should
  be concise, explicit about maps/domains and algebra, at an advanced
  mathematical-physics/Lie-group level. Use `c(phi)=cos(phi)`, `s(phi)=sin(phi)`
  when expressions grow. Keep the established direction-cosine index spacing.

## Completed artifacts

- [Mono-notes](documents/mono-notes/screen.pdf): Goldstein §§4.2–4.7,
  conjugation/Ad, Euler conventions and dictionary, TikZ diagram,
  Cayley–Klein/quaternions, Euler theorem, Rodrigues/exponential,
  Cayley–Hamilton, axis-angle endpoint cases and Chasles/screw motion.
- New [six-DOF theory source](documents/mono-notes/topics/free-rigid-body-simulation.tex):
  explicit 13-component constrained state, COM equations, invariants,
  midpoint/Cayley update and engineering reading. Screen PDF page 148.
  Last build: 493 portrait / 267 screen pages, zero LaTeX warnings;
  existing global layout overflows documented in the verification receipt.
- [Python dynamics](Python/rotations/DYNAMICS.md),
  `Python/rotations/rigid_body.py`, and `simulator_server.py`.
- [C++ library](CppMath/LieGroups/README.md): `include/rotations.hpp` and
  `include/rigid_body.hpp`.
- [Rust library](Rust/lie-groups/README.md): `src/goldstein.rs`,
  `src/multitrig.rs`, `src/rigid_body.rs` and matching tests/fixtures.
- [Browser frontend and controls](JavaScript/rigid-body-viewer/README.md).
  Three.js procedural WebGL, local assets, no runtime CDN. Four presets,
  keyboard/touch thrust and torque, orbit/chase camera, mass/dimensions/rates,
  pause/step/reset/time scale, conservation telemetry.

The numerical model has constant mass and full physical symmetric positive
COM inertia. It advances fixed `1/240 s` steps using implicit-midpoint body
momentum with a matching Cayley quaternion increment. External wrench is zero
by default. Releasing controls removes the wrench and preserves free coasting;
there is no artificial damping. The viewer's body is a uniform cuboid.
The ideal wrench controls do not implement actual thruster allocation or fuel.

## Verification and practical limits

Canonical receipt: [VERIFICATION.md](JavaScript/rigid-body-viewer/VERIFICATION.md).
Previous rotation receipt:
[Rust/lie-groups/docs/GOLDSTEIN-VERIFICATION.md](Rust/lie-groups/docs/GOLDSTEIN-VERIFICATION.md).
Durable JSON/text receipts, hashes and desktop/mobile screenshots:
`~/.cache/mathphysics-sixdof/runs/20261021/`.
Temporary integrated run: `/tmp/mathphysics-sixdof-final-20261021`.

- 26 Python/HTTP tests and 5 frontend command tests pass.
- Real Chrome desktop/mobile WebGL and keyboard/touch controls pass, including
  pause/exact single step, coast, torque, blur release and editable setup.
- C++ four CTests pass Debug/Release/sanitizers; Rust 48 integration, 3 unit,
  4 doctests pass Debug/Release, plus formatting and Clippy.
- Eight dynamics fixtures each from C++ and Rust agree with Python to
  `2.34e-14`; independent DOP853 and closed full-attitude solutions verify
  trajectory accuracy and second-order convergence. Full intermediate-axis
  flip, nondiagonal inertia, forcing and explicit solver failure are tested.
- Exact SymPy checks midpoint conservation/work and Cayley identities;
  existing symbolic rotation and Cadabra positive/negative checks pass.
- Default 60-second tumble at 240 Hz: attitude error `0.00127532 degrees`
  versus independent DOP853; energy/spin relative drift about `1e-12`.
  Conservation alone does not establish attitude phase accuracy.

This is a verified ideal rigid-body baseline. Measured mass properties,
environment, actuators, flexibility and slosh are separate physical modeling
extensions. No NASA check-case dataset has yet been executed.
SymPy/SciPy use the existing uv environment; Cadabra uses the installed
network-disabled Docker image because its Python binding differs from the
uv interpreter. CAS tools are not required to run the interactive simulator.

## Where to pick up

The requested baseline is complete; no unfinished implementation is implied.
First review/fly the current lab, then connect the equations to Goldstein
§§5.1, 5.3–5.6 and Wie §5.5, §§6.1–6.7. Wie is already available locally:
`Data/Public/books/EngineeringPhysics/Wie-SpaceVehicleDynamicsControl-2e/INDEX.md`
(relative to workspace2).
For an engineering continuation, Wie §§7.3–7.7 cover feedback/momentum
management/jet selection; Chapters 8–9 cover flexibility and slosh. Choose a
specific next physical model or controller with Ernest before extending scope.

A future agent should read this file and the verification receipt, preserve
all existing changes, check whether the server is still available, and resume
from Ernest's next instruction. No need to repeat all completed tests without
a source change or a new concern. No new work was started for this handoff.
