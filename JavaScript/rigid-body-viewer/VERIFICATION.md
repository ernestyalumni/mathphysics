# Six-DOF lab verification — 2026-10-21

The local production viewer and authoritative Python backend are implemented.
Run `./run_rigid_body_lab.sh` from the repository root, then open
<http://127.0.0.1:8765>. See [README.md](README.md) for commands and controls.
The Python, C++ and Rust libraries share the same body-to-inertial, scalar-first
Hamilton quaternion convention. The browser uses the Python implementation.

## Executed evidence

- Python: 26 tests, including existing rotation regressions, 10 dynamics groups
  and 6 server/HTTP groups. HTTP tests use a real ephemeral local server; batch
  solver failure leaves the previous state intact.
- Native dynamics: eight fixture cases each for C++ and Rust, compared with the
  independently written Python implementation. Maximum state/energy/inertial
  spin discrepancy `2.34e-14` (threshold `2e-10`).
- Independent SciPy DOP853: eight free/forced cases, including nondiagonal
  inertia and a nonidentity initial attitude. Closed spherical/principal and
  full symmetric-top attitude solutions; full intermediate-axis-flip test.
  Halving timestep divides symmetric-top attitude errors by `3.9984, 3.9996`.
- Exact SymPy: midpoint energy, torque work and Cayley/momentum inverse
  identities, plus the existing 34 rotation gates and 15 exact-ring records.
  Cadabra: three zero identities and a wrong-sign nonzero negative control.
- C++: four CTests pass Debug, Release and ASan/UBSan/LeakSanitizer.
  Rust: full Debug/Release suite (48 integration, 3 unit, 4 doctests), final
  eight dynamics groups, formatting and all-targets Clippy pass.
- Frontend: five Node command tests pass. Playwright with actual local Chrome
  renders WebGL and exercises keyboard thrust/torque, coasting, blur release,
  form typing, initial conditions, pause/exact single step, presets, speed,
  cameras and mobile touch. No page errors or API errors. Desktop and mobile
  screenshots were visually inspected. Browser tests begin from an explicit
  reset because the server intentionally has one shared pilot state.
- Both notes PDFs rebuilt: 493 portrait / 267 screen pages; zero LaTeX warnings.
  The new simulation material fits screen page 148. Global TOC number-width
  overflows remain (97/53 boxes overall); no new simulation-equation overflow.

The default cuboid (`m=100 kg`, dimensions `3×2×1 m`, initial body angular
velocity `[.12,.32,.5] rad/s`) was propagated for 60 seconds at `dt=1/240 s`.
Independent DOP853 comparison: attitude error `2.22585e-5 rad` (`0.00127532°`),
angular-velocity error norm `1.34902e-6 rad/s`, relative energy drift `9.50e-13`
and relative inertial-spin drift `5.69e-13`. Runtime was 6.39 seconds on this
machine; this is one measurement, not a performance guarantee.

## Reproduce

```sh
ROTATIONS_RUN_DIR=/tmp/mathphysics-sixdof-check bash Python/rotations/validation/run.sh
npm --prefix JavaScript/rigid-body-viewer test
./run_rigid_body_lab.sh
# In a second terminal, with the real simulator running:
npm --prefix JavaScript/rigid-body-viewer run test:browser
```

The rotation validation runner uses the repository uv environment and installed
local Cadabra Docker image. See [Python dynamics](../../Python/rotations/DYNAMICS.md)
and the [rotation runner](../../Python/rotations/README.md) for dependencies.
No Cadabra or symbolic CAS is required to run the simulator.

The final integrated receipt is `/tmp/mathphysics-sixdof-final-20261021`.
Durable text/JSON evidence and screenshots are copied under
`~/.cache/mathphysics-sixdof/runs/20261021/`; source and binary hashes identify
what was tested. Compiled binaries and node_modules remain outside versioned
source.

## Scope

This verifies an ideal constant-mass rigid body, with fixed physical COM
inertia and explicit net external wrenches. The integrator is second order;
conservation does not certify arbitrary long-horizon phase accuracy. The
viewer models a uniform cuboid; surface details are visual decoration.
Gravity, contacts, flexibility, fuel/slosh, actuator allocation and navigation
are not part of this baseline. No NASA check-case dataset has yet been run;
those are a later verification target, not evidence claimed here.
