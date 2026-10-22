# INERTIA — local rigid body flight lab

Interactive WebGL view of the authoritative Python 6-DOF rigid-body solver.
The browser renders state and sends body-frame wrenches; it never integrates
physics or replaces the numerical solver with an animation.

From the `mathphysics` repository root:

```sh
./run_rigid_body_lab.sh
```

Open the printed loopback URL, normally **http://127.0.0.1:8765**. The launcher
builds this frontend into `/tmp/mathphysics-rigid-body-viewer-dist` and serves
it with the solver on the same origin. `--no-build` reuses an existing build;
`--port 8766` selects another local port. This is a single-pilot local sandbox;
multiple browser tabs share its server state.

Dependencies are pinned in `package-lock.json`. After installation, the built
app uses no CDN, external fonts, assets, or runtime network services.

```sh
cd JavaScript/rigid-body-viewer
npm ci
npm test
npm run build
```

For frontend development, run the solver and `npm run dev`; Vite proxies `/api`
to `http://127.0.0.1:8765`. Override that target with `RIGID_BODY_API` when using
another solver port. `npm run build` writes outside the repository.

## Fly

Body axes are **+X nose, +Y starboard/right, +Z belly/down**. All frames remain
right-handed; the camera uses inertial -Z as its initial up direction.

| Action | Keys | Body-frame command |
|---|---|---|
| Forward / reverse | W / S | force ±X |
| Starboard / port | D / A | force ±Y |
| Descend / rise | F / R | force ±Z |
| Nose up / down | ↑ / ↓ | torque ±Y |
| Yaw right / left | → / ← | torque ±Z |
| Roll right / left | E / Q | torque ±X |
| Translation boost | Shift | 4× force, unchanged torque |
| Pause / resume | Space | preserves physical state |
| One solver step while paused | . | one fixed step, 1/240 s |
| Orbit / chase camera | C | changes the view only |

On-screen flight buttons apply the same momentary commands. Drag to orbit,
scroll to zoom; the chase camera follows body orientation. The camera follows
the body position while the inertial grid, angular momentum vector and trail
provide reference. The blue-green velocity arrow shows translational motion;
the amber momentum arrow shows inertial angular momentum. Body axis colors
are coral X, mint Y, and periwinkle Z. Amber wrench arrows and cyan exhaust
appear only while the server reports an applied wrench.

**Hold to fire; release to coast.** There is no drag, gravity, velocity damping
or attitude stabilization. Holding torque increases spin; release removes the
external torque and leaves free rigid-body dynamics. With no input, external
force and torque are exactly zero. All held controls are released on blur or
page hiding, and a zero-wrench request is sent. Flight keys are ignored while
editing input fields or selecting simulation speed.

The four presets demonstrate free tumble, stable principal-axis spin,
intermediate-axis instability and initially stationary thruster flight.
Changing mass, cuboid dimensions or initial body angular velocity and choosing
**Apply & restart** starts a fresh initial condition. Rate changes how much
simulation time is advanced per wall-clock interval; it does not change the
fixed solver step. Pause and single-step do not introduce damping.

Energy and momentum diagnostics use the server's coast-segment reference.
Under external thrust, those quantities can physically change. The displayed
relative drifts use denominator `max(reference magnitude, 1 SI unit)`, making
an initially zero reference meaningful rather than dividing by zero.

## Local API and browser verification

All requests are same-origin JSON. The backend owns validation and state.

- `GET /api/state`: current snapshot.
- `POST /api/reset`: `{preset,mass,dimensions,omega_B?}`; dimensions in metres,
  mass in kg, angular velocity in rad/s. The preset supplies initial rates if
  `omega_B` is omitted.
- `POST /api/advance`: `{wall_dt,force_B,torque_B,speed,paused,step?}`. Wrench is
  body-frame N and N·m. A single step advances exactly one backend fixed step.

Snapshots include `r_I`, `v_I`, scalar-first Hamilton `q=[w,x,y,z]`, `omega_B`,
`sim_time`, model configuration, `control`, and conservation `diagnostics`.
Rendering assigns Three.js quaternion `(x,y,z,w)` without changing handedness.

`window.__simTest` exposes read-only `state`, `keys`, `wrench`, `requestCount`,
`lastRequest`, `errors`, `connected`, `webgl`, `cameraMode`, and `paused` for
browser checks. It has no physics mutation or stepping hooks.

The five Node unit tests exercise every command sign, zero wrench, opposing
commands, boost, input validation and independent pointer/keyboard holds.
The `tests/browser/` Playwright tests use the running real backend and
installed Chrome:

```sh
SIM_URL=http://127.0.0.1:8765 npm run test:browser
```

No browser download is required. Set `CHROME_PATH` if Chrome is installed
elsewhere. The simulator must be running; E2E tests change its shared state.
