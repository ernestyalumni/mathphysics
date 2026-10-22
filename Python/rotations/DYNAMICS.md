# Ideal rigid-body simulation

`rigid_body.py` implements COM-based translation and Hamilton quaternion attitude:

\[
(r_I,v_I,q,\omega_B)\in\mathbb R^3\times\mathbb R^3\times S^3\times\mathbb R^3,
\quad R(q):B\to I,\quad m>0,\quad J_B=J_B^T>0,
\]
\[
\dot r_I=v_I,\quad m\dot v_I=R(q)F_B+F_I,\quad
\dot q=\tfrac12q\otimes(0,\omega_B),\quad
J_B\dot\omega_B=\tau_B-\omega_B\times J_B\omega_B.
\]

`J_B` is about the center of mass and may be nondiagonal. Imports also require
physical principal-moment triangle inequalities. Roundoff-sized asymmetry is
replaced by the symmetric part. Imported quaternion norm-squared must be within
`1e-10` of unity; accepted roundoff is normalized without flipping its sign.
There is no implicit gravity, drag, contact, mass variation or actuator model.

```python
from rigid_body import RigidBodyModel, RigidBodyState, step, diagnostics
model = RigidBodyModel.from_box(2.0, (1.0, 0.6, 0.4))
state = RigidBodyState(omega_body=(0.01, 1.0, 0.01))
for _ in range(240):
    state = step(model, state, 1/240)
print(diagnostics(model, state))
```

The discrete angular update solves in body angular momentum:

\[
h_k=J_B\omega_k,\quad h_m=\tfrac12(h_0+h_1),\quad
\omega_m=J_B^{-1}h_m,\qquad
h_1-h_0=\Delta t\,[\tau_B-\omega_m\times h_m],
\]
\[
\delta q=\frac{(1,\frac{\Delta t}{2}\omega_m)}
 {\sqrt{1+\frac{\Delta t^2}{4}|\omega_m|^2}},\quad
q_1=q_0\otimes\delta q,\quad q_m=q_0\otimes\sqrt{\delta q},
\]
\[
v_1=v_0+\frac{\Delta t}{m}[R(q_m)F_B+F_I],\qquad
r_1=r_0+\tfrac{\Delta t}{2}(v_0+v_1).
\]

Newton uses an analytic Jacobian and checks its residual explicitly. The default
bound is `1e-12*(1+||h0||inf+dt*||tau||inf)` in the chosen momentum units; the
leading `1` gives an absolute as well as relative tolerance. `step_with_info`
returns the state and `StepInfo(iterations,residual,residual_bound)`.
Nonconvergence, singular linear solves and nonfinite updates fail explicitly.
There is no arbitrary upper time-step cap.

For zero torque, the midpoint/Cayley pairing preserves the complete inertial
spin `H_I=R(q)J_B omega_B` and rotational energy to nonlinear solve tolerance.
With torque its discrete work identity is
`E1-E0=dt*omega_mid.dot(torque_body)`. The integrator is second order, and
**conserved energy does not imply accurate attitude phase**: even a spherical
top accumulates Cayley phase error. The reference tests measure full attitude,
angular velocity and convergence, rather than using invariants alone.

`angular_momentum_inertial` means spin about the COM. Diagnostics separately
report total angular momentum about the inertial origin, including `r_I × P_I`.
All forces/torques are constant during an individual step. Controls can vary
between steps. Body forces are rotated at the midpoint attitude; inertial forces
are entered separately, e.g. `force_inertial=(0,0,-mass*9.81)`.

## Independent checks

From the repository root, with the existing uv environment:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=Python/rotations .venv/bin/python -m unittest discover -s Python/rotations/tests -v
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=Python/rotations .venv/bin/python Python/rotations/validation/dynamics_oracle.py
```

The reference ODE is independently written and integrated by SciPy DOP853
(`rtol=2e-13`, `atol=2e-15`). Closed solutions cover spherical/principal rotations
and the complete symmetric-top attitude, including both precession and spin.
Other checks cover offdiagonal inertia, forced translation/torque, midpoint body
forces, a full intermediate-axis flip, Newton failures and quaternion continuity.
Exact SymPy checks verify the Cayley/momentum inverse relation and discrete
energy/work identities.

Actual tested environment: uv CPython 3.13.3, NumPy 2.4.4, SciPy 1.18.1 and
SymPy 1.14.0. A baseline is ideal mechanics, not an aerodynamic flight model or a
claim of reliable long-horizon phase at arbitrary step size. Native-language
fixture comparisons reconstruct the same body-to-inertial attitude convention.
