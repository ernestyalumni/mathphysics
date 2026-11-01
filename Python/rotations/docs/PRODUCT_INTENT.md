# Product Intent (PRD)

This current product-intent record is maintained by `pdd intent apply`.
Each accepted change links to an immutable intent event.

<!-- pdd-intent-entry:rotation-library-conventions-as-pdd-source-22448aa3:start -->
## Rotation library conventions as PDD source

- Intent event: [`docs/intents/intent__rotation-library-conventions-as-pdd-source-22448aa3.md`](intents/intent__rotation-library-conventions-as-pdd-source-22448aa3.md)
- Change kind: `add`
- Supersedes: none
- Scope: `existing_subproject_adoption`
- Technology: `python`

> As a researcher who uses the rotation library in Python/rotations to check Goldstein, Frankel and engineering attitude conventions, I want the library's conventions written down as the source it is regenerated from, so that a regenerated version gives the same rotations, angles and quaternions as the current one.
>
> The library MUST keep right-handed column vectors with x_I = R x_B, and Goldstein's coordinate matrix MUST be the transpose A_G = R.T.
> Euler angles MUST be given and returned in intrinsic product order: ZXZ and ZYZ as (phi, theta, psi), and ZYX as (heading, pitch, roll) meaning Rz(heading) Ry(pitch) Rx(roll).
> Extraction MUST return the middle angle in [0, pi] for ZXZ/ZYZ and in [-pi/2, pi/2] for ZYX; at gimbal lock it MUST set the last angle to zero and report singular=True.
> Quaternions MUST be scalar-first Hamilton quaternions, and matrix-to-quaternion extraction MUST return the representative with w >= 0.
> Axis-angle extraction MUST return an angle in [0, pi] and MUST use the axis (1, 0, 0) for the identity.
> Cayley-Klein parameters MUST satisfy gamma = -conj(beta) and delta = conj(alpha), with alpha = w + i z and beta = y + i x.
> The library MUST reject non-finite values, reflections, non-unit quaternions, non-unit axes and CK parameters that break the conjugacy constraints.
> The library MUST NOT silently project accepted data onto SO(3); the only normalization allowed is the documented one inside matrix-to-quaternion extraction.
> The rotation functions MUST NOT depend on anything outside the Python standard library.
> Never return a ZYX tuple in (roll, pitch, yaw) order.
>
> For example, from_euler('zyx', 0.4, 0.7, -0.3) followed by to_euler(R, 'zyx') returns (0.4, 0.7, -0.3) with singular=False.
> For example, a 90 degree turn about z maps the body x axis to the inertial y axis, and to_quaternion returns (cos(pi/4), 0, 0, sin(pi/4)).
<!-- pdd-intent-entry:rotation-library-conventions-as-pdd-source-22448aa3:end -->
