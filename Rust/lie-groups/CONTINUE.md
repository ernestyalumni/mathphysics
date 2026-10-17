# Continuation — Frankel rotation-group toolkit

The requested bounded implementation is complete on feat/lie-groups-so3-su2.
Read README.md, docs/SOURCE-AUDIT.md and docs/VERIFICATION.md first.

Preserve conventions:
- Active right-handed column-vector rotations; epsilon012=+1.
- E_i=hat(e_i), T_i=-i sigma_i/2, F_i=-i sigma_i; T maps to E, F maps to 2E.
- SU2 is the actual unit-quaternion group, not quaternions modulo sign.
- Distinct mutable general matrices versus private, validated group values.
- Exact symbolic residuals, numerical tolerance checks and finite witnesses
  are different evidence classes.

Optional future units (not completed and no background job running):
1. Arbitrary-precision Rust coefficient backend and uniform fallible arithmetic;
   coordinate with ../qft-symbolics/CONTINUE.md rather than introduce two CASes.
2. Explicit logarithm/branch APIs (pi rotations and antipodal spinors need a
   documented branch choice); additional independent stress/property cases.
3. Multiple symbolic angles or an expression-DAG layer, only if physics tasks
   need it. Current Trig has one common angle and retains half-angle metadata.
4. Higher-spin/Wigner representations with a documented Cartesian/spherical
   basis change; do not conflate j=1 matrices with the standard Cartesian basis.
5. If required, introduce a separate Euclidean3 index space when integrating
   with the typed QFT tensor kernel. No such integration is presently claimed.
6. Optional native half backend and performance measurements; existing Half
   is binary16 storage with promoted software operations, not hardware support.
7. Portable CI/container recipes using locally authorized dependencies/images.

Run bash validation/run.sh, keeping artifacts external, before new claims.
Never weaken determinant/domain checks merely because a CAS class name says SO.
No new physics, whole-CAS parity or arbitrary-input agent verifier is claimed.
