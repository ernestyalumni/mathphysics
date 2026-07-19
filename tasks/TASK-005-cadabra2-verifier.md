# TASK-005 — Cadabra2 verifier wrapper + known-answer tests

Status: TODO
Priority: P1
Depends on: — (parallel to the numeric track)
Estimated size: one session (3–5 h)
Branch suggestion: `feat/task-005-cadabra2-verifier`

## Goal

A thin, uniform wrapper so agent loops can submit a symbolic claim to Cadabra2 and get back
a machine-readable verdict. This is pillar-2 infrastructure from `../agent/RESEARCH-AGENDA.md`
(thread 2): the discovery loop is only credible if the verifier is boring and reliable.

## Deliverables

1. `Cadabra2/_verifier.py` — `verify(claim_script: path) -> {status: pass|fail|inconclusive,
   residual, notes, runtime_s}`. A claim script is a Cadabra2/Python file that computes a
   difference expression and reports whether it simplifies to zero.
2. `Cadabra2/_test/` — at least 3 known-answer tests:
   - a spinor identity from Srednicki Ch. 48–50 already computed in
     `amplitudes/Srednicki/individual_chapters/`,
   - a Fierz/color identity (port from `amplitudes/scripts/01_color_fierz_identity.py`),
   - one deliberately-wrong claim that must return `fail` (guards against a verifier that
     rubber-stamps everything).
3. `Cadabra2/README.md` — how to run (native install or the `cadabra2-ubuntu:24.04` Docker
   image per `../RESEARCH-PLAN.md`), plus the template ported from
   `Monoclaw/Python/Cadabra2/Srednicki/` (`mhv_parke_taylor.py` + `mhv_export_latex.py` pair).

## Definition of done

A single command (document it, e.g. `python Cadabra2/_verifier.py --run-tests`) runs all
known-answer tests: the true claims return `pass`, the wrong claim returns `fail`, and the
output is valid JSON. If Cadabra2 cannot be installed in this environment, document the
exact blocker in `Cadabra2/README.md` and implement the same wrapper interface backed by
SymPy so downstream tasks are unblocked — but say clearly which backend ran.
