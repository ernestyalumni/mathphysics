# Verification demo -- latest run

- Timestamp (UTC): 2026-07-27T20:46:39.258824+00:00
- Git SHA: 946f2ffdeff51598f8c8c79de97e39c1ed8e3861
- Branch: master
- Backend: python 3.13.3, numpy 2.5.1
- Reference: arXiv:2602.12176 (Guevara, Lupsasca, Skinner, Strominger, Weil)

## Pass/fail criteria

- Residual tolerance: 1e-10
- A candidate PASSES a case if the max residual across trials is below tolerance.
- A candidate FAILS a case if any trial's residual is at or above tolerance.

## Results

| Candidate | n | trials | seed | max residual | status |
|---|---|---|---|---|---|
| correct closed form | 3 | 10 | 1003 | 0.00e+00 | PASS |
| correct closed form | 4 | 10 | 1004 | 0.00e+00 | PASS |
| correct closed form | 5 | 10 | 1005 | 0.00e+00 | PASS |
| correct closed form | 6 | 10 | 1006 | 0.00e+00 | PASS |
| correct closed form | 7 | 10 | 1007 | 0.00e+00 | PASS |
| correct closed form | 8 | 3 | 1008 | 0.00e+00 | PASS |
| correct closed form | 9 | 3 | 1009 | 0.00e+00 | PASS |
| correct closed form | 10 | 3 | 1010 | 0.00e+00 | PASS |
| broken closed form (2^(n-1) bug) | 4 | 10 | 2004 | 5.00e-01 | FAIL |
| broken closed form (2^(n-1) bug) | 5 | 10 | 2005 | 5.00e-01 | FAIL |
| broken closed form (2^(n-1) bug) | 6 | 10 | 2006 | 5.00e-01 | FAIL |

## Summary

- Correct candidate: PASS on all n
- Broken candidate (deliberate off-by-one exponent bug): correctly caught (FAIL) on all cases

**LABEL: verified reproduction of a known result; no novelty claim.**
