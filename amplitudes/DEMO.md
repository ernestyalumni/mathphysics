# Two-minute verification demo

Shows the core product behavior: a candidate physics formula is checked against an
independently implemented verifier, and the run ends in an explicit PASS/FAIL ledger --
including a deliberately wrong candidate, so the demo shows the verifier actually rejecting
something, not just displaying green checks.

## Run it

From the repo root, in a fresh checkout:

```bash
uv run --project amplitudes/MHVamplitudes python amplitudes/scripts/09_verification_demo.py
```

No manual environment setup, no Docker, no Cadabra2. Total runtime is well under a second on
a development Mac; interactive/reading time for a reviewer is under two minutes.

## What it does

1. Prints the candidate formula (SMGA closed form, arXiv:2602.12176 eq. 16) and its
   provenance, and names the independent verifier (Berends-Giele recursion, same paper).
2. Compares the two independently implemented code paths across `n = 3..10` on random,
   momentum-conserving kinematic configurations (seeded, reproducible) -- the known-correct
   case. All cases PASS with residual `0.00e+00`.
3. Runs a second, deliberately wrong candidate (a plausible-looking bug: the normalization
   exponent is `2**(n-1)` instead of the correct `2**(n-2)`, built by recomposing the same
   real `sg_ij`/`sg_i_set` primitives incorrectly -- not a toy failure injected outside the
   physics). All cases correctly FAIL.
4. Writes a structured report to [`amplitudes/results/demo-latest.md`](results/demo-latest.md)
   with timestamp, git SHA, branch, backend versions, seeds, trial counts, and a
   residual-summary table.
5. Exits 0 only if the correct candidate passes everywhere and the broken candidate is caught
   everywhere -- i.e. the script itself fails loudly if the demo stops demonstrating what it
   claims to.
6. Ends with the explicit label: **verified reproduction of a known result; no novelty
   claim.**

## What a reviewer should look for

- The correct-candidate table: all residuals exactly `0.00e+00` through `n = 10`.
- The broken-candidate table: nonzero residuals, all marked `FAIL` -- proof the verifier
  rejects a wrong answer rather than rubber-stamping anything plausible-looking.
- The report file (`amplitudes/results/demo-latest.md`) is regenerated on every run and is
  not hand-edited -- it's evidence, not marketing copy.

## Truth label

This reproduces a known result from arXiv:2602.12176 through `n = 10`. It is a
reproducibility/engineering milestone, not a new physics result, and not evidence beyond the
paper's own all-`n` formula and proof.
