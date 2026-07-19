# Research Plan: Verified MHV Amplitudes as the First Discovery-Engine Wedge

**Last updated:** 2026-07-19
**Owner:** Ernest Yeung + TMPAI
**Active branch:** `feat/yc-mhv-execution-2026-07`

## Outcome

Build a theoretical-physics research agent whose claims are checked by independent numerical or symbolic verifiers. The first narrow wedge is single-minus gluon amplitudes in the half-collinear region studied in arXiv:2602.12176v2 (SMGA).

The immediate objective is not to claim a new result. It is to make the existing reproduction legible and runnable in under two minutes, then use that artifact as the foundation for more ambitious discovery loops.

## Verified baseline on 2026-07-19

- `amplitudes/scripts/08_smga_stripped_amplitudes.py` runs successfully and checks explicit low-multiplicity formulas, discrete structure, the soft theorem, and cyclicity.
- The focused MHV suite passes: **43 tests passed** across the amplitude unit tests and the SMGA integration test.
- The integration suite compares two separately implemented paths—the specialized SMGA Berends–Giele recurrence and the closed form—for multiplicities through `n = 10`.
- The complete package suite currently reports **62 passed and 2 failed**. Both failures are in pre-existing, uncommitted `massless_spinors.py` work and are not on the MHV demo path.
- The default `uv run pytest` entry point is not yet clean because package metadata/dependencies need hardening. See `tasks/TASK-000-baseline-reproducibility.md`.

These are reproduction and verification results. They are not evidence of original physics yet.

## July 19–25 critical path

1. **Reproducible entry point** — one documented command from a fresh checkout.
2. **Two-minute demo** — show a candidate formula, an independent recurrence, verifier output, and an honest result ledger.
3. **Application-facing explanation** — state the product as verifier-grounded scientific work, not “a chatbot for physics.”
4. **Record founder and demo videos** — founder video is a separate, founder-only one-minute recording; the demo may show the product.
5. **Submit by July 25** — July 26–27 are contingency, not planned workdays.

## Research sequence after the application

1. Generate a permanent Eq. (16) verification report for `n = 3..10` with seeds and run metadata.
2. Add field-theory consistency checks: cyclicity, reflection, U(1) decoupling, selected KK relations, soft and collinear behavior.
3. Stabilize a machine-readable Cadabra2/SymPy verifier interface with known-true and known-false tests.
4. Explore regions adjacent to `R_1` with a proposer → critic → verifier → logger loop.
5. Only after fresh-point verification and literature review, promote a surviving conjecture into a paper section.

## Deliberate non-goals before July 25

- No NMHV or graviton expansion.
- No general-purpose physics-agent platform rewrite.
- No elaborate dashboard if a CLI plus static report communicates the loop clearly.
- No arXiv novelty claim.
- No cleanup of unrelated Srednicki build artifacts or the uncommitted massless-spinor work.

## Source-of-truth files

- Agent-sized work: `tasks/README.md`
- Current execution queue: `TASKS.md`
- Demo handoff: `AGENT-HANDOFF-2026-07-19.md`
- Product charter: `TMP-MISSION.md`
- Loop and eval rules: `agent/DISCOVERY-LOOPS.md`, `agent/EVAL-AND-TOOLING.md`
- Public research material: `amplitudes/`

## Truth standard

Every public statement must be one of:

- **reproduced** — matches a cited known result;
- **verified** — a named verifier returned a recorded result;
- **conjectured** — proposed but not proved;
- **inconclusive** — the verifier could not decide.

“Discovered” is reserved for a result that survives independent verification and a real literature search.
