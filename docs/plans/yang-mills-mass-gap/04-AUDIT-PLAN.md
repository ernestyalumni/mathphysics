# Audit plan: prove or disprove parts of the construction

A Lean kernel has accepted the proof terms. So "disprove" cannot mean finding a wrong step inside a proof. It can mean any of the following, and each is a real way for the result to be less than advertised.

| # | Attack surface | What failure would look like | Who has covered it | Our task |
|---|---|---|---|---|
| 1 | The build itself | Sources do not rebuild from clean; a hash mismatch; a hidden `sorry` or axiom | Group: clean rebuild on Apple Silicon, axiom audit, 100,168 declarations replayed | YM-A01: independent clean rebuild on Linux x86-64 |
| 2 | Process isolation | Build code that tampers with the check | **Not covered.** Group ran Comparator through `fake-landrun.sh` on macOS: "do not test isolation of potentially hostile build code" | YM-A02: Comparator on Linux with real Landrun |
| 3 | Kernel soundness | A Lean kernel bug accepts a wrong term | **Not covered.** "An additional external kernel implementation was not used" | YM-A02b: replay the export in a second kernel if one supports this toolchain |
| 4 | The statement is weaker than it reads | A cheap object satisfies `JointExistence` | Partly: OpenAI itself said the public target omits lattice connections; the group wrote a stronger target | YM-A03: red-team both targets for degenerate models |
| 5 | Definitions do not mean what the names suggest | "Wightman axioms" that are not Wightman's; a "Wilson measure" that is not Wilson's | Group: human review of a glossary, bridge to independent lattice and OS formalizations (A1, A3) | YM-A04: independent transcription of the public statement, diffed against Armstrong's |
| 6 | The informal accounts are wrong | Armstrong's or Hairer's prose claims something the Lean does not say | Their own review cycles | YM-A05: spot-check units we depend on (III.9, III.12, IV.1, V.1, V.3) |
| 7 | Numbers disagree with physics | A coefficient that is not the textbook one | Group: leading trace coefficient (A4) | YM-T01–T03: recompute κ_G, C_𝔤, the F² two-point constant in Rust |
| 8 | Vacuity of hypotheses | A theorem quantified over an empty or trivial class | Unknown | Part of YM-A03 |

## What a finding is

A finding is a sentence of the form "*declaration X at path:line, read with its definitions unfolded to depth k, says S; the document D at location ℓ says S′; S ≠ S′ because …*" plus the commands that show it. Anything vaguer is a question, not a finding.

Findings go to `ledger/findings/YYYY-MM-DD-<slug>.md`. Before anything is reported upstream (issue on Armstrong's or Hairer's repo, mail to the group) it needs a review by a different model family and Ernest's explicit go-ahead. Upstream reports are external actions.

## Red-team method for attack surface 4

For each clause of `JointExistence` (use the clause names of `statement.md` §3: Reg, P1–P4, W0–W3, R1–R15, O1–O2, T1–T4):

1. State what the clause is supposed to exclude.
2. Try to satisfy the clause, together with all others, by a construction that is not Yang–Mills: zero insertions, constant insertions, a free field, a generalized free field with continuous mass, a theory at a different slope b_*, insertions that ignore the gauge field, a phase space X with one point and a trivial theory.
3. Record which clause blocks each attempt and how. If nothing blocks it, that is a finding.

The published non-claims already concede some of this (generalized free fields are not excluded; b_* is only positive in the public theorem). The task is to find what is **not** already conceded.

## What we do not audit

- The internals of the 70% of the proof that is the ultraviolet renormalization group, line by line. That is Hairer's project; we contribute there if we choose to (task type `upstream`), we do not duplicate it.
- Mathlib.

## Contributing upstream

Hairer's repository takes agent contributions under a self-claim protocol on GitHub issues (`hairer/humanise-yang-mills/docs/coordination.md`). Units we digest for our own purposes (IV.1, V.1, V.3 and the unit-scale interface) are natural contributions. This needs Ernest's GitHub account and his decision; no agent posts there on its own.
