# YM-A03 — Red-team the formal statement for degenerate models

Tier T3 · Workstream A · Depends on: none to start (reading); YM-O01 to test anything in Lean · Blocks: G1 · Review: T3 different family, plus strategist · Budget: two sessions. **Run twice independently** by two model families; do not read the other thread before both `SUBMIT`.

## Goal
Find out what the public theorem `JointExistence` and the group's `combined_acceptance` do **not** pin down, beyond the published non-claims.

## Inputs
`LEAN/ChallengeSpec.lean`, `PublicWilson.lean`, `PublicTheory.lean`, `PublicView.lean`, `Construction.lean`; `ARM/statement.md` §3–§4; `YMG/submissions/openai-20261028/checks/` for the group's target; method in `PLAN/04-AUDIT-PLAN.md` §"Red-team method".

## Do
1. Clause table: every conjunct of `JointExistence`, its clause name from `statement.md`, what it is meant to exclude.
2. For each candidate cheap model in the method section, walk the clauses and record the first one that blocks it, quoting the Lean.
3. Try to weaken: is there a clause whose removal would not change which of your candidate models are excluded? That clause is doing less than it appears.
4. Check quantifier traps: an existential that can be satisfied by an empty type, a "for all test functions" over a space with no elements for some n, a limit along a filter that is trivially eventually constant, insertions allowed to depend on (N, p) in a way that makes "lattice origin" vacuous.
5. Same exercise, shorter, for `combined_acceptance`: which of your survivors does it kill?

## Deliverable
`PLAN/notes/YM-A03-spec-redteam-<session>.md`. Anything that survives all clauses goes also to `PLAN/ledger/findings/`.

## Acceptance
Every candidate has a blocking clause with `path:line`, or is a finding. Already-conceded gaps (`PLAN/01-WHAT-OPENAI-PROVED.md` §2) are listed separately and not counted as findings.
