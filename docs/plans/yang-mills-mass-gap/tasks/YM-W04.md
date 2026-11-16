# YM-W04 — Read the Faizal–Shabir mass-gap claim and locate its crossover step

Tier T3 · Workstream W (watch) · Depends on: nothing · Review: T3 of a different family · Budget: half a day

## Why
arXiv:2606.19362 (M. Faizal, A. Shabir, 9 June 2026, hep-lat/hep-th), "Reflection-Positive Construction of a Four-Dimensional SU(N) Yang-Mills Theory with Mass Gap and Confinement", with four companion papers in Int. J. Geom. Meth. Mod. Phys. (2026) 2650111–2650114 included at its end, claims: a uniform transfer-operator gap and an area law from convergent strong-coupling character expansions; a multiscale analysis carrying them to the continuum; universality; and agreement with the weak-coupling continuum. If correct, it would settle what the program attacks. Clay listed the problem as unsolved on 2026-11-16. The strategy (`03-MASS-GAP-STRATEGY.md` §9) requires reading any such claim promptly and red-teaming it.

## Inputs
The public arXiv PDF (network read allowed). `03-MASS-GAP-STRATEGY.md` §2 (why the crossover from strong to weak coupling has no small parameter) and §5 route R3.

## Do
1. Fetch the PDF into `REPOS/yang-mills-work/watch/2606.19362/` (outside git). Record its SHA-256.
2. Write the theorem(s) actually stated, with hypotheses, in our notation. Separate stated theorems from prose claims.
3. Find the step that passes from the strong-coupling regime (small β, where cluster expansions converge) to the continuum limit (β → ∞ along asymptotic freedom). State exactly which inequality carries the gap across, what constants it depends on, and whether the paper proves it or cites it.
4. Check the claim of uniformity: is the gap bound uniform in β over the whole range, or only for β below a constant?
5. Note whether the paper uses or mentions the OpenAI construction.

## Deliverable
`notes/YM-W04-faizal-shabir.md`: five sections of protocol §5. Label every statement of the paper `source-claims`. If a gap in the argument is found, a `FINDING` with page and equation.

## Do not
Contact the authors or post anything outside the local network. Call the claim proved or refuted on the basis of the abstract. Spend more than the budget: if the crossover step cannot be located in half a day, say so.
