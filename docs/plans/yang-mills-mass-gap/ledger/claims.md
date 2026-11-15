# Claims ledger

Every claim the program relies on. One row each. Labels as in `README.md` §"Truth standard". Add rows; never silently change a label: strike the old one and date the change.

| # | Claim | Label | Evidence | Date |
|---|---|---|---|---|
| C1 | `YangMills.jointExistence` proves `JointExistence`; 100,168 project declarations replayed; axioms propext, Classical.choice, Quot.sound only | source-claims (kernel-checked by the group, not by us) | `LEAN/README.md` §Verification; assessment §9.1 | 2026-11-14 |
| C2 | The construction does not prove a mass gap; low-energy bounds are polynomial (u², s⁻²) | source-claims | assessment p.1, §2 item 4; `ARM/statement.md` §4 item 5; `LEAN/YM/EveryPhase/Decay.lean:173–192` (pointer from provenance, not opened by us) | 2026-11-14 |
| C3 | The infinite-volume limit is subsequential; phases not shown equal | source-claims | assessment §5, §6.2 | 2026-11-14 |
| C4 | Angular average of h(k)\|k\|⁴ over S³ is −(11/3)C; shell integral of \|k\|⁻⁴ is 2π² log L; κ_G = 11C/(24π²) | reproduced | `ledger/evidence/check_oneloop.py` (SymPy), run 2026-11-14 | 2026-11-14 |
| C5 | C_𝔤 = 8 for su(2), 18 for su(3) with ⟨X,Y⟩ = −(1/N) Re tr(XY); Killing form = −C × that form | reproduced | same script | 2026-11-14 |
| C6 | The RG stops at unit scale with a fixed large coupling β̄; partition function preserved exactly at every level | source-claims | `ARM/summary.md` §"The strategy", §"The block-spin step" | 2026-11-14 |
| C7 | The group's Comparator runs lacked Linux process isolation and used no second kernel | source-claims | assessment §9.2 | 2026-11-14 |
| C8 | OpenAI's N1–N3 prove 2D O(n≥3) exponential decay at all β, an O(4) full gap, and a unique massive O(3) continuum limit | source-claims (their preprints; formalization status per `openai/math/lean/docs/215.md`, not opened by us) | `openai/math/CONTENTS.md` family 215; N1 abstract and §1.3; N2 Thm 1.1, Cor 1.2 | 2026-11-14 |
| C9 | For positive trace-class T: (Tr Tⁿ)²/Tr T²ⁿ − 1 ≥ Σ_{i≥1}(λ_i/λ_0)ⁿ | conjectured (three-line derivation in `03` §4; to be checked by YM-G02) | — | 2026-11-14 |
| C10 | A mass gap for the continuum theory reduces to an infrared estimate for the unit-scale lattice system, uniform in volume | conjectured (strategist's reading; M0–M2 make it precise) | `03-MASS-GAP-STRATEGY.md` §2, §4 | 2026-11-14 |
| C11 | The Clay prize is for existence and mass gap together; publication, two-year wait, general acceptance required | source-claims (web search; primary text pending YM-W02) | `02-MILLENNIUM-PRIZE-STATUS.md` sources | 2026-11-14 |
