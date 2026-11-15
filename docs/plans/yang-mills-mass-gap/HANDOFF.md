# HANDOFF — Yang–Mills program

Last updated: 2026-11-14, end of session 1 (strategist: Fable 5.2). Branch `feat/yang-mills-mass-gap-program`, not pushed, not merged. Nothing is running in the background. No agents have been launched.

## State

**Done this session**

- Read at overview depth: the review group's assessment, Armstrong's statement / summary / provenance / FAQ, Hairer's coordination rules, the Lean repo's README, and the abstract and strategy sections of OpenAI's N1 and N2.
- Written and committed: charter (`README.md`), digest (`01`), prize answer (`02`), mass-gap strategy (`03`), audit plan (`04`), tooling roadmap with PDD review card (`05`), agent protocol (`06`), board, 32 task packets, four launch prompts, claims ledger, tracker wrapper `bin/ym`.
- Tracker: wave-1 packets registered as local work items 1–8 (`bin/ym list`). The tracker lives in `MP/.pdd/work_items/` and is not versioned.
- Reproduced with SymPy: the one-loop arithmetic (κ_G = 11C/(24π²)) and C_𝔤 = 2N² for su(2), su(3). Script: `ledger/evidence/check_oneloop.py`.
- `pdd intent plan` run (read-only). Intent ID `yang-mills-construction-audit-rust-lattice-gauge-d3c6da77`.

**Not done**

- `pdd intent apply`: not run. Needs Ernest to approve the review card in `05-TOOLING-ROADMAP.md`. No Rust code written.
- Lean toolchain: not installed. YM-O01 needs approval (disk).
- No packet has been executed. No deliverable exists under `notes/`.
- The strategist has **not** read: any Lean file beyond the top-level README; the unit outlines; N1 beyond §1 and the start of §2; N2 beyond §1; N3 at all; any of the classical literature (all citations in `03` §5 and §10 are from memory and flagged as such).

## Waiting on Ernest

1. Approve (or correct) the PDD review card in `05-TOOLING-ROADMAP.md` → unblocks all `T` and `N` packets.
2. Say whether to install Lean and the 8.8 GB compiled archive on the Expansion drive → unblocks YM-O01, A01, A02.
3. Say whether agents may fetch public papers (arXiv) and `git fetch` the upstream repos → unblocks R06, R07, W01, W02.
4. Say which harnesses and models to use for T1–T3, and whether the strategist should launch wave 1 itself (it can spawn sub-agents on request) or he will start them with `prompts/worker.md`.
5. Later: whether to contribute to Hairer's public repository under his GitHub account.

## Next actions, in order

1. Launch wave 1 (`BOARD.md` §"Suggested first launch"): eight packets, all reading and analysis, no approvals needed. Strongest model on YM-R01.
2. Coordinator pass; assign cross-family reviewers.
3. When R01 is accepted: open YM-G03 (twice, independently) and YM-G04. These decide the primary route.
4. When R02, R03, R05, R08 are accepted: strategist writes YM-G01, the problem sheet (milestone G0).
5. In parallel, once approved: YM-O01 → A01 → A02; YM-T01, T02, T05, T04.

## The strategic picture to preserve

- OpenAI's proof is a complete **ultraviolet** construction. It stops at unit scale with a weakly coupled lattice gauge system and reaches infinite volume by compactness alone.
- The mass gap is an **infrared** statement about that unit-scale system, uniform in volume. It is the classic open problem of lattice gauge theory at weak coupling, now in a precise setting.
- OpenAI's own two-dimensional papers (N1–N3) show the shape of a full solution: one new infrared idea (N1: fourth variation in two non-commuting directions ⇒ twist ratio → 1 ⇒ exponential decay) plus transport by the renormalization group. The Yang–Mills Lean has the transport machinery and none of the infrared idea.
- Our primary route is the gauge analogue of N1 (centre twists, vortex containers, embedded Z₂ gauge theory). Its known obstacles are frustration (Z₂ monopoles), loss of planarity, and the step from Wilson loops to local observables.
- Independent of whether that works, the conditional theorem "infrared estimate ⇒ gap and uniqueness", the Linux isolation audit, and the Rust tools are each worth having.

## Log

- 2026-11-14 — Session 1. Program created. See "State".
