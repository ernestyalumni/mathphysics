# YM-O01 — Lean toolchain and the pinned compiled build on this machine

Tier T1–T2 · Workstream A · Depends on: **Ernest's approval** (toolchain install, ~2.2 GB download, ~15 GB disk on the Expansion drive, which is 97% full with 449 GB free) · Blocks: YM-A01, A02, T04 phase 2, all Lean work · Review: T1 · Budget: half a day wall time

## Do
1. Install `elan` with `ELAN_HOME` on the Expansion drive, not on the root NVMe. No system packages.
2. Work in a **copy** of `LEAN` at `REPOS/yang-mills-work/OpenAI-YM-20261028-build` (`git worktree` or `cp -a`); the pinned checkout stays untouched.
3. Follow `ARM/README.md` §"Compiled Lean files": `lake exe cache get`, download the two archive parts and `SHA256SUMS` from the release `lean-build-3456e5e`, verify the checksum, `lake unpack`, `lake build WilsonReference WilsonConstruction`.
4. Run the three verification commands from `LEAN/README.md` §Verification (`VerifyJoint.lean`, `verification/Usage.lean`, `verification/CompactCheck.lean`) and save logs.
5. Confirm an editor-free query works: `lake env lean` on a one-line file that `#check`s `YangMills.jointExistence` and `#print axioms` it.

## Deliverable
`PLAN/notes/YM-O01-lean-setup.md` with exact commands, versions, disk use, wall times; logs in `PLAN/ledger/evidence/YM-O01/`.

## Acceptance
`#print axioms YangMills.jointExistence` shows exactly `propext`, `Classical.choice`, `Quot.sound`. Checksums match. State clearly that this loads **precompiled** files and is not an independent rebuild.

## Do not
Download with `gh` if it would require logging in; use the `curl` route in the release notes. Do not run `gh auth login`.
