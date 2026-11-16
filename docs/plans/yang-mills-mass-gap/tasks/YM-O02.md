# YM-O02 — First Mac Mini connection and Fable continuity

Tier T2 · Operations · Current priority requested by Ernest · Review: independent T2.

## Do

1. Read `07-CONTINUATION-ROADMAP.md`, protocol §10 and the Mac bootstrap prompt.
2. Keep machine-specific addresses/usernames/fingerprints in private workspace notes or untracked `.pdd/ym-lan.env` on the worker. Use the Linux canonical host for PDD/tracker writes.
3. Establish discovery → ping → TCP/22 → verified host identity → key authentication → same tracker read. Either direction can bootstrap; Mac-to-Linux is enough for the canonical tracker. Linux-to-Mac needs Remote Login enabled on the Mac.
4. Characterize and repair `bin/ym` portability without changing its tracker interface. No local Linux tools required on a remote Mac worker.
5. Once SSH works, clone the existing program branch over the LAN into a separate directory; preserve existing clones. Never clone or upload home directories or private material. Copy any uncommitted onboarding files explicitly and mark them as a local overlay, or wait for a coordinator-approved commit; a clone transfers committed content only.
6. Post PROGRESS hello and ANSWER on this thread from distinct machines. Record receipt time and stages honestly. Independent review is still required before `done`.
7. Update `CODEX-HANDOFF-2026-11-15.md` for Fable. Do not inject terminal keystrokes into an unidentified live Claude session or launch a billed model to deliver a report.

## Acceptance

Launcher regressions pass; real authenticated SSH identity matches the local console fingerprint; Mac reads an existing canonical thread and exchanges a reply. Report uncompleted stages. Router changes, port forwarding, opening model APIs and wide network scans are outside this packet.

## Deliverable

Public, sanitized handoff plus private connection receipt. Public plan/manuscript files contain no LAN identifiers. Once connected, R09 is a suitable first context task; mathematical red-team tasks require a strong reviewer.
