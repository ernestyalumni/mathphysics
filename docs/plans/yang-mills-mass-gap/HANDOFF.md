# HANDOFF — Yang–Mills program

Last updated: 2026-11-16 UTC, Codex follow-up. Shared checkout was master at 16e44a3; follow-up sources are on ym/YM-O02-codex-b7f4. One read-only research/review subagent ran at Ernest's request. The Mac daily source collector is scheduled; it does not launch a model.

## Session 3 (2026-11-15/16): local models — read `notes/local-models.md`

- **GPU containers on the desktop:** `nvidia-container-toolkit` was purged on 2026-11-11, so claw-work's launcher and every other GPU container fail. Ernest must reinstall it (sudo). **Workaround that works today, no sudo:** `bin/ym-local-llm start <model> <alias> <context>` passes the RTX 3060 into the llama.cpp container by hand and serves `http://127.0.0.1:8081/v1`.
- **Probe with known answers** of the three 9B profiles (`ledger/evidence/local-model-probe/probe.py`; results in the note, sections 4 and 4a). All three read a whole 20-page paper and answer factual questions correctly in about 20 seconds, catch one-line errors, and never flagged a correct proof. All miscalculate. Only Qwen3.8 distill Q5 sometimes caught a multi-step bookkeeping error.
- **Assignments:** local 9B models as *flagging* referees (union of the three; only INCORRECT flags count) and for verbatim indexing; Qwythos for tool use; nothing with calculation; no digests; no closing reviews. Section 6 of the note.
- **More capability:** a 27B-class model on the Mac (M4 Pro, 64 GB); Qwen3.8-27B per web search, unverified. Run the probe on it first.
- Protocol `06` §12: admission test before any local model gets a packet. YM-T14 extended with `referee` and `index` task types; that packet is now the way to put the local models to work.
- 2026-11-16: Ernest reinstalled the container toolkit; `launch.sh` works again. Tool-use and loop tests of Qwen3.8 distill Q5 (`agent_loop.py`, note section 4b): single tool calls and short plans 21/21, valid JSON always; long loops stop early with a false completion report about one run in three, so the harness must verify completion.
- The strategist's containers are stopped; Ernest's `llama-cpp-server` on port 8080 was left as he started it.

## Ernest's decisions of 2026-11-15 (all open questions from session 1 are closed)

1. PDD review card approved, with the amendment that agents may talk to other local agents, including on other machines of the local network (`05` review card, `06` §10).
2. Lean runs in Docker; installing on the Expansion drive is fine.
3. An agent with GitHub access may contribute to Hairer's public repository (`06` §11). Not possible today: the `gh` token on this machine is invalid.
4. Harnesses: Pi, Hermes Agent, Codex, Claude Code, Grok Build; local Qwen3.5 models when he is away (`06` §12, `prompts/launch.md`).
5. Direction: win with new tools (Rust symbolic and exact computation), do not rely completely on Lean (`03` §11, `05` §"Tools as the discovery loop").

## State

**Done in session 2**

- **Mathematical write-up** in the mono-notes: `documents/mono-notes/topics/yang-mills-lattice.tex` (Section "Lattice Yang–Mills theory: the scaling limit and the mass gap") and `topics/yang-mills-source-map.tex` (map from each labelled statement to a Lean declaration, file and line). Both PDF editions rebuilt (607 and 325 pages, no LaTeX warnings) and published in place.
  - Derived: C_𝔤 = 2N² for su(N); classical limit of the Wilson energy, β = 1/g₀²; slope b_* = 11C/(24π²) from the textbook coefficient with the normalization conversion; s⁻² decay from the u² spectral bound; gap ⇔ exponential decay; the massless one-particle rate u².
  - New, with proofs (unreviewed): purity bounds; the four-dimensional hypercubic defect and the gap bound on the doubled torus; a conditional lattice gap m ≥ log 2/n; the twist lemmas; the exact 2D Z₂ twist ratio.
  - Stated as conjectures: defect recurrence; small defect at a bounded physical scale; vanishing twist free energy.
- **Lean in Docker (YM-O01, in review).** Image `ym-lean:4.32.0`; working copy and compiled build under `REPOS/yang-mills-work/` (19 GB). Offline: `Build completed successfully (18126 jobs)`; `#print axioms YangMills.jointExistence` gives `[propext, Classical.choice, Quot.sound]`. This loads precompiled files; it is not an independent rebuild. Note: `notes/YM-O01-lean-setup.md`.
- **Rust, first tooling slice (in review).** 28 new tests, all passing, offline and locked:
  - `Rust/qft-symbolics/src/sphere.rs` (T01): exact sphere averages, sphere areas, shell coefficient 1/(8π²), slope 11C/(24π²).
  - `Rust/lie-groups/src/special_unitary.rs` (T02, su(N) only): Killing form and Casimir constant from the definition, exact integers.
  - `Rust/lattice-gauge` (new crate, no dependencies; T05 and the start of T10): torus geometry, Z₂ cochains with GF(2) elimination, exact Z₂ gauge density of states with twists, purity bounds. Its tests machine-check the two twist lemmas and the purity proposition of the notes.
- **PDD.** Both intents recorded by `pdd intent apply` in `docs/PRODUCT_INTENT.md` and `docs/intents/`. The provider-backed architecture workflow that follows is a generic 19-step generator; it reached step 4 in 25 minutes and was stopped by a time limit (first intent), and was stopped at step 2 after 4 minutes (amendment). **It generated no prompts, stories or code.** The Rust units were written by hand, prompt first, with prompts and `architecture.json` entries in each crate, as `Rust/qft-symbolics` has always been done. No regeneration from prompts was attempted. No user story was written for this program.
- **Coordination.** `SOURCES.md` and `bin/ym-sources` (pins verified OK for all five repositories); `bin/ym` works over SSH from other machines and locks writers; `bin/ym-lean`; protocol §10–§12; `prompts/launch.md`; new packet YM-T14.
- **Corrections to session 1.** `03` §2 had the expected mass as e^{−β̄/(2κ_G)}; the exponent is −β̄/κ_G. Only N1 of OpenAI's two-dimensional papers is formalized. Both fixed in place and recorded as claims C19, C20.

**Not done**

- The session-2 mathematical arguments and Rust tooling still await the reviews listed below. The later O02 launcher/transport and W03 collector have separate independent reviews recorded in notes/YM-O02-W03-independent-review.md.
- YM-O01 step 4 (the three verification scripts) not run. No clean rebuild (YM-A01), no Comparator (YM-A02).
- No SU(2) code, no Monte Carlo, no CUDA. T03, T04, T06–T11 not started.
- `bin/ym-micro` (runner for tool-less local models) not written: packet YM-T14.
- First Mac connection now reproduced through strict SSH and canonical tracker readback, independently checked. The original Pi model acknowledgement is not observed; no four-machine rollout is claimed.
- Nothing contributed upstream.
- The strategist has still **not** read: N1 beyond §2's start, N2 beyond §1 (its Eq. (1.6) and the claims around it are quoted from §1), N3, the unit outlines, or the classical lattice literature. The statement that the Wilson transfer operator is nonnegative in the form used by Prop. `ym:prop:defect` is a hypothesis in the notes, not something we established.

## Next actions, in order

1. **Reviews first.** Give YM-G02 (the new proofs) to a T3 model of a different family with a red-team instruction; give the three Rust units and YM-O01 to T2 reviewers. Tracker items 9–13.
2. **Wave 1** reading packets (items 1–8), strongest model on YM-R01, two families on YM-A03.
3. **Tools on the critical path:** YM-N01 (SU(2) heat bath on `Rust/lattice-gauge`) → YM-N02 (twist ratios); YM-T10 (extend exact enumeration: transfer matrix, Z_N) → YM-N03 (defect in solvable models: compute Δ_β(n) exactly for 2D and 3D Z₂ gauge theory and test the recurrence numerically); YM-T11 after R01 → YM-G03.
4. **Audit:** finish YM-O01, then YM-A01 in Docker (many hours of CPU; check disk first: 430 GB free at 97%).
5. YM-T14 so that the local 9B models have something safe to do.

## Things a successor must know

- The mono-notes are a **public** manuscript with a strict style (`documents/mono-notes/STYLE.md`): algebra in displayed equations, terse prose, no implementation names in the exposition. The Lean map is a separate section for that reason. Extend the topic; do not append corrections.
- `mathphysics` is a public repository. No local-network addresses, host names or credentials in it. They go in the untracked `MP/.pdd/ym-lan.env`.
- The local model server and GPU Monte Carlo share the RTX 3060.
- Do not rerun `pdd intent apply` expecting code: see "PDD" above. If Ernest wants the generator to run to the end, it needs a terminal session without a time limit and his say-so on provider usage.

## The strategic picture to preserve

- OpenAI's proof is a complete **ultraviolet** construction. It stops at unit scale with a weakly coupled lattice gauge system and reaches infinite volume by compactness alone.
- The mass gap is an **infrared** statement, uniform in volume. In the notes it is now reduced to two conjectures about one explicit quantity, the hypercubic defect Δ_β(n) = 16 log Z(n⁴) − log Z((2n)⁴): a recurrence, and smallness on a box of bounded physical size. Then three further steps: one coupling on both boxes, passage to the limit, uniqueness.
- OpenAI's two-dimensional papers show the shape of a full solution; only the first of them is formalized.
- Our route to the infrared input is the gauge analogue of N1 (centre twists). Its obstacles: frustration, loss of planarity, Wilson loops versus local observables.
- Our way of competing is computation before proof: exact small-lattice tests and Monte Carlo decide which mechanism is true before anyone tries to prove it. Lean certifies at the end.

## Log

- 2026-11-14 — Session 1. Program created: plan, strategy, protocol, 32 packets, tracker.
- 2026-11-15 — Session 2. Decisions received. Mono-notes topic with proofs and Lean map; Lean in Docker with the compiled build; first Rust slice (28 tests); LAN protocol; PDD intents recorded, generator not completed.
- 2026-11-15 — Session 3 started. Ernest asked which local models (claw-work llama.cpp profiles, or others) can take which packets. Working note: `notes/local-models.md` (written incrementally; read it first if this session was cut off).

- 2026-11-15 — Codex continuation requested by Ernest: first Mac Mini connection and complete reporting to Fable while its usage renews. Read [CODEX-HANDOFF-2026-11-15.md](CODEX-HANDOFF-2026-11-15.md) and [07-CONTINUATION-ROADMAP.md](07-CONTINUATION-ROADMAP.md). Existing local-model work preserved. Mac discovered and ping succeeds, but TCP/22 refuses; no authenticated LAN connection yet. Added O02/R09 packets, manual stories, Mac bootstrap prompt and portable tracker repair; no provider-backed regeneration.

## Codex follow-up, 16 November 2026 UTC

See `notes/YM-O02-hostname-trust-followup.md`. The Mac's host identity matches its agent's fingerprint; hostname trust is now present and strict SSH works through both address forms. Native Mac remote cases4/4 and Linux7/7 pass. Comment15 is a real Mac process readback of the strategist ANSWER, not original Pi model acknowledgement. O02 is in independent review, not automatically done.

Attribution correction: e478dc4 swept Codex's continuation plan/captures, Mac bootstrap prompt, O02/R09 packets, launcher repair/spec/tests and Codex handoff into Fable's local-model commit. Those artifacts were authored by Codex from Ernest's requests; Fable authored the session3 local-model/protocol material. The commit is left unchanged. Three original user-story sources are preserved on the follow-up feature branch; original master working-tree copies remain intact.

The shared master checkout was left on its existing branch; follow-up source work uses `ym/YM-O02-codex-b7f4`. Both mono-note editions already included the lattice/source-map topics. A dated front status/reading page has been added and both editions are being rebuilt. New daily watch sources/prompts are in `watch/`; the collector never calls an LLM, and queued AI review must not be reported as completed. Final deployment/verification outcomes follow below.

### Final follow-up verification

- O02 technical acceptance passed: strict hostname/IP SSH, native Mac remote cases4/4, Linux7/7, and live same-thread readback independently repeated by the reviewer. No original Pi model acknowledgement is asserted. Its next research suggestion remains R09(local:15), still unclaimed at the last check.
- Daily watch installed and crontab read back on the Mac: hourly local guard, once per Los Angeles date at/after09:00. Source collector only; queued AI-review prompt, no automatic model/provider usage. Initial live report fetched8/10 sources in18.09s, honestly partial: OpenAI announcement403, arXiv429. No auth bypass or repeated rate-limit retries. Second same-date invocation returned skipped-already-attempted. All11 watcher contracts pass on Linux and native Mac; independent code review accepts. Reports and cron backup stay private on Mac.
- Both reading PDFs published at the requested canonical mono-notes paths:608/327 pages, zero LaTeX citation/reference warnings. Existing overfull boxes47/3 unchanged. Rendered status pages inspected; status page10/6, full lattice topic530/287, source map539/292. screen.tex already inputs master.tex; no duplicate include added.
- The three original stories are saved unchanged with the feature work; the new daily-watch story is a manual acceptance source. No generated PDD contracts, provider-backed apply or whole-project regeneration claimed. Exact new request and read-only PDD plan are preserved.
