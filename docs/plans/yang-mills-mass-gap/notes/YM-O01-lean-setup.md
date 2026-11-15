# YM-O01 — Lean in Docker with the pinned compiled build

Session: `claude-fable52-strat` (strategist), 2026-11-15. Status: submitted, **not yet reviewed**.

## 1. Result

- The pinned proof loads in Lean 4.32.0 inside the Docker image `ym-lean:4.32.0`, with networking disabled, on this Linux x86-64 machine.
- `lake build WilsonReference WilsonConstruction` reports `Build completed successfully (18126 jobs)`.
- `#print axioms YangMills.jointExistence` prints `[propext, Classical.choice, Quot.sound]`.
- The compiled archive's SHA-256 equals the published one.
- Nothing was installed on the host. Disk use: 19 GB under `REPOS/yang-mills-work/` (17 GB working copy with Mathlib cache and project build, 2.1 GB downloads) plus the image.

## 2. Labels

| Statement | Label |
|---|---|
| The axioms line above was printed by Lean 4.32.0 here | **reproduced** (from precompiled files) |
| The proof is correct | **source-claims** for us. This run loaded `.olean` files compiled by someone else. It is **not** an independent rebuild and not a kernel re-check of the proof terms. That is YM-A01. |
| Checksum match | **reproduced** |

## 3. Evidence

Receipt: `ledger/evidence/YM-O01/receipt.txt`. Full logs: `REPOS/yang-mills-work/logs/{download,cache-get,unpack-build}.log` (not in git).

    PLAN/bin/ym-lean build                       # image from PLAN/docker/lean/Dockerfile
    PLAN/bin/ym-lean prepare                     # clone of LEAN at 3456e5e into REPOS/yang-mills-work/OpenAI-YM-20261028-build
    PLAN/bin/ym-lean run lake exe cache get      # Mathlib cache, 8639 files (network)
    curl -fL -O https://github.com/scottnarmstrong/yang-mills-informalization/releases/download/lean-build-3456e5e/{SHA256SUMS,ym-build-3456e5e.tgz.part0,ym-build-3456e5e.tgz.part1}
    cat ym-build-3456e5e.tgz.part0 ym-build-3456e5e.tgz.part1 > ym-build-3456e5e.tgz && sha256sum ym-build-3456e5e.tgz
    PLAN/bin/ym-lean offline lake unpack ym-build-3456e5e.tgz
    PLAN/bin/ym-lean offline lake build WilsonReference WilsonConstruction
    PLAN/bin/ym-lean offline lake env lean YmAxiomCheck.lean     # import Solution; #check; #print axioms

Wall time for checksum, unpack, build and axiom print together: 44 minutes (15:07–15:51 UTC).

## 4. Not done / not checked

- Whether `lake build` recompiled any module or only validated traces. The 44 minutes were not split by step. A reviewer should rerun `lake build` (it should now be immediate) and time `lake unpack` separately.
- The three verification scripts of `LEAN/README.md` §Verification (`VerifyJoint.lean`, `verification/Usage.lean`, `verification/CompactCheck.lean`): **not run**. Packet step 4 is open.
- `YmAxiomCheck.lean` was written into the working copy (untracked there). It is three lines; see the commands above.
- The Dockerfile fetches `elan-init.sh` from GitHub at build time without pinning a hash. For an audit-grade image, pin it.
- No independent rebuild (YM-A01), no Comparator under Landrun (YM-A02).

## 5. Questions for the strategist

None.

## Resume here

Run packet step 4 (the three verification scripts) with `PLAN/bin/ym-lean offline lake env lean <file>`, save the logs under `ledger/evidence/YM-O01/`, then hand to a reviewer. After that YM-A01 can start: a clean clone, no `lake unpack`, full `lake build`; expect many hours of CPU (the group reports 184 minutes on Apple Silicon).
