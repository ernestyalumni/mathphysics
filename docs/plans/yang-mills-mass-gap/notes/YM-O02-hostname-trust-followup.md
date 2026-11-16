# YM-O02 follow-up: hostname trust and native Mac round trip

## Result

The hostname/IP discrepancy is resolved. The desktop had an Ed25519 trust entry for the Mac IP, but none for its local hostname. The Mac-side fingerprint in tracker comment12 matched both live presentations, the existing IP trust entry and the public host-key file read through strict authenticated SSH. The hostname aliases were added after a backup; strict SSH through IP and hostname now succeeds. A Mac process has read Fable ANSWER(comment13) and posted receipt(comment15) without impersonating the Pi model.

## Labels

Connection and machine readback: reproduced by named commands. Original Pi-agent acknowledgement: not observed. Original launcher repair and this setup await independent review; author is not its own reviewer. No mathematical result is claimed.

## Evidence

- Linux `python3 docs/plans/yang-mills-mass-gap/tests/test_ym_launcher.py`:7/7 pass.
- Native Mac Bash3.2/Python3.9.6 remote-worker cases:4/4 pass (remote routing, trusted env routing, rejected remote new, literal metacharacters). The Mac does not need local PDD/flock.
- Launcher SHA256 c77230e560661f57970fcb695d0442d2b0adee7a006bfaf5d48dc918f3b0cb50; test SHA256 3bbe545ed738daeb429c34d5f2dbe01241525c47b6280664c38441b95bb6ce52 match on both machines.
- Strict SSH options: `-T -o BatchMode=yes -o StrictHostKeyChecking=yes -o HostKeyAlgorithms=ssh-ed25519 -o ConnectTimeout=6`; identity command: `hostname; id -un; uname -s; ssh-keygen -lf /etc/ssh/ssh_host_ed25519_key.pub`. Machine values, exact fingerprints and known_hosts backup location stay in the private connection receipt, not this public file.
- Mac executes its own launcher `show 14 </dev/null`, confirms Fable's session identifier and ANSWER, then posts as codex-mac-probe-b7f4. stdin is closed so the nested SSH command cannot consume later lines of the remote orchestration script.
- Tracker comments14–16 record the claim, machine readback and resolution.

## Not done / not checked

No original Pi model acknowledgement fabricated; no host-key bypass, private-key read or firewall/power change. This is one Mac, not a four-machine deployment. Reviewer verdict is a separate tracker entry.

## Questions for the strategist

Commit e478dc4 included the Codex-authored continuation roadmap, original request/plan captures, Mac prompt, O02/R09 packets, launcher/source specification/tests and Codex handoff along with Fable's session3 local-model work. Its message credits Fable alone. Correct attribution is Ernest's intent, Codex's listed continuation artifacts and Fable's local-model/protocol changes. This note corrects the record without rewriting e478dc4. The three original story sources are preserved unchanged in this feature worktree; the canonical master checkout's original untracked copies are left intact.
