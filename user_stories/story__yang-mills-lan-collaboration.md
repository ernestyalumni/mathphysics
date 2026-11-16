# Collaborate on Yang–Mills across the local network

## Story

As Ernest running research agents on four local computers, I want the first Mac Mini and then the other machines to share task claims, questions, results, reviews and restart notes with the canonical Linux checkout, so that work continues while a particular model's usage is unavailable.

Source: continuation outcome f and `request-2026-11-15-fable-reporting.txt`; approved LAN amendment `fb139368`. Manually authored acceptance source; no generated story contract or fresh PDD apply.

## Acceptance

1. Record name resolution, ping, TCP reachability, host identity, authentication and tracker round trip as separate stages. Ping MUST NOT count as authenticated SSH or working collaboration.
2. Verify host fingerprints through the machine's local console/agent before accepting keys. No passwords/private keys in prompts or reports; no disabling host-key verification or reconfiguring public routing/firewalls.
3. A Mac worker can read the canonical tracker through SSH without installing Linux PDD, GNU coreutils or flock on the Mac. Canonical writes remain serialized there. A remote worker MUST NOT create tracker items via `ym new`.
4. Different machines read the same task thread and exchange a bounded hello/answer, then use the established CLAIM/PROGRESS/SUBMIT/REVIEW protocol in separate branches/clones.
5. A sleeping/disconnected machine or missing credentials produces an explicit blocked stage; it MUST NOT be reported connected or done.
6. LAN addresses, hostnames and usernames stay in untracked local configuration/private workspace notes. Public prompts use placeholders. Never copy the main workspace memory, provider credentials or entire home directories to peers.
7. Every session records files changed, checks, limitations and next steps in the shared handoff and tracker for Claude Code/Fable. A file being written does not imply a live session read it; renewed model usage is not consumed automatically.

## Planned evidence

YM-O02 and launcher regression tests, followed by a real authenticated SSH probe, same-item read and hello/reply. Shell fixtures are regression evidence, not a Mac/network end-to-end pass.
