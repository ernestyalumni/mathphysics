# Deployment receipt — 16 November 2026 UTC

The reviewed collector is installed in a dedicated private directory on the Mac Mini. A single tagged user-crontab entry invokes an hourly guard; the collector makes at most one attempt per Los Angeles date at/after09:00. Previous crontab content was backed up and non-job entries preserved. Cron does not wake the Mac or backfill missed dates. Machine addresses, paths and backup location are in the private connection receipt.

Native Mac watcher contracts11/11 pass; SHA256 of deployed script matches the independently reviewed8c9d95f02df5019eb3ace032ea220ac78994d210e498c407ab38a6e61d7f8551. Linux contracts11/11 pass. First live collection finished in18.09s with8/10 successful sources. OpenAI's announcement returned403; arXiv returned429. Both remain unavailable/unknown, not unchanged. A second same-date runner invocation returned skipped-already-attempted. No bypass, credentials, provider calls or same-day retry loop.

The repo-versioned daily-watch-agent.md is copied into each run as a queued review prompt. This does not run an AI daily. The current source-status baseline was checked by an explicitly requested research subagent; future source changes need a separately assigned reviewer. This default preserves the user's limited model usage. Raw source snapshots, reports, state and cron logs remain outside Git.
