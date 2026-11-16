# Daily source watch and queued AI review

The collector is a Python standard-library operations utility. Symbolic research remains Rust first. It reads ten fixed public primary-source feeds/pages and writes private snapshots, hashes and a local report. It never runs a model, spends provider usage, updates source pins or promotes a mathematical claim. Dynamic HTML changes are review candidates, not reliable semantic diffs.

`daily-watch-agent.md` is the versioned context for an AI reviewer. Every collection queues that prompt with its snapshot directory. This session's research subagent performed a current primary-source audit; future reports are **not AI-reviewed** until an explicitly selected agent actually reviews them. No silent paid fallback.

Run manually:

```sh
python3 -B docs/plans/yang-mills-mass-gap/watch/ym_watch.py \
  --repo "$PWD" --state-dir /absolute/private/ym-watch-state
```

Use `--daily-at 09:00 --timezone America/Los_Angeles` for a scheduled invocation. An hourly cron guard is used so time-zone/DST handling is explicit and a sleeping Mac gets a later opportunity after waking. It makes at most one collection attempt per Los Angeles date, at/after 09:00. Plain cron does not wake or power on a Mac; missed dates are not backfilled. Each source has a hard 15s wall-clock fetch deadline, 4 MiB body cap and the run has a 180s collection-loop budget. Errors remain visible; there are no automatic repeated network retries on the same date. The next date attempts again.

Installer instructions and machine-local paths belong in the private connection receipt, not this public file. Preserve existing crontab entries, back up the old crontab, and add only the identifiable mathphysics-ym-watch entry. Uninstall by removing only that entry. Do not enable public ports or change sleep/power settings. Reports and downloaded source bodies stay outside Git, default state directory permissions 700 and file umask 077.

Validation: `python3 -B docs/plans/yang-mills-mass-gap/watch/test_ym_watch.py`. Tests check change detection, fetch-failure honesty, preserved hashes, scheduling without redundant requests, public URL/redirect restrictions and corrupt-state preservation. A live collection and a second same-day no-network invocation are additional deployment checks.

The fetch runs in a bounded child process so a slow-drip server cannot evade the deadline. Environment proxy settings are disabled rather than inheriting proxy credentials. Nonempty source/schema/unique-ID and state/hash validation happen before daily skip decisions. Output paths beneath the selected repository are rejected.
