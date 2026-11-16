#!/usr/bin/env python3
"""Bounded public-source collection; never a mathematical verdict or AI call."""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from zoneinfo import ZoneInfo

ALLOWED_HOSTS = {'raw.githubusercontent.com', 'api.github.com', 'www.claymath.org',
                 'claymath.org', 'openai.com', 'www.openai.com', 'export.arxiv.org'}
FETCH_DEADLINE_SECONDS = 15
LIMIT_BYTES = 4 * 1024 * 1024
CONTEXT_FILES = ['docs/plans/yang-mills-mass-gap/SOURCES.md',
                 'docs/plans/yang-mills-mass-gap/BOARD.md',
                 'docs/plans/yang-mills-mass-gap/HANDOFF.md',
                 'docs/plans/yang-mills-mass-gap/ledger/claims.md',
                 'documents/mono-notes/topics/yang-mills-source-map.tex']


def validate_url(url):
    parsed = urllib.parse.urlsplit(url)
    if (parsed.scheme != 'https' or parsed.hostname not in ALLOWED_HOSTS or
            parsed.username or parsed.password or parsed.port not in (None, 443)):
        raise ValueError('Source URL must be HTTPS on the fixed public allowlist')
    return url


class PublicRedirects(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        validate_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def fetch_url_in_process(url):
    validate_url(url)
    request = urllib.request.Request(url, headers={
        'User-Agent': 'mathphysics-public-source-watch/1.0',
        'Accept': 'application/json,text/html,text/plain,application/atom+xml,*/*'})
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}), PublicRedirects())
    with opener.open(request, timeout=15) as response:
        data = response.read(LIMIT_BYTES + 1)
    if len(data) > LIMIT_BYTES:
        raise ValueError('Source exceeds the 4 MiB collection limit')
    return data


def fetch_url(url):
    validate_url(url)
    try:
        result = subprocess.run(
            [sys.executable, '-B', str(Path(__file__).resolve()), '--fetch-url', url],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=FETCH_DEADLINE_SECONDS, check=False)
    except subprocess.TimeoutExpired as error:
        raise TimeoutError('Hard fetch wall-clock deadline exceeded') from error
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors='replace')[-1000:])
    if len(result.stdout) > LIMIT_BYTES:
        raise ValueError('Source exceeds the 4 MiB collection limit')
    return result.stdout


def validate_sources(sources):
    if not isinstance(sources, list) or not sources or len(sources) > 32:
        raise ValueError('Sources must be a nonempty list of at most 32 entries')
    seen = set()
    for source in sources:
        if not isinstance(source, dict):
            raise ValueError('Each source must be an object')
        identifier, url = source.get('id'), source.get('url')
        if (not isinstance(identifier, str) or not identifier or not identifier.isascii()
                or not identifier.replace('-', '').replace('_', '').isalnum()
                or identifier in seen):
            raise ValueError('Source identifiers must be valid and unique')
        if not isinstance(url, str):
            raise ValueError('Source URL must be a string')
        validate_url(url)
        seen.add(identifier)


def validate_state(state):
    if not isinstance(state, dict):
        raise ValueError('State must be an object')
    if not state:
        return
    if not isinstance(state.get('attempt_date'), str):
        raise ValueError('State date must be a string')
    dt.date.fromisoformat(state['attempt_date'])
    hashes = state.get('source_hashes')
    if not isinstance(hashes, dict):
        raise ValueError('State source hashes must be an object')
    for key, value in hashes.items():
        if (not isinstance(key, str) or not isinstance(value, str) or len(value) != 64
                or any(c not in '0123456789abcdef' for c in value)):
            raise ValueError('State contains an invalid source hash')


def write_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def context(repo):
    result = {'files': {}}
    for name in CONTEXT_FILES:
        path = repo / name
        result['files'][name] = (hashlib.sha256(path.read_bytes()).hexdigest()
                                  if path.is_file() else None)
    for field, args in [('revision', ['rev-parse', 'HEAD']),
                        ('branch', ['branch', '--show-current'])]:
        p = subprocess.run(['git', '-C', str(repo), *args], capture_output=True,
                           text=True, timeout=5)
        result[field] = p.stdout.strip() if p.returncode == 0 else 'unavailable'
    return result


def run_watch(repo, state_dir, sources, fetcher=fetch_url, now=None,
              daily_at=None, timezone='America/Los_Angeles', max_seconds=180):
    repo = repo.resolve()
    state_dir = state_dir.resolve()
    try:
        state_dir.relative_to(repo)
    except ValueError:
        pass
    else:
        raise ValueError('Private collection output must be outside the repository')
    validate_sources(sources)
    now = now or dt.datetime.now(dt.timezone.utc)
    if now.tzinfo is None:
        raise ValueError('Timezone-aware time required')
    local = now.astimezone(ZoneInfo(timezone))
    today = local.date().isoformat()
    state_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (state_dir / 'run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return {'status': 'skipped-running'}
        state_path = state_dir / 'state.json'
        state = json.loads(state_path.read_text()) if state_path.exists() else {}
        validate_state(state)
        if daily_at:
            hour, minute = map(int, daily_at.split(':'))
            if not (0 <= hour < 24 and 0 <= minute < 60):
                raise ValueError('Invalid daily time')
            if (local.hour, local.minute) < (hour, minute):
                return {'status': 'skipped-before-due'}
            if state.get('attempt_date') == today:
                return {'status': 'skipped-already-attempted'}
        run_id = now.astimezone(dt.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
        run_dir = state_dir / 'runs' / run_id
        run_dir.mkdir(parents=True, exist_ok=False, mode=0o700)
        hashes = dict(state.get('source_hashes', {}))
        rows = []
        started = time.monotonic()
        for source in sources:
            row = {'id': source['id'], 'url': source['url'],
                   'previous_sha256': hashes.get(source['id'])}
            try:
                if time.monotonic() - started >= max_seconds:
                    raise TimeoutError('Run time budget exhausted; source not checked')
                data = fetcher(source['url'])
                if not isinstance(data, bytes) or len(data) > LIMIT_BYTES:
                    raise ValueError('Invalid or oversized fetch result')
                digest = hashlib.sha256(data).hexdigest()
                row.update(sha256=digest, bytes=len(data), status=(
                    'first-seen' if row['previous_sha256'] is None else
                    'unchanged' if digest == row['previous_sha256'] else 'changed'))
                (run_dir / (source['id'] + '.source')).write_bytes(data)
                hashes[source['id']] = digest
            except Exception as error:  # A missing source must remain visible.
                row.update(status='unavailable', error=f'{type(error).__name__}: {error}')
            rows.append(row)
        result = {'schema': 'ym.source-watch.v1', 'run_utc': now.isoformat(),
                  'local_date': today, 'timezone': timezone,
                  'status': 'partial' if any(r['status'] == 'unavailable' for r in rows) else 'collected',
                  'semantic_review': 'not-run', 'sources': rows,
                  'repo_context': context(repo),
                  'elapsed_seconds': round(time.monotonic() - started, 2)}
        write_json(run_dir / 'report.json', result)
        lines = [f'# Yang–Mills source watch — {today}', '',
                 f'Collection: **{result["status"]}**. Semantic AI review: **not run**.', '',
                 'Changes are review candidates, not mathematical findings. Fetch failures are unknown, never no-change.',
                 'The last reviewed baseline (2026-11-16 UTC) supports the October existence construction; it does not establish a mass gap or a unique infinite-volume limit.',
                 'This collector cannot decide whether a new paper proves a theorem, or assign a percentage toward a prize.', '',
                 '| Source | Observation |', '|---|---|']
        for row in rows:
            lines.append(f'| [{row["id"]}]({row["url"]}) | {row["status"]} |')
        lines += ['', 'Read the versioned daily-watch-agent.md prompt before reviewing source snapshots.',
                  'Keep this report local. Do not automatically update source pins, manuscript text, task completion or mathematical claim labels.', '']
        (run_dir / 'report.md').write_text('\n'.join(lines))
        prompt = Path(__file__).with_name('daily-watch-agent.md')
        if prompt.is_file():
            (run_dir / 'review-prompt.md').write_text(prompt.read_text() +
                '\n\nLocal source-collection directory: ' + str(run_dir) + '\n')
        write_json(state_path, {'attempt_date': today, 'source_hashes': hashes,
                                'latest_run': run_id, 'latest_status': result['status']})
        write_json(state_dir / 'latest.json', {'run': run_id, 'status': result['status'],
                                             'report': str(run_dir / 'report.md')})
        return result


def main():
    if len(sys.argv) == 3 and sys.argv[1] == '--fetch-url':
        sys.stdout.buffer.write(fetch_url_in_process(sys.argv[2]))
        return 0
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', required=True, type=Path)
    parser.add_argument('--state-dir', required=True, type=Path)
    parser.add_argument('--sources', type=Path, default=Path(__file__).with_name('sources.json'))
    parser.add_argument('--daily-at', help='Once per local date at or after HH:MM; no network before due')
    parser.add_argument('--timezone', default='America/Los_Angeles')
    args = parser.parse_args()
    os.umask(0o077)
    result = run_watch(args.repo, args.state_dir,
                       json.loads(args.sources.read_text()), daily_at=args.daily_at,
                       timezone=args.timezone)
    print(json.dumps({key: result[key] for key in ('status', 'semantic_review') if key in result}))
    return 1 if result['status'] == 'partial' else 0


if __name__ == '__main__':
    raise SystemExit(main())
