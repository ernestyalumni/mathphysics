import datetime as dt
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import subprocess
import sys
import time

spec = importlib.util.spec_from_file_location('ym_watch', Path(__file__).with_name('ym_watch.py'))
watch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(watch)
SOURCES = [{'id': 'clay', 'url': 'https://www.claymath.org/millennium/yang-mills-the-maths-gap/'}]
NOW = dt.datetime(2026, 11, 16, 18, tzinfo=dt.timezone.utc)


class WatchContracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.repo = self.root / "repo"
        self.repo.mkdir()

    def run_watch(self, fetcher, **kwargs):
        return watch.run_watch(self.repo, self.root/'state', SOURCES, fetcher, now=kwargs.pop('now', NOW), **kwargs)

    def test_first_seen_change_and_no_change_are_not_semantic_verdicts(self):
        for payload, expected in [(b'a', 'first-seen'), (b'a', 'unchanged'), (b'b', 'changed')]:
            now = NOW + dt.timedelta(seconds=len(list((self.root/'state/runs').glob('*'))) if (self.root/'state/runs').exists() else 0)
            result = self.run_watch(lambda _: payload, now=now)
            self.assertEqual(result['sources'][0]['status'], expected)
            self.assertEqual(result['semantic_review'], 'not-run')

    def test_failure_preserves_last_good_hash_and_is_never_no_change(self):
        self.run_watch(lambda _: b'good')
        before = json.loads((self.root/'state/state.json').read_text())['source_hashes']
        def fail(_): raise TimeoutError('fixture failure')
        result = self.run_watch(fail, now=NOW+dt.timedelta(seconds=1))
        self.assertEqual(result['status'], 'partial')
        self.assertEqual(result['sources'][0]['status'], 'unavailable')
        after = json.loads((self.root/'state/state.json').read_text())['source_hashes']
        self.assertEqual(after, before)

    def test_schedule_does_no_network_before_due_or_twice_per_date(self):
        calls = []
        fetch = lambda url: calls.append(url) or b'ok'
        early = dt.datetime(2026,11,16,16,tzinfo=dt.timezone.utc)  # 08:00 PST
        self.assertEqual(self.run_watch(fetch, daily_at='09:00', now=early)['status'], 'skipped-before-due')
        self.assertEqual(calls, [])
        self.assertEqual(self.run_watch(fetch, daily_at='09:00')['status'], 'collected')
        self.assertEqual(self.run_watch(fetch, daily_at='09:00')['status'], 'skipped-already-attempted')
        self.assertEqual(len(calls), 1)

    def test_source_urls_and_redirects_cannot_target_private_or_unlisted_hosts(self):
        for url in ['http://www.claymath.org/', 'https://127.0.0.1/',
                    'https://example.com/', 'https://token@www.claymath.org/',
                    'https://www.claymath.org:8443/']:
            with self.assertRaises(ValueError): watch.validate_url(url)
        with self.assertRaises(ValueError):
            watch.PublicRedirects().redirect_request(None, None, 302, '', {}, 'http://127.0.0.1/')

    def test_rejects_empty_duplicate_or_malformed_sources_before_state_creation(self):
        for sources in [[], [SOURCES[0], SOURCES[0]], ['bad'], [{'id': 'x', 'url': 3}]]:
            with self.assertRaises(ValueError):
                watch.run_watch(self.repo, self.root/'new-state', sources, lambda _: b'x')
        self.assertFalse((self.root/'new-state').exists())

    def test_structurally_bad_state_is_not_a_successful_daily_skip(self):
        state = self.root/'state'
        state.mkdir()
        text = json.dumps({'attempt_date': '2026-11-16', 'source_hashes': None})
        (state/'state.json').write_text(text)
        with self.assertRaises(ValueError): self.run_watch(lambda _: b'no', daily_at='09:00')
        self.assertEqual((state/'state.json').read_text(), text)

    def test_output_inside_repo_is_rejected(self):
        with self.assertRaises(ValueError):
            watch.run_watch(self.repo, self.repo/'private', SOURCES, lambda _: b'x')
        self.assertFalse((self.repo/'private').exists())

    def test_slow_drip_worker_is_killed_by_wall_clock_deadline(self):
        real_run = subprocess.run
        def slow_worker(command, **kwargs):
            return real_run([sys.executable, '-c',
                "import sys,time\nwhile True: sys.stdout.write('x'); sys.stdout.flush(); time.sleep(.01)"], **kwargs)
        started = time.monotonic()
        with mock.patch.object(watch, 'FETCH_DEADLINE_SECONDS', .15), mock.patch.object(watch.subprocess, 'run', side_effect=slow_worker):
            with self.assertRaises(TimeoutError): watch.fetch_url(SOURCES[0]['url'])
        self.assertLess(time.monotonic()-started, 3)

    def test_fetch_does_not_inherit_credentialed_proxy(self):
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = b'public'
        opener = mock.Mock()
        opener.open.return_value = response
        with mock.patch.dict(watch.os.environ, {'https_proxy': 'http://user:secret@proxy.invalid:8888'}), mock.patch.object(watch.urllib.request, 'build_opener', return_value=opener) as build:
            self.assertEqual(watch.fetch_url_in_process(SOURCES[0]['url']), b'public')
        handlers = build.call_args.args
        proxy = next(h for h in handlers if isinstance(h, watch.urllib.request.ProxyHandler))
        self.assertEqual(proxy.proxies, {})

    def test_dst_and_next_local_date_get_one_attempt_each(self):
        calls=[]
        fetch=lambda url: calls.append(url) or b'public'
        day1=dt.datetime(2026,10,31,16,tzinfo=dt.timezone.utc)
        early=dt.datetime(2026,11,1,16,30,tzinfo=dt.timezone.utc)
        due=dt.datetime(2026,11,1,17,tzinfo=dt.timezone.utc)
        self.assertEqual(self.run_watch(fetch,daily_at='09:00',now=day1)['status'],'collected')
        self.assertEqual(self.run_watch(fetch,daily_at='09:00',now=early)['status'],'skipped-before-due')
        self.assertEqual(self.run_watch(fetch,daily_at='09:00',now=due)['status'],'collected')
        self.assertEqual(self.run_watch(fetch,daily_at='09:00',now=due)['status'],'skipped-already-attempted')
        self.assertEqual(len(calls),2)

    def test_corrupt_state_is_not_silently_overwritten(self):
        state = self.root/'state'
        state.mkdir()
        (state/'state.json').write_text('bad state')
        with self.assertRaises(json.JSONDecodeError): self.run_watch(lambda _: b'no')
        self.assertEqual((state/'state.json').read_text(), 'bad state')


if __name__ == '__main__':
    unittest.main(verbosity=2)
