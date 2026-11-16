import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

SCRIPT = Path(os.environ.get('YM_SCRIPT_UNDER_TEST', str(Path(__file__).resolve().parents[1] / 'bin/ym')))


class TrackerLauncher(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='ym launcher ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bin = self.root / 'tools'
        self.bin.mkdir()
        self.launcher = self.root / 'docs/plans/yang-mills-mass-gap/bin/ym'
        self.launcher.parent.mkdir(parents=True)
        self.launcher.write_bytes(SCRIPT.read_bytes())
        (self.root / '.pdd').mkdir()
        self.capture = self.root / 'capture.json'
        self.env = dict(os.environ, PATH=str(self.bin)+':'+os.environ['PATH'],
                        CAPTURE=str(self.capture), YM_PDD=str(self.bin / 'pdd'))
        for key in ('YM_HOST', 'YM_REMOTE_YM', 'YM_MATHPHYSICS_ROOT'):
            self.env.pop(key, None)
        recorder = '#!/usr/bin/env python3\nimport json,os,sys\nfrom pathlib import Path\nPath(os.environ["CAPTURE"]).write_text(json.dumps({"args":sys.argv[1:],"cwd":os.getcwd(),"local":os.getenv("PDD_LOCAL_ONLY")}))\nprint("fixture-ok")\n'
        for name in ('pdd', 'ssh'):
            p = self.bin / name
            p.write_text(recorder)
            p.chmod(0o755)
        # Model macOS/BSD readlink and the absence of GNU realpath.
        p = self.bin / 'readlink'
        p.write_text('#!/bin/sh\nif [ "$1" = -f ]; then echo "readlink: illegal option -- f" >&2; exit 64; fi\nexec /usr/bin/readlink "$@"\n')
        p.chmod(0o755)
        p = self.bin / 'realpath'
        p.write_text('#!/bin/sh\necho "realpath unavailable" >&2\nexit 127\n')
        p.chmod(0o755)

    def run_ym(self, *args, launcher=None):
        return subprocess.run(['/bin/bash', str(launcher or self.launcher), *args],
                              env=self.env, cwd=self.root, text=True,
                              capture_output=True, timeout=10)

    def captured(self):
        return json.loads(self.capture.read_text())

    def test_local_list_uses_canonical_root_and_local_only(self):
        r = self.run_ym('list')
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.captured()['cwd'], str(self.root))
        self.assertEqual(self.captured()['args'], ['--no-core-dump', 'work', 'list', '--state', 'open'])
        self.assertEqual(self.captured()['local'], '1')

    def test_remote_worker_needs_neither_local_pdd_nor_flock(self):
        self.env.update(YM_HOST='worker@canonical.invalid', YM_REMOTE_YM='/research checkout/bin/ym', YM_PDD='/missing')
        r = self.run_ym('show', '4')
        self.assertEqual(r.returncode, 0, r.stderr)
        args = self.captured()['args']
        self.assertIn('BatchMode=yes', args)
        self.assertEqual(args[-2], 'worker@canonical.invalid')
        self.assertEqual(shlex.split(args[-1]), ['YM_HOST=', '/research checkout/bin/ym', 'show', '4'])

    def test_env_file_routes_remote(self):
        (self.root / '.pdd/ym-lan.env').write_text('YM_HOST=worker@canonical.invalid\nYM_REMOTE_YM=/canonical/bin/ym\n')
        r = self.run_ym('list')
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.captured()['args'][-2], 'worker@canonical.invalid')

    def test_remote_new_is_rejected_without_ssh(self):
        self.env['YM_HOST'] = 'worker@canonical.invalid'
        r = self.run_ym('new', 'title', 'packet.md')
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertFalse(self.capture.exists())

    def test_remote_comment_preserves_metacharacters(self):
        self.env.update(YM_HOST='worker@canonical.invalid', YM_REMOTE_YM='/research/bin/ym')
        payload = "QUESTION a'b ; $(touch SHOULD_NOT_EXIST) `id`"
        r = self.run_ym('say', '4', 'codex-fixture', payload)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(shlex.split(self.captured()['args'][-1])[-1], payload)
        self.assertFalse((self.root / 'SHOULD_NOT_EXIST').exists())

    def test_local_new_resolves_relative_packet_without_realpath(self):
        packet = self.root / 'a packet.md'
        packet.write_text('task')
        r = self.run_ym('new', 'title', packet.name)
        self.assertEqual(r.returncode, 0, r.stderr)
        args = self.captured()['args']
        self.assertEqual(args[args.index('--body-file')+1], str(packet))

    def test_symlink_launcher_keeps_canonical_root(self):
        link = self.root / 'shortcut'
        link.symlink_to(self.launcher.relative_to(self.root))
        r = self.run_ym('list', launcher=link)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.captured()['cwd'], str(self.root))


if __name__ == '__main__':
    unittest.main(verbosity=2)
