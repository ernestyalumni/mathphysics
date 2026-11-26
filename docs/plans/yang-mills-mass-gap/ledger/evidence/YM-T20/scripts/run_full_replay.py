#!/usr/bin/env python3
"""Resumable, receipt-preserving hosted replay for the frozen YM referee set."""
import concurrent.futures
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import runpy
import signal
import subprocess
import sys
import threading
import time

ROOT = Path(os.environ['YM_REFEREE_ROOT'])
BASE = Path(os.environ['YM_REFEREE_BASE'])
HERE = BASE / 'full-replay'
RECEIPTS = HERE / 'receipts'
CASE_FILE = BASE / 'frozen-test.json'
CASE_SHA = 'b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe'
PARSER_FILE = ROOT / 'repos/mathphysics/docs/plans/yang-mills-mass-gap/bin/ym-micro'
MODELS = {
    'deepseek': ('deepseek/deepseek-v4.1-flash-1017-free', 8),
    'kimi': ('moonshotai/kimi-k3.1-pro', 2),
    'flashx': ('z-ai/glm-5.3-flashx', 2),
    'glm53': ('z-ai/glm-5.3', 2),
    'minimax': ('MiniMax-M3', 2),
    'qwen38': ('qwen/qwen3.8-max', 2),
}
TIMEOUT = 240
MAX_ATTEMPTS = 3
SYSTEM_PROMPT = 'Review the mathematical argument and follow the requested answer format.'
PRINT_LOCK = threading.Lock()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.tmp-{os.getpid()}-{threading.get_ident()}')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    os.replace(temp, path)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def redact(value, limit=2000):
    text = str(value or '')
    text = re.sub(r'Bearer\s+[^\s,;]+', 'Bearer [redacted]', text, flags=re.I)
    text = re.sub(r'\b(sk|key|token)-[A-Za-z0-9_.-]{8,}', r'\1-[redacted]', text)
    return text[-limit:]


def load_receipts(route, sample, case_id):
    folder = RECEIPTS / route / f'sample-{sample}'
    rows = []
    for path in sorted(folder.glob(case_id + '.attempt-*.json')):
        rows.append(json.loads(path.read_text()))
    rows.sort(key=lambda row: row['attempt'])
    return rows


def is_terminal(receipt):
    if receipt.get('transport_failure'):
        return False
    return bool(receipt.get('has_assistant_message')) and receipt.get('stop') != 'error'


def invoke(route, sample, case, attempt):
    model = MODELS[route][0]
    cmd = [
        'pi', '-p', '--mode', 'json', '--no-session', '--no-tools', '--no-mcp',
        '--no-extensions', '--no-skills', '--no-context-files', '--system-prompt',
        SYSTEM_PROMPT, '--provider', 'tokenrouter', '--model', model, '--thinking',
        'low', case['prompt'],
    ]
    started = time.time()
    proc = subprocess.Popen(
        cmd, cwd='/tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, start_new_session=True,
    )
    timed_out = False
    try:
        stdout, stderr = proc.communicate(timeout=TIMEOUT)
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        stdout, stderr = proc.communicate()
        returncode = None
    msg = {}
    event_types = []
    for line in stdout.splitlines():
        try:
            event = json.loads(line)
        except ValueError:
            continue
        if isinstance(event, dict):
            event_types.append(event.get('type'))
            if event.get('type') == 'message_end' and event.get('message', {}).get('role') == 'assistant':
                msg = event['message']
    answer = '\n'.join(
        block.get('text', '') for block in msg.get('content', [])
        if isinstance(block, dict) and block.get('type') == 'text'
    )
    has_message = bool(msg)
    stop = msg.get('stopReason')
    transport_failure = timed_out or returncode != 0 or not has_message or stop == 'error'
    parsed = PARSER(answer)
    expected = case['expected']
    planted = case.get('step')
    located_or_passed = (
        parsed['verdict'] == expected
        and (planted is None or parsed.get('step') == str(planted))
    )
    out_dir = RECEIPTS / route / f'sample-{sample}'
    dest = out_dir / f"{case['id']}.attempt-{attempt}.json"
    receipt = {
        'schema': 1,
        'route_name': route,
        'route': model,
        'reported_model': msg.get('model'),
        'case': case['id'],
        'group': case['group'],
        'sample': sample,
        'attempt': attempt,
        'prompt_sha256': hashlib.sha256(case['prompt'].encode()).hexdigest(),
        'expected': expected,
        'step': planted,
        'parsed': parsed,
        'answer': answer,
        'has_assistant_message': has_message,
        'stop': stop,
        'exit': 'timeout' if timed_out else returncode,
        'transport_failure': transport_failure,
        'terminal_response': not transport_failure,
        'located_or_passed': located_or_passed if not transport_failure else None,
        'usage': msg.get('usage'),
        'error': redact(msg.get('errorMessage'), 600),
        'stderr_tail': redact(stderr, 1200),
        'stdout_tail_if_no_message': '' if has_message else redact(stdout, 1200),
        'event_types': event_types[-20:],
        'seconds': round(time.time() - started, 1),
        'timeout_seconds': TIMEOUT,
        'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'command_options': {
            'thinking': 'low', 'tools': False, 'session': False,
            'system_prompt': SYSTEM_PROMPT, 'provider': 'tokenrouter',
        },
    }
    atomic_json(dest, receipt)
    with PRINT_LOCK:
        print(
            f"{route} s{sample} a{attempt} {case['id']} "
            f"{parsed['verdict']} step={parsed.get('step')} "
            f"located_or_passed={receipt['located_or_passed']} "
            f"exit={receipt['exit']} stop={stop} seconds={receipt['seconds']}",
            flush=True,
        )
    return receipt


def stop_requested():
    return (HERE / 'STOP').exists()


def run_sample(route, sample, case):
    if stop_requested():
        return 'stopped'
    while True:
        rows = load_receipts(route, sample, case['id'])
        if any(is_terminal(row) for row in rows) or len(rows) >= MAX_ATTEMPTS:
            return 'done'
        if stop_requested():
            return 'stopped'
        attempt = len(rows) + 1
        receipt = invoke(route, sample, case, attempt)
        if is_terminal(receipt) or attempt >= MAX_ATTEMPTS:
            return 'done'
        delay = attempt * 3 + int(hashlib.sha256(case['id'].encode()).hexdigest()[:4], 16) % 7
        for _ in range(delay):
            if stop_requested():
                return 'stopped'
            time.sleep(1)


def run_route(route, cases):
    model, workers = MODELS[route]
    atomic_json(HERE / f'status-{route}.json', {
        'state': 'running', 'model': model, 'concurrency': workers,
        'scheduled_samples': 2 * len(cases), 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    })
    jobs = [(sample, case) for case in cases for sample in (1, 2)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        outcomes = list(pool.map(lambda job: run_sample(route, job[0], job[1]), jobs))
    state = 'stopped' if 'stopped' in outcomes else 'complete'
    all_receipts = [row for sample, case in jobs
                    for row in load_receipts(route, sample, case['id'])]
    atomic_json(HERE / f'status-{route}.json', {
        'state': state, 'model': model, 'concurrency': workers,
        'scheduled_samples': len(jobs), 'request_receipts': len(all_receipts),
        'terminal_request_receipts': sum(is_terminal(row) for row in all_receipts),
        'transport_request_receipts': sum(not is_terminal(row) for row in all_receipts),
        'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    })
    print(f'ROUTE {route} {state}: {len(all_receipts)} request receipts', flush=True)


def main():
    os.umask(0o077)
    if not os.path.ismount(os.environ.get('YM_REFEREE_MOUNT', '/')):
        raise SystemExit('Expansion is not mounted; refusing to write elsewhere')
    if sha256(CASE_FILE) != CASE_SHA:
        raise SystemExit('Frozen case SHA-256 mismatch')
    cases = json.loads(CASE_FILE.read_text())
    if len(cases) != 43 or sum(c['expected'] == 'CORRECT' for c in cases) != 14:
        raise SystemExit('Unexpected frozen benchmark composition')
    requested = sys.argv[1:] or list(MODELS)
    unknown = set(requested) - set(MODELS)
    if unknown:
        raise SystemExit(f'Unknown routes: {sorted(unknown)}')
    HERE.mkdir(parents=True, exist_ok=True)
    special = requested[0] if len(requested) == 1 and requested[0] in {'glm53', 'minimax', 'qwen38'} else None
    metadata_name = f'run-metadata-{special}.json' if special else 'run-metadata.json'
    prereg_name = f'PREREGISTRATION-{special.upper()}.md' if special else 'PREREGISTRATION.md'
    metadata_file = HERE / metadata_name
    preregistration_file = HERE / prereg_name
    if special == 'qwen38':
        # Local version checks only; execute no mathematical request in this preflight.
        assert subprocess.check_output(['node', '--version'], text=True).strip() == 'v24.18.1'
        assert subprocess.check_output([str(Path.home() / '.local/bin/pi'), '--version'], text=True).strip() == '1.6.0'
    if special == 'qwen38' and metadata_file.exists():
        prior = json.loads(metadata_file.read_text())
        for key, expected in [('case_sha256', CASE_SHA), ('samples_per_case', 2),
                              ('timeout_seconds', TIMEOUT),
                              ('max_attempts_per_sample', MAX_ATTEMPTS),
                              ('parser_sha256', sha256(PARSER_FILE)),
                              ('runner_sha256', sha256(Path(__file__))),
                              ('preregistration_sha256', sha256(preregistration_file))]:
            assert prior[key] == expected, f'Qwen frozen field changed: {key}'
        assert prior['requested_routes'] == ['qwen38']
        assert prior['models']['qwen38'] == MODELS['qwen38'][0]
        assert prior['concurrency']['qwen38'] == 2
    atomic_json(metadata_file, {
        'case_file': str(CASE_FILE), 'case_sha256': CASE_SHA,
        'case_count': len(cases), 'proof_groups': len({c['group'] for c in cases}),
        'samples_per_case': 2, 'models': {k: v[0] for k, v in MODELS.items()},
        'concurrency': {k: v[1] for k, v in MODELS.items()},
        'timeout_seconds': TIMEOUT, 'max_attempts_per_sample': MAX_ATTEMPTS,
        'parser_sha256': sha256(PARSER_FILE),
        'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'preregistration_file': str(preregistration_file),
        'preregistration_sha256': sha256(preregistration_file),
        'started_or_resumed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'pid': os.getpid(), 'requested_routes': requested,
        **({'runtime_node': 'v24.18.1', 'runtime_pi': '1.6.0',
            'thinking': 'low', 'output_cap': 'provider_default',
            'provider_catalog_context': '128K', 'provider_catalog_displayed_output': '16.4K'} if special == 'qwen38' else {}),
    })
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(requested)) as pool:
        list(pool.map(lambda route: run_route(route, cases), requested))


if __name__ == '__main__':
    lib = runpy.run_path(str(PARSER_FILE))
    global PARSER
    PARSER = lib['parse_verdict']
    main()
