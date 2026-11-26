#!/usr/bin/env python3
"""Read-only interim checks; finalize a separate Qwen report after full resolution."""
import argparse
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import time

HERE = Path('<private-run>/full-replay')
BASE = HERE.parent
QWEN_UNIT = 'ym-referee-qwen38-20261126.service'
ROUTE = 'qwen38'
PRIMARY = ('kimi', 'deepseek', 'glm53')
ALL = (ROUTE,) + PRIMARY


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def library():
    lib = runpy.run_path(str(HERE / 'summarize_replay.py'))
    # Runtime-only alias; neither the frozen primary metadata nor its reporter is edited.
    lib['MODELS'].setdefault(ROUTE, 'qwen/qwen3.8-max')
    return lib


def resolved(summary):
    return summary['terminal_samples'] + summary['transport_failure_samples'] == 86


def main(wait):
    lib = library()
    primary_receipt_check = subprocess.run(
        ['sha256sum', '--quiet', '-c', 'PRIMARY-RECEIPTS.sha256'],
        cwd=HERE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert primary_receipt_check.returncode == 0, primary_receipt_check.stderr
    assert sha(BASE / 'frozen-test.json') == 'b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe'
    frozen_primary = json.loads((HERE / 'primary-summary.json').read_text())
    assert frozen_primary['complete']
    expected = {'kimi': (86, 78), 'deepseek': (74, 60), 'glm53': (86, 27)}
    for name, pair in expected.items():
        assert (frozen_primary['models'][name]['terminal_samples'],
                frozen_primary['models'][name]['primary_success_terminal']) == pair
    last = None
    while True:
        summary, samples = lib['route_summary'](ROUTE)
        assert summary['scheduled_samples'] == 86
        if resolved(summary):
            break
        last = summary
        active = subprocess.run(['systemctl', '--user', 'is-active', '--quiet', QWEN_UNIT])
        if active.returncode != 0 and not wait:
            return {'state': 'in_progress_or_not_running', 'summary': summary,
                    'final_comparison_written': False}
        if active.returncode != 0:
            raise SystemExit('Qwen runner is inactive before all 86 samples resolve; preserving partial receipts.')
        time.sleep(20)
    rows = {name: lib['route_summary'](name) for name in ALL}
    assert frozen_primary['complete']
    for name in PRIMARY:
        assert rows[name][0]['terminal_samples'] == frozen_primary['models'][name]['terminal_samples']
        assert rows[name][0]['primary_success_terminal'] == frozen_primary['models'][name]['primary_success_terminal']
    checks = ['sha256sum', '--quiet', '-c', 'PRIMARY-RECEIPTS.sha256']
    subprocess.run(checks, cwd=HERE, check=True, stdout=subprocess.PIPE)
    def pad(name):
        stat, cases = rows[name]
        returned = [x for x in cases if x['terminal']]
        parseable = [x for x in returned if x['verdict'] in ('CORRECT', 'INCORRECT')]
        return {'summary': stat, 'samples': cases, 'parseable': len(parseable),
                'success_given_parseable': sum(bool(x['success']) for x in parseable),
                'resolved': resolved(stat)}
    summaries = {name: pad(name) for name in ALL}
    indexed = {name: {(x['case'], x['sample']): x for x in summaries[name]['samples']} for name in ALL}
    matched = {}
    for title, condition in [('all_terminal', lambda x: x['terminal']),
                             ('all_parseable', lambda x: x['verdict'] in ('CORRECT', 'INCORRECT'))]:
        keys = [key for key in indexed[ROUTE] if all(condition(indexed[name][key]) for name in ALL)]
        matched[title] = {'observations': len(keys),
                          'success_by_route': {name: sum(bool(indexed[name][key]['success']) for key in keys) for name in ALL},
                          'warning': 'Post-hoc available-response subset; not an independent-case superiority test.'}
    receipt_root = HERE / 'receipts' / ROUTE
    manifest = ''
    receipt_count = 0
    allowed_cases = {c['id']: c for c in lib['CASES']}
    for path in sorted(receipt_root.glob('sample-*/*.attempt-*.json')):
        row = json.loads(path.read_text())
        assert row['schema'] == 1 and row['route_name'] == ROUTE
        assert row['route'] == 'qwen/qwen3.8-max'
        assert row['case'] in allowed_cases and row['sample'] in (1, 2)
        assert row['expected'] == allowed_cases[row['case']]['expected']
        assert row['group'] == allowed_cases[row['case']]['group']
        assert 1 <= row['attempt'] <= 3 and row['timeout_seconds'] == 240
        assert row['command_options']['thinking'] == 'low'
        assert row['command_options']['tools'] is False and row['command_options']['session'] is False
        manifest += sha(path) + '  ' + str(path.relative_to(HERE)) + '\n'
        receipt_count += 1
    assert receipt_count == summaries[ROUTE]['summary']['request_attempts']
    manifest_path = HERE / 'QWEN38-RECEIPTS.sha256'
    if manifest_path.exists():
        assert manifest_path.read_text() == manifest, 'Existing Qwen receipt manifest changed'
    else:
        tmp = manifest_path.with_name(manifest_path.name + '.tmp-' + str(os.getpid()))
        tmp.write_text(manifest)
        os.replace(tmp, manifest_path)
    csv_path = HERE / 'qwen38-case-results.csv'
    fields = ['case', 'group', 'expected', 'planted_step', 'sample']
    for name in ALL:
        fields += [name + '_' + key for key in ('category', 'verdict', 'returned_step', 'terminal', 'stop')]
    q_fields = ['requests_seconds', 'input_tokens', 'output_tokens', 'reasoning_tokens', 'cost_estimate_usd']
    fields += ['qwen38_' + key for key in q_fields]
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for case in lib['CASES']:
            for number in (1, 2):
                out = {'case': case['id'], 'group': case['group'], 'expected': case['expected'],
                       'planted_step': case.get('step'), 'sample': number}
                for name in ALL:
                    row = indexed[name][(case['id'], number)]
                    for key in ('category', 'verdict', 'returned_step', 'terminal', 'stop'):
                        out[name + '_' + key] = row[key]
                for key in q_fields:
                    out['qwen38_' + key] = indexed[ROUTE][(case['id'], number)][key]
                writer.writerow(out)
    generated = datetime.datetime.now(datetime.timezone.utc).isoformat()
    payload = {'generated_utc': generated, 'complete': True, 'route': 'qwen/qwen3.8-max',
               'scheduled_samples': 86, 'receipt_manifest': str(manifest_path),
               'receipt_manifest_sha256': sha(manifest_path),
               'primary_receipt_manifest_sha256': sha(HERE / 'PRIMARY-RECEIPTS.sha256'),
               'primary_summary_frozen_artifact_sha256': sha(HERE / 'primary-summary.json'),
               'preregistration_sha256': sha(HERE / 'PREREGISTRATION-QWEN38.md'),
               'runner_sha256': sha(HERE / 'run_full_replay.py'),
               'comparison_script_sha256': sha(Path(__file__)),
               'models': summaries, 'matched_available_response_sensitivity': matched,
               'limitations': ['Different launch/provider load; repeated samples are not independent cases.',
                               'Terminal, all-scheduled and post-hoc parseable denominators are distinct.',
                               'Mechanical labels and model agreement are not independent human certification.']}
    (HERE / 'qwen38-comparison.json').write_text(json.dumps(payload, indent=2) + '\n')
    def cell(name, key):
        return summaries[name]['summary'][key]
    def display(route):
        return {'qwen38': 'Qwen3.8-max', 'kimi': 'Kimi3.1-pro',
                'deepseek': 'DeepSeek4.1-free', 'glm53': 'GLM5.3'}[route]
    lines = ['# Qwen 3.8 Max versus immutable primary models', '',
             f'Generated UTC: {generated}. All Qwen observations are fully resolved against the same 43-case set, two samples per case.',
             'No expected label was supplied to any model. This is a separate comparison; the completed primary report was not overwritten.', '',
             '| Metric | ' + ' | '.join(display(n) for n in ALL) + ' |',
             '|---|' + '---:|' * len(ALL)]
    metrics = [('Terminal responses /86', 'terminal_samples'),
               ('Transport failures /86', 'transport_failure_samples'),
               ('Primary success / terminal', 'primary_success_terminal'),
               ('Primary success / all 86', 'primary_success_all_scheduled_sensitivity'),
               ('Parseable verdicts', None),
               ('Success given parseable', None),
               ('Originals passed /28', 'originals_passed'),
               ('Original false alarms', 'original_false_alarms'),
               ('Exact planted steps /58', 'mutants_exact_step'),
               ('Mutants passed', 'mutants_passed_as_correct'),
               ('Wrong/missing step', 'mutants_wrong_or_missing_step'),
               ('Unparsed', 'unparsed_terminal_samples'),
               ('Token truncations', 'token_truncated_terminal_samples'),
               ('Requests', 'request_attempts'), ('Retries', 'retry_attempts'),
               ('Timeouts', 'timeout_requests'), ('Input tokens', 'input_tokens_all_requests'),
               ('Output tokens', 'output_tokens_all_requests'), ('Reasoning tokens', 'reasoning_tokens_all_requests'),
               ('Mean seconds per terminal sample', 'terminal_sample_seconds_mean'),
               ('Estimated cost USD', 'cost_estimate_usd_all_requests')]
    for label, key in metrics:
        vals = []
        for name in ALL:
            if key == 'primary_success_terminal':
                val = f"{cell(name, key)}/{cell(name, 'terminal_samples')}"
            elif label == 'Parseable verdicts':
                val = summaries[name]['parseable']
            elif label == 'Success given parseable':
                val = f"{summaries[name]['success_given_parseable']}/{summaries[name]['parseable']}"
            else:
                val = cell(name, key)
            vals.append(f'${val:.6f}' if isinstance(val, float) and 'cost' in label.lower() else str(val))
        lines.append('| ' + label + ' | ' + ' | '.join(vals) + ' |')
    lines += ['', '## Secondary matched-observation views', '',
              'These subsets condition on all four routes producing terminal or parseable responses; they are post-hoc, selected views, not independent-case superiority estimates.']
    for name, m in matched.items():
        lines += ['', f'- **{name}:** {m["observations"]} observations. Successes: ' +
                  ', '.join(f'{display(n)} {score}' for n, score in m['success_by_route'].items()) + '.']
    lines += ['', '## Evidence boundary', '',
              'Exact per-case results are in `qwen38-case-results.csv` and `qwen38-comparison.json`. All raw Qwen requests and failures remain separately hashed in `QWEN38-RECEIPTS.sha256`. Usage and costs are Pi/TokenRouter estimates, not invoices; requests lacking usage cannot contribute unavailable token totals. Truncated or unparsed outputs, transport failures and mathematical failures remain distinct.',
              'These mechanical historical labels are not independent human certification. No agreement among models, successful compilation or finite-volume numerical result establishes a continuum Yang–Mills mass-gap proof.', '']
    (HERE / 'QWEN38-COMPARISON.md').write_text('\n'.join(lines))
    return {'state': 'complete', 'qwen_terminal': cell(ROUTE, 'terminal_samples'),
            'qwen_success_terminal': cell(ROUTE, 'primary_success_terminal'),
            'qwen_success_all_scheduled': cell(ROUTE, 'primary_success_all_scheduled_sensitivity'),
            'qwen_receipt_manifest_sha256': sha(manifest_path),
            'primary_receipts_verified_again': True}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--wait', action='store_true')
    ap.add_argument('--self-test', action='store_true')
    args = ap.parse_args()
    if args.self_test:
        lib = library()
        runner = runpy.run_path(str(HERE / 'run_full_replay.py'), run_name='offline_qwen_route_check')
        assert runner['MODELS'][ROUTE] == ('qwen/qwen3.8-max', 2)
        assert runner['TIMEOUT'] == 240 and runner['MAX_ATTEMPTS'] == 3
        s, samples = lib['route_summary'](ROUTE)
        assert s['scheduled_samples'] == 86
        print(json.dumps({'offline_checks': True, 'model_calls': 0,
                          'scheduled_samples': s['scheduled_samples'],
                          'resolved_at_check': s['terminal_samples'] + s['transport_failure_samples'],
                          'final_report_written': False}, indent=2))
    else:
        result = main(args.wait)
        if result:
            print(json.dumps(result, indent=2), flush=True)
        if not args.wait and result['state'] != 'complete':
            raise SystemExit(1)
