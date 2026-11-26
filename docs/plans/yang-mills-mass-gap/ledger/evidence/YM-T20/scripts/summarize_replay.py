#!/usr/bin/env python3
"""Deterministically summarize immutable full-replay receipts without judging mathematics."""
import csv
import datetime
import json
import os
from pathlib import Path
import statistics

ROOT = Path(os.environ['YM_REFEREE_ROOT'])
HERE = Path(os.environ['YM_REFEREE_BASE']) / 'full-replay'
CASES = json.loads((HERE.parent / 'frozen-test.json').read_text())
META = json.loads((HERE / 'run-metadata.json').read_text())
MODELS = dict(META['models'])
MODELS['glm53'] = 'z-ai/glm-5.3'
MODELS['minimax'] = 'MiniMax-M3'


def usage_sum(receipts, key):
    return sum((row.get('usage') or {}).get(key) or 0 for row in receipts)


def cost_sum(receipts):
    return sum(((row.get('usage') or {}).get('cost') or {}).get('total') or 0 for row in receipts)


def route_summary(route):
    requests = []
    for path in (HERE / 'receipts' / route).glob('sample-*/*.attempt-*.json'):
        requests.append(json.loads(path.read_text()))
    requests.sort(key=lambda row: (row['case'], row['sample'], row['attempt']))
    samples = []
    for case in CASES:
        for number in (1, 2):
            attempts = [r for r in requests if r['case'] == case['id'] and r['sample'] == number]
            final = attempts[-1] if attempts else None
            terminal = bool(final and final['terminal_response'])
            parsed = (final or {}).get('parsed') or {}
            expected = case['expected']
            planted = case.get('step')
            verdict = parsed.get('verdict') if terminal else None
            step = parsed.get('step') if terminal else None
            if not terminal:
                if final is None:
                    category = 'pending'
                elif len(attempts) < META['max_attempts_per_sample'] and not (HERE / 'STOP').exists():
                    category = 'retry_pending'
                else:
                    category = 'transport_failure'
                success = False if category == 'transport_failure' else None
            elif expected == 'CORRECT' and verdict == 'CORRECT':
                category, success = 'original_passed', True
            elif expected == 'CORRECT' and verdict == 'INCORRECT':
                category, success = 'original_false_alarm', False
            elif expected == 'INCORRECT' and verdict == 'INCORRECT' and step == str(planted):
                category, success = 'mutant_exact_step', True
            elif expected == 'INCORRECT' and verdict == 'CORRECT':
                category, success = 'mutant_passed', False
            elif expected == 'INCORRECT' and verdict == 'INCORRECT':
                category, success = 'mutant_wrong_or_missing_step', False
            else:
                category, success = 'unparsed', False
            samples.append({
                'route': route, 'case': case['id'], 'group': case['group'],
                'sample': number, 'expected': expected, 'planted_step': planted,
                'attempts': len(attempts), 'terminal': terminal,
                'category': category, 'success': success, 'verdict': verdict,
                'returned_step': step, 'stop': (final or {}).get('stop'),
                'requests_seconds': round(sum(r['seconds'] for r in attempts), 1),
                'input_tokens': usage_sum(attempts, 'input'),
                'output_tokens': usage_sum(attempts, 'output'),
                'reasoning_tokens': usage_sum(attempts, 'reasoning'),
                'cost_estimate_usd': cost_sum(attempts),
            })
    complete = [row for row in samples if row['attempts']]
    observed = [row for row in complete if row['category'] == 'transport_failure']
    terminal = [row for row in samples if row['terminal']]
    count = lambda category: sum(row['category'] == category for row in samples)
    pairs = []
    for case in CASES:
        pair = [row for row in samples if row['case'] == case['id']]
        if all(row['terminal'] for row in pair):
            pairs.append({
                'case': case['id'],
                'same_verdict': pair[0]['verdict'] == pair[1]['verdict'],
                'same_verdict_and_step': (pair[0]['verdict'], pair[0]['returned_step']) == (pair[1]['verdict'], pair[1]['returned_step']),
                'both_success': all(row['success'] for row in pair),
                'exactly_one_success': sum(bool(row['success']) for row in pair) == 1,
            })
    groups = []
    for name in sorted({case['group'] for case in CASES}):
        subset = [row for row in terminal if row['group'] == name]
        groups.append({
            'group': name, 'terminal': len(subset),
            'originals_passed': sum(row['category'] == 'original_passed' for row in subset),
            'false_alarms': sum(row['category'] == 'original_false_alarm' for row in subset),
            'mutants_exact': sum(row['category'] == 'mutant_exact_step' for row in subset),
            'mutants_passed': sum(row['category'] == 'mutant_passed' for row in subset),
            'mutants_wrong_step': sum(row['category'] == 'mutant_wrong_or_missing_step' for row in subset),
            'unparsed': sum(row['category'] == 'unparsed' for row in subset),
        })
    summary = {
        'route': MODELS[route], 'scheduled_samples': 2 * len(CASES),
        'samples_with_any_request': len(complete),
        'terminal_samples': len(terminal), 'transport_failure_samples': len(observed),
        'pending_samples': sum(row['category'] == 'pending' for row in samples),
        'retry_pending_samples': sum(row['category'] == 'retry_pending' for row in samples),
        'request_attempts': len(requests),
        'retry_attempts': sum(row['attempt'] > 1 for row in requests),
        'timeout_requests': sum(row['exit'] == 'timeout' for row in requests),
        'no_assistant_message_requests': sum(not row['has_assistant_message'] for row in requests),
        'assistant_error_requests': sum(row['stop'] == 'error' for row in requests),
        'nonzero_exit_requests': sum(isinstance(row['exit'], int) and row['exit'] != 0 for row in requests),
        'token_truncated_terminal_samples': sum(row['stop'] == 'length' for row in terminal),
        'primary_success_terminal': sum(row['success'] is True for row in terminal),
        'primary_success_all_scheduled_sensitivity': sum(row['success'] is True for row in samples),
        'originals_passed': count('original_passed'),
        'original_false_alarms': count('original_false_alarm'),
        'mutants_exact_step': count('mutant_exact_step'),
        'mutants_passed_as_correct': count('mutant_passed'),
        'mutants_wrong_or_missing_step': count('mutant_wrong_or_missing_step'),
        'unparsed_terminal_samples': count('unparsed'),
        'original_asserted_steps': [
            {'case': row['case'], 'sample': row['sample'], 'step': row['returned_step']}
            for row in terminal if row['category'] == 'original_false_alarm'
        ],
        'fully_observed_case_pairs': len(pairs),
        'same_verdict_pairs': sum(row['same_verdict'] for row in pairs),
        'same_verdict_and_step_pairs': sum(row['same_verdict_and_step'] for row in pairs),
        'both_samples_successful_cases': sum(row['both_success'] for row in pairs),
        'split_success_cases': sum(row['exactly_one_success'] for row in pairs),
        'request_seconds_sum': round(sum(row['seconds'] for row in requests), 1),
        'terminal_sample_seconds_mean': round(statistics.mean([row['requests_seconds'] for row in terminal]), 2) if terminal else None,
        'input_tokens_all_requests': usage_sum(requests, 'input'),
        'output_tokens_all_requests': usage_sum(requests, 'output'),
        'reasoning_tokens_all_requests': usage_sum(requests, 'reasoning'),
        'total_tokens_all_requests': usage_sum(requests, 'totalTokens'),
        'samples_missing_usage': sum(not row.get('usage') for row in requests),
        'cost_estimate_usd_all_requests': cost_sum(requests),
        'provider_reported_cost_total_requests': sum(bool(((row.get('usage') or {}).get('cost') or {}).get('totalFromProvider')) for row in requests),
        'reported_models': sorted({row.get('reported_model') for row in requests if row.get('reported_model')}),
        'by_group': groups,
    }
    return summary, samples


def main():
    generated = datetime.datetime.now(datetime.timezone.utc).isoformat()
    all_samples = []
    summaries = {}
    for route in MODELS:
        summary, samples = route_summary(route)
        summaries[route] = summary
        all_samples.extend(samples)
    complete = all(s['scheduled_samples'] == s['terminal_samples'] + s['transport_failure_samples'] for s in summaries.values())
    with (HERE / 'samples.csv').open('w', newline='') as f:
        fields = list(all_samples[0])
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_samples)
    def pct(n, d):
        return f'{100*n/d:.1f}%' if d else 'n/a'
    lines = [
        '# Frozen hosted-model referee replay: deterministic report', '',
        f'Generated UTC: {generated}. **Status: {"COMPLETE" if complete else "IN PROGRESS"}.**', '',
        'Primary requested comparison: Kimi, free DeepSeek, and `z-ai/glm-5.3`. FlashX receipts are supplementary and are not substituted for the larger GLM model.', '',
        'Primary success means an original is passed or a planted error is flagged at exactly the planted step. Mathematical categories below count only terminal responses; transport failures and pending samples are separate.', '',
        '| Model | Scheduled | Terminal | Transport-failed | Pending/retrying | Primary success / terminal | Pass+exact / all scheduled (sensitivity) | Unparsed |',
        '|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for route, s in summaries.items():
        lines.append(f"| {s['route']} | {s['scheduled_samples']} | {s['terminal_samples']} | {s['transport_failure_samples']} | {s['pending_samples'] + s['retry_pending_samples']} | {s['primary_success_terminal']}/{s['terminal_samples']} ({pct(s['primary_success_terminal'],s['terminal_samples'])}) | {s['primary_success_all_scheduled_sensitivity']}/{s['scheduled_samples']} | {s['unparsed_terminal_samples']} |")
    lines += ['', '## Mathematical outcomes (terminal responses)', '',
              '| Model | Original passed | Original false alarm | Mutant exact step | Mutant passed | Mutant wrong/missing step | Unparsed |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for route, s in summaries.items():
        lines.append(f"| {s['route']} | {s['originals_passed']} | {s['original_false_alarms']} | {s['mutants_exact_step']} | {s['mutants_passed_as_correct']} | {s['mutants_wrong_or_missing_step']} | {s['unparsed_terminal_samples']} |")
    lines += ['', '## Reliability, repeatability, usage, and cost', '',
              '| Model | Requests | Retries | Timeouts | Token truncations | Fully observed pairs | Same verdict+step | Both samples successful | Input tokens | Output tokens | Reasoning tokens | Cost USD (Pi estimate) |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for route, s in summaries.items():
        lines.append(f"| {s['route']} | {s['request_attempts']} | {s['retry_attempts']} | {s['timeout_requests']} | {s['token_truncated_terminal_samples']} | {s['fully_observed_case_pairs']} | {s['same_verdict_and_step_pairs']} | {s['both_samples_successful_cases']} | {s['input_tokens_all_requests']} | {s['output_tokens_all_requests']} | {s['reasoning_tokens_all_requests']} | ${s['cost_estimate_usd_all_requests']:.6f} |")
    for route, s in summaries.items():
        lines += ['', f"## {s['route']}: proof groups", '',
                  '| Group | Terminal | Original passed | False alarms | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |',
                  '|---|---:|---:|---:|---:|---:|---:|---:|']
        for row in s['by_group']:
            lines.append(f"| {row['group']} | {row['terminal']} | {row['originals_passed']} | {row['false_alarms']} | {row['mutants_exact']} | {row['mutants_passed']} | {row['mutants_wrong_step']} | {row['unparsed']} |")
        lines += ['', f"### {s['route']}: exact per-case outcomes", '',
                  '| Case | Expected | Planted step | Sample 1: category / verdict / returned step | Sample 2: category / verdict / returned step |',
                  '|---|---|---:|---|---|']
        for case in CASES:
            pair = sorted([row for row in all_samples if row['route'] == route and row['case'] == case['id']], key=lambda row: row['sample'])
            cells = []
            for row in pair:
                verdict = row['verdict'] or '—'
                returned = row['returned_step'] or '—'
                cells.append(f"{row['category']} / {verdict} / {returned}")
            lines.append(f"| {case['id']} | {case['expected']} | {case.get('step') or '—'} | {cells[0]} | {cells[1]} |")
    (HERE / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    (HERE / 'summary.json').write_text(json.dumps({'generated_utc': generated, 'complete': complete, 'models': summaries}, indent=2) + '\n')
    primary_routes = ['kimi', 'deepseek', 'glm53']
    primary_complete = all(summaries[route]['scheduled_samples'] == summaries[route]['terminal_samples'] + summaries[route]['transport_failure_samples'] for route in primary_routes)
    frozen_primary_path = HERE / 'primary-summary.json'
    frozen_primary = json.loads(frozen_primary_path.read_text()) if frozen_primary_path.exists() else None
    if primary_complete and frozen_primary and frozen_primary.get('complete'):
        current = {route: summaries[route] for route in primary_routes}
        assert frozen_primary['models'] == current, 'Completed primary results differ from frozen publication'
        # Do not refresh timestamps or let a later model rewrite this completed report.
        print(json.dumps({route: {k: v for k, v in summaries[route].items() if k not in ('by_group', 'original_asserted_steps', 'reported_models')} for route in MODELS}, indent=2))
        return
    (HERE / 'primary-summary.json').write_text(json.dumps({
        'generated_utc': generated, 'complete': primary_complete,
        'models': {route: summaries[route] for route in primary_routes},
    }, indent=2) + '\n')
    primary_lines = [
        '# Requested three-model planted-error referee comparison', '',
        f'Generated UTC: {generated}. **Status: {"COMPLETE" if primary_complete else "IN PROGRESS"}.**', '',
        'Routes: Kimi, free DeepSeek, and larger `z-ai/glm-5.3` (not FlashX). Exact per-case planted/returned steps are in `samples.csv` and `REPORT.md`; all raw responses and failures are in `receipts/`. Mathematical outcomes count terminal responses only.', '',
        '| Model | Terminal /86 | Transport-failed | Pending/retrying | Primary success | Original passed | Original false alarm | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for route in primary_routes:
        s = summaries[route]
        primary_lines.append(f"| {s['route']} | {s['terminal_samples']} | {s['transport_failure_samples']} | {s['pending_samples'] + s['retry_pending_samples']} | {s['primary_success_terminal']}/{s['terminal_samples']} ({pct(s['primary_success_terminal'], s['terminal_samples'])}) | {s['originals_passed']} | {s['original_false_alarms']} | {s['mutants_exact_step']} | {s['mutants_passed_as_correct']} | {s['mutants_wrong_or_missing_step']} | {s['unparsed_terminal_samples']} |")
    primary_lines += ['', '| Model | Requests | Retries | Timeouts | Token truncations | Input tokens | Output tokens | Reasoning tokens | Cost USD (Pi estimate) |',
                      '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for route in primary_routes:
        s = summaries[route]
        primary_lines.append(f"| {s['route']} | {s['request_attempts']} | {s['retry_attempts']} | {s['timeout_requests']} | {s['token_truncated_terminal_samples']} | {s['input_tokens_all_requests']} | {s['output_tokens_all_requests']} | {s['reasoning_tokens_all_requests']} | ${s['cost_estimate_usd_all_requests']:.6f} |")
    (HERE / 'PRIMARY-REPORT.md').write_text('\n'.join(primary_lines) + '\n')
    print(json.dumps({route: {k: v for k, v in s.items() if k not in ('by_group', 'original_asserted_steps', 'reported_models')} for route, s in summaries.items()}, indent=2))


if __name__ == '__main__':
    main()
