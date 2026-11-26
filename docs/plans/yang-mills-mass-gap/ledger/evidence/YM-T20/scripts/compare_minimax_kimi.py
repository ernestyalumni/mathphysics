#!/usr/bin/env python3
"""Compare frozen Kimi receipts with the separately launched MiniMax-M3 run."""
import datetime
import json
from pathlib import Path
import runpy

ROOT = Path('<workspace>')
HERE = ROOT/'<private-run>/full-replay'
LIB = runpy.run_path(str(HERE/'summarize_replay.py'))
FIELDS = ['scheduled_samples','terminal_samples','transport_failure_samples','pending_samples','retry_pending_samples','request_attempts','retry_attempts','timeout_requests','assistant_error_requests','token_truncated_terminal_samples','primary_success_terminal','primary_success_all_scheduled_sensitivity','originals_passed','original_false_alarms','mutants_exact_step','mutants_passed_as_correct','mutants_wrong_or_missing_step','unparsed_terminal_samples','fully_observed_case_pairs','same_verdict_and_step_pairs','both_samples_successful_cases','split_success_cases','terminal_sample_seconds_mean','input_tokens_all_requests','output_tokens_all_requests','reasoning_tokens_all_requests','samples_missing_usage','cost_estimate_usd_all_requests','reported_models']


def pad(summary, samples):
    out = {key: summary[key] for key in FIELDS}
    returned = [x for x in samples if x['terminal']]
    parseable = [x for x in returned if x['verdict'] in ('CORRECT','INCORRECT')]
    out['parseable_samples'] = len(parseable)
    out['success_given_parseable'] = sum(bool(x['success']) for x in parseable)
    out['success_per_all_terminal_percent'] = round(100*out['primary_success_terminal']/len(returned),2) if returned else None
    out['success_per_all_scheduled_percent'] = round(100*out['primary_success_all_scheduled_sensitivity']/out['scheduled_samples'],2)
    out['parseable_percent'] = round(100*len(parseable)/len(returned),2) if returned else None
    out['success_given_parseable_percent'] = round(100*out['success_given_parseable']/len(parseable),2) if parseable else None
    out['cost_per_primary_success_usd'] = out['cost_estimate_usd_all_requests']/out['primary_success_all_scheduled_sensitivity'] if out['primary_success_all_scheduled_sensitivity'] else None
    out['cases'] = []
    for case in LIB['CASES']:
        pair = sorted([x for x in samples if x['case']==case['id']], key=lambda x: x['sample'])
        out['cases'].append({'case':case['id'],'group':case['group'],'expected':case['expected'],'planted_step':case.get('step'),'samples':pair})
    return out


def main():
    summaries = {}
    for route in ('kimi','minimax'):
        summary, samples = LIB['route_summary'](route)
        if summary['scheduled_samples'] != summary['terminal_samples']+summary['transport_failure_samples']:
            raise SystemExit(f'{route} is not fully resolved; not writing final comparison')
        summaries[route] = pad(summary, samples)
    indexed = {route:{(s['case'],s['sample']):s for s in d['cases'] for s in s['samples']} for route,d in summaries.items()}
    keys = sorted(set(indexed['kimi']) & set(indexed['minimax']))
    matched = {}
    for name, predicate in [
        ('both_terminal', lambda m,k: m['terminal'] and k['terminal']),
        ('both_parseable', lambda m,k: m['verdict'] in ('CORRECT','INCORRECT') and k['verdict'] in ('CORRECT','INCORRECT')),
    ]:
        chosen = [key for key in keys if predicate(indexed['minimax'][key],indexed['kimi'][key])]
        matched[name] = {'observations':len(chosen),'minimax_success':sum(bool(indexed['minimax'][key]['success']) for key in chosen),'kimi_success':sum(bool(indexed['kimi'][key]['success']) for key in chosen),'discordant_minimax_only':[key for key in chosen if indexed['minimax'][key]['success'] and not indexed['kimi'][key]['success']],'discordant_kimi_only':[key for key in chosen if indexed['kimi'][key]['success'] and not indexed['minimax'][key]['success']], 'selection_warning':'post-hoc available-response subset; not an independent-case superiority test'}
    generated = datetime.datetime.now(datetime.timezone.utc).isoformat()
    output = {'generated_utc':generated,'complete':True,'models':summaries,'matched_observation_sensitivity':matched,'scope':'MiniMax-M3 (the only available MiniMax route) versus Kimi on the frozen 43 cases, two samples each. These are not fresh independently human-adjudicated labels and model agreement is not theorem certification.'}
    (HERE/'minimax-comparison.json').write_text(json.dumps(output,indent=2)+'\n')
    m,k = summaries['minimax'],summaries['kimi']
    lines = ['# MiniMax-M3 versus Kimi — frozen 43-case replay','',f'Generated UTC: {generated}. Each route has 86 scheduled samples; no expected label is supplied to the model. MiniMax is the only configured MiniMax route. MiniMax started later, so provider load and exact provider-side output caps are not controlled.','', 'Primary success means an original passed or a mutant flagged at exactly its planted step. Mathematical rates count terminal responses; all-scheduled sensitivity also retains transport failures.','','| Metric | MiniMax-M3 | Kimi |','|---|---:|---:|']
    metrics = [('Terminal samples','terminal_samples'),('Transport-failed samples','transport_failure_samples'),('Requests / retries',None),('Timeouts','timeout_requests'),('Unparsed terminal answers','unparsed_terminal_samples'),('Token truncations','token_truncated_terminal_samples'),('Primary success / terminal','primary_success_terminal'),('Success / all 86 scheduled','primary_success_all_scheduled_sensitivity'),('Parseable verdicts','parseable_samples'),('Success given parseable','success_given_parseable'),('Originals passed /28','originals_passed'),('Original false alarms','original_false_alarms'),('Exact planted step /58','mutants_exact_step'),('Mutants passed','mutants_passed_as_correct'),('Wrong/missing step','mutants_wrong_or_missing_step'),('Same verdict+step pairs','same_verdict_and_step_pairs'),('Both-sample-success cases','both_samples_successful_cases'),('Mean seconds per terminal sample','terminal_sample_seconds_mean'),('Input tokens, all requests','input_tokens_all_requests'),('Output tokens, all requests','output_tokens_all_requests'),('Reasoning tokens, all requests','reasoning_tokens_all_requests'),('Cost estimate USD','cost_estimate_usd_all_requests'),('Cost estimate per success USD','cost_per_primary_success_usd')]
    for label,key in metrics:
        if key is None:
            a,b = f"{m['request_attempts']} / {m['retry_attempts']}",f"{k['request_attempts']} / {k['retry_attempts']}"
        else:
            a,b=m[key],k[key]
            if key in ('cost_estimate_usd_all_requests','cost_per_primary_success_usd'):
                a,b=f'${a:.6f}' if a is not None else 'n/a',f'${b:.6f}' if b is not None else 'n/a'
            if key=='primary_success_terminal': a,b=f'{a}/{m["terminal_samples"]}',f'{b}/{k["terminal_samples"]}'
        lines.append(f'| {label} | {a} | {b} |')
    for title,key in [('Matched terminal observations','both_terminal'),('Matched parseable observations','both_parseable')]:
        x=matched[key]
        lines += ['',f'## {title}', '',f"Secondary, post-hoc subset: {x['observations']} observations; MiniMax {x['minimax_success']}/{x['observations']}, Kimi {x['kimi_success']}/{x['observations']}.", 'Discordances are recorded in `minimax-comparison.json`; these post-hoc common-response samples do not support an independent-case superiority claim.', '']
    lines += ['## Error and research boundary','', 'Receipt-level mathematical errors remain planted-label outcomes, not proof that each historical original has independent human certification. Exact per-case categories and planted/returned steps are reproduced below and in `samples.csv`. Different elapsed times are cumulative request time, not guaranteed wall-clock model latency. Usage/cost is Pi/TokenRouter accounting, not a provider invoice.','', 'A double model pass is not a theorem certificate. Use every actual Yang–Mills flag only as input to independent strategist/human or Lean adjudication.','', '## Exact MiniMax per-case comparison', '', '| Case | Expected | Planted step | MiniMax sample 1: category / verdict / step | MiniMax sample 2 | Kimi sample 1 | Kimi sample 2 |','|---|---|---:|---|---|---|---|']
    for case in LIB['CASES']:
        def cell(route,sample): 
            x=indexed[route][(case['id'],sample)]
            return f"{x['category']} / {x['verdict'] or '—'} / {x['returned_step'] or '—'}"
        lines.append(f"| {case['id']} | {case['expected']} | {case.get('step') or '—'} | {cell('minimax',1)} | {cell('minimax',2)} | {cell('kimi',1)} | {cell('kimi',2)} |")
    lines.append('')
    (HERE/'MINIMAX-COMPARISON.md').write_text('\n'.join(lines))
    print(json.dumps({'complete':True,'terminal_and_success':{route:{'terminal':d['terminal_samples'],'primary':d['primary_success_terminal'],'all_scheduled':d['primary_success_all_scheduled_sensitivity'],'parseable':d['parseable_samples'],'unparsed':d['unparsed_terminal_samples']} for route,d in summaries.items()}},indent=2))


if __name__=='__main__':
    main()
