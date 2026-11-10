#!/usr/bin/env python3
"""Huzaifa Paper 1: auditable orchestration and standalone HTML defense report.

Standard-library-only. Does not reconstruct unknown original experiment code.
Unmapped or missing stages BLOCK a claim of complete reproducibility.
"""
from __future__ import annotations

import argparse
import base64
import csv
import datetime as dt
import hashlib
import html
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import traceback
import zipfile

VERSION = '1.0.0'
REPORT_CSS = """
:root{--navy:#102b44;--teal:#087f8c;--ink:#203349;--light:#f2f6fa;--border:#d8e2ec;--red:#a32929;--gold:#8b5900}
*{box-sizing:border-box}body{font:15px/1.58 system-ui,-apple-system,Segoe UI,Arial,sans-serif;color:var(--ink);margin:0;background:#f5f7fa}
header{background:var(--navy);color:white;padding:32px max(25px,calc((100% - 1120px)/2))}header h1{font-size:28px;margin:4px 0 8px}header p{margin:0;color:#d7e7f2}
main{max-width:1120px;margin:26px auto;padding:0 18px 45px}section{background:white;border:1px solid var(--border);border-radius:10px;margin:18px 0;padding:23px;box-shadow:0 2px 8px #132a4209}
h2{font-size:20px;color:var(--navy);margin:0 0 13px}h3{font-size:16px;color:var(--navy)}.status{display:inline-block;padding:5px 12px;border-radius:16px;font-weight:700;font-size:13px;background:#e3edf2;color:var(--navy)}.ok{background:#e4f3e9;color:#196136}.bad{background:#ffe6e5;color:var(--red)}.warn{background:#fff1d6;color:var(--gold)}
table{border-collapse:collapse;width:100%;font-size:13px;margin:13px 0}th,td{border-bottom:1px solid var(--border);padding:9px 11px;text-align:left;vertical-align:top}th{background:#edf4f8;color:var(--navy)}tr:nth-child(even) td{background:#fafcfe}
small,.muted{color:#5b6d7d}.metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(185px,1fr));gap:12px}.metric{border:1px solid var(--border);border-radius:8px;padding:12px}.metric strong{font-size:23px;display:block;color:var(--navy)}.metric span{font-size:12px;color:#536a7a}
code,pre{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:12px;overflow-wrap:anywhere}pre{white-space:pre-wrap;background:#f3f6f8;padding:12px;border-radius:6px}figure{margin:20px auto;text-align:center}figure img{max-width:100%;max-height:560px;object-fit:contain}figcaption{font-size:12px;color:#526778}.alert{padding:12px;border-left:4px solid #d4931b;background:#fff9e8}.foot{font-size:12px;color:#617181}
@media print{body{background:white}section{box-shadow:none;break-inside:avoid}header{padding:20px;background:white;color:#102b44;border-bottom:2px solid #102b44}header p{color:#536a7a}}
"""


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def dump_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=False, ensure_ascii=False, allow_nan=False, default=str), encoding='utf-8')


def safe_relative(project, relative):
    if not isinstance(relative, str) or not relative.strip():
        raise ValueError('Empty path')
    p = Path(relative)
    if p.is_absolute() or '..' in p.parts or re.match(r'^[A-Za-z]:', relative):
        raise ValueError('Only project-relative paths are allowed: ' + relative)
    target = (project / p).resolve()
    if not target.is_relative_to(project.resolve()):
        raise ValueError('Path escapes project: ' + relative)
    return target


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def file_record(path, project, hash_limit=400 * 1024 * 1024):
    try:
        stat = path.stat()
        rel = str(path.relative_to(project)).replace('\\', '/')
        rec = {'path': rel, 'bytes': stat.st_size, 'mtime_utc': dt.datetime.fromtimestamp(stat.st_mtime, dt.timezone.utc).isoformat()}
        if stat.st_size <= hash_limit:
            rec['sha256'] = sha256(path)
        else:
            rec['sha256'] = None
            rec['hash_note'] = 'Large file: metadata only; not cryptographically verified'
        return rec
    except Exception as exc:
        return {'path': str(path), 'error': str(exc)}


def collect_inventory(project, pipeline):
    src = project / 'src'
    scripts = sorted(src.glob('*.py')) if src.is_dir() else []
    stages = []
    for stage in pipeline['stages']:
        hint = stage.get('hints', [])
        scored = []
        for script in scripts:
            name = script.stem.lower()
            score = sum(2 if h.lower() in name else 0 for h in hint)
            try:
                content = script.read_text(encoding='utf-8', errors='replace')
            except OSError:
                content = ''
            for output in stage.get('expected_outputs', []):
                if Path(output).name in content:
                    score += 3
                if Path(output).parent.name in content:
                    score += 2
            if score:
                scored.append({'script': str(script.relative_to(project)).replace('\\', '/'), 'score': score})
        scored.sort(key=lambda x: (-x['score'], x['script']))
        configured = stage.get('script', '')
        stages.append({
            'id': stage['id'], 'label': stage['label'],
            'configured_script': configured,
            'configured_exists': bool(configured and safe_relative(project, configured).is_file()),
            'suggestions': scored[:7],
            'output_presence': [{'path': x, 'exists': safe_relative(project, x).is_file()} for x in stage.get('expected_outputs', [])]
        })
    inputs = []
    for folder in ('data/raw', 'data/processed', 'results/step31_genomator', 'results/step31_genomator_1024'):
        p = project / folder
        if p.exists():
            items = sorted(x for x in p.iterdir() if x.is_file())
            inputs.append({'folder': folder, 'file_count': len(items), 'sample': [{'name': x.name, 'bytes': x.stat().st_size} for x in items[:35]]})
        else:
            inputs.append({'folder': folder, 'missing': True})
    return {'project': str(project), 'source_scripts': [file_record(x, project, 20 * 1024 * 1024) for x in scripts], 'stages': stages, 'inputs': inputs}


def environment(project, pipeline):
    info = {'captured_utc': utc_now(), 'runner_version': VERSION, 'os': platform.platform(), 'python': sys.version,
            'project': str(project), 'cwd': os.getcwd(), 'wsl_distro': 'Ubuntu', 'pipeline_schema': pipeline.get('schema_version')}
    packages = {}
    for name in ('numpy', 'pandas', 'scipy', 'scikit-learn', 'matplotlib', 'cyvcf2', 'vcfpy', 'msprime', 'tskit', 'genomator'):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = 'not available in runner environment'
    info['runner_packages'] = packages
    for label, cmd in (
        ('windows_python', [str(project / '.venv' / 'Scripts' / 'python.exe'), '--version']),
        ('wsl_python', ['wsl', '-d', 'Ubuntu', '--', '/home/huzaifa/genomator-venv/bin/python', '--version']),
        ('wsl_genomator', ['wsl', '-d', 'Ubuntu', '--', '/home/huzaifa/genomator-venv/bin/genomator', '--help']),
    ):
        try:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=12)
            info[label] = {'returncode': p.returncode, 'output': (p.stdout + p.stderr)[:1200]}
        except Exception as exc:
            info[label] = {'error': str(exc)}
    # Capture package versions from the actual WSL Genomator environment, not
    # only the Windows orchestrator. No pip executable is required.
    wsl_code = (
        'import importlib.metadata as m, json\n'
        'names = ["numpy","pandas","scipy","scikit-learn","matplotlib","cyvcf2","vcfpy","msprime","tskit","genomator"]\n'
        'out = {}\n'
        'for n in names:\n'
        '    try: out[n] = m.version(n)\n'
        '    except m.PackageNotFoundError: out[n] = "not installed"\n'
        'print(json.dumps(out))\n'
    )
    try:
        p = subprocess.run(['wsl','-d','Ubuntu','--','/home/huzaifa/genomator-venv/bin/python','-c',wsl_code],
                           capture_output=True, text=True, timeout=30)
        info['wsl_package_versions'] = json.loads(p.stdout) if p.returncode == 0 else {'error': p.stderr[:1500]}
    except Exception as exc:
        info['wsl_package_versions'] = {'error': str(exc)}
    try:
        p = subprocess.run(['git','-C',str(project),'rev-parse','HEAD'],capture_output=True,text=True,timeout=8)
        info['project_git_commit'] = p.stdout.strip() if p.returncode == 0 else None
    except Exception:
        info['project_git_commit'] = None
    return info


def to_wsl_path(path):
    s = str(path)
    match = re.match(r'^([A-Za-z]):[\\/](.*)$', s)
    if not match:
        raise RuntimeError('WSL executor requires Windows drive path: ' + s)
    return '/mnt/' + match.group(1).lower() + '/' + match.group(2).replace('\\', '/')


def command_for(stage, project, python_win, python_wsl, distro):
    script = safe_relative(project, stage['script'])
    executor = stage.get('executor', 'wsl_python')
    if executor == 'wsl_python':
        return ['wsl', '-d', distro, '--', python_wsl, to_wsl_path(script)]
    if executor == 'windows_python':
        return [str(python_win), str(script)]
    raise ValueError('Unsupported executor: ' + executor)


def validate_pipeline(project, pipeline):
    errors = []
    for stage in pipeline['stages']:
        script = stage.get('script', '').strip()
        if not script:
            errors.append(f"{stage['id']}: original script NOT MAPPED")
        else:
            try:
                if not safe_relative(project, script).is_file():
                    errors.append(f"{stage['id']}: script missing: {script}")
            except Exception as exc:
                errors.append(f"{stage['id']}: invalid script: {exc}")
        if not stage.get('expected_outputs'):
            errors.append(f"{stage['id']}: no expected outputs configured")
        for output in stage.get('expected_outputs', []):
            try:
                safe_relative(project, output)
            except Exception as exc:
                errors.append(f"{stage['id']}: invalid output path: {exc}")
    if not (project / 'data').exists():
        errors.append('data/ missing: original genomic input data unavailable')
    return errors


def run_stage(stage, project, out, args):
    script = safe_relative(project, stage['script'])
    log_path = out / 'logs' / (stage['id'] + '.log')
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = command_for(stage, project, args.python_win, args.python_wsl, args.distro)
    before = {}
    for item in stage['expected_outputs']:
        p = safe_relative(project, item)
        before[item] = (p.stat().st_mtime_ns, p.stat().st_size) if p.is_file() else None
    start_ns = time.time_ns()
    start = utc_now()
    started = time.monotonic()
    rc = None
    error = None
    with log_path.open('w', encoding='utf-8') as log:
        log.write('STAGE: ' + stage['id'] + '\nSTART: ' + start + '\nSCRIPT SHA256: ' + sha256(script) + '\n')
        log.write('COMMAND: ' + ' '.join(cmd) + '\n\n')
        log.flush()
        try:
            env = os.environ.copy()
            env.update({'PYTHONUNBUFFERED': '1', 'PYTHONHASHSEED': '0', 'MPLBACKEND': 'Agg', 'REPRO_RUN_ID': out.name})
            with subprocess.Popen(cmd, cwd=str(project), stdout=log, stderr=subprocess.STDOUT, env=env) as proc:
                try:
                    rc = proc.wait(timeout=int(stage.get('timeout_seconds', 7200)))
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                    rc = -9
                    error = 'Timed out; process killed'
        except Exception as exc:
            error = str(exc)
            log.write('\nRUNNER ERROR: ' + error + '\n')
    changed = []
    missing = []
    unchanged = []
    output_records = []
    for item in stage['expected_outputs']:
        p = safe_relative(project, item)
        if not p.is_file() or p.stat().st_size == 0:
            missing.append(item)
        else:
            output_records.append(file_record(p, project))
            current = (p.stat().st_mtime_ns, p.stat().st_size)
            if before[item] != current and p.stat().st_mtime_ns >= start_ns - 3_000_000_000:
                changed.append(item)
            else:
                unchanged.append(item)
    ok = rc == 0 and not missing and not unchanged
    return {'id': stage['id'], 'label': stage['label'], 'script': stage['script'], 'executor': stage['executor'],
            'command': cmd, 'started_utc': start, 'finished_utc': utc_now(), 'duration_seconds': round(time.monotonic() - started, 2),
            'returncode': rc, 'error': error, 'status': 'PASS' if ok else 'FAIL', 'missing_outputs': missing,
            'unchanged_outputs': unchanged, 'changed_outputs': changed, 'outputs': output_records,
            'log': str(log_path.relative_to(out)).replace('\\', '/')}


def as_number(v):
    try:
        value = float(v)
        if value != value or abs(value) == float('inf'):
            return None
        return value
    except (ValueError, TypeError):
        return None


def load_metrics(project):
    p = project / 'results' / 'step36_final_validation' / 'final_validation.json'
    if not p.is_file():
        return {'source': str(p), 'available': False, 'experiments': [], 'sensitivity': []}
    try:
        data = json.loads(p.read_text(encoding='utf-8'))
        experiments = []
        for row in data.get('experiments', []):
            boot = row.get('bootstrap') or {}
            experiments.append({
                'experiment': row.get('experiment'),
                'actual_auc': as_number(row.get('actual_auc')),
                'shadow_auc': as_number(row.get('shadow_auc')),
                'auc_gap': as_number(row.get('auc_gap')),
                'calibrated_auc': as_number(row.get('calibrated_auc')),
                'gap_ci': [as_number(x) for x in boot.get('gap_ci', [])],
                'snps': row.get('snps'), 'training_n': row.get('training_n'), 'heldout_n': row.get('heldout_n')
            })
        sensitivity = []
        for row in data.get('relatedness_sensitivity', []):
            sensitivity.append({k: as_number(v) for k, v in row.items()})
        result = {'source': str(p), 'available': True, 'experiments': experiments, 'sensitivity': sensitivity}
        repeats = project / 'results' / 'step33_independent_validation' / 'repeat_results.csv'
        if repeats.exists():
            with repeats.open(newline='', encoding='utf-8-sig') as f:
                rows = list(csv.DictReader(f))
            result['repeat_file'] = str(repeats)
            result['repeat_count'] = len(rows)
            result['repeat_columns'] = list(rows[0]) if rows else []
            actual_key = next((k for k in (result['repeat_columns']) if 'actual' in k.lower() and 'auc' in k.lower()), None)
            shadow_key = next((k for k in (result['repeat_columns']) if 'shadow' in k.lower() and 'auc' in k.lower()), None)
            if actual_key and shadow_key:
                a = [as_number(x[actual_key]) for x in rows]
                b = [as_number(x[shadow_key]) for x in rows]
                if all(v is not None for v in a + b) and a:
                    result['repeats'] = {'n_repeats': len(a), 'actual_auc_mean': sum(a) / len(a),
                                         'shadow_auc_mean': sum(b) / len(b), 'auc_gap_mean': (sum(a) - sum(b)) / len(a),
                                         'actual_aucs': a, 'shadow_aucs': b}
        return result
    except Exception as exc:
        return {'source': str(p), 'available': False, 'error': str(exc), 'experiments': [], 'sensitivity': []}


def compare(metrics, pipeline, is_fresh):
    comparisons = []
    reference = pipeline.get('reference_results', {})
    tolerance = float(pipeline.get('baseline_tolerance_auc', 0.00015))
    for exp in metrics.get('experiments', []):
        key = exp.get('experiment')
        ref = reference.get(key, {})
        for field in ('actual_auc', 'shadow_auc', 'auc_gap', 'calibrated_auc'):
            observed = exp.get(field)
            expected = ref.get(field)
            if observed is None or expected is None:
                continue
            delta = observed - expected
            comparisons.append({'experiment': key, 'metric': field, 'reference': expected, 'observed': observed,
                                'difference': delta, 'within_tolerance': abs(delta) <= tolerance, 'fresh_rerun': is_fresh})
        ci = exp.get('gap_ci') or []
        for index, field in enumerate(('gap_ci_lower', 'gap_ci_upper')):
            if len(ci) > index and ci[index] is not None and field in ref:
                delta = ci[index] - ref[field]
                comparisons.append({'experiment': key, 'metric': field, 'reference': ref[field], 'observed': ci[index],
                                    'difference': delta, 'within_tolerance': abs(delta) <= tolerance, 'fresh_rerun': is_fresh})
    rep = metrics.get('repeats', {})
    ref = reference.get('256_repeats', {})
    for field in ('n_repeats', 'actual_auc_mean', 'shadow_auc_mean', 'auc_gap_mean'):
        if field in rep and field in ref:
            delta = rep[field] - ref[field]
            comparisons.append({'experiment': '256_repeats', 'metric': field, 'reference': ref[field], 'observed': rep[field],
                                'difference': delta, 'within_tolerance': abs(delta) <= (0 if field == 'n_repeats' else tolerance),
                                'fresh_rerun': is_fresh})
    return comparisons


def fingerprint_inputs(project, pipeline):
    files = []
    for stage in pipeline['stages']:
        script = stage.get('script')
        if script:
            p = safe_relative(project, script)
            if p.is_file():
                files.append(p)
    for folder in ('data/raw', 'data/processed'):
        p = project / folder
        if p.is_dir():
            files.extend(x for x in p.rglob('*') if x.is_file() and not x.is_symlink())
    unique = sorted(set(files), key=lambda x: str(x))
    return [file_record(x, project) for x in unique[:1000]]


def esc(v):
    return html.escape(str(v), quote=True)


def format_float(x, digits=5):
    return f'{x:.{digits}f}' if isinstance(x, (int, float)) else '—'


def table(headers, rows):
    out = '<table><thead><tr>' + ''.join('<th>' + esc(x) + '</th>' for x in headers) + '</tr></thead><tbody>'
    for row in rows:
        out += '<tr>' + ''.join('<td>' + esc(x) + '</td>' for x in row) + '</tr>'
    return out + '</tbody></table>'


def figures_html(project, full_fresh):
    folder = project / 'results' / 'step35_publication' / 'figures'
    if not folder.is_dir():
        return '<p>No figure directory found.</p>'
    items = sorted(folder.glob('*.png'))[:8]
    if not items:
        return '<p>No PNG figures found.</p>'
    tag = 'Regenerated in this run' if full_fresh else 'Existing figures — not verified as regenerated'
    result = '<p class="muted">' + esc(tag) + '</p>'
    for item in items:
        if item.stat().st_size > 8 * 1024 * 1024:
            result += '<p>Image omitted because it exceeds 8 MB: ' + esc(item.name) + '</p>'
            continue
        b64 = base64.b64encode(item.read_bytes()).decode('ascii')
        result += f'<figure><img src="data:image/png;base64,{b64}" alt="{esc(item.name)}"><figcaption>{esc(item.name)}</figcaption></figure>'
    return result


def build_html(project, run_id, mode, status, reason, stage_results, inventory, metrics, comparisons, env, manifest, fresh):
    total = len(stage_results)
    passed = sum(s['status'] == 'PASS' for s in stage_results)
    cls = 'ok' if status == 'REPRODUCED' else ('bad' if status == 'FAILED' else 'warn')
    summary = []
    for exp in metrics.get('experiments', []):
        summary.append([exp.get('experiment'), format_float(exp.get('actual_auc')), format_float(exp.get('shadow_auc')),
                        format_float(exp.get('auc_gap')), format_float(exp.get('calibrated_auc')),
                        ' to '.join(format_float(x) for x in exp.get('gap_ci', []))])
    stages = []
    for s in stage_results:
        stages.append([s['id'], s['status'], s.get('script', ''), str(s.get('returncode')), format_float(s.get('duration_seconds'), 1),
                       '; '.join(s.get('missing_outputs', []) + s.get('unchanged_outputs', [])) or 'all refreshed'])
    cmp_rows = [[c['experiment'], c['metric'], format_float(c['reference']), format_float(c['observed']),
                 format_float(c['difference'], 6), 'PASS' if c['within_tolerance'] else 'MISMATCH'] for c in comparisons]
    stage_suggestions = []
    for item in inventory['stages']:
        best = ', '.join(f"{s['script']} ({s['score']})" for s in item['suggestions'][:3])
        stage_suggestions.append([item['id'], item['configured_script'] or 'UNMAPPED', 'yes' if item['configured_exists'] else 'no', best or 'none'])
    sensitivity = []
    for s in metrics.get('sensitivity', []):
        sensitivity.append([format_float(s.get('kinship_threshold'), 4), format_float(s.get('flagged_heldout'), 0),
                            format_float(s.get('remaining_heldout'), 0), format_float(s.get('actual_auc')),
                            format_float(s.get('shadow_auc')), format_float(s.get('auc_gap'))])
    warnings = [
        'The attack is nearest-synthetic-neighbor genotype mismatch, not DOMIAS.',
        'A shadow-control AUC near 0.5 does not by itself prove that ancestry has been calibrated away.',
        'These results do not establish direct copying, differential privacy, or statistically significant ancestry disparities.',
        'Kinship sensitivity used common variants on chromosomes 20–22, not genome-wide pedigree confirmation.',
        'The 1024-SNP result has one cohort split; the 256- and 1024-SNP panels are not a causal size experiment.',
        'Rerun reproducibility requires exact inputs, original code, dependencies, fixed seeds, and comparable Genomator settings.'
    ]
    content = [f'''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Huzaifa Paper 1 — Reproducibility Report</title><style>{REPORT_CSS}</style></head><body>
<header><p>RESEARCH DEFENSE · COMPUTATIONAL REPRODUCIBILITY</p><h1>When Ancestry Looks Like Memorization</h1><p>Ancestry-matched shadow controls for synthetic-genome membership inference · Run {esc(run_id)}</p></header><main>
<section><h2>Reproducibility verdict</h2><p><span class="status {cls}">{esc(status)}</span></p><p class="alert">{esc(reason)}</p><div class="metric-grid"><div class="metric"><strong>{passed}/{total}</strong><span>Executed stages passed</span></div><div class="metric"><strong>{len(comparisons)}</strong><span>Numeric checks</span></div><div class="metric"><strong>{'YES' if fresh else 'NO'}</strong><span>Fresh complete execution</span></div></div><p class="muted">Mode: {esc(mode)} · Generated UTC: {esc(utc_now())} · Project: {esc(project)}</p></section>
<section><h2>Membership-inference results</h2><p>Metrics are parsed from <code>results/step36_final_validation/final_validation.json</code>. Values from prior outputs are marked historical when no full rerun occurred.</p>
{table(['Panel','Actual AUC','Shadow AUC','ΔAUC','Adjusted AUC','95% bootstrap CI for ΔAUC'],summary) if summary else '<p>Metrics unavailable. No values have been fabricated.</p>'}
<h3>Kinship sensitivity (1024-SNP panel)</h3>{table(['Threshold','Flagged','Held-out remaining','Actual AUC','Shadow AUC','ΔAUC'],sensitivity) if sensitivity else '<p>Not available.</p>'}</section>
<section><h2>Original-result comparison</h2><p>Reference values are fixed from the earlier experimental logs. A match is meaningful as a fresh reproduction only if all required stages were executed successfully.</p>{table(['Panel','Metric','Reference','Observed','Difference','Tolerance check'],cmp_rows) if cmp_rows else '<p>No comparable metrics.</p>'}</section>
<section><h2>Pipeline execution and audit trail</h2>{table(['Stage','Status','Original script','Exit','Seconds','Output freshness / issues'],stages) if stages else '<p>No stages executed in this mode.</p>'}</section>
<section><h2>Script mapping and discovery</h2><p>Suggestions are heuristic; review before entering scripts in <code>reproduce/pipeline.json</code>. Do not execute unreviewed candidates.</p>{table(['Stage','Configured script','Exists','Candidate scripts (ranked)'],stage_suggestions)}</section>
<section><h2>Figures</h2>{figures_html(project, fresh)}</section>
<section><h2>Interpretation and boundaries</h2><ul>{''.join('<li>'+esc(w)+'</li>' for w in warnings)}</ul><p><strong>Defensible claim:</strong> In evaluated Genomator experiments, the actual-generator nearest-neighbor attack discriminated membership better than an independent ancestry-matched shadow-generator control. This is evidence of a membership-associated signal under the tested setup, not proof of its precise mechanism.</p></section>
<section><h2>Environment and input integrity</h2><p>Complete machine-readable records are saved beside this HTML report as <code>environment.json</code>, <code>input_manifest.json</code>, <code>pipeline_used.json</code>, <code>stage_results.json</code>, and <code>comparison.csv</code>.</p><p>Runner Python: {esc(env.get('python','').splitlines()[0])} · Platform: {esc(env.get('os'))} · Files fingerprinted: {len(manifest)}</p><p>Files over 400 MB are metadata-only and cannot be claimed cryptographically verified. Scripts are SHA-256 hashed.</p></section>
<section><h2>Reproduction protocol</h2><ol><li>Keep the original project data and original scripts in the project directory.</li><li>Map every required stage to its exact original source script.</li><li>Run <code>--mode full --allow-overwrite</code> using the configured environment. Original experiment scripts may overwrite existing results.</li><li>Inspect every stage log, output freshness check, metric comparison, and input manifest.</li><li>Retain this HTML together with the JSON and CSV audit trail for your defense.</li></ol></section>
<p class="foot">Automated verification is not a substitute for code review, ethical data governance, independent replication, or peer review. Generated with Huzaifa Paper 1 Reproduction Runner {esc(VERSION)}.</p></main></body></html>''']
    return ''.join(content)


def snapshot_code(project, out, pipeline):
    """Archive available research source code, never genotype data."""
    folder = out / 'code_snapshot'
    count = 0
    src = project / 'src'
    if src.is_dir():
        for path in sorted(src.rglob('*.py')):
            if path.is_symlink() or path.stat().st_size > 5 * 1024 * 1024:
                continue
            target = folder / path.relative_to(project)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            count += 1
    dump_json(folder / 'source_snapshot_manifest.json', {'scripts_copied': count, 'note': 'Data and Genomator installed source not included; only project src/*.py files.'})
    return count


def defense_zip(out):
    """Portable code + report + log bundle; never include patient/genotype data."""
    destination = out / 'defense_bundle.zip'
    with zipfile.ZipFile(destination, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.rglob('*')):
            if p.is_file() and p != destination and p.stat().st_size < 30 * 1024 * 1024:
                z.write(p, arcname=p.relative_to(out))
    return destination


def main():
    parser = argparse.ArgumentParser(description='Huzaifa Paper 1 reproduction and defense report')
    parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--pipeline', type=Path, default=None)
    parser.add_argument('--mode', choices=('inventory', 'full', 'report'), default='inventory')
    parser.add_argument('--allow-overwrite', action='store_true', help='Required for full mode: original scripts may overwrite outputs')
    parser.add_argument('--python-win', type=Path, default=None)
    parser.add_argument('--python-wsl', default='/home/huzaifa/genomator-venv/bin/python')
    parser.add_argument('--distro', default='Ubuntu')
    args = parser.parse_args()
    project = args.project.resolve()
    args.python_win = args.python_win or (project / '.venv' / 'Scripts' / 'python.exe')
    config_path = args.pipeline or (Path(__file__).resolve().parent / 'pipeline.json')
    pipeline = json.loads(config_path.read_text(encoding='utf-8'))
    if pipeline.get('schema_version') != 1:
        raise ValueError('Unsupported pipeline schema')
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = project / 'results' / 'repro_runs' / (stamp + '_' + args.mode)
    out.mkdir(parents=True, exist_ok=True)
    inventory = collect_inventory(project, pipeline)
    dump_json(out / 'inventory.json', inventory)
    dump_json(out / 'pipeline_used.json', pipeline)
    env = environment(project, pipeline)
    dump_json(out / 'environment.json', env)
    manifest = fingerprint_inputs(project, pipeline)
    dump_json(out / 'input_manifest.json', manifest)
    snapshot_code(project, out, pipeline)
    stage_results = []
    errors = []
    status = 'INVENTORY ONLY'
    reason = 'This run inventories the original scripts and data; it does not execute experiments or establish reproducibility.'
    if args.mode == 'full':
        errors = validate_pipeline(project, pipeline)
        if not args.allow_overwrite:
            errors.append('Full mode requires --allow-overwrite because original scripts may overwrite prior results.')
        if errors:
            status = 'BLOCKED'
            reason = 'Full rerun blocked before executing scripts: ' + '; '.join(errors)
        else:
            status = 'REPRODUCED'
            reason = 'All required stages ran, outputs refreshed, and comparisons will be evaluated.'
            for stage in pipeline['stages']:
                print('RUNNING', stage['id'], stage['label'], flush=True)
                result = run_stage(stage, project, out, args)
                stage_results.append(result)
                dump_json(out / 'stage_results.json', stage_results)
                print(' ', result['status'], f"({result['duration_seconds']} s)", flush=True)
                if result['status'] != 'PASS':
                    status = 'FAILED'
                    reason = 'Execution stopped at ' + stage['id'] + ': ' + (result.get('error') or ', '.join(result['missing_outputs'] + result['unchanged_outputs']) or f"exit code {result['returncode']}")
                    break
    elif args.mode == 'report':
        status = 'HISTORICAL REPORT ONLY'
        reason = 'Existing output files were read; no generator, attack, or validation stage was rerun.'
    dump_json(out / 'stage_results.json', stage_results)
    metrics = load_metrics(project)
    fresh = status == 'REPRODUCED' and len(stage_results) == len(pipeline['stages']) and all(x['status'] == 'PASS' for x in stage_results)
    comparisons = compare(metrics, pipeline, fresh)
    if fresh and (not metrics.get('available') or not comparisons):
        status = 'FAILED'
        reason = 'Fresh execution finished but numeric results could not be parsed and compared.'
        fresh = False
    if fresh and any(not c['within_tolerance'] for c in comparisons):
        status = 'MISMATCH'
        reason = 'All stages reran but one or more metrics differ from the locked reference values beyond tolerance. Review randomness and logs.'
        fresh = False
    if fresh and not any(c['experiment'] == '256_repeats' for c in comparisons):
        status = 'INCOMPLETE VERIFICATION'
        reason = 'Stage execution succeeded, but five-repeat summary could not be independently parsed from the repeat CSV.'
        fresh = False
    dump_json(out / 'metrics.json', metrics)
    with (out / 'comparison.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['experiment','metric','reference','observed','difference','within_tolerance','fresh_rerun'])
        writer.writeheader()
        writer.writerows(comparisons)
    dump_json(out / 'verdict.json', {'status': status, 'reason': reason, 'fresh_full_reproduction': fresh, 'mode': args.mode,
                                    'run_id': out.name, 'stages_executed': len(stage_results), 'stage_passes': sum(s['status'] == 'PASS' for s in stage_results),
                                    'comparison_count': len(comparisons), 'preflight_errors': errors})
    page = build_html(project, out.name, args.mode, status, reason, stage_results, inventory, metrics, comparisons, env, manifest, fresh)
    report_path = out / 'report.html'
    report_path.write_text(page, encoding='utf-8')
    latest = project / 'results' / 'repro_runs' / 'latest.html'
    shutil.copy2(report_path, latest)
    bundle = defense_zip(out)
    print('\nVERDICT:', status)
    print('REASON:', reason)
    print('HTML REPORT:', report_path)
    print('LATEST HTML:', latest)
    print('INVENTORY:', out / 'inventory.json')
    print('DEFENSE BUNDLE:', bundle)
    if args.mode == 'inventory':
        print('\nSTAGE MAPPING SUGGESTIONS (review before using):')
        for item in inventory['stages']:
            top = ', '.join(s['script'] for s in item['suggestions'][:2]) or 'NO CANDIDATE'
            print(f" - {item['id']}: configured={item['configured_script'] or 'UNMAPPED'}; candidates={top}")
    if errors:
        print('PRE-FLIGHT ISSUES:')
        for error in errors:
            print(' -', error)
    return 0 if args.mode != 'full' or fresh else 2


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(1)
