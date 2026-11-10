#!/usr/bin/env python3
"""Create the 2504-sample metadata from the 1000G Phase 3 panel (no genotypes)."""
from pathlib import Path
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
a = parser.parse_args()
panel = a.project / 'data/raw/1000g_phase3_samples.panel'
output = a.project / 'data/processed/chr22_pilot_sample_metadata.csv'
with panel.open(encoding='utf-8') as f:
    rows = [line.split() for line in f if line.strip()]
fields = ('sample', 'pop', 'super_pop', 'gender')
if not rows or tuple(rows[0][:4]) != fields:
    raise SystemExit('Unexpected panel columns; refusing to generate metadata')
data = rows[1:]
if len(data) != 2504 or any(len(row) < 4 for row in data) or len({row[0] for row in data}) != 2504:
    raise SystemExit('Unexpected panel sample count, missing columns, or duplicate samples')
output.parent.mkdir(parents=True, exist_ok=True)
with output.open('w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow([*fields, 'Unnamed: 4', 'Unnamed: 5'])
    writer.writerows([*row[:4], '', ''] for row in data)
print(f'Wrote {len(data)} ordered samples: {output}')
