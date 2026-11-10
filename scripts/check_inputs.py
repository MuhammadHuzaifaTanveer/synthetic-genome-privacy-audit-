#!/usr/bin/env python3
"""Validate existence and sample metadata consistency; not a full scientific validation."""
from pathlib import Path
import csv
import argparse

p = argparse.ArgumentParser()
p.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
a = p.parse_args()
root = a.project
required = [root / 'data/raw/1000g_phase3_samples.panel', root / 'data/processed/chr22_pilot_sample_metadata.csv']
for chrom in (20,21,22):
    required.append(root / f'data/raw/1000g_phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz')
    required.append(root / f'data/processed/chr{chrom}_expanded_genotypes.npz')
missing = [str(p) for p in required if not p.is_file()]
if missing:
    raise SystemExit('MISSING INPUTS:\n' + '\n'.join(missing))
with required[0].open(encoding='utf-8') as f:
    panel = [x.split() for x in f if x.strip()]
with required[1].open(encoding='utf-8', newline='') as f:
    meta = list(csv.DictReader(f))
fields=('sample','pop','super_pop','gender')
if panel[0][:4] != list(fields) or len(panel)-1 != len(meta) or len(meta) != 2504 or not all(all(meta[i][k] == row[j] for j,k in enumerate(fields)) for i,row in enumerate(panel[1:])):
    raise SystemExit('Metadata mismatch (row order or fields)')
print('Input files present and 2504 metadata records match panel order.')
