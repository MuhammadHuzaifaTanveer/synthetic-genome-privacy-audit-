
from pathlib import Path
import csv
import json
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from cyvcf2 import VCF

ROOT = Path("/mnt/d/research/Huzaifa/Paper1")
DATA = ROOT / "data"
OUT = ROOT / "results" / "step31_genomator"
OUT.mkdir(parents=True, exist_ok=True)

META = DATA / "processed" / "chr22_pilot_sample_metadata.csv"
NPZ = DATA / "processed" / "chr22_expanded_genotypes.npz"
VCF_FILES = sorted((DATA / "raw" / "1000g_phase3").glob("ALL.chr22*.vcf.gz"))
GENOMATOR = Path("/home/huzaifa/genomator-venv/bin/genomator")

SEED = 3101
TRAIN_PER_ANCESTRY = 20
TEST_PER_ANCESTRY = 10
N_SNPS = 256
N_SYNTHETIC = 50

rng = np.random.default_rng(SEED)

if not VCF_FILES:
    raise FileNotFoundError("Original chromosome 22 VCF not found")
if not GENOMATOR.exists():
    raise FileNotFoundError(f"Genomator entry point not found: {GENOMATOR}")

source = VCF_FILES[0]
metadata = pd.read_csv(META)
metadata["sample"] = metadata["sample"].astype(str)

train_ids = []
test_ids = []

for ancestry in sorted(metadata["super_pop"].dropna().unique()):
    candidates = metadata.loc[
        metadata["super_pop"] == ancestry, "sample"
    ].drop_duplicates().to_numpy()

    needed = TRAIN_PER_ANCESTRY + TEST_PER_ANCESTRY
    if len(candidates) < needed:
        raise ValueError(f"Insufficient individuals for {ancestry}")

    chosen = rng.choice(candidates, needed, replace=False)
    train_ids.extend(chosen[:TRAIN_PER_ANCESTRY])
    test_ids.extend(chosen[TRAIN_PER_ANCESTRY:])

train_ids = list(map(str, train_ids))
test_ids = list(map(str, test_ids))

if set(train_ids) & set(test_ids):
    raise RuntimeError("Training and held-out individuals overlap")

with np.load(NPZ, allow_pickle=True) as z:
    positions = np.asarray(z["positions"], dtype=np.int64)

positions = np.unique(positions)
target_positions = set(
    map(int, positions[np.linspace(0, len(positions) - 1, N_SNPS, dtype=int)])
)

reader = VCF(str(source))
index = {name: i for i, name in enumerate(reader.samples)}

missing = (set(train_ids) | set(test_ids)) - set(index)
if missing:
    raise ValueError(f"Individuals absent from source VCF: {sorted(missing)[:5]}")

train_idx = [index[x] for x in train_ids]
test_idx = [index[x] for x in test_ids]

train_vcf = OUT / "pilot_train.vcf"
test_vcf = OUT / "pilot_heldout.vcf"
synthetic_vcf = OUT / "pilot_synthetic.vcf"

def write_header(handle, samples):
    handle.write("##fileformat=VCFv4.2\n")
    handle.write("##contig=<ID=22>\n")
    handle.write(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    )
    handle.write(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples) + "\n"
    )

def genotype_text(genotype):
    a, b = int(genotype[0]), int(genotype[1])
    return f"{a}|{b}"

selected = 0

with train_vcf.open("w") as tr, test_vcf.open("w") as te:
    write_header(tr, train_ids)
    write_header(te, test_ids)

    for variant in reader:
        if variant.POS not in target_positions:
            continue
        if len(variant.REF) != 1 or len(variant.ALT) != 1:
            continue
        if len(variant.ALT[0]) != 1:
            continue

        genotypes = variant.genotypes
        selected_indices = train_idx + test_idx

        if any(
            genotypes[i][0] not in (0, 1)
            or genotypes[i][1] not in (0, 1)
            for i in selected_indices
        ):
            continue

        prefix = (
            f"{variant.CHROM}\t{variant.POS}\t.\t"
            f"{variant.REF}\t{variant.ALT[0]}\t.\tPASS\t.\tGT\t"
        )

        tr.write(
            prefix + "\t".join(genotype_text(genotypes[i]) for i in train_idx) + "\n"
        )
        te.write(
            prefix + "\t".join(genotype_text(genotypes[i]) for i in test_idx) + "\n"
        )
        selected += 1

reader.close()

if selected < 40:
    raise RuntimeError(f"Only {selected} usable SNPs; pilot requires at least 40")

print(f"Training individuals: {len(train_ids)}", flush=True)
print(f"Held-out individuals: {len(test_ids)}", flush=True)
print(f"SNPs: {selected}", flush=True)
print("Starting Genomator...", flush=True)

started = time.time()
command = [

    str(GENOMATOR),
    str(train_vcf),
    str(synthetic_vcf),
    str(N_SYNTHETIC),
    "1",
    "1",
]

try:
    run = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=900,
        cwd=str(GENOMATOR.parent),
    )
except subprocess.TimeoutExpired as exc:
    (OUT / "status.json").write_text(
        json.dumps({"status": "timeout", "seconds": 900}, indent=2)
    )
    raise RuntimeError("Genomator exceeded the 15-minute pilot limit") from exc

(OUT / "genomator_stdout.log").write_text(run.stdout)
(OUT / "genomator_stderr.log").write_text(run.stderr)

if run.returncode != 0:
    raise RuntimeError(
        f"Genomator failed with exit code {run.returncode}. "
        "Inspect genomator_stderr.log."
    )

def load_matrix(path):
    vcf = VCF(str(path))
    rows = []
    keys = []
    for variant in vcf:
        keys.append((variant.CHROM, variant.POS, variant.REF, variant.ALT[0]))
        rows.append([
            g[0] + g[1]
            if g[0] in (0, 1) and g[1] in (0, 1)
            else -1
            for g in variant.genotypes
        ])
    samples = list(vcf.samples)
    vcf.close()
    return keys, np.asarray(rows, dtype=np.int8).T, samples

train_keys, train_matrix, _ = load_matrix(train_vcf)
test_keys, test_matrix, _ = load_matrix(test_vcf)
synth_keys, synth_matrix, synth_samples = load_matrix(synthetic_vcf)

if train_keys != test_keys or train_keys != synth_keys:
    raise RuntimeError("VCF variant order or identity differs across datasets")

if np.any(synth_matrix < 0):
    raise RuntimeError("Generated VCF contains missing or invalid genotypes")

train_af = train_matrix.mean(axis=0) / 2
synth_af = synth_matrix.mean(axis=0) / 2
allele_frequency_mae = float(np.mean(np.abs(train_af - synth_af)))

nearest_distances = [
    int(np.min(np.sum(train_matrix != sample, axis=1)))
    for sample in synth_matrix
]

results = {
    "status": "completed",
    "generator": "Genomator",
    "seed": SEED,
    "training_samples": len(train_ids),
    "heldout_samples": len(test_ids),
    "synthetic_samples": len(synth_samples),
    "snps": len(train_keys),
    "ancestry_groups": sorted(metadata["super_pop"].dropna().unique().tolist()),
    "allele_frequency_mae": allele_frequency_mae,
    "synthetic_to_training_nearest_genotype_hamming": nearest_distances,
    "runtime_seconds": round(time.time() - started, 2),
    "privacy_conclusion": "Not established by this pilot",
}

(OUT / "pilot_results.json").write_text(json.dumps(results, indent=2))
pd.DataFrame({
    "chromosome": [x[0] for x in train_keys],
    "position": [x[1] for x in train_keys],
    "training_allele_frequency": train_af,
    "synthetic_allele_frequency": synth_af,
}).to_csv(OUT / "allele_frequency_comparison.csv", index=False)

pd.DataFrame({
    "sample": train_ids + test_ids,
    "split": ["train"] * len(train_ids) + ["heldout"] * len(test_ids),
}).merge(metadata[["sample", "super_pop"]], on="sample", how="left").to_csv(
    OUT / "cohort_manifest.csv", index=False
)

print("\nGENOMATOR PILOT COMPLETE")
print(json.dumps(results, indent=2))
