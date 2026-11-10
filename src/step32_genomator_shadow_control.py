import shutil
import os

from pathlib import Path
import json
import subprocess

import numpy as np
import pandas as pd
from cyvcf2 import VCF

ROOT = Path(os.environ.get("PAPER1_ROOT", Path(__file__).resolve().parents[1])).resolve()
BASE = ROOT / "results" / "step31_genomator"
OUT = ROOT / "results" / "step32_shadow_control"
OUT.mkdir(parents=True, exist_ok=True)

META = ROOT / "data/processed/chr22_pilot_sample_metadata.csv"
SOURCE_DIR = ROOT / "data/raw/1000g_phase3"
GENOMATOR = os.environ.get("GENOMATOR_BIN", shutil.which("genomator") or "/usr/local/bin/genomator")

rng = np.random.default_rng(3202)


def read_matrix(path):
    reader = VCF(str(path))
    ids = list(reader.samples)
    keys, rows = [], []

    for v in reader:
        keys.append((v.CHROM, v.POS, v.REF, v.ALT[0]))
        rows.append([
            g[0] + g[1] if g[0] in (0, 1) and g[1] in (0, 1)
            else -1
            for g in v.genotypes
        ])

    reader.close()
    return ids, keys, np.asarray(rows, dtype=np.int8).T


def auc(labels, scores):
    labels = np.asarray(labels)
    ranks = pd.Series(scores).rank(method="average").to_numpy()
    positive = labels == 1
    n1 = int(positive.sum())
    n0 = len(labels) - n1
    return float(
        (ranks[positive].sum() - n1 * (n1 + 1) / 2)
        / (n1 * n0)
    )


def attack_auc(train, heldout, synthetic):
    real = np.vstack([train, heldout])
    distances = np.min(
        np.count_nonzero(
            real[:, None, :] != synthetic[None, :, :],
            axis=2
        ),
        axis=1
    )
    labels = np.r_[
        np.ones(len(train), dtype=int),
        np.zeros(len(heldout), dtype=int)
    ]
    return auc(labels, -distances)


train_ids, keys, train = read_matrix(BASE / "pilot_train.vcf")
test_ids, test_keys, heldout = read_matrix(BASE / "pilot_heldout.vcf")
_, synth_keys, actual_synthetic = read_matrix(BASE / "pilot_synthetic.vcf")

if keys != test_keys or keys != synth_keys:
    raise RuntimeError("Variant mismatch in original experiment")

metadata = pd.read_csv(META)
metadata["sample"] = metadata["sample"].astype(str)

excluded = set(train_ids) | set(test_ids)
shadow_ids = []

for ancestry in sorted(metadata["super_pop"].dropna().unique()):
    available = metadata.loc[
        (metadata["super_pop"] == ancestry)
        & (~metadata["sample"].isin(excluded)),
        "sample"
    ].drop_duplicates().to_numpy()

    if len(available) < 20:
        raise RuntimeError(f"Insufficient shadow samples: {ancestry}")

    shadow_ids.extend(
        rng.choice(available, size=20, replace=False).tolist()
    )

if set(shadow_ids) & excluded:
    raise RuntimeError("Shadow cohort overlaps evaluation cohort")

sources = sorted(SOURCE_DIR.glob("ALL.chr22*.vcf.gz"))
if not sources:
    raise FileNotFoundError("Original chromosome 22 VCF missing")

source = VCF(str(sources[0]))
index = {s: i for i, s in enumerate(source.samples)}
missing = set(shadow_ids) - set(index)

if missing:
    raise RuntimeError(f"Shadow individuals missing: {sorted(missing)[:5]}")

indices = [index[s] for s in shadow_ids]
wanted = {key[1]: key for key in keys}
found = {}

for variant in source:
    if variant.POS not in wanted:
        continue

    identity = (
        variant.CHROM,
        variant.POS,
        variant.REF,
        variant.ALT[0]
    )

    if identity != wanted[variant.POS]:
        continue

    genotypes = [variant.genotypes[i] for i in indices]
    if any(g[0] not in (0, 1) or g[1] not in (0, 1)
           for g in genotypes):
        raise RuntimeError(f"Invalid shadow genotype: {variant.POS}")

    found[identity] = [
        f"{g[0]}{'|' if g[2] else '/'}{g[1]}"
        for g in genotypes
    ]

source.close()

if len(found) != len(keys):
    raise RuntimeError(
        f"Only recovered {len(found)} of {len(keys)} SNPs"
    )

shadow_input = OUT / "shadow_train.vcf"
shadow_output = OUT / "shadow_synthetic.vcf"

with shadow_input.open("w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##contig=<ID=22>\n")
    f.write(
        '##FORMAT=<ID=GT,Number=1,Type=String,'
        'Description="Genotype">\n'
    )
    f.write(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        "\tFORMAT\t" + "\t".join(shadow_ids) + "\n"
    )

    for chrom, pos, ref, alt in keys:
        f.write(
            f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"
            + "\t".join(found[(chrom, pos, ref, alt)]) + "\n"
        )

print("Shadow cohort:", len(shadow_ids), flush=True)
print("SNPs:", len(keys), flush=True)
print("Generating 50 shadow synthetic genomes...", flush=True)

try:
    process = subprocess.run(
        [
            GENOMATOR,
            str(shadow_input),
            str(shadow_output),
            "50", "1", "1"
        ],
        capture_output=True,
        text=True,
        timeout=900
    )
except subprocess.TimeoutExpired as exc:
    raise RuntimeError("Shadow generation timed out") from exc

(OUT / "generator_stdout.log").write_text(process.stdout)
(OUT / "generator_stderr.log").write_text(process.stderr)

if process.returncode != 0:
    raise RuntimeError(
        "Shadow Genomator failed; inspect generator_stderr.log"
    )

_, shadow_keys, shadow_synthetic = read_matrix(shadow_output)

if shadow_keys != keys:
    raise RuntimeError("Shadow variant mismatch")

if np.any(shadow_synthetic < 0):
    raise RuntimeError("Invalid shadow synthetic genotypes")

actual_auc = attack_auc(train, heldout, actual_synthetic)
shadow_auc = attack_auc(train, heldout, shadow_synthetic)

results = {
    "status": "completed",
    "snps": len(keys),
    "shadow_training_individuals": len(shadow_ids),
    "shadow_synthetic_individuals": len(shadow_synthetic),
    "evaluation_training_individuals": len(train_ids),
    "evaluation_heldout_individuals": len(test_ids),
    "actual_generator_auc": round(actual_auc, 4),
    "shadow_generator_auc": round(shadow_auc, 4),
    "actual_minus_shadow_auc": round(actual_auc - shadow_auc, 4),
    "interpretation_limitations": [
        "Single independent shadow generator",
        "No relatedness exclusion",
        "No split-to-split uncertainty",
        "Not a DOMIAS implementation"
    ]
}

(OUT / "shadow_control_results.json").write_text(
    json.dumps(results, indent=2)
)

pd.DataFrame({
    "sample": shadow_ids,
    "role": "shadow_generator_training"
}).to_csv(OUT / "shadow_cohort_manifest.csv", index=False)

print("\nSTEP 32 SHADOW CONTROL RESULTS")
print(json.dumps(results, indent=2))
