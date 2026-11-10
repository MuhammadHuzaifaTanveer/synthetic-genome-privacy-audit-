import shutil
import os

from pathlib import Path
import json
import subprocess
import time

import numpy as np
import pandas as pd
from cyvcf2 import VCF

ROOT = Path(os.environ.get("PAPER1_ROOT", Path(__file__).resolve().parents[1])).resolve()
PILOT = ROOT / "results/step31_genomator"
OUT = ROOT / "results/step33_independent_validation"
OUT.mkdir(parents=True, exist_ok=True)

META = ROOT / "data/processed/chr22_pilot_sample_metadata.csv"
SOURCE_DIR = ROOT / "data/raw/1000g_phase3"
GENOMATOR = os.environ.get("GENOMATOR_BIN", shutil.which("genomator") or "/usr/local/bin/genomator")

SEEDS = [3301, 3302, 3303, 3304, 3305]
TRAIN_PER_GROUP = 20
TEST_PER_GROUP = 10
SHADOW_PER_GROUP = 20
SYNTHETIC_N = 50
TIMEOUT_SECONDS = 180

def read_vcf(path):
    vcf = VCF(str(path))
    samples = list(vcf.samples)
    keys, rows = [], []

    for v in vcf:
        keys.append((v.CHROM, int(v.POS), v.REF, v.ALT[0]))
        rows.append([
            g[0] + g[1]
            if g[0] in (0, 1) and g[1] in (0, 1)
            else -1
            for g in v.genotypes
        ])

    vcf.close()

    if not rows:
        raise RuntimeError(f"No variants in {path}")

    return samples, keys, np.asarray(rows, dtype=np.int8).T


def calculate_auc(labels, scores):
    labels = np.asarray(labels)
    scores = np.asarray(scores)

    positive = labels == 1
    n1 = int(positive.sum())
    n0 = len(labels) - n1

    if n1 == 0 or n0 == 0:
        return float("nan")

    ranks = pd.Series(scores).rank(method="average").to_numpy()

    return float(
        (ranks[positive].sum() - n1 * (n1 + 1) / 2)
        / (n1 * n0)
    )


def membership_scores(train, test, synthetic, columns=None):
    if columns is not None:
        train = train[:, columns]
        test = test[:, columns]
        synthetic = synthetic[:, columns]

    real = np.vstack([train, test])

    distances = np.min(
        np.count_nonzero(
            real[:, None, :] != synthetic[None, :, :],
            axis=2
        ),
        axis=1
    )

    return -distances


def evaluate(train, test, synthetic, train_groups, test_groups):
    labels = np.r_[
        np.ones(len(train), dtype=int),
        np.zeros(len(test), dtype=int)
    ]

    groups = np.asarray(train_groups + test_groups)
    scores = membership_scores(train, test, synthetic)

    result = {
        "overall_auc": calculate_auc(labels, scores)
    }

    for group in sorted(set(groups)):
        mask = groups == group
        result[f"auc_{group}"] = calculate_auc(
            labels[mask], scores[mask]
        )

    # Minor allele frequency estimated from training samples only.
    af = train.mean(axis=0) / 2
    maf = np.minimum(af, 1 - af)

    rare = (maf >= 0.005) & (maf < 0.05)
    common = maf >= 0.05

    result["rare_snp_count"] = int(rare.sum())
    result["common_snp_count"] = int(common.sum())

    for name, mask in [("rare", rare), ("common", common)]:
        # Avoid interpreting tiny SNP subsets as reliable attacks.
        if mask.sum() >= 10:
            subset_scores = membership_scores(
                train, test, synthetic, mask
            )
            result[f"{name}_auc"] = calculate_auc(
                labels, subset_scores
            )
        else:
            result[f"{name}_auc"] = float("nan")

    return result


def minimum_cross_distance(left, right):
    distances = np.count_nonzero(
        left[:, None, :] != right[None, :, :],
        axis=2
    )
    return distances.min(axis=1)


def write_vcf(path, ids, sample_indices, keys, genotype_rows):
    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=22>\n")
        f.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,'
            'Description="Genotype">\n'
        )
        f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
            "\tINFO\tFORMAT\t" + "\t".join(ids) + "\n"
        )

        for key, row in zip(keys, genotype_rows):
            chrom, pos, ref, alt = key
            genotypes = [row[i] for i in sample_indices]

            f.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}"
                "\t.\tPASS\t.\tGT\t"
                + "\t".join(genotypes) + "\n"
            )


def generate(input_path, output_path, log_prefix):
    command = [
        GENOMATOR,
        str(input_path),
        str(output_path),
        str(SYNTHETIC_N),
        "1",
        "1"
    ]

    try:
        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SECONDS
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Genomator timed out: {input_path}"
        ) from exc

    Path(str(log_prefix) + "_stdout.log").write_text(
        process.stdout
    )
    Path(str(log_prefix) + "_stderr.log").write_text(
        process.stderr
    )

    if process.returncode != 0:
        raise RuntimeError(
            f"Genomator failed: {input_path}. "
            f"Inspect {log_prefix}_stderr.log"
        )


print("STEP 33: Loading chromosome 22 data...", flush=True)

_, selected_keys, _ = read_vcf(PILOT / "pilot_train.vcf")
wanted = set(selected_keys)

sources = sorted(SOURCE_DIR.glob("ALL.chr22*.vcf.gz"))
if not sources:
    raise FileNotFoundError("Chromosome 22 source VCF missing")

source = VCF(str(sources[0]))
source_samples = list(source.samples)
source_index = {s: i for i, s in enumerate(source_samples)}

found = {}

for variant in source:
    if not variant.ALT or len(variant.ALT) != 1:
        continue

    key = (
        variant.CHROM,
        int(variant.POS),
        variant.REF,
        variant.ALT[0]
    )

    if key not in wanted:
        continue

    gt = variant.genotypes

    dosage = np.array([
        g[0] + g[1]
        if g[0] in (0, 1) and g[1] in (0, 1)
        else -1
        for g in gt
    ], dtype=np.int8)

    strings = [
        f"{g[0]}{'|' if g[2] else '/'}{g[1]}"
        if g[0] in (0, 1) and g[1] in (0, 1)
        else "./."
        for g in gt
    ]

    found[key] = (dosage, strings)

    if len(found) == len(wanted):
        break

source.close()

if len(found) != len(selected_keys):
    raise RuntimeError(
        f"Recovered {len(found)} of {len(selected_keys)} SNPs"
    )

genotype_rows = [found[k][1] for k in selected_keys]
matrix = np.column_stack([found[k][0] for k in selected_keys])

metadata = pd.read_csv(META)
metadata["sample"] = metadata["sample"].astype(str)

# Only include individuals with complete genotypes at these SNPs.
complete = np.all(matrix >= 0, axis=1)
eligible_samples = {
    sample for sample, valid in zip(source_samples, complete)
    if valid
}

metadata = metadata[
    metadata["sample"].isin(eligible_samples)
].drop_duplicates("sample")

groups_available = sorted(
    metadata["super_pop"].dropna().unique()
)

if len(groups_available) != 5:
    raise RuntimeError("Expected five ancestry superpopulations")

results = []
cohort_records = []

started = time.time()

for repeat, seed in enumerate(SEEDS, start=1):
    rng = np.random.default_rng(seed)

    train_ids = []
    test_ids = []
    shadow_ids = []

    for group in groups_available:
        candidates = metadata.loc[
            metadata["super_pop"] == group,
            "sample"
        ].to_numpy()

        needed = (
            TRAIN_PER_GROUP
            + TEST_PER_GROUP
            + SHADOW_PER_GROUP
        )

        if len(candidates) < needed:
            raise RuntimeError(
                f"Insufficient eligible individuals in {group}"
            )

        selected = rng.choice(
            candidates, size=needed, replace=False
        ).tolist()

        train_ids.extend(selected[:TRAIN_PER_GROUP])
        test_ids.extend(
            selected[
                TRAIN_PER_GROUP:
                TRAIN_PER_GROUP + TEST_PER_GROUP
            ]
        )
        shadow_ids.extend(selected[-SHADOW_PER_GROUP:])

    if len(set(train_ids + test_ids + shadow_ids)) != (
        len(train_ids) + len(test_ids) + len(shadow_ids)
    ):
        raise RuntimeError("Overlapping cohorts within repeat")

    train_idx = [source_index[s] for s in train_ids]
    test_idx = [source_index[s] for s in test_ids]
    shadow_idx = [source_index[s] for s in shadow_ids]

    train = matrix[train_idx]
    test = matrix[test_idx]
    shadow = matrix[shadow_idx]

    repeat_dir = OUT / f"repeat_{repeat:02d}"
    repeat_dir.mkdir(exist_ok=True)

    actual_input = repeat_dir / "actual_train.vcf"
    shadow_input = repeat_dir / "shadow_train.vcf"
    actual_output = repeat_dir / "actual_synthetic.vcf"
    shadow_output = repeat_dir / "shadow_synthetic.vcf"

    write_vcf(
        actual_input, train_ids, train_idx,
        selected_keys, genotype_rows
    )
    write_vcf(
        shadow_input, shadow_ids, shadow_idx,
        selected_keys, genotype_rows
    )

    print(
        f"Repeat {repeat}/{len(SEEDS)}: generating actual...",
        flush=True
    )
    generate(
        actual_input, actual_output,
        repeat_dir / "actual"
    )

    print(
        f"Repeat {repeat}/{len(SEEDS)}: generating shadow...",
        flush=True
    )
    generate(
        shadow_input, shadow_output,
        repeat_dir / "shadow"
    )

    actual_samples, actual_keys, actual_synth = read_vcf(
        actual_output
    )
    shadow_samples, shadow_keys, shadow_synth = read_vcf(
        shadow_output
    )

    if actual_keys != selected_keys or shadow_keys != selected_keys:
        raise RuntimeError("Generated SNP identity mismatch")

    if (
        len(actual_samples) != SYNTHETIC_N
        or len(shadow_samples) != SYNTHETIC_N
    ):
        raise RuntimeError("Unexpected synthetic sample count")

    if np.any(actual_synth < 0) or np.any(shadow_synth < 0):
        raise RuntimeError("Missing generated genotypes")

    ancestry_map = dict(
        zip(metadata["sample"], metadata["super_pop"])
    )
    train_groups = [ancestry_map[s] for s in train_ids]
    test_groups = [ancestry_map[s] for s in test_ids]

    actual_metrics = evaluate(
        train, test, actual_synth,
        train_groups, test_groups
    )
    shadow_metrics = evaluate(
        train, test, shadow_synth,
        train_groups, test_groups
    )

    # Diagnostic only: not a formal kinship estimate.
    test_to_train = minimum_cross_distance(test, train)
    shadow_to_train = minimum_cross_distance(shadow, train)

    row = {
        "repeat": repeat,
        "seed": seed,
        "snps": len(selected_keys),
        "actual_auc": actual_metrics["overall_auc"],
        "shadow_auc": shadow_metrics["overall_auc"],
        "auc_gap": (
            actual_metrics["overall_auc"]
            - shadow_metrics["overall_auc"]
        ),
        "test_min_distance_to_train": int(
            test_to_train.min()
        ),
        "shadow_min_distance_to_train": int(
            shadow_to_train.min()
        ),
        "test_near_duplicate_count": int(
            (test_to_train <= 2).sum()
        ),
        "shadow_near_duplicate_count": int(
            (shadow_to_train <= 2).sum()
        ),
        "rare_snp_count": actual_metrics["rare_snp_count"],
        "common_snp_count": actual_metrics["common_snp_count"],
    }

    for group in groups_available:
        row[f"actual_{group}_auc"] = actual_metrics[f"auc_{group}"]
        row[f"shadow_{group}_auc"] = shadow_metrics[f"auc_{group}"]

    for category in ("rare", "common"):
        row[f"actual_{category}_auc"] = actual_metrics[
            f"{category}_auc"
        ]
        row[f"shadow_{category}_auc"] = shadow_metrics[
            f"{category}_auc"
        ]

    results.append(row)

    for role, ids in [
        ("train", train_ids),
        ("heldout", test_ids),
        ("shadow_train", shadow_ids)
    ]:
        for sample in ids:
            cohort_records.append({
                "repeat": repeat,
                "seed": seed,
                "sample": sample,
                "role": role,
                "ancestry": ancestry_map[sample]
            })

    pd.DataFrame(results).to_csv(
        OUT / "repeat_results.csv", index=False
    )
    pd.DataFrame(cohort_records).to_csv(
        OUT / "cohort_manifests.csv", index=False
    )

    print(
        f"  Actual AUC: {row['actual_auc']:.4f} | "
        f"Shadow AUC: {row['shadow_auc']:.4f} | "
        f"Gap: {row['auc_gap']:+.4f}",
        flush=True
    )

df = pd.DataFrame(results)

summary = {
    "status": "completed",
    "generator": "Genomator",
    "repeats": len(df),
    "snps": len(selected_keys),
    "train_per_repeat": TRAIN_PER_GROUP * 5,
    "heldout_per_repeat": TEST_PER_GROUP * 5,
    "shadow_train_per_repeat": SHADOW_PER_GROUP * 5,
    "synthetic_per_generator": SYNTHETIC_N,
    "actual_auc_mean": float(df["actual_auc"].mean()),
    "actual_auc_sd": float(df["actual_auc"].std()),
    "shadow_auc_mean": float(df["shadow_auc"].mean()),
    "shadow_auc_sd": float(df["shadow_auc"].std()),
    "auc_gap_mean": float(df["auc_gap"].mean()),
    "auc_gap_sd": float(df["auc_gap"].std()),
    "all_repeats_positive_gap": bool(
        (df["auc_gap"] > 0).all()
    ),
    "runtime_seconds": round(time.time() - started, 2),
    "limitations": [
        "256-SNP pilot panel, not genome-wide",
        "Five independently sampled cohort splits",
        "Individuals may recur across different repeats",
        "Genomator random state not explicitly controlled",
        "Near-duplicate diagnostic is not genome-wide kinship",
        "Rare SNP strata may contain too few variants",
        "Published DOMIAS baseline not implemented",
        "No formal inference about ancestry disparities"
    ]
}

(OUT / "validation_summary.json").write_text(
    json.dumps(summary, indent=2)
)

print("\nSTEP 33 INDEPENDENT VALIDATION COMPLETE")
print(json.dumps(summary, indent=2))
