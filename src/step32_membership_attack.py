import os

from pathlib import Path
import json

import numpy as np
import pandas as pd
from cyvcf2 import VCF

ROOT = Path(os.environ.get("PAPER1_ROOT", Path(__file__).resolve().parents[1])).resolve()
OUT = ROOT / "results" / "step31_genomator"
RNG = np.random.default_rng(3201)


def load_vcf(path):
    vcf = VCF(str(path))
    samples = list(vcf.samples)
    keys = []
    rows = []

    for variant in vcf:
        keys.append(
            (variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
        )
        rows.append([
            g[0] + g[1]
            if g[0] in (0, 1) and g[1] in (0, 1)
            else -1
            for g in variant.genotypes
        ])

    vcf.close()
    return samples, keys, np.asarray(rows, dtype=np.int8).T


train_ids, train_keys, train = load_vcf(OUT / "pilot_train.vcf")
test_ids, test_keys, test = load_vcf(OUT / "pilot_heldout.vcf")
synth_ids, synth_keys, synth = load_vcf(OUT / "pilot_synthetic.vcf")

if not (train_keys == test_keys == synth_keys):
    raise RuntimeError("Variant sets differ between VCF files")

if any(np.any(x < 0) for x in (train, test, synth)):
    raise RuntimeError("Missing or invalid genotype detected")

if set(train_ids) & set(test_ids):
    raise RuntimeError("Training/held-out sample overlap")

manifest = pd.read_csv(OUT / "cohort_manifest.csv")
ancestry_map = dict(zip(manifest["sample"], manifest["super_pop"]))


def nearest_distance(real, synthetic):
    return np.min(
        np.count_nonzero(
            real[:, None, :] != synthetic[None, :, :],
            axis=2
        ),
        axis=1
    )


def auc(labels, scores):
    labels = np.asarray(labels)
    scores = np.asarray(scores)

    positives = labels == 1
    negatives = labels == 0

    if positives.sum() == 0 or negatives.sum() == 0:
        return float("nan")

    ranks = pd.Series(scores).rank(method="average").to_numpy()
    n_pos = positives.sum()
    n_neg = negatives.sum()

    return float(
        (ranks[positives].sum() - n_pos * (n_pos + 1) / 2)
        / (n_pos * n_neg)
    )


train_dist = nearest_distance(train, synth)
test_dist = nearest_distance(test, synth)

# Lower distance indicates stronger membership evidence.
labels = np.r_[
    np.ones(len(train_ids), dtype=int),
    np.zeros(len(test_ids), dtype=int)
]
distances = np.r_[train_dist, test_dist]
scores = -distances

sample_ids = train_ids + test_ids
ancestries = np.array([ancestry_map[x] for x in sample_ids])

overall_auc = auc(labels, scores)

# Bootstrap both classes independently.
boot_aucs = []
positive_idx = np.where(labels == 1)[0]
negative_idx = np.where(labels == 0)[0]

for _ in range(2000):
    chosen = np.r_[
        RNG.choice(positive_idx, len(positive_idx), replace=True),
        RNG.choice(negative_idx, len(negative_idx), replace=True)
    ]
    boot_aucs.append(auc(labels[chosen], scores[chosen]))

ci_lower, ci_upper = np.quantile(boot_aucs, [0.025, 0.975])

group_results = {}

for group in sorted(set(ancestries)):
    mask = ancestries == group

    group_results[group] = {
        "train_n": int(np.sum(mask & (labels == 1))),
        "heldout_n": int(np.sum(mask & (labels == 0))),
        "auc": round(auc(labels[mask], scores[mask]), 4),
        "train_mean_distance": round(
            float(distances[mask & (labels == 1)].mean()), 3
        ),
        "heldout_mean_distance": round(
            float(distances[mask & (labels == 0)].mean()), 3
        ),
    }

# Permute labels within ancestry, preserving group composition.
permuted_aucs = []

for _ in range(2000):
    shuffled = labels.copy()

    for group in set(ancestries):
        indices = np.where(ancestries == group)[0]
        shuffled[indices] = RNG.permutation(shuffled[indices])

    permuted_aucs.append(auc(shuffled, scores))

permutation_p = (
    1 + sum(x >= overall_auc for x in permuted_aucs)
) / (len(permuted_aucs) + 1)

results = {
    "attack": "Nearest-neighbor genotype-distance membership inference",
    "generator": "Genomator",
    "snps": len(train_keys),
    "training_samples": len(train_ids),
    "heldout_samples": len(test_ids),
    "synthetic_samples": len(synth_ids),
    "overall_auc": round(overall_auc, 4),
    "bootstrap_95_percent_ci": [
        round(float(ci_lower), 4),
        round(float(ci_upper), 4)
    ],
    "ancestry_stratified_permutation_p": round(
        float(permutation_p), 5
    ),
    "ancestry_results": group_results,
    "limitations": [
        "Single train/heldout split",
        "Only 256 SNPs",
        "Nearest-neighbor attack, not DOMIAS",
        "No relatedness control yet",
        "Bootstrap does not capture split-to-split variation"
    ]
}

(OUT / "step32_attack_results.json").write_text(
    json.dumps(results, indent=2)
)

pd.DataFrame({
    "sample": sample_ids,
    "membership": labels,
    "ancestry": ancestries,
    "nearest_synthetic_distance": distances,
    "attack_score": scores
}).to_csv(OUT / "step32_attack_scores.csv", index=False)

print("\nSTEP 32 MEMBERSHIP-INFERENCE RESULTS")
print(json.dumps(results, indent=2))
