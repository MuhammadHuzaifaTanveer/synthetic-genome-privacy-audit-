
from pathlib import Path
import json

import numpy as np
import pandas as pd
from cyvcf2 import VCF

ROOT = Path("/mnt/d/research/Huzaifa/Paper1")
RESULTS = ROOT / "results"
OUT = RESULTS / "step36_final_validation"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(3601)
BOOTSTRAPS = 2000
KINSHIP_THRESHOLDS = (0.0442, 0.0884)


def read_vcf(path):
    reader = VCF(str(path))
    samples = list(reader.samples)
    keys = []
    rows = []

    for variant in reader:
        keys.append((
            variant.CHROM,
            int(variant.POS),
            variant.REF,
            variant.ALT[0]
        ))
        rows.append([
            g[0] + g[1]
            if g[0] in (0, 1) and g[1] in (0, 1)
            else -1
            for g in variant.genotypes
        ])

    reader.close()

    if not rows:
        raise RuntimeError(f"Empty VCF: {path}")

    return samples, keys, np.asarray(rows, dtype=np.int8).T


def auc(labels, scores):
    labels = np.asarray(labels)
    scores = np.asarray(scores)

    pos = labels == 1
    n1 = int(pos.sum())
    n0 = len(labels) - n1

    if not n1 or not n0:
        return float("nan")

    ranks = pd.Series(scores).rank(
        method="average"
    ).to_numpy()

    return float(
        (ranks[pos].sum() - n1 * (n1 + 1) / 2)
        / (n1 * n0)
    )


def attack_scores(real, synthetic):
    if np.any(real < 0) or np.any(synthetic < 0):
        raise RuntimeError("Invalid genotype encountered")

    return -np.min(
        np.count_nonzero(
            real[:, None, :] != synthetic[None, :, :],
            axis=2
        ),
        axis=1
    ).astype(float)


def paired_bootstrap(labels, groups, actual, shadow):
    strata = [
        np.where((labels == label) & (groups == group))[0]
        for group in sorted(set(groups))
        for label in (0, 1)
    ]

    gaps = []
    actual_aucs = []
    shadow_aucs = []

    for _ in range(BOOTSTRAPS):
        indices = np.concatenate([
            RNG.choice(s, size=len(s), replace=True)
            for s in strata if len(s)
        ])

        a = auc(labels[indices], actual[indices])
        b = auc(labels[indices], shadow[indices])

        actual_aucs.append(a)
        shadow_aucs.append(b)
        gaps.append(a - b)

    return {
        "actual_auc_ci": np.quantile(
            actual_aucs, [0.025, 0.975]
        ).tolist(),
        "shadow_auc_ci": np.quantile(
            shadow_aucs, [0.025, 0.975]
        ).tolist(),
        "gap_ci": np.quantile(
            gaps, [0.025, 0.975]
        ).tolist()
    }


def load_experiment(name, base, shadow_dir):
    train_ids, keys, train = read_vcf(
        base / "pilot_train.vcf"
    )
    test_ids, test_keys, test = read_vcf(
        base / "pilot_heldout.vcf"
    )
    _, actual_keys, actual_synth = read_vcf(
        base / "pilot_synthetic.vcf"
    )
    _, shadow_keys, shadow_synth = read_vcf(
        shadow_dir / "shadow_synthetic.vcf"
    )

    if not (keys == test_keys == actual_keys == shadow_keys):
        raise RuntimeError(f"Variant mismatch: {name}")

    if set(train_ids) & set(test_ids):
        raise RuntimeError(f"Overlapping real cohorts: {name}")

    manifest = pd.read_csv(base / "cohort_manifest.csv")
    ancestry_map = dict(zip(
        manifest["sample"].astype(str),
        manifest["super_pop"].astype(str)
    ))

    sample_ids = train_ids + test_ids
    groups = np.array([ancestry_map[s] for s in sample_ids])
    labels = np.r_[
        np.ones(len(train_ids), dtype=int),
        np.zeros(len(test_ids), dtype=int)
    ]

    real = np.vstack([train, test])
    actual_scores = attack_scores(real, actual_synth)
    shadow_scores = attack_scores(real, shadow_synth)

    # The shadow generator is a matched negative control.
    # Subtracting its score is an exploratory calibration.
    calibrated = actual_scores - shadow_scores

    result = {
        "experiment": name,
        "snps": len(keys),
        "training_n": len(train_ids),
        "heldout_n": len(test_ids),
        "actual_auc": auc(labels, actual_scores),
        "shadow_auc": auc(labels, shadow_scores),
        "calibrated_auc": auc(labels, calibrated),
    }

    result["auc_gap"] = (
        result["actual_auc"] - result["shadow_auc"]
    )

    result["bootstrap"] = paired_bootstrap(
        labels, groups, actual_scores, shadow_scores
    )

    ancestry_rows = []

    for group in sorted(set(groups)):
        mask = groups == group

        ancestry_rows.append({
            "experiment": name,
            "ancestry": group,
            "training_n": int(
                np.sum(mask & (labels == 1))
            ),
            "heldout_n": int(
                np.sum(mask & (labels == 0))
            ),
            "actual_auc": auc(
                labels[mask], actual_scores[mask]
            ),
            "shadow_auc": auc(
                labels[mask], shadow_scores[mask]
            ),
            "calibrated_auc": auc(
                labels[mask], calibrated[mask]
            )
        })

    scores_df = pd.DataFrame({
        "experiment": name,
        "sample": sample_ids,
        "ancestry": groups,
        "membership": labels,
        "actual_score": actual_scores,
        "shadow_score": shadow_scores,
        "calibrated_score": calibrated
    })

    return result, ancestry_rows, scores_df


print("STEP 36: Loading experiments...", flush=True)

experiments = [
    (
        "256_SNP",
        RESULTS / "step31_genomator",
        RESULTS / "step32_shadow_control"
    ),
    (
        "1024_SNP",
        RESULTS / "step31_genomator_1024",
        RESULTS / "step34_shadow_control_1024"
    )
]

summaries = []
ancestry_rows = []
score_frames = []

for name, base, shadow in experiments:
    result, group_rows, scores = load_experiment(
        name, base, shadow
    )

    summaries.append(result)
    ancestry_rows.extend(group_rows)
    score_frames.append(scores)

    print(
        f"{name}: actual={result['actual_auc']:.4f}, "
        f"shadow={result['shadow_auc']:.4f}, "
        f"calibrated={result['calibrated_auc']:.4f}",
        flush=True
    )

scores_df = pd.concat(score_frames, ignore_index=True)
scores_df.to_csv(
    OUT / "membership_scores.csv", index=False
)

pd.DataFrame(ancestry_rows).to_csv(
    OUT / "ancestry_calibration.csv", index=False
)


# ---------------------------------------------------------
# Relatedness sensitivity using processed chromosome data
# ---------------------------------------------------------

print(
    "STEP 36: Computing multi-chromosome kinship estimates...",
    flush=True
)

# Evaluate the larger-panel cohort.
target = scores_df[
    scores_df["experiment"] == "1024_SNP"
].copy()

train_meta = target[target["membership"] == 1]
test_meta = target[target["membership"] == 0]

pairs = []

for ancestry in sorted(target["ancestry"].unique()):
    train_samples = train_meta.loc[
        train_meta["ancestry"] == ancestry, "sample"
    ].tolist()

    test_samples = test_meta.loc[
        test_meta["ancestry"] == ancestry, "sample"
    ].tolist()

    for a in train_samples:
        for b in test_samples:
            pairs.append((ancestry, a, b))

pair_df = pd.DataFrame(
    pairs, columns=["ancestry", "train_sample", "test_sample"]
)

pair_df["het_both"] = 0
pair_df["ibs0"] = 0
pair_df["het_train"] = 0
pair_df["het_test"] = 0
pair_df["valid_sites"] = 0

all_ids = sorted(
    set(pair_df["train_sample"])
    | set(pair_df["test_sample"])
)

for chromosome in ("20", "21", "22"):
    path = (
        ROOT / "data/processed"
        / f"chr{chromosome}_expanded_genotypes.npz"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Relatedness analysis requires {path}"
        )

    with np.load(path, allow_pickle=True) as data:
        sample_names = [
            str(x) for x in data["samples"]
        ]

        index = {
            name: i for i, name in enumerate(sample_names)
        }

        missing = set(all_ids) - set(index)

        if missing:
            raise RuntimeError(
                f"chr{chromosome}: missing samples "
                f"{sorted(missing)[:5]}"
            )

        genotypes = data["genotypes"]
        maf = np.asarray(data["maf"])

        if genotypes.shape[0] == len(sample_names):
            selected = genotypes[
                [index[s] for s in all_ids], :
            ]
        elif genotypes.shape[1] == len(sample_names):
            selected = genotypes[
                :, [index[s] for s in all_ids]
            ].T
        else:
            raise RuntimeError(
                f"Unexpected chr{chromosome} genotype shape"
            )

        if selected.shape[1] != len(maf):
            raise RuntimeError(
                f"chr{chromosome}: MAF dimension mismatch"
            )

        # Restrict to common variants to reduce instability.
        keep = np.isfinite(maf) & (maf >= 0.05)
        selected = np.asarray(
            selected[:, keep], dtype=np.int8
        )

    local_index = {
        sample: i for i, sample in enumerate(all_ids)
    }

    train_idx = np.array([
        local_index[s] for s in pair_df["train_sample"]
    ])
    test_idx = np.array([
        local_index[s] for s in pair_df["test_sample"]
    ])

    # Process blocks to keep memory usage manageable.
    for start in range(0, selected.shape[1], 512):
        block = selected[:, start:start + 512]

        a = block[train_idx]
        b = block[test_idx]

        valid = (
            (a >= 0) & (a <= 2)
            & (b >= 0) & (b <= 2)
        )

        pair_df["het_both"] += np.sum(
            valid & (a == 1) & (b == 1), axis=1
        )
        pair_df["ibs0"] += np.sum(
            valid & (
                ((a == 0) & (b == 2))
                | ((a == 2) & (b == 0))
            ), axis=1
        )
        pair_df["het_train"] += np.sum(
            valid & (a == 1), axis=1
        )
        pair_df["het_test"] += np.sum(
            valid & (b == 1), axis=1
        )
        pair_df["valid_sites"] += valid.sum(axis=1)

    print(
        f"  chr{chromosome}: "
        f"{selected.shape[1]} common variants",
        flush=True
    )

denominator = (
    pair_df["het_train"] + pair_df["het_test"]
)

pair_df["kinship_estimate"] = np.where(
    denominator > 0,
    (
        pair_df["het_both"] - 2 * pair_df["ibs0"]
    ) / denominator,
    np.nan
)

pair_df.to_csv(
    OUT / "kinship_train_vs_heldout.csv",
    index=False
)

maximum_kinship = pair_df.groupby(
    "test_sample"
)["kinship_estimate"].max()

sensitivity = []

actual_result = next(
    r for r in summaries if r["experiment"] == "1024_SNP"
)

for threshold in KINSHIP_THRESHOLDS:
    flagged = set(
        maximum_kinship[
            maximum_kinship >= threshold
        ].index
    )

    retain = (
        (target["membership"] == 1)
        | (~target["sample"].isin(flagged))
    )

    filtered = target.loc[retain]

    n_heldout = int(
        (filtered["membership"] == 0).sum()
    )

    if n_heldout:
        a = auc(
            filtered["membership"],
            filtered["actual_score"]
        )
        b = auc(
            filtered["membership"],
            filtered["shadow_score"]
        )
    else:
        a = float("nan")
        b = float("nan")

    sensitivity.append({
        "kinship_threshold": threshold,
        "flagged_heldout": len(flagged),
        "remaining_heldout": n_heldout,
        "actual_auc": a,
        "shadow_auc": b,
        "auc_gap": a - b
    })

pd.DataFrame(sensitivity).to_csv(
    OUT / "relatedness_sensitivity.csv",
    index=False
)

report = {
    "status": "completed",
    "experiments": summaries,
    "relatedness_sensitivity": sensitivity,
    "relatedness_method": (
        "KING-style pairwise kinship estimator using "
        "common SNPs from chromosomes 20, 21, and 22; "
        "within-ancestry train-versus-heldout pairs"
    ),
    "interpretation": (
        "Exploratory privacy audit; ancestry-calibrated "
        "score and relatedness-filtered results are "
        "sensitivity analyses, not independent validation."
    ),
    "remaining_limitations": [
        "DOMIAS or another independently implemented "
        "published membership-inference baseline is missing",
        "Relatedness screening covers only chromosomes "
        "20, 21, and 22 and is not genome-wide",
        "No LD pruning or pedigree-confirmed kinship",
        "Only one 1024-SNP cohort split",
        "No formal statistical test establishes "
        "differences in privacy exposure by ancestry",
        "Calibration uses the same evaluation cohort "
        "and is not externally validated",
        "Synthetic-genome results are specific to "
        "the evaluated Genomator settings"
    ]
}

(OUT / "final_validation.json").write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print("\nSTEP 36 FINAL VALIDATION COMPLETE")
for row in summaries:
    print(
        f"{row['experiment']}: "
        f"actual AUC={row['actual_auc']:.4f}, "
        f"shadow AUC={row['shadow_auc']:.4f}, "
        f"gap={row['auc_gap']:+.4f}, "
        f"calibrated AUC={row['calibrated_auc']:.4f}"
    )
    print(
        "  Paired bootstrap gap 95% CI:",
        [
            round(v, 4)
            for v in row["bootstrap"]["gap_ci"]
        ]
    )

print("\nRelatedness sensitivity:")
print(pd.DataFrame(sensitivity).to_string(index=False))
print("\nResults:", OUT)
print("STATUS: Final sensitivity analysis complete.")
print("Publication readiness still requires scientific review.")
