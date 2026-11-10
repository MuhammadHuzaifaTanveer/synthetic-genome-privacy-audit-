import os

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("PAPER1_ROOT", Path(__file__).resolve().parents[1])).resolve()
RESULTS = ROOT / "results"
OUT = RESULTS / "step35_publication"
FIGURES = OUT / "figures"
TABLES = OUT / "tables"
MANUSCRIPT = ROOT / "manuscript"

for folder in (OUT, FIGURES, TABLES, MANUSCRIPT):
    folder.mkdir(parents=True, exist_ok=True)


def load_json(path):
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


print("STEP 35: Loading verified experiment results...", flush=True)

attack_256 = load_json(
    RESULTS / "step31_genomator" / "step32_attack_results.json"
)
shadow_256 = load_json(
    RESULTS / "step32_shadow_control" / "shadow_control_results.json"
)
attack_1024 = load_json(
    RESULTS / "step31_genomator_1024" / "step32_attack_results.json"
)
shadow_1024 = load_json(
    RESULTS / "step34_shadow_control_1024" / "shadow_control_results.json"
)
repeated = pd.read_csv(
    RESULTS / "step33_independent_validation" / "repeat_results.csv"
)

# Validate that the source files refer to the expected experiments.
assert attack_256["snps"] == 256
assert attack_1024["snps"] == 1024
assert shadow_256["snps"] == 256
assert shadow_1024["snps"] == 1024
assert len(repeated) == 5

summary_rows = [
    {
        "experiment": "256 SNPs: original",
        "snps": 256,
        "splits": 1,
        "actual_auc": attack_256["overall_auc"],
        "shadow_auc": shadow_256["shadow_generator_auc"],
        "auc_gap": shadow_256["actual_minus_shadow_auc"],
        "ci_lower": attack_256["bootstrap_95_percent_ci"][0],
        "ci_upper": attack_256["bootstrap_95_percent_ci"][1],
    },
    {
        "experiment": "256 SNPs: five-split mean",
        "snps": 256,
        "splits": 5,
        "actual_auc": repeated["actual_auc"].mean(),
        "shadow_auc": repeated["shadow_auc"].mean(),
        "auc_gap": repeated["auc_gap"].mean(),
        "ci_lower": np.nan,
        "ci_upper": np.nan,
    },
    {
        "experiment": "1024 SNPs: original",
        "snps": 1024,
        "splits": 1,
        "actual_auc": attack_1024["overall_auc"],
        "shadow_auc": shadow_1024["shadow_generator_auc"],
        "auc_gap": shadow_1024["actual_minus_shadow_auc"],
        "ci_lower": attack_1024["bootstrap_95_percent_ci"][0],
        "ci_upper": attack_1024["bootstrap_95_percent_ci"][1],
    }
]

summary = pd.DataFrame(summary_rows)
summary.to_csv(TABLES / "table1_overall_results.csv", index=False)

# Table 2: ancestry-specific results at both SNP scales.
ancestry_rows = []

for snps, attack in [(256, attack_256), (1024, attack_1024)]:
    for ancestry, metrics in attack["ancestry_results"].items():
        ancestry_rows.append({
            "snps": snps,
            "ancestry": ancestry,
            "train_n": metrics["train_n"],
            "heldout_n": metrics["heldout_n"],
            "auc": metrics["auc"],
            "train_mean_distance": metrics["train_mean_distance"],
            "heldout_mean_distance": metrics["heldout_mean_distance"],
        })

ancestry = pd.DataFrame(ancestry_rows)
ancestry.to_csv(TABLES / "table2_ancestry_results.csv", index=False)

# Table 3: repeated split and rare/common variant results.
repeated.to_csv(
    TABLES / "table3_repeated_split_results.csv",
    index=False
)

# Figure 1: Actual vs shadow control.
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(summary))
width = 0.36

ax.bar(
    x - width / 2,
    summary["actual_auc"],
    width,
    label="Actual generator"
)
ax.bar(
    x + width / 2,
    summary["shadow_auc"],
    width,
    label="Shadow generator"
)
ax.axhline(
    0.5, linestyle="--", linewidth=1,
    color="gray", label="Chance"
)
ax.set_xticks(x)
ax.set_xticklabels(summary["experiment"], rotation=12)
ax.set_ylabel("Membership-inference ROC-AUC")
ax.set_ylim(0, 1)
ax.set_title("Membership Inference: Actual vs Shadow Genomator")
ax.legend()
fig.tight_layout()
fig.savefig(FIGURES / "figure1_actual_vs_shadow.png", dpi=300)
fig.savefig(FIGURES / "figure1_actual_vs_shadow.pdf")
plt.close(fig)

# Figure 2: five repeated splits.
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(
    repeated["repeat"],
    repeated["actual_auc"],
    marker="o", label="Actual generator"
)
ax.plot(
    repeated["repeat"],
    repeated["shadow_auc"],
    marker="s", label="Shadow generator"
)
ax.axhline(0.5, linestyle="--", color="gray", linewidth=1)
ax.set_xlabel("Cohort split")
ax.set_ylabel("ROC-AUC")
ax.set_xticks(repeated["repeat"])
ax.set_ylim(0.4, 0.9)
ax.set_title("256-SNP Membership Inference Across Five Splits")
ax.legend()
fig.tight_layout()
fig.savefig(FIGURES / "figure2_repeated_splits.png", dpi=300)
fig.savefig(FIGURES / "figure2_repeated_splits.pdf")
plt.close(fig)

# Figure 3: ancestry-specific results.
pivot = ancestry.pivot(
    index="ancestry", columns="snps", values="auc"
)
fig, ax = plt.subplots(figsize=(8, 5))
pivot.plot(kind="bar", ax=ax)
ax.set_xlabel("1000 Genomes superpopulation")
ax.set_ylabel("ROC-AUC")
ax.set_ylim(0, 1)
ax.set_title("Ancestry-Stratified Attack Performance")
ax.legend(title="SNP count")
fig.tight_layout()
fig.savefig(FIGURES / "figure3_ancestry_auc.png", dpi=300)
fig.savefig(FIGURES / "figure3_ancestry_auc.pdf")
plt.close(fig)

# Figure 4: rare vs common SNP attack.
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(
    repeated["repeat"],
    repeated["actual_rare_auc"],
    marker="o", label="Rare SNPs"
)
ax.plot(
    repeated["repeat"],
    repeated["actual_common_auc"],
    marker="s", label="Common SNPs"
)
ax.axhline(0.5, linestyle="--", color="gray", linewidth=1)
ax.set_xlabel("Cohort split")
ax.set_ylabel("ROC-AUC")
ax.set_xticks(repeated["repeat"])
ax.set_ylim(0.4, 0.9)
ax.set_title("Rare vs Common SNP Membership Inference")
ax.legend()
fig.tight_layout()
fig.savefig(FIGURES / "figure4_variant_frequency.png", dpi=300)
fig.savefig(FIGURES / "figure4_variant_frequency.pdf")
plt.close(fig)

actual_mean = repeated["actual_auc"].mean()
shadow_mean = repeated["shadow_auc"].mean()
gap_mean = repeated["auc_gap"].mean()

manuscript = f"""# When Ancestry Looks Like Memorization:
# Demography-Calibrated Privacy Auditing of Synthetic Human Genomes

## Status
Exploratory research manuscript draft. Not yet submission-ready.

## Abstract (working draft)
Synthetic human genomes may preserve statistical properties of
training data while exposing information about the individuals used
to generate them. We investigated membership inference against
synthetic genotypes generated by Genomator using 1000 Genomes
Phase 3 data. A nearest-neighbor genotype-distance attack was
evaluated using ancestry-balanced training and held-out cohorts
across 256-SNP and 1024-SNP panels. Independently trained
shadow generators provided negative controls.

In an initial 256-SNP experiment, attack ROC-AUC was
{attack_256["overall_auc"]:.4f}, compared with
{shadow_256["shadow_generator_auc"]:.4f} for its shadow control.
Across five separately sampled 256-SNP cohort splits, mean
actual-generator ROC-AUC was {actual_mean:.4f}, versus
{shadow_mean:.4f} for shadow generators, yielding a mean
difference of {gap_mean:.4f}. In a single 1024-SNP experiment,
ROC-AUC was {attack_1024["overall_auc"]:.4f}, compared with
{shadow_1024["shadow_generator_auc"]:.4f} for its shadow control.

These exploratory findings suggest detectable membership
information under the evaluated generator, attack, and cohort
conditions. They do not establish ancestry-specific privacy
disparities, memorization mechanisms, or genome-wide privacy risk.

## 1. Introduction
- Motivation: synthetic genomic data sharing and privacy risks.
- Why population structure complicates membership inference.
- Distinguishing demographic similarity from training membership.
- Research question: can matched demographic shadow controls
  distinguish generator-specific membership signals?

## 2. Materials and Methods

### 2.1 Data
1000 Genomes Phase 3; chromosome 22 SNP panels.
Five superpopulations: AFR, AMR, EAS, EUR, SAS.

### 2.2 Cohorts
Ancestry-balanced training and held-out groups.
Original experiment: 100 training and 50 held-out genomes.
Shadow generator: 100 separate ancestry-matched genomes.
Fifty synthetic genomes generated per Genomator run.

### 2.3 Synthetic genome generation
Genomator; preserve exact command-line parameters and software
versions in final reproducibility documentation.

### 2.4 Membership-inference attack
Nearest-synthetic-neighbor genotype mismatch distance.
Lower distance is interpreted as stronger membership evidence.
ROC-AUC is the primary discrimination metric.

### 2.5 Negative controls
Independent ancestry-balanced shadow-generator training cohorts.

### 2.6 Uncertainty
Bootstrap confidence intervals for single-split AUC.
Five separately sampled 256-SNP cohort splits.
Individuals may recur across splits, so these are not fully
independent observations.

### 2.7 Variant-frequency stratification
Rare versus common SNP subsets defined from training allele
frequencies; counts and classification rules must be reported.

## 3. Results

### 3.1 Main membership-inference findings
See table1_overall_results.csv and Figure 1.

### 3.2 Repeated cohort validation
See table3_repeated_split_results.csv and Figure 2.

### 3.3 Ancestry-stratified results
See table2_ancestry_results.csv and Figure 3.
Differences are descriptive; no statistical evidence of
between-ancestry disparities has yet been established.

### 3.4 Variant-frequency analysis
See Figure 4.
Rare and common SNP subsets differ in size, limiting direct
interpretation of their relative attack performance.

## 4. Discussion
- Actual-generator AUC consistently exceeds shadow-control AUC.
- The result is compatible with generator-specific membership
  information but does not prove a memorization mechanism.
- Population structure and genotype similarity require
  careful control.
- Results are specific to Genomator, selected chromosome 22
  panels, cohort construction, and the evaluated attack.

## 5. Limitations and Required Validation
1. Published attack baseline such as DOMIAS has not been evaluated.
2. No genome-wide kinship/relatedness exclusion has been performed.
3. The 1024-SNP analysis uses a single cohort split.
4. Genomator stochasticity has not been fully controlled.
5. Small ancestry subgroups limit inference about disparities.
6. Genomic panels are limited to chromosome 22.
7. SNP frequency strata differ in variant count.
8. Shadow controls are generator-matched but are not a formal
   demographic-calibration model.
9. The original 256-SNP and 1024-SNP panels may differ in SNP
   composition, so their AUC difference is not a causal SNP-count
   effect.

## 6. Conclusion
The evaluated nearest-neighbor attack consistently distinguishes
training members from held-out genomes better than matched
shadow-generator controls in these exploratory experiments.
Further validation is needed before claiming demographic
privacy disparities or broadly generalizable privacy leakage.

## References
To be populated with verified bibliographic records and DOIs.
"""

manuscript_path = MANUSCRIPT / "huzaifa_paper1_draft.md"
manuscript_path.write_text(manuscript, encoding="utf-8")

report = {
    "status": "publication_package_generated",
    "research_status": "exploratory_not_submission_ready",
    "overall_results": summary_rows,
    "files": {
        "tables": str(TABLES),
        "figures": str(FIGURES),
        "manuscript": str(manuscript_path),
    },
    "remaining_critical_work": [
        "Published attack baseline",
        "Formal relatedness/kinship sensitivity analysis",
        "Ancestry-calibration evaluation",
        "Verify bibliography and methods",
        "Final manuscript scientific review",
    ],
}

(OUT / "publication_summary.json").write_text(
    json.dumps(report, indent=2, default=float),
    encoding="utf-8"
)

print("\nSTEP 35 PUBLICATION PACKAGE GENERATED")
print(f"Tables: {TABLES}")
print(f"Figures: {FIGURES}")
print(f"Manuscript: {manuscript_path}")
print("\nOverall results:")
print(summary.to_string(index=False))
print("\nSTATUS: Exploratory manuscript draft; not submission-ready.")
