import os

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("PAPER1_ROOT", Path(__file__).resolve().parents[1])).resolve()
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

metadata = pd.read_csv(
    PROCESSED / "chr22_pilot_sample_metadata.csv"
)
expected_samples = metadata["sample"].astype(str).to_numpy()

summary = []
total_snps = 0
total_missing = 0
total_calls = 0

print("=== STEP 30: GENOTYPE VALIDATION ===")

for chrom in (20, 21, 22):
    path = PROCESSED / f"chr{chrom}_expanded_genotypes.npz"

    with np.load(path, allow_pickle=True) as data:
        G = data["genotypes"]
        samples = data["samples"].astype(str)
        positions = data["positions"]
        maf = data["maf"]
        classes = data["variant_class"].astype(str)

    n_samples, n_snps = G.shape

    assert n_samples == 2504, "Unexpected sample count"
    assert np.array_equal(
        samples, expected_samples
    ), "Sample order mismatch"

    assert len(positions) == n_snps
    assert len(maf) == n_snps
    assert len(classes) == n_snps

    assert np.all(np.diff(positions) > 0), (
        "Positions are not strictly increasing"
    )

    assert np.all((G >= -1) & (G <= 2)), (
        "Invalid genotype values"
    )

    valid = G >= 0
    called = valid.sum(axis=0)

    assert np.all(called > 0), "Uncalled SNP detected"

    call_rates = called / n_samples

    assert np.all(call_rates >= 0.98), (
        "Call-rate threshold violated"
    )

    calculated_af = np.where(
        valid, G, 0
    ).sum(axis=0) / (2 * called)

    calculated_maf = np.minimum(
        calculated_af, 1 - calculated_af
    )

    assert np.all(calculated_maf >= 0.005 - 1e-8), (
        "MAF threshold violated"
    )

    max_maf_difference = float(
        np.max(np.abs(calculated_maf - maf))
    )

    assert max_maf_difference < 1e-5, (
        "Stored MAF values disagree with genotypes"
    )

    expected_classes = np.where(
        calculated_maf >= 0.05,
        "common",
        "low_frequency"
    )

    assert np.array_equal(
        classes, expected_classes
    ), "Variant-class mismatch"

    common_positions = positions[
        classes == "common"
    ]
    low_positions = positions[
        classes == "low_frequency"
    ]

    assert np.all(
        np.diff(common_positions) >= 3500
    ), "Common-SNP spacing violated"

    assert np.all(
        np.diff(low_positions) >= 6000
    ), "Lower-frequency spacing violated"

    missing_calls = int((~valid).sum())
    number_calls = int(G.size)

    record = {
        "chromosome": chrom,
        "individuals": int(n_samples),
        "snps": int(n_snps),
        "common_snps": int(
            (classes == "common").sum()
        ),
        "low_frequency_snps": int(
            (classes == "low_frequency").sum()
        ),
        "missing_call_fraction": (
            missing_calls / number_calls
        ),
        "minimum_maf": float(calculated_maf.min()),
        "maximum_maf": float(calculated_maf.max()),
        "maximum_maf_difference": max_maf_difference,
        "validation": "PASS"
    }

    summary.append(record)
    total_snps += n_snps
    total_missing += missing_calls
    total_calls += number_calls

    print(
        f"chr{chrom}: PASS | "
        f"{n_samples:,} individuals | "
        f"{n_snps:,} SNPs | "
        f"missing={100 * missing_calls / number_calls:.4f}%"
    )

assert total_snps >= 30000, (
    "Full-scale SNP target not reached"
)

report = {
    "status": "PASS",
    "individuals": 2504,
    "total_snps": int(total_snps),
    "total_genotype_calls": int(total_calls),
    "overall_missing_fraction": (
        total_missing / total_calls
    ),
    "chromosomes": summary,
    "limitations": [
        "Three chromosomes, not whole-genome coverage",
        "Position-thinned, not LD-pruned",
        "Ultra-rare variants excluded",
        "Relatedness not yet controlled",
        "Genotype dosages, not phased haplotypes"
    ]
}

output = RESULTS / "step30_genotype_validation.json"
output.write_text(
    json.dumps(report, indent=2),
    encoding="utf-8"
)

print("\n=== VALIDATION SUMMARY ===")
print(f"Individuals: {report['individuals']:,}")
print(f"SNPs: {total_snps:,}")
print(f"Genotype calls: {total_calls:,}")
print(
    f"Overall missing: "
    f"{100 * report['overall_missing_fraction']:.4f}%"
)
print(f"Report: {output}")
print("\nSTEP 30 VALIDATION PASSED")
