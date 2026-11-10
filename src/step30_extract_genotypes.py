import os

from pathlib import Path
import gzip
import json
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("PAPER1_ROOT", Path(__file__).resolve().parents[1])).resolve()
RAW = ROOT / "data" / "raw" / "1000g_phase3"
PROCESSED = ROOT / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

CHROMOSOMES = (20, 21, 22)

SUFFIX = (
    ".phase3_shapeit2_mvncall_integrated_v5b"
    ".20130502.genotypes.vcf.gz"
)

COMMON_MIN_MAF = 0.05
RARE_MIN_MAF = 0.005

COMMON_SPACING = 3500
RARE_SPACING = 6000

MIN_CALL_RATE = 0.98

metadata = pd.read_csv(
    PROCESSED / "chr22_pilot_sample_metadata.csv"
)

expected_samples = (
    metadata["sample"].astype(str).to_numpy()
)

if len(expected_samples) != 2504:
    raise ValueError(
        "Expected 2504 original individuals"
    )

if len(set(expected_samples)) != len(expected_samples):
    raise ValueError("Duplicate sample IDs in metadata")


def decode_genotype(field):
    """
    Decode diploid biallelic genotype dosage.

    Valid examples:
    0|0 -> 0
    0|1 -> 1
    1|1 -> 2

    Missing or unsupported calls -> -1
    """
    gt = field.split(b":", 1)[0]

    if (
        len(gt) != 3
        or gt[1] not in (ord("|"), ord("/"))
        or gt[0] not in (48, 49)
        or gt[2] not in (48, 49)
    ):
        return -1

    return (gt[0] - 48) + (gt[2] - 48)


def get_info_af(info):
    """
    Read global alternate-allele frequency
    from the VCF INFO field.
    """
    marker = b";AF="

    if info.startswith(b"AF="):
        value = info[3:].split(b";", 1)[0]
    elif marker in info:
        value = info.split(marker, 1)[1].split(
            b";", 1
        )[0]
    else:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def process_chromosome(chromosome):
    filename = f"ALL.chr{chromosome}{SUFFIX}"
    source = RAW / filename

    if not source.exists():
        raise FileNotFoundError(source)

    positions = []
    genotypes = []
    frequencies = []
    variant_classes = []

    last_common = -COMMON_SPACING
    last_rare = -RARE_SPACING

    counters = {
        "variant_lines": 0,
        "candidate_lines": 0,
        "selected_common": 0,
        "selected_rare": 0,
        "missing_af": 0,
        "low_call_rate": 0
    }

    sample_order = None
    previous_position = -1

    print(
        f"\nProcessing chromosome {chromosome}",
        flush=True
    )

    # Iterating to EOF also checks the gzip CRC.
    with gzip.open(source, "rb") as handle:
        for line in handle:
            if line.startswith(b"##"):
                continue

            if line.startswith(b"#CHROM"):
                header = line.rstrip(b"\r\n").split(
                    b"\t"
                )

                vcf_samples = [
                    value.decode("utf-8")
                    for value in header[9:]
                ]

                if (
                    len(vcf_samples) != len(expected_samples)
                    or set(vcf_samples) != set(expected_samples)
                ):
                    raise ValueError(
                        f"Sample mismatch on chr{chromosome}"
                    )

                lookup = {
                    sample: index
                    for index, sample in enumerate(vcf_samples)
                }

                sample_order = [
                    lookup[sample]
                    for sample in expected_samples
                ]

                print(
                    f"Verified {len(sample_order)} samples",
                    flush=True
                )
                continue

            if line.startswith(b"#"):
                continue

            if sample_order is None:
                raise ValueError("VCF header missing")

            counters["variant_lines"] += 1

            fields = line.rstrip(b"\r\n").split(
                b"\t", 9
            )

            if len(fields) != 10:
                raise ValueError("Malformed VCF record")

            chrom, pos, _, ref, alt, _, _, info, fmt, rest = fields

            position = int(pos)

            if position < previous_position:
                raise ValueError(
                    "VCF positions are not sorted"
                )

            previous_position = position

            # Only simple biallelic A/C/G/T SNPs.
            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in (b"A", b"C", b"G", b"T")
                or alt not in (b"A", b"C", b"G", b"T")
            ):
                continue

            if fmt.split(b":", 1)[0] != b"GT":
                continue

            if (
                position < last_common + COMMON_SPACING
                and position < last_rare + RARE_SPACING
            ):
                continue

            af = get_info_af(info)

            if af is None:
                counters["missing_af"] += 1
                continue

            approximate_maf = min(af, 1.0 - af)

            if approximate_maf < RARE_MIN_MAF:
                continue

            if approximate_maf >= COMMON_MIN_MAF:
                category = "common"

                if position < last_common + COMMON_SPACING:
                    continue
            else:
                category = "low_frequency"

                if position < last_rare + RARE_SPACING:
                    continue

            counters["candidate_lines"] += 1

            sample_fields = rest.split(b"\t")

            if len(sample_fields) != len(vcf_samples):
                raise ValueError(
                    "Unexpected genotype field count"
                )

            dosage = np.fromiter(
                (
                    decode_genotype(sample_fields[index])
                    for index in sample_order
                ),
                dtype=np.int8,
                count=len(expected_samples)
            )

            valid = dosage >= 0
            call_rate = float(valid.mean())

            if call_rate < MIN_CALL_RATE:
                counters["low_call_rate"] += 1
                continue

            actual_af = float(
                dosage[valid].mean() / 2
            )

            actual_maf = min(
                actual_af,
                1.0 - actual_af
            )

            if actual_maf < RARE_MIN_MAF:
                continue

            # Reclassify using actual genotypes.
            actual_category = (
                "common"
                if actual_maf >= COMMON_MIN_MAF
                else "low_frequency"
            )

            if actual_category == "common":
                if position < last_common + COMMON_SPACING:
                    continue
                last_common = position
                counters["selected_common"] += 1
            else:
                if position < last_rare + RARE_SPACING:
                    continue
                last_rare = position
                counters["selected_rare"] += 1

            positions.append(position)
            genotypes.append(dosage)
            frequencies.append(actual_maf)
            variant_classes.append(actual_category)

            if len(positions) % 5000 == 0:
                print(
                    f"  Selected {len(positions):,} SNPs",
                    flush=True
                )

    if not genotypes:
        raise ValueError(
            f"No SNPs selected on chr{chromosome}"
        )

    # SNP-major -> sample-major matrix.
    matrix = np.stack(
        genotypes,
        axis=0
    ).T.copy()

    output = (
        PROCESSED
        / f"chr{chromosome}_expanded_genotypes.npz"
    )

    np.savez_compressed(
        output,
        genotypes=matrix,
        samples=expected_samples,
        positions=np.array(
            positions,
            dtype=np.int32
        ),
        maf=np.array(
            frequencies,
            dtype=np.float32
        ),
        variant_class=np.array(
            variant_classes,
            dtype="U16"
        ),
        chromosome=np.array(
            chromosome,
            dtype=np.int16
        )
    )

    counters["total_selected"] = len(positions)
    counters["output_file"] = output.name

    print(
        f"Chromosome {chromosome}: "
        f"{len(positions):,} SNPs selected",
        flush=True
    )

    print(
        f"Common: {counters['selected_common']:,} | "
        f"Lower-frequency: {counters['selected_rare']:,}",
        flush=True
    )

    print(
        "Compressed VCF read completed "
        "(gzip integrity checked).",
        flush=True
    )

    return counters


def main():
    print("=== STEP 30: GENOTYPE EXTRACTION ===")
    print("Assembly: GRCh37")
    print("Individuals:", len(expected_samples))
    print("Chromosomes:", CHROMOSOMES)
    print("Minimum call rate:", MIN_CALL_RATE)
    print("Minimum MAF:", RARE_MIN_MAF)

    results = {}
    total = 0

    for chromosome in CHROMOSOMES:
        result = process_chromosome(chromosome)
        results[str(chromosome)] = result
        total += result["total_selected"]

    report = {
        "assembly": "GRCh37",
        "individuals": len(expected_samples),
        "total_snps": total,
        "chromosomes": results,
        "selection": {
            "common_min_maf": COMMON_MIN_MAF,
            "rare_min_maf": RARE_MIN_MAF,
            "common_spacing_bp": COMMON_SPACING,
            "rare_spacing_bp": RARE_SPACING,
            "minimum_call_rate": MIN_CALL_RATE,
            "ld_pruned": False
        }
    }

    report_path = (
        PROCESSED / "step30_extraction_report.json"
    )

    report_path.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8"
    )

    print("\n=== EXTRACTION SUMMARY ===")

    for chromosome, result in results.items():
        print(
            f"chr{chromosome}: "
            f"{result['total_selected']:,} SNPs"
        )

    print(f"\nTotal SNPs: {total:,}")
    print(f"Individuals: {len(expected_samples):,}")
    print(f"Report: {report_path}")

    print("\nIMPORTANT:")
    print(
        "Lower-frequency means MAF 0.5%-5%; "
        "ultra-rare variants are excluded."
    )
    print(
        "Variants are position-thinned, "
        "not linkage-disequilibrium-pruned."
    )
    print(
        "No relatedness filtering has "
        "been performed yet."
    )

    print("\nSTEP 30 EXTRACTION COMPLETE")


if __name__ == "__main__":
    main()
