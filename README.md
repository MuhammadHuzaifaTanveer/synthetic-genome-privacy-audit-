# Demography-Calibrated Privacy Auditing of Synthetic Human Genomes

Research-code snapshot for Huzaifa's Paper 1. This repository intentionally excludes the unpublished manuscript, individual-level genotype data, generated VCFs, local Python environments, experiment logs, defense bundles, and private files.

## Reproducibility status

**Not independently verified as fully reproducible.** Nine stage scripts are mapped in `reproduce/pipeline.json`, and the original runner's three unit tests passed on the researcher's machine. This is *not* evidence that the complete experiment reruns successfully on a clean machine.

The original scripts contain hardcoded Windows and WSL paths (e.g. `D:\research\Huzaifa\Paper1`, `/mnt/d/research/Huzaifa/Paper1`, and a Genomator executable under `/home/huzaifa/`). Adapt these paths locally before attempting a rerun. Some stages overwrite outputs. The original source-to-processed-data chain and exact third-party installation/version are not fully established.

## Software and data

Python 3.12 was used in the original environment. Observed dependencies include NumPy, pandas, and cyvcf2, plus Genomator for generation. **Dependency versions are not pinned.** Genomator is not redistributed here; install and cite the upstream version separately, respecting its license.

Original analyses used public 1000 Genomes Phase 3 chromosome 20, 21 and 22 VCFs and a 2,504-sample panel. Obtain data directly from the legitimate upstream dataset; check applicable data access and usage conditions. None of the individual-level genetic datasets are included here. Preprocessing scripts are under `src/`, but additional manual setup and validation are required.

## Run basic orchestrator tests

```bash
python -m unittest discover -s reproduce/tests -v
```

## Inspect the available pipeline (no experimental rerun)

```bash
python reproduce/run_all.py --mode inventory
```

Do **not** run `--mode full` until root paths, Genomator installation, data and stage prerequisites have been fixed. See `reproduce/pipeline.json` and the scripts for exact expected outputs and parameters. The original runner can generate manuscript-related files *locally* at Step 35, but manuscripts are not committed to this repository.

## Publication and privacy policy

This repository is code-only. Do not commit unpublished manuscripts, participant-level genotypes, source credentials, local logs, or HTML defense bundles. Figures and tables should be released separately only after appropriate review. Original third-party source and license details belong with their upstream authors.


## Data, selected results and environment

- [Data acquisition and reproduction instructions](docs/DATA_AND_REPRODUCTION.md)
- [Original input SHA-256 fingerprints](docs/input_fingerprints.json)
- [Selected historical aggregate results](public_results/MANIFEST.json)
- [Observed software versions](environment/observed_versions.json)
- `requirements.in` documents direct Python dependencies; Genomator is installed separately.

**Status:** source code and selected aggregate artifacts are public; independent end-to-end reproduction is **not yet verified**. All manuscript drafts and genotype-level data are excluded. The reproducer still relies on a Windows-to-WSL workflow with local environment configuration.
