# Data provenance and reproducibility

## Data source

The input is **1000 Genomes Project Phase 3 (GRCh37), 20130502**, 2,504 sample panel, chromosomes 20–22. Source landing page: https://www.internationalgenome.org/data-portal/data-collection/phase3

The source release files are hosted at `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`:

- `integrated_call_samples_v3.20130502.ALL.panel`
- `ALL.chr{20,21,22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` and corresponding `.tbi`

The downloader in `src/step30_download_genomes.py` fetches the chromosome VCFs and indexes. The sample panel must also be fetched from the same source and stored at `data/raw/1000g_phase3_samples.panel`. No genotypes are included in this repository. Input SHA-256 fingerprints for the original author's files are in `docs/input_fingerprints.json`.

## Preparation (Windows + WSL Ubuntu, from repository root)

1. Install Python 3.12 and the dependencies; configure Genomator per the upstream instructions: https://github.com/aehrc/genomator . Review its **non-commercial** license before reuse. Use a documented Genomator commit and record it for a publishable replication.
2. Obtain the panel at the release URL above. Download the three chromosome VCFs using `python src/step30_download_genomes.py` with `PAPER1_ROOT` set to this clone.
3. Run `python scripts/prepare_sample_metadata.py` to generate the ordered 2,504-sample CSV. Two empty trailing columns are included to preserve the historical schema.
4. Run `python src/step30_extract_genotypes.py` from an environment that supports this script, then `python scripts/check_inputs.py`. Run `src/step30_validate_genotypes.py` as an additional check.
5. Put the actual Genomator CLI on WSL `PATH` or set `GENOMATOR_BIN` in the WSL environment. For the Windows orchestrator, supply `--python-wsl` pointing to the chosen WSL Python, e.g. `python reproduce/run_all.py --mode inventory --python-wsl /path/to/venv/bin/python`.
6. Only after backing up any preexisting outputs and reviewing stage outputs, use `python reproduce/run_all.py --mode full --allow-overwrite --python-wsl /path/to/venv/bin/python`. **This may overwrite local results.** The runner uses Windows `wsl.exe`; it is not a native Linux driver.

**Environment caveats:** `requirements.in` provides the direct dependencies; `environment/observed_versions.json` records versions discoverable on the preparation machine. These are not a fully verified lockfile. Install the externally licensed Genomator separately. WSL distribution and executable paths need local configuration.

**Method caveats:** Genomator may not be deterministic solely because NumPy random seeds are fixed. The runner's reported verdict must be checked against actual fresh outputs. A claim of end-to-end reproduction requires running the full pipeline in an independent environment, not just unit tests or an inventory.

## Output provenance

`public_results/` contains **previously generated aggregate results** copied from the author's local experiment outputs, not newly reproduced values. Manifest hashes allow checking the exact published files. Individual-level scores, cohort manifests, raw/derived genotypes, logs, and manuscripts are excluded.
