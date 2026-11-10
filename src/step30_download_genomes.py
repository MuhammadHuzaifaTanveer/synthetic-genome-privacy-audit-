
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
import shutil
import time

ROOT = Path(r"D:\research\Huzaifa\Paper1")
DEST = ROOT / "data" / "raw" / "1000g_phase3"
DEST.mkdir(parents=True, exist_ok=True)

BASE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/"
    "vol1/ftp/release/20130502/"
)

CHROMOSOMES = (20, 21, 22)
SUFFIX = (
    ".phase3_shapeit2_mvncall_integrated_v5b"
    ".20130502.genotypes.vcf.gz"
)

CHUNK = 4 * 1024 * 1024


def download(url, destination):
    temporary = destination.with_name(
        destination.name + ".part"
    )

    if destination.exists():
        print(
            f"Already downloaded: {destination.name}",
            flush=True
        )
        return

    for attempt in range(1, 4):
        try:
            existing = (
                temporary.stat().st_size
                if temporary.exists()
                else 0
            )

            headers = {}
            if existing:
                headers["Range"] = f"bytes={existing}-"

            request = Request(
                url,
                headers=headers,
                method="GET"
            )

            with urlopen(request, timeout=120) as response:
                status = response.status

                if existing and status != 206:
                    existing = 0
                    mode = "wb"
                else:
                    mode = "ab" if existing else "wb"

                total = response.headers.get(
                    "Content-Length"
                )

                expected = (
                    existing + int(total)
                    if total is not None
                    else None
                )

                downloaded = existing
                last_report = downloaded

                print(
                    f"Downloading {destination.name}",
                    flush=True
                )

                with temporary.open(mode) as output:
                    while True:
                        block = response.read(CHUNK)

                        if not block:
                            break

                        output.write(block)
                        downloaded += len(block)

                        if (
                            downloaded - last_report
                            >= 25 * 1024 * 1024
                        ):
                            print(
                                f"  {downloaded / 1048576:.1f} MB",
                                flush=True
                            )
                            last_report = downloaded

            if expected is not None:
                if downloaded != expected:
                    raise IOError(
                        "Incomplete download: "
                        f"{downloaded} / {expected}"
                    )

            temporary.replace(destination)

            print(
                f"Saved: {destination.name} "
                f"({downloaded / 1048576:.1f} MB)",
                flush=True
            )
            return

        except (HTTPError, URLError, OSError) as error:
            print(
                f"Attempt {attempt}/3 failed: {error}",
                flush=True
            )

            if attempt == 3:
                raise

            time.sleep(attempt * 3)


def main():
    print("=== STEP 30: GENOMIC DATA ACQUISITION ===")
    print("Source: 1000 Genomes Phase 3")
    print("Genome assembly: GRCh37")
    print("Chromosomes: 20, 21, 22")
    print(f"Destination: {DEST}")

    for chromosome in CHROMOSOMES:
        filename = f"ALL.chr{chromosome}{SUFFIX}"

        for name in (filename, filename + ".tbi"):
            download(
                BASE_URL + name,
                DEST / name
            )

    print("\nFiles acquired:")
    for file in sorted(DEST.glob("*vcf.gz*")):
        print(
            f"{file.name}: "
            f"{file.stat().st_size / 1048576:.1f} MB"
        )

    print("\nSTEP 30 DOWNLOAD COMPLETE")
    print(
        "Next: validate files and extract "
        "a multi-chromosome genotype matrix."
    )


if __name__ == "__main__":
    main()
