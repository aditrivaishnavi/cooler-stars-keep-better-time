"""Download the TESS light curves used by the notebooks from MAST.

    python get_data.py          everything: TARS light curves for both samples (8,333 files, 0.73 GB)
                                and QLP light curves for the 60 comparison stars (660 files, 0.28 GB)
    python get_data.py --tars   only the TARS files
    python get_data.py --qlp    only the QLP files

Only sectors that the analysis uses are downloaded. Files already on disk are skipped, so the script
can be stopped and started again. Files go to data/lightcurves/, which git ignores.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests
from astroquery.mast import Observations

import wander as w

DOWNLOAD_URL = "https://mast.stsci.edu/api/v0.1/Download/file"


def tars_wanted() -> dict[int, set[int]]:
    wanted = {}
    for name in ("original_450", "replication_422"):
        for row in w.read_star_list(name):
            wanted[row["tic_id"]] = row["sectors_allowed"]
    return wanted


def qlp_wanted() -> dict[int, set[int]]:
    with (w.DATA / "qlp_check_stars.csv").open() as stream:
        return {int(r["tic_id"]): {int(s) for s in r["overlap_sectors"].split(",")} for r in csv.DictReader(stream)}


def find_files(provenance: str, suffix: str, wanted: dict[int, set[int]]) -> list[tuple[str, str]]:
    """(file name, MAST data URI) for every wanted sector of every wanted star."""
    found = set()
    ids = sorted(wanted)
    for start in range(0, len(ids), 100):
        chunk = ids[start:start + 100]
        observations = Observations.query_criteria(provenance_name=provenance, target_name=[str(t) for t in chunk])
        if len(observations):
            for row in Observations.get_product_list(observations):
                name = str(row["productFilename"])
                if not name.endswith(suffix):
                    continue
                tic = int(name.split("-")[1].split("_")[0])
                if tic in wanted and w.sector_number(name) in wanted[tic]:
                    found.add((name, str(row["dataURI"])))
        print(f"  {provenance}: looked up {min(start + 100, len(ids))} of {len(ids)} stars", flush=True)
    return sorted(found)


def download(files: list[tuple[str, str]], folder: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)

    def fetch(item: tuple[str, str]) -> bool:
        name, uri = item
        path = folder / name
        if path.exists():
            return True
        for _ in range(3):
            try:
                response = requests.get(DOWNLOAD_URL, params={"uri": uri}, timeout=180)
                response.raise_for_status()
            except requests.RequestException:
                continue
            partial = path.with_name(name + ".part")
            partial.write_bytes(response.content)
            partial.replace(path)
            return True
        return False

    with ThreadPoolExecutor(max_workers=16) as pool:
        done = list(pool.map(fetch, files))
    print(f"  {sum(done)} of {len(files)} files in {folder.relative_to(w.ROOT)}", flush=True)
    if not all(done):
        print("  some downloads failed; run the script again to retry them")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tars", action="store_true", help="only the TARS light curves")
    parser.add_argument("--qlp", action="store_true", help="only the QLP light curves")
    args = parser.parse_args()
    both = not (args.tars or args.qlp)
    if args.tars or both:
        print("TARS light curves (both samples)")
        download(find_files("TARS", "_lc.fits", tars_wanted()), w.TARS_DIR)
    if args.qlp or both:
        print("QLP light curves (60 comparison stars)")
        download(find_files("QLP", "_llc.fits", qlp_wanted()), w.QLP_DIR)


if __name__ == "__main__":
    main()
