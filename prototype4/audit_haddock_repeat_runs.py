import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path("prototype4/results")
OUT = ROOT / "haddock_run_reproducibility_20261101"

RUNS = [
    "haddock3_chain_A_pilot_20261031",
    "haddock3_chain_A_rep_seed20261031",
    "haddock3_chain_B_pilot_20261031",
    "haddock3_chain_B_rep_seed20261031",
]

STAGES = [
    "1_rigidbody",
    "2_seletop",
    "3_flexref",
    "4_emref",
    "5_clustfcc",
    "6_seletopclusts",
]

INTERESTING = re.compile(
    r"seed|restraint|ambig|unambig|pdb|topology|"
    r"sampling|nmodel|niter|water|cutoff|center|"
    r"radius|rigidbody|flexref|seletop|clust|molecule",
    re.I,
)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def scalar_text(value):
    if isinstance(value, (str, int, float, bool)) or value is None:
        return repr(value)[:180]
    return None


def collect(obj, found):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if INTERESTING.search(str(key)):
                scalar = scalar_text(value)
                if scalar is not None:
                    found[str(key)].add(scalar)
                elif isinstance(value, list):
                    for item in value[:5]:
                        scalar = scalar_text(item)
                        if scalar is not None:
                            found[str(key)].add(scalar)
            collect(value, found)
    elif isinstance(obj, list):
        for value in obj:
            collect(value, found)


def inspect_io(path):
    try:
        data = json.loads(path.read_text())
    except Exception as exc:
        print("  JSON ERROR:", exc)
        return

    found = defaultdict(set)
    collect(data, found)

    print("  Top-level keys:",
          list(data.keys())[:30] if isinstance(data, dict)
          else type(data).__name__)

    for key in sorted(found, key=str.lower):
        values = sorted(found[key])
        print(f"  {key}: {values[:8]}")
        if len(values) > 8:
            print(f"    ... {len(values)} unique values")


def inspect_cluster(path):
    if not path.exists():
        print("  cluster.out.gz MISSING")
        return

    print("  cluster.out.gz first 18 lines:")
    try:
        with gzip.open(path, "rt", errors="replace") as f:
            for _, line in zip(range(18), f):
                print("   ", line.rstrip()[:240])
    except Exception as exc:
        print("  READ ERROR:", exc)


def main():
    # Existing pose fingerprint data
    csv_path = OUT / "pose_contact_fingerprints.csv"

    print("\n========== MODEL HASH / DISTANCE AUDIT ==========")

    with csv_path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    for run in RUNS:
        subset = [r for r in rows if r["run"] == run]
        hashes = [r["sha256"] for r in subset]
        distances = [
            float(r["minimum_interface_distance_A"])
            for r in subset
            if r["minimum_interface_distance_A"] not in ("", "None")
        ]

        print(f"\n{run}")
        print("  models:", len(subset))
        print("  unique model hashes:", len(set(hashes)))
        print("  minimum distances:", sorted(set(distances)))
        print("  contact counts:",
              [int(r["contact_pair_count"]) for r in subset])

    for run in RUNS:
        run_dir = ROOT / run

        print("\n" + "=" * 80)
        print("RUN:", run)

        # Inspect stage parameters and seeds.
        for stage in STAGES:
            io_path = run_dir / stage / "io.json"
            if io_path.exists():
                print(f"\nSTAGE: {stage}")
                inspect_io(io_path)

        print("\nCLUSTER OUTPUT:")
        inspect_cluster(run_dir / "5_clustfcc" / "cluster.out.gz")

    print("\nAudit complete.")


if __name__ == "__main__":
    main()
