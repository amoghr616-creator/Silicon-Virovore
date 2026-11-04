from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path.cwd()
RESULTS = ROOT / "prototype4" / "results"
OUT = RESULTS / "haddock_run_reproducibility_20261101"

CSV_PATH = OUT / "pose_contact_fingerprints.csv"

RUN_DIRS = {
    "A_pilot": RESULTS / "haddock3_chain_A_pilot_20261031",
    "A_repeat": RESULTS / "haddock3_chain_A_rep_seed20261031",
    "B_pilot": RESULTS / "haddock3_chain_B_pilot_20261031",
    "B_repeat": RESULTS / "haddock3_chain_B_rep_seed20261031",
}

STAGES = [
    "1_rigidbody",
    "2_seletop",
    "3_flexref",
    "4_emref",
    "5_clustfcc",
    "6_seletopclusts",
]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def flatten(obj, prefix=""):
    """Flatten JSON scalar values without expanding huge structures."""
    result = {}

    if isinstance(obj, dict):
        for key, value in obj.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            result.update(flatten(value, path))

    elif isinstance(obj, list):
        if len(obj) <= 30 and all(
            isinstance(x, (str, int, float, bool)) or x is None
            for x in obj
        ):
            result[prefix] = json.dumps(obj, sort_keys=True)
        else:
            result[f"{prefix}.__length__"] = str(len(obj))
            for i, value in enumerate(obj[:5]):
                result.update(flatten(value, f"{prefix}[{i}]"))

    else:
        result[prefix] = repr(obj)[:250]

    return result


def read_io(run_dir, stage):
    path = run_dir / stage / "io.json"
    if not path.exists():
        return None

    try:
        return json.loads(path.read_text())
    except Exception as exc:
        print("JSON ERROR:", path, exc)
        return None


def show_input_metadata(run_name, run_dir, stage):
    data = read_io(run_dir, stage)

    if data is None:
        return None

    inp = data.get("input", {})
    out = data.get("output", {})

    print(f"\n  {stage}")
    print("    input type:", type(inp).__name__)
    print("    output type:", type(out).__name__)

    if isinstance(inp, dict):
        print("    input keys:", sorted(inp.keys()))

    if isinstance(out, dict):
        print("    output keys:", sorted(out.keys()))

    # Print relevant flattened input values.
    values = flatten(inp)
    selected = {}

    for key, value in values.items():
        low = key.lower()
        if any(term in low for term in (
            "seed", "restraint", "ambig", "unambig",
            "active", "passive", "sampling", "pdb",
            "topol", "mol", "parameter", "cutoff",
            "nmodel", "niter", "radius", "center",
            "water", "rigidbody", "flexref",
        )):
            selected[key] = value

    for key in sorted(selected)[:50]:
        print(f"    {key}: {selected[key]}")

    return {
        "input": inp,
        "output": out,
        "path": str(run_dir / stage / "io.json"),
    }


def compare_stage_inputs(label1, dir1, label2, dir2, stage):
    a = read_io(dir1, stage)
    b = read_io(dir2, stage)

    if a is None or b is None:
        print(f"\n{stage}: missing metadata")
        return

    fa = flatten(a.get("input", {}))
    fb = flatten(b.get("input", {}))

    common = sorted(set(fa) & set(fb))
    only_a = sorted(set(fa) - set(fb))
    only_b = sorted(set(fb) - set(fa))

    changed = [
        key for key in common
        if fa[key] != fb[key]
    ]

    print(f"\n{stage}: {label1} vs {label2}")
    print("  common input fields:", len(common))
    print("  pilot-only fields:", only_a[:15])
    print("  repeat-only fields:", only_b[:15])
    print("  changed common fields:", len(changed))

    for key in changed[:35]:
        print(f"  DIFF {key}")
        print(f"       pilot:  {fa[key]}")
        print(f"       repeat: {fb[key]}")


def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError(CSV_PATH)

    with CSV_PATH.open(newline="") as f:
        rows = list(csv.DictReader(f))

    grouped = defaultdict(list)
    for row in rows:
        grouped[row["run"]].append(row)

    print("\n========== CORRECTED MODEL AUDIT ==========")

    for run_name, run_dir in RUN_DIRS.items():
        models = grouped.get(run_name, [])

        print(f"\n{run_name}")
        print("  directory exists:", run_dir.exists())
        print("  CSV models:", len(models))

        hashes = []
        hash_mismatches = []
        contact_counts = []
        distances = []

        for row in models:
            recorded = row["sha256"]
            hashes.append(recorded)
            contact_counts.append(int(row["contact_pair_count"]))

            actual_path = ROOT / row["path"]

            if not actual_path.exists():
                hash_mismatches.append(
                    f"missing: {row['path']}"
                )
            elif sha256(actual_path) != recorded:
                hash_mismatches.append(
                    f"hash mismatch: {row['path']}"
                )

            try:
                distances.append(
                    float(row["minimum_interface_distance_A"])
                )
            except (ValueError, TypeError):
                pass

        print("  unique recorded hashes:", len(set(hashes)))
        print("  hash verification errors:", hash_mismatches[:10])
        print("  contact counts:", contact_counts)
        print("  unique minimum distances:", sorted(set(distances)))

    print("\n========== STAGE INPUT METADATA ==========")

    for run_name, run_dir in RUN_DIRS.items():
        print("\n" + "=" * 70)
        print(run_name)

        for stage in STAGES:
            show_input_metadata(run_name, run_dir, stage)

    print("\n========== PILOT VS REPEAT INPUT DIFFERENCES ==========")

    for protomer in ("A", "B"):
        pilot = RUN_DIRS[f"{protomer}_pilot"]
        repeat = RUN_DIRS[f"{protomer}_repeat"]

        print("\n" + "=" * 70)
        print(f"PROTOMER {protomer}")

        for stage in STAGES:
            compare_stage_inputs(
                f"{protomer}_pilot",
                pilot,
                f"{protomer}_repeat",
                repeat,
                stage,
            )

    print("\nAudit complete.")
    print("Output source:", CSV_PATH)


if __name__ == "__main__":
    main()
