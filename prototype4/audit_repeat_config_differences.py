from pathlib import Path
import re

ROOT = Path.cwd()
RUN_ROOT = (
    ROOT / "prototype4" / "results"
    / "p4_reliability_repeats_20261103"
)

RUNS = {
    "1917": RUN_ROOT / "HADDOCK3_chain_A_seed1917",
    "2917": RUN_ROOT / "HADDOCK3_chain_A_seed2917",
}

STAGES = [
    "0_topoaa",
    "1_rigidbody",
    "2_seletop",
    "3_flexref",
    "4_emref",
    "5_clustfcc",
    "6_seletopclusts",
]

# Values likely to differ only because each run has its own directory.
# Do not automatically ignore any scientific parameter.
PATH_KEYS = {
    "run_dir",
    "output_dir",
    "working_dir",
}

def parse_cfg(path):
    result = {}
    section = ""
    for raw in path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line
            continue
        if "=" not in line:
            continue

        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()

        full_key = f"{section}.{key}"
        if key.lower() == "iniseed":
            value = "<SEED>"
        elif key.lower() in PATH_KEYS:
            value = "<RUN_SPECIFIC_PATH>"

        result[full_key] = value
    return result

for stage in STAGES:
    path_a = RUNS["1917"] / stage / "params.cfg"
    path_b = RUNS["2917"] / stage / "params.cfg"

    print("\n" + "=" * 80)
    print(f"STAGE: {stage}")

    if not path_a.exists() or not path_b.exists():
        print("MISSING CONFIG")
        continue

    cfg_a = parse_cfg(path_a)
    cfg_b = parse_cfg(path_b)

    keys = sorted(set(cfg_a) | set(cfg_b))
    differences = [
        (key, cfg_a.get(key, "<MISSING>"), cfg_b.get(key, "<MISSING>"))
        for key in keys
        if cfg_a.get(key) != cfg_b.get(key)
    ]

    if not differences:
        print("No remaining parameter differences.")
    else:
        print(f"{len(differences)} differences:")
        for key, value_a, value_b in differences[:80]:
            print(f"  {key}")
            print(f"    1917: {value_a}")
            print(f"    2917: {value_b}")
        if len(differences) > 80:
            print("  ... output truncated")

