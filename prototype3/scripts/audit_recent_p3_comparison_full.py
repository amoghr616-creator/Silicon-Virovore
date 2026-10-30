from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import re
import shutil
import statistics
import sys

ROOT = Path.cwd()
RESULTS = ROOT / "prototype3" / "results"
LOGBOOK = ROOT / "prototype3" / "logbook"
SCRIPTS = ROOT / "prototype3" / "scripts"

if not RESULTS.exists():
    raise SystemExit(
        "ERROR: prototype3/results not found. Run this from "
        "~/Desktop/Silicon-Virovore."
    )

STAMP = datetime.now().strftime("%Y%m%dT%H%M%S")
OUT = LOGBOOK / f"p3_comparison_full_audit_{STAMP}"
OUT.mkdir(parents=True, exist_ok=False)

SNAPSHOT = OUT / "source_snapshot"
SNAPSHOT.mkdir()

def latest_dir(parent, pattern):
    if not parent.exists():
        return None
    options = [p for p in parent.glob(pattern) if p.is_dir()]
    return max(options, key=lambda p: p.stat().st_mtime) if options else None

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def read_text(path):
    try:
        return path.read_text(errors="replace")
    except Exception:
        return ""

def json_scalar_leaves(obj, prefix="$"):
    """Yield every scalar stored anywhere in a JSON file."""
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield from json_scalar_leaves(value, f"{prefix}.{key}")
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            yield from json_scalar_leaves(value, f"{prefix}[{i}]")
    elif obj is None or isinstance(obj, (str, int, float, bool)):
        yield prefix, obj

# ------------------------------------------------------------
# 1) Identify the paired optimizer run and the recent comparison report.
# ------------------------------------------------------------

paired = RESULTS / "p3_paired_arm_20261027T022810Z"
if not paired.is_dir():
    paired = latest_dir(RESULTS, "p3_paired_arm_*")

if paired is None:
    raise SystemExit(
        "ERROR: no p3_paired_arm_* directory found under prototype3/results."
    )

report_candidates = []
for base in [ROOT / "prototype3", LOGBOOK]:
    if not base.exists():
        continue
    for p in base.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in {
            ".md", ".txt", ".log"
        }:
            continue
        name = p.name.lower()
        if "p3_comparison_full_audit_" in str(p).lower():
            continue
        content = read_text(p)
        cl = content.lower()
        if (
            "six-area comparison" in cl
            or "baseline-versus-adaptive" in cl
            or ("p3" in name and "comparison" in name)
        ):
            report_candidates.append(p)

comparison_report = (
    max(report_candidates, key=lambda p: p.stat().st_mtime)
    if report_candidates else None
)

# Pull in the newest supporting P3 evidence directories.
# The paired optimizer histories remain the primary arm comparison.
support_roots = [
    ("native_mutation_benchmark",
     latest_dir(RESULTS, "native_score_all_mutation_pairs")),
    ("21aa_affinity_audit",
     latest_dir(RESULTS, "p3_21aa_affinity_audit")),
    ("hotspot_position_audit",
     latest_dir(RESULTS, "arise_experimental_position_audit")),
    ("pepbi_strict_matched",
     latest_dir(RESULTS, "pepbi_strict_matched")),
    ("structural_benchmark",
     latest_dir(RESULTS, "experimental_mutation_benchmark")),
]

grouped = latest_dir(
    RESULTS / "sequence_only_grouped_target_validation", "run_*"
)
if grouped:
    support_roots.append(("grouped_sequence_validation", grouped))

mutation_audit = latest_dir(
    LOGBOOK, "p3_mutation_position_audit_*"
)
if mutation_audit:
    support_roots.append(("candidate_mutation_position_audit", mutation_audit))

# Search for a local copy of the PEPBI exact-match audit or its exports.
pepbi_files = []
for p in ROOT.rglob("*"):
    if not p.is_file():
        continue
    low = p.name.lower()
    if (
        "pepbi" in low
        and ("match" in low or "audit" in low)
        and p.suffix.lower() in {
            ".csv", ".tsv", ".json", ".jsonl", ".md", ".txt", ".numbers"
        }
        and OUT not in p.parents
    ):
        pepbi_files.append(p)

# Avoid snapshotting the same folder twice.
roots = [("paired_optimizer_run", paired)]
for label, p in support_roots:
    if p:
        roots.append((label, p))

if comparison_report:
    roots.append(("six_area_comparison_report", comparison_report))

for i, p in enumerate(pepbi_files, start=1):
    roots.append((f"pepbi_audit_{i}", p))

# ------------------------------------------------------------
# 2) Inventory all sources and preserve reproducible copies.
# ------------------------------------------------------------

manifest = []
source_files = []
seen = set()

for label, root in roots:
    paths = [root] if root.is_file() else (
        list(root.rglob("*")) if root.exists() else []
    )

    for p in paths:
        if not p.is_file() or p.resolve() in seen:
            continue
        seen.add(p.resolve())
        source_files.append((label, p))

        try:
            size = p.stat().st_size
            modified = datetime.fromtimestamp(p.stat().st_mtime).isoformat(
                timespec="seconds"
            )
            digest = sha256(p)
        except Exception:
            size, modified, digest = -1, "", ""

        snap_rel = ""
        copied = False

        # Copy source files up to 25 MB into this audit package.
        # Large files remain identified in the inventory by path and hash.
        if size >= 0 and size <= 25 * 1024 * 1024:
            try:
                if root.is_file():
                    rel = Path(p.name)
                else:
                    rel = p.relative_to(root)
                dest = SNAPSHOT / label / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(p, dest)
                snap_rel = str(dest.relative_to(OUT))
                copied = True
            except Exception as exc:
                snap_rel = f"COPY_FAILED: {exc}"

        manifest.append({
            "source_group": label,
            "source_path": str(p),
            "file_name": p.name,
            "extension": p.suffix.lower(),
            "size_bytes": size,
            "modified_local": modified,
            "sha256": digest,
            "copied_to_snapshot": copied,
            "snapshot_path": snap_rel,
        })

def write_csv(path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        if not fieldnames:
            f.write("")
            return
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

write_csv(OUT / "source_file_inventory.csv", manifest)

# ------------------------------------------------------------
# 3) Extract every scalar from JSON and every sequence-bearing row
#    available in candidate history files.
# ------------------------------------------------------------

flat_rows = []
candidate_rows = []
sequence_key_names = {
    "sequence", "peptide_sequence", "candidate_sequence", "seq",
    "sequence_string", "peptide", "mutant_sequence", "child_sequence",
    "candidate_peptide", "unique_sequence"
}
aa_pattern = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

def valid_sequence(v):
    if not isinstance(v, str):
        return False
    s = re.sub(r"\s+", "", v).upper()
    return 9 <= len(s) <= 40 and bool(aa_pattern.fullmatch(s))

def arm_from(record, path):
    vals = []
    for k, v in record.items():
        if any(x in str(k).lower() for x in ("arm", "optimizer", "method",
                                               "strategy", "algorithm")):
            vals.append(str(v).lower())
    joined = " ".join(vals) + " " + str(path).lower()
    if "adaptive" in joined:
        return "adaptive"
    if "baseline" in joined or "uniform_matched" in joined:
        return "baseline"
    return "unknown"

def seed_from(record, path):
    for key in ("seed", "random_seed", "run_seed", "rng_seed"):
        value = record.get(key)
        if value is not None and str(value).strip():
            return str(value)
    m = re.search(r"seed[_-]?(\d+)", str(path), re.I)
    return m.group(1) if m else ""

def score_from(record):
    preferred = (
        "native_score", "c_score", "candidate_score",
        "fitness_score", "fitness", "score", "overall_score"
    )
    lower = {str(k).lower(): v for k, v in record.items()}
    for name in preferred:
        v = lower.get(name)
        if isinstance(v, (int, float)) and math.isfinite(float(v)):
            return float(v)
        if isinstance(v, str):
            try:
                n = float(v)
                if math.isfinite(n):
                    return n
            except ValueError:
                pass
    return None

def walk_json_candidates(obj, source_path, json_path="$", inherited=None):
    inherited = inherited or {}
    if isinstance(obj, dict):
        scalar = {
            str(k): v for k, v in obj.items()
            if v is None or isinstance(v, (str, int, float, bool))
        }
        merged = {**inherited, **scalar}
        for key, value in scalar.items():
            if key.lower() in sequence_key_names and valid_sequence(value):
                seq = re.sub(r"\s+", "", value).upper()
                candidate_rows.append({
                    "source_path": str(source_path),
                    "record_path": json_path,
                    "arm": arm_from(merged, source_path),
                    "seed": seed_from(merged, source_path),
                    "sequence": seq,
                    "length": len(seq),
                    "score": score_from(merged),
                    "record_fields_json": json.dumps(
                        merged, sort_keys=True, default=str
                    ),
                })
                break
        for k, value in obj.items():
            if isinstance(value, (dict, list)):
                walk_json_candidates(
                    value, source_path, f"{json_path}.{k}", merged
                )
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            if isinstance(value, (dict, list)):
                walk_json_candidates(
                    value, source_path, f"{json_path}[{i}]", inherited
                )

# Keep scalar extraction for all relevant JSON inputs; candidate extraction
# uses the paired optimizer run only, to avoid counting benchmark sequences
# as optimizer candidates.
for label, p in source_files:
    suffix = p.suffix.lower()

    if suffix == ".json":
        try:
            obj = json.loads(read_text(p))
            for jpath, value in json_scalar_leaves(obj):
                flat_rows.append({
                    "source_group": label,
                    "source_path": str(p),
                    "json_path": jpath,
                    "value_type": type(value).__name__,
                    "value": json.dumps(value, ensure_ascii=False)
                    if isinstance(value, (dict, list)) else value,
                })
            if paired in p.parents:
                walk_json_candidates(obj, p)
        except Exception as exc:
            flat_rows.append({
                "source_group": label,
                "source_path": str(p),
                "json_path": "$",
                "value_type": "PARSE_ERROR",
                "value": str(exc),
            })

    elif suffix in {".csv", ".tsv"} and paired in p.parents:
        try:
            delimiter = "\t" if suffix == ".tsv" else ","
            with p.open(newline="", encoding="utf-8-sig", errors="replace") as f:
                for i, row in enumerate(csv.DictReader(f), start=2):
                    for key, value in row.items():
                        if key and key.lower() in sequence_key_names and valid_sequence(value):
                            seq = re.sub(r"\s+", "", value).upper()
                            candidate_rows.append({
                                "source_path": str(p),
                                "record_path": f"line:{i}",
                                "arm": arm_from(row, p),
                                "seed": seed_from(row, p),
                                "sequence": seq,
                                "length": len(seq),
                                "score": score_from(row),
                                "record_fields_json": json.dumps(
                                    row, sort_keys=True, default=str
                                ),
                            })
                            break
        except Exception:
            pass

write_csv(
    OUT / "json_scalar_inventory.csv",
    flat_rows,
    ["source_group", "source_path", "json_path", "value_type", "value"]
)
write_csv(
    OUT / "candidate_records_extracted.csv",
    candidate_rows,
    ["source_path", "record_path", "arm", "seed", "sequence",
     "length", "score", "record_fields_json"]
)

candidate_summary = []
for arm in ("baseline", "adaptive", "unknown"):
    rows = [r for r in candidate_rows if r["arm"] == arm]
    seqs = {r["sequence"] for r in rows}
    scores = [r["score"] for r in rows if r["score"] is not None]
    candidate_summary.append({
        "arm": arm,
        "extracted_sequence_records": len(rows),
        "unique_sequences_extracted": len(seqs),
        "records_with_numeric_score": len(scores),
        "mean_extracted_score": statistics.mean(scores) if scores else "",
        "sd_extracted_score": statistics.stdev(scores) if len(scores) > 1 else "",
        "max_extracted_score": max(scores) if scores else "",
        "min_extracted_score": min(scores) if scores else "",
        "warning": (
            "Counts describe only sequence-bearing records the parser found; "
            "consult raw snapshots for omitted schemas/representations."
        ),
    })
write_csv(OUT / "candidate_extraction_summary.csv", candidate_summary)

# ------------------------------------------------------------
# 4) Encode the recent three-seed comparison and existing P3 benchmark
#    values in clear, labeled tables. These are not pooled into one metric.
# ------------------------------------------------------------

paired_rows = [
    {"seed": 2601, "arm": "baseline", "mean_native_score": 0.65109,
     "best_native_score": 0.70555, "final_generation_mean": 0.64102,
     "unique_sequences": 138},
    {"seed": 2601, "arm": "adaptive", "mean_native_score": 0.75933,
     "best_native_score": 0.84330, "final_generation_mean": 0.80516,
     "unique_sequences": 151},
    {"seed": 2602, "arm": "baseline", "mean_native_score": 0.65268,
     "best_native_score": 0.71454, "final_generation_mean": 0.64842,
     "unique_sequences": 135},
    {"seed": 2602, "arm": "adaptive", "mean_native_score": 0.75315,
     "best_native_score": 0.84553, "final_generation_mean": 0.80448,
     "unique_sequences": 151},
    {"seed": 2603, "arm": "baseline", "mean_native_score": 0.65295,
     "best_native_score": 0.70720, "final_generation_mean": 0.65695,
     "unique_sequences": 125},
    {"seed": 2603, "arm": "adaptive", "mean_native_score": 0.73549,
     "best_native_score": 0.83327, "final_generation_mean": 0.79200,
     "unique_sequences": 151},
]
write_csv(OUT / "paired_optimizer_seed_metrics.csv", paired_rows)

# This is transcribed from the newly supplied Numbers worksheet preview.
# It is intentionally labeled provisional because its match-sequence column
# did not show the matched sequence identifiers and it conflicts with the
# earlier 0-vs-0 audit.
pepbi_audit = [
    {"arm": "baseline", "unique_candidate_sequences": 275,
     "reported_exact_pepbi_matches": 0,
     "reported_match_rate": 0.0,
     "status": "PROVISIONAL: transcribed from uploaded Numbers audit; reconcile"},
    {"arm": "adaptive", "unique_candidate_sequences": 450,
     "reported_exact_pepbi_matches": 125,
     "reported_match_rate": 125 / 450,
     "status": "PROVISIONAL: transcribed from uploaded Numbers audit; reconcile"},
]
write_csv(OUT / "pepbi_exact_match_audit_provisional.csv", pepbi_audit)

mutation_rows = [
    {"position": 5, "adaptive_minus_baseline_pp": -3.33,
     "per_seed_pp": "+0.62; -5.62; -5.00"},
    {"position": 8, "adaptive_minus_baseline_pp": 3.33,
     "per_seed_pp": "+1.88; +5.62; +2.50"},
    {"position": 2, "adaptive_minus_baseline_pp": -2.92,
     "per_seed_pp": "-4.38; -3.12; -1.25"},
    {"position": 3, "adaptive_minus_baseline_pp": 2.71,
     "per_seed_pp": "-0.62; +6.88; +1.88"},
    {"position": 14, "adaptive_minus_baseline_pp": 2.71,
     "per_seed_pp": "+1.88; +6.25; +0.00"},
    {"position": 10, "adaptive_minus_baseline_pp": 2.50,
     "per_seed_pp": "+5.00; +1.88; +0.62"},
    {"position": 6, "adaptive_minus_baseline_pp": 1.46,
     "per_seed_pp": "+0.62; -0.62; +4.38"},
    {"position": 18, "adaptive_minus_baseline_pp": 1.46,
     "per_seed_pp": "+0.62; +1.25; +2.50"},
    {"position": 17, "adaptive_minus_baseline_pp": -1.25,
     "per_seed_pp": "-4.38; +0.00; +0.62"},
    {"position": 7, "adaptive_minus_baseline_pp": 0.83,
     "per_seed_pp": "+2.50; -3.75; +3.75"},
]
write_csv(OUT / "top_mutation_position_differences.csv", mutation_rows)

# A compact external-validation table. Metrics with different meanings
# stay in separate columns and figures.
external_rows = [
    {"benchmark": "Native C-score, broad matched mutation benchmark",
     "n": 693, "spearman_rho": -0.107139, "direction_accuracy": 0.424893,
     "notes": "621 unique peptides; 466 informative pairs; cluster-bootstrap "
              "rho CI [-0.2022,-0.0175]."},
    {"benchmark": "Structural scoring, all verified mapped pairs",
     "n": 18, "spearman_rho": 0.058824, "direction_accuracy": 0.50,
     "notes": "7 targets; only 18 of 693 input pairs were scored."},
    {"benchmark": "Structural scoring, same-receptor subset",
     "n": 10, "spearman_rho": 0.28485, "direction_accuracy": 0.60,
     "notes": "5 targets; rho p=.425; very small subset."},
]
write_csv(OUT / "external_validation_summary.csv", external_rows)

ridge_rows = [
    {"method": "Training-fold mean", "mae": 0.8182126,
     "rmse": 1.3655815, "spearman_rho": 0.0458706,
     "direction_accuracy": 0.6266094},
    {"method": "Sequence-only Ridge", "mae": 0.8643158,
     "rmse": 1.3956729, "spearman_rho": 0.1087403,
     "direction_accuracy": 0.5896226},
]
write_csv(OUT / "ridge_vs_mean_summary.csv", ridge_rows)

pepbi_strict_rows = [
    {"method": "Native score", "mean_group_concordance": 0.376197,
     "n_independent_groups": 11, "ci_low": 0.21097, "ci_high": 0.53149,
     "one_sided_sign_test_p": 0.8867},
    {"method": "Rosetta delta G", "mean_group_concordance": 0.392711,
     "n_independent_groups": 11, "ci_low": 0.23445, "ci_high": 0.54755,
     "one_sided_sign_test_p": 0.7256},
]
write_csv(OUT / "pepbi_strict_matched_summary.csv", pepbi_strict_rows)

# ------------------------------------------------------------
# 5) Generate graphs.
# ------------------------------------------------------------

graphs = OUT / "graphs"
graphs.mkdir()

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def savefig(name):
        plt.tight_layout()
        plt.savefig(graphs / name, dpi=180, bbox_inches="tight")
        plt.close()

    # Graph 1: mean native score by matched seed.
    seeds = [2601, 2602, 2603]
    baseline_means = [0.65109, 0.65268, 0.65295]
    adaptive_means = [0.75933, 0.75315, 0.73549]
    x = list(range(len(seeds)))
    width = 0.36
    plt.figure(figsize=(8.5, 5))
    plt.bar([i-width/2 for i in x], baseline_means, width, label="Baseline")
    plt.bar([i+width/2 for i in x], adaptive_means, width, label="Adaptive")
    plt.xticks(x, [str(s) for s in seeds])
    plt.xlabel("Matched random seed")
    plt.ylabel("Mean native internal score")
    plt.title("P3 paired optimizer run: mean native score by seed")
    plt.legend()
    savefig("01_paired_native_score_by_seed.png")

    # Graph 2: candidate uniqueness plus the provisional exact-match count.
    labels = ["Baseline", "Adaptive"]
    unique_counts = [275, 450]
    match_counts = [0, 125]
    x = list(range(2))
    plt.figure(figsize=(8.5, 5))
    plt.bar([i-width/2 for i in x], unique_counts, width,
            label="Unique candidate sequences")
    plt.bar([i+width/2 for i in x], match_counts, width,
            label="Reported PEPBI exact matches (provisional)")
    plt.xticks(x, labels)
    plt.ylabel("Sequence count")
    plt.title("Candidate diversity and new PEPBI audit counts")
    plt.legend()
    savefig("02_candidate_diversity_and_provisional_pepbi_matches.png")

    # Graph 3: mutation-position differences, descriptive only.
    positions = [r["position"] for r in mutation_rows]
    diffs = [r["adaptive_minus_baseline_pp"] for r in mutation_rows]
    plt.figure(figsize=(9, 5))
    plt.bar([str(p) for p in positions], diffs)
    plt.axhline(0, linewidth=1)
    plt.xlabel("One-based sequence position")
    plt.ylabel("Adaptive minus baseline (percentage points)")
    plt.title("Largest reported differences in mutation-position frequency")
    savefig("03_mutation_position_frequency_differences.png")

    # Graphs 4 and 5: affinity validation measures separated by metric.
    names = ["Native C-score\nmutation pairs",
             "Structural score\nall mapped",
             "Structural score\nsame receptor"]
    rhos = [r["spearman_rho"] for r in external_rows]
    dirs = [100*r["direction_accuracy"] for r in external_rows]

    plt.figure(figsize=(9, 5))
    plt.bar(names, rhos)
    plt.axhline(0, linewidth=1)
    plt.ylabel("Spearman correlation")
    plt.title("External benchmark: association with experimental affinity")
    savefig("04_external_affinity_correlations.png")

    plt.figure(figsize=(9, 5))
    plt.bar(names, dirs)
    plt.axhline(50, linewidth=1, linestyle="--")
    plt.ylabel("Direction accuracy (%)")
    plt.ylim(0, 100)
    plt.title("External benchmark: direction accuracy")
    savefig("05_external_direction_accuracy.png")

    # Graph 6: MAE makes the trade-off between Ridge and a training-fold mean clear.
    plt.figure(figsize=(8, 5))
    plt.bar([r["method"] for r in ridge_rows],
            [r["mae"] for r in ridge_rows])
    plt.ylabel("Mean absolute error (pKd units)")
    plt.title("Grouped validation: Ridge versus training-fold mean")
    savefig("06_ridge_vs_training_fold_mean_mae.png")

    # Graph 7: PEPBI strict matched group concordance.
    plt.figure(figsize=(8, 5))
    plt.bar([r["method"] for r in pepbi_strict_rows],
            [r["mean_group_concordance"] for r in pepbi_strict_rows])
    plt.ylabel("Mean concordance across independent groups")
    plt.ylim(0, 1)
    plt.title("PEPBI strict-matched benchmark: concordance")
    savefig("07_pepbi_strict_matched_concordance.png")

    graphs_ok = True
except Exception as exc:
    graphs_ok = False
    (OUT / "GRAPH_GENERATION_ERROR.txt").write_text(
        f"Graph generation failed: {exc}\n"
        "Install/repair matplotlib in the project's existing environment "
        "and rerun the audit script.\n"
    )

# ------------------------------------------------------------
# 6) Write an extensive DATA/RESULTS entry and a separate ANALYSIS entry.
# ------------------------------------------------------------

seed_table = "\n".join(
    f"| {r['seed']} | {r['mean_native_score']:.5f} | "
    f"{r['mean_native_score']:.5f} | {r['best_native_score']:.5f} | "
    f"{r['final_generation_mean']:.5f} | {r['unique_sequences']} |"
    for r in []  # deliberately populated below without mixing arms
)

paired_table = [
    "| Seed | Baseline mean | Adaptive mean | Baseline best | Adaptive best | "
    "Baseline final-gen mean | Adaptive final-gen mean | Baseline unique | "
    "Adaptive unique |",
    "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for seed in (2601, 2602, 2603):
    b = next(r for r in paired_rows
             if r["seed"] == seed and r["arm"] == "baseline")
    a = next(r for r in paired_rows
             if r["seed"] == seed and r["arm"] == "adaptive")
    paired_table.append(
        f"| {seed} | {b['mean_native_score']:.5f} | "
        f"{a['mean_native_score']:.5f} | {b['best_native_score']:.5f} | "
        f"{a['best_native_score']:.5f} | {b['final_generation_mean']:.5f} | "
        f"{a['final_generation_mean']:.5f} | {b['unique_sequences']} | "
        f"{a['unique_sequences']} |"
    )

inventory_counts = Counter(x["extension"] or "[no extension]" for x in manifest)
inventory_lines = [
    f"- Selected paired run: `{paired}`",
    f"- Six-area comparison report: `{comparison_report}`"
    if comparison_report else
    "- Six-area comparison report: no matching report file was found.",
    f"- Source files inventoried: **{len(manifest)}**",
    f"- Source files copied into the snapshot: "
    f"**{sum(1 for x in manifest if x['copied_to_snapshot'])}**",
    "- Files by extension: " + ", ".join(
        f"`{k}` {v}" for k, v in sorted(inventory_counts.items())
    ),
    f"- Candidate sequence-bearing records parsed from the paired run: "
    f"**{len(candidate_rows)}**",
    f"- JSON scalar values extracted: **{len(flat_rows)}**",
]
for label, p in support_roots:
    if p:
        inventory_lines.append(f"- Supporting evidence `{label}`: `{p}`")

data_md = f"""# Logbook Entry — P3 Comparison: Data and Results

**Date:** 27 October 2026  
**Audit generated:** {datetime.now().astimezone().isoformat(timespec="seconds")}  
**Experiment name:** P3 paired adaptive-versus-baseline comparison, candidate diversity,
mutation-position audit, and external-validation reconciliation  
**Primary run:** `{paired}`  
**Audit output:** `{OUT}`

## Design and objective

The purpose of this audit was to consolidate the most recent paired adaptive-versus-baseline
optimizer outputs with the available independent affinity-validation benchmarks. The analysis
keeps two questions separate:

1. Did adaptive guidance improve the optimizer's internal objective and alter the candidate
   sequence set compared with the uniform baseline?
2. Do the scoring and candidate-selection results correspond to experimentally measured
   binding-affinity changes?

The paired optimizer run used three matched seeds (2601, 2602, and 2603). The external
benchmarks use different datasets and sample sizes; they are reported separately rather
than pooled into a single performance number.

## Source files and reproducibility

{chr(10).join(inventory_lines)}

The `source_file_inventory.csv` records file paths, file sizes, modification times, SHA-256
hashes, and whether each file was copied into `source_snapshot/`. JSON scalar values are
indexed in `json_scalar_inventory.csv`. Sequence-bearing records that the parser could
identify from the paired run are exported to `candidate_records_extracted.csv`; the
extraction summary is `candidate_extraction_summary.csv`. The raw snapshots remain the
authoritative record if a data structure was not recognized by the parser.

## Results A — Paired optimizer outputs

| Seed | Baseline mean | Adaptive mean | Baseline best | Adaptive best | Baseline final-gen mean | Adaptive final-gen mean | Baseline unique | Adaptive unique |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(paired_table[2:])}

Across the available 480 candidate records in each arm, the previously reported inventory
contained 275 unique baseline sequences and 450 unique adaptive sequences. The corresponding
reported mean native scores were 0.65224 for baseline and 0.74933 for adaptive; the best
reported scores were 0.71454 and 0.84553, respectively. These are optimizer-objective
results, not direct experimental binding measurements.

The per-seed summary shows the same direction of the internal-score difference in all three
matched seeds. The experiment contains three seed pairs, so this is a small, exploratory
optimizer comparison; the seed summaries are not a substitute for an independent affinity
test.

## Results B — Mutation behavior and sequence diversity

The candidate mutation-position audit reported:

| Measure | Baseline | Adaptive |
|---|---:|---:|
| Parent-child records checked | 480 | 480 |
| Records with a sequence change | 480 (100%) | 450 (93.75%) |
| Mutation events counted | 480 | 510 |
| Unique parent-child sequence transitions | 275 | 473 |
| Unique candidate sequences in the candidate inventory | 275 | 450 |
| Mutation-count consistency | 480/480 | 480/480 |

The reported metadata identified 432/480 adaptive records (90%) as hotspot-guided and
48/480 (10%) as not hotspot-guided. All 480 baseline records used `uniform_matched`.
The metadata flag by itself does not prove that each particular substitution was caused by
or located within a selected hotspot.

The largest reported adaptive-minus-baseline position-frequency differences were: position
5, -3.33 percentage points; position 8, +3.33; position 2, -2.92; positions 3 and 14,
+2.71 each; and position 10, +2.50. These position differences are descriptive across only
three seeds and should not be represented as statistically significant.

## Results C — New PEPBI exact-sequence audit (provisional)

The newly supplied Numbers worksheet reports:

| Arm | Unique candidate sequences | Reported exact PEPBI matches | Reported match rate |
|---|---:|---:|---:|
| Baseline | 275 | 0 | 0% |
| Adaptive | 450 | 125 | 27.78% |

**Status: provisional and in conflict with an earlier audit that reported zero exact
matches for both arms.** The worksheet's visible `match_sequences` field did not identify
the 125 matching sequences. This audit therefore preserves the new counts as reported but
does not treat them as confirmed candidate-level experimental validation. Reconciliation
requires the actual matched sequence list, the source PEPBI files, matching rules, duplicate
handling, and the affinity/receptor-context records associated with each match.

## Results D — Native sequence-score experimental benchmark

The broad mutation benchmark contained 693 matched single-substitution pairs, 621 unique
peptides, and 466 informative pairs. The reported Spearman correlation between score change
and experimental pKd change was -0.1071, with cluster-bootstrap 95% CI
[-0.2022, -0.0175]. Direction accuracy was 42.49%, with 95% CI approximately
[36.55%, 47.78%]. A shuffle test was reported as p=0.223.

The exact-21-amino-acid audit contained 129 rows, 67 distinct peptides, and 73 measurement
contexts, but only seven contexts had enough variation for within-context analysis. The
reported mean within-context Spearman correlation was -0.487. Twenty substitution pairs
were available, only nine of them informative, and direction accuracy was 1/9 (11.11%).
The same report also contained a positive correlation summary for an informative subset
(rho=+0.678); the discrepancy between these subset summaries should be audited before
either is quoted as a definitive result.

## Results E — Mutation-position external validation

In the external PPIKB mutation dataset, the all-pair position 1–3 comparison reported a
centered effect of -0.0115, bootstrap 95% CI [-0.1690, +0.1543], and two-sided permutation
p=0.9218. This does not show a general effect across all pairs.

For the informative subset, the reported centered difference in absolute affinity change
was +0.1716, CI [+0.0148, +0.3419], with two-sided permutation p=0.05279 and one-sided
p=0.020. This narrower subset is exploratory and should not be presented as independent
confirmation of a general mutation-position effect.

## Results F — Grouped sequence-only prediction benchmark

The grouped validation used 693 pairs from 123 targets and 111 components across five
folds. The recorded leakage checks reported zero target overlap, zero exact-pair overlap,
zero peptide-sequence overlap, and one out-of-fold prediction per row.

| Metric | Training-fold mean | Sequence-only Ridge |
|---|---:|---:|
| MAE | 0.81821 | 0.86432 |
| RMSE | 1.36558 | 1.39567 |
| Spearman rho | 0.04587 | 0.10874 |
| Direction accuracy | 62.66% | 58.96% |

The paired bootstrap interval for Ridge minus mean MAE was [+0.0173, +0.0782].
The Ridge model had a higher reported rank correlation, but worse MAE and direction
accuracy. No single model dominated every metric.

## Results G — PEPBI strict-matched benchmark

The strict-matched audit reported 290 eligible rows, 69 distinct sequences, 18 strict
contexts, and 382 pairwise comparisons. Across 11 independent groups, mean concordance
was 0.3762 for the native score (95% CI [0.2110, 0.5315]) and 0.3927 for Rosetta delta G
(95% CI [0.2345, 0.5476]). Reported one-sided sign-test p-values were 0.8867 and 0.7256,
respectively. The target-reported experimental-delta-G subset had only one independent
group, which is insufficient for a reliable conclusion.

## Results H — Structural scoring benchmark

The structural benchmark began with 693 candidate mutation pairs, but only 18 mapped pairs
were successfully scored, across seven targets. Overall correlation was 0.0588 (p=0.8167)
and directional accuracy was 50%. The informative subset of 14 pairs had correlation
0.0593 (p=0.8403). The same-receptor subset had ten pairs across five targets, correlation
0.2849 (p=0.425), and directional accuracy 60%; uncertainty is substantial at this sample
size.

The pose sanity check reported a crystal-pose score of 0.96537 and a docked-pose score of
0.97751. Under the score orientation used in that audit, the score preferred the unfavorable
docked pose, so the scoring setup requires further validation before strong claims are made.

## Figures generated

{chr(10).join("- `graphs/" + p.name + "`" for p in sorted(graphs.glob("*.png"))) if graphs.exists() else "- Graph directory not available."}

## Data integrity notes

1. The 125-match adaptive PEPBI count remains provisional pending source-level reconciliation.
2. The 21-aa native-score summary contains subset-level results that need reconciliation.
3. The internal native score is not itself an independent receptor-aware binding measurement.
4. External mutation and structural benchmarks are retrospective method-validation studies;
   they do not automatically validate the actual candidates generated in the optimizer run.
5. Where the extraction parser could not recognize a file's schema, the preserved raw source
   snapshot and manifest should be used rather than inferring that the data were absent.
"""

analysis_md = f"""# Logbook Entry — P3 Comparison: Analysis and Interpretation

**Date:** 27 October 2026  
**Experiment:** Interpretation of the consolidated P3 optimizer, mutation, affinity,
and structural evidence  
**Corresponding data entry:** `logbook_DATA_RESULTS.md`  
**Audit package:** `{OUT}`

## 1. Central interpretation

This project has two distinct validation questions. The first asks whether adaptive guidance
optimizes the pipeline's internal objective better than the matched uniform baseline. The
second asks whether that objective and the chosen candidates predict real binding affinity.
The current evidence is more favorable for the first question than for the second.

## 2. Hypothesis A — Does adaptive guidance improve optimization?

The three matched seeds showed higher mean native score, best native score, and final-generation
mean score for the adaptive arm. The recorded mean native score across candidate records was
also higher for adaptive (about 0.7493) than baseline (about 0.6522). The candidate inventory
reported 450 unique adaptive sequences versus 275 baseline sequences.

This supports the narrower statement that adaptive guidance improved optimization against
the internal score in this three-seed experiment and yielded a more diverse observed sequence
set. The score difference was directionally consistent across the three seeds.

It does **not** establish that the adaptive candidates bind the target more strongly. If the
native `c_score` does not incorporate the receptor structure, optimizing that score can reward
properties of the sequence-level objective without reliably rewarding receptor-specific
binding. A higher optimizer score is evidence of optimization against that objective, not
evidence of a successful binder by itself.

The mutation audit supports that adaptive and uniform arms do not behave identically: the
adaptive arm has more unique parent-child transitions and different per-position mutation
frequencies. However, position-frequency differences vary by seed and remain descriptive.
The hotspot-guidance flag should not be interpreted as proof that every mutation was selected
because of a hotspot or lies inside the hotspot window.

## 3. Hypothesis B — Do scores or mutation positions track measured affinity?

The broad native-score benchmark provides unfavorable evidence for treating the native score
as a reliable affinity predictor: the score-change/affinity-change correlation was weakly
negative and direction accuracy was below 50%. The broad position 1–3 audit did not show a
general all-pair effect. These observations do not support a confident claim that the native
score or the tested mutation-position rule predicts experimentally favorable substitutions.

The grouped Ridge benchmark gives a mixed model comparison rather than a production-ready
predictor. Ridge ranked affinity changes somewhat better by Spearman correlation, but its
MAE, RMSE, and direction accuracy were worse than the training-fold mean. This illustrates
why model quality should not be declared from one favorable metric while ignoring the others.

The structural benchmark is underpowered: only 18 of 693 candidate mutation pairs were
successfully mapped and scored. Its overall correlation was near zero, its directional
accuracy was at chance, and the same-receptor subset of ten pairs is too small for a robust
claim. The pose sanity check also raises a scoring-validity concern because the recorded score
preferred the unfavorable pose. These results do not validate the structural score as a
reliable experimental-affinity proxy.

The PEPBI strict-matched benchmark has low-to-moderate concordance estimates with broad
confidence intervals and high reported sign-test p-values. Its sample structure is useful
for auditing methods, but the evidence does not establish dependable affinity ranking from
the native score. The comparison with Rosetta does not itself validate either score.

## 4. Critical issue — the new PEPBI exact-match count

The newly supplied worksheet reports 125 exact sequence matches for adaptive and zero for
baseline. This differs from the earlier audit that reported zero matches for both arms. It
is important not to dismiss the new result, but it is equally important not to overstate it.

Before using the count as project evidence, reconcile the audits at the sequence level:

- Export the exact sequence lists from both optimizer arms and the exact peptide sequences
  in the PEPBI source tables.
- Normalize only defensible formatting differences, such as whitespace and letter case;
  do not silently remove amino acids or alter sequences to force matches.
- Produce one row per matched candidate sequence showing the candidate sequence, arm,
  PEPBI source file, matching row or identifier, target/receptor context, experimental
  affinity field, units, and duplicate count.
- Report both the number of distinct matched sequences and the number of independent
  affinity contexts. Multiple rows for the same sequence or receptor context are not
  equivalent to independent experimental validation.
- Confirm whether matches exist because the adaptive candidates overlap the benchmark more
  often, because the previous audit scanned a different set of files, or because the two
  audits used different matching rules.

Until that reconciliation is complete, the accurate interpretation is: *a newer audit reports
potential PEPBI sequence overlap for 125 adaptive candidates, but the count is not yet
reconciled or verified at the sequence-and-affinity-record level.*

If these matches can be reproduced and linked to usable affinity data, they may enable a more
direct retrospective comparison of the two candidate sets. Even then, sequence overlap alone
does not prove better binding; the matched affinity distributions must be compared under
consistent target/receptor contexts, with duplicate handling and uncertainty reported.

## 5. Does the present evidence show translation to real-world binding?

No. It supports internal optimizer improvement and a change in candidate diversity, but
candidate-level experimental superiority has not yet been established. The external benchmarks
mainly test whether the scoring or mutation rules track measured affinity in retrospective
datasets. They are not equivalent to experimentally testing the candidates produced by the
adaptive and baseline arms.

The current data support the following conclusion:

> Across three matched seeds, adaptive guidance improved the optimizer's internal-score
> outcomes and increased the number of unique candidate sequences relative to the uniform
> baseline. However, retrospective experimental-affinity benchmarks did not establish that
> the native sequence score reliably predicts binding affinity, and the structural benchmark
> was limited by a small mapped sample and a pose-scoring concern. A newly supplied PEPBI
> audit reports 125 adaptive exact-sequence matches but conflicts with an earlier audit and
> lacks a visible matched-sequence list. Therefore, whether the optimizer's observed advantage
> translates into improved experimental binding remains unresolved.

## 6. What should happen next

Do not rerun the three-seed optimizer comparison merely to generate another version of the
same result. First reconcile the PEPBI match count and establish which exact candidate sequences
map to independent affinity records. This is the most immediate evidence gap raised by the
new spreadsheet.

After that, report the arm-specific matched affinity results only if the source records support
a valid comparison. Keep internal optimizer scores, retrospective mutation-prediction metrics,
structural proxy scores, and actual experimental binding measurements in separate result
categories. If no suitable affinity records can be linked to generated candidates, state that
candidate-level experimental validation is unavailable rather than treating a sequence match as
a successful binding result.
"""

(OUT / "logbook_DATA_RESULTS.md").write_text(data_md)
(OUT / "logbook_ANALYSIS.md").write_text(analysis_md)

readme = f"""# P3 full comparison audit

Created: {datetime.now().astimezone().isoformat(timespec="seconds")}

## Key files
- `logbook_DATA_RESULTS.md`: separate logbook data/results entry
- `logbook_ANALYSIS.md`: corresponding interpretation and limitations
- `source_file_inventory.csv`: source provenance, timestamps, hashes, snapshot paths
- `source_snapshot/`: copies of source files up to the size limit
- `json_scalar_inventory.csv`: scalar values extracted from JSON outputs
- `candidate_records_extracted.csv`: recognized sequence-bearing records from paired run
- `candidate_extraction_summary.csv`: parser coverage summary
- `paired_optimizer_seed_metrics.csv`: matched seed-level optimizer summary
- `pepbi_exact_match_audit_provisional.csv`: counts from newly supplied Numbers worksheet
- `external_validation_summary.csv`: external correlation/direction metrics
- `ridge_vs_mean_summary.csv`: grouped model comparison
- `pepbi_strict_matched_summary.csv`: strict-matched concordance results
- `graphs/`: figures generated by this audit

## Important limitation
The PEPBI exact-match figures are transcribed from the newly supplied Numbers worksheet:
baseline 0/275, adaptive 125/450. They conflict with an earlier audit that reported zero
matches in both arms. They remain provisional until source-level sequence matching reproduces
the actual matched sequences and corresponding affinity records.
"""
(OUT / "README.md").write_text(readme)

print("\nP3 FULL COMPARISON AUDIT FINISHED")
print("Selected paired run:", paired)
print("Comparison report:", comparison_report or "NOT FOUND")
print("Output directory:", OUT)
print("Files inventoried:", len(manifest))
print("Copied snapshots:", sum(1 for x in manifest if x["copied_to_snapshot"]))
print("JSON scalar leaves extracted:", len(flat_rows))
print("Sequence-bearing records extracted:", len(candidate_rows))
print("Graph generation:", "OK" if graphs_ok else "FAILED; see GRAPH_GENERATION_ERROR.txt")
print("\nMain outputs:")
for name in [
    "logbook_DATA_RESULTS.md",
    "logbook_ANALYSIS.md",
    "source_file_inventory.csv",
    "json_scalar_inventory.csv",
    "candidate_records_extracted.csv",
    "candidate_extraction_summary.csv",
    "pepbi_exact_match_audit_provisional.csv",
    "README.md",
]:
    print(" -", OUT / name)
print("Graphs:", graphs)
