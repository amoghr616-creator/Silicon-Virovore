from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json
from datetime import datetime, timezone

SOURCE = Path("prototype4/results/all_prototypes_csv_archive.zip")
OUTPUT = Path(
    "prototype4/results/AmphiMize_Teacher_Raw_Data_P1_P2_P3.zip"
)

P2_PREFIX = (
    "prototype2/p2_uniform_random_reset_vs_v2_100_20261023/"
)
P3_PREFIXES = [
    "prototype3/results/pepbi_corrected_benchmark/",
    "prototype3/results/pepbi_native_score_exploratory/",
    "prototype3/results/pepbi_strict_matched/",
]
P1_FILE = "p1graphdata/P1_ablation_summary.csv"

if not SOURCE.exists():
    raise SystemExit(f"Source archive not found: {SOURCE}")

readme = """AmphiMize — Teacher-Facing Raw Data Archive
================================================

P1 — BASELINE
-------------
Includes the available P1 ablation summary.
The P1 section contains the available baseline ablation summary.
The package documents the scope of the included P1 evidence.

P2 — ADAPTIVE OPTIMIZATION
-------------------------
Includes the 100-seed adaptive-versus-uniform comparison:
- 100 adaptive candidate-run tables
- 100 uniform-baseline candidate-run tables
- Paired seed differences
- Recalculated 95% confidence intervals
- Run summary

The candidate-run tables contain computational outputs.
They are not experimental binding measurements.

P3 — EXTERNAL PEPBI VALIDATION
------------------------------
Includes the relevant PEPBI benchmark records:
- Experimental binding measurements and Delta-G values
- Native-score association results
- Group-level Spearman correlation summaries
- Per-group correlation results
- Strictly matched validation tables
- Score repeatability data

SCOPE AND LIMITATIONS
---------------------
PEPBI is an external peptide-protein benchmark. It is not
direct experimental validation of HERV-K Env binding.

This archive contains selected research data, not the complete
repository or software environment.

Generated UTC: {created}
""".format(created=datetime.now(timezone.utc).isoformat())


with ZipFile(SOURCE, "r") as zin:
    available = set(zin.namelist())

    p2_files = sorted(
        n for n in available
        if n.startswith(P2_PREFIX)
        and n.lower().endswith((".csv", ".tsv"))
    )

    p3_files = sorted(
        n for n in available
        if any(n.startswith(prefix) for prefix in P3_PREFIXES)
        and n.lower().endswith((".csv", ".tsv"))
    )

    p1_exists = P1_FILE in available

    required_p2 = [
        P2_PREFIX + "paired_seed_differences.csv",
        P2_PREFIX + "paired_stats_recalculated_95ci.csv",
        P2_PREFIX + "run_summary.csv",
    ]

    required_p3 = [
        "prototype3/results/pepbi_corrected_benchmark/correlation_results.csv",
        "prototype3/results/pepbi_corrected_benchmark/pepbi_all_rows_with_native_scores.csv",
        "prototype3/results/pepbi_corrected_benchmark/group_level_correlation_audit/per_group_spearman.csv",
    ]

    missing = [
        p for p in required_p2 + required_p3
        if p not in available
    ]
    if missing:
        raise SystemExit("Required source files missing:\n" +
                         "\n".join(missing))

    if not p2_files or not p3_files:
        raise SystemExit("P2 or P3 source group is unexpectedly empty.")

    selected = []

    for name in p2_files:
        relative = name[len(P2_PREFIX):]
        selected.append((
            name,
            "P2_100_seed_comparison/" + relative
        ))

    for name in p3_files:
        relative = name.split(
            "prototype3/results/", 1
        )[1]
        selected.append((
            name,
            "P3_PEPBI_external_validation/" + relative
        ))

    if p1_exists:
        selected.append((
            P1_FILE,
            "P1_available_ablation_summary/P1_ablation_summary.csv"
        ))

    manifest_files = []

    with ZipFile(OUTPUT, "w", ZIP_DEFLATED) as zout:
        for source_name, archive_name in selected:
            data = zin.read(source_name)
            zout.writestr(archive_name, data)
            manifest_files.append({
                "archive_path": archive_name,
                "source_path": source_name,
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            })

        zout.writestr("README.txt", readme)

        manifest = {
            "archive_name": OUTPUT.name,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "source_archive": str(SOURCE),
            "p1_median_experimental_delta_g_found": False,
            "p1_ablation_summary_included": p1_exists,
            "p2_csv_tsv_count": len(p2_files),
            "p3_pepbi_csv_tsv_count": len(p3_files),
            "data_files": manifest_files,
        }

        zout.writestr(
            "SHA256_MANIFEST.json",
            json.dumps(manifest, indent=2)
        )

with ZipFile(OUTPUT, "r") as check:
    bad = check.testzip()
    if bad:
        raise SystemExit(f"ZIP integrity check failed: {bad}")
    total = len(check.namelist())

print("\nARCHIVE CREATED AND VERIFIED")
print("Path:", OUTPUT.resolve())
print("P1 ablation summary included:", p1_exists)
print("P1 median experimental Delta-G: NOT FOUND")
print("P2 files:", len(p2_files))
print("P3 PEPBI files:", len(p3_files))
print("Total archive entries:", total)
print("ZIP size: %.2f MB" % (OUTPUT.stat().st_size / 1_000_000))
print("Integrity check: PASS")
