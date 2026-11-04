import csv
import json
from pathlib import Path

ROOT = Path("results/unified_ablation_20261102")
PILOT = ROOT / "pilot2"
POOL = ROOT / "initial_pools/seed_631.fasta"
OUTPUT = ROOT / "pilot2_summary.txt"

def read_fasta(path):
    sequences = []
    current = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
            current = []
        else:
            current.append(line.upper())
    if current:
        sequences.append("".join(current))
    return sequences

def walk(obj, path=""):
    if isinstance(obj, dict):
        for key, value in obj.items():
            current = f"{path}.{key}" if path else key
            if isinstance(value, (dict, list)):
                yield from walk(value, current)
            else:
                yield current, key, value
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            yield from walk(value, f"{path}[{index}]")

def print_selected(section, title):
    wanted = {
        "run_id", "seed", "random_seed", "population_size",
        "generations", "generations_completed", "generation_count",
        "mutation_guidance_mode", "requested_mutation_guidance_mode",
        "adaptive_enabled", "arise_enabled", "ale_enabled",
        "hotspot_enabled", "cache_enabled", "initial_pool_sha256",
        "candidate_count", "total_candidates", "evaluated",
        "unique_sequences", "best_score", "best_overall",
        "mean_score", "mean_overall", "runtime_seconds",
        "target_evidence_coverage", "target_evidence_count",
        "pareto_reordered_candidates", "pareto_requested",
        "status", "run_status", "actual_evaluations",
    }
    print(f"\n{title}")
    found = False
    for full_path, key, value in walk(section):
        if key in wanted:
            print(f"  {full_path}: {value}")
            found = True
    if not found:
        print("  No matching scalar fields found.")

def main():
    initial = read_fasta(POOL)
    initial_set = set(initial)

    lines = []
    import contextlib
    import io

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        print("UNIFIED P1-P4 PILOT SUMMARY")
        print("Initial FASTA records:", len(initial))
        print("Initial unique sequences:", len(initial_set))

        run_dirs = sorted(
            p for p in PILOT.iterdir()
            if p.is_dir() and p.name.startswith("unified_pilot2_s631_")
        )

        population_results = {}

        for run_dir in run_dirs:
            arm = run_dir.name.rsplit("_", 1)[-1]
            print("\n" + "=" * 72)
            print("ARM:", arm)
            print("DIRECTORY:", run_dir)

            audit_path = run_dir / "experiment_audit.json"
            if audit_path.exists():
                audit = json.loads(audit_path.read_text())
                for section in ("configuration", "execution",
                                "optimization", "evidence"):
                    if section in audit:
                        print_selected(audit[section], section.upper())
            else:
                print("Missing experiment audit.")

            for generation in (1, 2):
                target_path = (
                    run_dir / "p4_target_evidence"
                    / f"generation_{generation:04d}.json"
                )
                if target_path.exists():
                    report = json.loads(target_path.read_text())
                    print(f"\nTARGET EVIDENCE G{generation}")
                    for key in (
                        "candidate_count",
                        "target_evidence_count",
                        "target_evidence_coverage",
                        "pareto_requested",
                        "pareto_reordered_candidates",
                        "reason",
                    ):
                        print(f"  {key}: {report.get(key)}")
                    print("  manifest:", report.get("manifest"))

            csv_path = (
                run_dir / "p4_integrated_assessment"
                / "generation_0001.csv"
            )
            if csv_path.exists():
                with csv_path.open(newline="") as handle:
                    reader = csv.DictReader(handle)
                    headers = reader.fieldnames or []
                    rows = list(reader)

                sequence_column = next(
                    (h for h in headers if h.lower() in {
                        "sequence", "candidate_sequence", "peptide_sequence"
                    }),
                    None,
                )

                print("\nINITIAL POPULATION CSV")
                print("  rows:", len(rows))
                print("  columns:", headers)

                if sequence_column:
                    observed = [
                        row[sequence_column].strip().upper()
                        for row in rows
                        if row.get(sequence_column)
                    ]
                    observed_set = set(observed)
                    matches = (
                        len(observed) == 100
                        and len(observed_set) == 100
                        and observed_set == initial_set
                    )
                    population_results[arm] = matches
                    print("  sequence column:", sequence_column)
                    print("  unique sequences:", len(observed_set))
                    print("  exact FASTA match:", matches)
                else:
                    print("  No recognizable sequence column.")
            else:
                print("Missing generation-1 CSV.")

        print("\n" + "=" * 72)
        print("SHARED INITIAL POOL VERIFICATION")
        for arm, result in sorted(population_results.items()):
            print(f"  {arm}: exact match = {result}")

        print(
            "ALL FOUR MATCHED:",
            len(population_results) == 4
            and all(population_results.values()),
        )

    text = buffer.getvalue()
    OUTPUT.write_text(text)
    print(text)
    print(f"\nSaved report: {OUTPUT.resolve()}")

if __name__ == "__main__":
    main()
