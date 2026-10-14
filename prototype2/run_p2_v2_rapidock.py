from pathlib import Path
import argparse
import csv
import subprocess
import sys

ROOT = Path.cwd()
PANEL = ROOT / "prototype2" / "p2_v2_panel.tsv"
RECEPTOR = ROOT / "prototype2" / "9SYA_chainA_final.pdb"
OUTROOT = ROOT / "prototype2" / "p2_v2_rapidock"
SITE = (-8.405, -5.823, 1.523)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int, default=100)
    parser.add_argument("--limit-per-method", type=int)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    if args.n_samples < 1:
        sys.exit("--n-samples must be positive.")

    if not PANEL.exists():
        sys.exit(f"Panel not found: {PANEL}")
    if not RECEPTOR.exists():
        sys.exit(f"Receptor not found: {RECEPTOR}")

    with PANEL.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    required = {
        "candidate_id", "method", "sequence",
        "generation_local_score", "run_id",
    }
    if not rows or not required.issubset(rows[0]):
        sys.exit(f"Panel is empty or missing columns: {required}")

    selected = []
    counts = {}
    for row in rows:
        method = row["method"]
        counts.setdefault(method, 0)

        if (
            args.limit_per_method is not None
            and counts[method] >= args.limit_per_method
        ):
            continue

        counts[method] += 1
        selected.append(row)

    succeeded, skipped, failed = [], [], []

    for row in selected:
        candidate_id = row["candidate_id"]
        sequence = row["sequence"].strip().upper()
        outdir = OUTROOT / candidate_id
        ranked_file = outdir / "ranked_poses.csv"

        if args.resume and ranked_file.exists():
            print(f"SKIP completed: {candidate_id}")
            skipped.append(candidate_id)
            continue

        cmd = [
            "hybridock-pep", "dock",
            "--peptide", sequence,
            "--receptor", str(RECEPTOR),
            "--site", *(str(x) for x in SITE),
            "--box", "52",
            "--n-samples", str(args.n_samples),
            "--output-dir", str(outdir),
        ]

        print("\n" + "=" * 72)
        print(f"{candidate_id} | {row['method']} | {sequence}")
        print(" ".join(cmd))

        if args.dry_run:
            continue

        outdir.mkdir(parents=True, exist_ok=True)

        try:
            result = subprocess.run(cmd, check=False)
        except FileNotFoundError:
            sys.exit(
                "Could not find hybridock-pep. Check the active environment "
                "and installation."
            )

        if result.returncode == 0 and ranked_file.exists():
            print(f"COMPLETED: {candidate_id}")
            succeeded.append(candidate_id)
        else:
            print(f"FAILED or missing ranked_poses.csv: {candidate_id}")
            failed.append(candidate_id)

    print("\nEvaluation batch finished.")
    print(f"Candidates selected: {len(selected)}")
    print(f"Completed: {len(succeeded)}")
    print(f"Skipped as completed: {len(skipped)}")
    print(f"Failed: {len(failed)}")
    if failed:
        print("Failed IDs:", ", ".join(failed))


if __name__ == "__main__":
    main()
