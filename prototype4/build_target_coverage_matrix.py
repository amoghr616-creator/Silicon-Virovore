import csv
import json
from collections import Counter
from pathlib import Path

REPORT = Path(
    "results/p4_live_gate_smoke_20261102/"
    "p4_live_gate_smoke_20261102/"
    "p4_integrated_assessment/generation_0001.json"
)

OUT = Path(
    "prototype4/results/"
    "live_target_coverage_audit_20261102"
)


def main():
    if not REPORT.is_file():
        raise FileNotFoundError(REPORT)

    report = json.loads(REPORT.read_text())

    manifest_path = Path(
        report["target_manifest"]["manifest_path"]
    )

    target_report_path = Path(
        report["source_provenance"]["target_report_path"]
    )

    if not target_report_path.is_absolute():
        target_report_path = Path.cwd() / target_report_path

    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)

    if not target_report_path.is_file():
        raise FileNotFoundError(target_report_path)

    manifest = json.loads(manifest_path.read_text())
    target_report = json.loads(target_report_path.read_text())

    complexes = manifest.get("complexes", [])
    reference_counts = Counter(
        item["sequence"].strip().upper()
        for item in complexes
        if isinstance(item.get("sequence"), str)
    )

    if not reference_counts:
        raise RuntimeError("No peptide sequences found in manifest.")

    target_candidates = target_report.get("candidates", [])
    integrated_candidates = report.get("candidates", [])

    if len(target_candidates) != len(integrated_candidates):
        raise RuntimeError(
            "Candidate count mismatch between target and integrated reports."
        )

    rows = []

    for target, integrated in zip(
        target_candidates, integrated_candidates
    ):
        sequence = target.get("sequence")

        if not isinstance(sequence, str):
            raise RuntimeError(
                f"Missing sequence for candidate "
                f"{integrated.get('candidate_id')}"
            )

        sequence = sequence.strip().upper()

        distances = {}

        for reference in reference_counts:
            if len(reference) != len(sequence):
                continue

            differences = [
                (i + 1, ref_aa, candidate_aa)
                for i, (ref_aa, candidate_aa)
                in enumerate(zip(reference, sequence))
                if ref_aa != candidate_aa
            ]

            distances[reference] = differences

        if not distances:
            raise RuntimeError(
                f"No equal-length reference sequence for "
                f"{integrated.get('candidate_id')}"
            )

        nearest_reference = min(
            distances,
            key=lambda ref: len(distances[ref])
        )

        substitutions = distances[nearest_reference]
        exact_match = len(substitutions) == 0

        rows.append({
            "candidate_id": integrated.get("candidate_id"),
            "candidate_index": integrated.get("candidate_index"),
            "sequence": sequence,
            "target_status": integrated.get("target_status"),
            "target_eligible": integrated.get("target_eligible"),
            "exact_manifest_match": exact_match,
            "nearest_reference_sequence": nearest_reference,
            "minimum_hamming_distance": len(substitutions),
            "substitutions": ";".join(
                f"{old}{position}{new}"
                for position, old, new in substitutions
            ) or "NONE",
            "matching_complex_count": (
                reference_counts[sequence]
                if exact_match else 0
            ),
        })

    exact_count = sum(
        row["exact_manifest_match"] for row in rows
    )

    OUT.mkdir(parents=True, exist_ok=True)

    csv_path = OUT / "candidate_sequence_coverage.csv"
    json_path = OUT / "candidate_sequence_coverage.json"

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys())
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "source_report": str(REPORT),
        "source_manifest": str(manifest_path),
        "candidate_count": len(rows),
        "distinct_manifest_sequences": len(reference_counts),
        "manifest_complex_count": len(complexes),
        "exact_sequence_matches": exact_count,
        "unmatched_candidates": len(rows) - exact_count,
        "exact_sequence_coverage": exact_count / len(rows),
        "interpretation": (
            "Exact-sequence coverage only. Near matches are not "
            "treated as structural evidence for variant peptides."
        ),
        "candidates": rows,
    }

    json_path.write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print("========== TARGET COVERAGE MATRIX ==========")
    print(f"Candidates: {len(rows)}")
    print(f"Distinct manifest sequences: {len(reference_counts)}")
    print(f"Manifest complexes: {len(complexes)}")
    print(f"Exact matches: {exact_count}")
    print(f"Unmatched: {len(rows) - exact_count}")
    print(f"Exact coverage: {exact_count / len(rows):.2%}")

    for row in rows:
        print(
            f"{row['candidate_id']}: "
            f"distance={row['minimum_hamming_distance']}, "
            f"substitutions={row['substitutions']}, "
            f"status={row['target_status']}"
        )

    print(f"\nCSV: {csv_path}")
    print(f"JSON: {json_path}")


if __name__ == "__main__":
    main()
