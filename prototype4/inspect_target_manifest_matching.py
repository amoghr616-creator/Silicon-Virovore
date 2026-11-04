import csv
import json
import re
from pathlib import Path

REPORT = Path(
    "results/p4_live_gate_smoke_20261102/"
    "p4_live_gate_smoke_20261102/"
    "p4_integrated_assessment/generation_0001.json"
)

AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
SEQUENCE_KEYS = {
    "sequence", "candidate_sequence", "peptide_sequence",
    "aa_sequence", "sequence_string",
}


def valid_sequence(value):
    if not isinstance(value, str):
        return None
    value = re.sub(r"\s+", "", value).upper()
    return value if AA.fullmatch(value) else None


def collect_sequences(obj, path="root"):
    found = []

    if isinstance(obj, dict):
        for key, value in obj.items():
            current = f"{path}.{key}"

            if (
                key.lower() in SEQUENCE_KEYS
                or "peptide_sequence" in key.lower()
                or "candidate_sequence" in key.lower()
            ):
                seq = valid_sequence(value)
                if seq:
                    found.append((current, seq))

            found.extend(collect_sequences(value, current))

    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            found.extend(collect_sequences(value, f"{path}[{i}]"))

    return found


def candidate_ids(obj):
    result = []

    for item in obj.get("candidates", []):
        result.append({
            "id": item.get("candidate_id"),
            "index": item.get("candidate_index"),
            "status": item.get("target_status"),
            "eligible": item.get("target_eligible"),
            "sequences": collect_sequences(item, "candidate"),
        })

    return result


def scan_json_file(path):
    try:
        if path.stat().st_size > 25_000_000:
            return []
        return collect_sequences(json.loads(path.read_text()))
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return []


def scan_csv_file(path):
    found = []

    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                return found

            seq_columns = [
                col for col in reader.fieldnames
                if "sequence" in col.lower()
            ]

            for row_number, row in enumerate(reader, start=2):
                for col in seq_columns:
                    seq = valid_sequence(row.get(col))
                    if seq:
                        found.append(
                            (f"{path.name}:row{row_number}.{col}", seq)
                        )

    except (OSError, UnicodeDecodeError, csv.Error):
        pass

    return found


def main():
    if not REPORT.is_file():
        raise FileNotFoundError(REPORT)

    report = json.loads(REPORT.read_text())
    manifest_path = Path(
        report["target_manifest"]["manifest_path"]
    )

    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)

    manifest = json.loads(manifest_path.read_text())
    manifest_sequences = collect_sequences(manifest, "manifest")

    print("========== CANDIDATE RECORDS ==========")

    records = candidate_ids(report)

    for record in records:
        print(
            f"\nIndex={record['index']} "
            f"ID={record['id']} "
            f"target_status={record['status']} "
            f"target_eligible={record['eligible']}"
        )

        for location, seq in record["sequences"]:
            print(f"  Sequence at {location}: {seq}")

    print("\n========== MANIFEST SEQUENCES ==========")
    print(f"Manifest: {manifest_path}")
    print(f"Sequence fields found: {len(manifest_sequences)}")

    manifest_set = set()

    for location, seq in manifest_sequences:
        manifest_set.add(seq)
        shown = seq if len(seq) <= 100 else seq[:100] + "..."
        print(f"  {location}: length={len(seq)} sequence={shown}")

    print("\n========== RUN FILE INVENTORY ==========")

    run_root = REPORT.parent.parent
    scan_roots = [run_root, run_root.parent]
    files = set()

    for root in scan_roots:
        if root.exists():
            for path in root.rglob("*"):
                if (
                    path.is_file()
                    and path.suffix.lower() in {".json", ".jsonl", ".csv"}
                ):
                    files.add(path)

    for path in sorted(files):
        print(path)

    print("\n========== CANDIDATE SEQUENCE SEARCH ==========")

    file_sequence_hits = []

    for path in sorted(files):
        if path.suffix.lower() == ".csv":
            hits = scan_csv_file(path)
        elif path.suffix.lower() == ".json":
            hits = scan_json_file(path)
        else:
            continue

        if hits:
            for location, seq in hits:
                file_sequence_hits.append((str(path), location, seq))
                match = seq in manifest_set
                print(
                    f"{path} | {location} | length={len(seq)} "
                    f"| exact_manifest_match={match} | {seq}"
                )

    if not file_sequence_hits:
        print("No candidate sequence fields found in scanned JSON/CSV files.")

    print("\n========== SUMMARY ==========")
    print(f"Candidates in report: {len(records)}")
    print(f"Distinct manifest sequences: {len(manifest_set)}")
    print(f"Candidate sequence records found: {len(file_sequence_hits)}")
    print("No files were modified.")


if __name__ == "__main__":
    main()
