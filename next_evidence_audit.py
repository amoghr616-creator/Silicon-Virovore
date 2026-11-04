#!/usr/bin/env python3
"""Read-only inventory of benchmark tables, HADDOCK configs, and pose archives."""

from pathlib import Path
import csv
import difflib
import hashlib
import json
import re
import tarfile

ROOT = Path.cwd()
SEARCH_ROOTS = [
    ROOT / "prototype3" / "results",
    ROOT / "prototype4" / "results",
    ROOT / "prototype4" / "configs",
]

KEYWORDS = (
    "ledger", "generation", "cycle", "skempi", "loco", "pepbi",
    "ppi", "manifest", "run_plan", "target_evidence", "reliability",
    "fingerprint", "pose", "experiment", "seed", "validation",
)

lines = []

def emit(message=""):
    lines.append(str(message))

def readable_text(path, max_bytes=5_000_000):
    try:
        if path.stat().st_size > max_bytes:
            return None
        return path.read_text(encoding="utf-8", errors="replace")
    except (OSError, UnicodeError):
        return None

def inspect_table(path):
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, [])
            row_count = sum(1 for _ in reader)
        return header, row_count
    except (OSError, csv.Error, UnicodeError) as exc:
        return [f"READ_ERROR: {exc}"], -1

def inspect_json(path):
    text = readable_text(path)
    if text is None:
        return "too large or unreadable"
    try:
        obj = json.loads(text)
    except json.JSONDecodeError as exc:
        return f"invalid JSON: {exc}"
    if isinstance(obj, dict):
        return f"dict keys={list(obj.keys())[:35]}"
    if isinstance(obj, list):
        first = obj[0] if obj else None
        if isinstance(first, dict):
            return f"list length={len(obj)}; first-row keys={list(first.keys())[:35]}"
        return f"list length={len(obj)}; first item type={type(first).__name__}"
    return f"root type={type(obj).__name__}"

def normalized_config(text):
    # Ignore only the HADDOCK random seed; retain every other setting.
    return re.sub(
        r"(?im)^(\s*iniseed\s*[:=]\s*)[^\r\n]+",
        r"\1<SEED>",
        text,
    )

emit("SILICON-VIROVORE NEXT EVIDENCE AUDIT")
emit(f"Repository: {ROOT}")
emit("=" * 78)

# 1. Data files and schemas
emit("\n[1] RELEVANT TABLES AND JSON FILES")
data_files = []
for search_root in SEARCH_ROOTS:
    if not search_root.exists():
        continue
    for path in search_root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".csv", ".tsv", ".json"}:
            continue
        searchable = str(path).lower()
        if any(word in searchable for word in KEYWORDS):
            data_files.append(path)

for path in sorted(set(data_files))[:100]:
    rel = path.relative_to(ROOT)
    if path.suffix.lower() in {".csv", ".tsv"}:
        header, count = inspect_table(path)
        emit(f"\nTABLE: {rel}")
        emit(f"  rows={count}; columns={len(header)}")
        emit(f"  header={header}")
    else:
        emit(f"\nJSON: {rel}")
        emit(f"  {inspect_json(path)}")

if not data_files:
    emit("No matching CSV, TSV, or JSON files found under the search roots.")

# 2. HADDOCK configs
emit("\n[2] HADDOCK CONFIGURATION COMPARISON")
configs = []
for search_root in SEARCH_ROOTS:
    if not search_root.exists():
        continue
    for path in search_root.rglob("*.cfg"):
        text = readable_text(path)
        if text is None:
            continue
        searchable = (str(path) + "\n" + text).lower()
        if "haddock" in searchable or "iniseed" in searchable:
            configs.append((path, text))

config_records = []
for path, text in sorted(configs, key=lambda item: str(item[0])):
    normalized = normalized_config(text)
    digest = hashlib.sha256(normalized.encode("utf-8")).hexdigest()
    config_records.append((path, text, digest))
    emit(f"\nCONFIG: {path.relative_to(ROOT)}")
    emit(f"  normalized_sha256={digest}")

unique_hashes = sorted({record[2] for record in config_records})
emit(f"\nUnique normalized config hashes: {len(unique_hashes)}")
if len(config_records) >= 2:
    baseline_path, baseline_text, baseline_hash = config_records[0]
    for path, text, digest in config_records[1:]:
        if digest == baseline_hash:
            continue
        diff = list(difflib.unified_diff(
            baseline_text.splitlines(),
            text.splitlines(),
            fromfile=str(baseline_path.relative_to(ROOT)),
            tofile=str(path.relative_to(ROOT)),
            lineterm="",
        ))
        emit(f"\nCONFIG DIFF versus {baseline_path.relative_to(ROOT)}:")
        emit("\n".join(diff[:100]))
        if len(diff) > 100:
            emit(f"... diff truncated; {len(diff)} total lines")

# 3. Loose final structures and archived structures
emit("\n[3] PDB STRUCTURE INVENTORY")
loose_pdbs = []
archives = []
for search_root in SEARCH_ROOTS:
    if not search_root.exists():
        continue
    loose_pdbs.extend(search_root.rglob("*.pdb"))
    for path in search_root.rglob("*"):
        if path.is_file() and (
            path.name.lower().endswith((".tar", ".tar.gz", ".tgz", ".tar.bz2"))
        ):
            archives.append(path)

emit(f"Loose PDB files: {len(loose_pdbs)}")
for path in sorted(loose_pdbs)[:30]:
    emit(f"  {path.relative_to(ROOT)}")

for archive in sorted(set(archives)):
    try:
        with tarfile.open(archive, "r:*") as tar:
            members = [
                member.name for member in tar.getmembers()
                if member.isfile() and member.name.lower().endswith(".pdb")
            ]
        if not members:
            continue
        emit(f"\nARCHIVE: {archive.relative_to(ROOT)}")
        emit(f"  PDB members={len(members)}")
        final_like = [
            name for name in members
            if any(token in name.lower() for token in (
                "emref", "clustfcc", "seletopclust", "selectopclust",
                "final", "cluster_1", "model",
            ))
        ]
        emit(f"  possible final/refined members={len(final_like)}")
        for name in (final_like or members)[:20]:
            emit(f"    {name}")
    except (OSError, tarfile.TarError) as exc:
        emit(f"\nARCHIVE ERROR: {archive.relative_to(ROOT)}: {exc}")

# 4. Confirm the expected repeat folders exist
emit("\n[4] EXPECTED HADDOCK REPEAT LOCATIONS")
expected = [
    ROOT / "prototype4" / "results" / "haddock3_chain_A_pilot_20261031",
    ROOT / "prototype4" / "results" / "p4_reliability_repeats_20261103",
]
for path in expected:
    emit(f"{path.relative_to(ROOT)}: {'FOUND' if path.exists() else 'NOT FOUND'}")
    if path.exists() and path.is_dir():
        children = sorted(p.name for p in path.iterdir())
        emit(f"  contents={children[:40]}")

report = "\n".join(lines) + "\n"
out_dir = ROOT / "prototype4" / "results"
out_dir.mkdir(parents=True, exist_ok=True)
report_path = out_dir / "next_evidence_audit.txt"
report_path.write_text(report, encoding="utf-8")

print(report)
print(f"\nSaved report: {report_path}")
