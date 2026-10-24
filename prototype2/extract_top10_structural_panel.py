#!/usr/bin/env python3
"""Freeze top 10 unique sequence-level candidates per method and make a 5-seed pose manifest."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re, sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "prototype2" / "p2_uniform_random_reset_vs_v2_100_20261023"
OUT = ROOT / "prototype2" / "structural_validation"
AA = set("ACDEFGHIKLMNPQRSTVWY")
METHODS = ("adaptive_substitution", "uniform_baseline")
DOCK_SEEDS = (81001, 81002, 81003, 81004, 81005)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="replace only this script's three output files")
    args = ap.parse_args()
    runs = EXPERIMENT / "runs"
    if not runs.is_dir():
        raise SystemExit(f"ERROR: missing final experiment runs: {runs}")
    output_files = [OUT / "selected_top10.csv", OUT / "selection_manifest.json", OUT / "pose_manifest.csv"]
    if not args.force and any(p.exists() for p in output_files):
        raise SystemExit("ERROR: structural-selection outputs already exist. Inspect them; rerun with --force only to intentionally replace these three files.")

    files = sorted(runs.rglob("top_candidates.csv"))
    if len(files) != 200:
        raise SystemExit(f"ERROR: expected 200 per-run top_candidates.csv files, found {len(files)}; not selecting from an incomplete experiment.")
    records = []
    source_rows = Counter()
    method_files = Counter()
    source_hashes = {}
    for path in files:
        parts = path.parts
        method = next((m for m in METHODS if any(m in part for part in parts)), None)
        if method is None:
            raise SystemExit(f"ERROR: cannot determine method from path: {path}")
        seed_match = re.search(r"seed(\d+)", "/".join(parts))
        if not seed_match:
            raise SystemExit(f"ERROR: cannot determine source seed from path: {path}")
        seed = int(seed_match.group(1))
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            required = {"Sequence", "Overall Score", "Generation"}
            if not required.issubset(set(reader.fieldnames or [])):
                raise SystemExit(f"ERROR: unexpected CSV schema in {path}: {reader.fieldnames}")
            file_count = 0
            for row in reader:
                seq = (row.get("Sequence") or "").strip().upper()
                try: score = float(row["Overall Score"])
                except (TypeError, ValueError): raise SystemExit(f"ERROR: invalid score in {path}: {row}")
                if not seq or not math.isfinite(score):
                    raise SystemExit(f"ERROR: empty sequence/nonfinite score in {path}: {row}")
                if len(seq) != 21 or not set(seq).issubset(AA):
                    raise SystemExit(f"ERROR: unexpected peptide sequence alphabet/length in {path}: {seq!r}")
                try: generation = int(row["Generation"])
                except (TypeError, ValueError): generation = -1
                records.append({"method": method, "sequence": seq, "sequence_score": score,
                    "source_seed": seed, "source_run_id": path.parent.name,
                    "generation": generation, "source_file": str(path.relative_to(ROOT))})
                file_count += 1
            source_rows[method] += file_count
            method_files[method] += 1
        source_hashes[str(path.relative_to(ROOT))] = sha256(path)

    if method_files != Counter({"adaptive_substitution": 100, "uniform_baseline": 100}):
        raise SystemExit(f"ERROR: not exactly 100 saved method runs each: {method_files}")
    OUT.mkdir(parents=True, exist_ok=True)
    selected = []
    stats = {}
    for method in METHODS:
        pool = [r for r in records if r["method"] == method]
        # Deterministic tie-breaking; highest score per sequence is retained globally.
        pool.sort(key=lambda r: (-r["sequence_score"], r["sequence"], r["source_seed"], r["generation"], r["source_file"]))
        best_by_sequence = {}
        for r in pool:
            best_by_sequence.setdefault(r["sequence"], r)
        unique_ranked = list(best_by_sequence.values())
        unique_ranked.sort(key=lambda r: (-r["sequence_score"], r["sequence"], r["source_seed"], r["generation"], r["source_file"]))
        if len(unique_ranked) < 10:
            raise SystemExit(f"ERROR: only {len(unique_ranked)} unique valid sequences for {method}")
        stats[method] = {"input_rows": len(pool), "unique_sequences": len(unique_ranked), "top10": []}
        for rank, rec in enumerate(unique_ranked[:10], start=1):
            out = dict(rec); out["rank"] = rank
            selected.append(out)
            stats[method]["top10"].append({"rank": rank, "sequence": rec["sequence"], "sequence_score": rec["sequence_score"], "source_seed": rec["source_seed"], "generation": rec["generation"], "source_run_id": rec["source_run_id"]})

    fields = ["method", "rank", "sequence", "sequence_score", "source_seed", "source_run_id", "generation", "source_file"]
    write_csv(output_files[0], selected, fields)
    pose_rows = []
    for cand in selected:
        for rep, seed in enumerate(DOCK_SEEDS, start=1):
            rel = Path("prototype2/structural_validation/poses") / cand["method"] / f"rank_{cand['rank']:02d}" / f"seed_{seed}"
            pose_rows.append({"method": cand["method"], "rank": cand["rank"], "sequence": cand["sequence"],
                "sequence_score": cand["sequence_score"], "source_seed": cand["source_seed"],
                "repeat": rep, "vina_seed": seed, "output_dir": rel.as_posix()})
    write_csv(output_files[2], pose_rows, ["method", "rank", "sequence", "sequence_score", "source_seed", "repeat", "vina_seed", "output_dir"])
    summary_file = EXPERIMENT / "run_summary.csv"
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "experiment": str(EXPERIMENT.relative_to(ROOT)),
        "selection_rule": "Global descending Overall Score across each method's 100 saved top_candidates.csv files; exact sequence deduplication retaining the highest-scoring record; deterministic lexical/run tie-breaks; no docking information used.",
        "score_field": "Overall Score",
        "sequence_field": "Sequence",
        "expected_sequence_length": 21,
        "valid_alphabet": "ACDEFGHIKLMNPQRSTVWY",
        "source_top_candidates_file_count": len(files),
        "source_rows_by_method": dict(source_rows),
        "files_by_method": dict(method_files),
        "run_summary_sha256": sha256(summary_file) if summary_file.is_file() else None,
        "source_csv_sha256": source_hashes,
        "methods": stats,
        "docking_plan": {"backend": "AutoDock Vina only (exploratory alternative; not RAPiDock/HybriDock)", "independent_seeds_per_candidate": list(DOCK_SEEDS), "planned_runs": len(pose_rows), "primary_panel_size_per_method": 10},
        "scope_warning": "This is a selected top-10 panel geometric screen, not the geometric success rate of all generated candidates or experimental binding evidence. The panel is frozen on sequence-level score before docking.",
    }
    output_files[1].write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Source files: {len(files)} (100 per method)")
    print(f"Candidate rows: {len(records)}; source rows per method: {dict(source_rows)}")
    print(f"Selected: {len(selected)} unique candidates; planned docking runs: {len(pose_rows)}")
    for method in METHODS:
        print(f"\n{method}: {stats[method]['unique_sequences']} unique sequences")
        for rec in stats[method]["top10"]:
            print(f"{rec['rank']:2d}  {rec['sequence_score']:.6f}  {rec['sequence']}  seed={rec['source_seed']} gen={rec['generation']}")
    print(f"\nWrote: {output_files[0]}\nWrote: {output_files[1]}\nWrote: {output_files[2]}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
