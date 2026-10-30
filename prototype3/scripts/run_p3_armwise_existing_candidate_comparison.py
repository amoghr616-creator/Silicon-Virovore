#!/usr/bin/env python3
"""Armwise retrospective comparison using EXISTING candidate histories only.

Does not run candidate generation, the optimizer, or docking.
Tests existing candidates against available measured mutation/peptide data.
"""
from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "prototype3" / "results"
STAMP = datetime.now().strftime("%Y%m%dT%H%M%S")
OUT = RESULTS / f"p3_armwise_direct_comparison_{STAMP}"
OUT.mkdir(parents=True, exist_ok=False)


def clean_seq(value):
    if value is None:
        return ""
    value = str(value).strip().upper()
    if not value or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", value):
        return ""
    return value


def write_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_csv(path):
    try:
        with path.open(newline="", errors="replace") as f:
            return list(csv.DictReader(f))
    except (OSError, UnicodeError, csv.Error):
        return []


def num(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def norm_header(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def find_header(row, aliases):
    lut = {norm_header(k): k for k in row}
    for a in aliases:
        if norm_header(a) in lut:
            return lut[norm_header(a)]
    return None


def score_of(row):
    # Prefer the recorded native sequence-fitness value, not a mixed/docking score.
    for key in ("native_fitness", "c_score", "arise_comparable_score"):
        x = num(row.get(key))
        if x is not None:
            return x
    objectives = row.get("raw_objectives")
    if isinstance(objectives, dict):
        x = num(objectives.get("fitness"))
        if x is not None:
            return x
    return None


def locate_existing_run():
    candidates = []
    for runs in RESULTS.glob("p3_paired_arm_*/runs"):
        arm_seeds = {}
        for arm in ("adaptive", "baseline"):
            histories = list(runs.glob(f"{arm}/seed*/*/candidate_history.json"))
            arm_seeds[arm] = {
                p.parent.parent.name.replace("seed", "")
                for p in histories
            }
        common = arm_seeds.get("adaptive", set()) & arm_seeds.get("baseline", set())
        if common:
            candidates.append((len(common), runs.stat().st_mtime, runs, common))
    if not candidates:
        raise SystemExit(
            "ABORT: no existing paired adaptive/baseline candidate histories found "
            "under prototype3/results/p3_paired_arm_*/runs. No optimizer was launched."
        )
    # Prefer the largest complete paired cohort, then its newest directory.
    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
    _, _, runs, common = candidates[0]
    return runs, sorted(common)


RUNS_ROOT, PAIRED_SEEDS = locate_existing_run()
records = []

for arm in ("adaptive", "baseline"):
    for seed in PAIRED_SEEDS:
        histories = sorted((RUNS_ROOT / arm / f"seed{seed}").glob("*/candidate_history.json"))
        for history in histories:
            try:
                obj = json.loads(history.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            if not isinstance(obj, list):
                continue
            for i, row in enumerate(obj):
                if not isinstance(row, dict):
                    continue
                seq = clean_seq(row.get("sequence"))
                parent = clean_seq(row.get("parent_sequence"))
                if not seq:
                    continue
                score = score_of(row)
                changed = []
                if parent and len(parent) == len(seq):
                    changed = [j + 1 for j, (a, b) in enumerate(zip(parent, seq)) if a != b]
                records.append({
                    "arm": arm,
                    "seed": str(seed),
                    "candidate_id": row.get("candidate_id", f"{history.name}:{i}"),
                    "sequence": seq,
                    "parent_sequence": parent,
                    "length": len(seq),
                    "native_score": score,
                    "recorded_mutation_count": row.get("mutation_count", len(changed)),
                    "changed_positions": ";".join(map(str, changed)),
                    "n_changed_positions": len(changed),
                    "source_file": str(history.relative_to(ROOT)),
                })

# De-duplicate repeated serialization of the same candidate ID within the same arm/seed.
dedup = {}
for r in records:
    key = (r["arm"], r["seed"], str(r["candidate_id"]), r["sequence"])
    dedup[key] = r
records = list(dedup.values())

if not records:
    raise SystemExit("ABORT: candidate histories were found but no sequences were extractable.")

write_csv(OUT / "candidate_level_existing_records.csv", records)

# Test 3/4 supporting rows: use exact parent-child sequence pairs where possible.
pair_files = set()
for pattern in (
    "**/single_substitution_pairs.csv",
    "**/corrected_external_mutation_pairs.csv",
    "**/corrected_mutation_pairs.csv",
):
    pair_files.update(RESULTS.glob(pattern))

external_pairs = []
for path in sorted(pair_files):
    for row in read_csv(path):
        a_key = find_header(row, ("sequence_a", "peptide_a", "sequence_wt", "parent_sequence"))
        b_key = find_header(row, ("sequence_b", "peptide_b", "sequence_mutant", "mutant_sequence"))
        if not a_key or not b_key:
            continue
        a, b = clean_seq(row.get(a_key)), clean_seq(row.get(b_key))
        if not a or not b:
            continue
        pka_key = find_header(row, ("pKd_a", "pkd_a", "pka", "affinity_a"))
        pkb_key = find_header(row, ("pKd_b", "pkd_b", "pkb", "affinity_b"))
        delta_key = find_header(row, ("delta_pKd", "delta_pkd", "delta_affinity"))
        pka, pkb = num(row.get(pka_key)) if pka_key else None, num(row.get(pkb_key)) if pkb_key else None
        delta = (pkb - pka) if pka is not None and pkb is not None else (
            num(row.get(delta_key)) if delta_key else None
        )
        external_pairs.append({
            "source_file": str(path.relative_to(ROOT)),
            "sequence_a": a, "sequence_b": b,
            "length_a": len(a), "length_b": len(b),
            "delta_pKd_as_b_minus_a": delta,
            "context_key": row.get("context_key", ""),
            "target_group_id": row.get("target_group_id", ""),
            "target_value": row.get("target_value", ""),
        })

# Preserve separate assay/target contexts instead of silently collapsing them.
pair_lookup = defaultdict(list)
for p in external_pairs:
    pair_lookup[frozenset((p["sequence_a"], p["sequence_b"]))].append(p)

pair_matches = []
for r in records:
    parent, child = r["parent_sequence"], r["sequence"]
    if not parent or parent == child:
        continue
    for p in pair_lookup.get(frozenset((parent, child)), []):
        d = p["delta_pKd_as_b_minus_a"]
        if d is not None and parent == p["sequence_a"]:
            oriented_delta = d
        elif d is not None:
            oriented_delta = -d
        else:
            oriented_delta = None
        pair_matches.append({
            "arm": r["arm"], "seed": r["seed"],
            "candidate_id": r["candidate_id"],
            "parent_sequence": parent, "candidate_sequence": child,
            "candidate_length": len(child),
            "native_score": r["native_score"],
            "measured_delta_pKd_parent_to_child": oriented_delta,
            "external_context_key": p["context_key"],
            "external_target_group_id": p["target_group_id"],
            "external_source_file": p["source_file"],
            "target_context_verified_for_9SYA": "NO—must verify source target",
        })

write_csv(OUT / "exact_candidate_parent_child_experimental_matches.csv", pair_matches)

# Direct arm-by-arm sequence overlap with existing PEPBI tables.
pepbi_files = [
    p for p in RESULTS.rglob("*.csv")
    if "pepbi" in str(p).lower()
    and OUT.name not in str(p)
    and not any(part.startswith(".") for part in p.parts)
]
pepbi_hits = []
for path in sorted(set(pepbi_files)):
    for row in read_csv(path):
        # Only columns likely to contain peptide/ligand sequences.
        seq_cols = [
            k for k in row
            if any(t in norm_header(k) for t in ("sequence", "peptide", "ligand"))
            and not any(t in norm_header(k) for t in ("target", "protein", "receptor"))
        ]
        for col in seq_cols:
            seq = clean_seq(row.get(col))
            if not seq:
                continue
            for arm in ("adaptive", "baseline"):
                arms_seqs = {r["sequence"] for r in records if r["arm"] == arm}
                if seq in arms_seqs:
                    pepbi_hits.append({
                        "arm": arm,
                        "candidate_sequence": seq,
                        "source_file": str(path.relative_to(ROOT)),
                        "matched_column": col,
                        "source_row_context": json.dumps(
                            {k: row[k] for k in row if any(
                                t in norm_header(k) for t in ("kd", "pkd", "affinity", "delta", "assay", "target", "context")
                            )},
                            ensure_ascii=False,
                        ),
                    })
write_csv(OUT / "pepbi_exact_sequence_matches_by_arm.csv", pepbi_hits)

# Does an experimental OOF validation record contain any generated candidate sequence?
oof_files = list(RESULTS.rglob("out_of_fold_predictions.csv"))
oof_hits = []
for path in sorted(oof_files):
    for row in read_csv(path):
        seq_cols = [
            k for k in row
            if any(t in norm_header(k) for t in ("sequence", "peptide", "ligand"))
            and not any(t in norm_header(k) for t in ("target", "protein", "receptor"))
        ]
        for col in seq_cols:
            seq = clean_seq(row.get(col))
            if not seq:
                continue
            for arm in ("adaptive", "baseline"):
                if seq in {r["sequence"] for r in records if r["arm"] == arm}:
                    oof_hits.append({
                        "arm": arm, "candidate_sequence": seq,
                        "source_file": str(path.relative_to(ROOT)),
                        "matched_column": col,
                        "prediction_row": json.dumps(row, ensure_ascii=False),
                    })
write_csv(OUT / "ridge_oof_candidate_sequence_matches.csv", oof_hits)

# Arm-specific candidate scores, diversity, and mutation positions.
arm_summaries = []
position_counts = Counter()
position_denominators = Counter()
seed_summaries = []

for arm in ("adaptive", "baseline"):
    arm_rows = [r for r in records if r["arm"] == arm]
    scores = [r["native_score"] for r in arm_rows if r["native_score"] is not None]
    arm_sequences = {r["sequence"] for r in arm_rows}
    arm_summaries.append({
        "arm": arm,
        "candidate_records": len(arm_rows),
        "unique_sequences": len(arm_sequences),
        "sequence_diversity_unique_over_records": len(arm_sequences) / len(arm_rows) if arm_rows else None,
        "records_with_native_score": len(scores),
        "mean_native_score": sum(scores) / len(scores) if scores else None,
        "median_native_score": sorted(scores)[len(scores)//2] if scores else None,
        "best_native_score": max(scores) if scores else None,
        "mean_changed_positions_per_record": (
            sum(r["n_changed_positions"] for r in arm_rows) / len(arm_rows)
            if arm_rows else None
        ),
        "candidate_parent_child_pairs_with_exact_external_measurements": sum(
            r["arm"] == arm for r in pair_matches
        ),
        "exact_PEPBI_candidate_sequences": len({
            r["candidate_sequence"] for r in pepbi_hits if r["arm"] == arm
        }),
        "candidate_sequences_found_in_Ridge_OOF_rows": len({
            r["candidate_sequence"] for r in oof_hits if r["arm"] == arm
        }),
    })
    by_seed = defaultdict(list)
    for r in arm_rows:
        by_seed[r["seed"]].append(r)
        for pos in range(1, (r["length"] or 0) + 1):
            position_denominators[(arm, pos)] += 1
        for pos in map(int, filter(None, r["changed_positions"].split(";"))):
            position_counts[(arm, pos)] += 1
    for seed, sr in sorted(by_seed.items()):
        ss = [r["native_score"] for r in sr if r["native_score"] is not None]
        seed_summaries.append({
            "arm": arm, "seed": seed, "records": len(sr),
            "unique_sequences": len({r["sequence"] for r in sr}),
            "mean_native_score": sum(ss) / len(ss) if ss else None,
            "best_native_score": max(ss) if ss else None,
        })

write_csv(OUT / "armwise_summary.csv", arm_summaries)
write_csv(OUT / "armwise_seed_summary.csv", seed_summaries)

position_rows = []
for arm in ("adaptive", "baseline"):
    positions = sorted(p for a, p in position_denominators if a == arm)
    for pos in positions:
        denom = position_denominators[(arm, pos)]
        count = position_counts[(arm, pos)]
        position_rows.append({
            "arm": arm, "position_1based": pos,
            "records_with_position_comparison": denom,
            "records_mutated_at_position": count,
            "percent_records_mutated_at_position": 100 * count / denom if denom else None,
        })
write_csv(OUT / "armwise_mutation_position_frequency.csv", position_rows)

# Build a paired-seed comparison from EXISTING candidate records.
idx = {(r["arm"], r["seed"]): r for r in seed_summaries}
paired = []
for seed in PAIRED_SEEDS:
    a, b = idx.get(("adaptive", seed)), idx.get(("baseline", seed))
    if not a or not b:
        continue
    paired.append({
        "seed": seed,
        "adaptive_mean_native_score": a["mean_native_score"],
        "baseline_mean_native_score": b["mean_native_score"],
        "adaptive_minus_baseline_mean_native_score": (
            a["mean_native_score"] - b["mean_native_score"]
            if a["mean_native_score"] is not None and b["mean_native_score"] is not None else None
        ),
        "adaptive_best_native_score": a["best_native_score"],
        "baseline_best_native_score": b["best_native_score"],
        "adaptive_minus_baseline_best_native_score": (
            a["best_native_score"] - b["best_native_score"]
            if a["best_native_score"] is not None and b["best_native_score"] is not None else None
        ),
        "adaptive_unique_sequences": a["unique_sequences"],
        "baseline_unique_sequences": b["unique_sequences"],
    })
write_csv(OUT / "paired_existing_seed_comparison.csv", paired)

# Shared historic metrics are context, not relabeled as arm-specific evidence.
shared_metrics = []
for path in [
    RESULTS / "native_score_all_mutation_pairs" / "metrics.json",
    RESULTS / "p3_21aa_affinity_audit" / "native_c_score_21aa" / "metrics.json",
    RESULTS / "arise_experimental_position_audit" / "metrics.json",
]:
    if path.exists():
        try:
            shared_metrics.append((str(path.relative_to(ROOT)), json.loads(path.read_text())))
        except (OSError, json.JSONDecodeError):
            pass

with (OUT / "shared_historical_metrics.json").open("w") as f:
    json.dump(
        {path: obj for path, obj in shared_metrics},
        f, indent=2, ensure_ascii=False
    )

# Report.
lines = [
    "# P3 direct adaptive-versus-baseline comparison using existing candidates",
    "",
    f"Generated: {STAMP}",
    f"Existing paired candidate directory: `{RUNS_ROOT.relative_to(ROOT)}`",
    f"Paired seeds used: {', '.join(PAIRED_SEEDS)}",
    "",
    "## What this run did",
    "",
    "- Reused stored candidate histories from both arms.",
    "- Compared native sequence-fitness values, sequence diversity, and mutation-position frequencies by arm.",
    "- Searched existing measured mutation-pair tables, PEPBI tables, and Ridge out-of-fold rows for exact candidate links.",
    "- Did not launch an optimizer, create new candidates, or run docking.",
    "",
    "## Direct armwise candidate comparison",
    "",
    "| Arm | Records | Unique sequences | Diversity | Mean native score | Best native score |",
    "|---|---:|---:|---:|---:|---:|",
]
for s in arm_summaries:
    diversity = s["sequence_diversity_unique_over_records"]
    lines.append(
        f'| {s["arm"]} | {s["candidate_records"]} | {s["unique_sequences"]} | '
        f'{diversity:.3f} | {s["mean_native_score"]:.6f} | {s["best_native_score"]:.6f} |'
        if s["mean_native_score"] is not None and s["best_native_score"] is not None
        else f'| {s["arm"]} | {s["candidate_records"]} | {s["unique_sequences"]} | — | — | — |'
    )
lines += [
    "",
    "Paired per-seed mean/best score differences are in `paired_existing_seed_comparison.csv`.",
    "These are comparisons of stored native sequence-fitness values, not measured binding affinities.",
    "",
    "## Tests 3 and 4 — experimental mutation-affinity linkage",
    "",
    "| Arm | Exact candidate parent-child pairs found in experimental mutation tables |",
    "|---|---:|",
]
for s in arm_summaries:
    lines.append(f'| {s["arm"]} | {s["candidate_parent_child_pairs_with_exact_external_measurements"]} |')
lines += [
    "",
    "Full matches, including assay and target-context fields, are in "
    "`exact_candidate_parent_child_experimental_matches.csv`.",
    "A match by peptide sequence alone does not prove the assay used the project's 9SYA target. "
    "An armwise affinity comparison is biologically interpretable only when sequence, measurement, "
    "and relevant target/assay context are verified. If there are no valid records, the armwise "
    "experimental outcome is **not estimable**; the shared historical score benchmarks remain supporting "
    "evidence only.",
    "",
    "## Test 5 — mutation positions",
    "",
    "Arm-specific position frequencies are in `armwise_mutation_position_frequency.csv`. "
    "They directly compare where the existing adaptive and baseline candidates mutate. This is distinct "
    "from proving that those positions improve experimentally measured affinity; that requires valid "
    "candidate-linked measurements in the correct experimental context.",
    "",
    "## Test 6 — Ridge versus training-fold mean",
    "",
    "| Arm | Candidate sequences found in existing Ridge out-of-fold rows |",
    "|---|---:|",
]
for s in arm_summaries:
    lines.append(f'| {s["arm"]} | {s["candidate_sequences_found_in_Ridge_OOF_rows"]} |')
lines += [
    "",
    "Exact row links are in `ridge_oof_candidate_sequence_matches.csv`. The overall grouped Ridge-versus-mean "
    "benchmark evaluates prediction methods on the historical dataset; it is not itself a comparison of the "
    "two optimizers. It becomes candidate-linked only for candidate sequences with valid measured labels.",
    "",
    "## Test 7 — PEPBI",
    "",
    "| Arm | Unique exact candidate sequences found in PEPBI tables |",
    "|---|---:|",
]
for s in arm_summaries:
    lines.append(f'| {s["arm"]} | {s["exact_PEPBI_candidate_sequences"]} |')
lines += [
    "",
    "Record-level matches are in `pepbi_exact_sequence_matches_by_arm.csv`. An exact sequence match is not "
    "a binding-success result unless its experimental affinity and target/assay context are verified. "
    "Zero matches means outcome accuracy is unavailable, not zero percent.",
    "",
    "## Test 9 — receptor-context control",
    "",
    "The native C-score is sequence-only and does not encode receptor-specific interactions; applying it to "
    "both candidate arms therefore cannot establish a receptor-dependent binding advantage. This limitation "
    "applies to adaptive and baseline alike. A genuinely receptor-aware comparison would require a separate "
    "receptor-conditioned score evaluated on both arms with the same protocol; it must not be relabeled as "
    "native C-score.",
    "",
    "## Bottom line",
    "",
    "This run produces arm-by-arm candidate comparisons without repeating optimization. It also attempts "
    "explicit candidate-to-experiment links instead of silently treating shared historical benchmarks as "
    "optimizer results. If measured candidate-linked affinity records are absent or fail target-context "
    "verification, the experiment cannot yet support a claim that the adaptive arm produces better real "
    "binders than the baseline. The internal score and diversity comparison remains direct, but it is not "
    "a measured-affinity result.",
    "",
    "## Output files",
    "",
]
for p in sorted(OUT.iterdir()):
    lines.append(f"- `{p.name}`")
(OUT / "P3_ARMWISE_DIRECT_COMPARISON.md").write_text("\n".join(lines) + "\n")

print("Armwise analysis complete. No optimizer or docking was run.")
print(f"Candidate source: {RUNS_ROOT.relative_to(ROOT)}")
print(f"Report: {OUT.relative_to(ROOT) / 'P3_ARMWISE_DIRECT_COMPARISON.md'}")
print(f"Output directory: {OUT.relative_to(ROOT)}")
for s in arm_summaries:
    print(
        f'{s["arm"]}: records={s["candidate_records"]}, '
        f'unique={s["unique_sequences"]}, mean_native_score={s["mean_native_score"]}, '
        f'exact_pair_hits={s["candidate_parent_child_pairs_with_exact_external_measurements"]}, '
        f'PEPBI_hits={s["exact_PEPBI_candidate_sequences"]}, '
        f'OOF_hits={s["candidate_sequences_found_in_Ridge_OOF_rows"]}'
    )
