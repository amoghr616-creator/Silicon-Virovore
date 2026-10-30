import csv
import json
import math
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path.cwd()
RESULTS = ROOT / "prototype3" / "results"
RUNS = RESULTS / "p3_paired_arm_20261027T022810Z" / "runs"

SOURCES = [
    (
        "PPIKB_KD_9_20",
        RESULTS / "ppikb_kd_benchmark" / "filtered_kd_measurements.csv",
        (9, 20),
    ),
    (
        "PPIKB_KD_EXACT21",
        RESULTS / "p3_21aa_affinity_audit" / "native_c_score_21aa"
        / "eligible_scored_measurements.csv",
        (21, 21),
    ),
]

STAMP = datetime.now().strftime("%Y%m%dT%H%M%S")
OUT = RESULTS / f"p3_candidate_ppikb_join_{STAMP}"
OUT.mkdir(parents=True, exist_ok=False)

AA = set("ACDEFGHIKLMNPQRSTVWY")


def norm(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def col(row, aliases):
    lookup = {norm(k): k for k in row}
    for name in aliases:
        if norm(name) in lookup:
            return lookup[norm(name)]
    return None


def get(row, aliases, default=""):
    key = col(row, aliases)
    return row.get(key, default) if key else default


def number(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def yes(x):
    return str(x).strip().lower() in {"true", "1", "yes"}


def no(x):
    return str(x).strip().lower() in {"false", "0", "no"}


def sequence(x):
    s = re.sub(r"\s+", "", str(x or "").upper())
    return s if s and set(s) <= AA else ""


def rows_from(path):
    if not path.exists():
        return None
    with path.open(newline="", errors="replace") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    return rows


def affinity_pkd(row):
    # Prefer a recorded pKd / p-affinity measurement.
    for alias in ("_p_affinity", "pKd", "p_affinity",
                  "measured_pKd", "pkd"):
        value = number(get(row, [alias], None))
        if value is not None:
            return value

    # Otherwise derive pKd from a measured Kd in molar units.
    molar = number(get(row, ["_molar", "Kd_molar", "molar"], None))
    if molar is not None and molar > 0:
        return -math.log10(molar)

    value = number(get(row, ["_value", "Kd_value"], None))
    unit = str(get(row, ["_unit", "Kd_unit"], "")).strip().lower()
    factors = {
        "m": 1.0, "mol": 1.0, "molar": 1.0,
        "mm": 1e-3, "millimolar": 1e-3,
        "um": 1e-6, "μm": 1e-6, "micromolar": 1e-6,
        "nm": 1e-9, "nanomolar": 1e-9,
        "pm": 1e-12, "picomolar": 1e-12,
    }
    if value is not None and value > 0 and unit in factors:
        return -math.log10(value * factors[unit])
    return None


def read_candidates():
    if not RUNS.exists():
        raise SystemExit(f"Saved candidate histories not found: {RUNS}")

    records = []
    for arm in ("adaptive", "baseline"):
        paths = sorted((RUNS / arm).rglob("candidate_history.json"))
        if not paths:
            raise SystemExit(f"No stored candidate histories for {arm}: {RUNS / arm}")

        for path in paths:
            try:
                history = json.loads(path.read_text())
            except Exception as exc:
                print("WARNING: skipped unreadable history:", path, repr(exc))
                continue

            seed_match = re.search(r"seed(\d+)", str(path))
            seed = seed_match.group(1) if seed_match else "unknown"

            if not isinstance(history, list):
                continue

            for item in history:
                if not isinstance(item, dict):
                    continue
                seq = sequence(item.get("sequence"))
                if not seq:
                    continue

                records.append({
                    "arm": arm,
                    "seed": seed,
                    "candidate_id": str(item.get("candidate_id", "")),
                    "sequence": seq,
                    "length": len(seq),
                    "parent_sequence": sequence(item.get("parent_sequence")),
                    "history_file": str(path.relative_to(ROOT)),
                })

    return records


def read_measurements(label, path, allowed_lengths):
    raw = rows_from(path)
    if raw is None:
        return [], {
            "source": str(path.relative_to(ROOT)),
            "status": "FILE_MISSING",
            "rows": 0,
        }

    result = []
    rejected = defaultdict(int)
    if not raw:
        return [], {
            "source": str(path.relative_to(ROOT)),
            "status": "EMPTY_FILE",
            "rows": 0,
        }

    for row_num, row in enumerate(raw, start=2):
        seq = sequence(get(row, [
            "_peptide", "Peptide_Sequence", "peptide_sequence",
            "canonical_peptide_sequence", "canonical_sequence",
            "peptide_seq", "sequence",
        ]))

        if not seq:
            rejected["invalid_or_missing_peptide_sequence"] += 1
            continue
        if not allowed_lengths[0] <= len(seq) <= allowed_lengths[1]:
            rejected["outside_expected_length_range"] += 1
            continue

        # Apply eligibility flags where they are available.
        valid_flags = [
            (("_eligible", "eligible"), yes),
            (("_target_valid", "target_valid"), yes),
            (("_canonical", "canonical"), yes),
            (("_length_ok", "length_ok"), yes),
            (("_linear", "linear"), yes),
            (("_unmodified", "unmodified"), yes),
        ]
        bad = False
        for aliases, predicate in valid_flags:
            key = col(row, aliases)
            if key and not predicate(row.get(key, "")):
                rejected["failed_" + aliases[0].lstrip("_")] += 1
                bad = True
                break
        if bad:
            continue

        censored_key = col(row, ("_censored", "censored"))
        if censored_key and not no(row.get(censored_key, "")):
            rejected["censored"] += 1
            continue

        type_value = str(get(row, [
            "_affinity_type", "affinity_type", "affinity_type_parsed"
        ])).strip().upper()
        if type_value and type_value != "KD":
            rejected["not_KD"] += 1
            continue

        pkd = affinity_pkd(row)
        if pkd is None:
            rejected["missing_usable_pKd"] += 1
            continue

        protein_seq = sequence(get(row, [
            "_protein", "Protein_Sequence", "protein_sequence",
            "target_sequence", "receptor_sequence",
        ]))
        protein_name = str(get(row, [
            "_name", "Protein_Name", "protein_name", "target_name"
        ])).strip()
        assay = str(get(row, [
            "_assay", "Affinity_Assay", "assay", "assay_name"
        ])).strip()
        uniprot = str(get(row, ["UniProt_ID", "uniprot_id"])).strip()
        pdb = str(get(row, ["PDB_ID", "pdb_id"])).strip()
        affinity_text = str(get(row, ["Affinity", "affinity"])).strip()
        mutations = str(get(row, [
            "protein_mutations", "target_mutations", "Protein_Mutations"
        ])).strip()
        site = str(get(row, ["binding_site", "Binding_Site"])).strip()
        temperature = str(get(row, [
            "temperature_K", "Temperature_K", "temperature"
        ])).strip()
        group_name = str(get(row, ["binding_group", "target_group"])).strip()

        # Prefer exact target sequence when supplied; retain assay and
        # any available variant/site/temperature fields.
        target_identity = protein_seq or protein_name or group_name
        context = json.dumps(
            [target_identity, assay, mutations, site, temperature],
            ensure_ascii=False,
        )

        target_text = " ".join([
            protein_name, group_name,
            str(get(row, ["Gene_Names", "gene_names", "gene"])),
            str(get(row, ["Title", "title"])),
        ]).lower()
        looks_hervk = bool(re.search(
            r"herv[\s-]*k|human endogenous retrovirus type k|"
            r"human endogenous retrovirus k",
            target_text,
        ))

        result.append({
            "dataset": label,
            "source_file": str(path.relative_to(ROOT)),
            "source_row": row_num,
            "source_id": str(get(row, ["ID", "id"])),
            "sequence": seq,
            "length": len(seq),
            "pKd": pkd,
            "protein_name": protein_name,
            "protein_sequence": protein_seq,
            "uniprot_id": uniprot,
            "assay": assay,
            "affinity_text": affinity_text,
            "PDB_ID": pdb,
            "mutations": mutations,
            "binding_site": site,
            "temperature_K": temperature,
            "target_context_key": context,
            "target_label_mentions_HERV_K": looks_hervk,
        })

    # Collapse repeated records for the same sequence and exact context
    # using median pKd; preserve record counts and source IDs.
    buckets = defaultdict(list)
    for row in result:
        buckets[(row["dataset"], row["target_context_key"], row["sequence"])].append(row)

    collapsed = []
    for (_, _, _), items in buckets.items():
        item = dict(items[0])
        values = sorted(x["pKd"] for x in items)
        n = len(values)
        item["pKd"] = (
            values[n // 2] if n % 2 else
            (values[n // 2 - 1] + values[n // 2]) / 2
        )
        item["source_rows_collapsed"] = n
        item["source_ids"] = ";".join(sorted({
            x["source_id"] for x in items if x["source_id"]
        }))
        collapsed.append(item)

    return collapsed, {
        "source": str(path.relative_to(ROOT)),
        "status": "loaded",
        "raw_rows": len(raw),
        "eligible_sequence_context_rows": len(collapsed),
        "distinct_sequences": len({r["sequence"] for r in collapsed}),
        "contexts": len({r["target_context_key"] for r in collapsed}),
        "rejected_rows": dict(rejected),
        "length_range": list(allowed_lengths),
    }


def write_csv(path, rows, fields):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    candidates = read_candidates()
    unique_arm_sequences = {
        arm: {r["sequence"] for r in candidates if r["arm"] == arm}
        for arm in ("adaptive", "baseline")
    }

    measurements, source_audit = [], []
    for label, path, lengths in SOURCES:
        loaded, audit = read_measurements(label, path, lengths)
        measurements.extend(loaded)
        source_audit.append(audit)

    # Match each unique generated sequence against every eligible
    # measured context, retaining source and target/assay metadata.
    candidate_arms = defaultdict(set)
    candidate_seeds = defaultdict(lambda: defaultdict(set))
    candidate_record_counts = defaultdict(lambda: defaultdict(int))
    for r in candidates:
        candidate_arms[r["sequence"]].add(r["arm"])
        candidate_seeds[r["sequence"]][r["arm"]].add(r["seed"])
        candidate_record_counts[r["sequence"]][r["arm"]] += 1

    measurements_by_sequence = defaultdict(list)
    for m in measurements:
        measurements_by_sequence[m["sequence"]].append(m)

    matched_rows = []
    for seq in sorted(candidate_arms):
        if seq not in measurements_by_sequence:
            continue

        for m in measurements_by_sequence[seq]:
            matched_rows.append({
                "candidate_sequence": seq,
                "candidate_length": len(seq),
                "arms_containing_sequence": ";".join(sorted(candidate_arms[seq])),
                "adaptive_seeds": ";".join(sorted(candidate_seeds[seq]["adaptive"])),
                "baseline_seeds": ";".join(sorted(candidate_seeds[seq]["baseline"])),
                "adaptive_candidate_records": candidate_record_counts[seq]["adaptive"],
                "baseline_candidate_records": candidate_record_counts[seq]["baseline"],
                "dataset": m["dataset"],
                "pKd": m["pKd"],
                "protein_name": m["protein_name"],
                "protein_sequence": m["protein_sequence"],
                "uniprot_id": m["uniprot_id"],
                "assay": m["assay"],
                "affinity_text": m["affinity_text"],
                "PDB_ID": m["PDB_ID"],
                "source_file": m["source_file"],
                "source_row": m["source_row"],
                "source_ids": m["source_ids"],
                "source_rows_collapsed": m["source_rows_collapsed"],
                "target_context_key": m["target_context_key"],
                "target_label_mentions_HERV_K": m["target_label_mentions_HERV_K"],
                "manual_target_and_assay_verification_required": "YES",
            })

    # Find exact target-assay contexts containing candidates unique to
    # both arms; do not count a peptide shared by both arms as evidence
    # that one arm did better.
    context_map = defaultdict(lambda: defaultdict(dict))
    meta_map = {}

    for match in matched_rows:
        seq = match["candidate_sequence"]
        context = match["target_context_key"]
        arms = candidate_arms[seq]
        meta_map[context] = match
        for arm in arms:
            context_map[context][arm][seq] = match["pKd"]

    paired_context_rows = []
    for context, arm_dict in context_map.items():
        adaptive = arm_dict.get("adaptive", {})
        baseline = arm_dict.get("baseline", {})
        # Only arm-exclusive sequences contribute to the arm contrast.
        a_only = sorted(set(adaptive) - set(baseline))
        b_only = sorted(set(baseline) - set(adaptive))
        if not a_only or not b_only:
            continue

        a_mean = sum(adaptive[s] for s in a_only) / len(a_only)
        b_mean = sum(baseline[s] for s in b_only) / len(b_only)
        meta = meta_map[context]
        paired_context_rows.append({
            "dataset": meta["dataset"],
            "protein_name": meta["protein_name"],
            "assay": meta["assay"],
            "context_key": context,
            "adaptive_unique_sequences_with_measurements": len(a_only),
            "baseline_unique_sequences_with_measurements": len(b_only),
            "adaptive_mean_pKd": a_mean,
            "baseline_mean_pKd": b_mean,
            "adaptive_minus_baseline_pKd": a_mean - b_mean,
            "target_label_mentions_HERV_K": meta["target_label_mentions_HERV_K"],
            "interpretation": "Descriptive matched-context comparison; manual biological context verification required.",
        })

    arm_summary = []
    for arm in ("adaptive", "baseline"):
        seqs = unique_arm_sequences[arm]
        matched = seqs & set(measurements_by_sequence)
        lengths21 = {s for s in matched if len(s) == 21}
        hervk = {
            s for s in matched
            if any(m["target_label_mentions_HERV_K"]
                   for m in measurements_by_sequence[s])
        }
        arm_summary.append({
            "arm": arm,
            "candidate_records": sum(r["arm"] == arm for r in candidates),
            "unique_candidate_sequences": len(seqs),
            "unique_sequences_with_any_eligible_exact_PPIKB_match": len(matched),
            "matched_unique_21aa_sequences": len(lengths21),
            "candidate_sequences_with_HERV_K_labelled_context": len(hervk),
            "matched_measurement_context_rows": sum(
                len(measurements_by_sequence[s]) for s in matched
            ),
            "same_context_armwise_comparison_estimable": any(
                r["adaptive_unique_sequences_with_measurements"] > 0
                and r["baseline_unique_sequences_with_measurements"] > 0
                for r in paired_context_rows
            ),
        })

    shared_seq = unique_arm_sequences["adaptive"] & unique_arm_sequences["baseline"]
    all_candidate_seq = unique_arm_sequences["adaptive"] | unique_arm_sequences["baseline"]
    summary = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "candidate_history_source": str(RUNS.relative_to(ROOT)),
        "candidate_record_count": len(candidates),
        "adaptive_unique_sequences": len(unique_arm_sequences["adaptive"]),
        "baseline_unique_sequences": len(unique_arm_sequences["baseline"]),
        "sequences_shared_by_both_arms": len(shared_seq),
        "all_unique_generated_sequences": len(all_candidate_seq),
        "measurement_source_audit": source_audit,
        "matched_candidate_measurement_rows": len(matched_rows),
        "candidate_sequences_with_any_exact_measurement_match": len({
            r["candidate_sequence"] for r in matched_rows
        }),
        "candidate_sequences_with_HERV_K_labelled_context": len({
            r["candidate_sequence"] for r in matched_rows
            if r["target_label_mentions_HERV_K"]
        }),
        "same_context_armwise_comparison_rows": len(paired_context_rows),
        "warning": (
            "Exact sequence matching does not prove HERV-K Env binding. "
            "Target, assay, measurement provenance and relevance must be verified. "
            "PPIKB has been examined previously in this project, so this is a "
            "retrospective data-linkage audit, not independent prospective validation."
        ),
    }

    candidate_fields = [
        "arm", "seed", "candidate_id", "sequence", "length",
        "parent_sequence", "history_file",
    ]
    match_fields = [
        "candidate_sequence", "candidate_length", "arms_containing_sequence",
        "adaptive_seeds", "baseline_seeds", "adaptive_candidate_records",
        "baseline_candidate_records", "dataset", "pKd", "protein_name",
        "protein_sequence", "uniprot_id", "assay", "affinity_text", "PDB_ID",
        "source_file", "source_row", "source_ids", "source_rows_collapsed",
        "target_context_key", "target_label_mentions_HERV_K",
        "manual_target_and_assay_verification_required",
    ]
    measurement_fields = [
        "dataset", "source_file", "source_row", "source_id", "sequence",
        "length", "pKd", "protein_name", "protein_sequence", "uniprot_id",
        "assay", "affinity_text", "PDB_ID", "mutations", "binding_site",
        "temperature_K", "target_context_key", "target_label_mentions_HERV_K",
        "source_rows_collapsed", "source_ids",
    ]

    write_csv(OUT / "saved_candidate_inventory.csv", candidates, candidate_fields)
    write_csv(OUT / "eligible_measurement_contexts.csv", measurements, measurement_fields)
    write_csv(OUT / "exact_candidate_measurement_matches.csv", matched_rows, match_fields)
    write_csv(
        OUT / "matched_same_context_arm_comparison.csv",
        paired_context_rows,
        [
            "dataset", "protein_name", "assay", "context_key",
            "adaptive_unique_sequences_with_measurements",
            "baseline_unique_sequences_with_measurements",
            "adaptive_mean_pKd", "baseline_mean_pKd",
            "adaptive_minus_baseline_pKd",
            "target_label_mentions_HERV_K", "interpretation",
        ],
    )
    write_csv(
        OUT / "armwise_join_summary.csv",
        arm_summary,
        [
            "arm", "candidate_records", "unique_candidate_sequences",
            "unique_sequences_with_any_eligible_exact_PPIKB_match",
            "matched_unique_21aa_sequences",
            "candidate_sequences_with_HERV_K_labelled_context",
            "matched_measurement_context_rows",
            "same_context_armwise_comparison_estimable",
        ],
    )
    write_csv(
        OUT / "measurement_source_audit.csv",
        source_audit,
        ["source", "status", "raw_rows", "eligible_sequence_context_rows",
         "distinct_sequences", "contexts", "length_range", "rejected_rows"],
    )
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    lines = [
        "# Saved-candidate PPIKB linkage audit",
        "",
        f"Generated: {summary['created']}",
        "",
        "## Purpose",
        "",
        "Determine whether any saved adaptive/baseline candidate sequence exactly matches an eligible PPIKB measurement, including the separately audited exact-21-aa records. No optimizer, docking, or scoring implementation was run or modified.",
        "",
        "## Data sources",
        "",
    ]
    for item in source_audit:
        lines.append(f"- `{item['source']}`: **{item['status']}**; rows loaded: {item.get('eligible_sequence_context_rows', 0)}.")
    lines += [
        "",
        "## Armwise linkage",
        "",
        "| Arm | Candidate records | Unique sequences | Exact measured sequence matches | Exact 21-aa matches | HERV-K-labelled contexts |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in arm_summary:
        lines.append(
            f"| {r['arm']} | {r['candidate_records']} | {r['unique_candidate_sequences']} | "
            f"{r['unique_sequences_with_any_eligible_exact_PPIKB_match']} | "
            f"{r['matched_unique_21aa_sequences']} | "
            f"{r['candidate_sequences_with_HERV_K_labelled_context']} |"
        )
    lines += [
        "",
        f"Sequence matches across both arms: {summary['sequences_shared_by_both_arms']}.",
        f"Matched same-context arm comparisons: {summary['same_context_armwise_comparison_rows']}.",
        "",
        "## Interpretation",
        "",
    ]
    if not matched_rows:
        lines += [
            "**No exact candidate sequence match was found in the eligible files searched.** The current PPIKB records therefore cannot supply measured affinity labels for the saved candidate lists in this audit. That does not mean the candidates have zero binding; their experimental outcomes are unmeasured here.",
            "",
        ]
    else:
        lines += [
            "Exact sequence matches were found. Inspect `exact_candidate_measurement_matches.csv` and verify the source record, target protein, assay, units, and experimental context before using a match as biological evidence.",
            "",
            "Matches against unrelated target proteins can support only a general peptide-data observation; they do not demonstrate HERV-K Env binding.",
            "",
        ]
    lines += [
        "An armwise measured-affinity comparison is estimable only in target-assay contexts that contain measured sequences unique to both arms. Any comparison across other targets remains descriptive and retrospective.",
        "",
        "## Conclusion",
        "",
        "This audit resolves whether existing eligible PPIKB files contain exact candidate links. It does not manufacture missing affinity measurements and it does not establish HERV-K-specific binding without verified target-relevant records.",
        "",
    ]
    (OUT / "PPIKB_LINKAGE_REPORT.md").write_text("\n".join(lines))

    print("\nCandidate-to-PPIKB audit completed.")
    print("Optimizer/docking: NOT RUN")
    print("Output:", OUT.relative_to(ROOT))
    print("\nArmwise summary:")
    for row in arm_summary:
        print(row)
    print("\nSource audit:")
    for row in source_audit:
        print(row)
    print("\nMatched candidate measurement rows:", len(matched_rows))
    print("Same-context armwise comparisons:", len(paired_context_rows))
    print("Report:", OUT.relative_to(ROOT) / "PPIKB_LINKAGE_REPORT.md")


if __name__ == "__main__":
    main()
