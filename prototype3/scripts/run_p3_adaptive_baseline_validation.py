from __future__ import annotations

import csv
import json
import math
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path.cwd()
RUNROOT = Path(sys.argv[1]).resolve()
OUT = RUNROOT / "p3_validation_comparison"
OUT.mkdir(parents=True, exist_ok=True)
AA = set("ACDEFGHIKLMNPQRSTVWY")
SEEDS = (2601, 2602, 2603)
ARMS = ("adaptive", "baseline")
rng = np.random.default_rng(20261027)


def normseq(x):
    if not isinstance(x, str):
        return None
    s = re.sub(r"\s+", "", x).upper()
    return s if 7 <= len(s) <= 60 and set(s) <= AA else None


def readjson(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return {}


def walk(obj):
    if isinstance(obj, dict):
        yield obj
        for v in obj.values():
            yield from walk(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk(v)


def direct_value(d, names):
    names = {n.lower() for n in names}
    return next((v for k, v in d.items() if k.lower() in names), None)


def all_value(obj, names):
    for d in walk(obj):
        v = direct_value(d, names)
        if v is not None and not isinstance(v, (dict, list)):
            return v
    return None


def load_candidates(folder):
    obj = readjson(folder / "candidate_history.json")
    records = []

    for d in walk(obj):
        seq = normseq(direct_value(d, {
            "sequence", "candidate_sequence", "peptide_sequence", "seq"
        }))
        score = direct_value(d, {
            "native_fitness", "c_score", "fitness_score", "native_score"
        })
        if not seq:
            continue
        try:
            score = float(score)
            if not math.isfinite(score):
                continue
        except (TypeError, ValueError):
            continue

        parent = normseq(direct_value(d, {"parent_sequence"}))
        gen = direct_value(d, {"generation", "generation_index", "gen"})
        try:
            gen = int(gen)
        except (TypeError, ValueError):
            gen = None

        records.append({
            "sequence": seq, "native_fitness": score,
            "parent_sequence": parent, "generation": gen
        })

    # If history records were nested in several ways, deduplicate exact records.
    unique = {}
    for r in records:
        unique[(r["sequence"], r["native_fitness"],
                r["parent_sequence"], r["generation"])] = r

    records = list(unique.values())
    if records:
        return records

    # Fallback to the run's top-candidate export.
    p = folder / "top_candidates.csv"
    if p.exists():
        for d in pd.read_csv(p).to_dict("records"):
            seq = normseq(direct_value(d, {
                "sequence", "candidate_sequence", "peptide_sequence", "seq"
            }))
            score = direct_value(d, {
                "native_fitness", "c_score", "fitness_score", "native_score"
            })
            try:
                score = float(score)
                if seq and math.isfinite(score):
                    records.append({
                        "sequence": seq, "native_fitness": score,
                        "parent_sequence": normseq(direct_value(
                            d, {"parent_sequence"}
                        )),
                        "generation": None,
                    })
            except (TypeError, ValueError):
                pass
    return records


def csvout(path, rows):
    if not rows:
        path.write_text("status\nno_eligible_records\n")
        return
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def ncol(df, aliases):
    def clean(x):
        return re.sub(r"[^a-z0-9]", "", str(x).lower())
    cols = {clean(c): c for c in df.columns}
    for a in aliases:
        if clean(a) in cols:
            return cols[clean(a)]
    return None


# 1. Load all six optimizer runs and their actual mutation events.
runs = {}
summaries = []
events = defaultdict(list)
selected = {}
hotspot_provenance = []

for arm in ARMS:
    for seed in SEEDS:
        folder = (
            RUNROOT / "runs" / arm / f"seed{seed}" / f"p3_{arm}_{seed}"
        )
        records = load_candidates(folder)
        runs[(arm, seed)] = {
            "folder": folder,
            "records": records,
            "seqs": {r["sequence"] for r in records},
        }

        for r in records:
            p, s = r["parent_sequence"], r["sequence"]
            if p and len(p) == len(s):
                events[(arm, seed)].extend(
                    i + 1 for i, (a, b) in enumerate(zip(p, s)) if a != b
                )

        scores = [r["native_fitness"] for r in records]
        gens = [r["generation"] for r in records if r["generation"] is not None]
        final_gen = max(gens) if gens else None
        final_scores = [
            r["native_fitness"] for r in records
            if final_gen is not None and r["generation"] == final_gen
        ]
        summaries.append({
            "arm": arm, "seed": seed,
            "candidate_records_read": len(records),
            "unique_sequences": len({r["sequence"] for r in records}),
            "mean_native_fitness": float(np.mean(scores)) if scores else None,
            "best_native_fitness": max(scores) if scores else None,
            "final_generation": final_gen,
            "final_generation_mean_fitness":
                float(np.mean(final_scores)) if final_scores else None,
            "parent_child_mutation_events": len(events[(arm, seed)]),
        })

        if arm != "adaptive":
            continue

        # Recover logged adaptive-selected positions. Do NOT substitute
        # frequently mutated positions and call them selected ARISE hotspots.
        for source_name, source_path in (
            ("experiment_audit", folder / "experiment_audit.json"),
            ("candidate_history", folder / "candidate_history.json"),
        ):
            obj = readjson(source_path)
            hits = []

            for d in walk(obj):
                for key, value in d.items():
                    k = key.lower()
                    if not (
                        ("hotspot" in k or "selected_positions" in k)
                        and "position" in k
                    ):
                        continue
                    if not isinstance(value, (list, tuple, set)):
                        continue
                    try:
                        vals = sorted({int(x) for x in value})
                    except (ValueError, TypeError):
                        continue
                    if not vals or len(vals) > 21:
                        continue
                    zero_based = (
                        "zero_based" in k or "0based" in k
                        or "index" in k or 0 in vals
                    )
                    pos = sorted({
                        v + 1 if zero_based else v
                        for v in vals
                        if 0 <= v <= 21
                    })
                    pos = [p for p in pos if 1 <= p <= 21]
                    if pos:
                        hits.append((key, pos, zero_based))

            if hits:
                counts = Counter(tuple(x[1]) for x in hits)
                win = list(counts.most_common(1)[0][0])
                key, _, zero_based = next(
                    h for h in reversed(hits) if h[1] == win
                )
                selected[seed] = win
                hotspot_provenance.append({
                    "seed": seed, "status": "logged positions found",
                    "source": source_name, "field": key,
                    "positions_1based": ",".join(map(str, win)),
                    "coordinate_conversion": (
                        "zero-based converted to one-based"
                        if zero_based else "treated as one-based"
                    ),
                })
                break

for seed in SEEDS:
    if seed not in selected:
        hotspot_provenance.append({
            "seed": seed,
            "status": "no explicit selected hotspot recovered",
            "source": "", "field": "",
            "positions_1based": "",
            "coordinate_conversion": "",
        })

csvout(OUT / "optimizer_runs.csv", summaries)
csvout(OUT / "adaptive_hotspot_provenance.csv", hotspot_provenance)

# 2. Load the already scored experimental single-mutation benchmark.
pair_path = ROOT / "prototype3/results/sequence_only_mutation_baseline/scored_mutation_pairs.csv"
pairrows = []
if pair_path.exists():
    df = pd.read_csv(pair_path)
    ca = ncol(df, ["sequence_a", "peptide_a", "seq_a"])
    cb = ncol(df, ["sequence_b", "peptide_b", "seq_b"])
    dp = ncol(df, ["delta_pKd", "delta_pKd_b_minus_a"])
    ds = ncol(df, ["delta_score", "delta_native_score", "delta_c_score"])
    pa = ncol(df, ["pKd_a", "pkd_a"])
    pb = ncol(df, ["pKd_b", "pkd_b"])
    poscol = ncol(df, ["mutation_position_1based", "mutation_position", "position"])
    contextcol = ncol(df, [
        "context_key", "target_group_id", "target_sequence",
        "target_id", "assay_id", "publication_id"
    ])

    if ca and cb:
        for d in df.to_dict("records"):
            a, b = normseq(str(d[ca])), normseq(str(d[cb]))
            if not a or not b or len(a) != len(b):
                continue
            changes = [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]
            if len(changes) != 1:
                continue

            try:
                exp = float(d[dp]) if dp else float(d[pb]) - float(d[pa])
                if not math.isfinite(exp):
                    continue
            except (TypeError, ValueError):
                continue

            try:
                score_delta = float(d[ds]) if ds else None
                if score_delta is not None and not math.isfinite(score_delta):
                    score_delta = None
            except (TypeError, ValueError):
                score_delta = None

            pairrows.append({
                "sequence_a": a, "sequence_b": b,
                "position": changes[0], "delta_pKd": exp,
                "abs_delta_pKd": abs(exp), "delta_score": score_delta,
                "context": str(d[contextcol]) if contextcol else "all",
                "pKd_a": float(d[pa]) if pa and pd.notna(d[pa]) else None,
                "pKd_b": float(d[pb]) if pb and pd.notna(d[pb]) else None,
            })

csvout(OUT / "experimental_pairs_reused.csv", pairrows)

# 3. Validate each adaptive-selected position set on the same external data.
#    The exact same selected positions are evaluated for adaptive and baseline.
posresults = []
for seed in SEEDS:
    positions = selected.get(seed)
    if not positions:
        posresults.append({
            "seed": seed, "arm": "both_same_test",
            "status": "blocked: no explicit adaptive-selected positions logged",
        })
        continue

    support = sorted({r["position"] for r in pairrows})
    used = sorted(set(positions) & set(support))
    data = [r for r in pairrows if r["position"] in support]

    if not used or len(used) == len(support) or not data:
        posresults.append({
            "seed": seed, "arm": "both_same_test",
            "status": "blocked: insufficient supported positional contrast",
            "positions": ",".join(map(str, positions)),
        })
        continue

    contexts = defaultdict(list)
    for r in data:
        contexts[r["context"]].append(r)
    centered = {}
    for ctx, rows in contexts.items():
        mean = float(np.mean([r["abs_delta_pKd"] for r in rows]))
        for r in rows:
            centered[id(r)] = r["abs_delta_pKd"] - mean

    inside = [r for r in data if r["position"] in used]
    outside = [r for r in data if r["position"] not in used]
    effect = float(np.mean([centered[id(r)] for r in inside]) -
                    np.mean([centered[id(r)] for r in outside]))

    perm = []
    for _ in range(3000):
        draw = set(rng.choice(support, size=len(used), replace=False).tolist())
        inn = [r for r in data if r["position"] in draw]
        out = [r for r in data if r["position"] not in draw]
        if inn and out:
            perm.append(float(
                np.mean([centered[id(r)] for r in inn]) -
                np.mean([centered[id(r)] for r in out])
            ))
    p = (1 + sum(abs(x) >= abs(effect) for x in perm)) / (1 + len(perm))

    for arm in ARMS:
        ev = events[(arm, seed)]
        posresults.append({
            "seed": seed, "arm": arm,
            "status": "computed; same frozen adaptive-selected positions",
            "positions": ",".join(map(str, used)),
            "pairs_inside": len(inside), "pairs_outside": len(outside),
            "context_centered_effect": effect,
            "random_position_permutation_p_two_sided": p,
            "arm_mutation_events": len(ev),
            "fraction_arm_mutation_events_inside_selected_positions":
                float(np.mean([v in used for v in ev])) if ev else None,
            "note": (
                "External window statistic is identical for both arms by design; "
                "only each arm's mutation-event hit rate differs."
            ),
        })
csvout(OUT / "mutation_position_validation.csv", posresults)

# 4. Native score versus measured affinity, reported separately by arm.
#    A matched candidate-pKd association is computed only for exact sequence
#    matches. A second exploratory test uses external measured mutation pairs
#    at the positions each arm actually mutated.
known_pkd = defaultdict(list)
for r in pairrows:
    if r["pKd_a"] is not None:
        known_pkd[r["sequence_a"]].append(r["pKd_a"])
    if r["pKd_b"] is not None:
        known_pkd[r["sequence_b"]].append(r["pKd_b"])

native_results = []
candidate_measurements = []
for arm in ARMS:
    seqscores = defaultdict(list)
    observed_positions = set()
    for seed in SEEDS:
        for r in runs[(arm, seed)]["records"]:
            seqscores[r["sequence"]].append(r["native_fitness"])
        observed_positions.update(events[(arm, seed)])

    for seq, scores in seqscores.items():
        vals = known_pkd.get(seq, [])
        if vals:
            candidate_measurements.append({
                "arm": arm, "sequence": seq,
                "mean_native_fitness": float(np.mean(scores)),
                "median_measured_pKd": float(np.median(vals)),
                "measurement_rows": len(vals),
            })

    matched = [r for r in candidate_measurements if r["arm"] == arm]
    rho_exact = None
    if len(matched) >= 3:
        x = [r["mean_native_fitness"] for r in matched]
        y = [r["median_measured_pKd"] for r in matched]
        if len(set(x)) > 1 and len(set(y)) > 1:
            rho_exact = float(spearmanr(x, y).statistic)

    subset = [r for r in pairrows if r["position"] in observed_positions]
    scored = [r for r in subset if r["delta_score"] is not None]
    rho_pos = None
    direction = None
    if len(scored) >= 3 and len({r["delta_score"] for r in scored}) > 1:
        rho_pos = float(spearmanr(
            [r["delta_score"] for r in scored],
            [r["delta_pKd"] for r in scored]
        ).statistic)
        informative = [r for r in scored if abs(r["delta_pKd"]) >= 0.3]
        if informative:
            direction = float(np.mean([
                np.sign(r["delta_score"]) == np.sign(r["delta_pKd"])
                for r in informative
            ]))

    native_results.append({
        "arm": arm,
        "unique_generated_sequences": len(seqscores),
        "exact_generated_sequences_with_measured_pKd": len(matched),
        "exact_candidate_native_score_vs_measured_pKd_spearman": rho_exact,
        "external_pairs_at_arm_mutated_positions": len(subset),
        "scored_external_pairs_at_arm_mutated_positions": len(scored),
        "external_delta_score_vs_delta_pKd_spearman": rho_pos,
        "external_informative_direction_accuracy": direction,
        "scope": (
            "Exact sequence matches measure generated-candidate association. "
            "Position-filtered external pairs are retrospective and exploratory, "
            "not measurements of the generated candidates."
        ),
    })

csvout(OUT / "native_score_affinity_by_arm.csv", native_results)
csvout(OUT / "exact_candidate_affinity_matches.csv", candidate_measurements)

# 5. PEPBI exact candidate-sequence overlap. Do not invent measurements.
pepbi_rows = []
workbook = ROOT / "prototype3/data/raw/PEPBI.xlsx"
pepbi_sequences = defaultdict(list)

if workbook.exists():
    try:
        sheets = pd.read_excel(workbook, sheet_name=None)
        candidate_seqs = set()
        for arm in ARMS:
            for seed in SEEDS:
                candidate_seqs.update(runs[(arm, seed)]["seqs"])

        for sheet, table in sheets.items():
            for col in table.columns:
                for idx, value in enumerate(table[col].tolist()):
                    seq = normseq(value)
                    if seq and seq in candidate_seqs:
                        pepbi_sequences[seq].append({
                            "sheet": sheet, "row": idx + 2,
                            "matched_column": str(col),
                        })

        for arm in ARMS:
            for seed in SEEDS:
                for seq in sorted(runs[(arm, seed)]["seqs"] & set(pepbi_sequences)):
                    for src in pepbi_sequences[seq]:
                        pepbi_rows.append({
                            "arm": arm, "seed": seed, "sequence": seq, **src
                        })
    except Exception as e:
        pepbi_rows.append({"status": f"workbook read failed: {e}"})
else:
    pepbi_rows.append({"status": f"workbook missing: {workbook}"})

csvout(OUT / "pepbi_exact_candidate_overlaps.csv", pepbi_rows)
pepbi_summary = []
for arm in ARMS:
    seqs = set().union(*(runs[(arm, seed)]["seqs"] for seed in SEEDS))
    matches = seqs & set(pepbi_sequences)
    pepbi_summary.append({
        "arm": arm, "unique_candidate_sequences": len(seqs),
        "exact_pepbi_sequence_matches": len(matches),
        "measured_affinity_validation": (
            "possible only for exact matches with a valid experimental outcome"
            if matches else "not estimable: no exact PEPBI candidate overlap"
        ),
    })
csvout(OUT / "pepbi_summary.csv", pepbi_summary)

# 6. New structural scoring only for actual complex PDBs in these six run folders.
AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}


def pdbchains(path):
    chains = defaultdict(dict)
    with Path(path).open(errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM  ") or len(line) < 27:
                continue
            if line[12:16].strip() != "CA":
                continue
            name = line[17:20].strip().upper()
            if name not in AA3:
                continue
            chain = line[21].strip() or "_"
            try:
                rn = int(line[22:26])
            except ValueError:
                continue
            chains[chain][(rn, line[26].strip())] = name
    return chains


def seqof(residues):
    return "".join(AA3[v] for _, v in sorted(residues.items()))


struct_inventory, struct_attempts = [], []
scorer = ROOT / "prototype3/run_target_aware_pilot.py"

for arm in ARMS:
    for seed in SEEDS:
        folder = runs[(arm, seed)]["folder"]
        candidates = runs[(arm, seed)]["seqs"]
        seen = 0
        for pdb in folder.rglob("*.pdb"):
            try:
                chains = pdbchains(pdb)
            except OSError:
                continue
            peptide_chains = [
                c for c, residues in chains.items()
                if seqof(residues) in candidates
            ]
            receptor_chains = [
                c for c, residues in chains.items()
                if len(residues) > 100
                and residues.get((203, ""), "") == "ARG"
                and residues.get((205, ""), "") == "ARG"
            ]
            eligible = bool(peptide_chains and receptor_chains)
            struct_inventory.append({
                "arm": arm, "seed": seed, "complex_pdb": str(pdb),
                "candidate_peptide_chains": ",".join(peptide_chains),
                "receptor_chains_with_Arg203_Arg205": ",".join(receptor_chains),
                "eligible": eligible,
            })
            if not eligible or seen >= 5:
                continue

            for pc in peptide_chains:
                for rc in receptor_chains:
                    dest = OUT / "structural_scores" / f"{arm}_{seed}_{pdb.stem}_{rc}_{pc}"
                    cmd = [
                        sys.executable, str(scorer),
                        "--complex-pdb", str(pdb),
                        "--receptor-chain", rc, "--peptide-chain", pc,
                        "--target-id", f"p3_{arm}_{seed}",
                        "--hotspot", f"{rc}:203:ARG",
                        "--hotspot", f"{rc}:205:ARG",
                        "--out-root", str(dest), "--skip-c-score",
                    ]
                    try:
                        result = subprocess.run(
                            cmd, cwd=ROOT, capture_output=True, text=True,
                            timeout=180
                        )
                        metrics = list(dest.rglob("benchmark_metrics.json"))
                        struct_attempts.append({
                            "arm": arm, "seed": seed, "complex_pdb": str(pdb),
                            "exit_code": result.returncode,
                            "metrics_file": str(metrics[-1]) if metrics else "",
                            "stderr_tail": result.stderr[-700:],
                        })
                    except Exception as e:
                        struct_attempts.append({
                            "arm": arm, "seed": seed, "complex_pdb": str(pdb),
                            "error": str(e),
                        })
                    seen += 1

csvout(OUT / "structural_complex_inventory.csv", struct_inventory)
csvout(OUT / "new_structural_scoring_attempts.csv", struct_attempts)

# 7. Locate the fresh strict grouped Ridge versus training-fold-mean run.
grouped_files = sorted(
    (RUNROOT / "grouped_sequence_validation").glob("run_*/metrics.json")
)
grouped_path = grouped_files[-1] if grouped_files else None
grouped_metrics = readjson(grouped_path) if grouped_path else {}

report = [
    "# P3 adaptive-versus-baseline validation",
    "",
    f"Optimizer runs: `{RUNROOT}`",
    f"Report directory: `{OUT}`",
    "",
    "## 1. Optimizer outputs",
    "",
    "| Arm | Seed | Candidates read | Unique sequences | Mean native fitness | Best | Final-generation mean |",
    "|---|---:|---:|---:|---:|---:|---:|",
]
for r in summaries:
    report.append(
        f"| {r['arm']} | {r['seed']} | {r['candidate_records_read']} | "
        f"{r['unique_sequences']} | {r['mean_native_fitness']} | "
        f"{r['best_native_fitness']} | {r['final_generation_mean_fitness']} |"
    )

report += [
    "",
    "## 2. Adaptive-selected mutation positions",
    "",
    "Source: `adaptive_hotspot_provenance.csv`. Position-window tests use the same "
    "adaptive-selected positions for the matched baseline. If positions were not logged, "
    "the test is marked blocked; observed mutation frequency is not relabeled as ARISE output.",
    "",
]
for r in hotspot_provenance:
    report.append(
        f"- Seed {r['seed']}: {r['status']}; positions={r.get('positions_1based', '')}; "
        f"source={r.get('source', '')}."
    )

report += [
    "",
    "Mutation statistics: `mutation_position_validation.csv`.",
    "",
    "## 3. Native score versus experimental affinity",
    "",
    "| Arm | Exact candidate/pKd matches | Exact-match rho | External pairs at observed mutation positions | External rho | External direction accuracy |",
    "|---|---:|---:|---:|---:|---:|",
]
for r in native_results:
    report.append(
        f"| {r['arm']} | {r['exact_generated_sequences_with_measured_pKd']} | "
        f"{r['exact_candidate_native_score_vs_measured_pKd_spearman']} | "
        f"{r['external_pairs_at_arm_mutated_positions']} | "
        f"{r['external_delta_score_vs_delta_pKd_spearman']} | "
        f"{r['external_informative_direction_accuracy']} |"
    )
report += [
    "",
    "Exact candidate matches are direct but may be too sparse to estimate accuracy. "
    "The external pair analysis filters a retrospective benchmark by positions observed "
    "in each arm; it is not experimental measurement of newly generated candidates.",
    "",
    "## 4. PEPBI",
    "",
]
for r in pepbi_summary:
    report.append(
        f"- {r['arm']}: {r['exact_pepbi_sequence_matches']} exact candidate-sequence "
        f"matches out of {r['unique_candidate_sequences']} unique sequences."
    )
report += [
    "",
    "Detailed matches: `pepbi_exact_candidate_overlaps.csv`. An overlap count alone is "
    "not evidence of predictive accuracy; measured outcomes must be available for the same sequence.",
    "",
    "## 5. New target-aware structural scoring",
    "",
    f"- Candidate-run PDBs inventoried: {len(struct_inventory)}",
    f"- Structural scoring attempts: {len(struct_attempts)}",
    f"- Eligible attempt metric files: {sum(bool(r.get('metrics_file')) and r.get('exit_code') == 0 for r in struct_attempts)}",
    "",
    "Files: `structural_complex_inventory.csv`, `new_structural_scoring_attempts.csv`. "
    "A generated candidate without a valid receptor–peptide complex pose cannot be scored structurally.",
    "",
    "## 6. Grouped Ridge versus training-fold mean",
    "",
    f"- Fresh metrics: `{grouped_path}`" if grouped_path else "- Fresh grouped metrics not found.",
    "- This is a separate sequence-only prediction benchmark, not a direct optimizer-arm comparison.",
    "",
    "## Interpretation",
    "",
    "These analyses distinguish internal optimization, retrospective position associations, "
    "exact-sequence affinity evidence, PEPBI overlap, and structural pose evidence. "
    "No experimental affinity is imputed for an unmeasured candidate. A higher native fitness "
    "score must not be described as stronger real-world binding without experimental validation.",
    "",
]

report_path = OUT / "P3_ADAPTIVE_BASELINE_VALIDATION.md"
report_path.write_text("\n".join(report), encoding="utf-8")

print("\n".join(report))
print(f"\nSaved final report: {report_path}")
print(f"CSV evidence tables: {OUT}")
