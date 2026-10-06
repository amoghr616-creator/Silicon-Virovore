import json
import glob
import math

REJECT = "MKLAVFALLVFFAESSDKARR"
H = set("AILMFWVY")
POS = set("KR")
NEG = set("DE")

def hydro_fraction(seq):
    return sum(a in H for a in seq) / len(seq)

def charge(seq):
    return sum(a in POS for a in seq) - sum(a in NEG for a in seq)

def max_hydrophobic_run(seq):
    best = cur = 0
    for a in seq:
        if a in H:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best

def seq_distance(a, b):
    return sum(x != y for x, y in zip(a, b))

records = {}

for path in glob.glob("../results/ablation-*-619/candidate_history.json"):
    with open(path) as f:
        data = json.load(f)

    def walk(x):
        if isinstance(x, dict):
            seq = x.get("sequence")
            if isinstance(seq, str) and len(seq) == 21 and seq != REJECT:
                records.setdefault(seq, {
                    "sequence": seq,
                    "source_files": set(),
                    "old_scores": []
                })
                records[seq]["source_files"].add(path)

                for key in (
                    "overall",
                    "score",
                    "fitness",
                    "generation_local_score",
                    "posthoc_run_normalized_score"
                ):
                    value = x.get(key)
                    if isinstance(value, (int, float)):
                        records[seq]["old_scores"].append(float(value))

            for v in x.values():
                walk(v)

        elif isinstance(x, list):
            for v in x:
                walk(v)

    walk(data)

rows = []

for r in records.values():
    seq = r["sequence"]
    scores = r["old_scores"]

    rows.append({
        "sequence": seq,
        "old_best_score": max(scores) if scores else float("nan"),
        "hydrophobic_fraction": hydro_fraction(seq),
        "net_charge": charge(seq),
        "max_hydrophobic_run": max_hydrophobic_run(seq),
        "distance_from_rejected": seq_distance(seq, REJECT),
        "runs_seen": len(r["source_files"])
    })

# First-pass computational sanity screen.
# These are TRIAGE heuristics, not biological cutoffs.
screened = [
    r for r in rows
    if r["hydrophobic_fraction"] <= 0.76
    and r["max_hydrophobic_run"] <= 8
]

# Sort primarily by diversity from the rejected candidate,
# not by the old optimization score.
screened.sort(
    key=lambda r: (
        r["distance_from_rejected"],
        -r["runs_seen"],
        -r["old_best_score"] if not math.isnan(r["old_best_score"]) else 0
    ),
    reverse=True
)

print("Unique candidates:", len(rows))
print("Candidates after sequence sanity screen:", len(screened))
print()

print(
    f"{'#':>3} {'Sequence':25s} {'Hydro':>6} "
    f"{'Charge':>7} {'MaxH':>5} {'Dist':>5} {'Runs':>5} {'OldScore':>8}"
)

for i, r in enumerate(screened[:30], 1):
    print(
        f"{i:3d} {r['sequence']:25s} "
        f"{r['hydrophobic_fraction']:6.3f} "
        f"{r['net_charge']:7.1f} "
        f"{r['max_hydrophobic_run']:5d} "
        f"{r['distance_from_rejected']:5d} "
        f"{r['runs_seen']:5d} "
        f"{r['old_best_score']:8.3f}"
    )

# Save the screened panel.
with open("p2_screened.tsv", "w") as f:
    f.write(
        "sequence\thydrophobic_fraction\tnet_charge\t"
        "max_hydrophobic_run\tdistance_from_rejected\t"
        "runs_seen\told_best_score\n"
    )
    for r in screened:
        f.write(
            f"{r['sequence']}\t"
            f"{r['hydrophobic_fraction']:.4f}\t"
            f"{r['net_charge']:.1f}\t"
            f"{r['max_hydrophobic_run']}\t"
            f"{r['distance_from_rejected']}\t"
            f"{r['runs_seen']}\t"
            f"{r['old_best_score']:.4f}\n"
        )
