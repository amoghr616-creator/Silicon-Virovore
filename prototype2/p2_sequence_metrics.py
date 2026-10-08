#!/usr/bin/env python3
"""Prototype 2 sequence/physicochemical metrics for the frozen six-candidate panel.

Computes only metrics that can be justified from sequence alone with the current
pipeline assumptions. Structure/docking/MD/ADME fields are intentionally not
fabricated; they are emitted as NA until measured/predicted by a validated module.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis

PANEL = Path(__file__).with_name("p2_panel.tsv")
SCREENED = Path(__file__).with_name("p2_screened.tsv")
OUT = Path(__file__).with_name("p2_sequence_metrics.tsv")

# Matches Silicon Virovore's current Eisenberg implementation.
EISENBERG = {
    "I":0.73,"F":0.61,"V":0.54,"L":0.53,"W":0.37,"M":0.26,"A":0.25,
    "G":0.16,"C":0.04,"Y":0.02,"P":-0.07,"T":-0.18,"S":-0.26,"H":-0.40,
    "E":-0.62,"N":-0.64,"Q":-0.85,"D":-0.72,"K":-1.10,"R":-1.76,
}
KD = {
    "A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,
    "G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,
    "P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2,
}
HYDROPHOBIC = set("AILMFWVY")
POS = set("KR")
NEG = set("DE")


def eisenberg_moment(seq: str) -> float:
    angle = math.radians(100.0)
    c = s = 0.0
    for i, aa in enumerate(seq):
        h = EISENBERG[aa]
        theta = (i + 1) * angle
        c += h * math.cos(theta)
        s += h * math.sin(theta)
    return math.hypot(c, s) / len(seq)


def max_hydrophobic_run(seq: str) -> int:
    best = cur = 0
    for aa in seq:
        if aa in HYDROPHOBIC:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def mean_kd(seq: str) -> float:
    return sum(KD[a] for a in seq) / len(seq)


def max_kd_window(seq: str, window: int = 7) -> float:
    vals = [sum(KD[a] for a in seq[i:i+window]) / window for i in range(len(seq)-window+1)]
    return max(vals)


def load_screened() -> dict[str, dict[str, str]]:
    if not SCREENED.exists():
        return {}
    with SCREENED.open() as f:
        rows = csv.DictReader(f, delimiter="\t")
        return {r["sequence"]: r for r in rows}

def load_panel() -> list[tuple[str, str]]:
    rows = []
    with PANEL.open() as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) >= 2 and row[0].startswith("P2-"):
                rows.append((row[0], row[1]))
    if len(rows) != 6:
        raise SystemExit(f"Expected exactly 6 P2 candidates in {PANEL}, found {len(rows)}")
    return rows


def main() -> None:
    screen = load_screened()
    rows = []
    for cid, seq in load_panel():
        old = screen.get(seq, {})
        pa = ProteinAnalysis(seq)
        rows.append({
            "candidate": cid,
            "sequence": seq,
            "length_aa": len(seq),
            "molecular_weight_da": round(pa.molecular_weight(), 3),
            "theoretical_pI": round(pa.isoelectric_point(), 3),
            "charge_pH_5": round(pa.charge_at_pH(5.0), 3),
            "charge_pH_6": round(pa.charge_at_pH(6.0), 3),
            "charge_pH_7": round(pa.charge_at_pH(7.0), 3),
            "charge_pH_7_4": round(pa.charge_at_pH(7.4), 3),
            "charge_pH_8": round(pa.charge_at_pH(8.0), 3),
            "charge_pH_9": round(pa.charge_at_pH(9.0), 3),
            "gravy_hydropathy": round(pa.gravy(), 3),
            "mean_kd_hydropathy": round(mean_kd(seq), 3),
            "max_kd_7mer": round(max_kd_window(seq, 7), 3),
            "hydrophobic_fraction": round(sum(a in HYDROPHOBIC for a in seq)/len(seq), 3),
            "max_hydrophobic_run": max_hydrophobic_run(seq),
            "net_charge_simple": sum(a in POS for a in seq) - sum(a in NEG for a in seq),
            "eisenberg_hydrophobic_moment": round(eisenberg_moment(seq), 4),
            "instability_index": round(pa.instability_index(), 3),
            "aromaticity_fraction": round(pa.aromaticity(), 3),
            "old_best_score": old.get("old_best_score", "NA"),
            "distance_from_rejected": old.get("distance_from_rejected", "NA"),
            "runs_seen": old.get("runs_seen", "NA"),
            "all_L_canonical": "yes",
            "structure_secondary_structure": "NA",
            "sasa": "NA",
            "rmsd": "NA",
            "rmsf": "NA",
            "conformational_entropy": "NA",
            "docking_affinity_kcal_mol": "NA",
            "md_stability": "NA",
            "proteolytic_half_life": "NA",
            "plasma_protein_binding": "NA",
            "clearance": "NA",
            "cytotoxicity": "NA",
            "hemolysis": "NA",
            "immunogenicity": "NA",
            "synthetic_feasibility": "NA",
            "cost_estimate": "NA",
        })

    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {OUT}")
    for r in rows:
        print(r["candidate"], r["sequence"], "MW=", r["molecular_weight_da"], "pI=", r["theoretical_pI"], "q7.4=", r["charge_pH_7_4"], "GRAVY=", r["gravy_hydropathy"], "muH=", r["eisenberg_hydrophobic_moment"])

if __name__ == "__main__":
    main()
