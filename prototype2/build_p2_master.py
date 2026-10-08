from pathlib import Path
import csv

CANDIDATES = {
    "P2-01": "MKQAVNALLVFFAGSSDAIRR",
    "P2-02": "DKLAVTALLVVFAESSDKLRR",
    "P2-03": "MKLAPFALLVVFAGASDWIRR",
    "P2-04": "MKLAVFALLVCFAESSDLVRR",
    "P2-05": "MKQAVFALLQVFAGSSDWIRR",
    "P2-06": "MKLAVFALRVFFAESSDAIRR",
}

FIELDS = [
    "candidate",
    "sequence",

    # Sequence / physicochemical
    "length",
    "molecular_weight_da",
    "pi",
    "charge_pH5",
    "charge_pH6",
    "charge_pH7",
    "charge_pH7_4",
    "charge_pH8",
    "charge_pH9",
    "gravy",
    "hydrophobic_moment",
    "hydrophobic_fraction",
    "max_hydrophobic_run",
    "amino_acid_composition",
    "stereochemistry",

    # Structure
    "structure_predicted",
    "structure_confidence_plddt",
    "secondary_structure",
    "structure_ca_rmsd_experimental_A",
    "sasa_A2",
    "radius_gyration_A",
    "structural_rmsd_A",
    "structural_rmsf_A",
    "conformational_variability",

    # Docking / interaction
    "docking_method",
    "docking_score_kcal_mol",
    "docking_pose_rmsd_experimental_A",
    "interface_contacts",
    "hydrogen_bonds",
    "salt_bridges",
    "contact_recall",
    "contact_jaccard",
    "pose_consistency",

    # MD
    "md_rmsd_A",
    "md_rmsf_A",
    "md_radius_gyration_A",
    "md_contact_persistence",
    "md_hbond_persistence",
    "md_stability_status",

    # Developability / safety predictions
    "protease_susceptibility",
    "predicted_half_life",
    "toxicity_prediction",
    "hemolysis_prediction",
    "immunogenicity_prediction",
    "synthetic_feasibility",
    "estimated_cost",

    # Algorithm / statistics
    "ga_fitness",
    "convergence_generation",
    "sequence_diversity",
    "arise_importance",
    "bootstrap_stability",
    "pareto_rank",

    # Evidence bookkeeping
    "evidence_status",
]

rows = []
for candidate, seq in CANDIDATES.items():
    row = {field: "NA" for field in FIELDS}
    row["candidate"] = candidate
    row["sequence"] = seq
    row["length"] = len(seq)
    row["stereochemistry"] = "L-canonical amino acids"
    row["evidence_status"] = "computational screening; downstream validation pending"
    rows.append(row)

out = Path("prototype2/p2_master_results.csv")
out.parent.mkdir(parents=True, exist_ok=True)

with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {out}")
print(f"Candidates: {len(rows)}")
print(f"Metrics per candidate: {len(FIELDS)}")
