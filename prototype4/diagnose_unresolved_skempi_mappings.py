import csv
import json
from pathlib import Path
from collections import OrderedDict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path.cwd()
FAMILY_MAP = Path(
    "prototype4/results/skempi_audit_20261029/family_map.csv"
)
STRUCTURES = Path("prototype4/data/structures/mmcif")

def normalize(seq):
    return "".join(str(seq).split()).upper()

def extract_chains(path):
    d = MMCIF2Dict(str(path))
    n = len(d["_atom_site.group_PDB"])

    def column(key, default="?"):
        values = d.get(key)
        if values is None:
            return [default] * n
        if not isinstance(values, list):
            return [values] * n
        if len(values) != n:
            raise ValueError(f"Column mismatch: {key}")
        return values

    groups = column("_atom_site.group_PDB")
    auth = column("_atom_site.auth_asym_id")
    comp = column("_atom_site.auth_comp_id")
    seqid = column("_atom_site.auth_seq_id")
    ins = column("_atom_site.pdbx_PDB_ins_code")
    model = column("_atom_site.pdbx_PDB_model_num", "1")
    alt = column("_atom_site.label_alt_id", ".")

    aa = {
        "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
        "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
        "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
        "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
        "MSE":"M", "SEC":"U", "PYL":"O", "HYP":"P",
    }

    chains = {}

    for i in range(n):
        if groups[i] != "ATOM" or str(model[i]) not in {"1", "1.0"}:
            continue
        if str(alt[i]) not in {".", "?", "A"}:
            continue

        chain = str(auth[i])
        residue = str(comp[i]).upper()
        if chain in {".", "?"} or residue not in aa:
            continue

        key = (str(seqid[i]), str(ins[i]))
        chains.setdefault(chain, OrderedDict())
        chains[chain][key] = aa[residue]

    return {
        chain: "".join(values.values())
        for chain, values in chains.items()
    }

def main():
    with FAMILY_MAP.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    diagnostics = []

    for row in rows:
        pdb = row["pdb_id"].strip().upper()
        path = STRUCTURES / f"{pdb.lower()}.cif"
        observed = extract_chains(path)

        for partner in ("protein1", "protein2"):
            chain_ids = json.loads(
                row[f"{partner}_auth_chain_ids_json"] or "[]"
            )
            stored = [
                normalize(x) for x in json.loads(
                    row[f"{partner}_sequences_json"] or "[]"
                )
            ]

            exact_matches = []
            for i, sequence in enumerate(stored):
                matches = [
                    chain for chain in chain_ids
                    if chain in observed and observed[chain] == sequence
                ]
                exact_matches.append({
                    "stored_sequence_index": i + 1,
                    "length": len(sequence),
                    "matching_chains": matches,
                })

            unmatched_chains = [
                chain for chain in chain_ids
                if chain not in observed
            ]

            unmatched_sequences = [
                x for x in exact_matches if not x["matching_chains"]
            ]

            # Show likely truncation or missing-residue cases by length only.
            # Length proximity is diagnostic, not proof of sequence identity.
            length_candidates = []
            for item in unmatched_sequences:
                candidates = []
                for chain in chain_ids:
                    if chain in observed:
                        delta = abs(
                            len(stored[item["stored_sequence_index"] - 1])
                            - len(observed[chain])
                        )
                        if delta <= 15:
                            candidates.append((chain, delta))
                length_candidates.append({
                    "stored_index": item["stored_sequence_index"],
                    "possible_chains_by_length_only": candidates,
                })

            diagnostics.append({
                "complex_id": row["complex_id"],
                "pdb_id": pdb,
                "partner": partner,
                "family_id": row["family_id"],
                "expected_chains": chain_ids,
                "observed_chains": sorted(observed),
                "stored_sequence_lengths": [len(x) for x in stored],
                "observed_chain_lengths": {
                    c: len(observed[c]) for c in chain_ids if c in observed
                },
                "exact_matches": exact_matches,
                "unmatched_chains": unmatched_chains,
                "length_candidates_not_identity_evidence": length_candidates,
            })

    unresolved = [
        x for x in diagnostics
        if x["unmatched_chains"] or any(
            not m["matching_chains"] for m in x["exact_matches"]
        )
    ]

    print("UNRESOLVED PARTNER GROUPS:", len(unresolved))

    for x in unresolved:
        print(
            f"\n{x['pdb_id']} | {x['partner']} | "
            f"family={x['family_id']}"
        )
        print("  Expected chains:", x["expected_chains"])
        print("  Observed chains:", x["observed_chains"])
        print("  Stored sequence lengths:", x["stored_sequence_lengths"])
        print("  Observed chain lengths:", x["observed_chain_lengths"])
        print("  Exact sequence matches:")
        for m in x["exact_matches"]:
            print(
                f"    stored#{m['stored_sequence_index']} "
                f"({m['length']} aa) -> {m['matching_chains']}"
            )
        if x["unmatched_chains"]:
            print("  Missing chains:", x["unmatched_chains"])

    out = Path(
        "prototype4/results/data_isolation_20261101/"
        "unresolved_mapping_diagnostics.json"
    )
    out.write_text(json.dumps(unresolved, indent=2) + "\n")
    print("\nFull diagnostic saved:", out)

if __name__ == "__main__":
    main()
