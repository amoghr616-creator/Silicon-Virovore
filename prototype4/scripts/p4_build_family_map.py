from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import re
import time
import urllib.request

from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from Bio import __version__ as biopython_version
from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


ROOT = Path.cwd()

AUDIT = ROOT / "prototype4/results/skempi_audit_20261029"
SOURCE = ROOT / "prototype4/data/source/skempi_v2.csv"
ELIGIBLE = AUDIT / "eligible_cycles.csv"
CACHE = ROOT / "prototype4/data/structures/mmcif"

OUT_MAP = AUDIT / "family_map.csv"
OUT_SEQUENCES = AUDIT / "sequence_inventory.csv"
OUT_EDGES = AUDIT / "family_similarity_edges.csv"
OUT_COUNTS = AUDIT / "family_group_counts.csv"
OUT_MUTATIONS = AUDIT / "mutation_chain_audit.csv"
OUT_MANIFEST = AUDIT / "family_map_manifest.json"

PRIMARY_IDENTITY = 0.30
COVERAGE_THRESHOLD = 0.70
MIN_LONG_ALIGNMENT = 30

IDENTITY_THRESHOLDS = [0.25, 0.30, 0.40]

CACHE.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def clean_sequence(value: str) -> str:
    value = str(value).upper()
    value = re.sub(r"\s+", "", value)
    value = re.sub(r"[^A-Z]", "", value)
    return value


def fetch_cif(pdb_id: str) -> Path:
    pdb_id = pdb_id.upper()
    path = CACHE / f"{pdb_id.lower()}.cif"

    if path.exists() and path.stat().st_size > 1000:
        return path

    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.cif.gz"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SiliconVirovore-P4/1.0 research"},
    )

    last_error = None

    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()

            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)

            if len(data) < 1000 or b"data_" not in data[:1000]:
                raise ValueError("Downloaded file does not look like mmCIF.")

            path.write_bytes(data)
            return path

        except Exception as exc:
            last_error = exc
            time.sleep(2 ** attempt)

    raise RuntimeError(
        f"Unable to download {pdb_id} from RCSB: {last_error}"
    )


def load_author_chain_sequences(pdb_id: str):
    cif_path = fetch_cif(pdb_id)
    data = MMCIF2Dict(str(cif_path))

    entity_ids = as_list(data.get("_entity_poly.entity_id"))
    canonical = as_list(
        data.get("_entity_poly.pdbx_seq_one_letter_code_can")
    )

    if not entity_ids or not canonical:
        raise ValueError(f"{pdb_id}: no polymer sequences in mmCIF")

    entity_sequences = {}

    for entity_id, sequence in zip(entity_ids, canonical):
        sequence = clean_sequence(sequence)
        if sequence:
            entity_sequences[str(entity_id)] = sequence

    label_ids = as_list(data.get("_struct_asym.id"))
    label_entities = as_list(data.get("_struct_asym.entity_id"))

    label_to_entity = {
        str(label): str(entity)
        for label, entity in zip(label_ids, label_entities)
    }

    atom_auth = as_list(data.get("_atom_site.auth_asym_id"))
    atom_label = as_list(data.get("_atom_site.label_asym_id"))

    if len(atom_auth) != len(atom_label):
        raise ValueError(f"{pdb_id}: inconsistent atom-site chain columns")

    auth_to_labels = defaultdict(set)

    for auth, label in zip(atom_auth, atom_label):
        auth = str(auth)
        label = str(label)

        if auth not in {"", ".", "?"}:
            if label in label_to_entity:
                auth_to_labels[auth].add(label)

    result = {}

    for auth, labels in auth_to_labels.items():
        found = set()

        for label in labels:
            entity = label_to_entity.get(label)
            sequence = entity_sequences.get(entity)

            if sequence:
                found.add(sequence)

        if len(found) == 1:
            result[auth] = {
                "sequence": next(iter(found)),
                "label_asym_ids": sorted(labels),
            }
        elif len(found) > 1:
            raise ValueError(
                f"{pdb_id}: author chain {auth} maps to "
                f"multiple distinct polymer sequences"
            )

    if not result:
        raise ValueError(f"{pdb_id}: no author-to-polymer sequence mappings")

    return result, sha256(cif_path)


def resolve_chain_group(group: str, available: set[str]):
    group = str(group).strip()

    if not group:
        raise ValueError("Empty chain group")

    if any(separator in group for separator in [",", ";"]):
        parts = [
            part.strip()
            for part in re.split(r"[,;]", group)
            if part.strip()
        ]

        if all(part in available for part in parts):
            return parts

        raise ValueError(f"Unresolved chain group {group!r}")

    exact = group in available

    if len(group) > 1:
        split_chars = list(group)
        split_valid = all(char in available for char in split_chars)

        if exact and split_valid:
            raise ValueError(
                f"Ambiguous chain group {group!r}: "
                "both a multi-character chain and multiple "
                "single-character chains are possible"
            )

        if split_valid:
            return split_chars

    if exact:
        return [group]

    raise ValueError(
        f"Chain group {group!r} not found. "
        f"Available author chain IDs: {sorted(available)}"
    )


def parse_complex_id(complex_id: str):
    parts = str(complex_id).split("_", 2)

    if len(parts) != 3:
        raise ValueError(f"Unexpected complex_id: {complex_id}")

    pdb_id, group1, group2 = parts

    if not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
        raise ValueError(f"Unexpected PDB identifier: {pdb_id}")

    return pdb_id.upper(), group1, group2


def parse_mutation_chain(mutation: str):
    parts = str(mutation).strip().split(":")

    # SKEMPI cleaned format:
    # original_residue : author_chain : residue_number : mutant_residue
    if len(parts) < 4:
        return None

    return parts[1].strip() or None


def similarity(seq1: str, seq2: str, aligner):
    if min(len(seq1), len(seq2)) < 8:
        return None

    alignment = aligner.align(seq1, seq2)[0]
    indices = alignment.indices

    valid = (indices[0] >= 0) & (indices[1] >= 0)

    aligned_length = int(valid.sum())

    if aligned_length == 0:
        return None

    i1 = indices[0][valid]
    i2 = indices[1][valid]

    matches = sum(
        seq1[int(a)] == seq2[int(b)]
        for a, b in zip(i1, i2)
    )

    identity = matches / aligned_length
    shorter_coverage = aligned_length / min(len(seq1), len(seq2))
    longer_coverage = aligned_length / max(len(seq1), len(seq2))

    return {
        "identity": identity,
        "aligned_length": aligned_length,
        "shorter_coverage": shorter_coverage,
        "longer_coverage": longer_coverage,
    }


def matches_family_rule(metrics, threshold):
    if metrics is None:
        return False

    # Short chains require near-identical, nearly complete matches.
    # This avoids calling a short conserved motif a protein family.
    if metrics["aligned_length"] < MIN_LONG_ALIGNMENT:
        return (
            metrics["aligned_length"] >= 8
            and metrics["identity"] >= 0.90
            and metrics["shorter_coverage"] >= 0.90
        )

    return (
        metrics["identity"] >= threshold
        and metrics["shorter_coverage"] >= COVERAGE_THRESHOLD
    )


class UnionFind:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def find(self, value):
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)

        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def connected_components(values, edges):
    uf = UnionFind(values)

    for edge in edges:
        uf.union(edge["complex_a"], edge["complex_b"])

    components = defaultdict(list)

    for value in values:
        components[uf.find(value)].append(value)

    return sorted(
        [sorted(members) for members in components.values()],
        key=lambda members: members[0],
    )


def main():
    for path in [SOURCE, ELIGIBLE]:
        if not path.is_file():
            raise FileNotFoundError(path)

    raw = pd.read_csv(SOURCE, sep=";", dtype=str).fillna("")
    cycles = pd.read_csv(ELIGIBLE)

    required_raw = {"#Pdb", "Protein 1", "Protein 2"}
    required_cycles = {
        "cycle_uid", "complex_id", "mutation_i", "mutation_j"
    }

    if not required_raw.issubset(raw.columns):
        raise ValueError("Missing required raw SKEMPI columns")

    if not required_cycles.issubset(cycles.columns):
        raise ValueError("Missing required eligible-cycle columns")

    complex_ids = sorted(cycles["complex_id"].astype(str).unique())

    print(f"Eligible cycles: {len(cycles)}")
    print(f"Distinct eligible complexes: {len(complex_ids)}")
    print(f"Biopython: {biopython_version}")

    # Protein names are descriptive source annotations, not chain identifiers.
    # Preserve every distinct label instead of assuming one label per #Pdb.

    def collect_labels(series):
        return sorted({
            str(value).strip()
            for value in series
            if str(value).strip()
        })

    raw_meta = raw.groupby("#Pdb", sort=False).agg(
        protein1_labels_json=(
            "Protein 1",
            lambda series: json.dumps(collect_labels(series))
        ),
        protein2_labels_json=(
            "Protein 2",
            lambda series: json.dumps(collect_labels(series))
        ),
    )

    # Explicit audit of changing descriptive labels.
    label_audit = raw_meta.reset_index().copy()

    label_audit["protein1_label_count"] = (
        label_audit["protein1_labels_json"]
        .map(lambda value: len(json.loads(value)))
    )

    label_audit["protein2_label_count"] = (
        label_audit["protein2_labels_json"]
        .map(lambda value: len(json.loads(value)))
    )

    label_audit.to_csv(
        AUDIT / "protein_label_variation_audit.csv",
        index=False,
    )

    variable_labels = label_audit[
        (label_audit["protein1_label_count"] > 1)
        | (label_audit["protein2_label_count"] > 1)
    ]

    print(
        "Complexes with multiple source protein labels:",
        len(variable_labels),
    )

    complexes = {}
    cif_hashes = {}
    unresolved = []

    # ---------------------------------------------------------
    # 1. Download structures and map author chain IDs to sequences.
    # ---------------------------------------------------------

    for complex_id in complex_ids:
        pdb_id, group1, group2 = parse_complex_id(complex_id)

        try:
            if complex_id not in raw_meta.index:
                raise ValueError("No matching #Pdb source record")

            author_sequences, cif_hash = load_author_chain_sequences(pdb_id)
            available = set(author_sequences)

            p1_chains = resolve_chain_group(group1, available)
            p2_chains = resolve_chain_group(group2, available)

            def make_partner(chains, partner):
                rows = []

                for chain in chains:
                    info = author_sequences[chain]

                    rows.append({
                        "auth_chain_id": chain,
                        "partner": partner,
                        "sequence": info["sequence"],
                        "length": len(info["sequence"]),
                        "label_asym_ids": info["label_asym_ids"],
                        "sequence_sha256": hashlib.sha256(
                            info["sequence"].encode()
                        ).hexdigest(),
                    })

                return rows

            p1 = make_partner(p1_chains, "protein1")
            p2 = make_partner(p2_chains, "protein2")

            if not p1 or not p2:
                raise ValueError("One interaction partner has no sequences")

            p1_labels = json.loads(
                raw_meta.loc[complex_id, "protein1_labels_json"]
            )
            p2_labels = json.loads(
                raw_meta.loc[complex_id, "protein2_labels_json"]
            )

            complexes[complex_id] = {
                "pdb_id": pdb_id,
                "group1": group1,
                "group2": group2,
                "protein1": " | ".join(p1_labels),
                "protein2": " | ".join(p2_labels),
                "protein1_name_labels": p1_labels,
                "protein2_name_labels": p2_labels,
                "p1": p1,
                "p2": p2,
                "all_chains": p1 + p2,
                "p1_ids": set(p1_chains),
                "p2_ids": set(p2_chains),
                "cif_sha256": cif_hash,
            }

            cif_hashes[pdb_id] = cif_hash

            print(
                f"Mapped {complex_id}: "
                f"P1={sum(x['length'] for x in p1)} aa, "
                f"P2={sum(x['length'] for x in p2)} aa"
            )

        except Exception as exc:
            unresolved.append({
                "complex_id": complex_id,
                "pdb_id": pdb_id,
                "reason": str(exc),
            })
            print(f"UNRESOLVED {complex_id}: {exc}")

    if unresolved:
        pd.DataFrame(unresolved).to_csv(
            AUDIT / "family_mapping_unresolved.csv",
            index=False,
        )

        raise RuntimeError(
            f"{len(unresolved)} complexes could not be mapped. "
            "Resolve these before family-held-out validation. "
            "See family_mapping_unresolved.csv"
        )

    # ---------------------------------------------------------
    # 2. Inventory all mapped chain sequences.
    # ---------------------------------------------------------

    sequence_rows = []

    for complex_id, record in complexes.items():
        for chain in record["all_chains"]:
            sequence_rows.append({
                "complex_id": complex_id,
                "pdb_id": record["pdb_id"],
                "partner": chain["partner"],
                "auth_chain_id": chain["auth_chain_id"],
                "label_asym_ids": json.dumps(chain["label_asym_ids"]),
                "sequence_length": chain["length"],
                "sequence_sha256": chain["sequence_sha256"],
                "sequence": chain["sequence"],
            })

    pd.DataFrame(sequence_rows).to_csv(
        OUT_SEQUENCES,
        index=False,
    )

    # ---------------------------------------------------------
    # 3. Audit mutation-chain membership.
    # ---------------------------------------------------------

    mutation_rows = []

    for row in cycles.itertuples(index=False):
        record = complexes[str(row.complex_id)]

        chain_i = parse_mutation_chain(row.mutation_i)
        chain_j = parse_mutation_chain(row.mutation_j)

        def classify(chain):
            if chain is None:
                return "unparsed"
            in_p1 = chain in record["p1_ids"]
            in_p2 = chain in record["p2_ids"]

            if in_p1 and in_p2:
                return "ambiguous_both_partners"
            if in_p1:
                return "protein1"
            if in_p2:
                return "protein2"
            return "not_in_declared_partner_groups"

        partner_i = classify(chain_i)
        partner_j = classify(chain_j)

        mutation_rows.append({
            "cycle_uid": row.cycle_uid,
            "complex_id": row.complex_id,
            "mutation_i": row.mutation_i,
            "mutation_i_chain": chain_i,
            "mutation_i_partner": partner_i,
            "mutation_j": row.mutation_j,
            "mutation_j_chain": chain_j,
            "mutation_j_partner": partner_j,
            "both_mutations_mapped": (
                partner_i in {"protein1", "protein2"}
                and partner_j in {"protein1", "protein2"}
            ),
            "cross_partner_cycle": (
                partner_i in {"protein1", "protein2"}
                and partner_j in {"protein1", "protein2"}
                and partner_i != partner_j
            ),
        })

    mutation_df = pd.DataFrame(mutation_rows)
    mutation_df.to_csv(OUT_MUTATIONS, index=False)

    # ---------------------------------------------------------
    # 4. Prepare local pairwise alignment.
    # ---------------------------------------------------------

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    similarity_cache = {}

    def get_similarity(a, b):
        key = tuple(sorted((a, b)))

        if key not in similarity_cache:
            similarity_cache[key] = similarity(a, b, aligner)

        return similarity_cache[key]

    # ---------------------------------------------------------
    # 5. Compare every pair of complexes.
    # ---------------------------------------------------------

    ids = sorted(complexes)
    pair_candidates = []

    for complex_a, complex_b in itertools.combinations(ids, 2):
        a = complexes[complex_a]
        b = complexes[complex_b]

        if a["pdb_id"] == b["pdb_id"]:
            pair_candidates.append({
                "complex_a": complex_a,
                "complex_b": complex_b,
                "reason": "same_pdb_entry",
                "matches": [],
            })
            continue

        matches = []

        for ca in a["all_chains"]:
            for cb in b["all_chains"]:
                metrics = get_similarity(
                    ca["sequence"],
                    cb["sequence"],
                )

                if metrics is None:
                    continue

                matches.append({
                    "chain_a": ca["auth_chain_id"],
                    "partner_a": ca["partner"],
                    "length_a": ca["length"],
                    "chain_b": cb["auth_chain_id"],
                    "partner_b": cb["partner"],
                    "length_b": cb["length"],
                    **metrics,
                })

        pair_candidates.append({
            "complex_a": complex_a,
            "complex_b": complex_b,
            "reason": "sequence_similarity",
            "matches": matches,
        })

    # ---------------------------------------------------------
    # 6. Build family groups at three identity thresholds.
    # ---------------------------------------------------------

    all_edge_rows = []
    group_count_rows = []
    primary_assignments = {}

    cycle_counts = cycles.groupby("complex_id").size().to_dict()

    for threshold in IDENTITY_THRESHOLDS:
        edges = []

        for pair in pair_candidates:
            if pair["reason"] == "same_pdb_entry":
                edges.append({
                    "complex_a": pair["complex_a"],
                    "complex_b": pair["complex_b"],
                    "reason": "same_pdb_entry",
                })
                continue

            passing = [
                match
                for match in pair["matches"]
                if matches_family_rule(match, threshold)
            ]

            if not passing:
                continue

            best = max(
                passing,
                key=lambda x: (
                    x["identity"],
                    x["shorter_coverage"],
                    x["aligned_length"],
                ),
            )

            edge = {
                "complex_a": pair["complex_a"],
                "complex_b": pair["complex_b"],
                "reason": "sequence_homology",
                **best,
            }

            edges.append(edge)

        components = connected_components(ids, edges)

        membership = {}

        for members in components:
            digest = hashlib.sha256(
                "|".join(members).encode()
            ).hexdigest()[:10]

            family_id = f"IF{int(threshold * 100)}_{digest}"

            for complex_id in members:
                membership[complex_id] = family_id

            family_cycles = sum(cycle_counts.get(x, 0) for x in members)

            group_count_rows.append({
                "identity_threshold": threshold,
                "shorter_coverage_threshold": COVERAGE_THRESHOLD,
                "n_families": len(components),
                "n_complexes": len(ids),
                "n_cycles": len(cycles),
                "largest_family_complexes": len(members),
                "largest_family_cycles": family_cycles,
                "family_id": family_id,
                "member_complexes_json": json.dumps(members),
            })

        for edge in edges:
            all_edge_rows.append({
                "identity_threshold": threshold,
                "shorter_coverage_threshold": COVERAGE_THRESHOLD,
                **edge,
            })

        if abs(threshold - PRIMARY_IDENTITY) < 1e-9:
            primary_assignments = membership

    # ---------------------------------------------------------
    # 7. Export completed family map.
    # ---------------------------------------------------------

    family_sizes = defaultdict(list)

    for complex_id, family_id in primary_assignments.items():
        family_sizes[family_id].append(complex_id)

    family_cycle_sizes = {
        family_id: sum(cycle_counts.get(x, 0) for x in members)
        for family_id, members in family_sizes.items()
    }

    map_rows = []

    for complex_id in ids:
        record = complexes[complex_id]
        family_id = primary_assignments[complex_id]

        p1_length = sum(x["length"] for x in record["p1"])
        p2_length = sum(x["length"] for x in record["p2"])

        map_rows.append({
            "complex_id": complex_id,
            "protein1": record["protein1"],
            "protein2": record["protein2"],
            "protein1_name_labels_json": json.dumps(
                record["protein1_name_labels"]
            ),
            "protein2_name_labels_json": json.dumps(
                record["protein2_name_labels"]
            ),
            "family_id": family_id,
            "family_notes": (
                "Interaction-family connected component; "
                "homology across either partner; "
                "30% identity over >=70% of shorter chain; "
                "short-chain matches require >=90% identity and coverage; "
                "same PDB entries grouped."
            ),
            "protein1_length": p1_length,
            "protein2_length": p2_length,
            "peptide_partner": "",
            "pdb_id": record["pdb_id"],
            "protein1_chain_group": record["group1"],
            "protein2_chain_group": record["group2"],
            "protein1_auth_chain_ids_json": json.dumps(
                [x["auth_chain_id"] for x in record["p1"]]
            ),
            "protein2_auth_chain_ids_json": json.dumps(
                [x["auth_chain_id"] for x in record["p2"]]
            ),
            "protein1_sequences_json": json.dumps(
                [x["sequence"] for x in record["p1"]]
            ),
            "protein2_sequences_json": json.dumps(
                [x["sequence"] for x in record["p2"]]
            ),
            "family_complex_count": len(family_sizes[family_id]),
            "family_cycle_count": family_cycle_sizes[family_id],
            "mapping_status": "sequence_mapping_verified",
            "peptide_partner_status": (
                "not_assigned; requires biological annotation"
            ),
        })

    family_map = pd.DataFrame(map_rows)
    family_map.to_csv(OUT_MAP, index=False)

    pd.DataFrame(all_edge_rows).to_csv(
        OUT_EDGES,
        index=False,
    )

    pd.DataFrame(group_count_rows).to_csv(
        OUT_COUNTS,
        index=False,
    )

    # ---------------------------------------------------------
    # 8. Manifest and summary.
    # ---------------------------------------------------------

    family_count = len(set(primary_assignments.values()))
    primary_largest = max(family_cycle_sizes.values())

    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": sha256(SOURCE),
        "eligible_cycles_sha256": sha256(ELIGIBLE),
        "script_sha256": sha256(Path(__file__)),
        "biopython_version": biopython_version,
        "primary_identity_threshold": PRIMARY_IDENTITY,
        "shorter_sequence_coverage_threshold": COVERAGE_THRESHOLD,
        "minimum_long_alignment": MIN_LONG_ALIGNMENT,
        "short_chain_rule": {
            "minimum_length": 8,
            "identity": 0.90,
            "shorter_coverage": 0.90,
        },
        "n_eligible_cycles": int(len(cycles)),
        "n_complexes": int(len(ids)),
        "n_mapped_complexes": int(len(complexes)),
        "n_family_groups_primary": int(family_count),
        "largest_family_cycle_count": int(primary_largest),
        "n_mutation_rows": int(len(mutation_df)),
        "n_mutations_both_mapped": int(
            mutation_df["both_mutations_mapped"].sum()
        ),
        "n_cross_partner_cycles": int(
            mutation_df["cross_partner_cycle"].sum()
        ),
        "n_unmapped_mutation_rows": int(
            (~mutation_df["both_mutations_mapped"]).sum()
        ),
        "protein_label_variation_audit": str(
            AUDIT / "protein_label_variation_audit.csv"
        ),
        "n_complexes_with_multiple_protein_labels": int(
            len(variable_labels)
        ),
        "pdb_mmcif_sha256": cif_hashes,
        "family_map": str(OUT_MAP),
        "sequence_inventory": str(OUT_SEQUENCES),
        "similarity_edges": str(OUT_EDGES),
        "threshold_sensitivity": str(OUT_COUNTS),
        "mutation_chain_audit": str(OUT_MUTATIONS),
        "interpretation": (
            "Sequence-based interaction-family grouping for SKEMPI "
            "protein-protein interaction cycles. Not a peptide-only "
            "or receptor-only family assignment."
        ),
    }

    OUT_MANIFEST.write_text(json.dumps(audit, indent=2))

    print("\nP4 FAMILY MAPPING COMPLETE")
    print("=" * 64)
    print(f"Cycles: {len(cycles)}")
    print(f"Complexes mapped: {len(complexes)}")
    print(f"Primary family groups: {family_count}")
    print(f"Largest family: {primary_largest} cycles")
    print(
        "Mutation rows with both chains mapped: "
        f"{audit['n_mutations_both_mapped']}/{len(mutation_df)}"
    )
    print(
        "Cross-partner mutation cycles: "
        f"{audit['n_cross_partner_cycles']}"
    )
    print("\nThreshold sensitivity:")
    print(pd.DataFrame(group_count_rows).groupby(
        "identity_threshold"
    )["n_families"].first().to_string())

    print("\nOutput files:")
    for path in [
        OUT_MAP,
        OUT_SEQUENCES,
        OUT_EDGES,
        OUT_COUNTS,
        OUT_MUTATIONS,
        OUT_MANIFEST,
    ]:
        print(path)


if __name__ == "__main__":
    main()
