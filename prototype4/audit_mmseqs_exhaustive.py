import csv
import json
import shutil
import subprocess
from pathlib import Path

OUT = Path(
    "prototype4/results/data_isolation_20261101/"
    "interface_window_mmseqs"
)

QUERY_FASTA = OUT / "skempi_interface_windows.fasta"
TARGET_FASTA = OUT / "hervk_9sya_target_windows.fasta"

QUERY_MAP = OUT / "skempi_interface_window_map.tsv"
TARGET_MAP = OUT / "hervk_target_window_map.tsv"

RAW = OUT / "mmseqs_exhaustive_raw.tsv"
TMP = OUT / "mmseqs_exhaustive_tmp"

FIELDS = [
    "query", "target", "pident", "alnlen",
    "qstart", "qend", "tstart", "tend",
    "qcov", "tcov", "evalue", "bits",
    "qaln", "taln",
]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def normalized_coverage(value):
    x = float(value)
    return x / 100.0 if x > 1.0 else x


def normalized_identity(value):
    x = float(value)
    return x * 100.0 if 0.0 <= x <= 1.0 else x


def aligned_interface_positions(hit, query_meta, target_meta):
    """
    Map interface positions to actual aligned residue columns.
    Gaps are handled explicitly; a simple start/end interval
    is not sufficient for gapped alignments.
    """
    qaln = hit["qaln"]
    taln = hit["taln"]

    if len(qaln) != len(taln):
        raise ValueError("Aligned query and target lengths differ")

    q_interface = set(
        json.loads(query_meta["interface_positions_relative_json"])
    )
    t_interface = set(
        json.loads(target_meta["interface_positions_relative_json"])
    )

    qpos = int(hit["qstart"])
    tpos = int(hit["tstart"])

    q_covered = set()
    t_covered = set()
    paired_interface_columns = []

    for qa, ta in zip(qaln, taln):
        q_has = qa != "-"
        t_has = ta != "-"

        if q_has and t_has:
            if qpos in q_interface:
                q_covered.add(qpos)
            if tpos in t_interface:
                t_covered.add(tpos)

            if qpos in q_interface and tpos in t_interface:
                paired_interface_columns.append((qpos, tpos))

        if q_has:
            qpos += 1
        if t_has:
            tpos += 1

    return {
        "query_interface_positions_aligned": len(q_covered),
        "target_interface_positions_aligned": len(t_covered),
        "paired_interface_alignment_columns": len(
            paired_interface_columns
        ),
        "query_interface_positions_json": json.dumps(
            sorted(q_covered)
        ),
        "target_interface_positions_json": json.dumps(
            sorted(t_covered)
        ),
        "paired_interface_columns_json": json.dumps(
            paired_interface_columns
        ),
    }


def main():
    for path in (QUERY_FASTA, TARGET_FASTA, QUERY_MAP, TARGET_MAP):
        if not path.exists():
            raise FileNotFoundError(path)

    mmseqs = shutil.which("mmseqs")
    if not mmseqs:
        raise RuntimeError("MMseqs2 not found in PATH")

    query_rows = read_tsv(QUERY_MAP)
    target_rows = read_tsv(TARGET_MAP)

    query_map = {r["window_id"]: r for r in query_rows}
    target_map = {r["window_id"]: r for r in target_rows}

    print("QUERY WINDOWS:", len(query_map))
    print("TARGET WINDOWS:", len(target_map))

    print("\nTARGET WINDOW DETAILS")
    for r in target_rows:
        print({
            "id": r["window_id"],
            "chain": r["auth_chain"],
            "length": r["sequence_length"],
            "interface_residues": r["interface_residue_count"],
            "interface_density": r["interface_density"],
            "sequence": r["sequence"],
        })

    if TMP.exists():
        shutil.rmtree(TMP)
    if RAW.exists():
        RAW.unlink()

    command = [
        mmseqs, "easy-search",
        str(QUERY_FASTA),
        str(TARGET_FASTA),
        str(RAW),
        str(TMP),

        # Exhaustive search: bypass normal k-mer prefiltering.
        "--exhaustive-search",

        # Capture weak local alignments for post-search filtering.
        "-e", "1000000",
        "--min-seq-id", "0.0",
        "-c", "0.0",
        "--cov-mode", "0",

        # Calculate real alignment identity.
        "--alignment-mode", "3",
        "--max-seqs", "100000",
        "--threads", "4",
        "--format-output", ",".join(FIELDS),
    ]

    print("\nRUNNING EXHAUSTIVE MMSEQS2 SEARCH")
    subprocess.run(command, check=True)

    enriched = []

    if RAW.exists():
        with RAW.open(encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f if line.strip()]
    else:
        lines = []

    for line in lines:
        parts = line.split("\t")
        if len(parts) != len(FIELDS):
            raise ValueError(
                f"Expected {len(FIELDS)} output fields, "
                f"received {len(parts)}"
            )

        hit = dict(zip(FIELDS, parts))

        qmeta = query_map[hit["query"]]
        tmeta = target_map[hit["target"]]

        identity = normalized_identity(hit["pident"])
        qcov = normalized_coverage(hit["qcov"])
        tcov = normalized_coverage(hit["tcov"])
        evalue = float(hit["evalue"])

        interface = aligned_interface_positions(
            hit, qmeta, tmeta
        )

        row = {
            **hit,
            "identity_percent_normalized": identity,
            "query_coverage_normalized": qcov,
            "target_coverage_normalized": tcov,
            "evalue_numeric": evalue,
            "query_pdb": qmeta["pdb_id"],
            "query_family": qmeta["family_id"],
            "query_partner": qmeta["partner"],
            "query_auth_chain": qmeta["auth_chain"],
            "target_auth_chain": tmeta["auth_chain"],
            **interface,
        }

        row["passes_original_thresholds"] = (
            identity >= 30.0
            and evalue <= 1e-3
            and qcov >= 0.70
            and tcov >= 0.70
        )

        row["passes_identity_evalue_only"] = (
            identity >= 30.0 and evalue <= 1e-3
        )

        row["passes_exploratory_coverage"] = (
            identity >= 20.0
            and evalue <= 1e-3
            and qcov >= 0.30
            and tcov >= 0.30
        )

        row["overlaps_interface_both_sides"] = (
            interface["query_interface_positions_aligned"] >= 1
            and interface["target_interface_positions_aligned"] >= 1
        )

        row["paired_interface_column_present"] = (
            interface["paired_interface_alignment_columns"] >= 1
        )

        enriched.append(row)

    enriched.sort(
        key=lambda r: (
            -float(r["bits"]),
            float(r["evalue_numeric"]),
        )
    )

    output_fields = (
        FIELDS[:-2]
        + [
            "qaln", "taln",
            "identity_percent_normalized",
            "query_coverage_normalized",
            "target_coverage_normalized",
            "evalue_numeric",
            "query_pdb", "query_family", "query_partner",
            "query_auth_chain", "target_auth_chain",
            "query_interface_positions_aligned",
            "target_interface_positions_aligned",
            "paired_interface_alignment_columns",
            "query_interface_positions_json",
            "target_interface_positions_json",
            "paired_interface_columns_json",
            "passes_original_thresholds",
            "passes_identity_evalue_only",
            "passes_exploratory_coverage",
            "overlaps_interface_both_sides",
            "paired_interface_column_present",
        ]
    )

    output_path = OUT / "mmseqs_exhaustive_audit.tsv"

    with output_path.open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=output_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(enriched)

    strict = [
        r for r in enriched
        if r["passes_original_thresholds"]
    ]

    identity_only = [
        r for r in enriched
        if r["passes_identity_evalue_only"]
    ]

    exploratory = [
        r for r in enriched
        if r["passes_exploratory_coverage"]
    ]

    interface_hits = [
        r for r in enriched
        if r["passes_original_thresholds"]
        and r["overlaps_interface_both_sides"]
    ]

    summary = {
        "query_windows": len(query_map),
        "target_windows": len(target_map),
        "all_pairs_possible": len(query_map) * len(target_map),
        "exhaustive_alignments_returned": len(enriched),

        "original_thresholds": {
            "identity_percent": 30.0,
            "evalue_max": 1e-3,
            "query_coverage_min": 0.70,
            "target_coverage_min": 0.70,
            "hits": len(strict),
        },

        "identity_and_evalue_only": {
            "identity_percent": 30.0,
            "evalue_max": 1e-3,
            "hits": len(identity_only),
        },

        "exploratory_relaxed_screen": {
            "identity_percent": 20.0,
            "evalue_max": 1e-3,
            "query_coverage_min": 0.30,
            "target_coverage_min": 0.30,
            "hits": len(exploratory),
        },

        "strict_hits_overlapping_interface_both_sides": len(
            interface_hits
        ),

        "strict_hits_with_paired_interface_alignment_column": sum(
            r["paired_interface_column_present"]
            for r in strict
        ),

        "unique_complexes_in_exploratory_hits": len({
            r["query_pdb"] for r in exploratory
        }),

        "unique_if30_groups_in_exploratory_hits": len({
            r["query_family"] for r in exploratory
        }),

        "outputs": [
            str(RAW),
            str(output_path),
        ],
    }

    summary_path = OUT / "mmseqs_exhaustive_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nEXHAUSTIVE SEARCH AUDIT COMPLETE")
    print(json.dumps(summary, indent=2))

    print("\nTOP 15 ALIGNMENTS")
    for r in enriched[:15]:
        print(
            r["query"], "->", r["target"],
            "identity=", round(
                r["identity_percent_normalized"], 2
            ),
            "qcov=", round(
                r["query_coverage_normalized"], 3
            ),
            "tcov=", round(
                r["target_coverage_normalized"], 3
            ),
            "E=", r["evalue"],
            "bits=", r["bits"],
            "interface_overlap=",
            r["overlaps_interface_both_sides"],
        )


if __name__ == "__main__":
    main()
