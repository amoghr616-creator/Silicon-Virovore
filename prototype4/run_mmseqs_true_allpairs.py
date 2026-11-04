import csv
import itertools
import json
import shutil
import subprocess
from pathlib import Path

OUT = Path(
    "prototype4/results/data_isolation_20261101/"
    "interface_window_mmseqs"
)

CORE = OUT / "true_allpairs"
QUERY_FASTA = OUT / "skempi_interface_windows.fasta"
TARGET_FASTA = OUT / "hervk_9sya_target_windows.fasta"
QUERY_MAP = OUT / "skempi_interface_window_map.tsv"
TARGET_MAP = OUT / "hervk_target_window_map.tsv"

FIELDS = [
    "query", "target", "pident", "alnlen",
    "qstart", "qend", "tstart", "tend",
    "qcov", "tcov", "evalue", "bits",
    "qaln", "taln",
]


def run(command, log):
    with log.open("w", encoding="utf-8") as f:
        subprocess.run(
            command,
            check=True,
            stdout=f,
            stderr=subprocess.STDOUT,
        )


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def norm_cov(value):
    value = float(value)
    return value / 100.0 if value > 1.0 else value


def interface_alignment(hit, qmeta, tmeta):
    q_positions = set(
        json.loads(qmeta["interface_positions_relative_json"])
    )
    t_positions = set(
        json.loads(tmeta["interface_positions_relative_json"])
    )

    qpos = int(hit["qstart"])
    tpos = int(hit["tstart"])

    q_overlap = set()
    t_overlap = set()
    paired = []

    for qa, ta in zip(hit["qaln"], hit["taln"]):
        q_has = qa != "-"
        t_has = ta != "-"

        if q_has and t_has:
            if qpos in q_positions:
                q_overlap.add(qpos)
            if tpos in t_positions:
                t_overlap.add(tpos)

            if qpos in q_positions and tpos in t_positions:
                paired.append((qpos, tpos))

        if q_has:
            qpos += 1
        if t_has:
            tpos += 1

    return {
        "query_interface_positions_aligned": len(q_overlap),
        "target_interface_positions_aligned": len(t_overlap),
        "paired_interface_alignment_columns": len(paired),
        "query_interface_positions_json": json.dumps(sorted(q_overlap)),
        "target_interface_positions_json": json.dumps(sorted(t_overlap)),
        "paired_interface_columns_json": json.dumps(paired),
    }


def main():
    mmseqs = shutil.which("mmseqs")
    if not mmseqs:
        raise RuntimeError("MMseqs2 executable not found")

    for path in (QUERY_FASTA, TARGET_FASTA, QUERY_MAP, TARGET_MAP):
        if not path.exists():
            raise FileNotFoundError(path)

    qrows = read_tsv(QUERY_MAP)
    trows = read_tsv(TARGET_MAP)

    qmap = {r["window_id"]: r for r in qrows}
    tmap = {r["window_id"]: r for r in trows}

    expected_pairs = len(qmap) * len(tmap)

    if expected_pairs != 402:
        raise RuntimeError(
            f"Expected 402 sequence pairs from this audit; "
            f"found {expected_pairs}. Inspect FASTA/maps."
        )

    if CORE.exists():
        shutil.rmtree(CORE)
    CORE.mkdir(parents=True)

    qdb = CORE / "queryDB"
    tdb = CORE / "targetDB"
    pref = CORE / "all_target_candidates"
    result = CORE / "allpairs_alignmentDB"
    tmp = CORE / "tmp"

    # 1. Create MMseqs sequence databases.
    run(
        [mmseqs, "createdb", str(QUERY_FASTA), str(qdb)],
        CORE / "createdb_query.log",
    )
    run(
        [mmseqs, "createdb", str(TARGET_FASTA), str(tdb)],
        CORE / "createdb_target.log",
    )

    # 2. Build a synthetic prefilter:
    # every query entry points to the entire target index.
    qindex = Path(str(qdb) + ".index")
    tindex = Path(str(tdb) + ".index")

    if not qindex.exists() or not tindex.exists():
        raise RuntimeError("MMseqs database index files missing")

    target_index_size = tindex.stat().st_size

    qindex_rows = [
        line.split("\t")
        for line in qindex.read_text().splitlines()
        if line.strip()
    ]

    if len(qindex_rows) != len(qmap):
        raise RuntimeError(
            f"Query index contains {len(qindex_rows)} entries, "
            f"expected {len(qmap)}"
        )

    # The MMseqs user-guide all-pairs method uses the target index
    # as the data file and points each query index entry to it.
    pref.symlink_to(tindex.resolve())

    with Path(str(pref) + ".index").open(
        "w", encoding="utf-8"
    ) as f:
        for fields in qindex_rows:
            f.write(f"{fields[0]}\t0\t{target_index_size}\n")

    Path(str(pref) + ".dbtype").write_bytes(
        bytes([7, 0, 0, 0])
    )

    # 3. Align the full target list against every query.
    # A very permissive E-value and zero identity/coverage filters
    # retain weak alignments for downstream, explicit filtering.
    align_command = [
        mmseqs, "align",
        str(qdb),
        str(tdb),
        str(pref),
        str(result),
        "-e", "1e100",
        "--min-seq-id", "0",
        "--min-aln-len", "0",
        "-c", "0",
        "--cov-mode", "0",
        "--alignment-mode", "3",
        "-a",
        "--threads", "4",
    ]

    print("Running direct MMseqs2 align on the full candidate list.")
    run(align_command, CORE / "align.log")

    # 4. Convert the alignment database to a readable table.
    raw = CORE / "allpairs_alignments.tsv"

    convert_command = [
        mmseqs, "convertalis",
        str(qdb),
        str(tdb),
        str(result),
        str(raw),
        "--format-output", ",".join(FIELDS),
        "--threads", "4",
    ]

    run(convert_command, CORE / "convertalis.log")

    # 5. Attach biological provenance and evaluate interface overlap.
    records = []
    seen_pairs = set()

    with raw.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) != len(FIELDS):
                raise ValueError(
                    f"Expected {len(FIELDS)} columns, got {len(parts)}"
                )

            hit = dict(zip(FIELDS, parts))
            qmeta = qmap[hit["query"]]
            tmeta = tmap[hit["target"]]

            pair = (hit["query"], hit["target"])
            if pair in seen_pairs:
                raise RuntimeError(
                    f"Duplicate query-target pair: {pair}"
                )
            seen_pairs.add(pair)

            identity_pct = float(hit["pident"])
            qcov = norm_cov(hit["qcov"])
            tcov = norm_cov(hit["tcov"])
            evalue = float(hit["evalue"])

            overlap = interface_alignment(hit, qmeta, tmeta)

            row = {
                **hit,
                "identity_percent": identity_pct,
                "query_coverage_fraction": qcov,
                "target_coverage_fraction": tcov,
                "evalue_numeric": evalue,
                "query_pdb": qmeta["pdb_id"],
                "query_family": qmeta["family_id"],
                "query_partner": qmeta["partner"],
                "query_auth_chain": qmeta["auth_chain"],
                "target_auth_chain": tmeta["auth_chain"],
                **overlap,
            }

            row["passes_original_thresholds"] = (
                identity_pct >= 30.0
                and evalue <= 1e-3
                and qcov >= 0.70
                and tcov >= 0.70
            )

            row["passes_exploratory_thresholds"] = (
                identity_pct >= 20.0
                and evalue <= 1e-3
                and qcov >= 0.30
                and tcov >= 0.30
            )

            row["overlaps_interface_both_sides"] = (
                overlap["query_interface_positions_aligned"] >= 1
                and overlap["target_interface_positions_aligned"] >= 1
            )

            row["has_paired_interface_column"] = (
                overlap["paired_interface_alignment_columns"] >= 1
            )

            records.append(row)

    records.sort(
        key=lambda r: (
            -float(r["bits"]),
            float(r["evalue_numeric"]),
        )
    )

    out_fields = FIELDS + [
        "identity_percent",
        "query_coverage_fraction",
        "target_coverage_fraction",
        "evalue_numeric",
        "query_pdb",
        "query_family",
        "query_partner",
        "query_auth_chain",
        "target_auth_chain",
        "query_interface_positions_aligned",
        "target_interface_positions_aligned",
        "paired_interface_alignment_columns",
        "query_interface_positions_json",
        "target_interface_positions_json",
        "paired_interface_columns_json",
        "passes_original_thresholds",
        "passes_exploratory_thresholds",
        "overlaps_interface_both_sides",
        "has_paired_interface_column",
    ]

    output_table = CORE / "allpairs_annotated.tsv"

    with output_table.open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f, fieldnames=out_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(records)

    strict = [
        r for r in records if r["passes_original_thresholds"]
    ]
    exploratory = [
        r for r in records if r["passes_exploratory_thresholds"]
    ]
    interface_hits = [
        r for r in strict if r["overlaps_interface_both_sides"]
    ]

    returned_pairs = {
        (r["query"], r["target"]) for r in records
    }

    missing_pairs = [
        pair
        for pair in itertools.product(qmap, tmap)
        if pair not in returned_pairs
    ]

    summary = {
        "method": (
            "MMseqs2 direct align using synthetic prefilter "
            "containing all target entries for every query"
        ),
        "query_windows": len(qmap),
        "target_windows": len(tmap),
        "expected_query_target_pairs": expected_pairs,
        "alignment_records_returned": len(records),
        "unique_pairs_returned": len(returned_pairs),
        "pairs_without_returned_alignment": len(missing_pairs),
        "missing_pair_examples": missing_pairs[:20],
        "original_thresholds": {
            "identity_percent_min": 30,
            "evalue_max": 1e-3,
            "query_coverage_min": 0.70,
            "target_coverage_min": 0.70,
            "hits": len(strict),
        },
        "exploratory_thresholds": {
            "identity_percent_min": 20,
            "evalue_max": 1e-3,
            "query_coverage_min": 0.30,
            "target_coverage_min": 0.30,
            "hits": len(exploratory),
        },
        "strict_hits_overlapping_interface_both_sides": len(
            interface_hits
        ),
        "strict_hits_with_paired_interface_column": sum(
            r["has_paired_interface_column"] for r in strict
        ),
        "unique_complexes_in_exploratory_hits": len({
            r["query_pdb"] for r in exploratory
        }),
        "unique_if30_groups_in_exploratory_hits": len({
            r["query_family"] for r in exploratory
        }),
        "outputs": [
            str(raw),
            str(output_table),
            str(CORE / "align.log"),
        ],
    }

    summary_path = CORE / "allpairs_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nTRUE ALL-PAIRS ALIGNMENT COMPLETE")
    print(json.dumps(summary, indent=2))

    print("\nTOP 15 ALIGNMENTS")
    for r in records[:15]:
        print(
            r["query"], "->", r["target"],
            "identity=", round(r["identity_percent"], 2),
            "qcov=", round(r["query_coverage_fraction"], 3),
            "tcov=", round(r["target_coverage_fraction"], 3),
            "E=", r["evalue"],
            "bits=", r["bits"],
            "interface_both_sides=",
            r["overlaps_interface_both_sides"],
        )


if __name__ == "__main__":
    main()
