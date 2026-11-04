import csv
import hashlib
import json
import urllib.request
from collections import defaultdict, Counter
from datetime import datetime, timezone
from pathlib import Path

OUT = Path("prototype4/results/data_isolation_20261101/cath")
DATA = Path("prototype4/data/cath_latest_release")

BASE = (
    "https://download.cathdb.info/cath/releases/"
    "latest-release/cath-classification-data/"
)

FILES = {
    "domain_list": "cath-domain-list.txt",
    "boundaries": "cath-domain-boundaries.txt",
}

INTERFACE_FILE = Path(
    "prototype4/results/data_isolation_20261101/"
    "interface_extraction/interface_residues.tsv"
)

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(name, filename):
    path = DATA / filename
    DATA.mkdir(parents=True, exist_ok=True)

    if path.exists() and path.stat().st_size > 0:
        print("Using cached:", path)
        return path, "cached"

    url = BASE + filename
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Silicon-Virovore/1.0"},
    )

    print("Downloading:", url)

    with urllib.request.urlopen(request, timeout=180) as response:
        last_modified = response.headers.get("Last-Modified")
        temporary = path.with_suffix(path.suffix + ".partial")

        with temporary.open("wb") as f:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)

        if temporary.stat().st_size == 0:
            raise RuntimeError(f"Empty download: {filename}")

        temporary.replace(path)

    return path, last_modified or "not-provided"


def parse_domain_list(path):
    classifications = {}

    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 12:
                continue

            domain_id = fields[0]

            # CATH levels: Class.Architecture.Topology.Homologous SF
            cath_h = ".".join(fields[1:5])

            classifications[domain_id] = {
                "domain_id": domain_id,
                "cath_class": fields[1],
                "cath_architecture": fields[2],
                "cath_topology": fields[3],
                "cath_homologous_superfamily": cath_h,
                "s35_cluster": fields[5],
                "domain_length": fields[10],
                "resolution_A": fields[11],
            }

    return classifications


def parse_boundaries(path):
    """
    Parse CATH Domall format.

    Each domain contains a number of continuous segments.
    Fragment records after the domain segments are not used.
    """
    domains = []

    with path.open(encoding="utf-8", errors="replace") as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3:
                continue

            chain_key = fields[0]
            domain_field = fields[1]
            fragment_field = fields[2]

            if not domain_field.startswith("D"):
                continue

            try:
                number_domains = int(domain_field[1:])
                number_fragments = int(fragment_field[1:])
            except ValueError:
                continue

            index = 3

            for domain_number in range(1, number_domains + 1):
                if index >= len(fields):
                    raise ValueError(
                        f"Truncated boundary line {line_number}"
                    )

                number_segments = int(fields[index])
                index += 1

                segments = []

                for _ in range(number_segments):
                    if index + 6 > len(fields):
                        raise ValueError(
                            f"Truncated segment on line {line_number}"
                        )

                    start_chain = fields[index]
                    start_residue = fields[index + 1]
                    start_insert = fields[index + 2]
                    end_chain = fields[index + 3]
                    end_residue = fields[index + 4]
                    end_insert = fields[index + 5]

                    index += 6

                    segments.append({
                        "start_chain": start_chain,
                        "start_residue": start_residue,
                        "start_insert": start_insert,
                        "end_chain": end_chain,
                        "end_residue": end_residue,
                        "end_insert": end_insert,
                    })

                domain_id = chain_key + f"{domain_number:02d}"

                domains.append({
                    "domain_id": domain_id,
                    "chain_key": chain_key,
                    "pdb_id": chain_key[:4].lower(),
                    "chain_character": chain_key[4:5],
                    "segments": segments,
                    "source_line": line_number,
                    "declared_fragment_count": number_fragments,
                })

    return domains


def normalize_insertion(value):
    if value in {"", "-", ".", "?", " "}:
        return ""
    return value.strip()


def residue_position(number, insertion):
    """
    Position ordering for author residue number + insertion code.
    """
    return (int(number), normalize_insertion(insertion))


def segment_contains(segment, chain, number, insertion):
    start_chain = segment["start_chain"]
    end_chain = segment["end_chain"]

    if start_chain != chain or end_chain != chain:
        return False

    try:
        position = residue_position(number, insertion)
        start = residue_position(
            segment["start_residue"],
            segment["start_insert"],
        )
        end = residue_position(
            segment["end_residue"],
            segment["end_insert"],
        )
    except ValueError:
        return False

    return start <= position <= end


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not INTERFACE_FILE.exists():
        raise FileNotFoundError(INTERFACE_FILE)

    domain_path, domain_modified = download_file(
        "domain_list", FILES["domain_list"]
    )
    boundary_path, boundary_modified = download_file(
        "boundaries", FILES["boundaries"]
    )

    classifications = parse_domain_list(domain_path)
    domains = parse_boundaries(boundary_path)

    print("\nCATH DATA LOADED")
    print("Classified domains:", len(classifications))
    print("Boundary entries:", len(domains))

    domain_index = defaultdict(list)

    for domain in domains:
        key = (
            domain["pdb_id"],
            domain["chain_character"],
        )
        domain_index[key].append(domain)

    with INTERFACE_FILE.open(
        newline="", encoding="utf-8"
    ) as f:
        interface_rows = list(csv.DictReader(f, delimiter="\t"))

    residue_results = []
    partner_groups = defaultdict(list)

    for row in interface_rows:
        pdb = row["pdb_id"].strip().lower()
        chain = row["auth_chain"].strip()

        residue_number = row["auth_seq_id"].strip()
        insertion = normalize_insertion(
            row["insertion_code"]
        )

        candidates = domain_index.get((pdb, chain), [])

        matching_domains = []

        for domain in candidates:
            if any(
                segment_contains(
                    segment,
                    chain,
                    residue_number,
                    insertion,
                )
                for segment in domain["segments"]
            ):
                classification = classifications.get(
                    domain["domain_id"]
                )

                matching_domains.append({
                    "domain_id": domain["domain_id"],
                    "classification": classification,
                })

        assigned_h = sorted({
            item["classification"]["cath_homologous_superfamily"]
            for item in matching_domains
            if item["classification"] is not None
        })

        if assigned_h:
            status = "assigned_cath_H"
        elif matching_domains:
            status = "domain_boundary_without_classification"
        elif not candidates:
            status = "chain_absent_from_cath_boundary_file"
        else:
            status = "interface_residue_outside_classified_domains"

        result = {
            "pdb_id": row["pdb_id"],
            "family_id_IF30": row["family_id"],
            "partner": row["partner"],
            "auth_chain": chain,
            "auth_seq_id": residue_number,
            "insertion_code": insertion,
            "resname": row["resname"],
            "matching_cath_domains_json": json.dumps([
                item["domain_id"] for item in matching_domains
            ]),
            "cath_H_superfamilies_json": json.dumps(assigned_h),
            "mapping_status": status,
        }

        residue_results.append(result)

        group_key = (
            row["pdb_id"],
            row["family_id"],
            row["partner"],
        )
        partner_groups[group_key].append(result)

    partner_results = []

    for (pdb, if30, partner), rows in sorted(partner_groups.items()):
        assigned = [
            row for row in rows
            if row["mapping_status"] == "assigned_cath_H"
        ]

        all_h = sorted({
            h
            for row in assigned
            for h in json.loads(row["cath_H_superfamilies_json"])
        })

        coverage = len(assigned) / len(rows) if rows else 0.0

        if not assigned:
            status = "no_CATH_H_assignment"
        elif coverage < 0.90:
            status = "partial_CATH_coverage"
        elif len(all_h) == 1:
            status = "single_CATH_H_superfamily"
        else:
            status = "multiple_CATH_H_superfamilies"

        partner_results.append({
            "pdb_id": pdb,
            "family_id_IF30": if30,
            "partner": partner,
            "interface_residue_count": len(rows),
            "CATH_H_assigned_residue_count": len(assigned),
            "CATH_H_coverage": round(coverage, 5),
            "CATH_H_superfamilies_json": json.dumps(all_h),
            "CATH_assignment_status": status,
        })

    def write_tsv(path, rows):
        if not rows:
            print("No rows for:", path)
            return
        with path.open(
            "w", newline="", encoding="utf-8"
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(rows[0].keys()),
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)

    residue_output = OUT / "interface_residue_cath_mapping.tsv"
    partner_output = OUT / "partner_interface_cath_assignments.tsv"

    write_tsv(residue_output, residue_results)
    write_tsv(partner_output, partner_results)

    status_counts = Counter(
        row["CATH_assignment_status"]
        for row in partner_results
    )

    residue_status_counts = Counter(
        row["mapping_status"]
        for row in residue_results
    )

    summary = {
        "retrieval_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "sources": {
            name: {
                "url": BASE + filename,
                "path": str(path),
                "last_modified": modified,
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for name, filename, path, modified in (
                (
                    "domain_list",
                    FILES["domain_list"],
                    domain_path,
                    domain_modified,
                ),
                (
                    "boundaries",
                    FILES["boundaries"],
                    boundary_path,
                    boundary_modified,
                ),
            )
        },
        "interface_residue_records": len(residue_results),
        "partner_groups": len(partner_results),
        "CATH_partner_assignment_status_counts": dict(status_counts),
        "residue_mapping_status_counts": dict(
            residue_status_counts
        ),
        "unique_CATH_H_superfamilies_observed": sorted({
            h
            for row in partner_results
            for h in json.loads(
                row["CATH_H_superfamilies_json"]
            )
        }),
        "outputs": [
            str(residue_output),
            str(partner_output),
        ],
        "interpretation_note": (
            "CATH assignments are domain-level classifications. "
            "Partner interfaces spanning multiple domains or with "
            "incomplete classification coverage are not forced into "
            "one superfamily."
        ),
    }

    (OUT / "cath_mapping_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nCATH INTERFACE MAPPING COMPLETE")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
