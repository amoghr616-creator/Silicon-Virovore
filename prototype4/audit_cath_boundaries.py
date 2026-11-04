from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import json
import re

BASE = Path("prototype4/data/cath_latest_release")

DOMAIN_LIST = BASE / "cath-domain-list.txt"
BOUNDARIES = BASE / "cath-domain-boundaries.txt"

OUTPUT = Path("prototype4/results/cath_boundary_domain_map.tsv")

TARGETS = {
    "1ao7A", "1ao7E",
    "1bj1H", "1danH",
    "1yy9A", "4cvwA",
    "4j2lA", "5m2oB",
}


def chain_key(value):
    value = value.strip()
    return value[:4].lower() + value[4:5]


def domain_key(value):
    value = value.strip()
    return value[:4].lower() + value[4:]


def read_records(path):
    if not path.exists():
        raise FileNotFoundError(path)

    return [
        (line_number, line.strip())
        for line_number, line in enumerate(
            path.read_text(
                encoding="utf-8",
                errors="replace"
            ).splitlines(),
            start=1
        )
        if line.strip()
        and not line.lstrip().startswith("#")
    ]


def read_domain_list():
    records = defaultdict(dict)

    for line_number, line in read_records(DOMAIN_LIST):
        fields = line.split()

        if len(fields) < 12:
            continue

        domain_id = domain_key(fields[0])
        chain = chain_key(domain_id[:5])

        if chain not in TARGETS:
            continue

        if domain_id in records[chain]:
            raise ValueError(
                f"Duplicate domain-list record: {domain_id}"
            )

        records[chain][domain_id] = {
            "cath_code": ".".join(fields[1:5]),
            "s35_cluster": fields[5],
            "length": fields[10],
            "resolution": fields[11],
            "source_line": line_number,
        }

    return records


def read_endpoint(tokens, pos):
    if pos + 3 > len(tokens):
        raise ValueError("Incomplete residue endpoint")

    return {
        "chain": tokens[pos],
        "residue": tokens[pos + 1],
        "insertion": tokens[pos + 2],
    }, pos + 3


def read_segment(tokens, pos):
    start, pos = read_endpoint(tokens, pos)
    end, pos = read_endpoint(tokens, pos)

    return {
        "start": start,
        "end": end,
    }, pos


def parse_boundary(line):
    fields = line.split()

    if len(fields) < 3:
        raise ValueError("Boundary record has fewer than 3 fields")

    chain = chain_key(fields[0])

    d_match = re.fullmatch(r"D(\d+)", fields[1])
    f_match = re.fullmatch(r"F(\d+)", fields[2])

    if not d_match or not f_match:
        raise ValueError(
            f"Invalid domain/fragment counts: {fields[:3]}"
        )

    n_domains = int(d_match.group(1))
    n_fragments = int(f_match.group(1))

    pos = 3
    domains = []

    for ordinal in range(1, n_domains + 1):
        if pos >= len(fields):
            raise ValueError("Missing domain segment count")

        n_segments = int(fields[pos])
        pos += 1

        if n_segments < 1:
            raise ValueError("Domain has no segments")

        segments = []

        for _ in range(n_segments):
            segment, pos = read_segment(fields, pos)
            segments.append(segment)

        # CATH convention: a single-domain chain has suffix 00.
        # Multi-domain chains use 01, 02, ... in boundary order.
        suffix = "00" if n_domains == 1 else f"{ordinal:02d}"

        domains.append({
            "domain_id": f"{chain}{suffix}",
            "ordinal": ordinal,
            "segments": segments,
        })

    fragments = []

    for fragment_number in range(1, n_fragments + 1):
        segment, pos = read_segment(fields, pos)

        if pos >= len(fields):
            raise ValueError("Fragment is missing residue count")

        count_match = re.fullmatch(
            r"\((\d+)\)",
            fields[pos]
        )

        if not count_match:
            raise ValueError(
                f"Invalid fragment residue count: {fields[pos]}"
            )

        fragments.append({
            "fragment_number": fragment_number,
            "segment": segment,
            "residue_count": int(count_match.group(1)),
        })

        pos += 1

    if pos != len(fields):
        raise ValueError(
            f"Unexpected trailing fields: {fields[pos:]}"
        )

    return {
        "chain": chain,
        "declared_domains": n_domains,
        "declared_fragments": n_fragments,
        "domains": domains,
        "fragments": fragments,
    }


def endpoint_text(endpoint):
    insertion = (
        "" if endpoint["insertion"] == "-"
        else endpoint["insertion"]
    )

    return (
        f"{endpoint['chain']}:"
        f"{endpoint['residue']}{insertion}"
    )


def segment_text(segment):
    return (
        endpoint_text(segment["start"])
        + "-"
        + endpoint_text(segment["end"])
    )


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024), b""
        ):
            digest.update(block)

    return digest.hexdigest()


def main():
    classifications = read_domain_list()

    boundaries = {}
    errors = []

    for line_number, line in read_records(BOUNDARIES):
        raw_chain = line.split()[0]
        chain = chain_key(raw_chain)

        if chain not in TARGETS:
            continue

        try:
            parsed = parse_boundary(line)

            if chain in boundaries:
                raise ValueError(
                    f"Duplicate boundary record for {chain}"
                )

            boundaries[chain] = parsed

        except Exception as exc:
            errors.append(
                f"{chain} line {line_number}: {exc}"
            )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "chain_id",
        "domain_id",
        "domain_ordinal",
        "cath_code",
        "s35_cluster",
        "domain_list_length",
        "resolution",
        "segments",
        "segment_count",
        "chain_fragment_count",
        "fragments",
        "classification_status",
    ]

    output_rows = []

    print("\n=== CATH BOUNDARY AUDIT ===")

    for chain in sorted(TARGETS):
        parsed = boundaries.get(chain)
        classified = classifications.get(chain, {})

        print(f"\n{chain}")

        if parsed is None:
            print("  NO VALID BOUNDARY RECORD")
            continue

        print(
            f"  domains={parsed['declared_domains']} "
            f"fragments={parsed['declared_fragments']}"
        )

        expected_ids = set()

        for domain in parsed["domains"]:
            domain_id = domain["domain_id"]
            expected_ids.add(domain_id)

            metadata = classified.get(domain_id)

            segments = domain["segments"]

            ranges = [
                segment_text(segment)
                for segment in segments
            ]

            status = (
                "CLASSIFIED"
                if metadata
                else "BOUNDARY_PRESENT_NO_DOMAIN_LIST_MATCH"
            )

            print(
                f"  {domain_id} | "
                f"segments={len(segments)} | "
                f"ranges={'; '.join(ranges)}"
            )

            if metadata:
                print(
                    f"    CATH={metadata['cath_code']} "
                    f"S35={metadata['s35_cluster']} "
                    f"length={metadata['length']}"
                )
            else:
                print(
                    "    NO EXACT CLASSIFICATION RECORD "
                    "IN THE INSPECTED DOMAIN LIST"
                )

            output_rows.append({
                "chain_id": chain,
                "domain_id": domain_id,
                "domain_ordinal": domain["ordinal"],
                "cath_code": (
                    metadata["cath_code"] if metadata else ""
                ),
                "s35_cluster": (
                    metadata["s35_cluster"] if metadata else ""
                ),
                "domain_list_length": (
                    metadata["length"] if metadata else ""
                ),
                "resolution": (
                    metadata["resolution"] if metadata else ""
                ),
                "segments": json.dumps(
                    segments, separators=(",", ":")
                ),
                "segment_count": len(segments),
                "chain_fragment_count": parsed[
                    "declared_fragments"
                ],
                "fragments": json.dumps(
                    parsed["fragments"], separators=(",", ":")
                ),
                "classification_status": status,
            })

        extra_ids = set(classified) - expected_ids

        if extra_ids:
            print(
                "  DOMAIN-LIST IDS NOT IN BOUNDARY RECORD:",
                sorted(extra_ids)
            )

        if parsed["fragments"]:
            for fragment in parsed["fragments"]:
                print(
                    "  EXCLUDED FRAGMENT "
                    f"{fragment['fragment_number']}: "
                    f"{segment_text(fragment['segment'])} "
                    f"(n={fragment['residue_count']})"
                )

    with OUTPUT.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print("\n=== FORMAT CHECKS ===")

    if "1ao7A" in boundaries:
        parsed = boundaries["1ao7A"]

        assert len(parsed["domains"]) == 2
        assert [
            len(d["segments"]) for d in parsed["domains"]
        ] == [1, 1]

        assert segment_text(
            parsed["domains"][0]["segments"][0]
        ) == "A:1-A:181"

        assert segment_text(
            parsed["domains"][1]["segments"][0]
        ) == "A:182-A:270"

        assert len(parsed["fragments"]) == 1

        print("PASS: 1ao7A two-domain and fragment parsing")

    if "1danH" in boundaries:
        parsed = boundaries["1danH"]

        assert len(parsed["domains"]) == 2
        assert [
            len(d["segments"]) for d in parsed["domains"]
        ] == [2, 2]

        assert len(parsed["fragments"]) == 1

        print("PASS: 1danH discontinuous-domain parsing")

    if "1yy9A" in boundaries:
        parsed = boundaries["1yy9A"]

        assert len(parsed["domains"]) == 4
        assert all(
            len(d["segments"]) == 1
            for d in parsed["domains"]
        )

        print("PASS: 1yy9A four-domain parsing")

    print("\n=== SOURCE HASHES ===")
    print("domain_list_sha256:", sha256(DOMAIN_LIST))
    print("boundaries_sha256:", sha256(BOUNDARIES))

    print("\n=== PARSE ERRORS ===")

    if errors:
        for error in errors:
            print("ERROR:", error)
    else:
        print("None")

    print("\n=== OUTPUT ===")
    print(OUTPUT)
    print("Rows:", len(output_rows))

    if errors:
        raise SystemExit(
            "\nFAILED: resolve parser errors before using this mapping."
        )


if __name__ == "__main__":
    main()
