from __future__ import annotations

import csv
import tarfile
from pathlib import Path

import numpy as np

ROOT = Path.cwd()
REFERENCE = ROOT / "prototype4/results/official_references/9SYA_reference.pdb"
SOURCE = ROOT / "prototype4/data/haddock3/receptors/9sya_source_chain_A.pdb"
RUN_ROOT = ROOT / "prototype4/results/p4_reliability_repeats_20261103"
OUT = ROOT / "prototype4/results/p4_target_mapping_20261103"

RUNS = {
    "1917": RUN_ROOT / "HADDOCK3_chain_A_seed1917",
    "2917": RUN_ROOT / "HADDOCK3_chain_A_seed2917",
}
POSES = [f"cluster_{i}_model_1.pdb" for i in range(1, 6)]

AA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H",
}


def parse_chain(data: bytes, chain_id: str):
    residues = []
    seen = set()

    for line in data.decode("utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 27:
            continue
        if line[21:22].strip() != chain_id:
            continue
        if line[16:17].strip() not in ("", "A"):
            continue

        resname = line[17:20].strip().upper()
        if resname not in AA:
            continue

        resseq = line[22:26].strip()
        icode = line[26:27].strip()
        key = (resseq, icode, resname)

        if key in seen:
            continue
        seen.add(key)

        residues.append({
            "resseq": resseq,
            "icode": icode,
            "resname": resname,
            "aa": AA[resname],
        })

    return residues


def extract_pose(run_dir: Path, filename: str) -> bytes:
    direct = list(run_dir.rglob(filename))
    if direct:
        return direct[0].read_bytes()

    for archive in sorted(run_dir.rglob("summary.tgz")):
        with tarfile.open(archive, "r:*") as tar:
            matches = [
                m for m in tar.getmembers()
                if m.isfile() and Path(m.name).name == filename
            ]
            if matches:
                handle = tar.extractfile(matches[0])
                if handle:
                    return handle.read()

    raise FileNotFoundError(f"{filename} not found in {run_dir}")


def align_sequences(seq_a: str, seq_b: str):
    """Global Needleman-Wunsch alignment; returns A-index -> B-index."""
    n, m = len(seq_a), len(seq_b)
    score = np.zeros((n + 1, m + 1), dtype=np.int32)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)

    for i in range(1, n + 1):
        score[i, 0] = -3 * i
        trace[i, 0] = 1  # up
    for j in range(1, m + 1):
        score[0, j] = -3 * j
        trace[0, j] = 2  # left

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = score[i - 1, j - 1] + (
                3 if seq_a[i - 1] == seq_b[j - 1] else -2
            )
            up = score[i - 1, j] - 3
            left = score[i, j - 1] - 3

            best = max(diagonal, up, left)
            score[i, j] = best
            trace[i, j] = (
                0 if best == diagonal
                else 1 if best == up
                else 2
            )

    mapping = {}
    i, j = n, m
    aligned = matches = 0

    while i > 0 or j > 0:
        move = trace[i, j]

        if i > 0 and j > 0 and move == 0:
            mapping[i - 1] = j - 1
            aligned += 1
            matches += seq_a[i - 1] == seq_b[j - 1]
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or move == 1):
            i -= 1
        else:
            j -= 1

    identity = matches / aligned if aligned else 0.0
    return mapping, identity, aligned


def find_reference_position(residues, position):
    matches = [
        i for i, residue in enumerate(residues)
        if residue["resseq"] == str(position)
        and residue["resname"] == "ARG"
    ]
    return matches


def describe(residue):
    if residue is None:
        return "UNMAPPED"
    return (
        f"{residue['resname']} "
        f"{residue['resseq']}{residue['icode']}"
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not REFERENCE.is_file() or not SOURCE.is_file():
        raise FileNotFoundError(
            f"Missing reference or source:\n{REFERENCE}\n{SOURCE}"
        )

    reference = parse_chain(REFERENCE.read_bytes(), "A")
    source = parse_chain(SOURCE.read_bytes(), "A")

    if not reference or not source:
        raise ValueError(
            f"Could not parse chain A: "
            f"reference={len(reference)}, source={len(source)}"
        )

    ref_seq = "".join(r["aa"] for r in reference)
    source_seq = "".join(r["aa"] for r in source)

    ref_to_source, identity, aligned = align_sequences(
        ref_seq, source_seq
    )

    print("REFERENCE / PREPARED SOURCE")
    print("  Reference chain A length:", len(reference))
    print("  Prepared source chain A length:", len(source))
    print(f"  Alignment identity: {identity:.4f}")
    print("  Aligned residue positions:", aligned)
    print("  Exact full-sequence match:", ref_seq == source_seq)

    results = []

    for position in (203, 205):
        ref_indices = find_reference_position(reference, position)

        print(f"\nREFERENCE ARG{position}")
        if not ref_indices:
            print("  NOT FOUND as ARG at that reference residue number.")
            results.append({
                "seed": "reference",
                "pose": "reference",
                "reference_position": position,
                "reference_residue": "NOT_FOUND",
                "source_residue": "NOT_MAPPED",
                "output_residue": "NOT_CHECKED",
                "sequence_identity_reference_source": identity,
                "sequence_identity_source_output": "",
                "mapping_status": "REFERENCE_POSITION_NOT_FOUND",
            })
            continue

        for ref_idx in ref_indices:
            source_idx = ref_to_source.get(ref_idx)
            source_res = source[source_idx] if source_idx is not None else None
            print(
                f"  Reference: index={ref_idx + 1}, "
                f"{describe(reference[ref_idx])}"
            )
            print(
                f"  Prepared source: "
                f"{describe(source_res)}"
            )

            for seed, run_dir in RUNS.items():
                # Also verify the receptor input supplied to this run.
                run_input = (
                    run_dir / "data/0_topoaa/9sya_source_chain_A.pdb"
                )
                run_input_residues = parse_chain(
                    run_input.read_bytes(), "A"
                )
                run_input_seq = "".join(
                    r["aa"] for r in run_input_residues
                )
                run_input_to_source, input_identity, _ = align_sequences(
                    source_seq, run_input_seq
                )

                for pose_name in POSES:
                    pose_data = extract_pose(run_dir, pose_name)
                    output = parse_chain(pose_data, "A")
                    output_seq = "".join(r["aa"] for r in output)

                    source_to_output, output_identity, _ = (
                        align_sequences(source_seq, output_seq)
                    )

                    output_idx = None
                    if source_idx is not None:
                        output_idx = source_to_output.get(source_idx)

                    output_res = (
                        output[output_idx]
                        if output_idx is not None else None
                    )

                    # Independent direct-numbering check.
                    direct = find_reference_position(output, position)
                    direct_desc = (
                        describe(output[direct[0]]) if direct
                        else "NOT_FOUND"
                    )

                    mapped_ok = (
                        output_res is not None
                        and output_res["resname"] == "ARG"
                    )

                    status = (
                        "MAPPED_ARG"
                        if mapped_ok
                        else "MAPPING_OR_IDENTITY_MISMATCH"
                    )

                    row = {
                        "seed": seed,
                        "pose": pose_name,
                        "reference_position": position,
                        "reference_residue": describe(reference[ref_idx]),
                        "source_residue": describe(source_res),
                        "output_residue": describe(output_res),
                        "direct_output_numbering": direct_desc,
                        "sequence_identity_reference_source": identity,
                        "sequence_identity_run_input_source": input_identity,
                        "sequence_identity_source_output": output_identity,
                        "mapping_status": status,
                    }
                    results.append(row)

                    print(
                        f"  {seed} {pose_name}: "
                        f"mapped={describe(output_res)}, "
                        f"direct-number={direct_desc}, "
                        f"identity(source/output)={output_identity:.4f}, "
                        f"status={status}"
                    )

    path = OUT / "target_residue_sequence_mapping.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(results[0].keys()),
        )
        writer.writeheader()
        writer.writerows(results)

    print("\nSaved:", path)


if __name__ == "__main__":
    main()
