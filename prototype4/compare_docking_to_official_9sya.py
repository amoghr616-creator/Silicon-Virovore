from pathlib import Path
from collections import OrderedDict
import hashlib

REFERENCE = Path(
    "prototype4/results/official_references/9SYA_reference.pdb"
)
DOCKING = Path("results/docking/model.003.00.pdb")

AA3 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D",
    "CYS":"C", "GLN":"Q", "GLU":"E", "GLY":"G",
    "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K",
    "MET":"M", "PHE":"F", "PRO":"P", "SER":"S",
    "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
    "MSE":"M",
}

HOTSPOTS = (141, 167, 185, 203, 205)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_chains(path):
    chains = OrderedDict()

    with path.open(errors="replace") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break

            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            if len(line) < 27 or line[16] not in (" ", "A"):
                continue

            residue = line[17:20].strip().upper()
            if residue not in AA3:
                continue

            chain = line[21].strip() or "_"
            residue_id = (
                line[22:26].strip(),
                line[26].strip(),
            )

            chains.setdefault(chain, OrderedDict())
            chains[chain].setdefault(residue_id, AA3[residue])

    return {
        chain: {
            "sequence": "".join(residues.values()),
            "residue_ids": list(residues.keys()),
        }
        for chain, residues in chains.items()
    }


def semiglobal_align(reference, target):
    """Align all of reference to its best-matching target region."""
    n, m = len(reference), len(target)
    gap = -2

    dp = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i * gap
        trace[i][0] = "U"

    for j in range(1, m + 1):
        trace[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            choices = [
                (
                    dp[i - 1][j - 1]
                    + (2 if reference[i - 1] == target[j - 1] else -1),
                    "D",
                ),
                (dp[i - 1][j] + gap, "U"),
                (dp[i][j - 1] + gap, "L"),
            ]
            dp[i][j], trace[i][j] = max(choices)

    end_j = max(range(m + 1), key=lambda j: dp[n][j])

    i, j = n, end_j
    mapping = {}
    paired = []
    gaps = 0

    while i > 0:
        direction = trace[i][j]

        if direction == "D":
            mapping[i - 1] = j - 1
            paired.append((reference[i - 1], target[j - 1]))
            i -= 1
            j -= 1

        elif direction == "U":
            gaps += 1
            i -= 1

        elif direction == "L":
            gaps += 1
            j -= 1

        else:
            raise RuntimeError("Invalid alignment traceback.")

    matches = sum(a == b for a, b in paired)

    return {
        "mapping": mapping,
        "identity": matches / len(paired) if paired else 0.0,
        "matches": matches,
        "paired": len(paired),
        "reference_coverage": len(paired) / n if n else 0.0,
        "gaps": gaps,
        "score": dp[n][end_j],
    }


def main():
    if not REFERENCE.is_file():
        raise FileNotFoundError(REFERENCE)

    if not DOCKING.is_file():
        raise FileNotFoundError(DOCKING)

    official = parse_chains(REFERENCE)
    modeled = parse_chains(DOCKING)

    if "A" not in modeled:
        raise RuntimeError("Docking file has no chain A.")

    target_sequence = modeled["A"]["sequence"]
    target_ids = modeled["A"]["residue_ids"]

    print("========== FILE PROVENANCE ==========")
    print(f"Official 9SYA SHA256: {sha256(REFERENCE)}")
    print(f"Docking model SHA256: {sha256(DOCKING)}")

    print("\n========== CHAIN INVENTORY ==========")
    for chain, data in official.items():
        print(
            f"Official chain {chain}: "
            f"{len(data['sequence'])} observed residues"
        )

    for chain, data in modeled.items():
        print(
            f"Docking chain {chain}: "
            f"{len(data['sequence'])} observed residues"
        )

    print("\n========== OFFICIAL 9SYA ALIGNMENT ==========")

    for chain, data in official.items():
        result = semiglobal_align(
            data["sequence"], target_sequence
        )

        print(f"\nOfficial chain: {chain}")
        print(f"Identity: {result['identity']:.3%}")
        print(
            f"Matched positions: "
            f"{result['matches']}/{result['paired']}"
        )
        print(
            f"Reference coverage: "
            f"{result['reference_coverage']:.3%}"
        )
        print(f"Alignment gaps: {result['gaps']}")
        print(f"Alignment score: {result['score']}")

        print("Hotspot residue mapping:")

        for hotspot in HOTSPOTS:
            indices = [
                i for i, residue_id
                in enumerate(data["residue_ids"])
                if residue_id[0] == str(hotspot)
            ]

            if not indices:
                print(f"  9SYA {hotspot}: not observed in coordinates")
                continue

            for ref_index in indices:
                model_index = result["mapping"].get(ref_index)

                if model_index is None:
                    print(
                        f"  9SYA {hotspot}: "
                        "not aligned to docking chain"
                    )
                    continue

                ref_aa = data["sequence"][ref_index]
                model_aa = target_sequence[model_index]
                model_id = target_ids[model_index]

                print(
                    f"  9SYA {hotspot}{ref_aa} -> "
                    f"docking A:{model_id[0]}{model_id[1]}{model_aa}"
                )

    print("\n========== END OF COMPARISON ==========")
    print("No original structure files were modified.")


if __name__ == "__main__":
    main()
