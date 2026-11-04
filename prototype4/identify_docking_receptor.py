from pathlib import Path
import hashlib

PDB = Path("results/docking/model.003.00.pdb")
REFERENCES = {
    "P61566_ER​​VK24": Path(
        "prototype4/results/official_references/P61566.fasta"
    ),
    "P61570_ER​​VK25": Path(
        "prototype4/results/official_references/P61570.fasta"
    ),
}

AA3 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D",
    "CYS":"C", "GLN":"Q", "GLU":"E", "GLY":"G",
    "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K",
    "MET":"M", "PHE":"F", "PRO":"P", "SER":"S",
    "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
    "MSE":"M",
}


def parse_fasta(path):
    if not path.is_file():
        raise FileNotFoundError(path)

    return "".join(
        line.strip().upper()
        for line in path.read_text().splitlines()
        if line and not line.startswith(">")
    )


def parse_pdb_chain(path, wanted_chain="A"):
    residues = {}

    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            if len(line) < 27:
                continue
            if line[21].strip() != wanted_chain:
                continue
            if line[16] not in (" ", "A"):
                continue

            aa = AA3.get(line[17:20].strip().upper())
            if aa is None:
                continue

            key = (line[22:26].strip(), line[26].strip())
            residues.setdefault(key, aa)

    return "".join(residues.values())


def needleman_wunsch(a, b):
    """Simple global alignment with linear gap penalties."""
    n, m = len(a), len(b)
    gap = -2

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap
        trace[i][0] = "U"

    for j in range(1, m + 1):
        score[0][j] = j * gap
        trace[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            choices = [
                (
                    score[i - 1][j - 1]
                    + (2 if a[i - 1] == b[j - 1] else -1),
                    "D",
                ),
                (score[i - 1][j] + gap, "U"),
                (score[i][j - 1] + gap, "L"),
            ]
            score[i][j], trace[i][j] = max(choices)

    i, j = n, m
    matches = 0
    paired = 0
    gaps = 0
    mismatches = []

    while i > 0 or j > 0:
        direction = trace[i][j]

        if direction == "D":
            paired += 1
            if a[i - 1] == b[j - 1]:
                matches += 1
            else:
                mismatches.append((i, a[i - 1], j, b[j - 1]))
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

    return {
        "score": score[n][m],
        "matches": matches,
        "paired": paired,
        "identity": matches / paired if paired else 0,
        "model_coverage": paired / n if n else 0,
        "reference_coverage": paired / m if m else 0,
        "gaps": gaps,
        "mismatches": list(reversed(mismatches)),
    }


def main():
    if not PDB.is_file():
        raise FileNotFoundError(PDB)

    model = parse_pdb_chain(PDB)

    print("========== DOCKING RECEPTOR ==========")
    print(f"Sequence length: {len(model)}")
    print(
        "SHA256:",
        hashlib.sha256(PDB.read_bytes()).hexdigest(),
    )

    results = {}

    for label, path in REFERENCES.items():
        reference = parse_fasta(path)

        result = needleman_wunsch(model, reference)
        results[label] = result

        print(f"\n========== {label} ==========")
        print(f"Reference length: {len(reference)}")
        print(f"Alignment score: {result['score']}")
        print(
            f"Identity: {result['identity']:.3%} "
            f"({result['matches']}/{result['paired']})"
        )
        print(f"Docking-chain coverage: {result['model_coverage']:.3%}")
        print(f"Reference coverage: {result['reference_coverage']:.3%}")
        print(f"Alignment gaps: {result['gaps']}")
        print(f"Substitutions: {len(result['mismatches'])}")

        print("First 20 aligned substitutions:")
        for item in result["mismatches"][:20]:
            model_pos, model_aa, ref_pos, ref_aa = item
            print(
                f"  Model {model_pos}{model_aa} "
                f"vs reference {ref_pos}{ref_aa}"
            )

    print("\n========== COMPARISON SUMMARY ==========")
    for label, result in results.items():
        print(
            f"{label}: identity={result['identity']:.3%}, "
            f"model_coverage={result['model_coverage']:.3%}, "
            f"reference_coverage={result['reference_coverage']:.3%}"
        )

    print("\nNo structural files or manifests were modified.")


if __name__ == "__main__":
    main()
