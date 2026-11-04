from pathlib import Path
import ast
import hashlib

SOURCE = Path("prototype4/p4.py")

EXPECTED_SHA256 = (
    "4ccb1f37b4be33b23636eab45584021b1df491178f7ac23b62f2609b7cca1181"
)

KEYWORDS = (
    "random_repetitions",
    "GroupKFold",
    "KFold",
    "ShuffleSplit",
    "train_test_split",
    "Ridge",
    "StandardScaler",
    "loco_predictions",
    "outer_group_metrics",
    "bootstrap",
    "predict",
    "feature",
    "family_id",
    "family",
    "split",
    "random",
    "validation_manifest",
    "calibrat",
    "uncertainty",
)

def main():
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)

    raw = SOURCE.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()

    print("SOURCE:", SOURCE)
    print("SHA256:", digest)

    if digest != EXPECTED_SHA256:
        raise SystemExit(
            "STOP: source hash differs from validation manifest."
        )

    source = raw.decode("utf-8")
    lines = source.splitlines()
    tree = ast.parse(source)

    print("\n=== SOURCE HASH VERIFIED ===")

    print("\n=== IMPORTS ===")
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            print(f"{node.lineno:4}: {ast.get_source_segment(source, node)}")

    print("\n=== TOP-LEVEL FUNCTION DEFINITIONS ===")
    functions = []

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            segment = ast.get_source_segment(source, node) or ""
            lower = segment.lower()
            score = sum(
                keyword.lower() in lower
                for keyword in KEYWORDS
            )

            functions.append((score, node, segment))

            print(
                f"{node.lineno:4}-{node.end_lineno:4} "
                f"{node.name:42s} relevance={score}"
            )

    print("\n=== RELEVANT FUNCTION IMPLEMENTATIONS ===")

    selected = [
        item for item in functions
        if item[0] >= 2
    ]

    selected.sort(
        key=lambda item: item[1].lineno
    )

    printed = 0
    MAX_LINES = 340

    for score, node, segment in selected:
        block = segment.splitlines()

        if printed + len(block) > MAX_LINES:
            print(
                "\nSOURCE OUTPUT LIMIT REACHED. "
                "Remaining functions listed above."
            )
            break

        print(
            f"\n--- {node.name} "
            f"(lines {node.lineno}-{node.end_lineno}) "
            f"relevance={score} ---"
        )

        for offset, line in enumerate(block):
            print(f"{node.lineno + offset:4}: {line}")

        printed += len(block)

    print("\n=== SPLIT / MODEL REFERENCE SEARCH ===")

    hits = []

    for i, line in enumerate(lines):
        low = line.lower()

        if any(keyword.lower() in low for keyword in KEYWORDS):
            start = max(1, i + 1 - 2)
            end = min(len(lines), i + 1 + 3)
            hits.append((start, end))

    hits.sort()

    merged = []

    for start, end in hits:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (
                merged[-1][0],
                max(merged[-1][1], end)
            )
        else:
            merged.append((start, end))

    search_printed = 0

    for start, end in merged:
        if search_printed >= 180:
            print("Additional reference lines omitted.")
            break

        print(f"\n--- lines {start}-{end} ---")

        for number in range(start, end + 1):
            print(f"{number:4}: {lines[number - 1]}")

        search_printed += end - start + 1

    print("\n=== END OF READ-ONLY SOURCE AUDIT ===")
    print("Production code was not modified.")


if __name__ == "__main__":
    main()
