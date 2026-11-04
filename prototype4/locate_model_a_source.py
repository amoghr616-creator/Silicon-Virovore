from pathlib import Path
import hashlib

TARGET = (
    "4ccb1f37b4be33b23636eab45584021b1df491178f7ac23b62f2609b7cca1181"
)

SKIP = {
    ".git", ".venv", "venv", "__pycache__",
    "node_modules", ".pytest_cache"
}

matches = []
candidates = []

for path in Path(".").rglob("*.py"):
    if any(part in SKIP for part in path.parts):
        continue

    try:
        raw = path.read_bytes()
        source = raw.decode("utf-8", errors="replace")
    except OSError:
        continue

    digest = hashlib.sha256(raw).hexdigest()

    if digest == TARGET:
        matches.append(path)

    terms = (
        "loco_predictions",
        "bootstrap_models",
        "GroupKFold",
        "outer_group_metrics",
        "ridge_alpha",
        "experimental_epistasis_kcal_mol",
    )

    count = sum(term in source for term in terms)

    if count >= 3:
        candidates.append((count, path, digest))

print("\n=== EXACT MANIFEST HASH MATCH ===")

if matches:
    for path in matches:
        print("EXACT MATCH:", path)
else:
    print("No exact match found in scanned Python files.")

print("\n=== LIKELY MODEL A SOURCE FILES ===")

for count, path, digest in sorted(
    candidates,
    key=lambda x: (-x[0], str(x[1]))
)[:20]:
    print(
        f"score={count}/6 | {path} | sha256={digest}"
    )

print("\n=== END ===")
