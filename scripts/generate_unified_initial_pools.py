from __future__ import annotations

import hashlib
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.config import DEFAULT_SEED_SEQUENCE

SEEDS = range(631, 641)
AA = "ACDEFGHIKLMNPQRSTVWY"

REFERENCE_MANIFEST = (
    ROOT / "prototype4/results/"
    "real_structure_integration_check_20261101/"
    "historical_pose_manifest.json"
)

OUT = ROOT / "results/unified_ablation_20261102/initial_pools"

seed_sequence = DEFAULT_SEED_SEQUENCE.upper()
if not seed_sequence or set(seed_sequence) - set(AA):
    raise ValueError("DEFAULT_SEED_SEQUENCE is not a canonical peptide.")

manifest = json.loads(REFERENCE_MANIFEST.read_text())
reference_sequences = {
    str(item["sequence"]).strip().upper()
    for item in manifest["complexes"]
}
if len(reference_sequences) != 1:
    raise ValueError(
        "Expected exactly one distinct historical structural sequence; "
        f"found {len(reference_sequences)}."
    )

reference_sequence = next(iter(reference_sequences))

if len(reference_sequence) != len(seed_sequence):
    raise ValueError("Historical structural sequence length differs from seed.")

if set(reference_sequence) - set(AA):
    raise ValueError("Historical structural sequence is not canonical.")

neighbors = []
for position, original in enumerate(seed_sequence):
    for replacement in AA:
        if replacement == original:
            continue
        variant = (
            seed_sequence[:position]
            + replacement
            + seed_sequence[position + 1:]
        )
        if variant != reference_sequence:
            neighbors.append(variant)

OUT.mkdir(parents=True, exist_ok=True)
summary = {
    "design": "shared matched initial sequence pools",
    "sequence_length": len(seed_sequence),
    "pool_size": 100,
    "pool_composition": {
        "configured_seed": 1,
        "historical_structure_reference_sequence": 1,
        "single_substitution_neighbors": 98,
    },
    "historical_structure_manifest": str(
        REFERENCE_MANIFEST.relative_to(ROOT)
    ),
    "historical_reference_sequence": reference_sequence,
    "pools": [],
}

for seed in SEEDS:
    rng = random.Random(seed)
    selected = rng.sample(neighbors, 98)
    pool = [seed_sequence, reference_sequence, *selected]

    if len(pool) != 100 or len(set(pool)) != 100:
        raise RuntimeError(f"Pool {seed} is not 100 unique sequences.")

    fasta = "".join(
        f">seed_{seed}_candidate_{index:03d}\n{sequence}\n"
        for index, sequence in enumerate(pool, start=1)
    )

    path = OUT / f"seed_{seed}.fasta"
    path.write_text(fasta, encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()

    summary["pools"].append({
        "seed": seed,
        "path": str(path.relative_to(ROOT)),
        "sequence_count": len(pool),
        "unique_sequences": len(set(pool)),
        "sha256": digest,
        "reference_sequence_present": reference_sequence in pool,
    })

    print(
        f"Seed {seed}: {len(pool)} unique sequences | "
        f"SHA-256 {digest}"
    )

manifest_out = OUT / "pool_manifest.json"
manifest_out.write_text(json.dumps(summary, indent=2) + "\n")

print(f"\nPool manifest: {manifest_out}")
