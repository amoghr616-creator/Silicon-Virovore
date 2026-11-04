import json
from pathlib import Path

REPORT = Path(
    "results/p4_live_gate_smoke_20261102/"
    "p4_live_gate_smoke_20261102/"
    "p4_integrated_assessment/generation_0001.json"
)

TERMS = (
    "target", "coverage", "candidate", "eligible", "abstain",
    "pareto", "population", "valid", "invalid", "missing",
    "reason", "status", "receptor", "pose", "complex",
    "source", "sequence_hash", "calibration",
)


def walk(obj, path="root"):
    if isinstance(obj, dict):
        for key, value in obj.items():
            current = f"{path}.{key}"
            key_lower = key.lower()

            if isinstance(value, list):
                if any(term in key_lower for term in TERMS):
                    print(f"\n{current}: list with {len(value)} entries")
                walk(value, current)

            elif isinstance(value, dict):
                if any(term in key_lower for term in TERMS):
                    print(f"\n{current}: object")
                walk(value, current)

            elif any(term in current.lower() for term in TERMS):
                print(f"{current}: {value!r}")

    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            walk(value, f"{path}[{index}]")


def main():
    if not REPORT.is_file():
        raise FileNotFoundError(f"Generation report not found: {REPORT}")

    data = json.loads(REPORT.read_text())

    print("========== LIVE TARGET COVERAGE AUDIT ==========")
    print(f"Report: {REPORT}")
    print(f"Top-level keys: {list(data.keys())}")

    walk(data)

    print("\n========== END OF READ-ONLY AUDIT ==========")


if __name__ == "__main__":
    main()
