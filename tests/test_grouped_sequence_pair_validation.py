import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "prototype3" / "scripts" / "run_grouped_sequence_pair_validation.py"
SPEC = importlib.util.spec_from_file_location("grouped_validation", SCRIPT)
mod = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(mod)


def test_pair_encoder_requires_exactly_one_canonical_substitution():
    a = "ACDEFGHIKL"
    b = "VCDEFGHIKL"
    x = mod.encode_pair(a, b)
    assert x.shape == (20 * 20 + 12 + 1 + 40,)
    assert x.sum() == 3  # substitution + length + position + available right flank
    try:
        mod.encode_pair(a, "VCDEYGHIKL")
    except ValueError as exc:
        assert "exactly one" in str(exc)
    else:
        raise AssertionError("Two substitutions should be rejected")


def test_leakage_components_join_shared_target_and_repeated_pair():
    frame = pd.DataFrame({
        "target_sequence": ["TARGET_A", "TARGET_A", "TARGET_B"],
        "pair_signature": ["PAIR_1", "PAIR_2", "PAIR_1"],
        "peptide_a": ["PEP_A", "PEP_B", "PEP_A"],
        "peptide_b": ["PEP_C", "PEP_D", "PEP_E"],
    })
    groups = mod.make_leakage_components(frame)
    assert groups[0] == groups[1] == groups[2]


def _synthetic_csv(path: Path) -> None:
    aa = mod.AA
    rows = []
    for i in range(6):
        target = "".join(aa[(i + j) % len(aa)] for j in range(25))
        base = "".join(aa[(i * 3 + j) % len(aa)] for j in range(10))
        # Keep each target's peptide-pair signatures distinct while creating
        # two exact single-substitution variants of the same base peptide.
        repl0 = aa[(mod.AA_INDEX[base[0]] + 7 + i) % len(aa)]
        repl5 = aa[(mod.AA_INDEX[base[5]] + 9 + i) % len(aa)]
        v0 = repl0 + base[1:]
        v5 = base[:5] + repl5 + base[6:]
        for sequence, pkd in [(base, 5.0 + 0.1 * i), (v0, 5.4 + 0.1 * i), (v5, 4.7 + 0.1 * i)]:
            rows.append({"_peptide": sequence, "_protein": target, "_assay": "assay_X", "_p_affinity": pkd})
    pd.DataFrame(rows).to_csv(path, index=False)


def test_end_to_end_grouped_validation_has_no_target_or_pair_leakage(tmp_path):
    inp = tmp_path / "synthetic.csv"
    _synthetic_csv(inp)
    out = mod.run_validation(inp, tmp_path / "results", n_splits=3, bootstrap_replicates=50, seed=42)
    metrics = __import__("json").loads((out / "metrics.json").read_text())
    assert metrics["single_substitution_pairs"] == 12
    assert metrics["unique_target_sequences"] == 6
    assert metrics["split_leakage_checks"]["all_folds_have_zero_target_overlap"] is True
    assert metrics["split_leakage_checks"]["all_folds_have_zero_exact_pair_overlap"] is True
    assert metrics["split_leakage_checks"]["all_folds_have_zero_peptide_sequence_overlap"] is True
    preds = pd.read_csv(out / "out_of_fold_predictions.csv")
    assert len(preds) == 12
    assert preds["prediction_sequence_only_ridge"].notna().all()
    audit = pd.read_json(out / "fold_split_audit.json")
    assert (audit["target_sequence_overlap"] == 0).all()
    assert (audit["exact_pair_signature_overlap"] == 0).all()
    assert (audit["peptide_sequence_overlap"] == 0).all()
