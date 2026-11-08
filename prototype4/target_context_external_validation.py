#!/usr/bin/env python3
"""Paired external test: sequence-only P3 versus target-context P4-TA.

Primary endpoint: leave-one-Binding-Group-out prediction of directly reported
PEPBI ITC Delta-G.  P3 uses the native sequence-only C score. P4-TA adds one
precomputed, complex-specific Rosetta Interface Analyzer feature (dG_separated).
All preprocessing and coefficient fitting occur inside each outer training fold.

This validates transfer of target-context structural evidence on PEPBI; it is
NOT a direct HERV-K Env binding experiment and does not establish 9SYA efficacy.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "prototype4/data/pepbi_target_context_external_validation.csv"
DEFAULT_OUT = ROOT / "prototype4/results/target_context_external_validation"
ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _native_scores(sequences: list[str]) -> dict[str, float]:
    # Import only after the CLI has validated input paths. The C API deliberately
    # evaluates target=NULL: this is the sequence-only P3 prior comparator.
    sys.path.insert(0, str(ROOT))
    import ctypes
    from c import bridge

    result: dict[str, float] = {}
    failures = []
    for sequence in sequences:
        metrics = bridge.NativeFitnessMetrics()
        status = bridge.lib.c_evaluate_sequence(
            sequence.encode("ascii"), ctypes.byref(metrics)
        )
        if status != 0:
            failures.append({"sequence": sequence, "status": int(status)})
            continue
        value = float(metrics.fitness_score)
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            failures.append({"sequence": sequence, "status": "invalid_score"})
            continue
        result[sequence] = value
    if failures:
        raise RuntimeError(
            "Native score failed for eligible benchmark sequences: "
            + json.dumps(failures[:10])
        )
    return result


def _read_rows(path: Path) -> list[dict[str, Any]]:
    required = [
        "Binding Group", "PEPBI Complex Name", "PDB ID", "Protein Sequence",
        "Peptide Sequence", "Peptide Length", "ΔG (kcal/mol)", "T (K)",
        "dG_separated", "dSASA_int", "nres_int", "hbonds_int", "packstat",
    ]
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        missing = [x for x in required if x not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Input is missing required columns: {missing}")
        raw = list(reader)

    rows = []
    for i, row in enumerate(raw, start=2):
        sequence = "".join((row.get("Peptide Sequence") or "").split()).upper()
        try:
            declared_length = int(float(row.get("Peptide Length", "")))
        except (TypeError, ValueError):
            continue
        if not (9 <= len(sequence) <= 21 and len(sequence) == declared_length):
            continue
        if not sequence or not set(sequence).issubset(ALPHABET):
            continue
        group = (row.get("Binding Group") or "").strip()
        if not group:
            continue
        experimental_dg = _finite(row.get("ΔG (kcal/mol)"))
        temperature = _finite(row.get("T (K)"))
        target_dg = _finite(row.get("dG_separated"))
        dsasa = _finite(row.get("dSASA_int"))
        nres = _finite(row.get("nres_int"))
        hbonds = _finite(row.get("hbonds_int"))
        packstat = _finite(row.get("packstat"))
        if None in (experimental_dg, temperature, target_dg):
            continue
        rows.append({
            "source_row": i,
            "binding_group": group,
            "complex_name": (row.get("PEPBI Complex Name") or "").strip(),
            "pdb_id": (row.get("PDB ID") or "").strip(),
            "protein_sequence": "".join((row.get("Protein Sequence") or "").split()).upper(),
            "peptide_sequence": sequence,
            "peptide_length": len(sequence),
            "experimental_dg_kcal_mol": experimental_dg,
            "temperature_K": temperature,
            "ria_dg_separated": target_dg,
            "ria_dSASA_int": dsasa,
            "ria_nres_int": nres,
            "ria_hbonds_int": hbonds,
            "ria_packstat": packstat,
        })
    if not rows:
        raise ValueError("No complete PEPBI rows remained after the predeclared 9–21 aa filter.")
    return rows


def _ridge_fit_predict(x_train: np.ndarray, y_train: np.ndarray,
                       x_test: np.ndarray, alpha: float) -> np.ndarray:
    """Ridge with training-fold-only standardization and unpenalized intercept."""
    x_mean = x_train.mean(axis=0)
    x_scale = x_train.std(axis=0, ddof=0)
    x_scale[x_scale < 1e-12] = 1.0
    z_train = (x_train - x_mean) / x_scale
    z_test = (x_test - x_mean) / x_scale
    y_mean = float(y_train.mean())
    yc = y_train - y_mean
    gram = z_train.T @ z_train
    coef = np.linalg.solve(gram + float(alpha) * np.eye(gram.shape[0]), z_train.T @ yc)
    return y_mean + z_test @ coef


def _mae(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean(np.abs(a - b)))


def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def _spearman(a: np.ndarray, b: np.ndarray) -> float | None:
    if len(a) < 3 or len(np.unique(a)) < 2 or len(np.unique(b)) < 2:
        return None
    value = float(spearmanr(a, b).statistic)
    return value if math.isfinite(value) else None


def _bootstrap_ci(values: np.ndarray, seed: int, reps: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    samples = np.empty(reps, dtype=float)
    for i in range(reps):
        samples[i] = rng.choice(values, size=n, replace=True).mean()
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def _exact_signflip_p(values: np.ndarray) -> float:
    """One-sided paired sign-flip p-value for mean benefit > 0."""
    n = len(values)
    observed = float(values.mean())
    if n > 20:
        raise ValueError("Exact sign-flip is limited to 20 independent groups.")
    count = total = 0
    for signs in itertools.product((-1.0, 1.0), repeat=n):
        statistic = float(np.mean(values * np.asarray(signs)))
        total += 1
        if statistic >= observed - 1e-15:
            count += 1
    return count / total


def run(input_path: Path, output_dir: Path, alpha: float = 1.0,
        bootstrap_reps: int = 20000, seed: int = 20261103) -> dict[str, Any]:
    rows = _read_rows(input_path)
    seqs = sorted({r["peptide_sequence"] for r in rows})
    scores = _native_scores(seqs)
    for row in rows:
        row["native_sequence_score"] = scores[row["peptide_sequence"]]

    groups = sorted({r["binding_group"] for r in rows})
    if len(groups) < 5:
        raise RuntimeError(
            f"Only {len(groups)} held-out binding groups; require >=5 for exploratory external validation."
        )

    # Keep an identical complete-case row cohort for P3 and every P4 variant.
    target_rows = np.asarray([r["experimental_dg_kcal_mol"] for r in rows], dtype=float)
    p3_x = np.asarray([[r["native_sequence_score"]] for r in rows], dtype=float)
    p4_x = np.asarray([
        [r["native_sequence_score"], r["ria_dg_separated"]]
        for r in rows
    ], dtype=float)
    p4_geometry_x = np.asarray([
        [r["native_sequence_score"], r["ria_dg_separated"],
         (r["ria_dSASA_int"] if r["ria_dSASA_int"] is not None else 0.0),
         (r["ria_nres_int"] if r["ria_nres_int"] is not None else 0.0)]
        for r in rows
    ], dtype=float)

    preds = {"P3_sequence_only": np.full(len(rows), np.nan),
             "P4_target_context_primary": np.full(len(rows), np.nan),
             "P4_target_context_geometry_sensitivity": np.full(len(rows), np.nan),
             "training_group_mean": np.full(len(rows), np.nan)}
    row_group = np.asarray([r["binding_group"] for r in rows], dtype=object)
    fold_records = []

    # Leave one full binding group out: same train/test rows for all models.
    for fold_no, heldout in enumerate(groups, start=1):
        test_idx = np.flatnonzero(row_group == heldout)
        train_idx = np.flatnonzero(row_group != heldout)
        if len(test_idx) == 0 or len(train_idx) < 5:
            continue
        ytr = target_rows[train_idx]
        preds["P3_sequence_only"][test_idx] = _ridge_fit_predict(
            p3_x[train_idx], ytr, p3_x[test_idx], alpha
        )
        preds["P4_target_context_primary"][test_idx] = _ridge_fit_predict(
            p4_x[train_idx], ytr, p4_x[test_idx], alpha
        )
        preds["P4_target_context_geometry_sensitivity"][test_idx] = _ridge_fit_predict(
            p4_geometry_x[train_idx], ytr, p4_geometry_x[test_idx], alpha
        )
        preds["training_group_mean"][test_idx] = float(ytr.mean())
        fold_records.append({
            "fold": fold_no,
            "heldout_binding_group": heldout,
            "n_train_rows": int(len(train_idx)),
            "n_test_rows": int(len(test_idx)),
        })

    if any(not np.isfinite(x).all() for x in preds.values()):
        raise RuntimeError("LOBO failed to produce predictions for every matched row.")

    output_dir.mkdir(parents=True, exist_ok=True)
    prediction_path = output_dir / "out_of_group_predictions.csv"
    with prediction_path.open("w", newline="", encoding="utf-8") as f:
        fields = list(rows[0].keys()) + list(preds.keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for i, row in enumerate(rows):
            w.writerow({**row, **{name: float(values[i]) for name, values in preds.items()}})

    # Independent inferential unit is Binding Group, not rows or pairwise combinations.
    group_rows = []
    for group in groups:
        idx = np.flatnonzero(row_group == group)
        record = {"binding_group": group, "n_rows": int(len(idx))}
        for name, yhat in preds.items():
            record[f"{name}_mae"] = _mae(target_rows[idx], yhat[idx])
            record[f"{name}_rmse"] = _rmse(target_rows[idx], yhat[idx])
            record[f"{name}_spearman"] = _spearman(target_rows[idx], yhat[idx])
        record["P3_minus_P4_primary_MAE"] = (
            record["P3_sequence_only_mae"]
            - record["P4_target_context_primary_mae"]
        )
        record["P3_minus_P4_geometry_MAE"] = (
            record["P3_sequence_only_mae"]
            - record["P4_target_context_geometry_sensitivity_mae"]
        )
        group_rows.append(record)

    group_path = output_dir / "binding_group_metrics.csv"
    with group_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(group_rows[0].keys()))
        w.writeheader(); w.writerows(group_rows)

    paired = np.asarray([x["P3_minus_P4_primary_MAE"] for x in group_rows], dtype=float)
    geom_paired = np.asarray([x["P3_minus_P4_geometry_MAE"] for x in group_rows], dtype=float)
    ci = _bootstrap_ci(paired, seed, bootstrap_reps)
    p = _exact_signflip_p(paired)
    # Directional and two-sided exact tests use binding group as the unit.
    p_worse = _exact_signflip_p(-paired)
    p_two_sided = min(1.0, 2.0 * min(p, p_worse))

    model_summary = []
    for name, yhat in preds.items():
        macro_mae = float(np.mean([x[f"{name}_mae"] for x in group_rows]))
        macro_rmse = float(np.mean([x[f"{name}_rmse"] for x in group_rows]))
        model_summary.append({
            "model": name,
            "n_rows": int(len(rows)),
            "n_binding_groups": len(groups),
            "pooled_mae_kcal_mol": _mae(target_rows, yhat),
            "pooled_rmse_kcal_mol": _rmse(target_rows, yhat),
            "pooled_spearman_rho": _spearman(target_rows, yhat),
            "macro_binding_group_mae_kcal_mol": macro_mae,
            "macro_binding_group_rmse_kcal_mol": macro_rmse,
        })
    summary_path = output_dir / "model_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(model_summary[0].keys()))
        w.writeheader(); w.writerows(model_summary)

    paired_path = output_dir / "paired_group_effects.csv"
    with paired_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["binding_group", "n_rows", "p3_mae", "p4_target_context_mae", "p3_minus_p4_mae_benefit"])
        w.writeheader()
        for x in group_rows:
            w.writerow({
                "binding_group": x["binding_group"], "n_rows": x["n_rows"],
                "p3_mae": x["P3_sequence_only_mae"],
                "p4_target_context_mae": x["P4_target_context_primary_mae"],
                "p3_minus_p4_mae_benefit": x["P3_minus_P4_primary_MAE"],
            })

    p3s = preds["P3_sequence_only"]
    p4s = preds["P4_target_context_primary"]
    direct_feature_correlations = {
        "native_sequence_score_vs_experimental_dG": _spearman(
            np.asarray([r["native_sequence_score"] for r in rows], dtype=float),
            target_rows,
        ),
        "ria_dG_separated_vs_experimental_dG": _spearman(
            np.asarray([r["ria_dg_separated"] for r in rows], dtype=float),
            target_rows,
        ),
        "note": (
            "Descriptive pooled correlations only; inferential comparison uses the "
            "leave-one-Binding-Group-out predictions and binding-group-level paired errors."
        ),
    }
    result = {
        "experiment": "P3 sequence-only versus P4 target-context structural augmentation",
        "source": str(input_path),
        "source_sha256": _sha256(input_path),
        "outcome": "directly reported PEPBI ITC ΔG (kcal/mol); lower is stronger binding",
        "target_context_feature": "Rosetta Interface Analyzer dG_separated from the corresponding peptide-protein complex; lower is more favorable",
        "P3_features": ["native_sequence_score"],
        "P4_primary_features": ["native_sequence_score", "ria_dg_separated"],
        "P4_geometry_sensitivity_features": ["native_sequence_score", "ria_dg_separated", "ria_dSASA_int", "ria_nres_int"],
        "validation": "leave-one-Binding-Group-out; no rows from a held-out binding group enter training",
        "preprocessing": "feature means/scales and Ridge coefficients are fitted on outer training fold only",
        "ridge_alpha_predeclared": alpha,
        "rows": len(rows),
        "unique_peptide_sequences": len(seqs),
        "independent_binding_groups": len(groups),
        "binding_groups": groups,
        "paired_primary_effect": {
            "definition": "P3 group MAE minus P4-TA group MAE; positive favors adding target-context structural evidence",
            "macro_mae_benefit_kcal_mol": float(paired.mean()),
            "group_bootstrap_95ci": [ci[0], ci[1]],
            "exact_one_sided_signflip_p": p,
            "exact_one_sided_signflip_p_for_p4_worsening": p_worse,
            "exact_two_sided_signflip_p": p_two_sided,
            "groups_improved": int((paired > 0).sum()),
            "groups_worsened": int((paired < 0).sum()),
            "groups_tied": int((np.abs(paired) <= 1e-12).sum()),
            "n_independent_groups": int(len(paired)),
        },
        "geometry_sensitivity_macro_mae_benefit_kcal_mol": float(geom_paired.mean()),
        "folds": fold_records,
        "models": model_summary,
        "direct_feature_correlations": direct_feature_correlations,
        "interpretation": (
            "This is an external transfer test of adding receptor-context RIA evidence to a sequence-only prior. "
            "It is not a direct 9SYA/HERV-K Env interaction measurement. P4 benefit requires lower held-out "
            "group-macro MAE with uncertainty supporting a positive paired difference; a positive point estimate "
            "alone is not proof. The PEPBI benchmark contains only the reported independent binding groups listed here."
        ),
        "artifacts": {
            "predictions": prediction_path.name,
            "group_metrics": group_path.name,
            "model_summary": summary_path.name,
            "paired_group_effects": paired_path.name,
        },
    }
    from src.p4_external_validation_gate import evaluate_external_validation_summary
    result["target_context_transfer_decision"] = evaluate_external_validation_summary(result)
    (output_dir / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (output_dir / "external_validation_gate.json").write_text(
        json.dumps(result["target_context_transfer_decision"], indent=2) + "\n",
        encoding="utf-8",
    )

    # Poster-ready primary plot: one point per independent binding group.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.2, 5.2))
        names = [x["binding_group"] for x in group_rows]
        x = np.arange(len(names))
        width = 0.38
        ax.bar(x - width / 2, [x["P3_sequence_only_mae"] for x in group_rows], width, label="P3: sequence-only")
        ax.bar(x + width / 2, [x["P4_target_context_primary_mae"] for x in group_rows], width, label="P4-TA: sequence + target-context RIA")
        ax.set_ylabel("Held-out experimental ΔG prediction MAE (kcal/mol)")
        ax.set_title("Target-context structural evidence: paired leave-one-group-out test")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "target_context_p3_vs_p4_group_mae.png", dpi=220)
        plt.close(fig)
    except Exception as exc:
        result["plot_warning"] = f"Could not render plot: {type(exc).__name__}: {exc}"
        (output_dir / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--bootstrap-reps", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20261103)
    args = parser.parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input file not found: {args.input}")
    if args.alpha <= 0 or args.bootstrap_reps < 1000:
        raise SystemExit("Require alpha > 0 and bootstrap-reps >= 1000.")
    result = run(args.input, args.out, args.alpha, args.bootstrap_reps, args.seed)
    print(json.dumps(result, indent=2))
    print(f"\nArtifacts: {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
