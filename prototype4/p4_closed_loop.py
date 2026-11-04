from __future__ import annotations

import argparse
import hashlib
import hmac
import json
import math
import secrets
import sqlite3
from pathlib import Path
from typing import Any

import numpy as np

VERSION = "p4_closed_loop_v1.0"
JACCARD_CUT_HEIGHT = 0.15
EXPLORATION_FRACTION = 0.15
EVSI_FUTURE_REPLICATES = 10

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS skempi_raw (
    raw_id INTEGER PRIMARY KEY,
    source_version TEXT NOT NULL,
    source_record_key TEXT,
    original_line TEXT NOT NULL,
    chain_partner_1 TEXT,
    chain_partner_2 TEXT,
    experimental_method TEXT,
    raw_sha256 TEXT NOT NULL UNIQUE,
    ingested_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TRIGGER IF NOT EXISTS skempi_raw_no_update
BEFORE UPDATE ON skempi_raw
BEGIN
    SELECT RAISE(ABORT, 'skempi_raw is immutable');
END;

CREATE TRIGGER IF NOT EXISTS skempi_raw_no_delete
BEFORE DELETE ON skempi_raw
BEGIN
    SELECT RAISE(ABORT, 'skempi_raw is immutable');
END;

CREATE TABLE IF NOT EXISTS cycle_nodes (
    node_id TEXT PRIMARY KEY,
    cycle_uid TEXT NOT NULL,
    state TEXT NOT NULL CHECK(state IN ('WT','M1','M2','M1M2')),
    raw_id INTEGER,
    mutation_partner_1 TEXT,
    mutation_partner_2 TEXT,
    mutation_1 TEXT,
    mutation_2 TEXT,
    topology TEXT CHECK(
        topology IN (
            'same_partner_1',
            'same_partner_2',
            'cross_partner_1x2',
            'unresolved'
        )
    ),
    interface_proximity TEXT CHECK(
        interface_proximity IN ('interface','distal','unresolved')
    ),
    FOREIGN KEY(raw_id) REFERENCES skempi_raw(raw_id),
    UNIQUE(cycle_uid, state)
);

CREATE INDEX IF NOT EXISTS idx_cycle_nodes_cycle
ON cycle_nodes(cycle_uid);

CREATE TABLE IF NOT EXISTS structure_runs (
    run_id TEXT PRIMARY KEY,
    cycle_uid TEXT,
    target_key TEXT,
    run_role TEXT,
    software TEXT NOT NULL,
    software_version TEXT NOT NULL,
    command_string TEXT NOT NULL,
    config_sha256 TEXT NOT NULL,
    random_seed INTEGER,
    input_sha256 TEXT NOT NULL,
    output_sha256 TEXT,
    trajectory_path TEXT,
    l_rmsd REAL,
    status TEXT NOT NULL,
    started_at TEXT,
    completed_at TEXT,
    metadata_json TEXT
);

CREATE INDEX IF NOT EXISTS idx_structure_runs_cycle
ON structure_runs(cycle_uid);

CREATE TABLE IF NOT EXISTS cycle_features (
    feature_id INTEGER PRIMARY KEY,
    cycle_uid TEXT NOT NULL,
    feature_name TEXT NOT NULL,
    feature_value REAL,
    feature_json TEXT,
    run_wt_id TEXT,
    run_m1_id TEXT,
    run_m2_id TEXT,
    run_m1m2_id TEXT,
    source_sha256 TEXT,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY(run_wt_id) REFERENCES structure_runs(run_id),
    FOREIGN KEY(run_m1_id) REFERENCES structure_runs(run_id),
    FOREIGN KEY(run_m2_id) REFERENCES structure_runs(run_id),
    FOREIGN KEY(run_m1m2_id) REFERENCES structure_runs(run_id),
    UNIQUE(cycle_uid, feature_name)
);

CREATE TABLE IF NOT EXISTS loco_predictions (
    prediction_id INTEGER PRIMARY KEY,
    model_version TEXT NOT NULL,
    cycle_uid TEXT NOT NULL,
    outer_fold INTEGER NOT NULL,
    heldout_group TEXT NOT NULL,
    y_true REAL NOT NULL,
    y_pred REAL NOT NULL,
    train_groups_sha256 TEXT NOT NULL,
    feature_set TEXT NOT NULL,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(model_version, cycle_uid, outer_fold)
);

CREATE TABLE IF NOT EXISTS validation_manifest (
    manifest_id TEXT PRIMARY KEY,
    manifest_json TEXT NOT NULL,
    content_sha256 TEXT NOT NULL,
    signature_hmac_sha256 TEXT,
    signature_key_id TEXT,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TRIGGER IF NOT EXISTS validation_manifest_no_update
BEFORE UPDATE ON validation_manifest
BEGIN
    SELECT RAISE(ABORT, 'validation_manifest is immutable');
END;

CREATE TRIGGER IF NOT EXISTS validation_manifest_no_delete
BEFORE DELETE ON validation_manifest
BEGIN
    SELECT RAISE(ABORT, 'validation_manifest is immutable');
END;
"""


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def open_registry(path: str | Path) -> sqlite3.Connection:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("PRAGMA foreign_keys = ON")
    con.executescript(SCHEMA)
    con.commit()
    return con


def add_raw_record(
    con: sqlite3.Connection,
    original_line: str,
    source_version: str,
    source_record_key: str | None = None,
    chain_partner_1: str | None = None,
    chain_partner_2: str | None = None,
    experimental_method: str | None = None,
) -> int:
    digest = sha256_text(original_line)
    con.execute(
        """
        INSERT OR IGNORE INTO skempi_raw
        (source_version, source_record_key, original_line,
         chain_partner_1, chain_partner_2, experimental_method, raw_sha256)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        (
            source_version, source_record_key, original_line,
            chain_partner_1, chain_partner_2, experimental_method, digest
        ),
    )
    con.commit()
    row = con.execute(
        "SELECT raw_id FROM skempi_raw WHERE raw_sha256 = ?", (digest,)
    ).fetchone()
    return int(row[0])


def add_cycle_node(
    con: sqlite3.Connection,
    node_id: str,
    cycle_uid: str,
    state: str,
    **metadata: Any,
) -> None:
    allowed = {
        "raw_id", "mutation_partner_1", "mutation_partner_2",
        "mutation_1", "mutation_2", "topology", "interface_proximity"
    }
    unknown = set(metadata) - allowed
    if unknown:
        raise ValueError(f"Unknown cycle-node fields: {sorted(unknown)}")

    fields = ["node_id", "cycle_uid", "state"] + list(metadata)
    values = [node_id, cycle_uid, state] + list(metadata.values())
    sql = (
        f"INSERT INTO cycle_nodes ({','.join(fields)}) "
        f"VALUES ({','.join('?' for _ in fields)})"
    )
    con.execute(sql, values)
    con.commit()


def add_structure_run(con: sqlite3.Connection, record: dict[str, Any]) -> None:
    required = {
        "run_id", "software", "software_version", "command_string",
        "config_sha256", "input_sha256", "status"
    }
    missing = required - set(record)
    if missing:
        raise ValueError(f"Missing structure-run fields: {sorted(missing)}")

    columns = [
        "run_id", "cycle_uid", "target_key", "run_role", "software",
        "software_version", "command_string", "config_sha256", "random_seed",
        "input_sha256", "output_sha256", "trajectory_path", "l_rmsd",
        "status", "started_at", "completed_at", "metadata_json"
    ]
    values = []
    for column in columns:
        value = record.get(column)
        if column == "metadata_json" and isinstance(value, dict):
            value = canonical_json(value)
        values.append(value)

    con.execute(
        f"INSERT INTO structure_runs ({','.join(columns)}) "
        f"VALUES ({','.join('?' for _ in columns)})",
        values,
    )
    con.commit()


def seal_manifest(
    con: sqlite3.Connection,
    manifest: dict[str, Any],
    key_path: str | Path,
) -> str:
    """
    HMAC-SHA256 authenticated integrity seal.
    This is a local symmetric signature, not a public-key signature.
    Keep the key private and back it up separately.
    """
    key_path = Path(key_path)
    key_path.parent.mkdir(parents=True, exist_ok=True)

    if not key_path.exists():
        key_path.write_bytes(secrets.token_bytes(32))
        key_path.chmod(0o600)

    key = key_path.read_bytes()
    body = canonical_json(manifest)
    digest = sha256_text(body)
    signature = hmac.new(key, body.encode(), hashlib.sha256).hexdigest()
    manifest_id = sha256_text(VERSION + body)[:24]
    key_id = hashlib.sha256(key).hexdigest()[:16]

    con.execute(
        """
        INSERT INTO validation_manifest
        (manifest_id, manifest_json, content_sha256,
         signature_hmac_sha256, signature_key_id)
        VALUES (?, ?, ?, ?, ?)
        """,
        (manifest_id, body, digest, signature, key_id),
    )
    con.commit()
    return manifest_id


def verify_manifest(
    con: sqlite3.Connection,
    manifest_id: str,
    key_path: str | Path,
) -> bool:
    row = con.execute(
        """
        SELECT manifest_json, content_sha256,
               signature_hmac_sha256, signature_key_id
        FROM validation_manifest WHERE manifest_id = ?
        """,
        (manifest_id,),
    ).fetchone()
    if row is None:
        return False

    body, digest, signature, key_id = row
    key = Path(key_path).read_bytes()
    expected = hmac.new(key, body.encode(), hashlib.sha256).hexdigest()

    return (
        sha256_text(body) == digest
        and hmac.compare_digest(expected, signature or "")
        and hashlib.sha256(key).hexdigest()[:16] == key_id
    )


# ---------------------------------------------------------------------
# MODEL C: CONTACT-FINGERPRINT ENTROPY AND CROSS-RUN REPRODUCIBILITY
# ---------------------------------------------------------------------

def jaccard_distance(a: set[str], b: set[str]) -> float:
    union = a | b
    if not union:
        return 0.0
    return 1.0 - len(a & b) / len(union)


def complete_linkage_jaccard(
    fingerprints: list[set[str]],
    cut_height: float = JACCARD_CUT_HEIGHT,
) -> list[list[int]]:
    """
    Agglomerative complete-linkage clustering using Jaccard distance.
    Clusters merge only while their maximum cross-cluster distance
    is <= cut_height.
    """
    n = len(fingerprints)
    if n == 0:
        raise ValueError("At least one contact fingerprint is required")

    distances = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            distances[i, j] = distances[j, i] = jaccard_distance(
                fingerprints[i], fingerprints[j]
            )

    clusters = [[i] for i in range(n)]

    while True:
        best = None
        best_distance = float("inf")

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = max(
                    distances[a, b]
                    for a in clusters[i]
                    for b in clusters[j]
                )
                if d < best_distance:
                    best_distance = d
                    best = (i, j)

        if best is None or best_distance > cut_height:
            break

        i, j = best
        merged = sorted(clusters[i] + clusters[j])
        clusters = [
            c for k, c in enumerate(clusters) if k not in (i, j)
        ]
        clusters.append(merged)
        clusters.sort(key=lambda c: c[0])

    return clusters


def contact_entropy(
    fingerprints: list[set[str]],
    cut_height: float = JACCARD_CUT_HEIGHT,
) -> dict[str, Any]:
    clusters = complete_linkage_jaccard(fingerprints, cut_height)
    counts = np.array([len(c) for c in clusters], dtype=float)
    probabilities = counts / counts.sum()

    entropy = float(-np.sum(probabilities * np.log(probabilities)))
    max_entropy = math.log(len(fingerprints)) if len(fingerprints) > 1 else 0.0
    normalized = entropy / max_entropy if max_entropy > 0 else 0.0

    labels = [None] * len(fingerprints)
    for cluster_id, cluster in enumerate(clusters):
        for frame_index in cluster:
            labels[frame_index] = cluster_id

    return {
        "n_frames": len(fingerprints),
        "n_clusters": len(clusters),
        "cluster_sizes": counts.astype(int).tolist(),
        "cluster_labels": labels,
        "entropy_nats": entropy,
        "normalized_entropy": float(normalized),
        "effective_cluster_count": float(math.exp(entropy)),
        "cut_height": float(cut_height),
        "distance": "Jaccard",
        "linkage": "complete",
    }


def run_consensus(
    fingerprints: list[set[str]],
    frequency_cutoff: float = 0.5,
) -> set[str]:
    if not fingerprints:
        raise ValueError("A run must contain at least one frame")
    counts: dict[str, int] = {}
    for frame in fingerprints:
        for contact in frame:
            counts[contact] = counts.get(contact, 0) + 1
    n = len(fingerprints)
    return {
        contact for contact, count in counts.items()
        if count / n >= frequency_cutoff
    }


def cross_run_reliability(
    runs: list[list[set[str]]],
) -> dict[str, Any]:
    if len(runs) < 2:
        raise ValueError("At least two independent runs are required")

    run_metrics = [contact_entropy(run) for run in runs]
    consensuses = [run_consensus(run) for run in runs]

    pairwise = []
    for i in range(len(consensuses)):
        for j in range(i + 1, len(consensuses)):
            pairwise.append(
                1.0 - jaccard_distance(consensuses[i], consensuses[j])
            )

    mean_entropy = float(np.mean([
        x["normalized_entropy"] for x in run_metrics
    ]))
    stability = float(np.mean(pairwise))

    # Uncalibrated diagnostic proxy only. It is not pose accuracy.
    proxy = 0.5 * (1.0 - mean_entropy) + 0.5 * stability

    return {
        "n_independent_runs": len(runs),
        "run_entropy": run_metrics,
        "mean_normalized_entropy": mean_entropy,
        "mean_cross_run_jaccard_similarity": stability,
        "reliability_proxy": proxy,
        "calibrated": False,
        "interpretation": (
            "Reproducibility diagnostic only; not experimental binding "
            "probability or structural accuracy."
        ),
    }


# ---------------------------------------------------------------------
# EXTERNAL CALIBRATION: STRICT GROUPED OUT-OF-FOLD RIDGE REGRESSION
# ---------------------------------------------------------------------

def _ridge_predict(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    alpha: float,
) -> np.ndarray:
    mean = x_train.mean(axis=0)
    scale = x_train.std(axis=0)
    scale[scale < 1e-12] = 1.0

    xt = (x_train - mean) / scale
    xv = (x_test - mean) / scale

    xt = np.column_stack([np.ones(len(xt)), xt])
    xv = np.column_stack([np.ones(len(xv)), xv])

    penalty = np.eye(xt.shape[1]) * alpha
    penalty[0, 0] = 0.0

    beta = np.linalg.solve(xt.T @ xt + penalty, xt.T @ y_train)
    return xv @ beta


def grouped_oof_ridge(
    features: np.ndarray,
    labels: np.ndarray,
    groups: np.ndarray,
    alpha: float = 1.0,
) -> np.ndarray:
    features = np.asarray(features, dtype=float)
    labels = np.asarray(labels, dtype=float)
    groups = np.asarray(groups)

    if features.ndim != 2 or len(features) != len(labels):
        raise ValueError("Feature matrix and labels have incompatible shapes")
    if len(groups) != len(labels):
        raise ValueError("Group vector length does not match labels")
    unique_groups = np.unique(groups)
    if len(unique_groups) < 3:
        raise ValueError("At least three independent families are required")

    predictions = np.full(len(labels), np.nan)

    for group in unique_groups:
        test = groups == group
        train = ~test
        if train.sum() <= features.shape[1]:
            raise ValueError(
                "Too few training rows for this grouped outer fold; "
                "reduce features or add independent families."
            )
        predictions[test] = _ridge_predict(
            features[train], labels[train], features[test], alpha
        )

    if np.isnan(predictions).any():
        raise RuntimeError("Some rows did not receive an outer-fold prediction")
    return predictions


def compare_entropy_incremental(
    baseline_features: np.ndarray,
    entropy_feature: np.ndarray,
    labels: np.ndarray,
    groups: np.ndarray,
    alpha: float = 1.0,
    bootstrap_replicates: int = 10000,
    seed: int = 20261031,
) -> dict[str, Any]:
    """
    Positive delta_rmse means entropy improves out-of-family prediction.
    CI is bootstrapped by family, not by individual frames/rows.
    """
    x0 = np.asarray(baseline_features, dtype=float)
    h = np.asarray(entropy_feature, dtype=float).reshape(-1, 1)
    y = np.asarray(labels, dtype=float)
    g = np.asarray(groups)

    p0 = grouped_oof_ridge(x0, y, g, alpha)
    p1 = grouped_oof_ridge(np.column_stack([x0, h]), y, g, alpha)

    rmse0 = float(np.sqrt(np.mean((y - p0) ** 2)))
    rmse1 = float(np.sqrt(np.mean((y - p1) ** 2)))
    delta = rmse0 - rmse1

    group_deltas = []
    unique_groups = np.unique(g)
    for group in unique_groups:
        idx = g == group
        r0 = np.sqrt(np.mean((y[idx] - p0[idx]) ** 2))
        r1 = np.sqrt(np.mean((y[idx] - p1[idx]) ** 2))
        group_deltas.append(float(r0 - r1))

    rng = np.random.default_rng(seed)
    group_deltas = np.asarray(group_deltas)
    boot = np.empty(bootstrap_replicates)

    for b in range(bootstrap_replicates):
        sample = rng.choice(
            group_deltas, size=len(group_deltas), replace=True
        )
        boot[b] = sample.mean()

    lo, hi = np.quantile(boot, [0.025, 0.975])

    return {
        "model_version": VERSION,
        "cv": "leave-one-family-out",
        "n_rows": len(y),
        "n_families": len(unique_groups),
        "baseline_oof_rmse": rmse0,
        "entropy_augmented_oof_rmse": rmse1,
        "delta_rmse_baseline_minus_augmented": float(delta),
        "family_bootstrap_ci95": [float(lo), float(hi)],
        "bootstrap_replicates": bootstrap_replicates,
        "entropy_incremental_validated": bool(lo > 0),
        "ci_excludes_zero": bool(lo > 0 or hi < 0),
        "entropy_improvement_supported": bool(lo > 0),
        "oof_predictions_baseline": p0.tolist(),
        "oof_predictions_entropy": p1.tolist(),
    }


# ---------------------------------------------------------------------
# MODEL B: VERSIONED TARGET-AWARE MULTI-OBJECTIVE RANKING ADAPTER
# ---------------------------------------------------------------------

def hotspot_proximity_score(d203: float, d205: float) -> float:
    """
    Geometric proximity diagnostic, not a binding energy.
    Distances are in Angstroms.
    """
    values = []
    for distance in (d203, d205):
        if not np.isfinite(distance) or distance < 0:
            raise ValueError("Invalid hotspot distance")
        values.append(math.exp(-max(distance - 4.0, 0.0) / 2.0))
    return float(np.mean(values))


def _pareto_fronts(
    objective_rows: list[dict[str, float]],
) -> list[int]:
    n = len(objective_rows)
    remaining = set(range(n))
    ranks = [-1] * n
    rank = 0

    while remaining:
        front = []
        for i in sorted(remaining):
            dominated = False
            for j in remaining:
                if i == j:
                    continue
                a = objective_rows[j]
                b = objective_rows[i]
                if all(a[k] >= b[k] for k in b) and any(
                    a[k] > b[k] for k in b
                ):
                    dominated = True
                    break
            if not dominated:
                front.append(i)

        if not front:
            raise RuntimeError("Pareto sorting failed")
        for i in front:
            ranks[i] = rank
            remaining.remove(i)
        rank += 1

    return ranks


def rank_target_candidates(
    candidates: list[dict[str, Any]],
    use_calibrated_reliability: bool = False,
) -> list[dict[str, Any]]:
    """
    Pareto ranking; does not alter c_score or historical fitness weights.
    Required keys: sequence_score, target_geometry_score,
                   arg203_min_A, arg205_min_A.
    """
    output = []
    objectives = []

    for candidate in candidates:
        row = dict(candidate)

        d203 = float(row["arg203_min_A"])
        d205 = float(row["arg205_min_A"])
        interface = float(row["target_geometry_score"])
        sequence = float(row["sequence_score"])

        row["target_objective_version"] = "HERVK_9SYA_target_pareto_v2"
        row["hotspot_proximity_score"] = hotspot_proximity_score(d203, d205)
        row["sulfate_displacement_flag"] = (
            row.get("sulfate_min_A") is not None
            and float(row["sulfate_min_A"]) < 4.0
        )
        row["target_eligible"] = not bool(row.get("severe_clash", False))

        objective = {
            "sequence_score": sequence,
            "target_geometry_score": interface,
            "hotspot_proximity_score": row["hotspot_proximity_score"],
        }

        if use_calibrated_reliability:
            if not row.get("reliability_calibrated", False):
                raise ValueError(
                    "Cannot use reliability in ranking until externally calibrated"
                )
            objective["reliability_score"] = float(
                row["reliability_score"]
            )

        # A severe-clash candidate is strictly dominated by a virtual
        # infeasible penalty in the final ranking.
        if not row["target_eligible"]:
            objective = {k: -1e12 for k in objective}

        objectives.append(objective)
        output.append(row)

    ranks = _pareto_fronts(objectives)
    for row, rank in zip(output, ranks):
        row["pareto_front"] = rank

    return sorted(
        output,
        key=lambda r: (not r["target_eligible"], r["pareto_front"])
    )


def attach_target_evidence(
    candidates: list[dict[str, Any]],
    evidence_by_sequence: dict[str, dict[str, Any]],
    use_calibrated_reliability: bool = False,
) -> list[dict[str, Any]]:
    enriched = []
    for candidate in candidates:
        row = dict(candidate)
        sequence = row.get("sequence")
        if not sequence or sequence not in evidence_by_sequence:
            raise KeyError(
                f"No target evidence for candidate sequence {sequence!r}"
            )
        row.update(evidence_by_sequence[sequence])
        enriched.append(row)
    return rank_target_candidates(enriched, use_calibrated_reliability)


# ---------------------------------------------------------------------
# CLOSED-LOOP FEEDBACK
# ---------------------------------------------------------------------

def parent_selection_weights(
    candidates: list[dict[str, Any]],
    exploration_fraction: float = EXPLORATION_FRACTION,
    use_calibrated_reliability: bool = False,
    temperature: float = 0.25,
) -> list[float]:
    if not candidates:
        raise ValueError("No candidates provided")
    if not 0 <= exploration_fraction < 1:
        raise ValueError("Exploration fraction must be in [0, 1)")

    values = np.asarray(
        [float(x["predicted_value"]) for x in candidates], dtype=float
    )
    vmin, vmax = float(values.min()), float(values.max())
    utility = (
        np.zeros_like(values)
        if vmax - vmin < 1e-12
        else (values - vmin) / (vmax - vmin)
    )

    if use_calibrated_reliability:
        for row in candidates:
            if not row.get("reliability_calibrated", False):
                raise ValueError(
                    "Reliability must be externally calibrated before "
                    "it can change parent-selection weights"
                )
        reliability = np.asarray([
            float(x["reliability_score"]) for x in candidates
        ])
        utility = utility * np.clip(reliability, 0.0, 1.0)

    logits = (utility - utility.max()) / max(temperature, 1e-6)
    exploit = np.exp(logits)
    exploit /= exploit.sum()

    n = len(candidates)
    weights = (
        (1.0 - exploration_fraction) * exploit
        + exploration_fraction / n
    )
    return weights.tolist()


def operator_gate(
    out_of_domain: bool,
    normalized_entropy: float,
    cross_run_similarity: float,
    reliability_calibrated: bool,
    entropy_limit: float = 0.50,
    stability_minimum: float = 0.70,
) -> dict[str, Any]:
    if not reliability_calibrated:
        return {
            "aggressive_operators_enabled": True,
            "local_operators_enabled": True,
            "gating_active": False,
            "reason": "Reliability model not externally calibrated",
        }

    reasons = []
    if out_of_domain:
        reasons.append("out_of_domain")
    if normalized_entropy > entropy_limit:
        reasons.append("high_contact_entropy")
    if cross_run_similarity < stability_minimum:
        reasons.append("low_cross_run_stability")

    gate = bool(reasons)
    return {
        "aggressive_operators_enabled": not gate,
        "local_operators_enabled": True,
        "gating_active": True,
        "reason": reasons or ["within_preregistered_limits"],
    }


def backbone_rmsd_radius_A(trajectory_backbones: Any) -> dict[str, Any]:
    """
    Calculate the mean Kabsch-aligned backbone RMSD to the medoid.

    Input shape:
        (n_independent_trajectories, n_consistently_mapped_backbone_atoms, 3)

    Supply one post-equilibration representative structure per independent
    trajectory. All structures must use the same atom selection and ordering.
    Docking poses from one HADDOCK execution are not independent trajectories.
    """
    xyz = np.asarray(trajectory_backbones, dtype=float)

    if xyz.ndim != 3 or xyz.shape[2] != 3:
        raise ValueError(
            "trajectory_backbones must have shape (n_trajectories, n_atoms, 3)"
        )

    n_trajectories, n_atoms, _ = xyz.shape
    if n_trajectories < 10:
        raise ValueError(
            f"Need at least 10 independent trajectories; received {n_trajectories}"
        )
    if n_atoms < 3:
        raise ValueError("At least three consistently mapped backbone atoms are required")
    if not np.isfinite(xyz).all():
        raise ValueError("Backbone coordinates contain non-finite values")

    def aligned_rmsd(reference, mobile):
        ref = reference - reference.mean(axis=0)
        mob = mobile - mobile.mean(axis=0)

        u, _, vt = np.linalg.svd(mob.T @ ref)
        correction = np.eye(3)
        correction[2, 2] = np.linalg.det(u @ vt)
        rotation = u @ correction @ vt

        fitted = mob @ rotation
        return float(np.sqrt(np.mean(np.sum((fitted - ref) ** 2, axis=1))))

    distances = np.zeros((n_trajectories, n_trajectories), dtype=float)
    for i in range(n_trajectories):
        for j in range(i + 1, n_trajectories):
            value = aligned_rmsd(xyz[i], xyz[j])
            distances[i, j] = value
            distances[j, i] = value

    mean_distance_by_structure = distances.mean(axis=1)
    medoid_index = int(np.argmin(mean_distance_by_structure))
    radius = float(mean_distance_by_structure[medoid_index])

    return {
        "n_independent_trajectories": int(n_trajectories),
        "n_backbone_atoms": int(n_atoms),
        "medoid_index": medoid_index,
        "mean_backbone_rmsd_radius_A": radius,
        "median_backbone_rmsd_to_medoid_A": float(
            np.median(distances[medoid_index])
        ),
        "rmsd_matrix_A": distances.tolist(),
        "alignment": "pairwise Kabsch least-squares superposition",
    }

def estimate_evsi(
    replicate_scores: dict[str, list[float]],
    future_replicates: int = EVSI_FUTURE_REPLICATES,
    simulations: int = 5000,
    seed: int = 20261031,
) -> dict[str, Any]:
    """
    Monte-Carlo normal-normal EVSI approximation for competing candidates.

    EVSI is expressed in the same utility units as replicate_scores. It is
    an estimate under this model, not an experimental result or an exact
    Bayesian decision value.
    """
    if future_replicates < 1:
        raise ValueError("future_replicates must be >=1")
    if simulations < 100:
        raise ValueError("simulations must be >=100")

    names = sorted(replicate_scores)
    if len(names) < 2:
        raise ValueError("EVSI requires at least two competing candidates")

    means, prior_vars, observation_vars = [], [], []
    for name in names:
        x = np.asarray(replicate_scores[name], dtype=float)
        if len(x) < 10:
            raise ValueError(f"{name} has fewer than 10 experimental replicates")
        if not np.isfinite(x).all():
            raise ValueError(f"{name} contains non-finite replicate scores")

        variance = float(x.var(ddof=1))
        means.append(float(x.mean()))
        prior_vars.append(variance / len(x))
        observation_vars.append(variance / future_replicates)

    means = np.asarray(means, dtype=float)
    prior_vars = np.asarray(prior_vars, dtype=float)
    observation_vars = np.asarray(observation_vars, dtype=float)

    rng = np.random.default_rng(seed)
    predictive_var = prior_vars + observation_vars
    predictive_draws = rng.normal(
        means,
        np.sqrt(np.maximum(predictive_var, 1e-15)),
        size=(simulations, len(names)),
    )

    gain = prior_vars / np.maximum(predictive_var, 1e-15)
    posterior_means = means + gain * (predictive_draws - means)
    posterior_maxima = posterior_means.max(axis=1)
    evsi = max(0.0, float(posterior_maxima.mean() - means.max()))

    posterior_samples = rng.normal(
        means,
        np.sqrt(np.maximum(prior_vars, 1e-15)),
        size=(simulations, len(names)),
    )
    winners = np.argmax(posterior_samples, axis=1)
    win_probabilities = {
        name: float(np.mean(winners == i))
        for i, name in enumerate(names)
    }

    return {
        "evsi_estimate": evsi,
        "future_replicates_per_candidate": int(future_replicates),
        "simulations": int(simulations),
        "candidate_win_probabilities": win_probabilities,
        "basis": "normal-normal approximation from experimental replicate scores",
        "utility_units": "same units as replicate_scores",
        "is_exact_zero_claim": False,
    }


def candidate_prune_decision(
    candidate_name: str,
    replicate_scores: dict[str, list[float]],
    evidence_type: str,
    reliability_calibrated: bool,
    evsi_tolerance: float = 1e-3,
    seed: int = 20261031,
    pilot_backbone_coordinates: Any | None = None,
    tau_cost: float | None = None,
    rmsd_radius_threshold_A: float = 1.2,
    additional_trajectories: int = 10,
) -> dict[str, Any]:
    """
    Conservative pruning gate.

    Pruning requires:
      - independent experimental replicate scores;
      - externally calibrated reliability;
      - >=10 experimental replicates per candidate;
      - >=10 independent post-equilibration structural trajectories;
      - mean backbone RMSD radius >= the configured threshold;
      - experimental interval dominance, low probability of being best,
        and EVSI <= tau_cost.

    pilot_backbone_coordinates contains one consistently atom-mapped
    backbone representative per independent trajectory. Missing trajectory
    evidence blocks pruning. A low RMSD radius requests more sampling.
    """
    if evidence_type != "experimental":
        return {
            "prune": False,
            "decision_eligible": False,
            "action": "collect_experimental_evidence",
            "reason": "Pruning requires independent experimental measurements",
        }

    if not reliability_calibrated:
        return {
            "prune": False,
            "decision_eligible": False,
            "action": "calibrate_reliability",
            "reason": "Reliability model is not externally calibrated",
        }

    if candidate_name not in replicate_scores:
        raise KeyError(candidate_name)
    if len(replicate_scores) < 2:
        raise ValueError("Pruning requires at least two competing candidates")
    if additional_trajectories < 1:
        raise ValueError("additional_trajectories must be >=1")
    if not np.isfinite(rmsd_radius_threshold_A) or rmsd_radius_threshold_A <= 0:
        raise ValueError("rmsd_radius_threshold_A must be positive and finite")

    for name, values in replicate_scores.items():
        x = np.asarray(values, dtype=float)
        if len(x) < 10:
            return {
                "prune": False,
                "decision_eligible": False,
                "action": "collect_experimental_replicates",
                "reason": f"{name} has fewer than 10 replicates",
            }
        if not np.isfinite(x).all():
            raise ValueError(f"{name} contains non-finite experimental scores")

    if pilot_backbone_coordinates is None:
        return {
            "prune": False,
            "decision_eligible": False,
            "action": "collect_initial_pilot_trajectories",
            "required_initial_trajectories": 10,
            "reason": "No independent trajectory backbone coordinates supplied",
        }

    coords = np.asarray(pilot_backbone_coordinates, dtype=float)
    if coords.ndim != 3 or coords.shape[2] != 3:
        raise ValueError(
            "pilot_backbone_coordinates must have shape (n_trajectories, n_atoms, 3)"
        )

    n_trajectories = int(coords.shape[0])
    if n_trajectories < 10:
        return {
            "prune": False,
            "decision_eligible": False,
            "action": "complete_initial_pilot_trajectories",
            "additional_trajectories": 10 - n_trajectories,
            "observed_trajectories": n_trajectories,
            "reason": "The initial RMSD gate requires 10 independent trajectories",
        }

    rmsd = backbone_rmsd_radius_A(coords)
    radius = rmsd["mean_backbone_rmsd_radius_A"]

    if radius < rmsd_radius_threshold_A:
        return {
            "prune": False,
            "decision_eligible": False,
            "action": "run_additional_sampling_trajectories",
            "additional_trajectories": int(additional_trajectories),
            "rmsd_radius": rmsd,
            "rmsd_radius_threshold_A": float(rmsd_radius_threshold_A),
            "reason": (
                "Backbone ensemble dispersion is below threshold; "
                "allocate additional trajectories before reevaluating"
            ),
        }

    effective_tau_cost = (
        float(evsi_tolerance) if tau_cost is None else float(tau_cost)
    )
    if not np.isfinite(effective_tau_cost) or effective_tau_cost < 0:
        raise ValueError("tau_cost must be finite and >=0")

    evsi = estimate_evsi(replicate_scores, seed=seed)
    candidate = np.asarray(replicate_scores[candidate_name], dtype=float)
    rivals = {
        k: np.asarray(v, dtype=float)
        for k, v in replicate_scores.items()
        if k != candidate_name
    }

    c_mean = float(candidate.mean())
    c_sem = float(candidate.std(ddof=1) / math.sqrt(len(candidate)))
    c_upper = c_mean + 1.96 * c_sem

    best_rival_lower = max(
        float(v.mean() - 1.96 * v.std(ddof=1) / math.sqrt(len(v)))
        for v in rivals.values()
    )

    win_probability = evsi["candidate_win_probabilities"][candidate_name]
    dominated = c_upper < best_rival_lower
    evsi_below_cost = evsi["evsi_estimate"] <= effective_tau_cost

    prune = bool(
        dominated
        and win_probability < 0.01
        and evsi_below_cost
    )

    return {
        "prune": prune,
        "decision_eligible": True,
        "action": "prune" if prune else "retain",
        "candidate_upper_95_approx": float(c_upper),
        "best_rival_lower_95_approx": float(best_rival_lower),
        "candidate_probability_of_being_best": float(win_probability),
        "evsi": evsi,
        "tau_cost": effective_tau_cost,
        "tau_cost_units": evsi["utility_units"],
        "evsi_below_tau_cost": bool(evsi_below_cost),
        "rmsd_radius": rmsd,
        "rmsd_radius_threshold_A": float(rmsd_radius_threshold_A),
        "reason": (
            "RMSD gate passed; experimental interval dominance, low "
            "probability of being best, and EVSI <= tau_cost"
            if prune else
            "RMSD gate passed, but experimental pruning criteria are not all met"
        ),
    }


# ---------------------------------------------------------------------
# PHASE A / PHASE B FAMILY SPLIT AND POWER ANALYSIS
# ---------------------------------------------------------------------

def validate_phase_split(
    phase_a: list[dict[str, Any]],
    phase_b: list[dict[str, Any]],
) -> dict[str, Any]:
    a_families = {str(x["family_id"]) for x in phase_a}
    b_families = {str(x["family_id"]) for x in phase_b}
    overlap = sorted(a_families & b_families)

    result = {
        "phase_a_complexes": len(phase_a),
        "phase_a_families": len(a_families),
        "phase_b_complexes": len(phase_b),
        "phase_b_families": len(b_families),
        "overlapping_families": overlap,
        "valid": (
            len(phase_a) >= 20
            and len(a_families) >= 10
            and len(phase_b) >= 40
            and len(b_families) >= 30
            and not overlap
        ),
    }
    return result


def family_bootstrap_power(
    family_effects: list[float],
    confirmatory_families: int = 30,
    bootstrap_replicates: int = 10000,
    alpha: float = 0.05,
    seed: int = 20261031,
    target_delta_rmse: float = 0.15,
    desired_power: float = 0.80,
    max_confirmatory_families: int = 300,
) -> dict[str, Any]:
    """
    Prospective power estimate for a paired family-level RMSE contrast.

    Each Phase A value must be:
        RMSE_baseline_for_family - RMSE_model_for_family

    Positive values mean the model reduces RMSE. The bootstrap resamples
    centered Phase A family-level contrasts, adds target_delta_rmse, and
    tests the simulated Phase B mean against zero using a two-sided
    one-sample t test.

    The Monte-Carlo estimate is not a guarantee. Families, not complexes
    or individual mutations, are the independent resampling units.
    """
    x = np.asarray(family_effects, dtype=float)

    if len(x) < 10 or not np.isfinite(x).all():
        raise ValueError(
            "Supply at least 10 finite Phase A family-level paired RMSE contrasts"
        )
    if confirmatory_families < 30:
        raise ValueError("Phase B requires at least 30 unseen families")
    if max_confirmatory_families < confirmatory_families:
        raise ValueError("max_confirmatory_families must be >= confirmatory_families")
    if bootstrap_replicates < 10000:
        raise ValueError("P4 confirmatory power analysis requires >=10,000 resamples")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    if not np.isfinite(target_delta_rmse) or target_delta_rmse <= 0:
        raise ValueError("target_delta_rmse must be positive and finite")
    if not 0 < desired_power < 1:
        raise ValueError("desired_power must be between 0 and 1")

    try:
        from scipy.stats import t as student_t
        from scipy.stats import nct
        scipy_available = True
    except ImportError:
        from statistics import NormalDist
        student_t = None
        nct = None
        scipy_available = False

    centered = x - x.mean()
    rng = np.random.default_rng(seed)

    # Phase B bootstrap distribution under the prespecified +0.15 RMSE effect.
    residual_draws = rng.choice(
        centered,
        size=(bootstrap_replicates, confirmatory_families),
        replace=True,
    )
    simulated = target_delta_rmse + residual_draws
    simulated_means = simulated.mean(axis=1)
    simulated_sds = simulated.std(axis=1, ddof=1)

    if scipy_available:
        critical = float(
            student_t.ppf(1 - alpha / 2, df=confirmatory_families - 1)
        )
    else:
        critical = float(NormalDist().inv_cdf(1 - alpha / 2))

    standard_errors = simulated_sds / math.sqrt(confirmatory_families)
    statistics = np.divide(
        np.abs(simulated_means),
        standard_errors,
        out=np.full(bootstrap_replicates, np.inf),
        where=standard_errors > 1e-12,
    )
    significant = statistics > critical
    power = float(significant.mean())

    # Wilson interval for Monte-Carlo uncertainty in estimated power.
    successes = int(significant.sum())
    z = 1.959963984540054
    denominator = 1 + z * z / bootstrap_replicates
    center = (
        power + z * z / (2 * bootstrap_replicates)
    ) / denominator
    half = (
        z * math.sqrt(
            power * (1 - power) / bootstrap_replicates
            + z * z / (4 * bootstrap_replicates ** 2)
        ) / denominator
    )
    power_ci = [
        float(max(0.0, center - half)),
        float(min(1.0, center + half)),
    ]

    # Approximate sample-size search using the noncentral-t distribution
    # and the Phase A family-level SD. This is a planning approximation,
    # separate from the bootstrap estimate at confirmatory_families.
    required_n = None
    pilot_sd = float(x.std(ddof=1))
    if pilot_sd <= 1e-12:
        required_n = 30
    else:
        for n in range(30, max_confirmatory_families + 1):
            if scipy_available:
                crit_n = float(student_t.ppf(1 - alpha / 2, df=n - 1))
                ncp = target_delta_rmse * math.sqrt(n) / pilot_sd
                approx_power = float(
                    nct.sf(crit_n, df=n - 1, nc=ncp)
                    + nct.cdf(-crit_n, df=n - 1, nc=ncp)
                )
            else:
                crit_n = float(NormalDist().inv_cdf(1 - alpha / 2))
                ncp = target_delta_rmse * math.sqrt(n) / pilot_sd
                approx_power = float(
                    NormalDist().cdf(-crit_n - ncp)
                    + 1 - NormalDist().cdf(crit_n - ncp)
                )
            if approx_power >= desired_power:
                required_n = n
                break

    return {
        "pilot_families": int(len(x)),
        "confirmatory_families": int(confirmatory_families),
        "target_delta_rmse_pKd": float(target_delta_rmse),
        "pilot_mean_paired_delta_rmse": float(x.mean()),
        "pilot_sd_paired_delta_rmse": pilot_sd,
        "estimated_power": power,
        "power_ci95_wilson": power_ci,
        "power_lower_ci_at_least_80_percent": bool(
            power_ci[0] >= desired_power
        ),
        "power_at_least_80_percent": bool(power >= desired_power),
        "alpha": float(alpha),
        "bootstrap_replicates": int(bootstrap_replicates),
        "bootstrap_successes": successes,
        "recommended_families_approx_for_target_power": required_n,
        "max_families_searched": int(max_confirmatory_families),
        "method": (
            "nonparametric resampling of centered Phase A family-level "
            "paired RMSE contrasts with target-effect shift; two-sided "
            "one-sample t test"
        ),
        "planning_sample_size_method": (
            "noncentral-t approximation using Phase A family-level SD"
        ),
        "scipy_student_t_available": scipy_available,
        "interpretation": (
            "Prospective power estimate, not a guarantee. Independent "
            "resampling unit is the CATH family."
        ),
    }


# ---------------------------------------------------------------------
# SELF-TESTS
# ---------------------------------------------------------------------

def self_test() -> None:
    frames = [
        {"A-B", "B-C"},
        {"A-B", "B-C"},
        {"A-B", "B-C"},
        {"X-Y"},
        {"X-Y"},
    ]
    result = contact_entropy(frames)
    assert result["n_clusters"] == 2
    assert result["cluster_sizes"] == [3, 2]

    identical = contact_entropy([{"A-B"}] * 5)
    assert identical["normalized_entropy"] == 0.0

    runs = [
        [{"A-B"}, {"A-B"}, {"A-B"}],
        [{"A-B"}, {"A-B"}, {"A-B"}],
    ]
    reliability = cross_run_reliability(runs)
    assert reliability["mean_cross_run_jaccard_similarity"] == 1.0
    assert reliability["calibrated"] is False

    candidates = [
        {
            "sequence": "AAAA",
            "sequence_score": 0.8,
            "target_geometry_score": 0.7,
            "arg203_min_A": 3.8,
            "arg205_min_A": 4.2,
            "severe_clash": False,
        },
        {
            "sequence": "BBBB",
            "sequence_score": 0.6,
            "target_geometry_score": 0.5,
            "arg203_min_A": 9.0,
            "arg205_min_A": 10.0,
            "severe_clash": False,
        },
    ]
    ranked = rank_target_candidates(candidates)
    assert len(ranked) == 2
    assert all("pareto_front" in x for x in ranked)

    weights = parent_selection_weights(
        [{"predicted_value": 0.9}, {"predicted_value": 0.1}]
    )
    assert abs(sum(weights) - 1.0) < 1e-10
    assert abs(min(weights) >= 0.075)

    gate = operator_gate(
        out_of_domain=True,
        normalized_entropy=0.9,
        cross_run_similarity=0.1,
        reliability_calibrated=False,
    )
    assert gate["aggressive_operators_enabled"] is True

    print("P4 ENGINE SELF-TESTS PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description="Silicon Virovore P4 engine")
    sub = parser.add_subparsers(dest="command", required=True)

    dbp = sub.add_parser("init-db")
    dbp.add_argument(
        "--db",
        default="prototype4/data/p4_closed_loop_registry.sqlite",
    )

    mp = sub.add_parser("seal-manifest")
    mp.add_argument("--db", default="prototype4/data/p4_closed_loop_registry.sqlite")
    mp.add_argument("--json", required=True)
    mp.add_argument(
        "--key",
        default="prototype4/configs/validation_manifest_hmac.key",
    )

    vp = sub.add_parser("verify-manifest")
    vp.add_argument("--db", default="prototype4/data/p4_closed_loop_registry.sqlite")
    vp.add_argument("--manifest-id", required=True)
    vp.add_argument(
        "--key",
        default="prototype4/configs/validation_manifest_hmac.key",
    )

    sub.add_parser("self-test")

    args = parser.parse_args()

    if args.command == "init-db":
        con = open_registry(args.db)
        tables = [
            row[0] for row in con.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            )
        ]
        con.close()
        print("Registry initialized:", args.db)
        print("Tables:", tables)

    elif args.command == "seal-manifest":
        manifest = json.loads(Path(args.json).read_text())
        con = open_registry(args.db)
        manifest_id = seal_manifest(con, manifest, args.key)
        con.close()
        print("Manifest sealed:", manifest_id)

    elif args.command == "verify-manifest":
        con = open_registry(args.db)
        valid = verify_manifest(con, args.manifest_id, args.key)
        con.close()
        print("Manifest valid:", valid)
        if not valid:
            raise SystemExit(1)

    elif args.command == "self-test":
        self_test()


if __name__ == "__main__":
    main()
