#!/usr/bin/env python3
"""Idempotent migration of the existing P4 Model A evidence into the registry."""

from __future__ import annotations

import hashlib
import json
import math
import sqlite3
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DB = ROOT / "prototype4/results/skempi_audit_featurefix_20261030/p4_evidence.sqlite"
MODEL_DB = ROOT / "prototype4/results/model_a_validation_20261030/p4_model_a_working.sqlite"
REGISTRY_DB = ROOT / "prototype4/data/p4_closed_loop_registry.sqlite"

STATE_MAP = {"P0": "WT", "Pi": "M1", "Pj": "M2", "Pij": "M1M2"}
MODEL_VERSION = "P4_ModelA_Ridge_group_bootstrap_v1"
BASELINE_VERSION = "P4_ModelA_fold_mean_baseline_v1"


def fail(message: str) -> None:
    raise RuntimeError(message)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def load_json(value: Any, label: str) -> Any:
    if isinstance(value, (dict, list)):
        return value
    try:
        return json.loads(value)
    except Exception as exc:
        fail(f"Cannot parse {label} as JSON: {exc}")


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    )


def table_columns(con: sqlite3.Connection, table: str) -> set[str]:
    return {
        row[1]
        for row in con.execute(f'PRAGMA table_info("{table}")').fetchall()
    }


def require_schema(con: sqlite3.Connection) -> None:
    required = {
        "skempi_raw": {
            "raw_id", "source_version", "source_record_key", "original_line",
            "chain_partner_1", "chain_partner_2", "experimental_method",
            "raw_sha256",
        },
        "cycle_nodes": {
            "node_id", "cycle_uid", "state", "raw_id",
            "mutation_partner_1", "mutation_partner_2",
            "mutation_1", "mutation_2", "topology", "interface_proximity",
        },
        "cycle_features": {
            "cycle_uid", "feature_name", "feature_value", "feature_json",
            "source_sha256",
        },
        "loco_predictions": {
            "model_version", "cycle_uid", "outer_fold", "heldout_group",
            "y_true", "y_pred", "train_groups_sha256", "feature_set",
        },
        "validation_manifest": {
            "manifest_id", "manifest_json", "content_sha256",
            "signature_hmac_sha256", "signature_key_id",
        },
        "structure_runs": {
            "run_id", "target_key", "run_role", "software",
            "software_version", "command_string", "config_sha256",
            "random_seed", "input_sha256", "output_sha256",
            "trajectory_path", "l_rmsd", "status", "started_at",
            "completed_at", "metadata_json",
        },
    }
    for table, columns in required.items():
        actual = table_columns(con, table)
        if not actual:
            fail(f"Registry is missing required table: {table}")
        missing = columns - actual
        if missing:
            fail(f"{table} is missing columns: {sorted(missing)}")


def parse_mutation(token: str) -> tuple[str, str]:
    parts = str(token).split(":")
    if len(parts) != 4:
        fail(f"Unexpected mutation token: {token!r}")
    partner, wild_type, position, mutant = parts
    if not partner or not wild_type or not position or not mutant:
        fail(f"Incomplete mutation token: {token!r}")
    return partner, f"{wild_type}{position}{mutant}"


def mutation_fields(node: sqlite3.Row) -> tuple:
    mutations = load_json(
        node["mutation_set_json"], f"mutation set for {node['node_uid']}"
    )
    if not isinstance(mutations, list):
        fail(f"Mutation set is not a list for {node['node_uid']}")
    if len(mutations) > 2:
        fail(
            f"{node['node_uid']} has {len(mutations)} mutations; "
            "registry schema only has two mutation slots. Refusing truncation."
        )

    parsed = [parse_mutation(x) for x in mutations]
    p1 = parsed[0][0] if len(parsed) >= 1 else None
    p2 = parsed[1][0] if len(parsed) >= 2 else None
    m1 = parsed[0][1] if len(parsed) >= 1 else None
    m2 = parsed[1][1] if len(parsed) >= 2 else None

    topology = "unresolved"
    if len(parsed) == 2:
        partners = {p1, p2}
        if partners == {"A"}:
            topology = "same_partner_1"
        elif partners == {"T"}:
            topology = "same_partner_2"
        elif partners == {"A", "T"}:
            topology = "cross_partner_1x2"

    # Interface proximity is not encoded robustly in this source dataset.
    # Preserve it as unresolved rather than inferring it from mutation labels.
    return p1, p2, m1, m2, topology, "unresolved"


def extract_feature_pairs(feature: sqlite3.Row) -> list[tuple[str, float]]:
    names = load_json(feature["feature_names_json"], "feature_names_json")
    values = load_json(feature["feature_values_json"], "feature_values_json")

    if not isinstance(names, list):
        fail(f"Feature names are not a list for {feature['cycle_uid']}")

    if isinstance(values, dict) and "values" in values:
        values = values["values"]

    if isinstance(values, dict):
        missing = [name for name in names if name not in values]
        if missing:
            fail(
                f"{feature['cycle_uid']} is missing {len(missing)} "
                f"feature values; first missing: {missing[:5]}"
            )
        pairs = [(name, values[name]) for name in names]
    elif isinstance(values, list):
        if len(names) != len(values):
            fail(
                f"Feature length mismatch for {feature['cycle_uid']}: "
                f"{len(names)} names vs {len(values)} values"
            )
        pairs = list(zip(names, values))
    else:
        fail(f"Unsupported feature-values format for {feature['cycle_uid']}")

    result = []
    for name, value in pairs:
        numeric = float(value)
        if not math.isfinite(numeric):
            fail(f"Non-finite feature {name} in {feature['cycle_uid']}")
        result.append((str(name), numeric))
    return result


def migrate() -> None:
    for path in (SOURCE_DB, MODEL_DB, REGISTRY_DB):
        if not path.is_file():
            fail(f"Required database not found: {path}")

    source = sqlite3.connect(f"file:{SOURCE_DB}?mode=ro", uri=True)
    model = sqlite3.connect(f"file:{MODEL_DB}?mode=ro", uri=True)
    source.row_factory = sqlite3.Row
    model.row_factory = sqlite3.Row

    registry = sqlite3.connect(REGISTRY_DB)
    registry.row_factory = sqlite3.Row
    registry.execute("PRAGMA foreign_keys=ON")
    registry.execute("PRAGMA busy_timeout=30000")

    try:
        require_schema(registry)

        source_raw = source.execute(
            "SELECT * FROM skempi_raw ORDER BY source_row_number"
        ).fetchall()
        source_nodes = source.execute(
            "SELECT * FROM cycle_nodes ORDER BY cycle_uid, node_state"
        ).fetchall()
        source_features = source.execute(
            "SELECT * FROM cycle_features ORDER BY cycle_uid"
        ).fetchall()
        predictions = model.execute(
            "SELECT * FROM loco_predictions ORDER BY fold_index, feature_uid"
        ).fetchall()
        source_manifest_row = model.execute(
            "SELECT * FROM validation_manifest LIMIT 1"
        ).fetchone()

        if not source_raw or not source_nodes or not source_features:
            fail("Source evidence database is missing required data")
        if not predictions or source_manifest_row is None:
            fail("Model A database lacks predictions or validation manifest")

        source_hashes = {row["source_sha256"] for row in source_raw}
        if len(source_hashes) != 1:
            fail(f"Expected one source dataset hash; found {len(source_hashes)}")
        dataset_sha256 = next(iter(source_hashes))

        source_manifest = dict(source_manifest_row)
        if source_manifest["input_dataset_sha256"] != dataset_sha256:
            fail("Manifest dataset hash does not match raw source rows")
        if int(source_manifest["cycle_count"]) != len(source_features):
            fail("Manifest cycle count does not match source feature rows")

        feature_to_cycle = {
            row["feature_uid"]: row["cycle_uid"] for row in source_features
        }

        # Each cycle must have exactly one of each of the four states.
        states_by_cycle: dict[str, set[str]] = {}
        for node in source_nodes:
            states_by_cycle.setdefault(node["cycle_uid"], set()).add(
                node["node_state"]
            )
        expected_states = set(STATE_MAP)
        bad_cycles = {
            cycle: sorted(states)
            for cycle, states in states_by_cycle.items()
            if states != expected_states
        }
        if bad_cycles:
            fail(
                f"{len(bad_cycles)} cycles do not have all four states; "
                f"examples: {list(bad_cycles.items())[:3]}"
            )
        if len(states_by_cycle) != len(source_features):
            fail("Node-cycle count does not equal feature-cycle count")

        # Verify source/model node and feature-cycle identity before writing.
        model_cycles = {
            row["cycle_uid"]
            for row in model.execute("SELECT cycle_uid FROM cycle_features")
        }
        if model_cycles != set(states_by_cycle):
            fail("Source and Model A working databases have different cycle sets")

        registry.execute("BEGIN IMMEDIATE")

        # 1. Raw records: preserve exact JSON payload and dataset-level SHA
        # in the imported manifest. Record hash is deterministic over
        # source_row_uid + newline + exact raw_json payload.
        raw_id_by_uid: dict[str, int] = {}
        for row in source_raw:
            uid = row["source_row_uid"]
            raw_json = row["raw_json"]
            row_hash = sha256_text(uid + "\n" + raw_json)

            existing = registry.execute(
                "SELECT raw_id, original_line, raw_sha256 "
                "FROM skempi_raw WHERE source_record_key=?",
                (uid,),
            ).fetchall()

            if len(existing) > 1:
                fail(f"Duplicate registry source_record_key: {uid}")

            if existing:
                record = existing[0]
                if (
                    record["original_line"] != raw_json
                    or record["raw_sha256"] != row_hash
                ):
                    fail(f"Existing raw record conflicts with source: {uid}")
                raw_id = int(record["raw_id"])
            else:
                raw = load_json(raw_json, f"raw_json {uid}")
                if not isinstance(raw, dict):
                    fail(f"Raw source row {uid} is not a JSON object")

                cur = registry.execute(
                    """
                    INSERT INTO skempi_raw (
                        source_version, source_record_key, original_line,
                        chain_partner_1, chain_partner_2, experimental_method,
                        raw_sha256
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        str(raw.get("SKEMPI version", "unknown")),
                        uid,
                        raw_json,
                        str(raw.get("Protein 1", "")),
                        str(raw.get("Protein 2", "")),
                        str(raw.get("Method", "")),
                        row_hash,
                    ),
                )
                raw_id = int(cur.lastrowid)

            raw_id_by_uid[uid] = raw_id

        # 2. Four-state nodes.
        for node in source_nodes:
            state = STATE_MAP.get(node["node_state"])
            if state is None:
                fail(f"Unknown source node state: {node['node_state']}")

            raw_id = raw_id_by_uid.get(node["source_row_uid"])
            if raw_id is None:
                fail(f"No imported raw row for node {node['node_uid']}")

            p1, p2, m1, m2, topology, proximity = mutation_fields(node)

            existing = registry.execute(
                "SELECT * FROM cycle_nodes WHERE node_id=?",
                (node["node_uid"],),
            ).fetchone()

            values = (
                node["cycle_uid"], state, raw_id, p1, p2, m1, m2,
                topology, proximity,
            )
            if existing:
                observed = (
                    existing["cycle_uid"], existing["state"], existing["raw_id"],
                    existing["mutation_partner_1"], existing["mutation_partner_2"],
                    existing["mutation_1"], existing["mutation_2"],
                    existing["topology"], existing["interface_proximity"],
                )
                if observed != values:
                    fail(f"Existing cycle node conflicts: {node['node_uid']}")
            else:
                registry.execute(
                    """
                    INSERT INTO cycle_nodes (
                        node_id, cycle_uid, state, raw_id,
                        mutation_partner_1, mutation_partner_2,
                        mutation_1, mutation_2, topology, interface_proximity
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (node["node_uid"], *values),
                )

        # 3. Per-cycle descriptors. Each descriptor is a separate feature row.
        feature_rows_expected = 0
        for feature in source_features:
            pairs = extract_feature_pairs(feature)
            metadata = {
                "source_feature_uid": feature["feature_uid"],
                "family_id": feature["family_id"],
                "complex_id": feature["complex_id"],
                "source_node_ids": {
                    "WT": feature["node_uid_p0"],
                    "M1": feature["node_uid_pi"],
                    "M2": feature["node_uid_pj"],
                    "M1M2": feature["node_uid_pij"],
                },
                "source_dataset_sha256": dataset_sha256,
                "feature_definition": "fixed amino-acid substitution/property descriptors",
            }

            for index, (name, value) in enumerate(pairs):
                feature_json = canonical_json({
                    **metadata,
                    "feature_index": index,
                })
                existing = registry.execute(
                    """
                    SELECT feature_value, feature_json, source_sha256
                    FROM cycle_features
                    WHERE cycle_uid=? AND feature_name=?
                    """,
                    (feature["cycle_uid"], name),
                ).fetchone()

                expected = (value, feature_json, dataset_sha256)
                if existing:
                    observed = (
                        float(existing["feature_value"]),
                        existing["feature_json"],
                        existing["source_sha256"],
                    )
                    if observed != expected:
                        fail(
                            f"Existing feature conflicts: "
                            f"{feature['cycle_uid']} / {name}"
                        )
                else:
                    registry.execute(
                        """
                        INSERT INTO cycle_features (
                            cycle_uid, feature_name, feature_value,
                            feature_json, source_sha256
                        ) VALUES (?, ?, ?, ?, ?)
                        """,
                        (feature["cycle_uid"], name, value,
                         feature_json, dataset_sha256),
                    )
                feature_rows_expected += 1

        # 4. Held-out predictions: primary Ridge and fold-mean baseline.
        # Store calibration/uncertainty metadata in feature_set as JSON.
        predictions_imported = 0
        for prediction in predictions:
            cycle_uid = feature_to_cycle.get(prediction["feature_uid"])
            if cycle_uid is None:
                fail(
                    "LOCO prediction references an unknown feature UID: "
                    f"{prediction['feature_uid']}"
                )

            common = (
                cycle_uid,
                int(prediction["fold_index"]),
                str(prediction["heldout_group"]),
                float(prediction["y_true"]),
                str(prediction["train_group_hash"]),
            )

            variants = [
                (
                    MODEL_VERSION,
                    float(prediction["y_pred"]),
                    {
                        "role": "primary_model_prediction",
                        "feature_definition": "410 fixed amino-acid substitution/property descriptors",
                        "calibration_status": prediction["calibration_status"],
                        "raw_uncertainty": float(prediction["raw_uncertainty"]),
                        "calibrated_abs_error": float(
                            prediction["calibrated_abs_error"]
                        ),
                    },
                ),
                (
                    BASELINE_VERSION,
                    float(prediction["baseline_pred"]),
                    {
                        "role": "training_fold_mean_baseline",
                        "calibration_status": "not_applicable",
                    },
                ),
            ]

            for version, y_pred, feature_set in variants:
                existing = registry.execute(
                    """
                    SELECT heldout_group, y_true, y_pred,
                           train_groups_sha256, feature_set
                    FROM loco_predictions
                    WHERE model_version=? AND cycle_uid=? AND outer_fold=?
                    """,
                    (version, cycle_uid, common[1]),
                ).fetchone()

                values = (
                    common[2], common[3], y_pred, common[4],
                    canonical_json(feature_set),
                )
                if existing:
                    observed = (
                        existing["heldout_group"],
                        float(existing["y_true"]),
                        float(existing["y_pred"]),
                        existing["train_groups_sha256"],
                        existing["feature_set"],
                    )
                    if observed != values:
                        fail(
                            f"Existing prediction conflicts: {version} / "
                            f"{cycle_uid}"
                        )
                else:
                    registry.execute(
                        """
                        INSERT INTO loco_predictions (
                            model_version, cycle_uid, outer_fold, heldout_group,
                            y_true, y_pred, train_groups_sha256, feature_set
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            version, common[0], common[1], common[2],
                            common[3], y_pred, common[4],
                            canonical_json(feature_set),
                        ),
                    )
                predictions_imported += 1

        # 5. Import the original manifest as provenance. Do not invent an
        # HMAC signature for an artifact that was not signed with this key.
        imported_manifest = {
            "source_validation_manifest": source_manifest,
            "migration": {
                "source_evidence_database": str(SOURCE_DB.relative_to(ROOT)),
                "source_model_database": str(MODEL_DB.relative_to(ROOT)),
                "source_dataset_sha256": dataset_sha256,
                "raw_record_hash_algorithm":
                    "sha256(source_row_uid + newline + exact raw_json)",
                "node_state_mapping": STATE_MAP,
                "feature_storage": "one registry row per cycle and descriptor",
                "prediction_versions": [MODEL_VERSION, BASELINE_VERSION],
                "signature_claim": "none; source manifest integrity hash preserved",
                "importer": "prototype4/import_p4_registry.py",
            },
        }
        manifest_json = canonical_json(imported_manifest)
        manifest_hash = sha256_text(manifest_json)
        manifest_id = "imported_model_a_" + str(source_manifest["manifest_uid"])

        existing_manifest = registry.execute(
            "SELECT manifest_json, content_sha256 FROM validation_manifest "
            "WHERE manifest_id=?",
            (manifest_id,),
        ).fetchone()
        if existing_manifest:
            if (
                existing_manifest["manifest_json"] != manifest_json
                or existing_manifest["content_sha256"] != manifest_hash
            ):
                fail("Imported Model A manifest conflicts with existing registry row")
        else:
            registry.execute(
                """
                INSERT INTO validation_manifest (
                    manifest_id, manifest_json, content_sha256,
                    signature_hmac_sha256, signature_key_id
                ) VALUES (?, ?, ?, NULL, NULL)
                """,
                (manifest_id, manifest_json, manifest_hash),
            )

        # Validate before committing.
        cycle_count = registry.execute(
            "SELECT COUNT(DISTINCT cycle_uid) FROM cycle_nodes"
        ).fetchone()[0]
        node_count = registry.execute(
            "SELECT COUNT(*) FROM cycle_nodes"
        ).fetchone()[0]
        raw_count = registry.execute(
            "SELECT COUNT(*) FROM skempi_raw"
        ).fetchone()[0]
        prediction_count = registry.execute(
            "SELECT COUNT(*) FROM loco_predictions "
            "WHERE model_version IN (?, ?)",
            (MODEL_VERSION, BASELINE_VERSION),
        ).fetchone()[0]
        feature_count = registry.execute(
            "SELECT COUNT(*) FROM cycle_features"
        ).fetchone()[0]

        if raw_count < len(source_raw):
            fail(f"Raw import incomplete: {raw_count} / {len(source_raw)}")
        if cycle_count < len(source_features):
            fail(
                f"Cycle import incomplete: {cycle_count} / {len(source_features)}"
            )
        if node_count < len(source_nodes):
            fail(f"Node import incomplete: {node_count} / {len(source_nodes)}")
        if prediction_count < 2 * len(predictions):
            fail(
                f"Prediction import incomplete: {prediction_count} / "
                f"{2 * len(predictions)}"
            )
        if feature_count < feature_rows_expected:
            fail(
                f"Feature import incomplete: {feature_count} / "
                f"{feature_rows_expected}"
            )

        fk_errors = registry.execute("PRAGMA foreign_key_check").fetchall()
        if fk_errors:
            fail(f"Registry foreign-key violations: {fk_errors[:5]}")

        registry.commit()

        print("\nP4 REGISTRY IMPORT SUCCESS")
        print(f"Source dataset SHA-256: {dataset_sha256}")
        print(f"Raw records in source: {len(source_raw)}")
        print(f"Cycle nodes in source: {len(source_nodes)}")
        print(f"Four-state cycles: {len(source_features)}")
        print(f"Descriptor rows expected/imported: {feature_rows_expected}")
        print(f"LOCO source predictions: {len(predictions)}")
        print(f"Model + baseline prediction rows: {2 * len(predictions)}")
        print(f"Registry raw records: {raw_count}")
        print(f"Registry cycle nodes: {node_count}")
        print(f"Registry descriptor rows: {feature_count}")
        print(f"Registry Model A prediction rows: {prediction_count}")
        print(f"Imported manifest: {manifest_id}")
        print("Foreign-key check: PASS")
        print("HMAC claim: none")

    except Exception:
        registry.rollback()
        raise
    finally:
        source.close()
        model.close()
        registry.close()


if __name__ == "__main__":
    migrate()
