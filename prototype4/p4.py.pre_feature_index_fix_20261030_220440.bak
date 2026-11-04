#!/usr/bin/env python3
"""Silicon Virovore P4: provenance-first SKEMPI mutation-cycle benchmark.

This module builds experimentally grounded double-mutant cycles from exact,
uncensored SKEMPI Kd measurements; stores immutable source provenance in SQLite;
performs leave-one-family-out prediction with all learned operations nested inside
each outer training split; and evaluates uncertainty-based abstention.

It does NOT run docking, infer HERV-K binding, or modify the ARISE/ALE optimizer.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import itertools
import json
import math
import re
import sqlite3
import statistics
import sys
import time
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold

PROJECT_ROOT = Path(__file__).resolve().parents[1]
R_KCAL = 0.00198720425864083
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA)}
NUMERIC_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
MUTATION_RE = re.compile(r"^([A-Z])([A-Za-z0-9])(-?\d+[A-Za-z]?)([A-Z])$")

# Fixed descriptors. Values are deterministic physicochemical lookup features,
# not learned from the experimental benchmark.
AA_PROP = {
    # hydropathy (KD), approximate side-chain volume, formal charge,
    # aromatic, polar, sulfur-containing, proline, glycine
    "A": (1.8,  88.6,  0, 0, 0, 0, 0, 0),
    "C": (2.5, 108.5,  0, 0, 1, 1, 0, 0),
    "D": (-3.5,111.1, -1, 0, 1, 0, 0, 0),
    "E": (-3.5,138.4, -1, 0, 1, 0, 0, 0),
    "F": (2.8, 189.9,  0, 1, 0, 0, 0, 0),
    "G": (-0.4, 60.1,  0, 0, 0, 0, 0, 1),
    "H": (-3.2,153.2,  0, 1, 1, 0, 0, 0),
    "I": (4.5, 166.7,  0, 0, 0, 0, 0, 0),
    "K": (-3.9,168.6,  1, 0, 1, 0, 0, 0),
    "L": (3.8, 166.7,  0, 0, 0, 0, 0, 0),
    "M": (1.9, 162.9,  0, 0, 0, 1, 0, 0),
    "N": (-3.5,114.1,  0, 0, 1, 0, 0, 0),
    "P": (-1.6,112.7,  0, 0, 0, 0, 1, 0),
    "Q": (-3.5,143.8,  0, 0, 1, 0, 0, 0),
    "R": (-4.5,173.4,  1, 0, 1, 0, 0, 0),
    "S": (-0.8, 89.0,  0, 0, 1, 0, 0, 0),
    "T": (-0.7,116.1,  0, 0, 1, 0, 0, 0),
    "V": (4.2, 140.0,  0, 0, 0, 0, 0, 0),
    "W": (-0.9,227.8,  0, 1, 1, 0, 0, 0),
    "Y": (-1.3,193.6,  0, 1, 1, 0, 0, 0),
}
LOCATION_NAMES = ("core", "rim", "support", "surface")


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_uid(*parts: Any, n: int = 24) -> str:
    payload = "\x1f".join(str(x) for x in parts)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:n]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def norm_text(value: Any) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def resolve_column(columns: Iterable[str], aliases: Iterable[str], label: str,
                   required: bool = True) -> str | None:
    by_norm = {re.sub(r"\s+", " ", str(c).strip().casefold()): str(c) for c in columns}
    for alias in aliases:
        key = re.sub(r"\s+", " ", alias.strip().casefold())
        if key in by_norm:
            return by_norm[key]
    if required:
        raise ValueError(
            f"Missing required SKEMPI column for {label}. Tried {list(aliases)}; "
            f"available columns: {list(columns)}"
        )
    return None


def detect_delimiter(path: str | Path) -> str:
    """Detect comma/semicolon/tab from the actual header; fail if ambiguous."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(65536)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t|")
        delim = dialect.delimiter
    except csv.Error:
        header = sample.splitlines()[0] if sample.splitlines() else ""
        counts = {d: header.count(d) for d in ",;\t|"}
        delim = max(counts, key=counts.get)
        if counts[delim] == 0:
            raise ValueError("Could not detect the source delimiter from the header.")
    # Cross-check that the first parsed line has many columns, not one giant cell.
    header_cols = next(csv.reader([sample.splitlines()[0]], delimiter=delim), [])
    if len(header_cols) < 10:
        raise ValueError(f"Delimiter detection produced only {len(header_cols)} columns ({delim!r}).")
    return delim


def parse_exact_positive_number(raw: Any, parsed_fallback: Any = "") -> float | None:
    """Return a positive exact numeric value. Inequalities/censoring are rejected."""
    raw_s = norm_text(raw)
    if raw_s:
        # Never accept a parsed numeric fallback when the original is censored.
        if not NUMERIC_RE.fullmatch(raw_s):
            return None
        try:
            value = float(raw_s)
        except ValueError:
            return None
    else:
        parsed_s = norm_text(parsed_fallback)
        if not parsed_s or not NUMERIC_RE.fullmatch(parsed_s):
            return None
        try:
            value = float(parsed_s)
        except ValueError:
            return None
    if not math.isfinite(value) or value <= 0:
        return None
    return value


def parse_mutation_set(raw: Any) -> tuple[tuple[str, str, str, str], ...] | None:
    """Parse SKEMPI cleaned/PDB notation e.g. YA52F or YA52F,GB17A."""
    text = norm_text(raw).replace(" ", "")
    if not text:
        return None
    chunks = [x for x in re.split(r"[,;]", text) if x]
    parsed: list[tuple[str, str, str, str]] = []
    for chunk in chunks:
        match = MUTATION_RE.fullmatch(chunk)
        if not match:
            return None
        wt, chain, residue, mutant = match.groups()
        if wt not in AA_INDEX or mutant not in AA_INDEX or wt == mutant:
            return None
        parsed.append((wt, chain, residue, mutant))
    # One physical site cannot be mutated twice in a state description.
    sites = [(m[1], m[2]) for m in parsed]
    if len(set(sites)) != len(sites):
        return None
    return tuple(sorted(parsed, key=lambda m: (m[1], m[2], m[0], m[3])))


def mutation_token(m: tuple[str, str, str, str]) -> str:
    return f"{m[0]}:{m[1]}:{m[2]}:{m[3]}"


def parse_location_list(raw: Any, n_mutations: int) -> list[str]:
    text = norm_text(raw).casefold()
    items = [norm_text(x).casefold() for x in re.split(r"[,;|]", text) if norm_text(x)]
    # The source may provide one class for a combined multi-mutant row or a
    # semicolon-delimited class per mutation. We only retain an aligned label
    # when the count matches; otherwise location is explicitly unknown.
    if len(items) != n_mutations:
        return ["unknown"] * n_mutations
    return [x if x in LOCATION_NAMES else "unknown" for x in items]


def read_skempi(path: str | Path) -> tuple[pd.DataFrame, str, str]:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"SKEMPI file not found: {path}")
    delimiter = detect_delimiter(path)
    df = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False,
                     encoding="utf-8-sig", na_filter=False)
    if len(df.columns) < 10:
        raise ValueError(f"Parsed only {len(df.columns)} columns; source delimiter may be wrong.")
    source_hash = sha256_file(path)
    return df, delimiter, source_hash


def build_column_map(df: pd.DataFrame) -> dict[str, str | None]:
    c = df.columns
    return {
        "pdb": resolve_column(c, ["#Pdb", "Pdb", "PDB"], "complex/PDB identifier"),
        "mutations": resolve_column(c, ["Mutation(s)_cleaned", "Mutation(s)_PDB", "Mutations"], "mutation string"),
        "aff_wt": resolve_column(c, ["Affinity_wt (M)", "Affinity_wt", "Kd_wt"], "wild-type Kd"),
        "aff_mut": resolve_column(c, ["Affinity_mut (M)", "Affinity_mut", "Kd_mut"], "mutant Kd"),
        "protein1": resolve_column(c, ["Protein 1", "Protein1"], "protein 1"),
        "protein2": resolve_column(c, ["Protein 2", "Protein2"], "protein 2"),
        "reference": resolve_column(c, ["Reference", "Ref"], "experimental reference"),
        "temperature": resolve_column(c, ["Temperature", "Temp (K)", "Temperature (K)"], "temperature in kelvin"),
        "method": resolve_column(c, ["Method", "Experimental method", "Affinity method"], "measurement method", required=False),
        "location": resolve_column(c, ["iMutation_Location(s)", "Mutation_Location(s)", "Mutation location"], "mutation location", required=False),
        "aff_wt_parsed": resolve_column(c, ["Affinity_wt_parsed", "Kd_wt_parsed"], "parsed wild-type Kd", required=False),
        "aff_mut_parsed": resolve_column(c, ["Affinity_mut_parsed", "Kd_mut_parsed"], "parsed mutant Kd", required=False),
        "holdout_type": resolve_column(c, ["Hold_out_type"], "SKEMPI holdout type", required=False),
        "holdout_proteins": resolve_column(c, ["Hold_out_proteins"], "SKEMPI holdout proteins", required=False),
    }


def normalized_context(row: dict[str, Any], cols: dict[str, str | None]) -> dict[str, Any] | None:
    pdb = norm_text(row.get(cols["pdb"], "")).upper()
    p1 = norm_text(row.get(cols["protein1"], ""))
    p2 = norm_text(row.get(cols["protein2"], ""))
    reference = norm_text(row.get(cols["reference"], ""))
    temp_s = norm_text(row.get(cols["temperature"], ""))
    try:
        temp_k = float(temp_s)
    except ValueError:
        return None
    # SKEMPI stores temperature in kelvin. Reject implausible/missing values.
    if not math.isfinite(temp_k) or not 250.0 <= temp_k <= 400.0:
        return None
    if not all([pdb, p1, p2, reference]):
        return None
    method_col = cols.get("method")
    method = norm_text(row.get(method_col, "")) if method_col else ""
    method = method or "NOT_REPORTED"
    # Exact context matching: do not silently combine papers, temperatures,
    # chain pairs, or methods. Numeric temperatures are canonicalized.
    temp_key = f"{temp_k:.4f}"
    key_parts = [pdb, p1, p2, reference, temp_key, method]
    context_uid = stable_uid(*key_parts)
    return {
        "context_uid": context_uid, "complex_id": pdb, "protein1": p1,
        "protein2": p2, "reference": reference, "temperature_k": temp_k,
        "method": method,
        "holdout_type": norm_text(row.get(cols.get("holdout_type") or "", "")),
        "holdout_proteins": norm_text(row.get(cols.get("holdout_proteins") or "", "")),
    }


def schema_sql() -> str:
    return """
    PRAGMA foreign_keys = ON;
    CREATE TABLE IF NOT EXISTS skempi_raw (
        source_row_uid TEXT PRIMARY KEY,
        source_sha256 TEXT NOT NULL,
        source_row_number INTEGER NOT NULL,
        raw_json TEXT NOT NULL,
        parsed_context_uid TEXT,
        parsed_mutations_json TEXT,
        kd_wt_m REAL,
        kd_mut_m REAL,
        import_status TEXT NOT NULL,
        import_reason TEXT NOT NULL DEFAULT ''
    );
    CREATE TABLE IF NOT EXISTS cycle_nodes (
        node_uid TEXT PRIMARY KEY,
        cycle_uid TEXT NOT NULL,
        node_state TEXT NOT NULL CHECK(node_state IN ('P0','Pi','Pj','Pij')),
        source_row_uid TEXT NOT NULL,
        corroborating_source_rows_json TEXT NOT NULL,
        mutation_set_json TEXT NOT NULL,
        complex_id TEXT NOT NULL,
        context_uid TEXT NOT NULL,
        reference TEXT NOT NULL,
        temperature_k REAL NOT NULL,
        measurement_method TEXT NOT NULL,
        kd_m REAL NOT NULL CHECK(kd_m > 0),
        pkd REAL NOT NULL,
        delta_g_kcal_mol REAL NOT NULL,
        FOREIGN KEY(source_row_uid) REFERENCES skempi_raw(source_row_uid)
    );
    CREATE TABLE IF NOT EXISTS structure_runs (
        run_uid TEXT PRIMARY KEY,
        cycle_uid TEXT NOT NULL,
        node_uid TEXT,
        engine TEXT NOT NULL,
        engine_version TEXT NOT NULL,
        command TEXT NOT NULL,
        config_hash TEXT NOT NULL,
        random_seed TEXT NOT NULL,
        input_structure_hash TEXT NOT NULL,
        output_artifact_hash TEXT NOT NULL,
        status TEXT NOT NULL,
        elapsed_seconds REAL,
        cpu_core_seconds REAL,
        gpu_seconds REAL,
        l_rmsd_angstrom REAL,
        details_json TEXT NOT NULL DEFAULT '{}',
        created_utc TEXT NOT NULL,
        FOREIGN KEY(node_uid) REFERENCES cycle_nodes(node_uid),
        FOREIGN KEY(cycle_uid) REFERENCES cycle_features(cycle_uid)
    );
    CREATE TABLE IF NOT EXISTS cycle_features (
        feature_uid TEXT PRIMARY KEY,
        cycle_uid TEXT NOT NULL UNIQUE,
        family_id TEXT,
        complex_id TEXT NOT NULL,
        node_uid_p0 TEXT NOT NULL,
        node_uid_pi TEXT NOT NULL,
        node_uid_pj TEXT NOT NULL,
        node_uid_pij TEXT NOT NULL,
        feature_names_json TEXT NOT NULL,
        feature_values_json TEXT NOT NULL,
        experimental_epistasis_kcal_mol REAL NOT NULL,
        FOREIGN KEY(node_uid_p0) REFERENCES cycle_nodes(node_uid),
        FOREIGN KEY(node_uid_pi) REFERENCES cycle_nodes(node_uid),
        FOREIGN KEY(node_uid_pj) REFERENCES cycle_nodes(node_uid),
        FOREIGN KEY(node_uid_pij) REFERENCES cycle_nodes(node_uid)
    );
    CREATE TABLE IF NOT EXISTS loco_predictions (
        prediction_uid TEXT PRIMARY KEY,
        feature_uid TEXT NOT NULL,
        heldout_group TEXT NOT NULL,
        split_unit TEXT NOT NULL,
        train_group_hash TEXT NOT NULL,
        train_groups_json TEXT NOT NULL,
        fold_index INTEGER NOT NULL,
        y_true REAL NOT NULL,
        y_pred REAL NOT NULL,
        baseline_pred REAL NOT NULL,
        raw_uncertainty REAL NOT NULL,
        calibrated_abs_error REAL NOT NULL,
        calibration_status TEXT NOT NULL,
        FOREIGN KEY(feature_uid) REFERENCES cycle_features(feature_uid)
    );
    CREATE TABLE IF NOT EXISTS validation_manifest (
        manifest_uid TEXT PRIMARY KEY,
        created_utc TEXT NOT NULL,
        input_dataset_sha256 TEXT NOT NULL,
        script_sha256 TEXT NOT NULL,
        split_unit TEXT NOT NULL,
        config_json TEXT NOT NULL,
        software_json TEXT NOT NULL,
        metrics_json TEXT NOT NULL,
        family_count INTEGER NOT NULL,
        cycle_count INTEGER NOT NULL,
        integrity_note TEXT NOT NULL
    );
    CREATE INDEX IF NOT EXISTS idx_nodes_cycle ON cycle_nodes(cycle_uid);
    CREATE INDEX IF NOT EXISTS idx_features_family ON cycle_features(family_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_group ON loco_predictions(heldout_group);
    CREATE INDEX IF NOT EXISTS idx_runs_cycle ON structure_runs(cycle_uid);
    """


def connect_db(path: str | Path) -> sqlite3.Connection:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("PRAGMA foreign_keys = ON")
    con.executescript(schema_sql())
    return con


def import_source_rows(con: sqlite3.Connection, df: pd.DataFrame, cols: dict[str, str | None],
                       source_hash: str) -> tuple[list[dict[str, Any]], Counter]:
    parsed_rows: list[dict[str, Any]] = []
    status_counts: Counter = Counter()
    con.execute("DELETE FROM skempi_raw")
    for ix, series in df.iterrows():
        # Preserve parsed cell strings exactly as pandas read them. Normalization
        # is applied only to derived fields, never to the immutable raw JSON.
        raw = {str(k): ("" if v is None else str(v)) for k, v in series.to_dict().items()}
        source_row_uid = stable_uid(source_hash, ix + 2)
        context = normalized_context(raw, cols)
        mutation_col = cols["mutations"]
        mutations = parse_mutation_set(raw.get(mutation_col, ""))
        wt_parsed_col, mut_parsed_col = cols.get("aff_wt_parsed"), cols.get("aff_mut_parsed")
        kd_wt = parse_exact_positive_number(
            raw.get(cols["aff_wt"], ""), raw.get(wt_parsed_col, "") if wt_parsed_col else "")
        kd_mut = parse_exact_positive_number(
            raw.get(cols["aff_mut"], ""), raw.get(mut_parsed_col, "") if mut_parsed_col else "")
        reasons = []
        if context is None:
            reasons.append("invalid_or_missing_context_or_temperature")
        if mutations is None:
            reasons.append("unparseable_mutation_string")
        if kd_wt is None:
            reasons.append("wild_type_affinity_missing_or_censored")
        if kd_mut is None:
            reasons.append("mutant_affinity_missing_or_censored")
        if context is None or mutations is None:
            status = "not_cycle_parseable"
        elif kd_mut is None:
            status = "mutant_measurement_unusable_wt_may_still_be_usable"
        elif kd_wt is None:
            status = "mutant_usable_wt_unusable"
        else:
            status = "row_measurements_usable"
        context_uid = context["context_uid"] if context else None
        mutation_json = json.dumps(mutations) if mutations is not None else None
        con.execute(
            "INSERT INTO skempi_raw VALUES (?,?,?,?,?,?,?,?,?,?)",
            (source_row_uid, source_hash, int(ix + 2), json.dumps(raw, sort_keys=True),
             context_uid, mutation_json, kd_wt, kd_mut, status, ";".join(reasons))
        )
        status_counts[status] += 1
        if context is not None and mutations is not None:
            parsed_rows.append({
                "source_row_uid": source_row_uid, "source_row_number": int(ix + 2),
                "raw": raw, "context": context, "mutations": mutations,
                "kd_wt": kd_wt, "kd_mut": kd_mut,
            })
    con.commit()
    return parsed_rows, status_counts


def near_identical(values: list[float], rtol: float = 1e-6) -> bool:
    if not values:
        return False
    return all(math.isclose(values[0], x, rel_tol=rtol, abs_tol=0.0) for x in values[1:])


def aa_delta(old: str, new: str) -> np.ndarray:
    return np.asarray(AA_PROP[new], dtype=float) - np.asarray(AA_PROP[old], dtype=float)


def feature_names() -> list[str]:
    names = [f"sub_{a}{b}" for a in AA for b in AA if a != b]
    props = ["hydropathy", "volume", "charge", "aromatic", "polar", "sulfur", "proline", "glycine"]
    for mode in ("sum", "absdiff", "product"):
        names.extend(f"{mode}_{p}" for p in props)
    names.extend([f"location_{x}_count" for x in LOCATION_NAMES])
    names.extend(["location_unknown_count", "same_chain"])
    return names


def substitution_feature_index(old: str, new: str) -> int:
    """Return the index matching feature_names()' custom amino-acid ordering."""
    if old not in AA_INDEX or new not in AA_INDEX:
        raise ValueError(f"Unknown amino acid substitution: {old}->{new}")
    if old == new:
        raise ValueError("A substitution feature requires old != new")

    i = AA_INDEX[old]
    j = AA_INDEX[new]
    return i * 19 + (j if j < i else j - 1)


def cycle_feature_vector(m1: tuple[str, str, str, str], m2: tuple[str, str, str, str],
                         loc1: str, loc2: str) -> np.ndarray:
    names = feature_names()
    x = np.zeros(len(names), dtype=float)
    # Mutation order is canonical/symmetric; count substitutions using
    # the same custom amino-acid order as feature_names().
    for old, new, _chain, _res in sorted([(m1[0], m1[3], m1[1], m1[2]),
                                          (m2[0], m2[3], m2[1], m2[2])]):
        x[substitution_feature_index(old, new)] += 1.0

    idx = 20 * 19
    d1, d2 = aa_delta(m1[0], m1[3]), aa_delta(m2[0], m2[3])
    x[idx:idx+8] = d1 + d2
    x[idx+8:idx+16] = np.abs(d1 - d2)
    x[idx+16:idx+24] = d1 * d2
    idx += 24
    locs = [loc1, loc2]
    for loc in LOCATION_NAMES:
        x[idx] = sum(v == loc for v in locs)
        idx += 1
    x[idx] = sum(v not in LOCATION_NAMES for v in locs)
    idx += 1
    x[idx] = float(m1[1] == m2[1])
    return x


def _best_exact_state(rows: list[dict[str, Any]]) -> tuple[float, list[str]] | None:
    usable = [r for r in rows if r["kd_mut"] is not None]
    if not usable:
        return None
    vals = [float(r["kd_mut"]) for r in usable]
    if not near_identical(vals):
        return None
    # No averaging of disagreeing measurements. Identical repeated source rows
    # are retained as corroboration, with one deterministic source row as FK.
    usable.sort(key=lambda r: r["source_row_uid"])
    return vals[0], [r["source_row_uid"] for r in usable]


def extract_cycles(con: sqlite3.Connection, parsed_rows: list[dict[str, Any]],
                   cols: dict[str, str | None]) -> tuple[pd.DataFrame, dict[str, Any]]:
    by_context: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in parsed_rows:
        by_context[row["context"]["context_uid"]].append(row)
    cycle_rows: list[dict[str, Any]] = []
    rejection_counts: Counter = Counter()
    raw_loc_col = cols.get("location")
    con.execute("DELETE FROM cycle_nodes")
    con.execute("DELETE FROM cycle_features")

    for context_uid, records in sorted(by_context.items()):
        context = records[0]["context"]
        wt_vals = [float(r["kd_wt"]) for r in records if r["kd_wt"] is not None]
        if not wt_vals:
            rejection_counts["no_exact_uncensored_wild_type_baseline"] += 1
            continue
        if not near_identical(wt_vals):
            rejection_counts["conflicting_wild_type_baselines_same_context"] += 1
            continue
        kd0 = wt_vals[0]
        wt_records = sorted([r for r in records if r["kd_wt"] is not None], key=lambda r: r["source_row_uid"])
        p0_source = wt_records[0]["source_row_uid"]
        p0_corroborating = [r["source_row_uid"] for r in wt_records]
        single_rows: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
        double_rows: dict[tuple[tuple[str, str, str, str], ...], list[dict[str, Any]]] = defaultdict(list)
        row_locations: dict[str, dict[tuple[str, str, str, str], str]] = {}
        for r in records:
            muts = r["mutations"]
            raw = r["raw"]
            locs = parse_location_list(raw.get(raw_loc_col, "") if raw_loc_col else "", len(muts))
            row_locations[r["source_row_uid"]] = {m: locs[i] for i, m in enumerate(muts)}
            if r["kd_mut"] is None:
                continue
            if len(muts) == 1:
                single_rows[muts[0]].append(r)
            elif len(muts) == 2:
                double_rows[muts].append(r)
            else:
                rejection_counts["mutation_states_not_single_or_double"] += 1

        single_states: dict[tuple[str, str, str, str], tuple[float, list[str], dict[str, Any]]] = {}
        for mut, entries in single_rows.items():
            best = _best_exact_state(entries)
            if best is None:
                rejection_counts["conflicting_or_missing_single_mutant_measurement"] += 1
                continue
            kd, sources = best
            source = next(x for x in entries if x["source_row_uid"] == sources[0])
            loc = row_locations.get(source["source_row_uid"], {}).get(mut, "unknown")
            single_states[mut] = (kd, sources, {"row": source, "location": loc})

        for m1, m2 in itertools.combinations(sorted(single_states), 2):
            if m1[1:3] == m2[1:3]:
                rejection_counts["two_mutations_same_physical_site"] += 1
                continue
            double_key = tuple(sorted((m1, m2), key=lambda m: (m[1], m[2], m[0], m[3])))
            if double_key not in double_rows:
                rejection_counts["missing_exact_double_mutant_state"] += 1
                continue
            best_double = _best_exact_state(double_rows[double_key])
            if best_double is None:
                rejection_counts["conflicting_or_missing_double_mutant_measurement"] += 1
                continue
            kdij, p_ij_sources = best_double
            record_pij = next(x for x in double_rows[double_key] if x["source_row_uid"] == p_ij_sources[0])
            # Experimental epistasis is computed on the binding free-energy scale.
            # G = RT ln(Kd / 1 M); epsilon = Gij - Gi - Gj + G0.
            t = float(context["temperature_k"])
            g0, gi, gj, gij = [R_KCAL * t * math.log(v) for v in
                               (kd0, single_states[m1][0], single_states[m2][0], kdij)]
            eps = gij - gi - gj + g0
            cycle_uid = stable_uid(context_uid, mutation_token(m1), mutation_token(m2))
            node_specs = [
                ("P0", (), kd0, p0_source, p0_corroborating),
                ("Pi", (m1,), single_states[m1][0], single_states[m1][1][0], single_states[m1][1]),
                ("Pj", (m2,), single_states[m2][0], single_states[m2][1][0], single_states[m2][1]),
                ("Pij", (m1, m2), kdij, p_ij_sources[0], p_ij_sources),
            ]
            node_uids: dict[str, str] = {}
            for state, muts, kd, source_uid, corroborating in node_specs:
                node_uid = stable_uid(cycle_uid, state)
                node_uids[state] = node_uid
                pkd = -math.log10(kd)
                dg = R_KCAL * t * math.log(kd)
                con.execute(
                    """INSERT OR REPLACE INTO cycle_nodes
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (node_uid, cycle_uid, state, source_uid, json.dumps(corroborating),
                     json.dumps([mutation_token(m) for m in muts]), context["complex_id"],
                     context_uid, context["reference"], t, context["method"], kd, pkd, dg)
                )
            loc1 = single_states[m1][2]["location"]
            loc2 = single_states[m2][2]["location"]
            x = cycle_feature_vector(m1, m2, loc1, loc2)
            names = feature_names()
            feature_uid = stable_uid("feature", cycle_uid)
            con.execute(
                """INSERT OR REPLACE INTO cycle_features
                (feature_uid,cycle_uid,family_id,complex_id,node_uid_p0,node_uid_pi,node_uid_pj,
                 node_uid_pij,feature_names_json,feature_values_json,experimental_epistasis_kcal_mol)
                VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                (feature_uid, cycle_uid, None, context["complex_id"], node_uids["P0"],
                 node_uids["Pi"], node_uids["Pj"], node_uids["Pij"], json.dumps(names),
                 json.dumps(x.tolist()), float(eps))
            )
            cycle_rows.append({
                "cycle_uid": cycle_uid, "feature_uid": feature_uid,
                "complex_id": context["complex_id"], "context_uid": context_uid,
                "protein1": context["protein1"], "protein2": context["protein2"],
                "reference": context["reference"], "temperature_k": t,
                "method": context["method"], "mutation_i": mutation_token(m1),
                "mutation_j": mutation_token(m2), "location_i": loc1, "location_j": loc2,
                "kd_p0_m": kd0, "kd_pi_m": single_states[m1][0],
                "kd_pj_m": single_states[m2][0], "kd_pij_m": kdij,
                "epsilon_exp_kcal_mol": eps,
                "node_uid_p0": node_uids["P0"], "node_uid_pi": node_uids["Pi"],
                "node_uid_pj": node_uids["Pj"], "node_uid_pij": node_uids["Pij"],
            })
    con.commit()
    summary = {
        "contexts_seen": len(by_context), "eligible_cycles": len(cycle_rows),
        "cycle_rejection_counts": dict(rejection_counts),
    }
    return pd.DataFrame(cycle_rows), summary


def write_family_template(cycles: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    if cycles.empty:
        # derive distinct complex ids from raw provenance table if no cycles found
        rows = []
    else:
        rows = (cycles[["complex_id", "protein1", "protein2"]]
                .drop_duplicates("complex_id").sort_values("complex_id").to_dict("records"))
    out = pd.DataFrame(rows, columns=["complex_id", "protein1", "protein2"])
    out["family_id"] = ""
    out["family_notes"] = ""
    # Fill these only after checking the exact chain sequences in the cited PDB.
    # They permit a separately reported short-peptide-partner subset.
    out["protein1_length"] = ""
    out["protein2_length"] = ""
    out["peptide_partner"] = ""  # protein1 / protein2 / none
    out.to_csv(path, index=False)
    return path


def audit_command(args: argparse.Namespace) -> None:
    df, delimiter, source_hash = read_skempi(args.input)
    cols = build_column_map(df)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    db_path = outdir / "p4_evidence.sqlite"
    if db_path.exists() and not args.overwrite:
        raise FileExistsError(f"{db_path} exists. Use --overwrite to rebuild this P4 audit database.")
    if db_path.exists():
        db_path.unlink()
    con = connect_db(db_path)
    parsed, row_counts = import_source_rows(con, df, cols, source_hash)
    cycles, cycle_summary = extract_cycles(con, parsed, cols)
    cycles.to_csv(outdir / "eligible_cycles.csv", index=False)
    family_template = write_family_template(cycles, outdir / "family_map_template.csv")
    summary = {
        "created_utc": utc_now(), "input_file": str(Path(args.input).resolve()),
        "input_sha256": source_hash, "delimiter": "TAB" if delimiter == "\t" else delimiter,
        "raw_rows": int(len(df)), "raw_columns": list(df.columns),
        "column_map": cols, "row_status_counts": dict(row_counts),
        **cycle_summary,
        "notes": [
            "Temperature is interpreted as kelvin, matching SKEMPI 2.0 documentation.",
            "Only exact positive numeric Kd values are eligible; censored/inequality values are rejected.",
            "Four states must match exact complex, chain pair, reference, temperature, and method.",
            "Conflicting alternative WT or state measurements are not averaged.",
            "SKEMPI Hold_out_type/Hold_out_proteins are preserved but are not treated as sequence-family labels.",
            "P0 provenance lists all corroborating source rows and links one canonical source row.",
        ],
        "database": str(db_path.resolve()), "family_map_template": str(family_template.resolve()),
    }
    (outdir / "p4_audit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    con.close()
    print(json.dumps(summary, indent=2))


def load_feature_matrix(con: sqlite3.Connection) -> tuple[pd.DataFrame, list[str], np.ndarray, np.ndarray]:
    df = pd.read_sql_query("SELECT * FROM cycle_features ORDER BY cycle_uid", con)
    if df.empty:
        raise ValueError("No eligible four-state cycles are in the database. Inspect p4_audit_summary.json.")
    names = json.loads(df.iloc[0]["feature_names_json"])
    X = np.asarray([json.loads(v) for v in df["feature_values_json"]], dtype=float)
    y = df["experimental_epistasis_kcal_mol"].to_numpy(dtype=float)
    if X.ndim != 2 or X.shape[1] != len(names) or not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("Feature matrix contains invalid dimensions or non-finite values.")
    return df, names, X, y


def attach_family_map(con: sqlite3.Connection, family_map_path: str | Path,
                      split_unit: str) -> pd.DataFrame:
    df = pd.read_csv(family_map_path, dtype=str, keep_default_na=False)
    required = {"complex_id", "family_id"}
    if not required.issubset(df.columns):
        raise ValueError(f"Family map must have columns {sorted(required)}. Found {list(df.columns)}")
    df["complex_id"] = df["complex_id"].map(norm_text)
    df["family_id"] = df["family_id"].map(norm_text)
    if df["complex_id"].duplicated().any():
        raise ValueError("Family map contains duplicate complex_id rows.")
    feat = pd.read_sql_query("SELECT feature_uid, complex_id FROM cycle_features", con)
    if split_unit == "complex":
        mapping = {x: x for x in feat["complex_id"].unique()}
    else:
        mapping = dict(zip(df["complex_id"], df["family_id"]))
        missing = sorted(set(feat["complex_id"]) - set(mapping))
        blank = sorted(k for k in set(feat["complex_id"]) if not mapping.get(k))
        if missing or blank:
            raise ValueError(
                f"Family map incomplete: {len(missing)} missing complexes and {len(blank)} blank family IDs. "
                "Fill family_id for every eligible complex; do not substitute SKEMPI Hold_out_type."
            )
    feat["family_id"] = feat["complex_id"].map(mapping)
    # Optional user-verified chain metadata. Missing values remain missing and
    # cannot silently be treated as peptide labels.
    for col in ("protein1_length", "protein2_length", "peptide_partner"):
        if col in df.columns:
            extra = dict(zip(df["complex_id"], df[col]))
            feat[col] = feat["complex_id"].map(extra)
        else:
            feat[col] = ""
    con.executemany("UPDATE cycle_features SET family_id=? WHERE feature_uid=?",
                    [(r.family_id, r.feature_uid) for r in feat.itertuples()])
    con.commit()
    return feat


def fit_ridge_ensemble(X_train: np.ndarray, y_train: np.ndarray, groups_train: np.ndarray,
                       X_test: np.ndarray, n_boot: int, seed: int,
                       alpha: float = 10.0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Group-bootstrap Ridge ensemble; returns mean, SD, transformed train/test."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups_train)
    if len(uniq) < 2:
        raise ValueError("At least two independent training groups are required for group bootstrap.")
    preds = []
    scaler = StandardScaler().fit(X_train)
    Xtr_s, Xte_s = scaler.transform(X_train), scaler.transform(X_test)
    group_indices = {g: np.flatnonzero(groups_train == g) for g in uniq}
    for _ in range(max(10, int(n_boot))):
        chosen = rng.choice(uniq, size=len(uniq), replace=True)
        indices = np.concatenate([group_indices[g] for g in chosen])
        # Ridge itself is scale-stable; scaler is fit per bootstrap replicate to
        # avoid reusing held-out feature statistics in inner folds.
        model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
        model.fit(X_train[indices], y_train[indices])
        preds.append(model.predict(X_test))
    pred_array = np.asarray(preds, dtype=float)
    mean = pred_array.mean(axis=0)
    sd = pred_array.std(axis=0, ddof=1)
    # Feature-space distance is calculated only against the training partition.
    distances = np.sqrt(((Xte_s[:, None, :] - Xtr_s[None, :, :]) ** 2).mean(axis=2)).min(axis=1)
    raw_u = np.log1p(np.maximum(sd, 0.0)) + np.log1p(np.maximum(distances, 0.0))
    return mean, sd, raw_u, distances


def nested_fold_predictions(X_train: np.ndarray, y_train: np.ndarray, groups_train: np.ndarray,
                            X_test: np.ndarray, n_boot: int, seed: int,
                            alpha: float = 10.0) -> tuple[np.ndarray, np.ndarray, np.ndarray, str, int]:
    """Inner group-OOF predictions fit the uncertainty-to-error calibration only."""
    uniq = np.unique(groups_train)
    n_splits = min(5, len(uniq))
    if n_splits < 3:
        raise ValueError("At least three training groups are required for nested calibration.")
    splitter = GroupKFold(n_splits=n_splits)
    oof_raw_u: list[float] = []
    oof_abs_err: list[float] = []
    oof_pred: list[float] = []
    oof_true: list[float] = []
    for fold, (itr, iva) in enumerate(splitter.split(X_train, y_train, groups_train)):
        if len(np.unique(groups_train[itr])) < 2:
            continue
        pred, _sd, raw_u, _dist = fit_ridge_ensemble(
            X_train[itr], y_train[itr], groups_train[itr], X_train[iva],
            n_boot=n_boot, seed=seed + 1009 * (fold + 1), alpha=alpha)
        oof_raw_u.extend(raw_u.tolist())
        oof_abs_err.extend(np.abs(y_train[iva] - pred).tolist())
        oof_pred.extend(pred.tolist())
        oof_true.extend(y_train[iva].tolist())
    test_pred, _sd, test_raw_u, _dist = fit_ridge_ensemble(
        X_train, y_train, groups_train, X_test, n_boot=n_boot,
        seed=seed + 900001, alpha=alpha)
    calibration_n = len(oof_raw_u)
    if calibration_n >= 12 and len(np.unique(oof_raw_u)) >= 4:
        calibrator = IsotonicRegression(increasing=True, out_of_bounds="clip")
        calibrator.fit(np.asarray(oof_raw_u), np.asarray(oof_abs_err))
        calibrated = calibrator.predict(test_raw_u)
        status = "nested_inner_group_oof_isotonic"
    elif calibration_n >= 1:
        # A constant training-only error estimate is safer than fitting a noisy
        # calibration curve with too few independent OOF observations.
        calibrated = np.full(len(test_pred), float(np.median(oof_abs_err)))
        status = "nested_inner_oof_constant_insufficient_for_ranked_abstention"
    else:
        raise ValueError("No inner out-of-fold predictions available for uncertainty calibration.")
    return test_pred, test_raw_u, calibrated, status, calibration_n


def metric_summary(y: np.ndarray, pred: np.ndarray) -> dict[str, Any]:
    y = np.asarray(y, dtype=float); pred = np.asarray(pred, dtype=float)
    mae = float(np.mean(np.abs(y - pred)))
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    if len(y) >= 3 and np.std(y) > 0 and np.std(pred) > 0:
        rho, p = spearmanr(y, pred)
        rho, p = float(rho), float(p)
    else:
        rho, p = None, None
    informative = np.abs(y) >= 0.5
    direction = (np.sign(y[informative]) == np.sign(pred[informative])) if informative.any() else np.asarray([])
    return {
        "n": int(len(y)), "mae_kcal_mol": mae, "rmse_kcal_mol": rmse,
        "spearman_rho": rho, "spearman_p_value": p,
        "n_informative_abs_epsilon_ge_0_5": int(informative.sum()),
        "direction_accuracy_informative": float(direction.mean()) if len(direction) else None,
    }


def group_bootstrap_ci(df: pd.DataFrame, value_col: str, reps: int = 2000,
                       seed: int = 20261029) -> dict[str, Any]:
    groups = df["heldout_group"].unique().tolist()
    if len(groups) < 5:
        return {"n_groups": len(groups), "ci95": None,
                "status": "not_reported_fewer_than_5_independent_groups"}
    rng = np.random.default_rng(seed)
    # Macro-average independent held-out groups so a large complex family does
    # not dominate the uncertainty interval merely because it has more cycles.
    by_group = {g: float(df.loc[df.heldout_group == g, value_col].mean()) for g in groups}
    observed = float(np.mean(list(by_group.values())))
    vals = []
    for _ in range(reps):
        sample = rng.choice(groups, size=len(groups), replace=True)
        vals.append(float(np.mean([by_group[g] for g in sample])))
    return {"n_groups": len(groups), "macro_mean": observed,
            "ci95": [float(np.quantile(vals, .025)), float(np.quantile(vals, .975))],
            "status": "group_bootstrap_percentile_ci_macro_averaged_conditional_on_fixed_oof_predictions"}


def risk_coverage(preds: pd.DataFrame, reps: int = 500, seed: int = 20261029) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    n = len(preds)
    if n < 10:
        return pd.DataFrame(), {"status": "not_reported_fewer_than_10_oof_cycles", "aurc": None}
    # The abstention decision uses the fold-calibrated expected absolute error,
    # which is in kcal/mol and is the intended cross-fold comparable quantity.
    score = preds["calibrated_abs_error"].to_numpy(float)
    rankable = len(np.unique(np.round(score, 10))) >= max(4, int(.1 * n))
    order = np.argsort(score, kind="mergesort")
    abs_err = np.abs(preds["y_true"].to_numpy(float) - preds["y_pred"].to_numpy(float))
    rows = []
    for coverage_pct in range(10, 101, 10):
        k = max(1, int(math.ceil(n * coverage_pct / 100)))
        selected = order[:k]
        selective_mae = float(abs_err[selected].mean())
        random_maes = []
        for _ in range(reps):
            ix = rng.choice(n, size=k, replace=False)
            random_maes.append(float(abs_err[ix].mean()))
        rows.append({
            "coverage": k / n, "retained_n": k,
            "uncertainty_selected_mae": selective_mae,
            "random_retention_mae_mean": float(np.mean(random_maes)),
            "random_retention_mae_ci95_low": float(np.quantile(random_maes, .025)),
            "random_retention_mae_ci95_high": float(np.quantile(random_maes, .975)),
            "improvement_vs_random_mae": float(np.mean(random_maes) - selective_mae),
            "uncertainty_rankable": bool(rankable),
        })
    tab = pd.DataFrame(rows).sort_values("coverage")
    aurc = float(np.trapezoid(tab["uncertainty_selected_mae"], tab["coverage"])) if rankable else None
    random_aurc = float(np.trapezoid(tab["random_retention_mae_mean"], tab["coverage"]))
    return tab, {
        "status": "rankable" if rankable else "uncertainty_not_sufficiently_variable_to_rank",
        "n_oof": n, "selection_score": "nested_calibrated_expected_absolute_error_kcal_mol",
        "coverage_range": [0.1, 1.0], "aurc_mae_kcal_mol": aurc,
        "random_retention_aurc_mae_kcal_mol": random_aurc,
        "aurc_improvement_vs_random": (random_aurc - aurc) if aurc is not None else None,
    }


def validate_command(args: argparse.Namespace) -> None:
    con = connect_db(args.db)
    family_frame = attach_family_map(con, args.family_map, args.split_unit)
    feat, names, X, y = load_feature_matrix(con)
    group_col = "family_id"
    groups_map = dict(zip(family_frame["complex_id"], family_frame["family_id"]))
    groups = feat["complex_id"].map(groups_map).to_numpy()
    n_groups = len(np.unique(groups))
    if args.peptide_only:
        needed = {"protein1_length", "protein2_length", "peptide_partner"}
        missing_cols = needed - set(family_frame.columns)
        if missing_cols:
            raise ValueError(f"--peptide-only requires filled family-map columns: {sorted(needed)}")
        keep = []
        for row in family_frame.itertuples():
            partner = norm_text(row.peptide_partner).casefold()
            try:
                l1, l2 = int(row.protein1_length), int(row.protein2_length)
            except (TypeError, ValueError):
                keep.append(False); continue
            valid = ((partner == "protein1" and 5 <= l1 <= args.max_peptide_length)
                     or (partner == "protein2" and 5 <= l2 <= args.max_peptide_length))
            keep.append(valid)
        family_frame = family_frame.loc[keep].copy()
        if family_frame.empty:
            raise ValueError("No cycles pass the independently verified peptide-partner length annotation.")
        keep_complexes = set(family_frame["complex_id"])
        feat = feat[feat["complex_id"].isin(keep_complexes)].reset_index(drop=True)
        X = np.asarray([json.loads(v) for v in feat["feature_values_json"]], dtype=float)
        y = feat["experimental_epistasis_kcal_mol"].to_numpy(dtype=float)
        groups_map = dict(zip(family_frame["complex_id"], family_frame["family_id"]))
        groups = feat["complex_id"].map(groups_map).to_numpy()
        n_groups = len(np.unique(groups))
    if n_groups < 5:
        raise ValueError(
            f"Only {n_groups} independent {args.split_unit} groups remain. Need >=5 for the P4 grouped-validation gate."
        )
    if len(y) < 20:
        raise ValueError(f"Only {len(y)} eligible cycles; do not train a predictive model with fewer than 20 cycles.")

    con.execute("DELETE FROM loco_predictions")
    records = []
    outer = GroupKFold(n_splits=n_groups)  # exactly one family/complex per outer fold
    for fold_idx, (train_idx, test_idx) in enumerate(outer.split(X, y, groups), start=1):
        train_groups = groups[train_idx]
        test_groups = np.unique(groups[test_idx])
        if len(test_groups) != 1:
            raise AssertionError("Leave-one-group-out split unexpectedly contains multiple held-out groups.")
        heldout = str(test_groups[0])
        train_hash = hashlib.sha256("\n".join(sorted(map(str, np.unique(train_groups)))).encode()).hexdigest()
        pred, raw_u, calibrated, cal_status, cal_n = nested_fold_predictions(
            X[train_idx], y[train_idx], train_groups, X[test_idx],
            n_boot=args.bootstrap_models, seed=args.seed + fold_idx * 7919, alpha=args.ridge_alpha)
        baseline = float(np.mean(y[train_idx]))
        for local_i, global_i in enumerate(test_idx):
            r = {
                "feature_uid": str(feat.iloc[global_i]["feature_uid"]),
                "heldout_group": heldout, "split_unit": args.split_unit,
                "train_group_hash": train_hash,
                "train_groups_json": json.dumps(sorted(map(str, np.unique(train_groups)))),
                "fold_index": fold_idx,
                "y_true": float(y[global_i]), "y_pred": float(pred[local_i]),
                "baseline_pred": baseline, "raw_uncertainty": float(raw_u[local_i]),
                "calibrated_abs_error": float(max(0.0, calibrated[local_i])),
                "calibration_status": f"{cal_status};inner_oof_n={cal_n}",
            }
            records.append(r)
    pred_df = pd.DataFrame(records)
    if pred_df["feature_uid"].duplicated().any() or set(pred_df["feature_uid"]) != set(feat["feature_uid"]):
        raise AssertionError("Outer predictions must contain every eligible cycle exactly once.")
    for r in pred_df.to_dict("records"):
        pred_uid = stable_uid("prediction", r["feature_uid"], args.split_unit, r["heldout_group"])
        con.execute(
            """INSERT INTO loco_predictions VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (pred_uid, r["feature_uid"], r["heldout_group"], r["split_unit"],
             r["train_group_hash"], r["train_groups_json"], int(r["fold_index"]), float(r["y_true"]),
             float(r["y_pred"]), float(r["baseline_pred"]), float(r["raw_uncertainty"]),
             float(r["calibrated_abs_error"]), r["calibration_status"])
        )
    con.commit()

    p4out = Path(args.out); p4out.mkdir(parents=True, exist_ok=True)
    pred_df.to_csv(p4out / "loco_oof_predictions.csv", index=False)
    fold_manifest = (pred_df.groupby(["fold_index", "heldout_group", "train_group_hash", "train_groups_json"], as_index=False)
                     .agg(n_test_cycles=("feature_uid", "size")))
    fold_manifest["n_train_groups"] = fold_manifest["train_groups_json"].map(lambda x: len(json.loads(x)))
    fold_manifest.to_csv(p4out / "outer_fold_manifest.csv", index=False)
    abs_error = np.abs(pred_df["y_true"] - pred_df["y_pred"])
    pred_df["abs_error"] = abs_error
    pred_df["baseline_abs_error"] = np.abs(pred_df["y_true"] - pred_df["baseline_pred"])
    calibration_bins = pred_df.copy()
    calibration_bins["calibration_bin"] = pd.qcut(
        calibration_bins["calibrated_abs_error"].rank(method="first"),
        q=min(5, max(1, len(calibration_bins) // 5)), labels=False, duplicates="drop")
    calibration_table = (calibration_bins.groupby("calibration_bin", dropna=False)
                         .agg(n=("y_true", "size"),
                              mean_predicted_abs_error=("calibrated_abs_error", "mean"),
                              mean_observed_abs_error=("abs_error", "mean"))
                         .reset_index())
    calibration_table["calibration_bias_pred_minus_observed"] = (
        calibration_table["mean_predicted_abs_error"] - calibration_table["mean_observed_abs_error"])
    calibration_table.to_csv(p4out / "uncertainty_calibration_bins.csv", index=False)
    outer_group_rows = []
    for group_name, block in pred_df.groupby("heldout_group", sort=True):
        group_model_mae = float(np.mean(np.abs(block.y_true - block.y_pred)))
        group_baseline_mae = float(np.mean(np.abs(block.y_true - block.baseline_pred)))
        outer_group_rows.append({"heldout_group": group_name, "n": int(len(block)),
                                 "model_mae": group_model_mae,
                                 "baseline_mae": group_baseline_mae,
                                 "model_minus_baseline_mae": group_model_mae - group_baseline_mae})
    pd.DataFrame(outer_group_rows).to_csv(p4out / "outer_group_metrics.csv", index=False)
    base_metrics = metric_summary(pred_df["y_true"].to_numpy(), pred_df["baseline_pred"].to_numpy())
    model_metrics = metric_summary(pred_df["y_true"].to_numpy(), pred_df["y_pred"].to_numpy())
    error_rho = (spearmanr(pred_df["raw_uncertainty"], abs_error)
                 if len(pred_df) >= 3 and pred_df["raw_uncertainty"].nunique() > 1 and np.std(abs_error) > 0
                 else None)
    calibrated_rho = (spearmanr(pred_df["calibrated_abs_error"], abs_error)
                      if len(pred_df) >= 3 and pred_df["calibrated_abs_error"].nunique() > 1 and np.std(abs_error) > 0
                      else None)
    if error_rho is not None and math.isfinite(float(error_rho.statistic)):
        uncertainty_association = {"raw_uncertainty_abs_error_spearman_rho": float(error_rho.statistic),
                                   "raw_uncertainty_abs_error_p": float(error_rho.pvalue)}
    else:
        uncertainty_association = {"raw_uncertainty_abs_error_spearman_rho": None, "raw_uncertainty_abs_error_p": None}
    if calibrated_rho is not None and math.isfinite(float(calibrated_rho.statistic)):
        uncertainty_association["calibrated_error_abs_error_spearman_rho"] = float(calibrated_rho.statistic)
        uncertainty_association["calibrated_error_abs_error_p"] = float(calibrated_rho.pvalue)
    else:
        uncertainty_association["calibrated_error_abs_error_spearman_rho"] = None
        uncertainty_association["calibrated_error_abs_error_p"] = None
    rc, rc_summary = risk_coverage(pred_df, reps=args.random_repetitions, seed=args.seed)
    if not rc.empty:
        rc.to_csv(p4out / "risk_coverage.csv", index=False)
    pred_df.to_csv(p4out / "loco_oof_predictions_with_errors.csv", index=False)
    # Family-cluster bootstrap CI for paired absolute-error difference.
    paired = pred_df.copy()
    paired["model_minus_baseline_abs_error"] = paired["abs_error"] - paired["baseline_abs_error"]
    model_ci = group_bootstrap_ci(paired.assign(model_abs_error=paired["abs_error"]), "model_abs_error", seed=args.seed)
    baseline_ci = group_bootstrap_ci(paired.assign(base_abs_error=paired["baseline_abs_error"]), "base_abs_error", seed=args.seed + 1)
    diff_ci = group_bootstrap_ci(paired.assign(diff=paired["model_minus_baseline_abs_error"]), "diff", seed=args.seed + 2)

    raw_hash_row = con.execute("SELECT source_sha256 FROM skempi_raw LIMIT 1").fetchone()
    data_hash = raw_hash_row[0] if raw_hash_row else ""
    software = {}
    for pkg in ("numpy", "pandas", "scipy", "scikit-learn"):
        try: software[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError: software[pkg] = "not-found"
    script_hash = sha256_file(__file__)
    metrics = {
        "n_cycles": int(len(pred_df)), "n_independent_groups": int(n_groups),
        "split_unit": args.split_unit, "peptide_only_subset": bool(args.peptide_only),
        "n_independent_groups_after_filters": int(n_groups),
        "model": "group-bootstrap ensemble of fixed-alpha Ridge regressions",
        "ridge_alpha": args.ridge_alpha, "bootstrap_models_per_fit": args.bootstrap_models,
        "baseline_training_fold_mean": base_metrics, "sequence_mutation_descriptor_model": model_metrics,
        "model_minus_baseline_mae": float(model_metrics["mae_kcal_mol"] - base_metrics["mae_kcal_mol"]),
        "group_bootstrap_ci_model_mae": model_ci,
        "group_bootstrap_ci_baseline_mae": baseline_ci,
        "group_bootstrap_ci_model_minus_baseline_mae": diff_ci,
        "uncertainty_error_association": uncertainty_association,
        "uncertainty_calibration": {
            "mean_predicted_abs_error_kcal_mol": float(pred_df["calibrated_abs_error"].mean()),
            "mean_observed_abs_error_kcal_mol": float(abs_error.mean()),
            "mean_calibration_bias_predicted_minus_observed_kcal_mol": float(pred_df["calibrated_abs_error"].mean() - abs_error.mean()),
            "calibration_bins_file": "uncertainty_calibration_bins.csv",
            "note": "All calibration mappings are fit inside each outer training partition using inner group-OOF residuals; calibration tables evaluate only outer held-out predictions.",
        },
        "outer_group_metrics_file": "outer_group_metrics.csv",
        "outer_fold_manifest_file": "outer_fold_manifest.csv",
        "validation_strength": "exploratory_lt_10_independent_groups" if n_groups < 10 else "family_count_at_least_10_still_requires_external_replication",
        "risk_coverage": rc_summary,
        "interpretation": "Experimental epistasis prediction on SKEMPI-like PPI cycles; not peptide-binding evidence for HERV-K Env.",
    }
    manifest_uid = str(uuid.uuid4())
    config = {"seed": args.seed, "bootstrap_models": args.bootstrap_models,
              "random_repetitions": args.random_repetitions, "ridge_alpha": args.ridge_alpha,
              "feature_names": names, "feature_definition": "fixed amino-acid substitution/property descriptors"}
    con.execute(
        """INSERT INTO validation_manifest VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
        (manifest_uid, utc_now(), data_hash, script_hash, args.split_unit,
         json.dumps(config, sort_keys=True), json.dumps(software, sort_keys=True),
         json.dumps(metrics, sort_keys=True), n_groups, len(pred_df),
         "SHA-256 hashes are integrity fingerprints; no digital signature is asserted.")
    )
    con.commit(); con.close()
    (p4out / "p4_validation_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


def require_sha256(value: Any, label: str, allow_blank: bool = False) -> str:
    value = norm_text(value).lower()
    if allow_blank and not value:
        return ""
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{label} must be a 64-character SHA-256 hex digest.")
    return value


def import_structures_command(args: argparse.Namespace) -> None:
    """Import already-completed structural runs; never launches an engine."""
    df = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    required = {"cycle_uid", "engine", "engine_version", "command", "config_hash",
                "random_seed", "input_structure_hash", "output_artifact_hash", "status"}
    if not required.issubset(df.columns):
        raise ValueError(f"Structure manifest missing columns: {sorted(required - set(df.columns))}")
    con = connect_db(args.db)
    known_cycles = {r[0] for r in con.execute("SELECT cycle_uid FROM cycle_features")}
    known_nodes = {r[0]: r[1] for r in con.execute("SELECT node_uid, cycle_uid FROM cycle_nodes")}
    inserted = 0
    for row in df.to_dict("records"):
        cycle_uid = norm_text(row["cycle_uid"])
        node_uid = norm_text(row.get("node_uid", "")) or None
        if cycle_uid not in known_cycles:
            raise ValueError(f"Unknown cycle_uid in structure manifest: {cycle_uid}")
        if node_uid and node_uid not in known_nodes:
            raise ValueError(f"Unknown node_uid in structure manifest: {node_uid}")
        if node_uid and known_nodes[node_uid] != cycle_uid:
            raise ValueError(f"node_uid {node_uid} does not belong to cycle_uid {cycle_uid}")
        def optional_float(key: str) -> float | None:
            v = norm_text(row.get(key, ""))
            if not v: return None
            x = float(v)
            if not math.isfinite(x) or x < 0: raise ValueError(f"Invalid {key}={v}")
            return x
        run_uid = norm_text(row.get("run_uid", "")) or str(uuid.uuid4())
        engine = norm_text(row["engine"])
        engine_version = norm_text(row["engine_version"])
        config_hash = require_sha256(row["config_hash"], "config_hash")
        input_hash = require_sha256(row["input_structure_hash"], "input_structure_hash")
        success_status = norm_text(row["status"]).casefold() in {"success", "completed", "ok"}
        output_hash = require_sha256(row.get("output_artifact_hash", ""), "output_artifact_hash",
                                     allow_blank=not success_status)
        if not all([engine, engine_version, norm_text(row["command"]), norm_text(row["random_seed"]) ]):
            raise ValueError("Structure manifest must include engine, version, command, and random_seed.")
        details = {k: v for k, v in row.items() if k not in required and k not in {
            "run_uid", "node_uid", "elapsed_seconds", "cpu_core_seconds", "gpu_seconds", "l_rmsd_angstrom"}}
        con.execute(
            """INSERT INTO structure_runs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (run_uid, cycle_uid, node_uid, engine, engine_version,
             norm_text(row["command"]), config_hash, norm_text(row["random_seed"]),
             input_hash, output_hash,
             norm_text(row["status"]), optional_float("elapsed_seconds"),
             optional_float("cpu_core_seconds"), optional_float("gpu_seconds"),
             optional_float("l_rmsd_angstrom"), json.dumps(details, sort_keys=True), utc_now())
        )
        inserted += 1
    con.commit(); con.close()
    print(f"Imported {inserted} completed structure-run records. No structural engine was executed.")


def cost_report_command(args: argparse.Namespace) -> None:
    con = connect_db(args.db)
    rows = con.execute("""SELECT status, COUNT(*), SUM(elapsed_seconds), SUM(cpu_core_seconds),
                          SUM(gpu_seconds), SUM(CASE WHEN l_rmsd_angstrom IS NOT NULL THEN 1 ELSE 0 END)
                          FROM structure_runs GROUP BY status ORDER BY status""").fetchall()
    total = con.execute("""SELECT COUNT(*), SUM(elapsed_seconds), SUM(cpu_core_seconds),
                           SUM(gpu_seconds), SUM(CASE WHEN status NOT IN ('success','completed','ok') THEN 1 ELSE 0 END)
                           FROM structure_runs""").fetchone()
    con.close()
    result = {
        "run_count": int(total[0] or 0), "elapsed_wall_hours": (total[1] or 0) / 3600.0,
        "cpu_core_hours": (total[2] or 0) / 3600.0, "gpu_hours": (total[3] or 0) / 3600.0,
        "failed_or_non_success_count": int(total[4] or 0),
        "by_status": [{"status": r[0], "n_runs": r[1],
                       "elapsed_seconds": r[2], "cpu_core_seconds": r[3],
                       "gpu_seconds": r[4], "n_runs_with_l_rmsd": r[5]} for r in rows],
        "accounting_note": "Wall time, CPU core-hours, and GPU-hours are reported separately and are not added together.",
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Silicon Virovore P4 experimental mutation-cycle benchmark")
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("audit", help="ingest raw SKEMPI file, extract strict four-state cycles, and create SQLite evidence DB")
    a.add_argument("--input", required=True, help="Original SKEMPI v2 CSV/TSV file")
    a.add_argument("--out", default="prototype4/results/audit")
    a.add_argument("--overwrite", action="store_true")
    a.set_defaults(func=audit_command)
    v = sub.add_parser("validate", help="nested leave-one-family-out validation and risk-coverage analysis")
    v.add_argument("--db", required=True)
    v.add_argument("--family-map", required=True, help="CSV with complex_id,family_id; family IDs must be independently sequence-clustered")
    v.add_argument("--split-unit", choices=["family", "complex"], default="family")
    v.add_argument("--out", default="prototype4/results/validation")
    v.add_argument("--bootstrap-models", type=int, default=30)
    v.add_argument("--random-repetitions", type=int, default=500)
    v.add_argument("--ridge-alpha", type=float, default=10.0)
    v.add_argument("--peptide-only", action="store_true", help="restrict to complexes with user-verified partner length 5..--max-peptide-length")
    v.add_argument("--max-peptide-length", type=int, default=30)
    v.add_argument("--seed", type=int, default=20261029)
    v.set_defaults(func=validate_command)
    s = sub.add_parser("import-structures", help="import externally completed structural runs; does not launch docking")
    s.add_argument("--db", required=True)
    s.add_argument("--manifest", required=True)
    s.set_defaults(func=import_structures_command)
    c = sub.add_parser("cost-report", help="report wall/CPU/GPU accounting separately")
    c.add_argument("--db", required=True)
    c.add_argument("--out", default="prototype4/results/structure_cost_report.json")
    c.set_defaults(func=cost_report_command)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        print(f"P4 ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
