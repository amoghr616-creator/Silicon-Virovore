from __future__ import annotations

import csv
import itertools
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from prototype4 import p4


class P4PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def synthetic_rows(self, n_families=5):
        columns = [
            "#Pdb", "Mutation(s)_cleaned", "iMutation_Location(s)",
            "Affinity_mut (M)", "Affinity_mut_parsed", "Affinity_wt (M)",
            "Affinity_wt_parsed", "Reference", "Protein 1", "Protein 2",
            "Temperature", "Method", "Hold_out_type", "Hold_out_proteins",
        ]
        rows = []
        aa_sets = [
            ["YA10F", "GB20A", "LC30D", "MD40E"],
            ["KA11E", "VB21G", "SC31N", "TD41K"],
            ["RA12W", "HG22F", "IC32M", "LE42Q"],
            ["FA13Y", "NG23S", "KC33R", "VD43A"],
            ["EA14D", "TG24V", "MC34W", "KL44I"],
        ]
        for fam in range(n_families):
            complex_id = f"{fam+1:04d}_A_B"
            muts = aa_sets[fam % len(aa_sets)]
            kd0 = 1.0e-8 * (1.0 + 0.1 * fam)
            for i, mut in enumerate(muts):
                kd = kd0 * (1.5 + 0.25 * i + 0.1 * fam)
                rows.append([
                    complex_id, mut, ("core" if i % 2 else "rim"),
                    f"{kd:.10g}", f"{kd:.10g}", f"{kd0:.10g}", f"{kd0:.10g}",
                    f"Reference {fam}", f"Receptor {fam}", f"Binder {fam}",
                    "298.15", "SPR", "AB/AG", "",
                ])
            for i, (m1, m2) in enumerate(itertools.combinations(muts, 2)):
                kd = kd0 * (1.1 + 0.13 * i + 0.05 * fam)
                rows.append([
                    complex_id, f"{m1},{m2}", "core,rim",
                    f"{kd:.10g}", f"{kd:.10g}", f"{kd0:.10g}", f"{kd0:.10g}",
                    f"Reference {fam}", f"Receptor {fam}", f"Binder {fam}",
                    "298.15", "SPR", "AB/AG", "",
                ])
        return columns, rows

    def write_synthetic(self, delimiter=";"):
        cols, rows = self.synthetic_rows()
        path = self.root / "skempi_synthetic.csv"
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(cols)
            writer.writerows(rows)
        return path

    def build_audit_db(self, delimiter=";"):
        source = self.write_synthetic(delimiter)
        df, detected, digest = p4.read_skempi(source)
        cols = p4.build_column_map(df)
        db = self.root / "p4.sqlite"
        con = p4.connect_db(db)
        parsed, statuses = p4.import_source_rows(con, df, cols, digest)
        cycles, summary = p4.extract_cycles(con, parsed, cols)
        con.close()
        return source, db, cycles, summary

    def test_delimiter_detection_semicolon_comma_tab(self):
        for delimiter, label in [(";", ";"), (",", ","), ("\t", "\t")]:
            path = self.write_synthetic(delimiter)
            self.assertEqual(p4.detect_delimiter(path), label)

    def test_exact_numeric_rejects_censored_values(self):
        self.assertEqual(p4.parse_exact_positive_number("1.0e-8"), 1.0e-8)
        self.assertIsNone(p4.parse_exact_positive_number(">1.0e-8", "1.0e-8"))
        self.assertIsNone(p4.parse_exact_positive_number("<5e-9"))
        self.assertIsNone(p4.parse_exact_positive_number("0"))

    def test_mutation_parser_and_canonical_order(self):
        a = p4.parse_mutation_set("YA52F,GB17A")
        b = p4.parse_mutation_set("GB17A,YA52F")
        self.assertEqual(a, b)
        self.assertIsNone(p4.parse_mutation_set("YA52F,YA52W"))
        self.assertIsNone(p4.parse_mutation_set("YA52Y"))

    def test_cycle_features_symmetric_under_mutation_order(self):
        m1 = p4.parse_mutation_set("YA52F")[0]
        m2 = p4.parse_mutation_set("GB17A")[0]
        x1 = p4.cycle_feature_vector(m1, m2, "core", "rim")
        x2 = p4.cycle_feature_vector(m2, m1, "rim", "core")
        self.assertTrue(np.array_equal(x1, x2))
        self.assertEqual(len(x1), len(p4.feature_names()))

    def test_extracts_complete_cycles_and_links_four_nodes(self):
        _source, db, cycles, summary = self.build_audit_db()
        self.assertEqual(len(cycles), 30)
        self.assertEqual(summary["eligible_cycles"], 30)
        con = sqlite3.connect(db)
        self.assertEqual(con.execute("SELECT COUNT(*) FROM cycle_nodes").fetchone()[0], 120)
        for row in con.execute("SELECT node_uid_p0,node_uid_pi,node_uid_pj,node_uid_pij FROM cycle_features"):
            self.assertEqual(len(set(row)), 4)
        con.close()

    def test_conflicting_wildtype_values_reject_entire_context(self):
        cols, rows = self.synthetic_rows(n_families=1)
        rows[1][5] = "2e-8"
        rows[1][6] = "2e-8"
        path = self.root / "conflict.csv"
        with path.open("w", newline="") as f:
            w = csv.writer(f, delimiter=";"); w.writerow(cols); w.writerows(rows)
        df, _, digest = p4.read_skempi(path)
        cmap = p4.build_column_map(df)
        con = p4.connect_db(self.root / "conflict.sqlite")
        parsed, _ = p4.import_source_rows(con, df, cmap, digest)
        cycles, summary = p4.extract_cycles(con, parsed, cmap)
        self.assertEqual(len(cycles), 0)
        self.assertEqual(summary["cycle_rejection_counts"]["conflicting_wild_type_baselines_same_context"], 1)
        con.close()

    def test_grouped_nested_validation_writes_one_prediction_per_cycle(self):
        _source, db, cycles, _summary = self.build_audit_db()
        family_map = self.root / "family_map.csv"
        mapping = cycles[["complex_id"]].drop_duplicates().sort_values("complex_id")
        mapping["family_id"] = [f"FAM{i+1}" for i in range(len(mapping))]
        mapping.to_csv(family_map, index=False)
        args = p4.build_parser().parse_args([
            "validate", "--db", str(db), "--family-map", str(family_map),
            "--out", str(self.root / "validation"), "--bootstrap-models", "10",
            "--random-repetitions", "20",
        ])
        p4.validate_command(args)
        result = pd.read_csv(self.root / "validation" / "loco_oof_predictions.csv")
        self.assertEqual(len(result), len(cycles))
        self.assertEqual(result["feature_uid"].nunique(), len(cycles))
        self.assertEqual(result["heldout_group"].nunique(), 5)
        self.assertTrue((result["calibration_status"].str.contains("inner_oof_n=")).all())
        con = sqlite3.connect(db)
        self.assertEqual(con.execute("SELECT COUNT(*) FROM loco_predictions").fetchone()[0], len(cycles))
        self.assertEqual(con.execute("SELECT COUNT(*) FROM validation_manifest").fetchone()[0], 1)
        con.close()


    def test_structure_run_import_requires_realistic_hash_shape_and_cost_fields(self):
        _source, db, cycles, _summary = self.build_audit_db()
        con = sqlite3.connect(db)
        cycle_uid = con.execute("SELECT cycle_uid FROM cycle_features LIMIT 1").fetchone()[0]
        node_uid = con.execute("SELECT node_uid FROM cycle_nodes WHERE cycle_uid=? LIMIT 1", (cycle_uid,)).fetchone()[0]
        con.close()
        manifest = self.root / "structure_runs.csv"
        pd.DataFrame([{
            "cycle_uid": cycle_uid, "node_uid": node_uid, "engine": "TEST_ENGINE",
            "engine_version": "0.0-test", "command": "test invocation",
            "config_hash": "c" * 64, "random_seed": "1",
            "input_structure_hash": "a" * 64, "output_artifact_hash": "b" * 64,
            "status": "success", "elapsed_seconds": "12.5",
            "cpu_core_seconds": "25", "gpu_seconds": "0", "l_rmsd_angstrom": "2.1",
        }]).to_csv(manifest, index=False)
        args = type("Args", (), {"db": str(db), "manifest": str(manifest)})()
        p4.import_structures_command(args)
        con = sqlite3.connect(db)
        row = con.execute("SELECT COUNT(*), SUM(cpu_core_seconds), SUM(gpu_seconds) FROM structure_runs").fetchone()
        con.close()
        self.assertEqual(row[0], 1)
        self.assertEqual(row[1], 25.0)
        self.assertEqual(row[2], 0.0)

    def test_peptide_only_filter_uses_explicit_partner_annotation(self):
        _source, db, cycles, _summary = self.build_audit_db()
        family_map = self.root / "family_map_peptide.csv"
        mapping = cycles[["complex_id"]].drop_duplicates().sort_values("complex_id")
        mapping["family_id"] = [f"FAM{i+1}" for i in range(len(mapping))]
        mapping["protein1_length"] = 300
        mapping["protein2_length"] = 21
        mapping["peptide_partner"] = "protein2"
        mapping.to_csv(family_map, index=False)
        args = p4.build_parser().parse_args([
            "validate", "--db", str(db), "--family-map", str(family_map),
            "--out", str(self.root / "peptide_validation"), "--bootstrap-models", "10",
            "--random-repetitions", "10", "--peptide-only",
        ])
        p4.validate_command(args)
        metrics = json.loads((self.root / "peptide_validation" / "p4_validation_metrics.json").read_text())
        self.assertTrue(metrics["peptide_only_subset"])
        self.assertEqual(metrics["n_cycles"], len(cycles))

    def test_family_map_missing_mapping_is_hard_error(self):
        _source, db, cycles, _summary = self.build_audit_db()
        family_map = self.root / "bad_family_map.csv"
        pd.DataFrame({"complex_id": [cycles.iloc[0].complex_id], "family_id": ["FAM1"]}).to_csv(family_map, index=False)
        con = p4.connect_db(db)
        with self.assertRaisesRegex(ValueError, "Family map incomplete"):
            p4.attach_family_map(con, family_map, "family")
        con.close()


if __name__ == "__main__":
    unittest.main()
