import tempfile
import unittest
from pathlib import Path

from prototype2.target_audit import (
    TargetMappingError,
    audit_residue_mapping,
    diagnose_pose,
    repair_from_reference,
)


def atom(serial, name, resname, chain, number, x, y, z, icode=""):
    element = name.strip()[0]
    return (
        f"ATOM  {serial:5d} {name:^4} {resname:>3} {chain}{number:4d}{icode:1}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2}\n"
    )


def write_structure(path, residues, missing_atoms=()):
    serial = 1
    lines = []
    for chain, number, icode, resname in residues:
        for name in ("N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"):
            if (chain, number, icode, name) in missing_atoms:
                continue
            lines.append(atom(serial, name, resname, chain, number, float(serial), 0.0, 0.0, icode))
            serial += 1
    path.write_text("".join(lines))


def write_peptide(path, sequence, start_x, chain="B", bad_coordinate=False):
    serial = 1
    lines = []
    names = {"A": "ALA", "G": "GLY", "R": "ARG"}
    for number, one_letter in enumerate(sequence, 1):
        for offset, name in enumerate(("N", "CA", "C", "O", "CB")):
            x = float("nan") if bad_coordinate and serial == 1 else start_x + number + offset * 0.2
            lines.append(atom(serial, name, names[one_letter], chain, number, x, 0.0, 0.0))
            serial += 1
    path.write_text("".join(lines))


class TargetAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.residues = [
            ("A", 201, "", "LEU"),
            ("A", 202, "", "SER"),
            ("A", 203, "", "ARG"),
            ("A", 204, "", "PRO"),
            ("A", 205, "", "ARG"),
            ("A", 206, "", "VAL"),
        ]

    def tearDown(self):
        self.temp.cleanup()

    def test_correct_mapping_and_fail_fast_success(self):
        reference = self.root / "reference.pdb"
        observed = self.root / "observed.pdb"
        write_structure(reference, self.residues)
        write_structure(observed, self.residues)
        result = audit_residue_mapping(observed, reference_path=reference, fail_fast=True)
        self.assertTrue(result["mapping_valid"])
        self.assertEqual(result["target_residues"][0]["pdb_author_residue_number"], 203)

    def test_renumbered_prepared_receptor_maps_by_sequence(self):
        reference = self.root / "reference.pdb"
        observed = self.root / "observed_renumbered.pdb"
        write_structure(reference, self.residues)
        renumbered = [(chain, index, icode, resname) for index, (chain, _, icode, resname) in enumerate(self.residues, 1)]
        write_structure(observed, renumbered)
        result = audit_residue_mapping(observed, reference_path=reference)
        self.assertTrue(result["mapping_valid"])
        self.assertEqual(result["target_residues"][0]["pdb_author_residue_number"], 3)
        self.assertEqual(result["target_residues"][1]["pdb_author_residue_number"], 5)

    def test_missing_anchor_is_reported(self):
        reference = self.root / "reference.pdb"
        observed = self.root / "observed.pdb"
        write_structure(reference, self.residues)
        write_structure(observed, [item for item in self.residues if item[1] != 205])
        result = audit_residue_mapping(observed, reference_path=reference)
        self.assertFalse(result["mapping_valid"])
        self.assertEqual(result["target_residues"][1]["mapping_status"], "not_present_in_observed_sequence")
        with self.assertRaises(TargetMappingError):
            audit_residue_mapping(observed, reference_path=reference, fail_fast=True)

    def test_wrong_residue_identity_is_not_substituted(self):
        reference = self.root / "reference.pdb"
        observed = self.root / "observed.pdb"
        write_structure(reference, self.residues)
        wrong = [item if item[1] != 203 else (item[0], item[1], item[2], "LYS") for item in self.residues]
        write_structure(observed, wrong)
        result = audit_residue_mapping(observed, reference_path=reference)
        self.assertFalse(result["mapping_valid"])
        self.assertEqual(result["target_residues"][0]["mapping_status"], "wrong_residue_identity")

    def test_insertion_code_is_recorded(self):
        reference = self.root / "reference.pdb"
        observed = self.root / "observed.pdb"
        residues = [(c, n, "A" if n == 203 else i, r) for c, n, i, r in self.residues]
        write_structure(reference, residues)
        write_structure(observed, residues)
        result = audit_residue_mapping(observed, reference_path=reference)
        target = result["target_residues"][0]
        self.assertTrue(result["mapping_valid"])
        self.assertEqual(target["pdb_author_insertion_code"], "A")

    def test_missing_guanidinium_atom_is_reported(self):
        reference = self.root / "reference.pdb"
        observed = self.root / "observed.pdb"
        write_structure(reference, self.residues)
        write_structure(observed, self.residues, missing_atoms={("A", 203, "", "NH2")})
        result = audit_residue_mapping(observed, reference_path=reference)
        self.assertFalse(result["mapping_valid"])
        self.assertEqual(result["target_residues"][0]["mapping_status"], "missing_required_side_chain_atoms")
        self.assertEqual(result["target_residues"][0]["missing_arginine_side_chain_atoms"], ["NH2"])

    def test_repair_preserves_residue_and_adds_only_missing_atoms(self):
        reference = self.root / "reference.pdb"
        observed = self.root / "observed.pdb"
        output = self.root / "repaired.pdb"
        write_structure(reference, self.residues)
        write_structure(observed, self.residues, missing_atoms={("A", 203, "", "NH2")})
        result = repair_from_reference(observed, reference, output)
        self.assertEqual(result["preparation_status"], "repaired_from_reference")
        self.assertEqual([item["atom_name"] for item in result["added_atoms"]], ["NH2"])
        self.assertTrue(result["target_mapping_valid"])
        self.assertIn("NH2", output.read_text())

    def test_pose_diagnostic_reports_per_anchor_metrics_and_search_volume(self):
        reference = self.root / "reference.pdb"
        receptor = self.root / "receptor.pdb"
        pose = self.root / "pose.pdb"
        write_structure(reference, self.residues)
        write_structure(receptor, self.residues)
        # The synthetic peptide is near ARG203/ARG205, but not overlapping the
        # receptor at the documented 2 Å severe-clash cutoff.
        write_peptide(pose, "AAAA", 34.0)
        result = diagnose_pose(
            receptor,
            pose,
            "AAAA",
            reference_path=reference,
            site_center=(31.0, 0.0, 0.0),
            box_size_A=20.0,
        )
        self.assertTrue(result["mapping_valid"])
        self.assertTrue(result["sequence_valid"])
        self.assertIn("203", result["target_metrics"])
        self.assertIn("minimum_guanidinium_distance_A", result["target_metrics"]["203"])
        self.assertIn("5.0", result["target_metrics"]["203"]["contact_cutoffs"])
        self.assertIn("peptide_residues", result["target_metrics"]["203"]["contact_cutoffs"]["8.0"])
        self.assertIn("search_volume", result)
        self.assertIn("severe_clash_count", result)

    def test_sequence_mismatch_fails_closed(self):
        reference = self.root / "reference.pdb"
        receptor = self.root / "receptor.pdb"
        pose = self.root / "pose.pdb"
        write_structure(reference, self.residues)
        write_structure(receptor, self.residues)
        write_peptide(pose, "AAAA", 34.0)
        result = diagnose_pose(receptor, pose, "AAGA", reference_path=reference)
        self.assertIn("sequence_mismatch", result["failure_categories"])
        self.assertFalse(result["sequence_valid"])

    def test_ambiguous_peptide_chain_fails_closed(self):
        reference = self.root / "reference.pdb"
        receptor = self.root / "receptor.pdb"
        pose_a = self.root / "pose_a.pdb"
        pose_b = self.root / "pose_b.pdb"
        pose = self.root / "pose.pdb"
        write_structure(reference, self.residues)
        write_structure(receptor, self.residues)
        write_peptide(pose_a, "AAAA", 34.0, chain="B")
        write_peptide(pose_b, "AAAA", 44.0, chain="C")
        pose.write_text(pose_a.read_text() + pose_b.read_text())
        result = diagnose_pose(receptor, pose, "AAAA", reference_path=reference)
        self.assertIn("ambiguous_chain", result["failure_categories"])
        self.assertFalse(result["sequence_valid"])

    def test_malformed_coordinate_record_fails_closed(self):
        reference = self.root / "reference.pdb"
        receptor = self.root / "receptor.pdb"
        pose = self.root / "pose.pdb"
        write_structure(reference, self.residues)
        write_structure(receptor, self.residues)
        write_peptide(pose, "AAAA", 34.0, bad_coordinate=True)
        result = diagnose_pose(receptor, pose, "AAAA", reference_path=reference)
        self.assertIn("malformed_pose", result["failure_categories"])
        self.assertFalse(result["coordinates_valid"])

    def test_severe_clash_is_independent_of_backend_labels(self):
        reference = self.root / "reference.pdb"
        receptor = self.root / "receptor.pdb"
        pose = self.root / "pose.pdb"
        write_structure(reference, self.residues)
        write_structure(receptor, self.residues)
        write_peptide(pose, "AAAA", 23.0)
        result = diagnose_pose(receptor, pose, "AAAA", reference_path=reference)
        self.assertGreater(result["severe_clash_count"], 0)
        self.assertIn("severe_clash", result["failure_categories"])


if __name__ == "__main__":
    unittest.main()
