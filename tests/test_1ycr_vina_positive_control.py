import unittest
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import numpy as np

from prototype2.run_1ycr_vina_positive_control import (
    BENCHMARK_ROOT,
    _contacts,
    _parse_index_map,
    _parse_pdb_sequences,
)


class YcrVinaPositiveControlTests(unittest.TestCase):
    def test_deposited_and_observed_peptide_residues_map_explicitly(self):
        deposited, observed, mapping, receptor_lines, ligand_lines = _parse_pdb_sequences(
            BENCHMARK_ROOT / "1YCR.pdb"
        )
        self.assertEqual(deposited, "SQETFSDLWKLLPEN")
        self.assertEqual(observed, "ETFSDLWKLLPEN")
        self.assertEqual(len(mapping), 13)
        self.assertEqual(
            [row["deposited_sequence_index_1based"] for row in mapping],
            list(range(3, 16)),
        )
        self.assertEqual(
            [row["author_residue_number"] for row in mapping],
            list(range(17, 30)),
        )
        self.assertEqual(
            "".join(row["observed_residue"] for row in mapping),
            observed,
        )
        self.assertEqual({line[21:22] for line in receptor_lines}, {"A"})
        self.assertEqual({line[21:22] for line in ligand_lines}, {"B"})
        self.assertEqual(len(ligand_lines), 113)

    def test_index_map_parser_handles_multiple_lines_and_duplicates(self):
        text = "\n".join([
            "REMARK INDEX MAP 44 1 45 2",
            "REMARK INDEX MAP 50 3 159 4",
        ])
        self.assertEqual(_parse_index_map(text), {44: 1, 45: 2, 50: 3, 159: 4})
        with self.assertRaises(ValueError):
            _parse_index_map("REMARK INDEX MAP 44 1\nREMARK INDEX MAP 44 2")
        with self.assertRaises(ValueError):
            _parse_index_map("REMARK INDEX MAP 44")
        with self.assertRaises(ValueError):
            _parse_index_map("REMARK INDEX MAP")

    def test_contact_keys_use_peptide_sequence_index(self):
        structure = Structure("synthetic")
        model = Model(0)
        structure.add(model)
        receptor = Chain("A")
        peptide = Chain("B")
        model.add(receptor)
        model.add(peptide)

        receptor_residue = Residue((" ", 25, " "), "ALA", " ")
        receptor_residue.add(Atom(
            "CA", np.array([0.0, 0.0, 0.0]), 1.0, 1.0, " ", " CA ", 1, element="C"
        ))
        receptor.add(receptor_residue)

        peptide_residue = Residue((" ", 17, " "), "GLU", " ")
        peptide_residue.add(Atom(
            "CA", np.array([3.0, 0.0, 0.0]), 1.0, 1.0, " ", " CA ", 2, element="C"
        ))
        peptide.add(peptide_residue)
        self.assertEqual(_contacts(receptor, peptide, 5.0), {(25, 1)})


if __name__ == "__main__":
    unittest.main()