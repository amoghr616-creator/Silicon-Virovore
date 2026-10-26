
from pathlib import Path

import pytest

from prototype3.run_target_aware_pilot import _far_separating_translations
from src.target_aware_scoring import score_peptide_receptor_complex, translate_peptide_chain_pdb


def _atom(serial, atom_name, resname, chain, resid, x, y, z, element="C"):
    return (f"ATOM  {serial:5d} {atom_name:4s} {resname:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}\n")


def _complex(path: Path) -> Path:
    rows = []
    serial = 1
    for resid, xyz in enumerate(((0, 0, 0), (10, 0, 0), (0, 10, 0)), 1):
        rows.append(_atom(serial, "CA", "ALA", "A", resid, *xyz)); serial += 1
    for resid, xyz in enumerate(((3, 0, 0), (3, 1, 0), (3, 0, 1)), 1):
        rows.append(_atom(serial, "CA", "GLY", "B", resid, *xyz)); serial += 1
    rows.append("END\n")
    path.write_text("".join(rows), encoding="utf-8")
    return path


def test_prespecified_far_controls_are_valid_no_contact_poses(tmp_path):
    source = _complex(tmp_path / "complex.pdb")
    translations = _far_separating_translations(source, "A", "B", margin_A=30.0)
    assert len(translations) == 3
    for index, (dx, dy, dz) in enumerate(translations, start=1):
        out = tmp_path / f"control_{index}.pdb"
        translate_peptide_chain_pdb(source, out, "B", dx, dy, dz)
        result = score_peptide_receptor_complex(out, receptor_chain="A", peptide_chain="B")
        assert result["score_valid"] is True
        assert result["structural_score_status"] == "success"
        assert result["target_aware_structural_score"] is not None
        assert result["structural_metrics"]["severe_clash_heavy_atom_pairs"] == 0
        assert result["structural_metrics"]["interchain_heavy_atom_contact_pairs"]["5A"] == 0
        assert result["structural_metrics"]["minimum_interchain_heavy_atom_distance_A"] >= 29.99
        assert any("peptide_far_from_receptor" in x for x in result["eligibility_reasons"])


def test_far_control_margin_must_be_finite_and_not_too_small(tmp_path):
    source = _complex(tmp_path / "complex.pdb")
    with pytest.raises(ValueError):
        _far_separating_translations(source, "A", "B", margin_A=2.0)
    with pytest.raises(ValueError):
        _far_separating_translations(source, "A", "B", margin_A=float("nan"))
