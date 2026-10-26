from __future__ import annotations

from pathlib import Path

import pytest

from src.target_aware_scoring import (
    HotspotSpec,
    StructuralInputError,
    StructuralScoreConfig,
    aligned_peptide_ca_rmsd,
    make_structural_cache_key,
    score_peptide_receptor_complex,
    translate_peptide_chain_pdb,
)


def _pdb_atom(serial: int, atom: str, res: str, chain: str, resid: int,
              x: float, y: float, z: float, element: str = "C") -> str:
    return (f"ATOM  {serial:5d} {atom:4s} {res:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}\n")


def make_complex(path: Path, *, peptide_dx: float = 3.0, receptor_chain: str = "A",
                 peptide_chain: str = "B", coords_override: str | None = None) -> Path:
    lines = []
    serial = 1
    for resid, (x, y, z) in enumerate(((0, 0, 0), (0, 10, 0), (0, 0, 10)), 1):
        lines.append(_pdb_atom(serial, "CA", "ALA", receptor_chain, resid, x, y, z)); serial += 1
    for resid, (x, y, z) in enumerate(((peptide_dx, 0, 0), (peptide_dx, 10, 0), (peptide_dx, 0, 10)), 1):
        lines.append(_pdb_atom(serial, "CA", "GLY", peptide_chain, resid, x, y, z)); serial += 1
    lines.append("END\n")
    content = "".join(lines)
    if coords_override is not None:
        content = content.replace("   3.000", coords_override, 1)
    path.write_text(content, encoding="utf-8")
    return path


def test_valid_complex_has_finite_repeatable_score_and_preserves_c_score_field(tmp_path):
    pdb = make_complex(tmp_path / "valid.pdb")
    a = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="B")
    b = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="B")
    assert a["structural_score_status"] == "success"
    assert a["score_valid"] is True
    assert 0.0 <= a["target_aware_structural_score"] <= 1.0
    assert a["target_aware_structural_score"] == b["target_aware_structural_score"]
    assert a["c_score"] is None
    assert a["c_score_status"] == "not_computed_by_structural_scorer"
    assert a["score_version"] == "target_aware_structural_score_v1.1_clash_gated"
    assert a["structural_metrics"]["buried_surface_area_A2"] is None


def test_missing_file_has_null_score_and_preparation_failure(tmp_path):
    result = score_peptide_receptor_complex(tmp_path / "missing.pdb", receptor_chain="A", peptide_chain="B")
    assert result["target_aware_structural_score"] is None
    assert result["score_valid"] is False
    assert result["structural_score_status"] == "preparation_failure"


def test_missing_peptide_or_wrong_chain_is_invalid(tmp_path):
    pdb = make_complex(tmp_path / "valid.pdb")
    missing = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="C")
    wrong = score_peptide_receptor_complex(pdb, receptor_chain="X", peptide_chain="B")
    assert missing["target_aware_structural_score"] is None
    assert missing["structural_score_status"] == "invalid_input"
    assert wrong["target_aware_structural_score"] is None
    assert wrong["structural_score_status"] == "invalid_input"


def test_same_chain_ids_and_separate_coordinate_files_are_rejected(tmp_path):
    pdb = make_complex(tmp_path / "valid.pdb")
    same_chain = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="A")
    separate = score_peptide_receptor_complex(pdb, tmp_path / "separate.pdb", receptor_chain="A", peptide_chain="B")
    assert same_chain["structural_score_status"] == "invalid_input"
    assert separate["structural_score_status"] == "invalid_input"
    assert separate["target_aware_structural_score"] is None


def test_nonfinite_coordinate_is_rejected(tmp_path):
    pdb = make_complex(tmp_path / "nan.pdb")
    text = pdb.read_text()
    # Set a coordinate field to NaN without changing the fixed-width PDB layout.
    line = text.splitlines()[0]
    line = line[:30] + f"{'nan':>8}" + line[38:]
    pdb.write_text("\n".join([line, *text.splitlines()[1:]]) + "\n")
    result = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="B")
    assert result["target_aware_structural_score"] is None
    assert result["structural_score_status"] == "invalid_input"
    assert "Non-finite" in result["eligibility_reasons"][0]


def test_severe_overlap_is_rejected_with_null_score_and_metrics_preserved(tmp_path):
    pdb = make_complex(tmp_path / "overlap.pdb", peptide_dx=0.0)
    result = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="B")
    assert result["score_valid"] is False
    assert result["structural_score_status"] == "invalid_complex"
    assert result["target_aware_structural_score"] is None
    assert result["structural_metrics"]["severe_clash_heavy_atom_pairs"] >= 3
    assert any("severe_interchain_overlaps" in reason for reason in result["eligibility_reasons"])
    assert result["structural_metrics"]["components"]["severe_clash_penalty_0_1"] > 0
    assert 0.0 <= result["structural_metrics"]["composite_score_diagnostic_before_clash_gate"] <= 1.0


def test_far_no_contact_pose_is_valid_low_score_and_not_zero_failure(tmp_path):
    pdb = make_complex(tmp_path / "far.pdb", peptide_dx=35.0)
    result = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="B")
    assert result["structural_score_status"] == "success"
    assert result["score_valid"] is True
    assert result["target_aware_structural_score"] is not None
    assert result["structural_metrics"]["interchain_heavy_atom_contact_pairs"]["5A"] == 0
    assert result["target_aware_structural_score"] < 0.05
    assert any("peptide_far_from_receptor" in reason for reason in result["eligibility_reasons"])


def test_valid_pose_can_have_zero_contacts_at_stricter_cutoff(tmp_path):
    pdb = make_complex(tmp_path / "near.pdb", peptide_dx=4.4)
    result = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="B")
    contacts = result["structural_metrics"]["interchain_heavy_atom_contact_pairs"]
    assert result["score_valid"] is True
    assert contacts["3.5A"] == 0
    assert contacts["5A"] > 0


def test_hotspot_must_match_author_chain_number_and_residue_name(tmp_path):
    pdb = make_complex(tmp_path / "hotspot.pdb")
    good = score_peptide_receptor_complex(
        pdb, config=StructuralScoreConfig(hotspot_residues=(HotspotSpec("A", 1, "ALA"),)),
        receptor_chain="A", peptide_chain="B")
    bad = score_peptide_receptor_complex(
        pdb, config=StructuralScoreConfig(hotspot_residues=(HotspotSpec("A", 203, "ARG"),)),
        receptor_chain="A", peptide_chain="B")
    assert good["score_valid"] is True
    assert good["structural_metrics"]["hotspot_metrics"][0]["hotspot"] == "A:1:ALA"
    assert bad["score_valid"] is False
    assert bad["target_aware_structural_score"] is None
    assert "hotspot_not_found_or_renumbered" in bad["eligibility_reasons"][0]


def test_controlled_translation_changes_score_and_cache_key(tmp_path):
    native = make_complex(tmp_path / "native.pdb")
    decoy = tmp_path / "decoy.pdb"
    translate_peptide_chain_pdb(native, decoy, "B", 20, 0, 0)
    a = score_peptide_receptor_complex(native, receptor_chain="A", peptide_chain="B")
    b = score_peptide_receptor_complex(decoy, receptor_chain="A", peptide_chain="B")
    assert a["target_aware_structural_score"] > b["target_aware_structural_score"]
    k1 = make_structural_cache_key("GGG", "receptor-hash", "pose-hash")
    k2 = make_structural_cache_key("GGG", "different-receptor-hash", "pose-hash")
    k3 = make_structural_cache_key("GGG", "receptor-hash", "different-pose-hash")
    k4 = make_structural_cache_key("GGG", "receptor-hash", "pose-hash", StructuralScoreConfig(geometry_optimum_A=3.0))
    assert len({k1, k2, k3, k4}) == 4


def test_reference_pose_rmsd_uses_receptor_alignment(tmp_path):
    ref = make_complex(tmp_path / "reference.pdb", peptide_dx=3.0)
    moved = tmp_path / "translated_peptide.pdb"
    translate_peptide_chain_pdb(ref, moved, "B", 0, 0, 2.0)
    same = aligned_peptide_ca_rmsd(ref, ref, receptor_chain="A", peptide_chain="B")
    offset = aligned_peptide_ca_rmsd(ref, moved, receptor_chain="A", peptide_chain="B")
    assert same["receptor_aligned_peptide_CA_RMSD_A"] == pytest.approx(0.0, abs=1e-7)
    assert offset["receptor_aligned_peptide_CA_RMSD_A"] == pytest.approx(2.0, abs=1e-6)
    assert offset["receptor_CA_correspondences"] == 3
    assert offset["peptide_CA_correspondences"] == 3


def test_reference_rmsd_matches_peptide_by_sequence_when_author_numbers_differ(tmp_path):
    ref = make_complex(tmp_path / "reference_numbering.pdb", peptide_dx=3.0)
    moving = tmp_path / "moving_peptide_renumbered.pdb"
    changed = []
    for line in ref.read_text().splitlines():
        if line[:6].strip() in {"ATOM", "HETATM"} and line[21:22] == "B":
            number = int(line[22:26]) + 16
            line = line[:22] + f"{number:4d}" + line[26:]
        changed.append(line)
    moving.write_text("\n".join(changed) + "\n", encoding="utf-8")

    result = aligned_peptide_ca_rmsd(
        ref, moving, receptor_chain="A", peptide_chain="B"
    )
    assert result["peptide_CA_correspondences"] == 3
    assert result["peptide_reference_sequence"] == "GGG"
    assert result["peptide_moving_sequence"] == "GGG"
    assert result["peptide_correspondence_method"] == (
        "exact_sequence_order_ignore_author_residue_number"
    )
    assert result["receptor_aligned_peptide_CA_RMSD_A"] == pytest.approx(
        0.0, abs=1e-7
    )


def test_reference_rmsd_rejects_different_peptide_sequences(tmp_path):
    ref = make_complex(tmp_path / "reference_sequence.pdb", peptide_dx=3.0)
    moving = tmp_path / "different_peptide_sequence.pdb"
    changed = []
    replaced = False
    for line in ref.read_text().splitlines():
        if not replaced and line[:6].strip() == "ATOM" and line[21:22] == "B":
            line = line[:17] + "ALA" + line[20:]
            replaced = True
        changed.append(line)
    moving.write_text("\n".join(changed) + "\n", encoding="utf-8")

    with pytest.raises(
        StructuralInputError, match="peptide_sequence_mismatch_for_RMSD"
    ):
        aligned_peptide_ca_rmsd(
            ref, moving, receptor_chain="A", peptide_chain="B"
        )


def test_malformed_coordinate_field_is_rejected(tmp_path):
    pdb = make_complex(tmp_path / "malformed.pdb")
    lines = pdb.read_text().splitlines()
    lines[0] = lines[0][:30] + " BADCOOR" + lines[0][38:]
    pdb.write_text("\n".join(lines) + "\n")
    result = score_peptide_receptor_complex(pdb, receptor_chain="A", peptide_chain="B")
    assert result["score_valid"] is False
    assert result["target_aware_structural_score"] is None
    assert result["structural_score_status"] == "invalid_input"


def test_reference_rmsd_corrects_global_rotation_and_translation(tmp_path):
    ref = make_complex(tmp_path / "reference.pdb", peptide_dx=3.0)
    moved = tmp_path / "rigid_body_moved.pdb"
    moved_lines = []
    for line in ref.read_text().splitlines():
        if line[:6].strip() in {"ATOM", "HETATM"}:
            x, y, z = (float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
            # Proper rigid transform: 90-degree z rotation + translation.
            xn, yn, zn = -y + 10.0, x - 5.0, z + 3.0
            line = line[:30] + f"{xn:8.3f}{yn:8.3f}{zn:8.3f}" + line[54:]
        moved_lines.append(line)
    moved.write_text("\n".join(moved_lines) + "\n")
    result = aligned_peptide_ca_rmsd(ref, moved, receptor_chain="A", peptide_chain="B")
    assert result["receptor_aligned_peptide_CA_RMSD_A"] == pytest.approx(0.0, abs=1e-6)



def test_reference_rmsd_supports_explicit_moving_chain_mapping(tmp_path):
    ref = make_complex(tmp_path / "reference.pdb", peptide_dx=3.0)
    moved = tmp_path / "different_chain_ids.pdb"
    lines = []
    for line in ref.read_text().splitlines():
        if line[:6].strip() in {"ATOM", "HETATM"}:
            chain = line[21:22]
            if chain == "A":
                line = line[:21] + "X" + line[22:]
            elif chain == "B":
                line = line[:21] + "Y" + line[22:]
        lines.append(line)
    moved.write_text("\n".join(lines) + "\n")
    result = aligned_peptide_ca_rmsd(
        ref, moved, receptor_chain="A", peptide_chain="B",
        moving_receptor_chain="X", moving_peptide_chain="Y")
    assert result["receptor_aligned_peptide_CA_RMSD_A"] == pytest.approx(0.0, abs=1e-7)
    assert result["moving_receptor_chain"] == "X"
    assert result["moving_peptide_chain"] == "Y"
