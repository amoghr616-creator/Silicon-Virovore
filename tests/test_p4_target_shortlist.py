from pathlib import Path

import pytest

from src.models import Candidate
from src.p4_candidate_haddock import (
    extract_receptor_chain_a,
    pdb_chain_sequences,
)
from src.p4_target_shortlist import (
    select_p4_structure_shortlist,
    write_p4_target_shortlist,
)


def test_structure_shortlist_is_fitness_seeded_then_max_min_diverse():
    # B is the fitness leader; D is farther from B than the second-highest
    # fitness option C and should receive the second structure slot.
    candidates = [
        Candidate("AAAAAAAAAAAAAAAAAAAAA", c_score=0.95),
        Candidate("AAAAAAAACAAAAAAAAAAAA", c_score=0.90),
        Candidate("RRRRRRRRRRRRRRRRRRRRR", c_score=0.70),
        Candidate("CCCCCCCCCCCCCCCCCCCCC", c_score=0.80),
    ]

    selected = select_p4_structure_shortlist(candidates, limit=2)

    assert selected == [candidates[0], candidates[3]]
    assert candidates[0].metadata["p4_structure_shortlist"]["rank"] == 1
    assert candidates[3].metadata["p4_structure_shortlist"]["rank"] == 2
    assert candidates[1].metadata["p4_structure_shortlist"]["selected"] is False
    assert candidates[2].metadata["p4_structure_shortlist"]["selected"] is False


def test_structure_shortlist_is_deterministic_and_rejects_bad_limit():
    first = [
        Candidate("A" * 21, c_score=0.5),
        Candidate("C" * 21, c_score=0.5),
        Candidate("R" * 21, c_score=0.5),
    ]
    second = [
        Candidate("A" * 21, c_score=0.5),
        Candidate("C" * 21, c_score=0.5),
        Candidate("R" * 21, c_score=0.5),
    ]

    left = [c.sequence for c in select_p4_structure_shortlist(first, 3)]
    right = [c.sequence for c in select_p4_structure_shortlist(second, 3)]
    assert left == right

    with pytest.raises(ValueError):
        select_p4_structure_shortlist(first, 0)


def test_target_shortlist_ranks_only_valid_evidence_and_never_reorders_population(tmp_path):
    better = Candidate("A" * 21, c_score=0.6)
    lower = Candidate("C" * 21, c_score=0.9)
    missing = Candidate("R" * 21, c_score=0.99)

    better.metadata.update({
        "candidate_id": "better-target",
        "p4_target_evidence": {
            "status": "success",
            "target_geometry_score_mean": 0.8,
            "hotspot_proximity_score": 0.7,
        },
        "p4_reliability": {"status": "diagnostic_only"},
    })
    lower.metadata.update({
        "candidate_id": "lower-target",
        "p4_target_evidence": {
            "status": "success",
            "target_geometry_score_mean": 0.7,
            "hotspot_proximity_score": 0.9,
        },
        "p4_reliability": {"status": "diagnostic_only"},
    })
    missing.metadata["p4_target_evidence"] = {"status": "missing_complex"}

    original_order = [better, lower, missing]
    before = [candidate.sequence for candidate in original_order]
    report = write_p4_target_shortlist(
        original_order,
        generation=3,
        output_directory=tmp_path,
    )

    assert [candidate.sequence for candidate in original_order] == before
    assert report["valid_target_evidence_count"] == 2
    assert report["target_evidence_coverage"] == pytest.approx(2 / 3)
    assert report["population_reordered"] is False
    assert report["parent_selection_changed"] is False
    assert report["reliability_calibrated"] is False
    assert report["ranked_candidates"][0]["candidate_id"] == "better-target"
    assert Path(report["report_path"]).is_file()


def test_extract_receptor_chain_a_fallback_is_sequence_verified(tmp_path):
    project_root = Path(__file__).resolve().parents[1]
    source = project_root / "data/receptor/9SYA.pdb"
    if not source.is_file():
        pytest.skip("Bundled 9SYA source PDB is not present in this checkout")

    destination = tmp_path / "9sya_chain_A_only.pdb"
    extracted = extract_receptor_chain_a(source, destination)
    sequences = pdb_chain_sequences(extracted.read_bytes())
    original = pdb_chain_sequences(source.read_bytes())

    assert sequences == {"A": original["A"]}
    assert len(sequences["A"]) == 324
