"""Report-only contact measurements for the HERV-K Env Arg203/Arg205 site."""

from __future__ import annotations

from Bio.PDB import PDBParser


TARGET_RESIDUES = (203, 205)
ARG_GUANIDINIUM_ATOMS = {"NE", "CZ", "NH1", "NH2"}
ARG_NITROGEN_ATOMS = {"NE", "NH1", "NH2"}

ACIDIC_SIDECHAIN_OXYGENS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}


def _heavy_atoms(residue):
    """Return non-hydrogen atoms, including when PDB element fields are blank."""
    atoms = []
    for atom in residue.get_atoms():
        element = (atom.element or atom.get_name()[0]).strip().upper()
        if element not in {"H", "D"}:
            atoms.append(atom)
    return atoms


def measure_arg203_arg205_contacts(
    pose_path: str,
    *,
    receptor_chain_id: str,
    peptide_chain_id: str,
    contact_cutoff_angstrom: float = 4.0,
    salt_bridge_cutoff_angstrom: float = 4.0,
) -> dict:
    """
    Measure contacts in a complex PDB containing receptor and peptide chains.

    Residue numbers are interpreted from the pose's PDB residue IDs.
    Chain IDs must be specified explicitly. This function does not infer
    remapped residue numbering or prove that a pose represents binding.
    """
    if receptor_chain_id == peptide_chain_id:
        raise ValueError("Receptor and peptide chains must be different")
    if contact_cutoff_angstrom <= 0 or salt_bridge_cutoff_angstrom <= 0:
        raise ValueError("Distance cutoffs must be positive")

    structure = PDBParser(QUIET=True).get_structure("complex", pose_path)
    model = structure[0]

    try:
        receptor = model[receptor_chain_id]
        peptide = model[peptide_chain_id]
    except KeyError as exc:
        raise ValueError(f"Required chain is missing: {exc}") from exc

    receptor_residues = {
        residue.id[1]: residue
        for residue in receptor
        if residue.id[0] == " "
    }

    peptide_residues = [
        residue for residue in peptide if residue.id[0] == " "
    ]
    peptide_atom_records = [
        (residue, atom)
        for residue in peptide_residues
        for atom in _heavy_atoms(residue)
    ]

    results = {}

    for residue_number in TARGET_RESIDUES:
        target_residue = receptor_residues.get(residue_number)

        if target_residue is None or target_residue.resname != "ARG":
            actual = (
                None if target_residue is None else target_residue.resname
            )
            raise ValueError(
                f"Expected ARG at receptor chain {receptor_chain_id}, "
                f"residue {residue_number}; found {actual}. "
                "Check the pose's chain IDs and residue numbering."
            )

        target_atoms = [
            atom
            for atom in target_residue
            if atom.get_name() in ARG_GUANIDINIUM_ATOMS
        ]
        target_nitrogens = [
            atom
            for atom in target_residue
            if atom.get_name() in ARG_NITROGEN_ATOMS
        ]

        if not target_atoms or not target_nitrogens:
            raise ValueError(
                f"Arg{residue_number} is missing required side-chain atoms"
            )

        minimum_distance = float("inf")
        contacting_peptide_residues = set()
        salt_bridge_peptide_residues = set()

        for peptide_residue, peptide_atom in peptide_atom_records:
            for target_atom in target_atoms:
                distance = peptide_atom - target_atom
                minimum_distance = min(minimum_distance, distance)

                if distance <= contact_cutoff_angstrom:
                    contacting_peptide_residues.add(
                        peptide_residue.id[1]
                    )

            allowed_acid_atoms = ACIDIC_SIDECHAIN_OXYGENS.get(
                peptide_residue.resname, set()
            )

            if peptide_atom.get_name() not in allowed_acid_atoms:
                continue

            for target_nitrogen in target_nitrogens:
                if (
                    peptide_atom - target_nitrogen
                    <= salt_bridge_cutoff_angstrom
                ):
                    salt_bridge_peptide_residues.add(
                        peptide_residue.id[1]
                    )

        results[f"arg{residue_number}_min_distance_angstrom"] = round(
            minimum_distance, 3
        )
        results[f"arg{residue_number}_contacting_peptide_residue_count"] = (
            len(contacting_peptide_residues)
        )
        results[f"arg{residue_number}_acidic_residue_contact_count"] = (
            len(salt_bridge_peptide_residues)
        )

    results["site_metric_status"] = "measured_not_binding_validation"
    results["receptor_chain_id"] = receptor_chain_id
    results["peptide_chain_id"] = peptide_chain_id
    results["contact_cutoff_angstrom"] = contact_cutoff_angstrom
    results["salt_bridge_cutoff_angstrom"] = salt_bridge_cutoff_angstrom

    return results