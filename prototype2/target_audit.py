"""Fail-closed structural audits for the Prototype 2 target workflow.

This module intentionally uses only the Python standard library.  The project
environment does not currently contain BioPython, PDBFixer, OpenMM, or gemmi,
so target identity and coordinate checks must not depend on an optional parser.

The important distinction made here is between a PDB author residue number
and the residue key in a prepared receptor.  Prepared files are mapped to the
author-numbered reference by sequence alignment; a numeric match in a second
file is never treated as proof of identity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA1_TO_3 = {value: key for key, value in AA3_TO_1.items()}
ARG_SIDECHAIN_ATOMS = ("NE", "CZ", "NH1", "NH2")
BACKBONE_ATOMS = ("N", "CA", "C", "O")


class TargetMappingError(ValueError):
    """Raised when the requested target site cannot be established."""


class PoseValidationError(ValueError):
    """Raised when a pose cannot be independently evaluated."""


@dataclass(frozen=True)
class Atom:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float | None
    bfactor: float | None
    element: str
    line_number: int

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)


@dataclass
class Residue:
    chain: str
    resseq: int
    icode: str
    resname: str
    atoms: list[Atom] = field(default_factory=list)

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)

    @property
    def one_letter(self) -> str | None:
        return AA3_TO_1.get(self.resname)

    @property
    def atom_names(self) -> set[str]:
        return {atom.name.strip() for atom in self.atoms}

    @property
    def coordinates_valid(self) -> bool:
        return all(math.isfinite(value) for atom in self.atoms for value in (atom.x, atom.y, atom.z))


def _float_or_none(text: str) -> float | None:
    text = text.strip()
    if not text or text in {".", "?"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_pdb(path: str | Path) -> tuple[list[Atom], list[str]]:
    """Parse coordinate records without normalizing author numbering."""

    path = Path(path)
    atoms: list[Atom] = []
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    for line_number, raw in enumerate(lines, 1):
        line = raw.rstrip("\r\n")
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            continue
        try:
            serial = int(line[6:11].strip() or 0)
            resseq = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = re.sub(r"[^A-Za-z]", "", line[12:16]).upper()[:1]
        atoms.append(
            Atom(
                record=line[0:6].strip(),
                serial=serial,
                name=line[12:16].strip(),
                altloc=line[16:17].strip(),
                resname=line[17:20].strip().upper(),
                chain=line[21:22].strip(),
                resseq=resseq,
                icode=line[26:27].strip(),
                x=x,
                y=y,
                z=z,
                occupancy=_float_or_none(line[54:60]) if len(line) >= 60 else None,
                bfactor=_float_or_none(line[60:66]) if len(line) >= 66 else None,
                element=element,
                line_number=line_number,
            )
        )
    return atoms, lines


def _pdb_record_parse_errors(lines: Sequence[str]) -> list[dict[str, object]]:
    """Report malformed coordinate records that ``parse_pdb`` must skip."""

    errors: list[dict[str, object]] = []
    for line_number, raw in enumerate(lines, 1):
        line = raw.rstrip("\r\n")
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            errors.append({"line": line_number, "reason": "coordinate_record_too_short"})
            continue
        fields = {
            "serial": line[6:11].strip(),
            "author_residue_number": line[22:26].strip(),
            "x": line[30:38].strip(),
            "y": line[38:46].strip(),
            "z": line[46:54].strip(),
        }
        invalid_fields = []
        for name in ("serial", "author_residue_number"):
            try:
                int(fields[name])
            except ValueError:
                invalid_fields.append(name)
        for name in ("x", "y", "z"):
            try:
                float(fields[name])
            except ValueError:
                invalid_fields.append(name)
        if invalid_fields:
            errors.append({"line": line_number, "reason": "invalid_coordinate_record", "fields": invalid_fields})
    return errors


def residues_from_atoms(atoms: Iterable[Atom], chain: str | None = None) -> list[Residue]:
    grouped: dict[tuple[str, int, str], Residue] = {}
    order: list[tuple[str, int, str]] = []
    for atom in atoms:
        if chain is not None and atom.chain != chain:
            continue
        key = atom.residue_key
        if key not in grouped:
            grouped[key] = Residue(atom.chain, atom.resseq, atom.icode, atom.resname)
            order.append(key)
        grouped[key].atoms.append(atom)
    return [grouped[key] for key in order]


def _residue_sort_key(residue: Residue) -> tuple[str, int, str]:
    return (residue.chain, residue.resseq, residue.icode)


def chain_residues(path: str | Path, chain: str) -> list[Residue]:
    atoms, _ = parse_pdb(path)
    return residues_from_atoms(atoms, chain)


def structure_diagnostics(path: str | Path, *, chain: str = "A") -> dict[str, object]:
    """Collect preparation checks without silently repairing anything."""

    path = Path(path)
    atoms, lines = parse_pdb(path)
    parse_errors = _pdb_record_parse_errors(lines)
    residues = residues_from_atoms(atoms, chain)
    duplicate_atom_names = []
    missing_backbone = []
    alternate_locations = []
    invalid_coordinate_atoms = []
    for residue in residues:
        names = [atom.name for atom in residue.atoms]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            duplicate_atom_names.append({"residue": f"{residue.resname}{residue.resseq}{residue.icode}", "atom_names": duplicates})
        missing = [name for name in BACKBONE_ATOMS if name not in residue.atom_names]
        if missing:
            missing_backbone.append({"residue": f"{residue.resname}{residue.resseq}{residue.icode}", "missing_atoms": missing})
        for atom in residue.atoms:
            if atom.altloc:
                alternate_locations.append({"residue": f"{residue.resname}{residue.resseq}{residue.icode}", "atom": atom.name, "altloc": atom.altloc, "occupancy": atom.occupancy})
    # Atom coordinates are scalar fields, so keep the finite-value check
    # explicit and independent of optional structural libraries.
    invalid_coordinate_atoms = [
        {"line": atom.line_number, "atom": atom.name, "residue": f"{atom.resname}{atom.resseq}{atom.icode}"}
        for atom in atoms
        if not all(math.isfinite(value) for value in (atom.x, atom.y, atom.z))
    ]
    atom_keys = {(atom.chain, atom.resseq, atom.icode, atom.name, atom.altloc) for atom in atoms}
    ssbond_records = [line.strip() for line in lines if line.startswith("SSBOND")]
    conect_records = [line.strip() for line in lines if line.startswith("CONECT")]
    residue_numbers = [res.resseq for res in residues]
    numeric_gaps = [
        {"from": left, "to": right, "missing_count": right - left - 1}
        for left, right in zip(residue_numbers, residue_numbers[1:])
        if right - left > 1
    ]
    return {
        "path": str(path),
        "chain": chain,
        "atom_count": len(atoms),
        "residue_count": len(residues),
        "sequence": "".join(res.one_letter or "X" for res in residues),
        "residue_number_range": [min(residue_numbers), max(residue_numbers)] if residue_numbers else None,
        "numeric_residue_gaps": numeric_gaps,
        "missing_backbone_atoms": missing_backbone,
        "alternate_location_records": alternate_locations,
        "insertion_codes": sorted({res.icode for res in residues if res.icode}),
        "hydrogen_atom_count": sum(1 for atom in atoms if _is_hydrogen(atom)),
        "disulfide_records": ssbond_records,
        "conect_records": conect_records,
        "duplicate_atom_keys": len(atoms) - len(atom_keys),
        "invalid_coordinate_atoms": invalid_coordinate_atoms,
        "malformed_coordinate_records": parse_errors,
        "internal_structural_problems": {
            "invalid_coordinates": bool(invalid_coordinate_atoms),
            "malformed_coordinate_records": bool(parse_errors),
            "missing_backbone_atoms": bool(missing_backbone),
            "duplicate_atom_keys": bool(len(atoms) - len(atom_keys)),
        },
    }


def _pdb_metadata(path: Path) -> dict[str, object]:
    metadata: dict[str, object] = {"structure_id": None, "header": None, "title": []}
    for raw in path.read_text(errors="replace").splitlines():
        if raw.startswith("HEADER"):
            metadata["header"] = raw.strip()
            metadata["structure_id"] = raw[62:66].strip() or None if len(raw) >= 66 else None
        elif raw.startswith("TITLE"):
            metadata["title"].append(raw[10:].strip())
    return metadata


def _parse_mmcif_atom_site(path: Path) -> dict[tuple[str, int, str, str], int]:
    """Return (author chain, author seq, insertion, atom) -> label seq.

    This is a deliberately small mmCIF loop reader for atom_site mapping.  It
    supports the ordinary whitespace/quoted values used by deposited files and
    returns an empty mapping for a file without a usable atom_site loop.
    """

    lines = path.read_text(errors="replace").splitlines()
    headers: list[str] = []
    rows: dict[tuple[str, int, str, str], int] = {}
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        j = i + 1
        loop_headers: list[str] = []
        while j < len(lines) and lines[j].strip().startswith("_"):
            loop_headers.append(lines[j].strip().split()[0])
            j += 1
        if not loop_headers or not any(item.startswith("_atom_site.") for item in loop_headers):
            i = j
            continue
        headers = loop_headers
        values: list[str] = []
        while j < len(lines):
            stripped = lines[j].strip()
            if not stripped or stripped.startswith("#"):
                if values:
                    break
                j += 1
                continue
            if stripped == "loop_" or stripped.startswith("_"):
                break
            values.extend(shlex.split(stripped, comments=False))
            j += 1
        width = len(headers)
        indices = {name: headers.index(name) for name in headers}
        required = {
            "_atom_site.auth_asym_id", "_atom_site.auth_seq_id",
            "_atom_site.label_seq_id", "_atom_site.auth_atom_id",
        }
        if not required.issubset(indices):
            return {}
        ins_index = indices.get("_atom_site.pdbx_PDB_ins_code")
        comp_index = indices.get("_atom_site.auth_comp_id")
        for offset in range(0, len(values) - width + 1, width):
            row = values[offset:offset + width]
            try:
                chain = row[indices["_atom_site.auth_asym_id"]]
                seq = int(row[indices["_atom_site.auth_seq_id"]])
                label_seq = int(row[indices["_atom_site.label_seq_id"]])
                atom = row[indices["_atom_site.auth_atom_id"]].strip()
            except (ValueError, IndexError):
                continue
            icode = row[ins_index].strip() if ins_index is not None else ""
            if icode in {".", "?"}:
                icode = ""
            if comp_index is not None and row[comp_index] in {".", "?"}:
                continue
            rows[(chain, seq, icode, atom)] = label_seq
        return rows
    return {}


def _global_alignment(reference: Sequence[str], observed: Sequence[str]) -> tuple[list[tuple[int | None, int | None]], int]:
    """Needleman-Wunsch alignment used only for residue identity mapping."""

    match, mismatch, gap = 2, -1, -2
    n, m = len(reference), len(observed)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [["M"] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * gap
        trace[i][0] = "U"
    for j in range(1, m + 1):
        score[0][j] = j * gap
        trace[0][j] = "L"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            candidates = (
                (score[i - 1][j - 1] + (match if reference[i - 1] == observed[j - 1] else mismatch), "M"),
                (score[i - 1][j] + gap, "U"),
                (score[i][j - 1] + gap, "L"),
            )
            score[i][j], trace[i][j] = max(candidates, key=lambda item: item[0])
    aligned: list[tuple[int | None, int | None]] = []
    i, j = n, m
    while i or j:
        direction = trace[i][j]
        if direction == "M":
            aligned.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif direction == "U":
            aligned.append((i - 1, None))
            i -= 1
        else:
            aligned.append((None, j - 1))
            j -= 1
    aligned.reverse()
    return aligned, score[n][m]


def _context(residues: Sequence[Residue], index: int, radius: int = 2) -> list[dict[str, object]]:
    lo = max(0, index - radius)
    hi = min(len(residues), index + radius + 1)
    return [
        {
            "chain": residues[i].chain,
            "author_residue_number": residues[i].resseq,
            "insertion_code": residues[i].icode,
            "residue_name": residues[i].resname,
            "one_letter": residues[i].one_letter,
        }
        for i in range(lo, hi)
    ]


def _residue_record(residue: Residue, context: list[dict[str, object]], label_seq: int | None = None) -> dict[str, object]:
    return {
        "chain": residue.chain,
        "pdb_author_residue_number": residue.resseq,
        "pdb_author_insertion_code": residue.icode,
        "mmcif_label_residue_number": label_seq,
        "residue_name": residue.resname,
        "one_letter": residue.one_letter,
        "local_sequence_context": context,
        "present": True,
        "coordinates_valid": residue.coordinates_valid,
        "atom_names": sorted(residue.atom_names),
        "required_arginine_side_chain_atoms_present": all(atom in residue.atom_names for atom in ARG_SIDECHAIN_ATOMS),
        "missing_arginine_side_chain_atoms": [atom for atom in ARG_SIDECHAIN_ATOMS if atom not in residue.atom_names],
        "is_arginine": residue.resname == "ARG",
    }


def audit_residue_mapping(
    structure_path: str | Path,
    *,
    chain: str = "A",
    target_author_residues: Sequence[int] = (203, 205),
    reference_path: str | Path | None = None,
    mmcif_path: str | Path | None = None,
    fail_fast: bool = False,
) -> dict[str, object]:
    """Audit target residue identity in one PDB representation.

    If ``reference_path`` is provided, target residues are identified in the
    reference sequence and mapped into the observed file by global sequence
    alignment.  This handles renumbered/cropped prepared receptors.
    """

    structure_path = Path(structure_path)
    atoms, _ = parse_pdb(structure_path)
    observed = [residue for residue in residues_from_atoms(atoms, chain) if residue.one_letter]
    reference = (
        [residue for residue in chain_residues(reference_path, chain) if residue.one_letter]
        if reference_path
        else observed
    )
    label_map = _parse_mmcif_atom_site(Path(mmcif_path)) if mmcif_path else {}
    alignment, alignment_score = _global_alignment(
        [res.one_letter or "X" for res in reference],
        [res.one_letter or "X" for res in observed],
    )
    ref_to_observed: dict[int, int] = {
        ref_i: obs_i
        for ref_i, obs_i in alignment
        if ref_i is not None and obs_i is not None
    }

    records: list[dict[str, object]] = []
    valid = True
    failure_reasons: list[str] = []
    for requested in target_author_residues:
        ref_indices = [i for i, residue in enumerate(reference) if residue.resseq == requested]
        if len(ref_indices) != 1:
            valid = False
            reason = f"reference chain {chain} does not contain exactly one author residue {requested}"
            failure_reasons.append(reason)
            records.append({
                "requested_author_residue_number": requested,
                "present": False,
                "mapping_status": "reference_missing_or_ambiguous",
                "diagnostic": reason,
            })
            continue
        ref_index = ref_indices[0]
        ref_residue = reference[ref_index]
        obs_index = ref_to_observed.get(ref_index)
        if obs_index is None:
            valid = False
            reason = f"author residue {requested} cannot be aligned into {structure_path.name}"
            failure_reasons.append(reason)
            records.append({
                "requested_author_residue_number": requested,
                "reference_residue_name": ref_residue.resname,
                "present": False,
                "mapping_status": "not_present_in_observed_sequence",
                "diagnostic": reason,
                "reference_local_sequence_context": _context(reference, ref_index),
            })
            continue
        observed_residue = observed[obs_index]
        label_seq = None
        for atom in observed_residue.atoms:
            label_seq = label_map.get((atom.chain, atom.resseq, atom.icode, atom.name))
            if label_seq is not None:
                break
        record = _residue_record(observed_residue, _context(observed, obs_index), label_seq)
        record.update({
            "requested_author_residue_number": requested,
            "reference_author_residue_number": ref_residue.resseq,
            "reference_residue_name": ref_residue.resname,
            "reference_sequence_index": ref_index,
            "observed_sequence_index": obs_index,
            "sequence_alignment_identity": ref_residue.one_letter == observed_residue.one_letter,
            "mapping_status": "mapped",
        })
        if ref_residue.resname != "ARG" or observed_residue.resname != "ARG":
            valid = False
            reason = (
                f"requested author residue {requested} mapped as "
                f"reference={ref_residue.resname}, observed={observed_residue.resname}"
            )
            failure_reasons.append(reason)
            record["mapping_status"] = "wrong_residue_identity"
            record["diagnostic"] = reason
        elif not record["required_arginine_side_chain_atoms_present"]:
            valid = False
            reason = (
                f"ARG{requested} is missing required guanidinium atoms: "
                f"{', '.join(record['missing_arginine_side_chain_atoms'])}"
            )
            failure_reasons.append(reason)
            record["mapping_status"] = "missing_required_side_chain_atoms"
            record["diagnostic"] = reason
        records.append(record)

    result: dict[str, object] = {
        "structure_identifier": structure_path.stem,
        "structure_path": str(structure_path),
        "structure_sha256": hashlib.sha256(structure_path.read_bytes()).hexdigest(),
        "source_metadata": _pdb_metadata(structure_path),
        "chain": chain,
        "reference_path": str(reference_path) if reference_path else None,
        "mmcif_path": str(mmcif_path) if mmcif_path else None,
        "author_numbering_is_observed_numbering": reference_path is None,
        "sequence_alignment_score": alignment_score,
        "observed_residue_count": len(observed),
        "target_residues": records,
        "mapping_valid": valid,
        "failure_reasons": failure_reasons,
    }
    if fail_fast and not valid:
        raise TargetMappingError(json.dumps(result, indent=2))
    return result


def compare_receptor_mappings(
    paths: Sequence[str | Path],
    *,
    reference_path: str | Path,
    chain: str = "A",
    mmcif_path: str | Path | None = None,
    target_author_residues: Sequence[int] = (203, 205),
) -> dict[str, object]:
    audits = [
        audit_residue_mapping(
            path,
            chain=chain,
            target_author_residues=target_author_residues,
            reference_path=reference_path,
            mmcif_path=mmcif_path,
        )
        for path in paths
    ]
    return {
        "reference_path": str(reference_path),
        "chain": chain,
        "target_author_residues": list(target_author_residues),
        "all_representations_valid": all(item["mapping_valid"] for item in audits),
        "representations": audits,
    }


def _format_atom_line(serial: int, source: Atom, *, altloc: str = "") -> str:
    element = source.element.rjust(2)
    occupancy = source.occupancy if source.occupancy is not None else 1.0
    bfactor = source.bfactor if source.bfactor is not None else 0.0
    return (
        f"ATOM  {serial:5d} {source.name:^4}{altloc:1}{source.resname:>3} "
        f"{source.chain:1}{source.resseq:4d}{source.icode:1}   "
        f"{source.x:8.3f}{source.y:8.3f}{source.z:8.3f}"
        f"{occupancy:6.2f}{bfactor:6.2f}          {element}\n"
    )


def _selected_reference_atoms(residue: Residue) -> dict[str, Atom]:
    by_altloc: dict[str, list[Atom]] = {}
    for atom in residue.atoms:
        by_altloc.setdefault(atom.altloc, []).append(atom)
    if "" in by_altloc:
        preferred = ""
    else:
        preferred = max(
            by_altloc,
            key=lambda alt: sum((atom.occupancy or 0.0) for atom in by_altloc[alt]),
        )
    selected = {atom.name: atom for atom in by_altloc.get("", [])}
    selected.update({atom.name: atom for atom in by_altloc.get(preferred, [])})
    return selected


def repair_from_reference(
    input_path: str | Path,
    reference_path: str | Path,
    output_path: str | Path,
    *,
    chain: str = "A",
) -> dict[str, object]:
    """Repair missing heavy atoms from a matching deposited conformer.

    This is intentionally conservative: it never deletes residues, changes
    existing coordinates, adds hydrogens, or renumbers atoms.  Only atom names
    present in the same residue/conformer of the reference but absent from the
    input are inserted.  The output is a new file.
    """

    input_path, reference_path, output_path = map(Path, (input_path, reference_path, output_path))
    input_atoms, input_lines = parse_pdb(input_path)
    reference_atoms, _ = parse_pdb(reference_path)
    input_residues = {res.key: res for res in residues_from_atoms(input_atoms, chain)}
    reference_residues = {res.key: res for res in residues_from_atoms(reference_atoms, chain)}
    input_names = {key: res.atom_names for key, res in input_residues.items()}
    missing: list[dict[str, object]] = []
    insertions: dict[int, list[str]] = {}
    next_serial = max((atom.serial for atom in input_atoms), default=0) + 1

    # Locate the last input line for each residue.  The source files are PDB
    # records in residue order; preserving that order keeps downstream tools
    # deterministic.
    last_line_for_key: dict[tuple[str, int, str], int] = {}
    for atom in input_atoms:
        last_line_for_key[atom.residue_key] = max(last_line_for_key.get(atom.residue_key, 0), atom.line_number)

    for key, input_residue in input_residues.items():
        reference_residue = reference_residues.get(key)
        if reference_residue is None or input_residue.resname != reference_residue.resname:
            continue
        selected = _selected_reference_atoms(reference_residue)
        for name, source_atom in selected.items():
            if source_atom.element == "H" or name in input_names[key]:
                continue
            serial = next_serial
            next_serial += 1
            line = _format_atom_line(serial, source_atom)
            insertion_line = last_line_for_key[key]
            insertions.setdefault(insertion_line, []).append(line)
            missing.append({
                "chain": key[0],
                "author_residue_number": key[1],
                "insertion_code": key[2],
                "residue_name": input_residue.resname,
                "atom_name": name,
                "reference_altloc": source_atom.altloc,
                "reference_serial": source_atom.serial,
                "coordinates": [source_atom.x, source_atom.y, source_atom.z],
            })

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for line_number, line in enumerate(input_lines, 1):
            handle.write(line)
            for addition in insertions.get(line_number, []):
                handle.write(addition)

    before = audit_residue_mapping(input_path, chain=chain, reference_path=reference_path)
    after = audit_residue_mapping(output_path, chain=chain, reference_path=reference_path)
    before_diagnostics = structure_diagnostics(input_path, chain=chain)
    after_diagnostics = structure_diagnostics(output_path, chain=chain)
    input_coordinate_map = {
        (atom.chain, atom.resseq, atom.icode, atom.name, atom.altloc): (atom.x, atom.y, atom.z)
        for atom in input_atoms
    }
    output_atoms, _ = parse_pdb(output_path)
    changed_existing_coordinates = []
    for atom in output_atoms:
        key = (atom.chain, atom.resseq, atom.icode, atom.name, atom.altloc)
        if key in input_coordinate_map and input_coordinate_map[key] != (atom.x, atom.y, atom.z):
            changed_existing_coordinates.append({"key": key, "before": input_coordinate_map[key], "after": (atom.x, atom.y, atom.z)})
    report = {
        "repair_method": "copy_missing_heavy_atoms_from_matching_reference_conformer",
        "input_path": str(input_path),
        "reference_path": str(reference_path),
        "output_path": str(output_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "reference_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        "output_sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "settings": {
            "chain": chain,
            "selected_altloc": "highest occupancy when no blank conformer exists",
            "hydrogens_added": False,
            "existing_coordinates_changed": bool(changed_existing_coordinates),
            "residues_deleted": False,
            "renumbering_performed": False,
        },
        "added_atoms": missing,
        "input_structure_diagnostics": before_diagnostics,
        "output_structure_diagnostics": after_diagnostics,
        "missing_reference_atoms_repaired": len(missing),
        "changed_existing_coordinates": changed_existing_coordinates,
        "pre_repair_mapping": before,
        "post_repair_mapping": after,
        "preparation_status": "repaired_from_reference" if missing else "no_repairs_needed",
        "target_mapping_valid": bool(after["mapping_valid"]),
    }
    return report


def _distance(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _is_hydrogen(atom: Atom) -> bool:
    return atom.element.upper() == "H" or atom.name.upper().startswith("H")


def _find_peptide_residues(atoms: Sequence[Atom], sequence: str) -> list[Residue]:
    matches: list[list[Residue]] = []
    for chain in sorted({atom.chain for atom in atoms}):
        residues = [res for res in residues_from_atoms(atoms, chain) if res.one_letter]
        if "".join(res.one_letter for res in residues) == sequence:
            matches.append(residues)
    if len(matches) != 1:
        raise PoseValidationError(
            f"expected exactly one chain with peptide sequence {sequence}, found {len(matches)}"
        )
    return matches[0]


def _mapped_target_atoms(
    receptor_atoms: Sequence[Atom],
    mapping: dict[str, object],
    target_author_residues: Sequence[int],
) -> dict[int, list[Atom]]:
    """Return target-residue atoms keyed by the requested author number."""

    records = {
        int(record["requested_author_residue_number"]): record
        for record in mapping["target_residues"]
    }
    result: dict[int, list[Atom]] = {}
    for requested in target_author_residues:
        record = records[requested]
        observed_key = (
            record["chain"],
            record["pdb_author_residue_number"],
            record["pdb_author_insertion_code"],
        )
        result[requested] = [
            atom for atom in receptor_atoms
            if atom.residue_key == observed_key and not _is_hydrogen(atom)
        ]
    return result


def _target_atom_coordinates(
    receptor_path: str | Path,
    *,
    reference_path: str | Path,
    target_author_residues: Sequence[int],
    chain: str,
) -> dict[tuple[int, str], tuple[float, float, float]]:
    """Get target guanidinium coordinates in a mapped representation."""

    atoms, _ = parse_pdb(receptor_path)
    mapping = audit_residue_mapping(
        receptor_path,
        chain=chain,
        target_author_residues=target_author_residues,
        reference_path=reference_path,
    )
    if not mapping["mapping_valid"]:
        raise TargetMappingError(str(mapping["failure_reasons"]))
    target_atoms = _mapped_target_atoms(atoms, mapping, target_author_residues)
    return {
        (requested, atom.name): (atom.x, atom.y, atom.z)
        for requested, residue_atoms in target_atoms.items()
        for atom in residue_atoms
        if atom.name in ARG_SIDECHAIN_ATOMS
    }


def coordinate_frame_audit(
    paths: Sequence[str | Path],
    *,
    reference_path: str | Path,
    chain: str = "A",
    target_author_residues: Sequence[int] = (203, 205),
    site_center: Sequence[float] | None = None,
    tolerance_A: float = 0.5,
) -> dict[str, object]:
    """Compare target coordinates across receptor representations.

    Sequence-aligned residue mapping is performed independently for every
    representation.  The reference coordinates are then compared by target
    author residue and guanidinium atom, so a prepared-file residue number is
    never confused with the biological numbering.
    """

    reference = Path(reference_path)
    canonical = _target_atom_coordinates(
        reference,
        reference_path=reference,
        target_author_residues=target_author_residues,
        chain=chain,
    )
    rows: list[dict[str, object]] = []
    for path_value in paths:
        path = Path(path_value)
        mapping = audit_residue_mapping(
            path,
            chain=chain,
            target_author_residues=target_author_residues,
            reference_path=reference,
        )
        row: dict[str, object] = {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "mapping_valid": mapping["mapping_valid"],
            "target_residues": mapping["target_residues"],
            "coordinate_transform_status": "unresolved",
            "coordinate_frame_comparable": False,
            "max_target_atom_delta_A": None,
            "rms_target_atom_delta_A": None,
            "target_site_center": list(site_center) if site_center is not None else None,
            "target_site_center_distances_A": {},
        }
        if not mapping["mapping_valid"]:
            row["failure_reasons"] = mapping["failure_reasons"]
            rows.append(row)
            continue
        observed = _target_atom_coordinates(
            path,
            reference_path=reference,
            target_author_residues=target_author_residues,
            chain=chain,
        )
        deltas = [
            math.sqrt(sum((observed[key][axis] - canonical[key][axis]) ** 2 for axis in range(3)))
            for key in canonical
            if key in observed
        ]
        row["max_target_atom_delta_A"] = round(max(deltas), 6) if deltas else None
        row["rms_target_atom_delta_A"] = round(math.sqrt(sum(delta * delta for delta in deltas) / len(deltas)), 6) if deltas else None
        comparable = bool(deltas) and max(deltas) <= tolerance_A
        row["coordinate_frame_comparable"] = comparable
        row["coordinate_transform_status"] = "no_coordinate_transform_detected" if comparable else "coordinate_mismatch"
        if site_center is not None:
            center = tuple(float(value) for value in site_center)
            for requested in target_author_residues:
                coords = [observed[(requested, atom)] for atom in ARG_SIDECHAIN_ATOMS if (requested, atom) in observed]
                if coords:
                    centroid = tuple(sum(point[axis] for point in coords) / len(coords) for axis in range(3))
                    row["target_site_center_distances_A"][str(requested)] = round(
                        math.sqrt(sum((centroid[axis] - center[axis]) ** 2 for axis in range(3))), 6
                    )
        rows.append(row)
    return {
        "reference_path": str(reference),
        "chain": chain,
        "target_author_residues": list(target_author_residues),
        "coordinate_comparison_tolerance_A": tolerance_A,
        "all_mappings_valid": all(row["mapping_valid"] for row in rows),
        "all_coordinate_frames_comparable": all(row["coordinate_frame_comparable"] for row in rows),
        "representations": rows,
    }


def _centroid(atoms: Sequence[Atom]) -> tuple[float, float, float] | None:
    if not atoms:
        return None
    return tuple(
        sum(getattr(atom, axis) for atom in atoms) / len(atoms)
        for axis in ("x", "y", "z")
    )


def _distance_to_point(point: Sequence[float], center: Sequence[float]) -> float:
    return math.sqrt(sum((float(point[i]) - float(center[i])) ** 2 for i in range(3)))


def _pairwise_min_and_clashes(
    receptor_atoms: Sequence[Atom],
    peptide_atoms: Sequence[Atom],
    severe_clash_cutoff: float,
) -> tuple[float, int]:
    """Return exact minimum distance and severe-overlap count.

    NumPy is used when available because the batch diagnostic can contain
    hundreds of poses.  The standard-library fallback keeps single-pose
    audits usable in minimal environments.
    """

    try:
        import numpy as np  # type: ignore

        receptor_coords = np.asarray([(atom.x, atom.y, atom.z) for atom in receptor_atoms], dtype=float)
        peptide_coords = np.asarray([(atom.x, atom.y, atom.z) for atom in peptide_atoms], dtype=float)
        distances_squared = ((receptor_coords[:, None, :] - peptide_coords[None, :, :]) ** 2).sum(axis=2)
        return (
            float(math.sqrt(float(distances_squared.min()))),
            int((distances_squared < severe_clash_cutoff ** 2).sum()),
        )
    except ImportError:
        minimum = math.inf
        clashes = 0
        cutoff_squared = severe_clash_cutoff ** 2
        for receptor_atom in receptor_atoms:
            for peptide_atom in peptide_atoms:
                distance_squared = sum(
                    (getattr(receptor_atom, axis) - getattr(peptide_atom, axis)) ** 2
                    for axis in ("x", "y", "z")
                )
                minimum = min(minimum, distance_squared)
                if distance_squared < cutoff_squared:
                    clashes += 1
        return math.sqrt(minimum), clashes


def diagnose_pose(
    receptor_path: str | Path,
    pose_path: str | Path,
    peptide_sequence: str,
    *,
    target_author_residues: Sequence[int] = (203, 205),
    chain: str = "A",
    reference_path: str | Path | None = None,
    site_center: Sequence[float] | None = None,
    box_size_A: float | None = None,
    contact_cutoffs_A: Sequence[float] = (3.0, 4.0, 5.0, 8.0),
    structural_contact_cutoff_A: float = 4.5,
    severe_clash_cutoff: float = 2.0,
    mapping_result: dict[str, object] | None = None,
    coordinate_frame_result: dict[str, object] | None = None,
) -> dict[str, object]:
    """Produce raw, per-anchor diagnostics for one parseable pose.

    This function is deliberately descriptive.  It does not tune a contact
    threshold or convert a docking score into a binding claim.
    """

    receptor_path, pose_path = Path(receptor_path), Path(pose_path)
    receptor_atoms, receptor_lines = parse_pdb(receptor_path)
    pose_atoms, pose_lines = parse_pdb(pose_path)
    parse_errors = {
        "receptor": _pdb_record_parse_errors(receptor_lines),
        "pose": _pdb_record_parse_errors(pose_lines),
    }
    result: dict[str, object] = {
        "receptor_path": str(receptor_path),
        "pose_path": str(pose_path),
        "peptide_sequence_requested": peptide_sequence,
        "site_center_A": list(site_center) if site_center is not None else None,
        "box_size_A": box_size_A,
        "severe_clash_cutoff_A": severe_clash_cutoff,
        "contact_cutoffs_A": list(contact_cutoffs_A),
        "parse_errors": parse_errors,
        "failure_categories": [],
    }
    categories: list[str] = result["failure_categories"]
    if any(parse_errors.values()):
        categories.append("malformed_pose")

    try:
        mapping = mapping_result or audit_residue_mapping(
            receptor_path,
            chain=chain,
            target_author_residues=target_author_residues,
            reference_path=reference_path,
        )
    except (OSError, ValueError) as exc:
        result.update({"mapping_valid": False, "diagnostic_error": str(exc)})
        categories.append("wrong_mapping")
        result["primary_failure_category"] = categories[0]
        return result
    result["mapping_valid"] = mapping["mapping_valid"]
    result["target_mapping"] = mapping["target_residues"]
    if not mapping["mapping_valid"]:
        categories.append("wrong_mapping")
        result["mapping_failure_reasons"] = mapping["failure_reasons"]
        result["primary_failure_category"] = categories[0]
        return result

    try:
        peptide_residues = _find_peptide_residues(pose_atoms, peptide_sequence)
    except PoseValidationError as exc:
        if "found 0" in str(exc):
            categories.append("sequence_mismatch")
        elif "found" in str(exc):
            categories.append("ambiguous_chain")
        else:
            categories.append("malformed_pose")
        result.update({"sequence_valid": False, "diagnostic_error": str(exc)})
        result["primary_failure_category"] = categories[0]
        return result
    peptide_atoms = [atom for residue in peptide_residues for atom in residue.atoms if not _is_hydrogen(atom)]
    receptor_heavy = [atom for atom in receptor_atoms if not _is_hydrogen(atom)]
    if not peptide_atoms or not receptor_heavy:
        categories.append("malformed_pose")
        result.update({"sequence_valid": True, "coordinates_valid": False, "diagnostic_error": "missing heavy-atom coordinates"})
        result["primary_failure_category"] = categories[0]
        return result
    coordinates_valid = all(
        math.isfinite(value)
        for atom in (*receptor_heavy, *peptide_atoms)
        for value in (atom.x, atom.y, atom.z)
    )
    if not coordinates_valid:
        categories.append("malformed_pose")
    result.update({
        "sequence_valid": True,
        "coordinates_valid": coordinates_valid,
        "peptide_chain": peptide_residues[0].chain,
        "peptide_atom_count_heavy": len(peptide_atoms),
        "peptide_residue_count": len(peptide_residues),
        "peptide_sequence_observed": "".join(residue.one_letter or "X" for residue in peptide_residues),
    })
    peptide_missing_backbone = [
        {
            "residue": f"{residue.resname}{residue.resseq}{residue.icode}",
            "missing_atoms": [name for name in BACKBONE_ATOMS if name not in residue.atom_names],
        }
        for residue in peptide_residues
        if any(name not in residue.atom_names for name in BACKBONE_ATOMS)
    ]
    result["peptide_missing_backbone_atoms"] = peptide_missing_backbone
    if peptide_missing_backbone:
        categories.append("malformed_pose")

    target_atoms_by_residue = _mapped_target_atoms(receptor_atoms, mapping, target_author_residues)
    target_metrics: dict[str, object] = {}
    for requested in target_author_residues:
        anchor_atoms = target_atoms_by_residue[requested]
        guanidinium_atoms = [atom for atom in anchor_atoms if atom.name in ARG_SIDECHAIN_ATOMS]
        distances_by_peptide_atom = {
            id(peptide_atom): min(_distance(target_atom, peptide_atom) for target_atom in anchor_atoms)
            for peptide_atom in peptide_atoms
        }
        guanidinium_distances = {
            id(peptide_atom): min(_distance(target_atom, peptide_atom) for target_atom in guanidinium_atoms)
            for peptide_atom in peptide_atoms
        }
        per_cutoff: dict[str, object] = {}
        for cutoff in contact_cutoffs_A:
            contacting_atoms = [atom for atom in peptide_atoms if distances_by_peptide_atom[id(atom)] <= cutoff]
            contacting_ids = {id(atom) for atom in contacting_atoms}
            guanidinium_contacting = [atom for atom in peptide_atoms if guanidinium_distances[id(atom)] <= cutoff]
            per_cutoff[str(cutoff)] = {
                "peptide_atom_count_within_A": len(contacting_atoms),
                "peptide_residues": sorted({
                    f"{atom.resname}{atom.resseq}{atom.icode}" for atom in contacting_atoms
                }),
                "peptide_residue_count": len({
                    (atom.chain, atom.resseq, atom.icode) for atom in contacting_atoms
                }),
                "guanidinium_peptide_atom_count_within_A": len(guanidinium_contacting),
                "guanidinium_peptide_residue_count": len({
                    (atom.chain, atom.resseq, atom.icode) for atom in guanidinium_contacting
                }),
                "peptide_atom_ids_within_A": len(contacting_ids),
            }
        target_metrics[str(requested)] = {
            "anchor_residue_number": requested,
            "anchor_atom_count": len(anchor_atoms),
            "guanidinium_atom_names": [atom.name for atom in guanidinium_atoms],
            "minimum_peptide_heavy_distance_A": round(min(distances_by_peptide_atom.values()), 3),
            "minimum_guanidinium_distance_A": round(min(guanidinium_distances.values()), 3),
            "contact_cutoffs": per_cutoff,
        }
    result["target_metrics"] = target_metrics
    result["target_contact_cutoff_A"] = structural_contact_cutoff_A
    result["target_contact_atom_count"] = sum(
        sum(
            1 for peptide_atom in peptide_atoms
            if min(_distance(target_atom, peptide_atom) for target_atom in target_atoms_by_residue[requested])
            <= structural_contact_cutoff_A
        )
        for requested in target_author_residues
    )
    minimum_distance, severe_clashes = _pairwise_min_and_clashes(
        receptor_heavy, peptide_atoms, severe_clash_cutoff,
    )
    result["minimum_receptor_peptide_heavy_atom_distance_A"] = round(minimum_distance, 3)
    result["severe_clash_count"] = severe_clashes
    result["pose_geometry_status"] = "valid" if coordinates_valid and severe_clashes == 0 else "severe_clashes_present"

    ca_atoms = [atom for atom in peptide_atoms if atom.name == "CA"]
    heavy_centroid = _centroid(peptide_atoms)
    ca_centroid = _centroid(ca_atoms)
    centroids = {"heavy_atom": heavy_centroid, "CA": ca_centroid}
    result["peptide_centroids_A"] = {key: list(value) if value is not None else None for key, value in centroids.items()}
    if site_center is not None:
        result["centroid_distances_to_site_center_A"] = {
            key: round(_distance_to_point(value, site_center), 3) if value is not None else None
            for key, value in centroids.items()
        }
    if box_size_A is not None and site_center is not None:
        half = float(box_size_A) / 2.0
        center = tuple(float(value) for value in site_center)
        inside = lambda atom: all(abs(getattr(atom, axis) - center[i]) <= half for i, axis in enumerate(("x", "y", "z")))
        inside_atoms = [atom for atom in peptide_atoms if inside(atom)]
        result["search_volume"] = {
            "definition": "axis-aligned cube centered at site_center_A",
            "half_edge_A": half,
            "any_peptide_heavy_atom_inside": bool(inside_atoms),
            "fraction_peptide_heavy_atoms_inside": round(len(inside_atoms) / len(peptide_atoms), 6),
            "all_peptide_heavy_atoms_inside": len(inside_atoms) == len(peptide_atoms),
            "CA_centroid_inside": bool(ca_centroid and all(abs(ca_centroid[i] - center[i]) <= half for i in range(3))),
            "heavy_atom_centroid_inside": bool(heavy_centroid and all(abs(heavy_centroid[i] - center[i]) <= half for i in range(3))),
        }
        # A few atoms can cross a large cube even when the pose centroid is
        # on a different receptor surface.  Use both centroids for the
        # failure classification and retain the atom-level inclusion metrics
        # above for inspection.
        if not result["search_volume"]["CA_centroid_inside"] and not result["search_volume"]["heavy_atom_centroid_inside"]:
            categories.append("search_region_miss")
    elif site_center is not None:
        result["search_volume"] = {"definition": "site center supplied without box size"}

    if coordinate_frame_result is not None:
        result["coordinate_frame"] = coordinate_frame_result
        if not coordinate_frame_result["coordinate_frame_comparable"]:
            categories.append("coordinate_mismatch")
    elif reference_path is not None:
        try:
            frame = coordinate_frame_audit(
                [receptor_path],
                reference_path=reference_path,
                chain=chain,
                target_author_residues=target_author_residues,
            )["representations"][0]
            result["coordinate_frame"] = frame
            if not frame["coordinate_frame_comparable"]:
                categories.append("coordinate_mismatch")
        except (OSError, ValueError, TargetMappingError) as exc:
            result["coordinate_frame"] = {"coordinate_frame_error": str(exc)}
            categories.append("coordinate_mismatch")
    if severe_clashes:
        categories.append("severe_clash")
    has_structural_contacts = any(
        metrics["minimum_peptide_heavy_distance_A"] <= structural_contact_cutoff_A
        for metrics in target_metrics.values()
    )
    result["target_contact_evidence_status"] = "contacts_present" if has_structural_contacts else "no_target_contacts"
    if not has_structural_contacts and not any(category in categories for category in ("wrong_mapping", "malformed_pose", "sequence_mismatch")):
        categories.append("geometrically_valid_but_misses_target_contacts")
    result["failure_categories"] = categories
    result["primary_failure_category"] = categories[0] if categories else "target_contacts_present"
    return result


def evaluate_pose(
    receptor_path: str | Path,
    pose_path: str | Path,
    peptide_sequence: str,
    *,
    target_author_residues: Sequence[int] = (203, 205),
    chain: str = "A",
    reference_path: str | Path | None = None,
    contact_cutoff: float = 4.5,
    severe_clash_cutoff: float = 2.0,
) -> dict[str, object]:
    """Independently validate coordinates, sequence, target contacts, clashes."""

    receptor_atoms, _ = parse_pdb(receptor_path)
    pose_atoms, _ = parse_pdb(pose_path)
    mapping = audit_residue_mapping(
        receptor_path,
        chain=chain,
        target_author_residues=target_author_residues,
        reference_path=reference_path,
    )
    if not mapping["mapping_valid"]:
        raise TargetMappingError(
            f"target mapping failed for {receptor_path}: {mapping['failure_reasons']}"
        )
    peptide_residues = _find_peptide_residues(pose_atoms, peptide_sequence)
    peptide_atoms = [atom for residue in peptide_residues for atom in residue.atoms if not _is_hydrogen(atom)]
    receptor_heavy = [atom for atom in receptor_atoms if not _is_hydrogen(atom)]
    if not peptide_atoms or not receptor_heavy:
        raise PoseValidationError("pose or receptor has no heavy-atom coordinates")
    target_atoms: list[Atom] = []
    target_record_by_number = {
        record["requested_author_residue_number"]: record
        for record in mapping["target_residues"]
    }
    for requested in target_author_residues:
        record = target_record_by_number[requested]
        observed_key = (
            record["chain"],
            record["pdb_author_residue_number"],
            record["pdb_author_insertion_code"],
        )
        target_atoms.extend(
            atom for atom in receptor_atoms
            if atom.residue_key == observed_key and not _is_hydrogen(atom)
        )
    contact_rows: list[dict[str, object]] = []
    all_contacts = 0
    for target_atom in target_atoms:
        distances = [(_distance(target_atom, peptide_atom), peptide_atom) for peptide_atom in peptide_atoms]
        distance, peptide_atom = min(distances, key=lambda item: item[0])
        if distance <= contact_cutoff:
            all_contacts += 1
            contact_rows.append({
                "target_residue_number": next(
                    key for key, record in target_record_by_number.items()
                    if record["pdb_author_residue_number"] == target_atom.resseq
                ),
                "target_atom": target_atom.name,
                "peptide_chain": peptide_atom.chain,
                "peptide_residue_number": peptide_atom.resseq,
                "peptide_residue_name": peptide_atom.resname,
                "peptide_atom": peptide_atom.name,
                "distance_A": round(distance, 3),
            })
    clash_rows: list[dict[str, object]] = []
    min_distance = math.inf
    for receptor_atom in receptor_heavy:
        for peptide_atom in peptide_atoms:
            distance = _distance(receptor_atom, peptide_atom)
            min_distance = min(min_distance, distance)
            if distance < severe_clash_cutoff:
                clash_rows.append({
                    "receptor_residue": f"{receptor_atom.resname}{receptor_atom.resseq}{receptor_atom.icode}",
                    "receptor_atom": receptor_atom.name,
                    "peptide_residue": f"{peptide_atom.resname}{peptide_atom.resseq}{peptide_atom.icode}",
                    "peptide_atom": peptide_atom.name,
                    "distance_A": round(distance, 3),
                })
    return {
        "receptor_path": str(receptor_path),
        "pose_path": str(pose_path),
        "peptide_sequence": peptide_sequence,
        "target_mapping_valid": True,
        "pose_coordinates_valid": all(
            math.isfinite(value)
            for atom in (*receptor_heavy, *peptide_atoms)
            for value in (atom.x, atom.y, atom.z)
        ),
        "peptide_sequence_valid": True,
        "peptide_chain": peptide_residues[0].chain,
        "peptide_atom_count_heavy": len(peptide_atoms),
        "target_contact_cutoff_A": contact_cutoff,
        "target_contact_atom_count": all_contacts,
        "target_contacts": contact_rows,
        "minimum_receptor_peptide_heavy_atom_distance_A": round(min_distance, 3),
        "severe_clash_cutoff_A": severe_clash_cutoff,
        "severe_clash_count": len(clash_rows),
        "severe_clashes": clash_rows,
        "pose_geometry_status": "valid" if not clash_rows else "severe_clashes_present",
        "target_contact_evidence_status": "contacts_present" if all_contacts else "no_target_contacts",
    }


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    mapping_parser = subparsers.add_parser("mapping")
    mapping_parser.add_argument("--reference", required=True, type=Path)
    mapping_parser.add_argument("--receptor", action="append", required=True, type=Path)
    mapping_parser.add_argument("--chain", default="A")
    mapping_parser.add_argument("--mmcif", type=Path)
    mapping_parser.add_argument("--output", required=True, type=Path)

    repair_parser = subparsers.add_parser("repair")
    repair_parser.add_argument("--input", required=True, type=Path)
    repair_parser.add_argument("--reference", required=True, type=Path)
    repair_parser.add_argument("--output", required=True, type=Path)
    repair_parser.add_argument("--report", required=True, type=Path)
    repair_parser.add_argument("--chain", default="A")

    pose_parser = subparsers.add_parser("pose")
    pose_parser.add_argument("--receptor", required=True, type=Path)
    pose_parser.add_argument("--pose", required=True, type=Path)
    pose_parser.add_argument("--sequence", required=True)
    pose_parser.add_argument("--reference", required=True, type=Path)
    pose_parser.add_argument("--output", required=True, type=Path)
    pose_parser.add_argument("--chain", default="A")

    args = parser.parse_args(argv)
    if args.command == "mapping":
        result = compare_receptor_mappings(
            args.receptor,
            reference_path=args.reference,
            chain=args.chain,
            mmcif_path=args.mmcif,
        )
        _write_json(args.output, result)
        print(json.dumps({"all_representations_valid": result["all_representations_valid"], "output": str(args.output)}, indent=2))
        return 0 if result["all_representations_valid"] else 2
    if args.command == "repair":
        result = repair_from_reference(args.input, args.reference, args.output, chain=args.chain)
        _write_json(args.report, result)
        print(json.dumps({"preparation_status": result["preparation_status"], "added_atoms": len(result["added_atoms"]), "output": str(args.output)}, indent=2))
        return 0 if result["target_mapping_valid"] else 2
    result = evaluate_pose(args.receptor, args.pose, args.sequence, reference_path=args.reference, chain=args.chain)
    _write_json(args.output, result)
    print(json.dumps({"pose_geometry_status": result["pose_geometry_status"], "target_contact_atom_count": result["target_contact_atom_count"], "output": str(args.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
