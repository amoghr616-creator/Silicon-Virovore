from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import tarfile
import json
import uuid
from datetime import datetime, timezone

from pathlib import Path
from typing import Any

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def pdb_chain_sequences(data: bytes) -> dict[str, str]:
    residues: dict[str, dict[tuple[int, str], str]] = {}

    for line in data.decode("utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 27:
            continue

        if line[16:17].strip() not in ("", "A"):
            continue

        chain = line[21:22].strip()
        resname = line[17:20].strip().upper()

        if resname not in AA3:
            continue

        try:
            number = int(line[22:26].strip())
        except ValueError:
            continue

        icode = line[26:27].strip()
        residues.setdefault(chain, {})[(number, icode)] = AA3[resname]

    return {
        chain: "".join(
            residue_map[key] for key in sorted(residue_map)
        )
        for chain, residue_map in residues.items()
    }



def extract_receptor_chain_a(
    source_pdb: Path,
    destination: Path,
) -> Path:
    """Create a chain-A-only receptor from the bundled official 9SYA PDB.

    The result is accepted only when its parsed chain-A sequence exactly
    matches chain A in the source. This prevents an arbitrary fallback file
    from silently becoming the HADDOCK receptor.
    """
    source_pdb = Path(source_pdb).resolve()
    destination = Path(destination).resolve()
    if not source_pdb.is_file():
        raise FileNotFoundError(source_pdb)

    source_bytes = source_pdb.read_bytes()
    source_sequences = pdb_chain_sequences(source_bytes)
    expected = source_sequences.get("A")
    if not expected:
        raise ValueError(
            f"Fallback 9SYA file has no parseable chain A: {source_pdb}"
        )

    lines = []
    for line in source_bytes.decode("utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") and len(line) >= 22 and line[21:22] == "A":
            lines.append(line)

    if not lines:
        raise ValueError(f"No ATOM records found for chain A in {source_pdb}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines + ["TER", "END", ""]) , encoding="utf-8")

    observed = pdb_chain_sequences(destination.read_bytes())
    if observed != {"A": expected}:
        destination.unlink(missing_ok=True)
        raise ValueError(
            "Extracted fallback receptor chain failed exact sequence validation."
        )
    return destination

def normalized_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper()




def _pdb_chain_segment_identifiers(pdb_path: Path):
    """Read chain and segment identifiers from PDB coordinate records."""
    chains, segments = set(), set()

    for line in Path(pdb_path).read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            padded = line.ljust(76)
            chain = padded[21].strip()
            segment = padded[72:76].strip()[:1]
            if chain:
                chains.add(chain)
            if segment:
                segments.add(segment)

    return chains, segments


def prepare_peptide_input_for_haddock(
    source_pdb: Path,
    destination: Path,
    *,
    expected_sequence: str,
) -> Path:
    """Write a separate chain-B/segment-B peptide PDB without changing coordinates."""
    source_pdb = Path(source_pdb).resolve()
    destination = Path(destination).resolve()

    if not source_pdb.is_file():
        raise FileNotFoundError(source_pdb)
    if source_pdb == destination:
        raise ValueError("Refusing to overwrite the original peptide structure.")

    observed = pdb_chain_sequences(source_pdb.read_bytes())
    if len(observed) != 1 or next(iter(observed.values())) != expected_sequence:
        raise ValueError(
            f"Peptide sequence/chain validation failed: {observed}; "
            f"expected one chain with {expected_sequence}"
        )

    output_lines = []
    atom_count = 0

    for raw in source_pdb.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        record = raw[:6].strip().upper()

        if record in {"ATOM", "HETATM", "ANISOU", "TER"}:
            chars = list(raw.ljust(76))
            chars[21] = "B"
            chars[72:76] = list("B   ")
            output_lines.append("".join(chars).rstrip())

            if record in {"ATOM", "HETATM"}:
                atom_count += 1
        else:
            output_lines.append(raw)

    if atom_count == 0:
        raise ValueError("Peptide structure contains no coordinate records.")

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        "\n".join(output_lines) + "\n",
        encoding="utf-8",
    )

    verified = pdb_chain_sequences(destination.read_bytes())
    if verified != {"B": expected_sequence}:
        destination.unlink(missing_ok=True)
        raise ValueError(f"Prepared peptide sequence validation failed: {verified}")

    chains, segments = _pdb_chain_segment_identifiers(destination)
    if chains != {"B"} or segments != {"B"}:
        destination.unlink(missing_ok=True)
        raise ValueError(
            f"Expected chain B/segment B; got chains={chains}, segments={segments}"
        )

    return destination


def validate_haddock_identifier_uniqueness(
    receptor_pdb: Path,
    peptide_pdb: Path,
) -> None:
    """Check the combined chain/segment identifier namespace before HADDOCK."""
    receptor_chains, receptor_segments = (
        _pdb_chain_segment_identifiers(receptor_pdb)
    )
    peptide_chains, peptide_segments = (
        _pdb_chain_segment_identifiers(peptide_pdb)
    )

    receptor_ids = receptor_chains | receptor_segments
    peptide_ids = peptide_chains | peptide_segments
    shared = receptor_ids & peptide_ids

    if not receptor_ids or not peptide_ids or shared:
        raise ValueError(
            "chain_identifier_collision: "
            f"receptor IDs={sorted(receptor_ids)}, "
            f"peptide IDs={sorted(peptide_ids)}, "
            f"shared={sorted(shared)}"
        )


def make_config(
    template: str,
    *,
    run_dir: Path,
    receptor: Path,
    peptide: Path,
    seed: int,
    ncores: int,
    peptide_length: int,
) -> str:
    text, n = re.subn(
        r'(?m)^run_dir\s*=.*$',
        f'run_dir = "{run_dir.resolve()}"',
        template,
        count=1,
    )
    if n != 1:
        raise ValueError("Template must contain exactly one run_dir.")

    text, n = re.subn(
        r"(?ms)^molecules\s*=\s*\[.*?^\]",
        (
            "molecules = [\n"
            f'    "{receptor.resolve()}",\n'
            f'    "{peptide.resolve()}",\n'
            "]"
        ),
        text,
        count=1,
    )
    if n != 1:
        raise ValueError("Could not replace template molecule list.")

    text = re.sub(
        r"(?m)^ncores\s*=.*$",
        f"ncores = {int(ncores)}",
        text,
    )

    text, seed_count = re.subn(
        r"(?m)^(\s*)iniseed\s*=\s*\d+\s*$",
        lambda m: f"{m.group(1)}iniseed = {int(seed)}",
        text,
    )
    if seed_count == 0:
        raise ValueError("Template contains no iniseed settings.")

    text, flex_count = re.subn(
        r"(?m)^(\s*fle_end_1\s*=\s*)\d+\s*$",
        rf"\g<1>{int(peptide_length)}",
        text,
    )
    if flex_count == 0:
        raise ValueError("Template contains no fle_end_1 settings.")

    return text


def extract_final_models(
    archive: Path,
    destination: Path,
    *,
    receptor_sequence: str,
    candidate_sequence: str,
) -> list[Path]:
    destination.mkdir(parents=True, exist_ok=True)
    valid = []

    with tarfile.open(archive, "r:*") as tar:
        members = sorted(
            (
                m for m in tar.getmembers()
                if m.isfile()
                and re.fullmatch(
                    r"cluster_\d+_model_1\.pdb",
                    Path(m.name).name,
                )
            ),
            key=lambda m: Path(m.name).name,
        )

        for member in members:
            handle = tar.extractfile(member)
            if handle is None:
                continue

            content = handle.read()
            sequences = pdb_chain_sequences(content)

            # Strict validation: chain A must be 9SYA receptor;
            # chain B must match this candidate peptide exactly.
            if sequences.get("A") != receptor_sequence:
                print(
                    "REJECTED receptor sequence mismatch:",
                    member.name,
                    sequences.get("A", "<missing>"),
                )
                continue

            if sequences.get("B") != candidate_sequence:
                print(
                    "REJECTED peptide sequence mismatch:",
                    member.name,
                    sequences.get("B", "<missing>"),
                )
                continue

            path = destination / Path(member.name).name
            path.write_bytes(content)
            valid.append(path)

    return valid


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _valid_execution_manifest(
    manifest_path: Path,
    *,
    seed: int,
    protocol_sha256: str,
    receptor_sha256: str,
    peptide_sha256: str,
    archive_path: Path,
) -> dict[str, Any] | None:
    """Return a prior execution only when its provenance is independently checkable."""
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        if payload.get("status") != "success" or payload.get("return_code") != 0:
            return None
        pose_hashes = payload.get("verified_pose_sha256")
        if not isinstance(pose_hashes, dict) or not pose_hashes:
            return None
        if any(
            not isinstance(value, str)
            or len(value) != 64
            or any(ch not in "0123456789abcdef" for ch in value.lower())
            for value in pose_hashes.values()
        ):
            return None
        expected = {
            "random_seed": int(seed),
            "protocol_sha256": protocol_sha256,
            "receptor_sha256": receptor_sha256,
            "peptide_input_sha256": peptide_sha256,
        }
        if any(payload.get(key) != value for key, value in expected.items()):
            return None
        if not payload.get("execution_id") or not payload.get("independent_run_id"):
            return None
        if not archive_path.is_file():
            return None
        if payload.get("archive_sha256") != _sha256_file(archive_path):
            return None
        return payload
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return None


def generate_candidate_target_complexes(
    candidates: list[Any],
    *,
    project_root: Path,
    template_config: Path,
    receptor_pdb: Path,
    output_root: Path,
    max_candidates: int = 1,
    seeds: tuple[int, ...] = (1917, 2917),
    ncores: int = 4,
) -> dict[str, Any]:
    """Generate exact-sequence 9SYA complexes with auditable independent runs.

    Each seed is a separate HADDOCK3 process. The per-seed execution manifest
    records the effective protocol hash, input hashes, command, archive hash,
    pose hashes, exit status, and unique execution ID. Existing output is only
    reused when this manifest verifies; an old/incomplete directory is never
    silently treated as an independent successful execution.

    Candidate-specific failures are recorded and do not abort evaluation of
    the remaining shortlist. Geometry/docking output is not binding affinity.
    """
    executable = shutil.which("haddock3")
    if executable is None:
        raise RuntimeError("haddock3 executable not found in the active environment.")
    if not template_config.is_file():
        raise FileNotFoundError(template_config)
    if not isinstance(max_candidates, int) or max_candidates < 1:
        raise ValueError("max_candidates must be an integer >= 1")
    if not isinstance(ncores, int) or ncores < 1:
        raise ValueError("ncores must be an integer >= 1")
    normalized_seeds = tuple(int(seed) for seed in seeds)
    if len(normalized_seeds) < 2 or len(set(normalized_seeds)) != len(normalized_seeds):
        raise ValueError("At least two distinct HADDOCK seeds are required for repeatability.")
    if any(seed < 0 for seed in normalized_seeds):
        raise ValueError("HADDOCK seeds must be nonnegative integers.")

    project_root = Path(project_root).resolve()
    template_config = Path(template_config).resolve()
    receptor_pdb = Path(receptor_pdb).expanduser().resolve()
    output_root = Path(output_root).expanduser().resolve()

    if not receptor_pdb.is_file():
        bundled_9sya = project_root / "data/receptor/9SYA.pdb"
        if not bundled_9sya.is_file():
            raise FileNotFoundError(
                "HADDOCK receptor input is missing. Expected either "
                f"{receptor_pdb} or bundled official structure {bundled_9sya}."
            )
        receptor_pdb = extract_receptor_chain_a(
            bundled_9sya,
            output_root / "_receptor/9sya_source_chain_A.pdb",
        )

    template = template_config.read_text(encoding="utf-8")
    receptor_sequences = pdb_chain_sequences(receptor_pdb.read_bytes())
    receptor_sequence = receptor_sequences.get("A")
    if not receptor_sequence:
        raise ValueError("Could not parse receptor chain A.")

    receptor_sha256 = _sha256_file(receptor_pdb)
    template_sha256 = _sha256_file(template_config)
    output_root.mkdir(parents=True, exist_ok=True)

    summary: dict[str, Any] = {
        "candidates_requested": min(len(candidates), max_candidates),
        "candidates_attempted": 0,
        "candidates_with_two_independent_runs": 0,
        "candidates_with_any_verified_pose": 0,
        "verified_pose_count": 0,
        "independent_run_count": 0,
        "failed_candidate_count": 0,
        "candidate_results": [],
        "reliability_calibrated": False,
        "binding_affinity_claim": False,
    }

    considered = 0
    for candidate in candidates:
        metadata = candidate.metadata if isinstance(candidate.metadata, dict) else {}
        candidate.metadata = metadata
        if metadata.get("target_complexes"):
            summary["candidate_results"].append({
                "candidate_id": metadata.get("candidate_id"),
                "sequence": normalized_sequence(candidate.sequence),
                "status": "existing_complexes_not_regenerated",
                "registered_complex_count": len(metadata.get("target_complexes") or []),
            })
            continue
        if considered >= max_candidates:
            break
        considered += 1
        summary["candidates_attempted"] += 1

        sequence = normalized_sequence(candidate.sequence)
        candidate_id = str(metadata.get("candidate_id") or sha256_bytes(sequence.encode())[:12])
        safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", candidate_id)
        candidate_result: dict[str, Any] = {
            "candidate_id": candidate_id,
            "sequence": sequence,
            "sequence_sha256": sha256_bytes(sequence.encode()),
            "status": "started",
            "seed_results": [],
        }
        candidate_entries: list[dict[str, Any]] = []
        candidate_runs = 0

        try:
            peptide_path_value = getattr(candidate, "structure_path", None)
            if not peptide_path_value or not Path(peptide_path_value).is_file():
                raise FileNotFoundError("Candidate-specific peptide structure is unavailable.")
            peptide_path = Path(peptide_path_value).expanduser().resolve()
            peptide_sequences = pdb_chain_sequences(peptide_path.read_bytes())
            if sequence not in peptide_sequences.values():
                raise ValueError(
                    f"Peptide structure sequence mismatch: observed={peptide_sequences}; "
                    f"expected_sha256={sha256_bytes(sequence.encode())}"
                )

            haddock_peptide_path = output_root / "_peptide_inputs" / f"{safe_id}_chain_B.pdb"
            peptide_path = prepare_peptide_input_for_haddock(
                peptide_path, haddock_peptide_path, expected_sequence=sequence
            )
            validate_haddock_identifier_uniqueness(receptor_pdb, peptide_path)
            peptide_sha256 = _sha256_file(peptide_path)

            for seed in normalized_seeds:
                seed_root = output_root / safe_id / f"seed_{seed}"
                run_dir = seed_root / "haddock_run"
                config_dir = output_root / "configs"
                config_dir.mkdir(parents=True, exist_ok=True)
                config_path = config_dir / f"{safe_id}_seed{seed}.cfg"
                log_path = config_dir / f"{safe_id}_seed{seed}.log"
                execution_manifest_path = seed_root / "execution_manifest.json"

                config = make_config(
                    template,
                    run_dir=run_dir,
                    receptor=receptor_pdb,
                    peptide=peptide_path,
                    seed=seed,
                    ncores=ncores,
                    peptide_length=len(sequence),
                )
                normalized_config = re.sub(
                    r'(?m)^run_dir\s*=.*$', 'run_dir = "<RUN_DIR>"', config
                )
                normalized_config = re.sub(
                    r"(?m)^(\s*iniseed\s*=\s*)\d+\s*$",
                    r"\g<1><SEED>", normalized_config,
                )
                protocol_sha256 = sha256_bytes(
                    (normalized_config + "\nTEMPLATE_SHA256=" + template_sha256).encode()
                )
                archive = run_dir / "analysis/6_seletopclusts_analysis/summary.tgz"

                prior = _valid_execution_manifest(
                    execution_manifest_path,
                    seed=seed,
                    protocol_sha256=protocol_sha256,
                    receptor_sha256=receptor_sha256,
                    peptide_sha256=peptide_sha256,
                    archive_path=archive,
                )
                if prior is not None:
                    models = extract_final_models(
                        archive,
                        seed_root / "complexes",
                        receptor_sequence=receptor_sequence,
                        candidate_sequence=sequence,
                    )
                    if not models:
                        candidate_result["seed_results"].append({
                            "seed": seed,
                            "status": "reused_archive_failed_sequence_validation",
                        })
                        continue
                    manifest_sha256 = _sha256_file(execution_manifest_path)
                    for model_path in models:
                        candidate_entries.append({
                            "path": str(model_path.resolve()),
                            "sha256": _sha256_file(model_path),
                            "target_structure_id": "9SYA",
                            "receptor_chain": "A",
                            "peptide_chain": "B",
                            "provenance": "candidate_specific_HADDOCK3_predicted_complex",
                            "independent_run_id": prior["independent_run_id"],
                            "random_seed": int(seed),
                            "protocol_sha256": protocol_sha256,
                            "execution_id": prior["execution_id"],
                            "execution_manifest_path": str(execution_manifest_path.resolve()),
                            "execution_manifest_sha256": manifest_sha256,
                        })
                    candidate_runs += 1
                    summary["independent_run_count"] += 1
                    candidate_result["seed_results"].append({
                        "seed": seed,
                        "status": "reused_verified_execution",
                        "execution_id": prior["execution_id"],
                        "pose_count": len(models),
                    })
                    continue

                # Never treat a pre-existing unmanifested run directory as new
                # evidence. It may be a partial or copied execution.
                if run_dir.exists() or execution_manifest_path.exists():
                    candidate_result["seed_results"].append({
                        "seed": seed,
                        "status": "unverified_existing_execution_directory",
                        "reason": "Existing HADDOCK output lacks a matching execution manifest; preserve for manual review.",
                    })
                    continue

                run_id = f"{safe_id}:HADDOCK3:seed{seed}"
                execution_id = uuid.uuid4().hex
                started_at = datetime.now(timezone.utc).isoformat()
                config_path.write_text(config, encoding="utf-8")
                config_sha256 = _sha256_file(config_path)
                command = [executable, str(config_path)]
                result = subprocess.run(
                    command,
                    cwd=project_root,
                    capture_output=True,
                    text=True,
                )
                log_text = result.stdout + "\n--- STDERR ---\n" + result.stderr
                log_path.write_text(log_text, encoding="utf-8")
                log_sha256 = _sha256_file(log_path)

                if result.returncode != 0 or not archive.is_file():
                    failed_manifest = {
                        "schema_version": 1,
                        "status": "failed",
                        "execution_id": execution_id,
                        "independent_run_id": run_id,
                        "candidate_id": candidate_id,
                        "sequence_sha256": sha256_bytes(sequence.encode()),
                        "random_seed": int(seed),
                        "protocol_sha256": protocol_sha256,
                        "template_sha256": template_sha256,
                        "config_sha256": config_sha256,
                        "receptor_sha256": receptor_sha256,
                        "peptide_input_sha256": peptide_sha256,
                        "command": command,
                        "started_at_utc": started_at,
                        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                        "return_code": int(result.returncode),
                        "log_path": str(log_path.resolve()),
                        "log_sha256": log_sha256,
                        "archive_path": str(archive.resolve()),
                        "failure": "HADDOCK process failed or final archive is missing.",
                    }
                    _atomic_json(execution_manifest_path, failed_manifest)
                    candidate_result["seed_results"].append({
                        "seed": seed,
                        "status": "failed",
                        "return_code": int(result.returncode),
                        "execution_manifest": str(execution_manifest_path),
                    })
                    continue

                models = extract_final_models(
                    archive,
                    seed_root / "complexes",
                    receptor_sequence=receptor_sequence,
                    candidate_sequence=sequence,
                )
                if not models:
                    failed_manifest = {
                        "schema_version": 1,
                        "status": "failed",
                        "execution_id": execution_id,
                        "independent_run_id": run_id,
                        "candidate_id": candidate_id,
                        "sequence_sha256": sha256_bytes(sequence.encode()),
                        "random_seed": int(seed),
                        "protocol_sha256": protocol_sha256,
                        "template_sha256": template_sha256,
                        "config_sha256": config_sha256,
                        "receptor_sha256": receptor_sha256,
                        "peptide_input_sha256": peptide_sha256,
                        "command": command,
                        "started_at_utc": started_at,
                        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                        "return_code": int(result.returncode),
                        "log_path": str(log_path.resolve()),
                        "log_sha256": log_sha256,
                        "archive_path": str(archive.resolve()),
                        "archive_sha256": _sha256_file(archive),
                        "failure": "Final archive contained no exact-sequence receptor/peptide models.",
                    }
                    _atomic_json(execution_manifest_path, failed_manifest)
                    candidate_result["seed_results"].append({
                        "seed": seed,
                        "status": "no_sequence_verified_models",
                        "execution_manifest": str(execution_manifest_path),
                    })
                    continue

                pose_hashes = {model_path.name: _sha256_file(model_path) for model_path in models}
                success_manifest = {
                    "schema_version": 1,
                    "status": "success",
                    "execution_id": execution_id,
                    "independent_run_id": run_id,
                    "candidate_id": candidate_id,
                    "sequence_sha256": sha256_bytes(sequence.encode()),
                    "random_seed": int(seed),
                    "protocol_sha256": protocol_sha256,
                    "template_sha256": template_sha256,
                    "config_sha256": config_sha256,
                    "receptor_sha256": receptor_sha256,
                    "peptide_input_sha256": peptide_sha256,
                    "command": command,
                    "started_at_utc": started_at,
                    "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                    "return_code": int(result.returncode),
                    "log_path": str(log_path.resolve()),
                    "log_sha256": log_sha256,
                    "archive_path": str(archive.resolve()),
                    "archive_sha256": _sha256_file(archive),
                    "verified_pose_sha256": pose_hashes,
                }
                _atomic_json(execution_manifest_path, success_manifest)
                manifest_sha256 = _sha256_file(execution_manifest_path)

                for model_path in models:
                    candidate_entries.append({
                        "path": str(model_path.resolve()),
                        "sha256": _sha256_file(model_path),
                        "target_structure_id": "9SYA",
                        "receptor_chain": "A",
                        "peptide_chain": "B",
                        "provenance": "candidate_specific_HADDOCK3_predicted_complex",
                        "independent_run_id": run_id,
                        "random_seed": int(seed),
                        "protocol_sha256": protocol_sha256,
                        "execution_id": execution_id,
                        "execution_manifest_path": str(execution_manifest_path.resolve()),
                        "execution_manifest_sha256": manifest_sha256,
                    })
                candidate_runs += 1
                summary["independent_run_count"] += 1
                summary["verified_pose_count"] += len(models)
                candidate_result["seed_results"].append({
                    "seed": seed,
                    "status": "success",
                    "execution_id": execution_id,
                    "pose_count": len(models),
                    "execution_manifest": str(execution_manifest_path),
                })

            if candidate_entries:
                metadata["target_complexes"] = candidate_entries
                metadata["target_structure_id"] = "9SYA"
                metadata["target_receptor_chain"] = "A"
                metadata["target_peptide_chain"] = "B"
                metadata["p4_complex_generation"] = {
                    "status": "success" if candidate_runs >= 2 else "partial_success",
                    "independent_run_count": candidate_runs,
                    "verified_pose_count": len(candidate_entries),
                    "sequence_verified": True,
                    "execution_manifests_verified": True,
                    "reliability_calibrated": False,
                }
                summary["candidates_with_any_verified_pose"] += 1
                if candidate_runs >= 2:
                    summary["candidates_with_two_independent_runs"] += 1
                    candidate_result["status"] = "two_independent_runs_verified"
                else:
                    candidate_result["status"] = "partial_target_evidence"
            else:
                metadata["p4_complex_generation"] = {
                    "status": "no_valid_complexes",
                    "independent_run_count": candidate_runs,
                    "sequence_verified": False,
                    "reliability_calibrated": False,
                }
                candidate_result["status"] = "no_valid_complexes"
        except Exception as exc:
            candidate_result["status"] = "candidate_generation_error"
            candidate_result["error"] = f"{type(exc).__name__}: {exc}"
            metadata["p4_complex_generation"] = {
                "status": "candidate_generation_error",
                "error": candidate_result["error"],
                "independent_run_count": candidate_runs,
                "reliability_calibrated": False,
            }

        if candidate_result["status"] not in {
            "two_independent_runs_verified",
            "partial_target_evidence",
            "no_valid_complexes",
        }:
            summary["failed_candidate_count"] += 1
        summary["candidate_results"].append(candidate_result)

    _atomic_json(output_root / "target_complex_generation_summary.json", summary)
    print("\nCandidate HADDOCK generation summary:", {
        key: value for key, value in summary.items() if key != "candidate_results"
    })
    return summary
