from pathlib import Path
from datetime import datetime
from collections import defaultdict, Counter
from difflib import SequenceMatcher
import os
import re
import csv
import io
import gzip
import json
import math
import hashlib
import random
import zipfile
import tarfile
import shutil

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "NumPy is required. Activate the project's .venv and install the "
        "project's existing scientific Python dependencies."
    ) from exc


ROOT = Path.cwd().resolve()
STAMP = datetime.now().strftime("%Y%m%dT%H%M%S")

PACKAGE = None
PACKAGE_CANDIDATES = sorted(
    p for p in (ROOT / "prototype3/logbook").glob("p3_final_defensible_package_*")
    if p.is_dir()
)

if PACKAGE_CANDIDATES:
    PACKAGE = PACKAGE_CANDIDATES[-1]
    OUT = PACKAGE / f"candidate_structural_translation_test_{STAMP}"
else:
    # The package may exist only as a ZIP; create a separate output folder.
    OUT = ROOT / "prototype3/logbook" / (
        f"p3_candidate_structural_translation_test_{STAMP}"
    )

TARGET = ROOT / "prototype2/target_aware_audit/9SYA_chainA_repaired.pdb"

if not TARGET.exists():
    matches = list(
        (ROOT / "prototype2").rglob("9SYA_chainA_repaired.pdb")
    )
    if matches:
        TARGET = matches[0]

# Predefined geometric success criteria.
MIN_INTERCHAIN_A = 2.0
MAX_ARG203_A = 5.0
MAX_ARG205_A = 5.0
MIN_TARGET_SEQUENCE_MATCH = 0.65

N_BOOTSTRAP = 5000
N_PERMUTATIONS = 10000
RANDOM_SEED = 20261027

SCAN_ROOTS = [
    ROOT / "prototype2",
    ROOT / "prototype3/results",
]

SKIP_DIRS = {
    ".git", ".venv", "venv", "__pycache__", "node_modules",
    "logbook", "source_archive", "source_snapshot",
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}
WATER_NAMES = {"HOH", "WAT", "DOD"}

MODEL_FIELDS = [
    "arm", "study_group", "cluster_id", "model_number",
    "source_file", "peptide_sequence", "peptide_length",
    "receptor_chain", "peptide_chain",
    "target_local_sequence_match_fraction",
    "min_interchain_heavy_atom_distance_A",
    "peptide_to_ARG203_min_distance_A",
    "peptide_to_ARG205_min_distance_A",
    "eligible_for_comparison", "passes_all_geometric_criteria", "status",
]


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def is_pdb_name(name):
    n = name.lower()
    return n.endswith((".pdb", ".ent", ".pdb.gz", ".ent.gz"))


def is_archive_name(name):
    n = name.lower()
    return n.endswith((".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2"))


def classify_arm(parts):
    labels = set()

    for part in parts:
        s = str(part).lower()
        adaptive = bool(
            re.search(r"adaptive|(?<![a-z])adapt(?![a-z])", s)
        )
        baseline = bool(re.search(r"baseline|uniform", s))

        # A shared directory named, for example, "adaptive_vs_baseline"
        # does not identify the arm of an individual model.
        if adaptive and baseline:
            continue
        if adaptive:
            labels.add("adaptive")
        elif baseline:
            labels.add("baseline")

    return next(iter(labels)) if len(labels) == 1 else "unknown"


def normalize_component(part):
    s = str(part)
    s = re.sub(
        r"uniform[_-]?baseline", "{ARM}", s, flags=re.IGNORECASE
    )
    s = re.sub(
        r"adaptive|baseline|uniform", "{ARM}", s, flags=re.IGNORECASE
    )
    return s


def study_group_from_parts(parts):
    # Drop the final PDB member/file name, then remove numbered cluster
    # directories so models from separate clusters remain in the same study.
    out = []
    for part in list(parts)[:-1]:
        if re.fullmatch(
            r"(?:cluster|clust)[_-]?\d+", str(part), flags=re.IGNORECASE
        ):
            continue
        out.append(normalize_component(part))

    return "/".join(out) or "unidentified_study_group"


def cluster_from_parts(parts):
    text = "/".join(map(str, parts))
    m = re.search(
        r"(?:cluster|clust)[ _-]*0*(\d+)", text, flags=re.IGNORECASE
    )
    return f"cluster_{int(m.group(1))}" if m else ""


def decode_pdb_bytes(data, member_name=""):
    if member_name.lower().endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


def read_path_text(path):
    if path.name.lower().endswith((".pdb.gz", ".ent.gz")):
        with gzip.open(path, "rt", errors="replace") as f:
            return f.read()
    return path.read_text(errors="replace")


def iter_arm_labeled_pdb_sources():
    """Yield existing arm-labeled PDB models; read archive members in place."""
    yielded = set()

    for base in SCAN_ROOTS:
        if not base.exists():
            continue

        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [
                d for d in dirnames
                if d not in SKIP_DIRS
                and not d.startswith("p3_final_defensible_package_")
                and not d.startswith("p3_comparison_full_audit_")
            ]

            parent = Path(dirpath)

            for filename in sorted(filenames):
                path = parent / filename
                rel = path.relative_to(ROOT)
                low = filename.lower()

                if is_pdb_name(filename):
                    parts = list(rel.parts)
                    arm = classify_arm(parts)
                    if arm == "unknown":
                        continue

                    try:
                        content = read_path_text(path)
                    except Exception:
                        continue

                    key = ("file", str(rel))
                    if key in yielded:
                        continue
                    yielded.add(key)

                    yield {
                        "source_file": str(rel),
                        "parts": parts,
                        "arm": arm,
                        "text": content,
                    }

                elif is_archive_name(filename):
                    # A combined archive may contain separate adaptive and
                    # baseline directories inside it, so classify each member.
                    try:
                        if low.endswith(".zip"):
                            with zipfile.ZipFile(path) as archive:
                                for info in archive.infolist():
                                    if info.is_dir() or not is_pdb_name(info.filename):
                                        continue
                                    try:
                                        raw = archive.read(info.filename)
                                        content = decode_pdb_bytes(
                                            raw, info.filename
                                        )
                                    except Exception:
                                        continue

                                    member_parts = info.filename.replace(
                                        "\\", "/"
                                    ).split("/")
                                    parts = list(rel.parts) + member_parts
                                    arm = classify_arm(parts)
                                    if arm == "unknown":
                                        continue

                                    label = f"{rel}::{info.filename}"
                                    key = ("archive_member", label)
                                    if key in yielded:
                                        continue
                                    yielded.add(key)

                                    yield {
                                        "source_file": label,
                                        "parts": parts,
                                        "arm": arm,
                                        "text": content,
                                    }

                        elif low.endswith((".tar", ".tar.gz", ".tgz", ".tar.bz2")):
                            with tarfile.open(path, "r:*") as archive:
                                for member in archive.getmembers():
                                    if not member.isfile() or not is_pdb_name(
                                        member.name
                                    ):
                                        continue
                                    try:
                                        stream = archive.extractfile(member)
                                        if stream is None:
                                            continue
                                        raw = stream.read()
                                        content = decode_pdb_bytes(
                                            raw, member.name
                                        )
                                    except Exception:
                                        continue

                                    member_parts = member.name.replace(
                                        "\\", "/"
                                    ).split("/")
                                    parts = list(rel.parts) + member_parts
                                    arm = classify_arm(parts)
                                    if arm == "unknown":
                                        continue

                                    label = f"{rel}::{member.name}"
                                    key = ("archive_member", label)
                                    if key in yielded:
                                        continue
                                    yielded.add(key)

                                    yield {
                                        "source_file": label,
                                        "parts": parts,
                                        "arm": arm,
                                        "text": content,
                                    }
                    except Exception:
                        # A corrupt/unreadable archive is excluded rather than
                        # silently interpreted as a negative structural result.
                        continue


def parse_models(text):
    models = []
    current = []

    for line in text.splitlines():
        record = line[:6].strip().upper()

        if record == "MODEL":
            if current:
                models.append(current)
                current = []
        elif record == "ENDMDL":
            if current:
                models.append(current)
                current = []
        elif record in {"ATOM", "HETATM"}:
            current.append(line)

    if current:
        models.append(current)

    return [m for m in models if m]


def parse_heavy_atoms(lines):
    atoms = []

    for line in lines:
        record = line[:6].strip().upper()
        if record not in {"ATOM", "HETATM"} or len(line) < 54:
            continue

        altloc = line[16:17].strip()
        if altloc not in {"", "A", "1"}:
            continue

        atom_name = line[12:16].strip()
        resname = line[17:20].strip().upper()
        chain = line[21:22].strip() or "_"
        icode = line[26:27].strip()

        if resname in WATER_NAMES:
            continue

        try:
            resseq = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue

        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            atom_letters = re.sub(r"^[0-9]+", "", atom_name)
            element = atom_letters[:1].upper()

        if element in {"H", "D"}:
            continue

        atoms.append({
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "resname": resname,
            "atom_name": atom_name,
            "element": element,
            "xyz": (x, y, z),
        })

    return atoms


def chain_sequence_and_map(atoms):
    residue_names = {}
    for a in atoms:
        if a["resname"] in AA3:
            key = (a["resseq"], a["icode"])
            residue_names[key] = a["resname"]

    keys = sorted(residue_names)
    sequence = "".join(AA3[residue_names[k]] for k in keys)
    return sequence, residue_names


def local_sequence_match_fraction(reference, candidate):
    if not reference or not candidate:
        return 0.0
    matcher = SequenceMatcher(None, reference, candidate, autojunk=False)
    matching_residues = sum(
        block.size for block in matcher.get_matching_blocks()
    )
    return matching_residues / min(len(reference), len(candidate))


def load_target_reference(path):
    if not path.exists():
        return "", None

    try:
        models = parse_models(path.read_text(errors="replace"))
        if not models:
            return "", None

        atoms = parse_heavy_atoms(models[0])
        chains = defaultdict(list)
        for atom in atoms:
            chains[atom["chain"]].append(atom)

        options = []
        for chain, chain_atoms in chains.items():
            seq, residue_map = chain_sequence_and_map(chain_atoms)
            if seq:
                options.append((len(seq), seq, residue_map, chain))

        if not options:
            return "", None

        _, seq, residue_map, _ = max(options, key=lambda item: item[0])
        return seq, residue_map

    except Exception:
        return "", None


def min_pair_distance(coords_a, coords_b):
    if len(coords_a) == 0 or len(coords_b) == 0:
        return float("nan")

    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    diff = a[:, None, :] - b[None, :, :]
    distances = np.sqrt(np.sum(diff * diff, axis=2))
    return float(np.min(distances))


def score_one_model(lines, target_sequence):
    atoms = parse_heavy_atoms(lines)
    if not atoms:
        return {
            "eligible": 0,
            "pass_geometry": "",
            "status": "no_heavy_atoms",
        }

    chains = defaultdict(list)
    for atom in atoms:
        chains[atom["chain"]].append(atom)

    # Select a receptor chain that contains both predefined ARG positions
    # and most closely matches the repaired 9SYA chain sequence.
    receptor_options = []
    for chain, chain_atoms in chains.items():
        arg203 = [
            a for a in chain_atoms
            if a["resseq"] == 203 and a["resname"] == "ARG"
        ]
        arg205 = [
            a for a in chain_atoms
            if a["resseq"] == 205 and a["resname"] == "ARG"
        ]
        if not arg203 or not arg205:
            continue

        seq, _ = chain_sequence_and_map(chain_atoms)
        match_fraction = local_sequence_match_fraction(
            target_sequence, seq
        ) if target_sequence else None

        receptor_options.append({
            "chain": chain,
            "atoms": chain_atoms,
            "arg203": arg203,
            "arg205": arg205,
            "sequence": seq,
            "match_fraction": match_fraction,
        })

    if not receptor_options:
        return {
            "eligible": 0,
            "pass_geometry": "",
            "status": "ARG203_or_ARG205_not_found_as_expected",
        }

    receptor_options.sort(
        key=lambda r: (
            -1.0 if r["match_fraction"] is None else r["match_fraction"],
            len(r["sequence"]),
        ),
        reverse=True,
    )
    receptor = receptor_options[0]

    peptide_options = []
    for chain, chain_atoms in chains.items():
        if chain == receptor["chain"]:
            continue

        seq, _ = chain_sequence_and_map(chain_atoms)
        if 5 <= len(seq) <= 50:
            peptide_options.append({
                "chain": chain,
                "atoms": chain_atoms,
                "sequence": seq,
            })

    if not peptide_options:
        return {
            "eligible": 0,
            "pass_geometry": "",
            "receptor_chain": receptor["chain"],
            "target_match_fraction": receptor["match_fraction"],
            "status": "no_plausible_peptide_chain_found",
        }

    # The target project uses 21-aa candidates. If multiple short chains
    # exist, prefer the chain whose canonical length is closest to 21 aa.
    peptide_options.sort(
        key=lambda p: (abs(len(p["sequence"]) - 21), len(p["sequence"]))
    )
    peptide = peptide_options[0]

    receptor_coords = [a["xyz"] for a in receptor["atoms"]]
    peptide_coords = [a["xyz"] for a in peptide["atoms"]]
    arg203_coords = [a["xyz"] for a in receptor["arg203"]]
    arg205_coords = [a["xyz"] for a in receptor["arg205"]]

    min_interchain = min_pair_distance(receptor_coords, peptide_coords)
    d203 = min_pair_distance(arg203_coords, peptide_coords)
    d205 = min_pair_distance(arg205_coords, peptide_coords)

    match_fraction = receptor["match_fraction"]
    if match_fraction is None:
        eligible = False
        status = "target_reference_missing"
    elif match_fraction < MIN_TARGET_SEQUENCE_MATCH:
        eligible = False
        status = "receptor_sequence_match_below_threshold"
    else:
        eligible = True
        status = "scored"

    passed = (
        min_interchain >= MIN_INTERCHAIN_A
        and d203 <= MAX_ARG203_A
        and d205 <= MAX_ARG205_A
    ) if eligible else ""

    return {
        "eligible": int(eligible),
        "pass_geometry": int(bool(passed)) if eligible else "",
        "status": status,
        "receptor_chain": receptor["chain"],
        "peptide_chain": peptide["chain"],
        "peptide_sequence": peptide["sequence"],
        "peptide_length": len(peptide["sequence"]),
        "target_match_fraction": match_fraction,
        "min_interchain": min_interchain,
        "d203": d203,
        "d205": d205,
    }


def bootstrap_ci(values, rng):
    if len(values) < 2:
        return None, None

    draws = []
    for _ in range(N_BOOTSTRAP):
        sample = [rng.choice(values) for _ in range(len(values))]
        draws.append(sum(sample) / len(sample))

    draws.sort()
    return (
        draws[int(0.025 * (len(draws) - 1))],
        draws[int(0.975 * (len(draws) - 1))],
    )


def cluster_rates_for(records):
    grouped = defaultdict(list)
    for row in records:
        if row["eligible_for_comparison"] != 1:
            continue
        if not row["cluster_id"]:
            continue
        grouped[row["cluster_id"]].append(
            int(row["passes_all_geometric_criteria"])
        )

    return {
        cluster: sum(values) / len(values)
        for cluster, values in grouped.items()
    }


def cluster_comparison(baseline_rates, adaptive_rates):
    if len(baseline_rates) < 3 or len(adaptive_rates) < 3:
        return {
            "difference": "",
            "ci_low": "",
            "ci_high": "",
            "p_two_sided": "",
            "status": "descriptive_only; fewer than 3 identified clusters per arm",
        }

    rng = random.Random(RANDOM_SEED)
    base = list(baseline_rates.values())
    adapt = list(adaptive_rates.values())

    observed = sum(adapt) / len(adapt) - sum(base) / len(base)

    boot_differences = []
    for _ in range(N_BOOTSTRAP):
        b = [rng.choice(base) for _ in base]
        a = [rng.choice(adapt) for _ in adapt]
        boot_differences.append(
            sum(a) / len(a) - sum(b) / len(b)
        )

    boot_differences.sort()
    lo = boot_differences[int(0.025 * (len(boot_differences) - 1))]
    hi = boot_differences[int(0.975 * (len(boot_differences) - 1))]

    pooled = base + adapt
    n_adapt = len(adapt)
    extreme = 0

    for _ in range(N_PERMUTATIONS):
        permuted = pooled[:]
        rng.shuffle(permuted)
        perm_a = permuted[:n_adapt]
        perm_b = permuted[n_adapt:]
        diff = sum(perm_a) / len(perm_a) - sum(perm_b) / len(perm_b)
        if abs(diff) >= abs(observed):
            extreme += 1

    p = (extreme + 1) / (N_PERMUTATIONS + 1)

    return {
        "difference": observed,
        "ci_low": lo,
        "ci_high": hi,
        "p_two_sided": p,
        "status": "cluster-level exploratory comparison",
    }


def main():
    target_sequence, _ = load_target_reference(TARGET)

    all_rows = []
    source_inventory = []
    seen_structures = set()
    source_count = 0
    parsed_model_count = 0
    duplicate_count = 0

    for source in iter_arm_labeled_pdb_sources():
        source_count += 1
        parts = source["parts"]
        study_group = study_group_from_parts(parts)
        cluster_id = cluster_from_parts(parts)
        models = parse_models(source["text"])
        parsed_model_count += len(models)

        source_valid = 0
        source_duplicate = 0

        for model_number, model_lines in enumerate(models, start=1):
            # Deduplicate identical coordinate models within the same study
            # and arm, but do not merge identical structures across arms.
            signature = hashlib.sha256(
                "\n".join(model_lines).encode("utf-8")
            ).hexdigest()
            duplicate_key = (
                signature, source["arm"], study_group, cluster_id
            )

            if duplicate_key in seen_structures:
                duplicate_count += 1
                source_duplicate += 1
                continue
            seen_structures.add(duplicate_key)

            metrics = score_one_model(model_lines, target_sequence)

            def fmt_num(value):
                if value is None:
                    return ""
                if isinstance(value, float) and not math.isfinite(value):
                    return ""
                return value

            row = {
                "arm": source["arm"],
                "study_group": study_group,
                "cluster_id": cluster_id,
                "model_number": model_number,
                "source_file": source["source_file"],
                "peptide_sequence": metrics.get("peptide_sequence", ""),
                "peptide_length": metrics.get("peptide_length", ""),
                "receptor_chain": metrics.get("receptor_chain", ""),
                "peptide_chain": metrics.get("peptide_chain", ""),
                "target_local_sequence_match_fraction": fmt_num(
                    metrics.get("target_match_fraction")
                ),
                "min_interchain_heavy_atom_distance_A": fmt_num(
                    metrics.get("min_interchain")
                ),
                "peptide_to_ARG203_min_distance_A": fmt_num(
                    metrics.get("d203")
                ),
                "peptide_to_ARG205_min_distance_A": fmt_num(
                    metrics.get("d205")
                ),
                "eligible_for_comparison": int(
                    metrics.get("eligible", 0)
                ),
                "passes_all_geometric_criteria": metrics.get(
                    "pass_geometry", ""
                ),
                "status": metrics.get("status", "unknown"),
            }
            all_rows.append(row)

            if row["eligible_for_comparison"] == 1:
                source_valid += 1

        source_inventory.append({
            "source_file": source["source_file"],
            "arm": source["arm"],
            "study_group": study_group,
            "cluster_id": cluster_id,
            "pdb_models_parsed": len(models),
            "eligible_models_scored": source_valid,
            "duplicate_models_skipped": source_duplicate,
        })

    write_csv(
        OUT / "models_scored.csv",
        all_rows,
        MODEL_FIELDS,
    )
    write_csv(
        OUT / "source_inventory.csv",
        source_inventory,
        [
            "source_file", "arm", "study_group", "cluster_id",
            "pdb_models_parsed", "eligible_models_scored",
            "duplicate_models_skipped",
        ],
    )

    eligible_rows = [
        r for r in all_rows if r["eligible_for_comparison"] == 1
    ]
    study_groups = sorted({r["study_group"] for r in eligible_rows})

    cluster_rows = []
    for group in study_groups:
        for arm in ("baseline", "adaptive"):
            subset = [
                r for r in eligible_rows
                if r["study_group"] == group and r["arm"] == arm
                and r["cluster_id"]
            ]
            by_cluster = defaultdict(list)
            for row in subset:
                by_cluster[row["cluster_id"]].append(
                    int(row["passes_all_geometric_criteria"])
                )

            for cluster, values in sorted(by_cluster.items()):
                cluster_rows.append({
                    "study_group": group,
                    "arm": arm,
                    "cluster_id": cluster,
                    "eligible_models": len(values),
                    "passing_models": sum(values),
                    "cluster_pass_rate": sum(values) / len(values),
                })

    write_csv(
        OUT / "cluster_summary.csv",
        cluster_rows,
        [
            "study_group", "arm", "cluster_id",
            "eligible_models", "passing_models", "cluster_pass_rate",
        ],
    )

    group_summaries = []
    cluster_rates_cache = {}

    for group in study_groups:
        arm_records = {
            arm: [
                r for r in eligible_rows
                if r["study_group"] == group and r["arm"] == arm
            ]
            for arm in ("baseline", "adaptive")
        }

        # Only a within-study comparison with both arms is interpretable.
        if not arm_records["baseline"] or not arm_records["adaptive"]:
            continue

        base = arm_records["baseline"]
        adapt = arm_records["adaptive"]

        base_pass = sum(int(r["passes_all_geometric_criteria"]) for r in base)
        adapt_pass = sum(int(r["passes_all_geometric_criteria"]) for r in adapt)

        br = cluster_rates_for(base)
        ar = cluster_rates_for(adapt)
        cluster_rates_cache[(group, "baseline")] = br
        cluster_rates_cache[(group, "adaptive")] = ar

        inference = cluster_comparison(br, ar)
        group_summaries.append({
            "study_group": group,
            "baseline_models": len(base),
            "baseline_passes": base_pass,
            "baseline_model_pass_rate": base_pass / len(base),
            "baseline_identified_clusters": len(br),
            "adaptive_models": len(adapt),
            "adaptive_passes": adapt_pass,
            "adaptive_model_pass_rate": adapt_pass / len(adapt),
            "adaptive_identified_clusters": len(ar),
            "adaptive_minus_baseline_cluster_rate_difference":
                inference["difference"],
            "cluster_bootstrap_95ci_low": inference["ci_low"],
            "cluster_bootstrap_95ci_high": inference["ci_high"],
            "cluster_permutation_p_two_sided": inference["p_two_sided"],
            "inference_status": inference["status"],
        })

    write_csv(
        OUT / "study_group_summary.csv",
        group_summaries,
        [
            "study_group",
            "baseline_models", "baseline_passes", "baseline_model_pass_rate",
            "baseline_identified_clusters",
            "adaptive_models", "adaptive_passes", "adaptive_model_pass_rate",
            "adaptive_identified_clusters",
            "adaptive_minus_baseline_cluster_rate_difference",
            "cluster_bootstrap_95ci_low", "cluster_bootstrap_95ci_high",
            "cluster_permutation_p_two_sided", "inference_status",
        ],
    )

    # The report explicitly distinguishes structural geometry from measured affinity.
    report = [
        "# Existing-candidate structural translation test",
        "",
        f"**Generated:** {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"**Target reference:** `{TARGET}`",
        f"**Target reference sequence loaded:** {'yes' if target_sequence else 'no'}",
        "",
        "## Question",
        "",
        "Does the adaptive arm's advantage remain visible when existing candidate "
        "structures are evaluated with the same target-aware geometric criteria?",
        "",
        "This is a structural-proxy test. It is not an experimental measurement of "
        "binding affinity and cannot by itself prove stronger real-world binding.",
        "",
        "## Predefined criteria",
        "",
        f"1. Minimum receptor–peptide heavy-atom distance >= {MIN_INTERCHAIN_A:.1f} Å.",
        f"2. Minimum peptide-to-ARG203 heavy-atom distance <= {MAX_ARG203_A:.1f} Å.",
        f"3. Minimum peptide-to-ARG205 heavy-atom distance <= {MAX_ARG205_A:.1f} Å.",
        "",
        "A model passes only when it meets all three criteria. Model-level pass rates "
        "are reported descriptively; where cluster labels can be identified, the "
        "exploratory comparison uses clusters as the resampling/permutation units.",
        "",
        "## Scan summary",
        "",
        f"- Arm-labeled PDB files/archive members discovered: {source_count}.",
        f"- PDB models parsed: {parsed_model_count}.",
        f"- Identical model duplicates skipped within the same arm/study/cluster: {duplicate_count}.",
        f"- Models eligible after receptor-sequence verification: {len(eligible_rows)}.",
        f"- Model rows recorded, including ineligible structures: {len(all_rows)}.",
        "",
        "## Comparable study groups",
        "",
    ]

    if not group_summaries:
        report += [
            "**No study group had eligible structures from both arms.** Therefore, "
            "this scan did not produce a valid adaptive-versus-baseline structural "
            "comparison. This is not evidence for either arm.",
            "",
            "Check `source_inventory.csv` and `models_scored.csv` for arm labels, "
            "cluster paths, ARG203/ARG205 availability, receptor sequence matches, "
            "and peptide-chain selection. If the structures are not present or "
            "cannot be linked to the actual candidate sets, a controlled structural "
            "evaluation would need to be carried out before drawing a conclusion.",
            "",
        ]
    else:
        report += [
            "| Existing study group | Baseline pass rate | Adaptive pass rate | Adaptive minus baseline cluster rate | 95% cluster-bootstrap CI | Cluster permutation p |",
            "|---|---:|---:|---:|---:|---:|",
        ]

        for s in group_summaries:
            diff = s["adaptive_minus_baseline_cluster_rate_difference"]
            lo = s["cluster_bootstrap_95ci_low"]
            hi = s["cluster_bootstrap_95ci_high"]
            p = s["cluster_permutation_p_two_sided"]

            diff_s = f"{diff:.3f}" if isinstance(diff, (int, float)) else "not estimated"
            ci_s = (
                f"[{lo:.3f}, {hi:.3f}]"
                if isinstance(lo, (int, float)) and isinstance(hi, (int, float))
                else "not estimated"
            )
            p_s = f"{p:.4f}" if isinstance(p, (int, float)) else "not estimated"

            report.append(
                f"| `{s['study_group']}` | "
                f"{s['baseline_model_pass_rate']:.1%} "
                f"({s['baseline_passes']}/{s['baseline_models']}) | "
                f"{s['adaptive_model_pass_rate']:.1%} "
                f"({s['adaptive_passes']}/{s['adaptive_models']}) | "
                f"{diff_s} | {ci_s} | {p_s} |"
            )

        report += [
            "",
            "## Interpretation limits",
            "",
            "- A higher geometric pass rate would support greater structural plausibility "
            "under these predefined criteria, not prove greater experimental affinity.",
            "- A similar or lower pass rate would indicate that the internal-score "
            "advantage was not mirrored by this structural proxy.",
            "- Docking models from the same cluster are related observations; model-level "
            "percentages must not be treated as independent biological replicates.",
            "- The source path and peptide sequence should be checked against the optimizer's "
            "actual candidate inventory before claiming that the tested poses represent "
            "the relevant candidate set.",
            "- The existing P3 paired optimizer run was previously reported to contain no "
            "candidate PDBs. If these structures come from a different P2 run, label the "
            "result as validation of that P2 structural set, not as a direct scoring of all "
            "480 current P3 records per arm.",
            "",
        ]

    report += [
        "## Files",
        "",
        "- `models_scored.csv`: per-model geometry measurements and eligibility.",
        "- `study_group_summary.csv`: same-study adaptive-versus-baseline comparisons.",
        "- `cluster_summary.csv`: per-cluster pass rates.",
        "- `source_inventory.csv`: exact source paths and scan coverage.",
        "- `structural_translation_pass_rates.png`: comparison graph, when comparable groups exist.",
        "",
    ]

    (OUT / "REPORT.md").write_text("\n".join(report), encoding="utf-8")

    # Separate logbook notes that can be incorporated after verifying the source
    # models correspond to the optimizer candidate sets.
    (OUT / "DATA_RESULTS_LOGBOOK.md").write_text(
        "# Logbook Supplement — Existing Candidate Structural Test\n\n"
        f"**Date:** {datetime.now().date().isoformat()}\n\n"
        "The analysis reused existing PDB structures and structural archives. "
        "No optimizer run, new seed comparison, or new docking job was launched "
        "by this script.\n\n"
        f"The test applied the predefined interchain-distance criterion "
        f"(minimum >= {MIN_INTERCHAIN_A:.1f} Å) and the ARG203 and ARG205 contact "
        f"criteria (each minimum <= {MAX_ARG203_A:.1f} Å and "
        f"{MAX_ARG205_A:.1f} Å). A model passed only if all three conditions held.\n\n"
        f"Arm-labeled sources discovered: {source_count}. PDB models parsed: "
        f"{parsed_model_count}. Models eligible after receptor-sequence verification: "
        f"{len(eligible_rows)}. Comparable study groups: {len(group_summaries)}.\n\n"
        "See `REPORT.md`, `models_scored.csv`, and `study_group_summary.csv` for exact "
        "measurements. These are geometric structural proxies, not experimental binding "
        "affinity measurements. Verify each model's peptide sequence and source lineage "
        "before attributing the result to a specific optimizer run.\n",
        encoding="utf-8",
    )

    if group_summaries:
        analysis_text = [
            "# Logbook Supplement — Structural Test Analysis",
            "",
            "This test asks whether adaptive-versus-baseline differences persist under "
            "a target-aware geometric evaluation of existing candidate structures.",
            "",
        ]
        for s in group_summaries:
            br = s["baseline_model_pass_rate"]
            ar = s["adaptive_model_pass_rate"]
            if ar > br:
                direction = "higher"
            elif ar < br:
                direction = "lower"
            else:
                direction = "equal"

            analysis_text.append(
                f"For `{s['study_group']}`, the adaptive model pass rate was "
                f"{ar:.1%} and the baseline pass rate was {br:.1%}; adaptive was "
                f"{direction} on this structural criterion."
            )

            if isinstance(
                s["adaptive_minus_baseline_cluster_rate_difference"], (int, float)
            ):
                d = s["adaptive_minus_baseline_cluster_rate_difference"]
                lo = s["cluster_bootstrap_95ci_low"]
                hi = s["cluster_bootstrap_95ci_high"]
                p = s["cluster_permutation_p_two_sided"]
                analysis_text.append(
                    f"The adaptive-minus-baseline difference in mean cluster pass rate "
                    f"was {d:.3f}, 95% cluster-bootstrap interval "
                    f"[{lo:.3f}, {hi:.3f}], exploratory two-sided cluster permutation "
                    f"p={p:.4f}."
                )
            else:
                analysis_text.append(
                    "Cluster-level inference was not estimated because the identified "
                    "cluster count was too small or cluster labels were unavailable."
                )

        analysis_text += [
            "",
            "These results can indicate whether an internal-score advantage is accompanied "
            "by a difference in structural plausibility. They do not establish stronger "
            "experimental binding. A direct real-world conclusion requires comparable "
            "experimental binding-affinity measurements for candidates from both arms. "
            "Candidate sequence and run lineage must also be verified before attributing "
            "the result to the current optimizer comparison.",
        ]
    else:
        analysis_text = [
            "# Logbook Supplement — Structural Test Analysis",
            "",
            "The existing structural scan did not find a study group with eligible "
            "structures from both arms. This test therefore produced no comparative "
            "structural evidence and does not favor either method.",
            "",
            "The next valid step is to locate or generate comparable target-aware structures "
            "for both arms using the same receptor and scoring procedure. Even then, structural "
            "plausibility would remain a computational proxy rather than proof of experimental "
            "binding.",
        ]

    (OUT / "ANALYSIS_LOGBOOK.md").write_text(
        "\n".join(analysis_text) + "\n", encoding="utf-8"
    )

    # Plot only comparable groups.
    if group_summaries:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            labels = [
                s["study_group"].split("/")[-2:]
                for s in group_summaries
            ]
            labels = ["/".join(x) for x in labels]

            baseline_rates = [
                s["baseline_model_pass_rate"] for s in group_summaries
            ]
            adaptive_rates = [
                s["adaptive_model_pass_rate"] for s in group_summaries
            ]

            x = np.arange(len(group_summaries))
            width = 0.36
            fig, ax = plt.subplots(figsize=(max(8, 1.8 * len(x)), 5.5))

            bars1 = ax.bar(
                x - width / 2, baseline_rates, width,
                label="Baseline",
            )
            bars2 = ax.bar(
                x + width / 2, adaptive_rates, width,
                label="Adaptive",
            )

            ax.set_ylabel("Models passing all geometric criteria")
            ax.set_ylim(0, 1.12)
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=25, ha="right")
            ax.set_title("Existing target-aware structural validation")
            ax.legend()
            ax.bar_label(bars1, fmt="%.0%%", padding=2)
            ax.bar_label(bars2, fmt="%.0%%", padding=2)

            fig.text(
                0.01, 0.01,
                "Structural plausibility proxy only; not measured binding affinity. "
                "Verify source run and peptide sequence lineage.",
                fontsize=8,
            )
            fig.tight_layout(rect=[0, 0.06, 1, 1])

            graph = OUT / "structural_translation_pass_rates.png"
            fig.savefig(graph, dpi=200, bbox_inches="tight")
            plt.close(fig)

            if PACKAGE is not None and PACKAGE.exists():
                graphs_dir = PACKAGE / "graphs"
                poster_dir = PACKAGE / "poster_graphs"
                graphs_dir.mkdir(exist_ok=True)
                poster_dir.mkdir(exist_ok=True)
                shutil.copy2(
                    graph,
                    graphs_dir / "10_target_aware_structural_translation_proxy.png",
                )
                shutil.copy2(
                    graph,
                    poster_dir / "04_target_aware_structural_translation_proxy.png",
                )

        except Exception as exc:
            (OUT / "GRAPH_GENERATION_ERROR.txt").write_text(
                repr(exc), encoding="utf-8"
            )

    print("\nEXISTING-CANDIDATE STRUCTURAL TEST FINISHED")
    print("Output folder:", OUT)
    print("PDB models parsed:", parsed_model_count)
    print("Models eligible after receptor sequence verification:", len(eligible_rows))
    print("Comparable study groups:", len(group_summaries))
    print("Model results:", OUT / "models_scored.csv")
    print("Group comparison:", OUT / "study_group_summary.csv")
    print("Report:", OUT / "REPORT.md")
    print("Data logbook supplement:", OUT / "DATA_RESULTS_LOGBOOK.md")
    print("Analysis logbook supplement:", OUT / "ANALYSIS_LOGBOOK.md")
    if group_summaries:
        print("Graph:", OUT / "structural_translation_pass_rates.png")
    else:
        print("No valid same-study structural comparison was found.")
    print("\nNo optimizer run, seed comparison, or docking job was launched.")


if __name__ == "__main__":
    main()
