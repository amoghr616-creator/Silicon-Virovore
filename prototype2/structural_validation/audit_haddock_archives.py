#!/usr/bin/env python3
"""Audit HADDOCK PDB summary archives: 1YCR native redocking RMSD and 9SYA geometry."""
import argparse, csv, math, re, shutil, statistics, tarfile
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.Align import PairwiseAligner

AA_CUSTOM = {"MSE":"M", "HID":"H", "HIE":"H", "HIP":"H", "UNK":"X"}
PARSER = PDBParser(QUIET=True)
ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.match_score = 2.0
ALIGNER.mismatch_score = -1.0
ALIGNER.open_gap_score = -8.0
ALIGNER.extend_gap_score = -0.5


def get_chain(structure, chain_id, path):
    # Accept either a Bio.PDB Structure or a Model already selected via structure[0].
    # A Model also exposes get_models(), so hasattr(..., "get_models") is not a safe type test.
    if isinstance(structure, Structure):
        models = list(structure.get_models())
        if not models:
            raise ValueError(f"{path}: structure has no models")
        model = models[0]
    else:
        model = structure
    if chain_id not in model:
        raise ValueError(f"{path}: chain {chain_id!r} missing; available chains: {list(model.child_dict)}")
    return model[chain_id]


def aa_residues(chain):
    return [r for r in chain.get_residues() if is_aa(r, standard=False) and r.has_id("CA")]


def sequence(residues):
    return "".join(seq1(r.resname, custom_map=AA_CUSTOM, undef_code="X") for r in residues)


def align_residues(native_residues, model_residues):
    s_native, s_model = sequence(native_residues), sequence(model_residues)
    if not s_native or not s_model:
        raise ValueError("Cannot align an empty protein/peptide chain")
    aln = ALIGNER.align(s_native, s_model)[0]
    pairs, model_to_native = [], {}
    for ni, mi in zip(aln.indices[0], aln.indices[1]):
        ni, mi = int(ni), int(mi)
        if ni >= 0 and mi >= 0 and s_native[ni] == s_model[mi]:
            pairs.append((ni, mi))
            model_to_native[mi] = ni
    return pairs, model_to_native, s_native, s_model


def heavy_atom_cloud(residues, local_to_mapped=None):
    coords, residue_ids = [], []
    for i, residue in enumerate(residues):
        if local_to_mapped is not None and i not in local_to_mapped:
            continue
        mapped_id = local_to_mapped[i] if local_to_mapped is not None else i
        for atom in residue.get_atoms():
            element = (atom.element or re.sub(r"^[0-9]+", "", atom.get_name())[:1]).strip().upper()
            if element in {"H", "D"}:
                continue
            coords.append(atom.coord.astype(float))
            residue_ids.append(mapped_id)
    if not coords:
        return np.empty((0, 3), dtype=float), np.array([], dtype=int)
    return np.asarray(coords, dtype=float), np.asarray(residue_ids, dtype=int)


def pairwise_distances(a, b):
    if len(a) == 0 or len(b) == 0:
        return np.empty((len(a), len(b)), dtype=float)
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)


def residue_contact_set(distances, residue_ids_a, residue_ids_b, cutoff=5.0):
    ai, bi = np.where(distances <= cutoff)
    return {(int(residue_ids_a[i]), int(residue_ids_b[j])) for i, j in zip(ai, bi)}


def residue_heavy_coords(residue):
    coords = []
    for atom in residue.get_atoms():
        element = (atom.element or re.sub(r"^[0-9]+", "", atom.get_name())[:1]).strip().upper()
        if element not in {"H", "D"}:
            coords.append(atom.coord.astype(float))
    return np.asarray(coords, dtype=float) if coords else np.empty((0, 3), dtype=float)


def target_residue_contact(chain, target_number, peptide_residues):
    target = next((r for r in chain.get_residues()
                   if r.id[1] == target_number and r.resname.upper() == "ARG"), None)
    if target is None:
        return None, None
    target_xyz = residue_heavy_coords(target)
    peptide_xyz, _ = heavy_atom_cloud(peptide_residues)
    if len(target_xyz) == 0 or len(peptide_xyz) == 0:
        return None, None
    d = pairwise_distances(target_xyz, peptide_xyz)
    return float(d.min()), int(np.count_nonzero(d <= 5.0))


def round_or_blank(value, digits=4):
    if value is None:
        return ""
    try:
        if not math.isfinite(float(value)):
            return ""
        return round(float(value), digits)
    except (TypeError, ValueError):
        return value


def extract_pdbs(tar_path, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:*") as tf:
        for member in tf.getmembers():
            if not member.isfile() or not member.name.lower().endswith(".pdb"):
                continue
            basename = Path(member.name).name
            if not basename:
                continue
            source = tf.extractfile(member)
            if source is None:
                continue
            with source, (out_dir / basename).open("wb") as dest:
                shutil.copyfileobj(source, dest)
    models = sorted(out_dir.glob("*.pdb"))
    if not models:
        raise ValueError(f"No PDB files extracted from {tar_path}")
    return models


def cluster_label(filename):
    m = re.search(r"cluster(\d+)", filename, flags=re.I)
    return m.group(1) if m else "unknown"


def analyze_control(model_path, native_model):
    model = PARSER.get_structure(model_path.stem, str(model_path))[0]
    nrec_chain, npep_chain = native_model["A"], native_model["B"]
    mrec_chain, mpep_chain = get_chain(model, "A", model_path), get_chain(model, "B", model_path)
    nrec, npep = aa_residues(nrec_chain), aa_residues(npep_chain)
    mrec, mpep = aa_residues(mrec_chain), aa_residues(mpep_chain)
    rec_pairs, rec_map, nrec_seq, mrec_seq = align_residues(nrec, mrec)
    pep_pairs, pep_map, npep_seq, mpep_seq = align_residues(npep, mpep)

    rec_fit = [(i, j) for i, j in rec_pairs if nrec[i].has_id("CA") and mrec[j].has_id("CA")]
    if len(rec_fit) < 3:
        raise ValueError(f"{model_path.name}: only {len(rec_fit)} matched receptor CA atoms")
    fixed_rec = np.asarray([nrec[i]["CA"].coord for i, j in rec_fit], dtype=float)
    moving_rec = np.asarray([mrec[j]["CA"].coord for i, j in rec_fit], dtype=float)
    sup_rec = SVDSuperimposer(); sup_rec.set(fixed_rec, moving_rec); sup_rec.run()
    rot, tran = sup_rec.get_rotran()
    rec_ca_rmsd = float(sup_rec.get_rms())

    fixed_bb, moving_bb = [], []
    for ni, mi in pep_pairs:
        for atom_name in ("N", "CA", "C"):
            if npep[ni].has_id(atom_name) and mpep[mi].has_id(atom_name):
                fixed_bb.append(npep[ni][atom_name].coord.astype(float))
                moving_bb.append(mpep[mi][atom_name].coord.astype(float))
    fixed_bb, moving_bb = np.asarray(fixed_bb), np.asarray(moving_bb)
    if len(fixed_bb) < 3:
        frame_rmsd = fitted_rmsd = None
    else:
        moving_transformed = moving_bb @ rot + tran
        frame_rmsd = float(np.sqrt(np.mean(np.sum((fixed_bb - moving_transformed) ** 2, axis=1))))
        sup_pep = SVDSuperimposer(); sup_pep.set(fixed_bb, moving_bb); sup_pep.run()
        fitted_rmsd = float(sup_pep.get_rms())

    nrxyz, nrids = heavy_atom_cloud(nrec, {i: i for i in range(len(nrec))})
    npxyz, npids = heavy_atom_cloud(npep, {i: i for i in range(len(npep))})
    native_d = pairwise_distances(nrxyz, npxyz)
    native_contacts = residue_contact_set(native_d, nrids, npids, 5.0)
    mrxyz, mrids = heavy_atom_cloud(mrec, rec_map)
    mpxyz, mpids = heavy_atom_cloud(mpep, pep_map)
    model_d_mapped = pairwise_distances(mrxyz, mpxyz)
    model_contacts = residue_contact_set(model_d_mapped, mrids, mpids, 5.0)
    intersection = len(native_contacts & model_contacts)
    precision = intersection / len(model_contacts) if model_contacts else None
    recall = intersection / len(native_contacts) if native_contacts else None
    f1 = (2 * precision * recall / (precision + recall)) if precision is not None and recall is not None and precision + recall else None

    # Geometry screen uses all heavy atoms, not just residues successfully sequence-mapped.
    all_rxyz, all_rids = heavy_atom_cloud(mrec)
    all_pxyz, all_pids = heavy_atom_cloud(mpep)
    geom_d = pairwise_distances(all_rxyz, all_pxyz)
    contacts_5 = residue_contact_set(geom_d, all_rids, all_pids, 5.0)
    return {
        "group": "control_1YCR", "model": model_path.name, "cluster": cluster_label(model_path.name),
        "receptor_sequence_native": nrec_seq, "receptor_sequence_model": mrec_seq,
        "peptide_sequence_native": npep_seq, "peptide_sequence_model": mpep_seq,
        "receptor_CA_fit_pairs": len(rec_fit), "peptide_backbone_atom_pairs": len(fixed_bb),
        "receptor_CA_RMSD_A": rec_ca_rmsd,
        "peptide_backbone_RMSD_receptor_aligned_A": frame_rmsd,
        "peptide_backbone_RMSD_peptide_fitted_A": fitted_rmsd,
        "native_residue_contacts_5A": len(native_contacts), "model_residue_contacts_5A": len(contacts_5),
        "native_contact_precision": precision, "native_contact_recall": recall, "native_contact_F1": f1,
        "interchain_min_distance_A": float(geom_d.min()) if geom_d.size else None,
        "atom_pairs_lt2A": int(np.count_nonzero(geom_d < 2.0)),
        "atom_pairs_lt2_5A": int(np.count_nonzero(geom_d < 2.5)),
        "atom_pairs_le5A": int(np.count_nonzero(geom_d <= 5.0)),
        "residue_pairs_le5A": len(contacts_5),
        "Arg203_min_A": None, "Arg203_atom_pairs_le5A": None,
        "Arg205_min_A": None, "Arg205_atom_pairs_le5A": None,
        "geometry_screen_pass": None,
    }


def analyze_candidate(model_path, group):
    model = PARSER.get_structure(model_path.stem, str(model_path))[0]
    rec_chain, pep_chain = get_chain(model, "A", model_path), get_chain(model, "B", model_path)
    rec, pep = aa_residues(rec_chain), aa_residues(pep_chain)
    rxyz, rids = heavy_atom_cloud(rec)
    pxyz, pids = heavy_atom_cloud(pep)
    d = pairwise_distances(rxyz, pxyz)
    contacts = residue_contact_set(d, rids, pids, 5.0)
    arg203_min, arg203_pairs = target_residue_contact(rec_chain, 203, pep)
    arg205_min, arg205_pairs = target_residue_contact(rec_chain, 205, pep)
    below2 = int(np.count_nonzero(d < 2.0))
    below25 = int(np.count_nonzero(d < 2.5))
    passes = bool(below2 == 0 and (arg203_pairs or 0) > 0 and (arg205_pairs or 0) > 0)
    return {
        "group": group, "model": model_path.name, "cluster": cluster_label(model_path.name),
        "receptor_sequence_native": "", "receptor_sequence_model": sequence(rec),
        "peptide_sequence_native": "", "peptide_sequence_model": sequence(pep),
        "receptor_CA_fit_pairs": None, "peptide_backbone_atom_pairs": None,
        "receptor_CA_RMSD_A": None, "peptide_backbone_RMSD_receptor_aligned_A": None,
        "peptide_backbone_RMSD_peptide_fitted_A": None,
        "native_residue_contacts_5A": None, "model_residue_contacts_5A": len(contacts),
        "native_contact_precision": None, "native_contact_recall": None, "native_contact_F1": None,
        "interchain_min_distance_A": float(d.min()) if d.size else None,
        "atom_pairs_lt2A": below2, "atom_pairs_lt2_5A": below25,
        "atom_pairs_le5A": int(np.count_nonzero(d <= 5.0)), "residue_pairs_le5A": len(contacts),
        "Arg203_min_A": arg203_min, "Arg203_atom_pairs_le5A": arg203_pairs,
        "Arg205_min_A": arg205_min, "Arg205_atom_pairs_le5A": arg205_pairs,
        "geometry_screen_pass": passes,
    }


def median_field(rows, field):
    vals = [float(r[field]) for r in rows if r.get(field) not in (None, "") and math.isfinite(float(r[field]))]
    return round(statistics.median(vals), 4) if vals else ""


def write_csv(path, rows):
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: round_or_blank(v, 6) if isinstance(v, (float, np.floating)) else v for k, v in row.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", required=True, help="Full native 1YCR.pdb containing receptor chain A and p53 chain B")
    ap.add_argument("--control-tar", required=True)
    ap.add_argument("--adaptive-tar", required=True)
    ap.add_argument("--baseline-tar", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    for p in [args.native, args.control_tar, args.adaptive_tar, args.baseline_tar]:
        if not Path(p).is_file():
            raise SystemExit(f"ERROR: file not found: {p}")
    out = Path(args.out).expanduser().resolve(); out.mkdir(parents=True, exist_ok=True)
    control_dir, adaptive_dir, baseline_dir = out / "control_models", out / "adaptive_models", out / "baseline_models"
    control_models = extract_pdbs(args.control_tar, control_dir)
    adaptive_models = extract_pdbs(args.adaptive_tar, adaptive_dir)
    baseline_models = extract_pdbs(args.baseline_tar, baseline_dir)
    native_struct = PARSER.get_structure("native_1YCR", str(Path(args.native).expanduser()))[0]
    if "A" not in native_struct or "B" not in native_struct:
        raise SystemExit(f"ERROR: native 1YCR must contain chain A and chain B; found {list(native_struct.child_dict)}")

    rows = []
    for path in control_models:
        rows.append(analyze_control(path, native_struct))
    for group, paths in [("adaptive_9SYA", adaptive_models), ("baseline_9SYA", baseline_models)]:
        for path in paths:
            rows.append(analyze_candidate(path, group))

    write_csv(out / "model_metrics.csv", rows)
    summaries = []
    for group in ["control_1YCR", "adaptive_9SYA", "baseline_9SYA"]:
        rr = [r for r in rows if r["group"] == group]
        summary = {
            "group": group, "n_models": len(rr),
            "median_interchain_min_distance_A": median_field(rr, "interchain_min_distance_A"),
            "median_atom_pairs_le5A": median_field(rr, "atom_pairs_le5A"),
            "median_residue_pairs_le5A": median_field(rr, "residue_pairs_le5A"),
            "models_with_any_pair_lt2A": sum((r["atom_pairs_lt2A"] or 0) > 0 for r in rr),
            "models_with_any_pair_lt2_5A": sum((r["atom_pairs_lt2_5A"] or 0) > 0 for r in rr),
            "median_peptide_backbone_RMSD_receptor_aligned_A": median_field(rr, "peptide_backbone_RMSD_receptor_aligned_A"),
            "n_control_RMSD_le2A": sum(r["peptide_backbone_RMSD_receptor_aligned_A"] is not None and r["peptide_backbone_RMSD_receptor_aligned_A"] <= 2.0 for r in rr),
            "median_native_contact_recall": median_field(rr, "native_contact_recall"),
            "median_native_contact_precision": median_field(rr, "native_contact_precision"),
            "models_passing_preliminary_9SYA_geometry_screen": sum(r["geometry_screen_pass"] is True for r in rr),
        }
        summaries.append(summary)
    write_csv(out / "summary.csv", summaries)

    cluster_summaries = []
    for group in ["control_1YCR", "adaptive_9SYA", "baseline_9SYA"]:
        groups = sorted({r["cluster"] for r in rows if r["group"] == group})
        for cluster in groups:
            cr = [r for r in rows if r["group"] == group and r["cluster"] == cluster]
            cluster_summaries.append({
                "group": group, "cluster": cluster, "n_models": len(cr),
                "median_residue_pairs_le5A": median_field(cr, "residue_pairs_le5A"),
                "median_atom_pairs_le5A": median_field(cr, "atom_pairs_le5A"),
                "median_receptor_aligned_peptide_RMSD_A": median_field(cr, "peptide_backbone_RMSD_receptor_aligned_A"),
                "median_native_contact_recall": median_field(cr, "native_contact_recall"),
                "geometry_pass_count": sum(r["geometry_screen_pass"] is True for r in cr),
            })
    write_csv(out / "cluster_summary.csv", cluster_summaries)

    print("\n=== HADDOCK STRUCTURAL AUDIT ===")
    print(f"Native reference: {Path(args.native).resolve()}")
    print(f"Models: control={len(control_models)}, adaptive={len(adaptive_models)}, baseline={len(baseline_models)}")
    for s in summaries:
        print(f"\n{s['group']} (n={s['n_models']})")
        print(f"  median minimum inter-chain distance: {s['median_interchain_min_distance_A']} Å")
        print(f"  median atom pairs <=5 Å: {s['median_atom_pairs_le5A']}")
        print(f"  median residue pairs <=5 Å: {s['median_residue_pairs_le5A']}")
        print(f"  models with any pair <2.0 Å: {s['models_with_any_pair_lt2A']}/{s['n_models']}")
        print(f"  models with any pair <2.5 Å: {s['models_with_any_pair_lt2_5A']}/{s['n_models']}")
        if s['group'] == 'control_1YCR':
            print(f"  median receptor-aligned p53 backbone RMSD: {s['median_peptide_backbone_RMSD_receptor_aligned_A']} Å")
            print(f"  models with receptor-aligned backbone RMSD <=2.0 Å: {s['n_control_RMSD_le2A']}/{s['n_models']}")
            print(f"  median native-contact recall / precision: {s['median_native_contact_recall']} / {s['median_native_contact_precision']}")
        else:
            print(f"  preliminary geometry-screen passes: {s['models_passing_preliminary_9SYA_geometry_screen']}/{s['n_models']}")
    print(f"\nWrote: {out / 'model_metrics.csv'}")
    print(f"Wrote: {out / 'summary.csv'}")
    print(f"Wrote: {out / 'cluster_summary.csv'}")
    print("Important: the 9SYA screen is geometric only; the existing adaptive/baseline HADDOCK jobs used different peptide active-residue restraints.")

if __name__ == "__main__":
    main()
