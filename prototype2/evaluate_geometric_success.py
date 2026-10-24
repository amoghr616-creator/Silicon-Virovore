#!/usr/bin/env python3
"""Evaluate the frozen Vina pose panel using hotspot contact + atom-type vdW clashes.

This is a geometric screening endpoint, not experimental binding evidence.
Missing or failed planned runs remain failures in all denominators.
"""
from __future__ import annotations
import csv, hashlib, json, math, statistics
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SV = ROOT / "prototype2" / "structural_validation"
RECEPTOR = ROOT / "prototype2" / "target_aware_audit" / "9SYA_chainA_repaired.pdb"
HOT_ATOMS = {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}
# Approximate Bondi-style van der Waals radii (angstroms) for a transparent first-pass clash screen.
VDW = {"C":1.70,"N":1.55,"O":1.52,"S":1.80,"P":1.80,"F":1.47,"CL":1.75,"BR":1.85,"I":1.98,"SE":1.90}
TYPE_TO_ELEMENT = {"A":"C","C":"C","N":"N","NA":"N","NS":"N","O":"O","OA":"O","OS":"O","S":"S","SA":"S","P":"P","F":"F","CL":"CL","BR":"BR","I":"I","H":"H","HD":"H","HS":"H"}

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()

def read_receptor(path: Path):
    all_atoms=[]; hot={203:[],205:[]}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ","HETATM")): continue
        name=line[12:16].strip(); chain=line[21:22].strip(); resname=line[17:20].strip(); resseq=int(line[22:26])
        element=(line[76:78].strip().upper() or name[:1].upper())
        xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
        if not all(math.isfinite(v) for v in xyz): raise ValueError(f"Non-finite receptor coordinate: {line}")
        if element in {"H","D"}: continue
        a={"name":name,"chain":chain,"resname":resname,"resseq":resseq,"element":element,"xyz":xyz}
        all_atoms.append(a)
        if chain=="A" and resseq in hot and resname=="ARG" and name in HOT_ATOMS:
            hot[resseq].append(a)
    if not all_atoms or not hot[203] or not hot[205]: raise ValueError("Receptor has no atoms or chain A Arg203/Arg205 mapping is invalid")
    for res in (203,205):
        present={a["name"] for a in hot[res]}
        missing=HOT_ATOMS-present
        if missing: raise ValueError(f"Arg A{res} side-chain atoms missing: {sorted(missing)}")
    return all_atoms,hot

def parse_pdbqt_models(path: Path):
    text=path.read_text(errors="replace")
    models=[]; current=[]; score=None; in_model=False; model_seen=False
    for line in text.splitlines():
        if line.startswith("MODEL"):
            model_seen=True
            if current: models.append((current,score))
            current=[]; score=None; in_model=True
        elif line.startswith("ENDMDL"):
            if current: models.append((current,score))
            current=[]; score=None; in_model=False
        elif line.startswith("REMARK VINA RESULT:"):
            try: score=float(line.split()[3])
            except (IndexError,ValueError): score=None
        elif line.startswith(("ATOM  ","HETATM")):
            if not model_seen or in_model:
                parts=line.split()
                atom_type=(parts[-1].upper() if parts else line[76:78].strip().upper())
                element=TYPE_TO_ELEMENT.get(atom_type)
                if element is None:
                    element=TYPE_TO_ELEMENT.get(line[76:78].strip().upper())
                if element is None: raise ValueError(f"Unknown PDBQT atom type {atom_type!r} in {path}: {line}")
                xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])],dtype=float)
                if not np.isfinite(xyz).all(): raise ValueError(f"Non-finite pose coordinate in {path}")
                if element not in {"H","D"}:
                    if element not in VDW: raise ValueError(f"No VDW radius for ligand element {element}")
                    current.append({"serial":int(line[6:11]),"atom_type":atom_type,"element":element,"xyz":xyz})
    if current: models.append((current,score))
    # PDBQT files occasionally have one model without MODEL/ENDMDL.
    if not models and current: models=[(current,score)]
    return models

def dist_matrix(xyz1, xyz2):
    a=np.asarray(xyz1,dtype=float); b=np.asarray(xyz2,dtype=float)
    return np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2)

def evaluate_pose(ligand_atoms, score, receptor_atoms, hotspots):
    if not ligand_atoms: raise ValueError("pose contains no peptide heavy atoms")
    lig_xyz=np.array([a["xyz"] for a in ligand_atoms],dtype=float)
    rec_xyz=np.array([a["xyz"] for a in receptor_atoms],dtype=float)
    d=dist_matrix(lig_xyz,rec_xyz)
    min_inter=float(d.min())
    lig_r=np.array([VDW[a["element"]] for a in ligand_atoms],dtype=float)
    rec_r=np.array([VDW[a["element"]] for a in receptor_atoms],dtype=float)
    overlap=lig_r[:,None]+rec_r[None,:]-d
    severe=int(np.count_nonzero(overlap>0.4))
    hotspot_mins={}
    hotspot_contact_atoms=set()
    for res in (203,205):
        hot_xyz=np.array([a["xyz"] for a in hotspots[res]],dtype=float)
        dh=dist_matrix(lig_xyz,hot_xyz)
        hotspot_mins[res]=float(dh.min())
        hotspot_contact_atoms.update(i for i in range(len(ligand_atoms)) if float(dh[i].min())<=3.5)
    min_hot=min(hotspot_mins.values())
    valid=math.isfinite(min_inter) and math.isfinite(min_hot) and score is not None and math.isfinite(float(score))
    pose_success=bool(valid and min_hot<=3.5 and severe==0)
    dual=bool(valid and hotspot_mins[203]<=3.5 and hotspot_mins[205]<=3.5 and severe==0)
    return {"pose_score_kcal_mol":score,"ligand_heavy_atom_count":len(ligand_atoms),"min_arg203_sidechain_A":hotspot_mins[203],"min_arg205_sidechain_A":hotspot_mins[205],"min_hotspot_A":min_hot,"hotspot_contacting_ligand_heavy_atom_count":len(hotspot_contact_atoms),"min_interchain_heavy_atom_A":min_inter,"severe_vdw_overlap_pair_count_gt0_4A":severe,"geometry_valid":int(valid),"pose_success":int(pose_success),"dual_arg_pose_success":int(dual),"failure_reason":"" if pose_success else ("invalid_geometry_or_missing_score" if not valid else "no_hotspot_contact" if min_hot>3.5 else "severe_vdw_overlap")}

def read_csv(path):
    with path.open(newline="") as f: return list(csv.DictReader(f))
def write_csv(path,rows,fields):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore"); w.writeheader(); w.writerows(rows)
def num(s):
    try:return float(s)
    except (TypeError,ValueError):return None

def main():
    if not RECEPTOR.is_file(): raise SystemExit(f"Missing receptor: {RECEPTOR}")
    selected_path=SV/"selected_top10.csv"; manifest_path=SV/"pose_manifest.csv"
    if not selected_path.is_file() or not manifest_path.is_file(): raise SystemExit("Missing selected_top10.csv or pose_manifest.csv; run extract_top10_structural_panel.py first.")
    receptor_atoms,hotspots=read_receptor(RECEPTOR)
    selected=read_csv(selected_path); planned=read_csv(manifest_path)
    if len(selected)!=20 or len(planned)!=100: raise SystemExit(f"Expected 20 selected candidates and 100 planned runs, got {len(selected)} and {len(planned)}")
    # Method/rank coverage must be exactly 10 per method; 5 docking repeats per candidate.
    counts=defaultdict(int)
    for r in planned: counts[(r["method"],int(r["rank"]))]+=1
    if len(counts)!=20 or any(v!=5 for v in counts.values()): raise SystemExit("Pose manifest does not contain exactly 5 repeats for every selected candidate.")
    pose_rows=[]; run_rows=[]; candidate_rows=[]
    pose_summary=defaultdict(list); run_summary={}
    for run in planned:
        out=ROOT/run["output_dir"]; status_path=out/"run_status.json"; pose_path=out/"vina_out.pdbqt"
        status=None; failure=""
        if not status_path.is_file(): failure="missing_run_status"
        else:
            try: status=json.loads(status_path.read_text())
            except Exception: failure="invalid_run_status_json"
        if status and not failure:
            if status.get("status")!="completed": failure="docking_not_completed"
            elif status.get("sequence")!=run["sequence"] or int(status.get("vina_seed",-1))!=int(run["vina_seed"]): failure="run_status_sequence_or_seed_mismatch"
            elif status.get("receptor_sha256")!=sha256(RECEPTOR): failure="receptor_hash_mismatch"
            elif not (SV/"receptor_meeko.pdbqt").is_file(): failure="prepared_receptor_missing"
            elif status.get("receptor_pdbqt_sha256")!=sha256(SV/"receptor_meeko.pdbqt"): failure="prepared_receptor_hash_mismatch"
            elif not pose_path.is_file(): failure="missing_pose_file"
            elif status.get("pose_sha256")!=sha256(pose_path): failure="pose_hash_mismatch"
        evaluations=[]
        if not failure:
            try:
                models=parse_pdbqt_models(pose_path)
                expected=int(status.get("ligand_heavy_atom_count",-1))
                if not models: raise ValueError("no pose models parsed")
                for pose_idx,(atoms,score) in enumerate(models, start=1):
                    result=evaluate_pose(atoms,score,receptor_atoms,hotspots)
                    result.update({"method":run["method"],"rank":run["rank"],"sequence":run["sequence"],"sequence_score":run["sequence_score"],"source_seed":run["source_seed"],"repeat":run["repeat"],"vina_seed":run["vina_seed"],"pose_index":pose_idx,"run_output_dir":run["output_dir"]})
                    if expected>0 and result["ligand_heavy_atom_count"]!=expected:
                        result["geometry_valid"]=0; result["pose_success"]=0; result["dual_arg_pose_success"]=0; result["failure_reason"]=f"ligand_heavy_atom_count_mismatch_expected_{expected}"
                    pose_rows.append(result); evaluations.append(result)
            except Exception as exc: failure=f"pose_evaluation_error:{type(exc).__name__}:{exc}"
        if failure:
            pose_rows.append({"method":run["method"],"rank":run["rank"],"sequence":run["sequence"],"sequence_score":run["sequence_score"],"source_seed":run["source_seed"],"repeat":run["repeat"],"vina_seed":run["vina_seed"],"pose_index":"","run_output_dir":run["output_dir"],"failure_reason":failure,"geometry_valid":0,"pose_success":0,"dual_arg_pose_success":0})
        success=any(int(p.get("pose_success",0))==1 for p in evaluations)
        dual=any(int(p.get("dual_arg_pose_success",0))==1 for p in evaluations)
        run_record={"method":run["method"],"rank":int(run["rank"]),"sequence":run["sequence"],"sequence_score":float(run["sequence_score"]),"source_seed":int(run["source_seed"]),"repeat":int(run["repeat"]),"vina_seed":int(run["vina_seed"]),"output_dir":run["output_dir"],"status":"completed" if not failure else "failed_or_missing","pose_count":len(evaluations),"run_success":int(success),"dual_arg_run_success":int(dual),"failure_reason":failure or ("" if success else "no_pose_met_geometric_success_criteria")}
        run_rows.append(run_record); run_summary[(run["method"],int(run["rank"]))] = run_summary.get((run["method"],int(run["rank"])),[])+[run_record]
    selected_by={ (r["method"],int(r["rank"])):r for r in selected }
    for method_rank,cand in sorted(selected_by.items()):
        runs=run_summary.get(method_rank,[])
        successes=sum(int(r["run_success"]) for r in runs); dual_successes=sum(int(r["dual_arg_run_success"]) for r in runs)
        candidate_success=int(len(runs)==5 and successes>=3); dual_candidate=int(len(runs)==5 and dual_successes>=3)
        candidate_rows.append({"method":cand["method"],"rank":cand["rank"],"sequence":cand["sequence"],"sequence_score":cand["sequence_score"],"source_seed":cand["source_seed"],"source_run_id":cand["source_run_id"],"generation":cand["generation"],"planned_runs":5,"present_run_records":len(runs),"geometric_successful_runs":successes,"dual_arg_successful_runs":dual_successes,"candidate_geometric_success":candidate_success,"candidate_dual_arg_success":dual_candidate})
    summary=[]
    for method in ("adaptive_substitution","uniform_baseline"):
        cr=[r for r in candidate_rows if r["method"]==method]; rr=[r for r in run_rows if r["method"]==method]
        prs=[r for r in pose_rows if r.get("method")==method and r.get("pose_score_kcal_mol") not in (None,"")]
        sr=sum(int(r["run_success"]) for r in rr); dsr=sum(int(r["dual_arg_run_success"]) for r in rr)
        cnum=sum(int(r["candidate_geometric_success"]) for r in cr); dnum=sum(int(r["candidate_dual_arg_success"]) for r in cr)
        vals=[float(r["pose_score_kcal_mol"]) for r in prs]
        mins=[float(r["min_interchain_heavy_atom_A"]) for r in prs if r.get("min_interchain_heavy_atom_A") not in (None,"")]
        clashes=[int(r["severe_vdw_overlap_pair_count_gt0_4A"]) for r in prs if r.get("severe_vdw_overlap_pair_count_gt0_4A") not in (None,"")]
        summary.append({"method":method,"candidates_in_frozen_panel":len(cr),"candidate_geometric_success_n":cnum,"candidate_geometric_success_denominator":10,"GSR_top10":cnum/10 if len(cr)==10 else None,"candidate_dual_arg_success_n":dnum,"candidate_dual_arg_denominator":10,"dual_arg_candidate_rate":dnum/10 if len(cr)==10 else None,"successful_runs_n":sr,"planned_runs_denominator":50,"SR_run":sr/50 if len(rr)==50 else None,"dual_arg_successful_runs_n":dsr,"mean_pose_vina_score_kcal_mol":statistics.mean(vals) if vals else None,"median_pose_vina_score_kcal_mol":statistics.median(vals) if vals else None,"pose_min_interchain_distance_median_A":statistics.median(mins) if mins else None,"pose_severe_clash_pairs_total":sum(clashes),"failed_or_missing_runs_n":sum(r["status"]!="completed" for r in rr),"interpretation":"Vina-only geometric screen; not experimental binding evidence"})
    write_csv(SV/"pose_results.csv",pose_rows,["method","rank","sequence","sequence_score","source_seed","repeat","vina_seed","pose_index","pose_score_kcal_mol","ligand_heavy_atom_count","min_arg203_sidechain_A","min_arg205_sidechain_A","min_hotspot_A","hotspot_contacting_ligand_heavy_atom_count","min_interchain_heavy_atom_A","severe_vdw_overlap_pair_count_gt0_4A","geometry_valid","pose_success","dual_arg_pose_success","failure_reason","run_output_dir"])
    write_csv(SV/"run_results.csv",run_rows,["method","rank","sequence","sequence_score","source_seed","repeat","vina_seed","output_dir","status","pose_count","run_success","dual_arg_run_success","failure_reason"])
    write_csv(SV/"candidate_results.csv",candidate_rows,["method","rank","sequence","sequence_score","source_seed","source_run_id","generation","planned_runs","present_run_records","geometric_successful_runs","dual_arg_successful_runs","candidate_geometric_success","candidate_dual_arg_success"])
    write_csv(SV/"method_summary.csv",summary,["method","candidates_in_frozen_panel","candidate_geometric_success_n","candidate_geometric_success_denominator","GSR_top10","candidate_dual_arg_success_n","candidate_dual_arg_denominator","dual_arg_candidate_rate","successful_runs_n","planned_runs_denominator","SR_run","dual_arg_successful_runs_n","mean_pose_vina_score_kcal_mol","median_pose_vina_score_kcal_mol","pose_min_interchain_distance_median_A","pose_severe_clash_pairs_total","failed_or_missing_runs_n","interpretation"])
    print((SV/"method_summary.csv").read_text())
    print(f"Wrote geometric metrics to {SV}; denominator is 10 candidates and 50 runs per method, including failures/missing runs.")

if __name__=="__main__": main()
