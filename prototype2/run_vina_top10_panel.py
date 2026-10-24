#!/usr/bin/env python3
"""Exploratory Vina-only docking for a frozen top-10-per-method peptide panel.

Requires RDKit, Meeko, NumPy and a Meeko mk_prepare_receptor.py script in one Python environment.
This does NOT emulate RAPiDock/HybriDock, nor does a Vina score establish binding.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, shutil, subprocess, sys, time, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SV = ROOT / "prototype2" / "structural_validation"
MANIFEST = SV / "pose_manifest.csv"
RECEPTOR_PDB = ROOT / "prototype2" / "target_aware_audit" / "9SYA_chainA_repaired.pdb"
VINA_DEFAULT = ROOT / "vina_1.2.7_mac_aarch64"

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""): h.update(chunk)
    return h.hexdigest()

def receptor_atoms(path: Path):
    atoms = []
    hotspot = {203: [], 205: []}
    sidechain = {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  "): continue
        name = line[12:16].strip()
        chain = line[21:22].strip()
        resname = line[17:20].strip()
        resseq = int(line[22:26])
        element = line[76:78].strip().upper() or name[:1].upper()
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if not all(math.isfinite(x) for x in xyz): raise ValueError("nonfinite receptor coordinate")
        if chain == "A" and resseq in hotspot and resname == "ARG" and name in sidechain:
            hotspot[resseq].append(xyz)
        if element not in {"H", "D"}:
            atoms.append({"name":name,"chain":chain,"resname":resname,"resseq":resseq,"element":element,"xyz":xyz})
    if not hotspot[203] or not hotspot[205]: raise ValueError("chain A Arg203/Arg205 side-chain mapping failed")
    center = [sum(x[k] for r in hotspot.values() for x in r) / sum(len(v) for v in hotspot.values()) for k in range(3)]
    return atoms, hotspot, center

def validate_receptor_pdbqt_frame(source_pdb: Path, receptor_pdbqt: Path) -> dict:
    source = {}
    wanted = {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}
    for line in source_pdb.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[21:22].strip() != "A": continue
        resseq = int(line[22:26]); name = line[12:16].strip()
        if resseq in {203, 205} and line[17:20].strip() == "ARG" and name in wanted:
            source[(resseq, name)] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    prepared = {}
    for line in receptor_pdbqt.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")): continue
        resname=line[17:20].strip(); chain=line[21:22].strip(); resseq=int(line[22:26]); name=line[12:16].strip()
        if chain=="A" and resseq in {203,205} and resname=="ARG" and name in wanted:
            prepared[(resseq,name)] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
    if not source or set(source) != set(prepared):
        raise RuntimeError(f"Cannot establish prepared receptor hotspot atom mapping; source={sorted(source)}, prepared={sorted(prepared)}")
    max_dev=max(math.dist(source[k],prepared[k]) for k in source)
    if max_dev > 0.20:
        raise RuntimeError(f"Prepared receptor coordinates appear shifted from source frame (max Arg203/205 atom deviation {max_dev:.3f} A)")
    return {"mapped_hotspot_atom_count":len(source),"max_hotspot_atom_deviation_A":max_dev}

def prepare_ligand(sequence: str, seed: int, output: Path) -> dict:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mol = Chem.MolFromSequence(sequence)
    if mol is None: raise RuntimeError(f"RDKit cannot construct peptide sequence: {sequence}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3(); params.randomSeed = int(seed); params.maxIterations = 100
    # Long linear peptides commonly fail distance-geometry embedding with ETKDG's default coordinate initialization.
    # Random-coordinate embedding succeeded on the tested 21-mer RDKit smoke case; keep this deterministic and logged.
    params.useRandomCoords = True
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        params.randomSeed = int(seed) + 1
        status = AllChem.EmbedMolecule(mol, params)
    if status != 0: raise RuntimeError(f"RDKit ETKDG random-coordinate embedding failed with status {status} on two seeds")
    if AllChem.MMFFHasAllMoleculeParams(mol): AllChem.MMFFOptimizeMolecule(mol, maxIters=400)
    setups = MoleculePreparation().prepare(mol)
    if len(setups) != 1: raise RuntimeError(f"Expected one Meeko setup, got {len(setups)}")
    pdbqt, success, error = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
    if not success: raise RuntimeError(f"Meeko ligand preparation failed: {error}")
    output.parent.mkdir(parents=True, exist_ok=True); output.write_text(pdbqt)
    # Ensure the emitted ligand includes a sensible heavy-atom count for a 21-residue peptide.
    heavy = sum(1 for line in pdbqt.splitlines() if line.startswith(("ATOM  ", "HETATM")) and line.split()[-1].upper() not in {"H", "HD", "HS", "D"})
    if heavy < len(sequence) * 3: raise RuntimeError(f"Prepared peptide has unexpectedly few heavy atoms: {heavy}")
    return {"sequence":sequence,"conformer_seed":seed,"ligand_heavy_atom_count":heavy,"ligand_pdbqt_sha256":sha256(output)}

def run_one(row: dict, root: Path, receptor_pdbqt: Path, receptor_hash: str, vina: Path, box_center: list[float], box_size: float, exhaustiveness: int, num_modes: int, energy_range: float, ligands: dict[str, tuple[Path,dict]], retry_failed: bool=False) -> dict:
    out = root / row["output_dir"]
    status_file = out / "run_status.json"
    if status_file.exists():
        try: old = json.loads(status_file.read_text())
        except Exception: old = {}
        if old.get("status") == "completed" and (out / "vina_out.pdbqt").is_file():
            return {"status":"skipped_completed","output_dir":str(out),"sequence":row["sequence"],"seed":int(row["vina_seed"])}
        if not retry_failed and old:
            return {"status":"skipped_existing_status_not_completed","output_dir":str(out),"sequence":row["sequence"],"seed":int(row["vina_seed"])}
    out.mkdir(parents=True, exist_ok=True)
    ligand_file, ligand_meta = ligands[row["method"] + "|" + row["sequence"]]
    pose = out / "vina_out.pdbqt"
    command = [str(vina), "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_file),
        "--center_x", f"{box_center[0]:.3f}", "--center_y", f"{box_center[1]:.3f}", "--center_z", f"{box_center[2]:.3f}",
        "--size_x", str(box_size), "--size_y", str(box_size), "--size_z", str(box_size),
        "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes), "--energy_range", str(energy_range),
        "--cpu", "1", "--seed", str(int(row["vina_seed"])), "--out", str(pose)]
    payload = {"status":"running","method":row["method"],"rank":int(row["rank"]),"sequence":row["sequence"],"sequence_score":float(row["sequence_score"]),"source_seed":int(row["source_seed"]),"repeat":int(row["repeat"]),"vina_seed":int(row["vina_seed"]),"receptor_pdb":str(RECEPTOR_PDB),"receptor_sha256":sha256(RECEPTOR_PDB),"receptor_pdbqt":str(receptor_pdbqt),"receptor_pdbqt_sha256":receptor_hash,"ligand_pdbqt":str(ligand_file),**ligand_meta,"vina_executable":str(vina),"vina_executable_sha256":sha256(vina),"center_A":box_center,"box_size_A":[box_size]*3,"exhaustiveness":exhaustiveness,"num_modes":num_modes,"energy_range_kcal_mol":energy_range,"command":command}
    status_file.write_text(json.dumps(payload, indent=2) + "\n")
    started=time.time()
    try:
        proc=subprocess.run(command, capture_output=True, text=True, check=False)
        (out/"vina.stdout.txt").write_text(proc.stdout); (out/"vina.stderr.txt").write_text(proc.stderr)
        payload.update({"returncode":proc.returncode,"runtime_seconds":round(time.time()-started,3),"stdout_sha256":sha256(out/"vina.stdout.txt"),"stderr_sha256":sha256(out/"vina.stderr.txt")})
        if proc.returncode == 0 and pose.is_file() and pose.stat().st_size > 100:
            text=pose.read_text(errors="replace")
            models=text.count("MODEL")
            if models == 0: models=1
            scores=[]
            for line in text.splitlines():
                if line.startswith("REMARK VINA RESULT:"):
                    try: scores.append(float(line.split()[3]))
                    except Exception: pass
            if not scores: raise RuntimeError("Vina returned a pose file but no REMARK VINA RESULT scores")
            payload.update({"status":"completed","pose_count":models,"vina_scores_kcal_mol":scores,"pose_sha256":sha256(pose)})
        else:
            payload.update({"status":"failed","failure_reason":f"Vina exit={proc.returncode}; pose exists={pose.is_file()}"})
    except Exception as exc:
        payload.update({"status":"failed","failure_reason":f"{type(exc).__name__}: {exc}","traceback":traceback.format_exc(),"runtime_seconds":round(time.time()-started,3)})
    status_file.write_text(json.dumps(payload, indent=2)+"\n")
    return {"status":payload["status"],"method":row["method"],"rank":int(row["rank"]),"seed":int(row["vina_seed"]),"sequence":row["sequence"],"output_dir":str(out),"runtime_seconds":payload.get("runtime_seconds"),"failure_reason":payload.get("failure_reason","")}

def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--python", type=Path, default=Path(sys.executable), help="Python in the RDKit/Meeko environment")
    ap.add_argument("--meeko-receptor-prep", type=Path, default=None)
    ap.add_argument("--vina", type=Path, default=VINA_DEFAULT)
    ap.add_argument("--smoke-test", action="store_true", help="run exactly the first manifest row; no others")
    ap.add_argument("--run-all", action="store_true", help="run the full manifest after a successful smoke test exists")
    ap.add_argument("--workers", type=int, default=4, help="parallel Vina processes; each Vina process uses one CPU")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--box-size", type=float, default=40.0)
    ap.add_argument("--retry-failed", action="store_true")
    args=ap.parse_args()
    if args.smoke_test == args.run_all:
        raise SystemExit("Choose exactly one of --smoke-test or --run-all")
    if args.workers<1 or args.exhaustiveness<1 or args.num_modes<1 or args.box_size<=0: raise SystemExit("Invalid worker/search parameters")
    for p in (RECEPTOR_PDB, MANIFEST, args.vina):
        if not p.is_file(): raise SystemExit(f"Required file missing: {p}")
    vina=args.vina.resolve()
    receptor_atoms(RECEPTOR_PDB)  # fail closed on hotspot residue identity
    try:
        from rdkit import Chem  # noqa: F401
        from meeko import MoleculePreparation, PDBQTWriterLegacy  # noqa: F401
        import numpy  # noqa: F401
    except Exception as exc: raise SystemExit(f"RDKit/Meeko/NumPy imports unavailable in {sys.executable}: {exc}")
    receptor_prep=args.meeko_receptor_prep
    if receptor_prep is None:
        candidate=Path(sys.executable).resolve().parent / "mk_prepare_receptor.py"
        if candidate.is_file(): receptor_prep=candidate
        else: raise SystemExit(f"Cannot find mk_prepare_receptor.py next to {sys.executable}; pass --meeko-receptor-prep /absolute/path/to/mk_prepare_receptor.py")
    if not receptor_prep.is_file(): raise SystemExit(f"Meeko receptor-prep script missing: {receptor_prep}")
    rows=list(csv.DictReader(MANIFEST.open(newline="")))
    if len(rows)!=100: raise SystemExit(f"Expected 100 pose-manifest rows (20 candidates x 5 seeds), found {len(rows)}")
    if args.run_all:
        smoke_rows=[r for r in rows if int(r["rank"])==1]
        if not smoke_rows: raise SystemExit("No rank-1 row to check for smoke completion")
        statuses=[]
        for r in smoke_rows:
            s=ROOT / r["output_dir"] / "run_status.json"
            if s.exists():
                try: statuses.append(json.loads(s.read_text()).get("status"))
                except Exception: statuses.append(None)
        if "completed" not in statuses:
            raise SystemExit("BLOCKED: no successful rank-1 smoke run found. Run --smoke-test and inspect diagnostics before --run-all.")
    SV.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt=SV/"receptor_meeko.pdbqt"
    prep_log=SV/"receptor_prepare.log"
    protocol_path=SV/"vina_protocol.json"
    source_receptor_hash=sha256(RECEPTOR_PDB)
    if receptor_pdbqt.is_file():
        if not protocol_path.is_file():
            raise SystemExit(f"A receptor PDBQT exists without a provenance file: {receptor_pdbqt}. Move it aside and rerun; do not reuse an unverified receptor.")
        old_protocol=json.loads(protocol_path.read_text())
        if old_protocol.get("receptor_sha256")!=source_receptor_hash or old_protocol.get("receptor_pdbqt_sha256")!=sha256(receptor_pdbqt):
            raise SystemExit("Existing receptor PDBQT/protocol hash mismatch. Preserve the old files, move them aside, and regenerate receptor prep.")
    else:
        cmd=[str(args.python.resolve()),str(receptor_prep.resolve()),"--read_pdb",str(RECEPTOR_PDB),"-o",str(SV/"receptor_meeko"),"-p",str(receptor_pdbqt)]
        proc=subprocess.run(cmd,capture_output=True,text=True,check=False)
        prep_log.write_text("COMMAND: "+" ".join(cmd)+"\n\nSTDOUT:\n"+proc.stdout+"\nSTDERR:\n"+proc.stderr)
        if proc.returncode!=0 or not receptor_pdbqt.is_file(): raise SystemExit(f"Meeko receptor prep failed; see {prep_log}")
    receptor_hash=sha256(receptor_pdbqt)
    rec_atoms, hot, center=receptor_atoms(RECEPTOR_PDB)
    frame_check=validate_receptor_pdbqt_frame(RECEPTOR_PDB,receptor_pdbqt)
    (SV/"vina_protocol.json").write_text(json.dumps({"scope":"Vina-only exploratory structural screen; separate from RAPiDock/HybriDock","receptor":str(RECEPTOR_PDB),"receptor_sha256":source_receptor_hash,"receptor_pdbqt_sha256":receptor_hash,"receptor_chain":"A","hotspot":"ARG A203/A205 side-chain heavy atoms","receptor_frame_check":frame_check,"box_center_A":center,"box_size_A":[args.box_size]*3,"exhaustiveness":args.exhaustiveness,"num_modes":args.num_modes,"energy_range_kcal_mol":5.0,"workers":args.workers,"per_vina_cpu":1,"seeds":"from pose_manifest.csv","interpretation":"docking scores are scoring-function outputs, not experimental affinity"},indent=2)+"\n")
    target_rows=rows[:1] if args.smoke_test else rows
    # Prepare one deterministic conformer for each sequence actually used in this invocation.
    selected_for_prep=list(csv.DictReader((SV/"selected_top10.csv").open(newline="")))
    if args.smoke_test:
        first=target_rows[0]
        selected_for_prep=[c for c in selected_for_prep if c["method"]==first["method"] and c["sequence"]==first["sequence"]]
    ligands={}
    for cand in selected_for_prep:
        key=cand["method"]+"|"+cand["sequence"]
        if key not in ligands:
            ligpath=SV/"ligands"/(f"{cand['method']}_rank{int(cand['rank']):02d}.pdbqt")
            ligmeta=prepare_ligand(cand["sequence"], 55100+int(cand["rank"])+(0 if cand["method"]=="adaptive_substitution" else 100), ligpath)
            ligands[key]=(ligpath,ligmeta)
    results=[]
    with ThreadPoolExecutor(max_workers=1 if args.smoke_test else args.workers) as ex:
        futures=[ex.submit(run_one,r,ROOT,receptor_pdbqt,receptor_hash,vina,center,args.box_size,args.exhaustiveness,args.num_modes,5.0,ligands,args.retry_failed) for r in target_rows]
        for fut in as_completed(futures):
            row=fut.result(); results.append(row)
            print(f"{row['status']:30s} {row.get('method')} rank={row.get('rank')} seed={row.get('seed')} {row.get('runtime_seconds')}s",flush=True)
    log_path=SV/("smoke_test_summary.json" if args.smoke_test else "vina_run_summary.json")
    log_path.write_text(json.dumps({"run_count":len(results),"results":results},indent=2)+"\n")
    successes=sum(r["status"] in {"completed","skipped_completed"} for r in results)
    print(f"Finished {len(results)} jobs; successful/completed-or-skipped={successes}. Summary: {log_path}")
    if args.smoke_test and successes!=1: return 2
    return 0 if successes==len(results) else 1

if __name__=="__main__": raise SystemExit(main())
