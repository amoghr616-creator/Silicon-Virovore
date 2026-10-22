from pathlib import Path
import argparse, csv, json, math
import numpy as np

WATER = {"HOH", "WAT", "DOD"}

def read_atoms(path, chain=None):
    atoms = []
    for line in path.read_text(errors="replace").splitlines():
        if line[:6].strip() not in {"ATOM", "HETATM"} or len(line) < 54:
            continue
        ch = line[21:22].strip()
        if chain is not None and ch != chain:
            continue
        residue = line[17:20].strip().upper()
        if residue in WATER:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = line[12:16].strip().upper().lstrip("0123456789")[:1]
        if element in {"", "H", "D"}:
            continue
        try:
            xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        except ValueError:
            continue
        if not all(math.isfinite(v) for v in xyz):
            continue
        atoms.append({
            "atom": line[12:16].strip(), "resname": residue,
            "chain": ch, "resid": line[22:26].strip(),
            "icode": line[26:27].strip(), "xyz": xyz,
        })
    return atoms

def label(a):
    return f'{a["chain"]}:{a["resname"]}{a["resid"]}{a["icode"].strip()}:{a["atom"]}'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", type=Path,
        default=Path("prototype2/target_aware_audit/9SYA_chainA_repaired.pdb"))
    ap.add_argument("--poses-root", type=Path,
        default=Path("prototype2/p2_rapidock"))
    ap.add_argument("--out", type=Path,
        default=Path("prototype2/p2_rapidock/clash_forensics.json"))
    ap.add_argument("--chain", default="A")
    ap.add_argument("--site-center", nargs=3, type=float,
        default=[-8.405, -5.823, 1.523])
    args = ap.parse_args()

    if not args.receptor.is_file():
        raise SystemExit(f"Receptor not found: {args.receptor}")
    if not args.poses_root.is_dir():
        raise SystemExit(f"Pose directory not found: {args.poses_root}")

    receptor = read_atoms(
        args.receptor, None if args.chain.upper() == "ALL" else args.chain
    )
    if not receptor:
        raise SystemExit("No receptor heavy atoms parsed; check chain/path.")

    files = sorted(p for p in args.poses_root.rglob("pose_*.pdb")
                   if "reverseprocess" not in p.name.lower())
    if not files:
        files = sorted(p for p in args.poses_root.rglob("rank[0-9]*.pdb")
                       if "reverseprocess" not in p.name.lower())
    if not files:
        raise SystemExit(f"No pose_*.pdb or rankN.pdb files under {args.poses_root}")

    R = np.asarray([a["xyz"] for a in receptor], dtype=float)
    site = np.asarray(args.site_center, dtype=float)
    summary, details = [], []

    for path in files:
        ligand = read_atoms(path)
        if not ligand:
            summary.append({"pose": str(path), "status": "NO_HEAVY_ATOMS"})
            continue

        L = np.asarray([a["xyz"] for a in ligand], dtype=float)
        D = np.linalg.norm(L[:, None, :] - R[None, :, :], axis=2)
        i_min, j_min = np.unravel_index(np.argmin(D), D.shape)
        clash_idx = np.argwhere(D < 2.0)
        clashes = sorted(
            (float(D[i, j]), int(i), int(j)) for i, j in clash_idx
        )

        def target_min(resid):
            indices = [
                j for j, a in enumerate(receptor)
                if a["resname"] == "ARG" and a["resid"] == str(resid)
            ]
            return round(float(D[:, indices].min()), 4) if indices else None

        def pair_record(distance, i, j):
            return {
                "distance_A": round(distance, 4),
                "peptide_atom": label(ligand[i]),
                "receptor_atom": label(receptor[j]),
                "peptide_xyz": ligand[i]["xyz"],
                "receptor_xyz": receptor[j]["xyz"],
            }

        centroid = L.mean(axis=0)
        row = {
            "pose": str(path), "status": "OK",
            "ligand_heavy_atoms": len(ligand),
            "minimum_distance_A": round(float(D[i_min, j_min]), 4),
            "pairs_lt_1_5A": int((D < 1.5).sum()),
            "pairs_lt_2A": int((D < 2.0).sum()),
            "pairs_lt_3A": int((D < 3.0).sum()),
            "pairs_le_4A": int((D <= 4.0).sum()),
            "ARG203_min_A": target_min(203),
            "ARG205_min_A": target_min(205),
            "peptide_centroid_A": [round(float(v), 3) for v in centroid],
            "centroid_to_site_A": round(float(np.linalg.norm(centroid-site)), 3),
        }
        closest_flat = np.argsort(D, axis=None)[:min(10, D.size)]
        closest = [
            pair_record(float(D.flat[k]), int(k // D.shape[1]), int(k % D.shape[1]))
            for k in closest_flat
        ]
        details.append({
            **row,
            "closest_10_pairs": closest,
            "clashes_lt_2A_examples": [
                pair_record(d, i, j) for d, i, j in clashes[:20]
            ],
        })
        summary.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "method": "Independent heavy-atom distances; no docking scores used",
        "receptor": str(args.receptor),
        "receptor_chain": args.chain,
        "site_center_A": args.site_center,
        "pose_files_found": len(files),
        "poses_with_pairs_lt_2A": sum(
            r.get("pairs_lt_2A", 0) > 0 for r in summary
        ),
        "poses": details,
    }, indent=2) + "\n")

    csv_path = args.out.with_suffix(".csv")
    columns = [
        "pose", "status", "ligand_heavy_atoms", "minimum_distance_A",
        "pairs_lt_1_5A", "pairs_lt_2A", "pairs_lt_3A", "pairs_le_4A",
        "ARG203_min_A", "ARG205_min_A", "centroid_to_site_A",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary)

    valid = [r for r in summary if r.get("status") == "OK"]
    print(f"Pose files: {len(files)} | Valid: {len(valid)}")
    if valid:
        print("Poses with heavy-atom pairs under 2 Å:",
              sum(r["pairs_lt_2A"] > 0 for r in valid))
        print("Median minimum distance (Å):",
              round(float(np.median([r["minimum_distance_A"] for r in valid])), 3))
        print("Median pairs under 2 Å:",
              float(np.median([r["pairs_lt_2A"] for r in valid])))
        print("Median ARG203 distance (Å):",
              round(float(np.nanmedian([
                  r["ARG203_min_A"] if r["ARG203_min_A"] is not None else np.nan
                  for r in valid
              ])), 3))
        print("Median ARG205 distance (Å):",
              round(float(np.nanmedian([
                  r["ARG205_min_A"] if r["ARG205_min_A"] is not None else np.nan
                  for r in valid
              ])), 3))
    print(f"JSON report: {args.out}")
    print(f"CSV report:  {csv_path}")

if __name__ == "__main__":
    main()
