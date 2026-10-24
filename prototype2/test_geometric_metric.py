#!/usr/bin/env python3
"""Small synthetic tests for the geometry predicates; these are not docking validation."""
import importlib.util
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("evaluate_geometric_success", HERE / "evaluate_geometric_success.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

def atom(name, element, xyz): return {"name":name,"element":element,"xyz":np.array(xyz,dtype=float)}
def run_tests():
    hot203=[atom(n,"N" if n in {"NE","NH1","NH2"} else "C", [0.0,0.0,0.0]) for n in ("CB","CG","CD","NE","CZ","NH1","NH2")]
    hot205=[atom(n,"N" if n in {"NE","NH1","NH2"} else "C", [20.0,0.0,0.0]) for n in ("CB","CG","CD","NE","CZ","NH1","NH2")]
    receptor=hot203+hot205+[atom("FAR","C",[100,100,100])]
    hotspots={203:hot203,205:hot205}
    # N-to-O distance of 3.1 A is a hotspot proximity pass without a severe VDW overlap.
    ligand=[atom("O","O",[3.1,0,0])]
    passed=mod.evaluate_pose(ligand,-5.0,receptor,hotspots)
    assert passed["pose_success"]==1, passed
    assert passed["dual_arg_pose_success"]==0, passed
    far=mod.evaluate_pose([atom("O","O",[50,50,50])],-5.0,receptor,hotspots)
    assert far["pose_success"]==0, far
    overlap=mod.evaluate_pose([atom("C","C",[0,0,0])],-5.0,receptor,hotspots)
    assert overlap["pose_success"]==0 and overlap["severe_vdw_overlap_pair_count_gt0_4A"]>0, overlap
    invalid=mod.evaluate_pose([atom("O","O",[3.1,0,0])],None,receptor,hotspots)
    assert invalid["geometry_valid"]==0 and invalid["pose_success"]==0, invalid
    print("PASS: contact/no-clash, no hotspot contact, severe-overlap rejection, missing score invalidation")

if __name__=="__main__": run_tests()
