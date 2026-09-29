# -*- coding: utf-8 -*-
"""
run_docking.py
Run AutoDock Vina 1.1.2 docking via subprocess (no Python bindings in 1.1.2).

Usage:
  python run_docking.py PTH1R   # dock the 8 PTH1R-target ligands
  python run_docking.py SLCO4C1 # dock the 2 SLCO4C1-target ligands

Box centers come from *_site.json; config is pure-ASCII; ligand filenames ASCII.
Results -> docking_results/<target>/<ligand>_out.pdbqt + <ligand>.log
"""
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import subprocess
import sys

VINA = r"G:\AT1-HF\Docking\vina.exe"
BASE = os.path.join(PROJECT_ROOT, "work/drug/docking")
REC_DIR = f"{BASE}/receptor"
LIG_DIR = f"{BASE}/ligand"
OUT_DIR = f"{BASE}/results"

# ligand -> (target, box_size)
LIGAND_MAP = {
    "PTH1R": {
        "Kamebanin": 22, "Anisaldehyde": 22, "Ginkgolide_NPC88890": 22,
        "Cucurbitacin_B": 22, "Chrysin": 22, "Roemerine": 22,
        "Tryptanthrin": 22, "Digoxigenin": 22,
    },
    "SLCO4C1": {
        "Digitoxin": 30, "Digoxigenin": 30,
    },
}

RECEPTOR = {
    "PTH1R": ("6NBF_rec.pdbqt", "6NBF_site.json"),
    "SLCO4C1": ("AF_Q6ZQN7_rec.pdbqt", "AF-Q6ZQN7_site.json"),
}


def dock_one(target, ligand, box_size, center, rec_pdbqt, out_dir, seed=42):
    config = (
        f"receptor = {REC_DIR}/{rec_pdbqt}\n"
        f"ligand = {LIG_DIR}/{ligand}.pdbqt\n"
        f"center_x = {center[0]:.3f}\n"
        f"center_y = {center[1]:.3f}\n"
        f"center_z = {center[2]:.3f}\n"
        f"size_x = {box_size}\n"
        f"size_y = {box_size}\n"
        f"size_z = {box_size}\n"
        f"exhaustiveness = 16\n"
        f"num_modes = 9\n"
        f"seed = {seed}\n"
    )
    config.encode("ascii")  # will raise if non-ASCII sneaks in
    cfg_path = os.path.join(out_dir, f"{ligand}.conf")
    with open(cfg_path, "w") as fh:
        fh.write(config)

    out_pdbqt = os.path.join(out_dir, f"{ligand}_out.pdbqt")
    log_path = os.path.join(out_dir, f"{ligand}.log")
    cmd = [VINA, "--config", cfg_path, "--out", out_pdbqt, "--log", log_path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r, out_pdbqt, log_path


def parse_scores(out_pdbqt):
    """parse best score per mode from Vina output pdbqt REMARK lines."""
    scores = []
    with open(out_pdbqt, encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                # format: REMARK VINA RESULT:  -7.6  0.000  0.000
                try:
                    scores.append((float(parts[3]), float(parts[4]), float(parts[5])))
                except (IndexError, ValueError):
                    pass
    return scores


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "PTH1R"
    lig_map = LIGAND_MAP[target]
    rec_pdbqt, site_json = RECEPTOR[target]
    site = json.load(open(f"{REC_DIR}/{site_json}", encoding="utf-8"))
    center = site["box_center"]

    out_dir = f"{OUT_DIR}/{target}"
    os.makedirs(out_dir, exist_ok=True)

    print(f"=== docking target {target}: {len(lig_map)} ligands ===")
    print(f"receptor {rec_pdbqt}; center {center}")
    results = []
    for ligand, box_size in lig_map.items():
        r, out_pdbqt, log_path = dock_one(target, ligand, box_size, center, rec_pdbqt, out_dir)
        if r.returncode != 0 or not os.path.exists(out_pdbqt):
            print(f"[FAIL] {ligand}: rc={r.returncode}")
            err = (r.stderr or "")[-400:]
            print("   ", err.replace("\n", " ")[:400])
            continue
        scores = parse_scores(out_pdbqt)
        best = scores[0][0] if scores else None
        results.append((ligand, best))
        print(f"[OK] {ligand}: best {best} kcal/mol ({len(scores)} modes)")
    print("\n=== ranking (best kcal/mol, more negative = stronger) ===")
    for lig, sc in sorted(results, key=lambda x: (x[1] is None, x[1] if x[1] is not None else 0)):
        print(f"  {lig:24s} {sc}")
    # save summary
    with open(f"{out_dir}/_scores.txt", "w") as fh:
        for lig, sc in sorted(results, key=lambda x: (x[1] is None, x[1] if x[1] is not None else 0)):
            fh.write(f"{lig}\t{sc}\n")
    print(f"scores -> {out_dir}/_scores.txt")


if __name__ == "__main__":
    main()
