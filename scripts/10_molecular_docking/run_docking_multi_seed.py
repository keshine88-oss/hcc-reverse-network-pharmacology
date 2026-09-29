# -*- coding: utf-8 -*-
"""
run_multi_seed.py
Multi-seed stability check: re-dock every ligand with 4 seeds (42, 1, 7, 2024)
and report per-ligand best-score range. Docking runs in a thread pool (parallelism
on the outer loop; Vina cpu=1 per task).
"""
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import subprocess
from concurrent.futures import ThreadPoolExecutor

VINA = r"G:\AT1-HF\Docking\vina.exe"
BASE = os.path.join(PROJECT_ROOT, "work/drug/docking")
REC_DIR = f"{BASE}/receptor"
LIG_DIR = f"{BASE}/ligand"
OUT_DIR = f"{BASE}/results"
SEEDS = [42, 1, 7, 2024]

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


def dock(target, ligand, box_size, center, rec_pdbqt, seed):
    out_dir = f"{OUT_DIR}/{target}"
    os.makedirs(out_dir, exist_ok=True)
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
    config.encode("ascii")
    cfg = os.path.join(out_dir, f"{ligand}_seed{seed}.conf")
    out_pdbqt = os.path.join(out_dir, f"{ligand}_seed{seed}_out.pdbqt")
    log = os.path.join(out_dir, f"{ligand}_seed{seed}.log")
    with open(cfg, "w") as fh:
        fh.write(config)
    r = subprocess.run([VINA, "--config", cfg, "--out", out_pdbqt, "--log", log],
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out_pdbqt):
        return None
    scores = []
    for line in open(out_pdbqt):
        if line.startswith("REMARK VINA RESULT:"):
            p = line.split()
            try:
                scores.append(float(p[3]))
            except (IndexError, ValueError):
                pass
    return scores[0] if scores else None


def run_ligand(args):
    target, ligand, box_size, center, rec_pdbqt = args
    vals = {}
    for seed in SEEDS:
        vals[seed] = dock(target, ligand, box_size, center, rec_pdbqt, seed)
    return ligand, vals


def main():
    tasks = []
    for target, ligs in LIGAND_MAP.items():
        rec_pdbqt, site_json = RECEPTOR[target]
        center = json.load(open(f"{REC_DIR}/{site_json}", encoding="utf-8"))["box_center"]
        for lig, box in ligs.items():
            tasks.append((target, lig, box, center, rec_pdbqt))

    results = {}
    with ThreadPoolExecutor(max_workers=4) as ex:
        for lig, vals in ex.map(run_ligand, tasks):
            results[lig] = vals
            sc = [v for v in vals.values() if v is not None]
            rng = (max(sc) - min(sc)) if sc else None
            print(f"{lig:24s} " + "  ".join(f"s{k}={v}" for k, v in vals.items()) + f"   range={rng}")

    # summary
    with open(f"{OUT_DIR}/multi_seed_summary.txt", "w") as fh:
        for lig, vals in sorted(results.items()):
            sc = [v for v in vals.values() if v is not None]
            rng = (max(sc) - min(sc)) if sc else None
            fh.write(f"{lig}\t" + "\t".join(str(vals[s]) for s in SEEDS) + f"\trange={rng}\n")
    print(f"\nsummary -> {OUT_DIR}/multi_seed_summary.txt")


if __name__ == "__main__":
    main()
