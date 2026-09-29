# -*- coding: utf-8 -*-
"""
analyze_hbonds.py
Strict hydrogen-bond analysis between receptor (Meeko pdbqt, polar H kept) and
the docked ligand poses (Vina output pdbqt, mode 1).

Criterion (angle-aware, not distance-only):
  donor   = N/O bearing polar H (within 1.3 A)
  acceptor= all O + N without H
  H-bond  = D...A <= 3.5 A  AND  H...A <= 2.6 A  AND  angle(D-H...A) >= 120 deg

Ligand donor capacity is counted first (O-H / N-H); ligands with no donor can
only accept H-bonds from the receptor side.
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
import numpy as np
import gemmi

RES_DIR = os.path.join(PROJECT_ROOT, "work/drug/docking/results")
REC_PDBQT = {
    "PTH1R": os.path.join(PROJECT_ROOT, "work/drug/docking/receptor/6NBF_rec.pdbqt"),
    "SLCO4C1": os.path.join(PROJECT_ROOT, "work/drug/docking/receptor/AF_Q6ZQN7_rec.pdbqt"),
}


def read_pdbqt_atoms(path):
    atoms = []  # (pos np.array, element, resid_label, atom_name, has_polar_h)
    model_count = 0
    with open(path, encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                model_count += 1
                if model_count > 1:
                    break  # only first model (best pose)
                continue
            if line.startswith("ENDMDL"):
                continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            # PDBQT columns: 30-38 x, 38-46 y, 46-54 z; 17-20 atom name; 23-26 resname; 21 chain; 22-26 resi
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            aname = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            resi = line[22:26].strip()
            elem = line[77:79].strip() or aname[0]
            atoms.append({
                "pos": np.array([x, y, z]),
                "elem": elem,
                "aname": aname,
                "label": f"{resname} {chain}{resi}",
            })
    return atoms


def find_polar_h(atom, all_atoms):
    """polar H (AutoDock type HD) within 1.3 A of a heavy atom"""
    for a in all_atoms:
        if a["elem"] == "HD":
            if np.linalg.norm(a["pos"] - atom["pos"]) <= 1.3:
                return True
    return False


def analyze(target, ligand):
    rec = read_pdbqt_atoms(REC_PDBQT[target])
    lig = read_pdbqt_atoms(f"{RES_DIR}/{target}/{ligand}_out.pdbqt")

    O_TYPES = {"OA", "O", "OS"}
    N_TYPES = {"N", "NA", "NS"}

    # ligand donor capacity
    lig_donors = 0
    lig_acceptors = []
    for a in lig:
        if a["elem"] in O_TYPES or a["elem"] in N_TYPES:
            has_h = find_polar_h(a, lig)
            if has_h:
                lig_donors += 1
            else:
                lig_acceptors.append(a)
    # receptor donors (N/O with polar H) and acceptors (O + N w/o H)
    rec_donors = [a for a in rec if (a["elem"] in O_TYPES or a["elem"] in N_TYPES) and find_polar_h(a, rec)]
    rec_acceptors = [a for a in rec if a["elem"] in O_TYPES or (a["elem"] in N_TYPES and not find_polar_h(a, rec))]

    hbonds = []
    # receptor donor -> ligand acceptor
    for D in rec_donors:
        H = None
        for a in rec:
            if a["elem"] == "HD":
                d = np.linalg.norm(a["pos"] - D["pos"])
                if d <= 1.3:
                    H = a; break
        if H is None:
            continue
        for A in lig_acceptors:
            dDA = np.linalg.norm(D["pos"] - A["pos"])
            dHA = np.linalg.norm(H["pos"] - A["pos"])
            if dDA > 3.5 or dHA > 2.6:
                continue
            v1 = D["pos"] - H["pos"]
            v2 = A["pos"] - H["pos"]
            ang = np.degrees(np.arccos(np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9), -1, 1)))
            if ang >= 120:
                hbonds.append((D["label"], D["aname"], A["aname"], round(dDA, 2), round(dHA, 2), round(ang, 0), "rec->lig"))

    return lig_donors, hbonds


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "PTH1R"
    ligands = {
        "PTH1R": ["Ginkgolide_NPC88890", "Cucurbitacin_B", "Digoxigenin", "Kamebanin",
                   "Roemerine", "Chrysin", "Tryptanthrin", "Anisaldehyde"],
        "SLCO4C1": ["Digitoxin", "Digoxigenin"],
    }[target]
    print(f"=== H-bond analysis: {target} ===\n")
    for lig in ligands:
        try:
            nd, hb = analyze(target, lig)
        except FileNotFoundError:
            print(f"{lig}: output missing"); continue
        print(f"{lig}: ligand donors={nd}, H-bonds={len(hb)}")
        for d, dan, aan, dDA, dHA, ang, dirn in hb:
            print(f"    {d} {dan} -> lig {aan}  D..A={dDA} A  H..A={dHA} A  angle={ang:.0f}")
    print("\ndone")
