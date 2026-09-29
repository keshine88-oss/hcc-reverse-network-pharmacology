# -*- coding: utf-8 -*-
"""
prep_annotations.py
Choose annotated residues per job (H-bond partners union hydrophobic contacts,
<=6 total) and export geometry (CA, outer atom, H-bond endpoints, ligand centroid)
as {name}_annot.json for render_pub.py / annotate_pub.py.

Inputs: figs/{rec}_H.pdb (polar H kept), figs/{target}_{lig}_pose.pdb (model-1 pose)
"""
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import gemmi

FIG = os.path.join(PROJECT_ROOT, "work/drug/docking/figs")
PUB = os.path.join(PROJECT_ROOT, "work/drug/docking/pub")
os.makedirs(PUB, exist_ok=True)

AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}
O_SET = {"OA", "O", "OS"}
N_SET = {"N", "NA", "NS"}
MAX_RES = 6


def read_pdb_atoms(path):
    """parse PDB; returns list of dicts (pos, elem, aname, resname, chain, resi, serial)"""
    atoms = []
    model = 0
    for line in open(path, encoding="ascii", errors="replace"):
        if line.startswith("MODEL"):
            model += 1
            if model > 1:
                break
            continue
        if line.startswith("ENDMDL"):
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            pos = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        atoms.append({
            "pos": pos,
            "elem": line[76:78].strip() or line[12:16].strip()[0],
            "aname": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "chain": line[21:22].strip(),
            "resi": line[22:26].strip(),
            "serial": int(line[6:11]),
        })
    return atoms


def find_polar_h(atom, all_atoms):
    return any(a["elem"] == "H" and np.linalg.norm(a["pos"] - atom["pos"]) <= 1.3
               for a in all_atoms)


def mass_center(atoms):
    sw = np.zeros(3)
    tw = 0.0
    for a in atoms:
        w = gemmi.Element(a["elem"]).weight
        sw += w * a["pos"]
        tw += w
    return sw / tw


def analyze(rec, lig):
    # ligand donors / acceptors
    lig_acc = []
    lig_don = 0
    for a in lig:
        if a["elem"] in O_SET | N_SET:
            if find_polar_h(a, lig):
                lig_don += 1
            else:
                lig_acc.append(a)
    rec_don = [a for a in rec if (a["elem"] in O_SET | N_SET) and find_polar_h(a, rec)]

    hbonds = []
    hb_res = {}
    for D in rec_don:
        H = next((a for a in rec if a["elem"] == "H"
                  and np.linalg.norm(a["pos"] - D["pos"]) <= 1.3), None)
        if H is None:
            continue
        for A in lig_acc:
            dDA = float(np.linalg.norm(D["pos"] - A["pos"]))
            dHA = float(np.linalg.norm(H["pos"] - A["pos"]))
            if dDA > 3.5 or dHA > 2.6:
                continue
            v1, v2 = D["pos"] - H["pos"], A["pos"] - H["pos"]
            ang = float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) /
                     (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9), -1, 1))))
            if ang >= 120:
                key = (D["chain"], D["resi"], D["aname"], A["serial"])
                if key not in {h["key"] for h in hbonds}:
                    hbonds.append({"key": key, "D": D, "A": A, "dist": dDA})
                    hb_res[(D["chain"], D["resi"])] = D["resname"]

    # hydrophobic contacts: rec carbon within 4.5 A of lig carbon
    lig_c = [a for a in lig if a["elem"] in ("C", "A")]
    rec_c = [a for a in rec if a["elem"] in ("C", "A")]
    contacts = {}
    for rc in rec_c:
        for lc in lig_c:
            if np.linalg.norm(rc["pos"] - lc["pos"]) <= 4.5:
                k = (rc["chain"], rc["resi"])
                contacts[k] = contacts.get(k, 0) + 1
                break
    hydro = sorted(contacts.items(), key=lambda kv: -kv[1])
    return hbonds, hb_res, hydro, lig_don


def main():
    jobs = [
        ("6NBF_rec_H.pdb", "PTH1R_Ginkgolide_NPC88890_pose.pdb", "PTH1R_ginkgolide"),
        ("6NBF_rec_H.pdb", "PTH1R_Cucurbitacin_B_pose.pdb", "PTH1R_cucurbitacinB"),
        ("6NBF_rec_H.pdb", "PTH1R_Digoxigenin_pose.pdb", "PTH1R_digoxigenin"),
        ("6NBF_rec_H.pdb", "PTH1R_Kamebanin_pose.pdb", "PTH1R_kamebanin"),
        ("6NBF_rec_H.pdb", "PTH1R_Roemerine_pose.pdb", "PTH1R_roemerine"),
        ("6NBF_rec_H.pdb", "PTH1R_Chrysin_pose.pdb", "PTH1R_chrysin"),
        ("6NBF_rec_H.pdb", "PTH1R_Tryptanthrin_pose.pdb", "PTH1R_tryptanthrin"),
        ("6NBF_rec_H.pdb", "PTH1R_Anisaldehyde_pose.pdb", "PTH1R_anisaldehyde"),
        ("AF_Q6ZQN7_rec_H.pdb", "SLCO4C1_Digitoxin_pose.pdb", "SLCO4C1_digitoxin"),
        ("AF_Q6ZQN7_rec_H.pdb", "SLCO4C1_Digoxigenin_pose.pdb", "SLCO4C1_digoxigenin"),
    ]
    for rec_f, lig_f, name in jobs:
        rec = read_pdb_atoms(f"{FIG}/{rec_f}")
        lig = read_pdb_atoms(f"{FIG}/{lig_f}")
        centroid = mass_center([a for a in lig if a["elem"] != "H"])

        hbonds, hb_res, hydro, lig_don = analyze(rec, lig)

        # residue map from receptor atoms
        def res_atoms(key):
            return [a for a in rec if (a["chain"], a["resi"]) == key and a["elem"] != "H"]

        residues = []
        used = set()
        # H-bond partners first
        for (ch, ri), rn in hb_res.items():
            used.add((ch, ri))
            atoms = res_atoms((ch, ri))
            ca = next((a for a in atoms if a["aname"] == "CA"), atoms[0])
            outer = max(atoms, key=lambda a: np.linalg.norm(a["pos"] - centroid))
            residues.append({
                "chain": ch, "resi": ri, "resname": rn,
                "label": f"{rn}-{ri.strip()}",
                "ca": list(np.round(ca["pos"], 3)),
                "outer": list(np.round(outer["pos"], 3)),
                "hbond": True,
            })
        # then hydrophobic contacts (exclude hbond residues), fill up to MAX_RES
        for (ch, ri), cnt in hydro:
            if len(residues) >= MAX_RES:
                break
            if (ch, ri) in used:
                continue
            atoms = res_atoms((ch, ri))
            if not atoms:
                continue
            rn = atoms[0]["resname"]
            ca = next((a for a in atoms if a["aname"] == "CA"), atoms[0])
            outer = max(atoms, key=lambda a: np.linalg.norm(a["pos"] - centroid))
            residues.append({
                "chain": ch, "resi": ri, "resname": rn,
                "label": f"{rn}-{ri.strip()}",
                "ca": list(np.round(ca["pos"], 3)),
                "outer": list(np.round(outer["pos"], 3)),
                "hbond": False,
            })
            used.add((ch, ri))

        hb_out = []
        for h in hbonds:
            D, A = h["D"], h["A"]
            hb_out.append({
                "d_chain": D["chain"], "d_resi": D["resi"], "d_name": D["aname"],
                "a_index": A["serial"],
                "dist": round(h["dist"], 2),
                "d_pos": list(np.round(D["pos"], 3)),
                "a_pos": list(np.round(A["pos"], 3)),
                "h_pos": None,
            })
        out = {
            "residues": residues,
            "hbonds": hb_out,
            "lig_donors": lig_don,
            "lig_centroid": list(np.round(centroid, 3)),
            "lig_atoms": [list(np.round(a["pos"], 3)) for a in lig if a["elem"] != "H"],
        }
        with open(f"{PUB}/{name}_annot.json", "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1, ensure_ascii=False)
        print(f"{name}: residues={[r['label'] for r in residues]} hbonds={len(hb_out)} lig_donors={lig_don}")


if __name__ == "__main__":
    main()
