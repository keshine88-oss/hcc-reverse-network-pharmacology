# -*- coding: utf-8 -*-
"""
convert_pdbqt_to_pdb.py
Convert Meeko receptor pdbqt (keeps polar H) and Vina output pose (model 1) to PDB
for PyMOL rendering. AutoDock type (cols 77-79) -> element (cols 77-78).

  receptor: 6NBF_rec.pdbqt -> 6NBF_rec_H.pdb   (polar H kept)
             AF_Q6ZQN7_rec.pdbqt -> AF_Q6ZQN7_rec_H.pdb
  ligands:   <lig>_out.pdbqt (model 1) -> <lig>_pose.pdb
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

BASE = os.path.join(PROJECT_ROOT, "work/drug/docking")
REC_DIR = f"{BASE}/receptor"
RES_DIR = f"{BASE}/results"
FIG_DIR = f"{BASE}/figs"

AD2ELEM = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
    "OA": "O", "OS": "O", "S": "S", "SA": "S", "H": "H", "HD": "H",
    "F": "F", "Cl": "CL", "Br": "BR", "I": "I", "P": "P",
}


def pdbqt_atom_to_pdb(line):
    if len(line) < 54:
        return None
    prefix = line[:66]  # PDB-compatible up to temperature factor (cols 1-66)
    # ensure element columns 77-78 (1-indexed) are filled
    ad = line[77:].strip() if len(line) >= 78 else ""
    elem = AD2ELEM.get(ad, "")
    if not elem:
        # infer from atom name
        aname = line[12:16].strip()
        elem = aname[0] if aname and aname[0].isalpha() else "C"
    # rebuild: prefix (66) + 10 spaces + element right-justified 2 cols
    out = prefix + " " * 10 + elem.rjust(2)
    return out


def convert_receptor(rec_pdbqt, out_pdb):
    with open(rec_pdbqt) as fh, open(out_pdb, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                p = pdbqt_atom_to_pdb(line)
                if p:
                    out.write(p + "\n")
            elif line.startswith(("END", "TER")):
                out.write(line)
    print(f"receptor -> {out_pdb}")


def convert_pose(lig_pdbqt, out_pdb):
    model = 0
    lines = []
    with open(lig_pdbqt) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                model += 1
                if model > 1:
                    break
                continue
            if line.startswith("ENDMDL"):
                continue
            if line.startswith(("ATOM", "HETATM")):
                p = pdbqt_atom_to_pdb(line)
                if p:
                    lines.append(p)
    with open(out_pdb, "w") as out:
        out.write("HETATM" + lines[0][6:] + "\n" if lines else "")
        for p in lines:
            out.write(p + "\n")
        out.write("END\n")
    print(f"pose -> {out_pdb}")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    convert_receptor(f"{REC_DIR}/6NBF_rec.pdbqt", f"{FIG_DIR}/6NBF_rec_H.pdb")
    convert_receptor(f"{REC_DIR}/AF_Q6ZQN7_rec.pdbqt", f"{FIG_DIR}/AF_Q6ZQN7_rec_H.pdb")

    poses = {
        "PTH1R": ["Ginkgolide_NPC88890", "Cucurbitacin_B", "Digoxigenin", "Kamebanin",
                   "Roemerine", "Chrysin", "Tryptanthrin", "Anisaldehyde"],
        "SLCO4C1": ["Digitoxin", "Digoxigenin"],
    }
    for target, ligs in poses.items():
        for lig in ligs:
            convert_pose(f"{RES_DIR}/{target}/{lig}_out.pdbqt", f"{FIG_DIR}/{target}_{lig}_pose.pdb")


if __name__ == "__main__":
    main()
