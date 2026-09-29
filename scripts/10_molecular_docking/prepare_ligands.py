# -*- coding: utf-8 -*-
"""
prepare_ligands.py
Prepare the 9 prioritized TCM active-ingredient ligands (+ optional controls)
as AutoDock Vina PDBQT files.

  SMILES -> RDKit 3D (ETKDG + MMFF) -> Meeko (rigid_macrocycles) -> PDBQT
  Post-check: scan atom types (cols 77-79) against the Vina-1.1.2 supported set.

Uses system Python 3.14 (RDKit + Meeko). No network.
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

LIG_DIR = os.path.join(PROJECT_ROOT, "docking/ligands")

# name (ASCII) -> SMILES
LIGANDS = {
    "Digitoxin": "O=C1OCC(=C1)[C@H]1CC[C@]2([C@]1(C)CC[C@H]1[C@H]2CC[C@H]2[C@]1(C)CC[C@@H](C2)O[C@H]1C[C@H](O)[C@@H]([C@H](O1)C)O[C@H]1C[C@H](O)[C@@H]([C@H](O1)C)O[C@H]1C[C@H](O)[C@@H]([C@H](O1)C)O)O",
    "Digoxigenin": "C[C@]12CC[C@@H](C[C@H]1CC[C@@H]1[C@@H]2C[C@H]([C@]2(C)[C@H](CC[C@]12O)C1=CC(=O)OC1)O)O",
    "Kamebanin": "C=C1[C@@H]2CC[C@H]3[C@@]4(C)[C@H](C[C@H]([C@@]3(C1=O)[C@@H]2O)O)C(C)(C)CC[C@@H]4O",
    "Anisaldehyde": "COc1ccc(cc1)C=O",
    "Ginkgolide_NPC88890": "O=C1O[C@@H]2[C@@]([C@@H]1C)(O)[C@@]13[C@]4([C@H]2O)[C@H](OC3=O)[C@H]([C@H](C24C(O1)OC(=O)[C@@H]2O)C(C)(C)C)O",
    "Cucurbitacin_B": "CC(=O)OC(C)(C)/C=C/C(=O)[C@@](C)([C@H]1[C@@H](C[C@@]2(C)[C@@H]3CC=C4[C@H](C[C@@H](C(=O)C4(C)C)O)[C@]3(C)C(=O)C[C@]12C)O)O",
    "Chrysin": "Oc1cc(O)c2c(c1)oc(cc2=O)c1ccccc1",
    "Roemerine": "CN1CCc2cc3c(c4-c5ccccc5CC1c24)OCO3",
    "Tryptanthrin": "c1ccc2c(c1)c(=O)n1-c3ccccc3C(=O)c1n2",
}

VINA_TYPES = {
    "C", "A", "N", "NA", "NS", "OA", "OS", "S", "SA", "H", "HD", "F",
    "Cl", "Br", "I", "P",
}


def prep_one(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[FAIL] {name}: SMILES parse error")
        return None
    mol = Chem.AddHs(mol)
    # ETKDG embed with seed retries, then MMFF
    ok = False
    for seed in (42, 1, 7, 2024):
        m = Chem.Mol(mol)
        cid = AllChem.EmbedMolecule(m, randomSeed=seed)
        if cid == 0:
            try:
                AllChem.MMFFOptimizeMolecule(m, maxIters=1000)
            except Exception:
                pass
            mol = m
            ok = True
            break
    if not ok:
        print(f"[FAIL] {name}: 3D embedding failed")
        return None

    preparator = MoleculePreparation(rigid_macrocycles=True)
    setups = preparator.prepare(mol)
    if not setups:
        print(f"[FAIL] {name}: empty setup list")
        return None
    pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        print(f"[FAIL] {name}: meeko write error: {err}")
        return None

    # atom-type precheck (cols 77-79)
    bad_types = set()
    for line in pdbqt_str.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            t = line[77:].strip()
            if t and t not in VINA_TYPES:
                bad_types.add(t)
    if bad_types:
        print(f"[WARN] {name}: non-Vina atom types {sorted(bad_types)}")

    path = os.path.join(LIG_DIR, f"{name}.pdbqt")
    with open(path, "w") as fh:
        fh.write(pdbqt_str)
    n_atoms = sum(1 for l in pdbqt_str.splitlines() if l.startswith(("ATOM", "HETATM")))
    print(f"[OK] {name}: {n_atoms} atoms -> {os.path.basename(path)}")
    return path


if __name__ == "__main__":
    os.makedirs(LIG_DIR, exist_ok=True)
    for nm, smi in LIGANDS.items():
        prep_one(nm, smi)
    print("done")
