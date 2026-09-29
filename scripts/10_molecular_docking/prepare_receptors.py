# -*- coding: utf-8 -*-
"""
prepare_receptors.py
Prepare receptor structures for AutoDock Vina docking.

  - 6NBF (PTH1R + LA-PTH + Gs + Nb35): keep receptor chain R and ligand peptide
    chain P; drop G-protein (A/B/G) and nanobody (N); drop lipids (CLR/PLM) and
    waters; write 6NBF_receptor.pdb (chain R, standard residues only).
    Box center = centroid of the LA-PTH N-terminal residues (the portion that
    inserts into the transmembrane orthosteric pocket), which defines the
    small-molecule orthosteric site.
  - AF-Q6ZQN7 (SLCO4C1 AlphaFold v6): report per-residue pLDDT distribution and
    write the receptor PDB (single chain, all standard residues).

Uses gemmi (system Python 3.14). No network. Pure stdlib + gemmi.
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json
import gemmi

REC_DIR = os.path.join(PROJECT_ROOT, "docking/receptors")

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # modified residues treat as standard-like: keep (they map to common AA)
    "MSE", "SEC", "PYL",
}


def is_standard_aa(res):
    return res.name in STD_AA


def centroid_of_residues(residues):
    """mass-weighted centroid (element weight) of a list of gemmi.Residue."""
    sx = sy = sz = w = 0.0
    for r in residues:
        for a in r:
            wt = gemmi.Element(a.element.name).weight
            p = a.pos
            sx += wt * p.x
            sy += wt * p.y
            sz += wt * p.z
            w += wt
    if w == 0:
        return None
    return (sx / w, sy / w, sz / w)


def process_6nbf():
    st = gemmi.read_structure(f"{REC_DIR}/6NBF.pdb")
    model = st[0]

    # identify chains
    chains = {ch.name: ch for ch in model}

    # --- ligand peptide N-terminal residues (chain P, first 12 residues) ---
    chain_p = chains.get("P")
    nterm = [r for r in chain_p if r.seqid.num <= 12 and is_standard_aa(r)]
    c = centroid_of_residues(nterm)
    print(f"LA-PTH N-terminal (resi<=12) centroid: ({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f})")
    # full peptide centroid for reference
    c_full = centroid_of_residues([r for r in chain_p if is_standard_aa(r)])
    print(f"LA-PTH full peptide centroid: ({c_full[0]:.3f}, {c_full[1]:.3f}, {c_full[2]:.3f})")

    # --- build receptor-only structure (chain R, standard residues) ---
    rec = gemmi.Structure()
    rec.name = "6NBF_PTH1R_receptor"
    rec.cell = st.cell
    rec.spacegroup_hm = st.spacegroup_hm
    m = gemmi.Model("1")
    ch = gemmi.Chain("R")
    n_kept = 0
    for r in chains["R"]:
        if is_standard_aa(r):
            # deep-copy the residue
            import copy
            ch.add_residue(copy.deepcopy(r))
            n_kept += 1
    m.add_chain(ch)
    rec.add_model(m)
    rec.write_pdb(f"{REC_DIR}/6NBF_receptor.pdb")
    print(f"6NBF_receptor.pdb written: chain R {n_kept} standard residues")

    # ligand peptide reference (full, standard residues) as separate PDB for box definition
    lig = gemmi.Structure()
    lig.name = "LA-PTH_ligand"
    lig.cell = st.cell
    lig.spacegroup_hm = st.spacegroup_hm
    m2 = gemmi.Model("1")
    ch2 = gemmi.Chain("P")
    for r in chain_p:
        if is_standard_aa(r):
            import copy
            ch2.add_residue(copy.deepcopy(r))
    m2.add_chain(ch2)
    lig.add_model(m2)
    lig.write_pdb(f"{REC_DIR}/6NBF_ligand_LA-PTH.pdb")

    site = {
        "protein": "PTH1R",
        "pdb": "6NBF",
        "box_center": list(c),
        "box_center_note": "centroid of LA-PTH N-terminal (resi<=12) in transmembrane orthosteric pocket",
        "box_center_full_peptide": list(c_full),
    }
    with open(f"{REC_DIR}/6NBF_site.json", "w", encoding="utf-8") as fh:
        json.dump(site, fh, indent=2, ensure_ascii=False)
    print("6NBF_site.json written")


def process_alphafold():
    st = gemmi.read_structure(f"{REC_DIR}/AF-Q6ZQN7.pdb")
    model = st[0]
    # pLDDT is stored in B-factor
    bf = []
    for ch in model:
        for r in ch:
            for a in r:
                bf.append(a.b_iso)
    n = len(bf)
    mean = sum(bf) / n
    very_low = sum(1 for b in bf if b < 50)
    low = sum(1 for b in bf if 50 <= b < 70)
    conf = sum(1 for b in bf if 70 <= b < 90)
    very_high = sum(1 for b in bf if b >= 90)
    print(f"AF-Q6ZQN7 pLDDT: n={n}, mean={mean:.1f}, "
          f"<50={very_low}({100*very_low/n:.1f}%) 50-70={low}({100*low/n:.1f}%) "
          f"70-90={conf}({100*conf/n:.1f}%) >=90={very_high}({100*very_high/n:.1f}%)")

    # write a receptor-only PDB (standard residues; AlphaFold may include non-std)
    rec = gemmi.Structure()
    rec.name = "AF_SLCO4C1"
    rec.cell = st.cell
    rec.spacegroup_hm = st.spacegroup_hm
    m = gemmi.Model("1")
    ch = gemmi.Chain("A")
    n_kept = 0
    for src in model:
        for r in src:
            if is_standard_aa(r):
                import copy
                ch.add_residue(copy.deepcopy(r))
                n_kept += 1
    m.add_chain(ch)
    rec.add_model(m)
    rec.write_pdb(f"{REC_DIR}/AF-Q6ZQN7_receptor.pdb")
    print(f"AF-Q6ZQN7_receptor.pdb written: {n_kept} residues")


if __name__ == "__main__":
    process_6nbf()
    print()
    process_alphafold()
