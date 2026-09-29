# -*- coding: utf-8 -*-
"""
render_combined.py  (PyMOL headless)
Per ligand, render:
  overview.png     whole receptor cartoon (bright cyan) + ligand sticks (purple),
                   orthoscopic, full view          + calibration frame
  closeup.png      pocket close-up: semi-transparent cyan cartoon (ligand never
                   occluded), annotated residues teal, ligand mauve, H-bond
                   yellow dashes                   + calibration frame
Colors follow the user's reference figure.
"""
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pymol import cmd

FIG = os.path.join(PROJECT_ROOT, "work/drug/docking/figs")
PUB = os.path.join(PROJECT_ROOT, "work/drug/docking/pub")
LOG = f"{PUB}/render_combined_log.txt"
os.makedirs(PUB, exist_ok=True)
logf = open(LOG, "w")


def log(*a):
    logf.write(" ".join(str(x) for x in a) + "\n")
    logf.flush()


SEP = 10.0
OV_W, OV_H = 1400, 1050
CL_W, CL_H = 1800, 1350

# reference-figure palette
CYAN = [0.30, 0.80, 0.80]
MAUVE = [0.71, 0.50, 0.75]
TEAL = [0.28, 0.53, 0.50]
YELLOW = [0.86, 0.80, 0.22]
RED = [0.80, 0.16, 0.16]
BLUE = [0.12, 0.10, 0.62]
HGRAY = [0.72, 0.72, 0.72]

JOBS = [
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

cmd.set("orthoscopic", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("ray_opaque_background", 1)
cmd.set("specular", 0)
cmd.bg_color("white")


def calib_frame(v, W, H, tag):
    """empty scene, 6 color spheres, flat light, same view; returns None"""
    cmd.delete("all")
    ctr = v[12:15]
    d = -v[11]
    axes = [((SEP, 0, 0), "S1", "red"), ((-SEP, 0, 0), "S2", "cyan"),
            ((0, SEP, 0), "S3", "green"), ((0, -SEP, 0), "S4", "magenta"),
            ((0, 0, SEP), "S5", "blue"), ((0, 0, -SEP), "S6", "yellow")]
    rows = []
    for i, (off, nm, _c) in enumerate(axes):
        x, y, z = ctr[0] + off[0], ctr[1] + off[1], ctr[2] + off[2]
        rows.append(f"HETATM{i+1:5d}  {nm:<4}LIG B   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C")
    pdb = "\n".join(rows) + "\nEND\n"
    cmd.read_pdbstr(pdb, "calib")
    cmd.hide("everything")
    cmd.show("spheres", "calib")
    cmd.set("sphere_scale", 0.9, "calib")
    for off, nm, color in axes:
        cmd.color(color, f"calib and name {nm}")
    cmd.set("ambient", 1.0)
    cmd.set("direct", 0.0)
    cmd.set("reflect", 0.0)
    cmd.set("specular", 0.0)
    v2 = list(v)
    v2[15] = d - 150
    v2[16] = d + 150
    cmd.set_view(v2)
    cmd.ray(W, H)
    cmd.png(f"{PUB}/{tag}.png", dpi=300)
    cmd.set("ambient", 0.5)
    cmd.set("direct", 0.5)
    log("calib", tag)


def style_ligand():
    cmd.show("sticks", "lig")
    cmd.set("stick_radius", 0.24, "lig")
    cmd.set_color("mauve_c", MAUVE)
    cmd.color("mauve_c", "lig and elem C")
    cmd.color("red", "lig and elem O")
    cmd.color("blue", "lig and elem N")
    cmd.hide("sticks", "lig and elem H")
    cmd.show("sticks", "lig and elem H and neighbor (elem O+N)")
    cmd.set_color("h_gray", HGRAY)
    cmd.color("h_gray", "lig and elem H")


for rec, lig, name in JOBS:
    try:
        # ---------- overview ----------
        cmd.delete("all")
        cmd.load(f"{FIG}/{rec}", "rec")
        cmd.load(f"{FIG}/{lig}", "lig")
        cmd.hide("everything")
        cmd.show("cartoon", "rec")
        cmd.set_color("cyan_b", CYAN)
        cmd.color("cyan_b", "rec")
        cmd.set("cartoon_transparency", 0.0, "rec")
        style_ligand()
        cmd.set("stick_radius", 0.20, "lig")
        cmd.orient("rec")
        cmd.ray(OV_W, OV_H)
        cmd.png(f"{PUB}/{name}_overview.png", dpi=300)
        log("overview", name)
        v_ov = list(cmd.get_view())
        calib_frame(v_ov, OV_W, OV_H, f"{name}_overview_calib")

        # ---------- closeup ----------
        cmd.delete("all")
        cmd.load(f"{FIG}/{rec}", "rec")
        cmd.load(f"{FIG}/{lig}", "lig")
        cmd.hide("everything")
        cmd.show("cartoon", "rec")
        cmd.color("cyan_b", "rec")
        cmd.set("cartoon_transparency", 0.55, "rec")  # never occlude the ligand

        ann = json.load(open(f"{PUB}/{name}_annot.json", encoding="utf-8"))
        res_sel = " or ".join(f"(chain {r['chain']} and resi {r['resi']})" for r in ann["residues"])
        cmd.show("sticks", "rec and (" + res_sel + ")")
        cmd.hide("sticks", "rec and not (" + res_sel + ")")
        cmd.set_color("teal_c", TEAL)
        cmd.color("teal_c", "rec and (" + res_sel + ") and elem C")
        cmd.color("red", "rec and (" + res_sel + ") and elem O")
        cmd.color("blue", "rec and (" + res_sel + ") and elem N")

        style_ligand()

        # H-bond dashes are drawn later in compose_combined.py (Pillow) from the
        # same projected coordinates as the distance labels, so dash and label
        # are guaranteed consistent (PyMOL `index` != PDB `serial` here because
        # the pose PDB repeats its first atom).
        cmd.hide("everything", "hb*")

        cmd.select("viewsel", "lig or (rec within 8 of lig)")
        cmd.orient("viewsel")
        cmd.ray(CL_W, CL_H)
        cmd.png(f"{PUB}/{name}_closeup.png", dpi=300)
        log("closeup", name)
        v_cl = list(cmd.get_view())
        calib_frame(v_cl, CL_W, CL_H, f"{name}_closeup_calib")

        with open(f"{PUB}/{name}_views.json", "w") as fh:
            json.dump({"overview": v_ov, "closeup": v_cl,
                       "ov_size": [OV_W, OV_H], "cl_size": [CL_W, CL_H], "SEP": SEP}, fh)
    except Exception as e:
        log("ERROR", name, repr(e))

logf.close()
