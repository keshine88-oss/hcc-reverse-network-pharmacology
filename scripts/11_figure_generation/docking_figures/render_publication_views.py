# -*- coding: utf-8 -*-
"""
render_pub.py  (PyMOL headless: pymol -cq render_pub.py)
Publication-grade rendering pass:
  main figure : cartoon (light gray) + annotated residues (teal sticks)
                + ligand (lavender sticks) + H-bond dashes (yellow)
  calib frame : separate empty scene, 6 mutually-exclusive color spheres at
                ctr +/- SEP along world XYZ, flat lighting, same view matrix
Logs to render_pub_log.txt (headless print does not reach stdout).
"""
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pymol import cmd

FIG = os.path.join(PROJECT_ROOT, "work/drug/docking/figs")
PUB = os.path.join(PROJECT_ROOT, "work/drug/docking/pub")
LOG = os.path.join(PROJECT_ROOT, "work/drug/docking/pub/render_pub_log.txt")
os.makedirs(PUB, exist_ok=True)
logf = open(LOG, "w")

def log(*a):
    logf.write(" ".join(str(x) for x in a) + "\n")
    logf.flush()

SEP = 10.0         # calibration sphere offset (A) -- large enough that spheres never occlude each other
W, H = 1800, 1350

JOBS = [
    ("6NBF_rec_H.pdb", "PTH1R_Ginkgolide_NPC88890_pose.pdb", "PTH1R_ginkgolide"),
    ("6NBF_rec_H.pdb", "PTH1R_Cucurbitacin_B_pose.pdb", "PTH1R_cucurbitacinB"),
    ("6NBF_rec_H.pdb", "PTH1R_Digoxigenin_pose.pdb", "PTH1R_digoxigenin"),
    ("AF_Q6ZQN7_rec_H.pdb", "SLCO4C1_Digitoxin_pose.pdb", "SLCO4C1_digitoxin"),
    ("AF_Q6ZQN7_rec_H.pdb", "SLCO4C1_Digoxigenin_pose.pdb", "SLCO4C1_digoxigenin"),
]

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("orthoscopic", 1)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("two_sided_lighting", 0)
cmd.set("specular", 0)
cmd.set("dash_color", [0.86, 0.80, 0.22])
cmd.set("dash_width", 3.4)
cmd.set("dash_gap", 0.26)
cmd.set("stick_radius", 0.15)

for rec, lig, name in JOBS:
    try:
        cmd.delete("all")
        cmd.load(f"{FIG}/{rec}", "rec")
        cmd.load(f"{FIG}/{lig}", "lig")

        # --- read annotated residues (chosen by prep_annotations.py) ---
        ann = json.load(open(f"{PUB}/{name}_annot.json", encoding="utf-8"))
        res_sel_parts = [f"(chain {r['chain']} and resi {r['resi']})" for r in ann["residues"]]
        res_sel = " or ".join(res_sel_parts)

        cmd.hide("everything")
        cmd.show("cartoon", "rec")
        cmd.color("gray60", "rec")
        cmd.set("cartoon_transparency", 0.25, "rec")

        # annotated residues only (hard rule: hide all other pocket sticks)
        cmd.show("sticks", "rec and (" + res_sel + ")")
        cmd.hide("sticks", "rec and not (" + res_sel + ")")
        cmd.set_color("teal_c", [0.28, 0.53, 0.50])
        cmd.color("teal_c", "rec and (" + res_sel + ") and elem C")
        cmd.color("red", "rec and (" + res_sel + ") and elem O")
        cmd.color("blue", "rec and (" + res_sel + ") and elem N")

        # ligand
        cmd.show("sticks", "lig")
        cmd.set("stick_radius", 0.24, "lig")
        cmd.set_color("lav_c", [0.66, 0.66, 0.85])
        cmd.color("lav_c", "lig and elem C")
        cmd.color("red", "lig and elem O")
        cmd.color("blue", "lig and elem N")
        # polar H only: hide all H then show H bonded to O/N
        cmd.hide("sticks", "lig and elem H")
        cmd.show("sticks", "lig and elem H and neighbor (elem O+N)")
        cmd.set_color("h_gray", [0.75, 0.75, 0.75])
        cmd.color("h_gray", "lig and elem H")

        # H-bond dashes (precomputed endpoints -> CGO-free: use distance between named atoms)
        hb = ann.get("hbonds", [])
        for i, h in enumerate(hb):
            dsel = f"rec and (chain {h['d_chain']} and resi {h['d_resi']} and name {h['d_name']})"
            asel = f"lig and index {h['a_index']}"
            cmd.distance(f"hb{i}", dsel, asel)
        cmd.hide("labels", "hb*")

        # view: ligand + one shell
        cmd.select("viewsel", "lig or (rec within 8 of lig)")
        cmd.orient("viewsel")

        cmd.ray(W, H)
        cmd.png(f"{PUB}/{name}.png", dpi=300)
        log("rendered", name)

        # ---- calibration frame in a fresh empty scene ----
        v = list(cmd.get_view())
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
        cmd.delete("all")
        cmd.read_pdbstr(pdb, "calib")
        cmd.show("spheres", "calib")
        cmd.set("sphere_scale", 0.9, "calib")
        for off, nm, color in axes:
            cmd.color(color, f"calib and name {nm}")
        # flat lighting for calibration
        cmd.set("ambient", 1.0)
        cmd.set("direct", 0.0)
        cmd.set("reflect", 0.0)
        cmd.set("specular", 0.0)
        cmd.set("ray_opaque_background", 1)
        cmd.bg_color("white")
        v2 = list(v)
        v2[15] = d - 150
        v2[16] = d + 150
        cmd.set_view(v2)
        cmd.ray(W, H)
        cmd.png(f"{PUB}/{name}_calib.png", dpi=300)
        cmd.set("ambient", 0.5)
        cmd.set("direct", 0.5)
        log("calib", name)
        # stash view for the annotator
        with open(f"{PUB}/{name}_view.json", "w") as fh:
            json.dump({"view": v, "W": W, "H": H, "SEP": SEP}, fh)
    except Exception as e:
        log("ERROR", name, repr(e))

logf.close()
