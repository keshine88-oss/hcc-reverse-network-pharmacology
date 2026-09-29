# -*- coding: utf-8 -*-
"""
render_poses.py  (run under PyMOL headless: pymol -cq render_poses.py)
Render binding-mode figures: receptor cartoon (light gray) + pocket residues
(light green sticks) + ligand (pink sticks) + H-bond dashes (yellow).
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pymol import cmd

FIG = os.path.join(PROJECT_ROOT, "work/drug/docking/figs")
OUT = os.path.join(PROJECT_ROOT, "docking/figures")

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("cartoon_fancy_helices", 1)

JOBS = [
    # (receptor, ligand, outname, show_hb)
    ("6NBF_rec_H.pdb", "PTH1R_Ginkgolide_NPC88890_pose.pdb", "PTH1R_ginkgolide.png", True),
    ("6NBF_rec_H.pdb", "PTH1R_Cucurbitacin_B_pose.pdb", "PTH1R_cucurbitacinB.png", True),
    ("6NBF_rec_H.pdb", "PTH1R_Digoxigenin_pose.pdb", "PTH1R_digoxigenin.png", True),
    ("AF_Q6ZQN7_rec_H.pdb", "SLCO4C1_Digitoxin_pose.pdb", "SLCO4C1_Digitoxin.png", False),
    ("AF_Q6ZQN7_rec_H.pdb", "SLCO4C1_Digoxigenin_pose.pdb", "SLCO4C1_Digoxigenin.png", False),
]

for rec, lig, outname, show_hb in JOBS:
    cmd.delete("all")
    cmd.load(f"{FIG}/{rec}", "rec")
    cmd.load(f"{FIG}/{lig}", "lig")

    cmd.hide("everything")
    cmd.show("cartoon", "rec")
    cmd.color("gray78", "rec")
    cmd.set("cartoon_transparency", 0.25, "rec")

    # pocket residues within 5 A of ligand (heavy atoms only)
    cmd.select("pocket", "rec within 5 of (lig and not elem H)")
    cmd.show("sticks", "pocket")
    cmd.color("palegreen", "pocket and elem C")

    # ligand
    cmd.show("sticks", "lig")
    cmd.color("hotpink", "lig and elem C")
    cmd.color("red", "lig and elem O")
    cmd.color("blue", "lig and elem N")
    cmd.set("stick_radius", 0.22, "lig")

    # H-bonds (receptor polar H kept from Meeko)
    cmd.hide("everything", "lig and elem H")  # hide ligand H (nonpolar show only)
    cmd.show("sticks", "lig and elem H and neighbor elem O")
    cmd.hide("sticks", "lig and elem H")
    if show_hb:
        cmd.distance("hb", "lig", "rec", mode=2)
        cmd.hide("labels", "hb")
        cmd.set("dash_color", "yellow", "hb")
        cmd.set("dash_width", 3.5, "hb")

    # orient on ligand + pocket
    cmd.select("viewsel", "lig or pocket")
    cmd.orient("viewsel")

    cmd.ray(1200, 900)
    cmd.png(f"{OUT}/{outname}", dpi=200)
    print("rendered", outname)

cmd.quit()
