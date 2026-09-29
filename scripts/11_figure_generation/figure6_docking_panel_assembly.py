# -*- coding: utf-8 -*-
"""
figure6_docking_panel_assembly.py - Figure 6: molecular docking
  a   : binding-energy bar chart
  b-k : binding-site close-ups of all 10 ligands (PTH1R x8 + SLCO4C1 x2)
Close-ups are cropped from the per-ligand *_combined.png (overview|closeup),
taking the already-annotated close-up half.
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from PIL import Image, ImageDraw, ImageFont

PUB = os.path.join(PROJECT_ROOT, "work/drug/docking/pub/combined")
CHART = os.path.join(PROJECT_ROOT, "figures")
BAR = os.path.join(PROJECT_ROOT, "work/val/fig6/Figure6a_binding_energy.png")
OUT = os.path.join(PROJECT_ROOT, "figures/Figure6_docking.png")

FONT_B = "C:/Windows/Fonts/arialbd.ttf"
FONT_R = "C:/Windows/Fonts/arial.ttf"

# close-up crop box inside each combined figure
X0, Y0, X1, Y1 = 1602, 46, 3068, 1146

PANELS = [
    ("b", "PTH1R_ginkgolide",     "Ginkgolide C \u00b7 PTH1R"),
    ("c", "PTH1R_cucurbitacinB",  "Cucurbitacin B \u00b7 PTH1R"),
    ("d", "PTH1R_digoxigenin",    "Digoxigenin \u00b7 PTH1R"),
    ("e", "PTH1R_kamebanin",      "Kamebanin \u00b7 PTH1R"),
    ("f", "PTH1R_roemerine",      "Roemerine \u00b7 PTH1R"),
    ("g", "PTH1R_chrysin",        "Chrysin \u00b7 PTH1R"),
    ("h", "PTH1R_tryptanthrin",   "Tryptanthrin \u00b7 PTH1R"),
    ("i", "PTH1R_anisaldehyde",   "Anisaldehyde \u00b7 PTH1R"),
    ("j", "SLCO4C1_digitoxin",    "Digitoxin \u00b7 SLCO4C1"),
    ("k", "SLCO4C1_digoxigenin",  "Digoxigenin \u00b7 SLCO4C1"),
]

NCOL = 5
CELL_W = 560          # close-up panel width (px)
GAP = 22
MARGIN = 34
BAR_W = NCOL * CELL_W + (NCOL - 1) * GAP
LAB_H = 46            # per-panel label strip height

font_tag = ImageFont.truetype(FONT_B, 30)
font_name = ImageFont.truetype(FONT_R, 24)

# ---- load panels ----
cells = []
for tag, name, title in PANELS:
    im = Image.open(f"{PUB}/{name}_combined.png").convert("RGB").crop((X0, Y0, X1, Y1))
    ch = int(im.height * CELL_W / im.width)
    cells.append((tag, title, im.resize((CELL_W, ch), Image.LANCZOS)))

cell_h = cells[0][2].height
# ---- bar chart (resize to grid width) ----
bar = Image.open(BAR).convert("RGB")
bar_h = int(bar.height * BAR_W / bar.width)
bar = bar.resize((BAR_W, bar_h), Image.LANCZOS)

nrow = (len(cells) + NCOL - 1) // NCOL
W = MARGIN * 2 + BAR_W
H = MARGIN + bar_h + GAP + nrow * (LAB_H + cell_h + GAP) + MARGIN
canvas = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(canvas)

canvas.paste(bar, (MARGIN, MARGIN))

y = MARGIN + bar_h + GAP
for i, (tag, title, im) in enumerate(cells):
    r, c = divmod(i, NCOL)
    x = MARGIN + c * (CELL_W + GAP)
    yy = y + r * (LAB_H + cell_h + GAP)
    d.text((x + 4, yy + 8), tag, fill="black", font=font_tag)
    d.text((x + 40, yy + 14), title, fill=(60, 60, 60), font=font_name)
    canvas.paste(im, (x, yy + LAB_H))

canvas.save(OUT)
print("saved", OUT, canvas.size)
