# -*- coding: utf-8 -*-
"""
prep_figures.py
Prepare PNG figures for manuscript embedding:
  fig7: convert the tcm screening overview tiff -> PNG (width 1800 px)
  fig8: composite the 5 binding-mode PNGs into a 2x3 panel with labels a-e
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from PIL import Image, ImageDraw, ImageFont

CHART = os.path.join(PROJECT_ROOT, "figures")
FIG7_TIFF = os.path.join(PROJECT_ROOT, "results/reverse_drug_mining/tcm_screening/tcm_screening_overview.tiff")
DOCK = os.path.join(PROJECT_ROOT, "docking/figures")

FONT = "C:/Windows/Fonts/msyh.ttc"


def make_fig7():
    im = Image.open(FIG7_TIFF)
    w, h = im.size
    nw = 1800
    nh = int(h * nw / w)
    im = im.convert("RGB").resize((nw, nh), Image.LANCZOS)
    im.save(f"{CHART}/tcm_screening_overview.png")
    print(f"fig7: {im.size} -> tcm_screening_overview.png")


def make_fig8():
    panels = [
        ("a", os.path.join(PROJECT_ROOT, "work/drug/docking/pub/combined/PTH1R_ginkgolide_combined.png"),
         "PTH1R · Ginkgolide C (−8.8 kcal/mol)"),
        ("b", os.path.join(PROJECT_ROOT, "work/drug/docking/pub/combined/SLCO4C1_digitoxin_combined.png"),
         "SLCO4C1 · Digitoxin (−12.3 kcal/mol)"),
    ]
    target_w = 2800
    imgs = []
    for lab, path, title in panels:
        im = Image.open(path).convert("RGB")
        nw = target_w
        nh = int(im.height * nw / im.width)
        imgs.append((lab, im.resize((nw, nh), Image.LANCZOS), title))
    lab_h = 70
    total_h = sum(im.height + lab_h for _, im, _ in imgs) + 20
    canvas = Image.new("RGB", (target_w + 40, total_h), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font_lab = ImageFont.truetype(FONT, 40)
        font_title = ImageFont.truetype(FONT, 30)
    except Exception:
        font_lab = font_title = ImageFont.load_default()
    y = 20
    for lab, im, title in imgs:
        draw.text((24, y + 8), lab, fill="black", font=font_lab)
        draw.text((70, y + 12), title, fill="#333333", font=font_title)
        canvas.paste(im, (20, y + lab_h))
        y += lab_h + im.height
    canvas.save(f"{CHART}/Figure7_docking_binding_modes.png")
    print(f"fig7: {canvas.size} -> Figure7_docking_binding_modes.png")


if __name__ == "__main__":
    make_fig7()
    make_fig8()
