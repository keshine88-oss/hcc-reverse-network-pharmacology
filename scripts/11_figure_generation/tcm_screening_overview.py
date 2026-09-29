# -*- coding: utf-8 -*-
"""
plot_tcm_figure.py
TCM active-ingredient screening summary figure (quantitative grid, 3 panels).

  a | screening funnel: target-page associations -> quantitative records ->
      curated-TCM-annotated compounds -> covered herbs
  b | hero: 14 most potent nameable PTH1R natural products (log10 nM),
      red = sub-micromolar; star = dual-target (also SLCO4C1)
  c | SLCO4C1 hits with ChEMBL cross-validation note

QA notes (exclusion rule): panel b shows the 14 most potent PTH1R hits that
have a resolvable common or class name; three unnamed InChIKey-only entries
(NPC602666 766.7 nM; NPC599838 / NPC606152 4466.8 nM) are excluded from the
panel but retained in npass_tcm_hits.csv. NPC88890 is shown under its class
name (ginkgolide-type cage lactone).
"""
import csv
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.insert(0, r"C:\Users\ZKX\.workbuddy\skills\nature-figure\scripts")
from audit_panel_alignment import require_matplotlib_panel_alignment

DATA = os.path.join(PROJECT_ROOT, "work/drug/npass_tcm_hits.csv")
PIVOT = os.path.join(PROJECT_ROOT, "work/drug/npass_tcm_herb_pivot.csv")
OUT = os.path.join(PROJECT_ROOT, "results/reverse_drug_mining/tcm_screening/tcm_screening_overview")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "Arial", "SimHei", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.unicode_minus": False,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

STEEL = "#6E8CA8"      # neutral signal family
RED = "#C0504D"        # accent: sub-micromolar
GRAYS = ["#A9BCCB", "#88A2B6", "#67869F", "#476A80"]  # sequential funnel

rows = list(csv.DictReader(open(DATA, encoding="utf-8-sig")))
pivot = list(csv.DictReader(open(PIVOT, encoding="utf-8-sig")))

# ---- funnel counts (computed from data, not hard-coded) ----
uniq = {r["npc_id"] for r in rows}
quant = {r["npc_id"] for r in rows if r["potency_nM"]}
herb = {r["npc_id"] for r in rows if r["tcm_herbs"]}
n_herbs = len(pivot)
funnel = [
    ("Natural products linked on the\nNPASS target page", len(uniq), "compounds"),
    ("with quantitative activity records", len(quant), "compounds"),
    ("curated with TCM source annotations", len(herb), "compounds"),
    ("covering TCM / medicinal sources", n_herbs, "medicinal materials"),
]

# ---- panel b: top-14 nameable PTH1R hits ----
SKIP_UNNAMED = {"NPC602666", "NPC599838", "NPC606152"}
DISPLAY = {"NPC88890": "unidentified ginkgolide-type compound*"}
pth = []
for r in rows:
    if r["gene"] != "PTH1R" or not r["potency_nM"] or r["npc_id"] in SKIP_UNNAMED:
        continue
    pth.append(r)
pth.sort(key=lambda r: float(r["potency_nM"]))
pth = pth[:14]

def disp_name(r):
    n = DISPLAY.get(r["npc_id"], r["compound_name"])
    if r["npc_id"] == "NPC189588":
        n += " ★"
    return n

WRAP = {
    "7-Methoxy-3-(4-Methoxyphenyl)Chromen-4-One": "7-Methoxy-3-(4-Methoxyphenyl)\nChromen-4-One",
    "2-Anilinonaphthalene-1,4-Dione": "2-Anilinonaphthalene-\n1,4-Dione",
}

def disp_herb(r):
    h = r["tcm_herbs"].split(";")[0].strip() if r["tcm_herbs"] else "—"
    h = h.replace(" (microbial source)", "")
    return h.split("（")[0]  # trim inner parenthetical for display

# ---- panel c: SLCO4C1 ----
slc = [r for r in rows if r["gene"] == "SLCO4C1"]
slc_q = sorted([r for r in slc if r["potency_nM"]], key=lambda r: float(r["potency_nM"]))
slc_assoc = [r for r in slc if not r["potency_nM"]]

MM = 1 / 25.4
fig = plt.figure(figsize=(180 * MM, 80 * MM))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.75, 0.95], wspace=0.70,
                      left=0.085, right=0.985, top=0.84, bottom=0.175)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[0, 2])

for ax, lab in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
    ax.text(-0.02, 1.10, lab, transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", va="top", ha="left")

# ---------- panel a: funnel ----------
y = list(range(len(funnel)))[::-1]
maxv = max(v for _, v, _ in funnel)
for (lab, v, unit), yi, c in zip(funnel, y, GRAYS):
    ax_a.barh(yi, v, height=0.62, color=c, edgecolor="none")
    ax_a.text(v * 0.04, yi, f"{v} {unit}", va="center", fontsize=6.5,
              fontweight="bold", color="white")
    ax_a.text(-maxv * 0.03, yi, lab, va="center", ha="right", fontsize=6.5)
ax_a.set_xlim(0, maxv * 1.05)
ax_a.set_ylim(-0.6, len(funnel) - 0.4)
ax_a.axis("off")
ax_a.set_title("Screening funnel (unique compounds)", fontsize=7.5, pad=2)

# ---------- panel b: hero potency bars ----------
names = [disp_name(r) for r in pth]
herbs = [disp_herb(r) for r in pth]
vals = [float(r["potency_nM"]) for r in pth]
cols = [RED if v < 1000 else STEEL for v in vals]
ylab = [WRAP.get(n, n) for n in names]

yi = list(range(len(pth)))[::-1]
ax_b.barh(yi, vals, height=0.66, color=cols, edgecolor="none")
for v, y0, h in zip(vals, yi, herbs):
    txt = f"{v:,.0f}" if v >= 1000 else f"{v:,.1f}"
    ax_b.text(v * 1.12, y0, f"{txt}｜{h}", va="center", fontsize=6)
ax_b.axvline(1000, color="#888888", lw=0.7, ls=(0, (3, 2)))
ax_b.text(1000, 13.95, "1 μM", fontsize=6, color="#555555", ha="center", va="bottom")
ax_b.set_xscale("log")
ax_b.set_xlim(500, 60000)
ax_b.set_ylim(-0.6, 14.6)
ax_b.set_yticks(yi)
ax_b.set_yticklabels(ylab, fontsize=6.5)
for tick, r in zip(ax_b.get_yticklabels(), pth):
    if r["npc_id"] == "NPC189588":
        tick.set_fontweight("bold")
ax_b.set_xlabel("NPASS potency (nM, log scale)", fontsize=7)
ax_b.set_title("High-potency natural products for PTH1R (top 14)", fontsize=7.5, pad=2)
ax_b.tick_params(axis="x", labelsize=6.5)
ax_b.set_xticks([1000, 10000])
ax_b.set_xticklabels(["1,000", "10,000"])

# ---------- panel c: SLCO4C1 ----------
names_c = [r["compound_name"] + (" ★" if r["npc_id"] == "NPC189588" else "") for r in slc_q]
vals_c = [float(r["potency_nM"]) for r in slc_q]
yi_c = list(range(len(slc_q)))[::-1]
ax_c.barh(yi_c, vals_c, height=0.5, color=[RED] * len(vals_c), edgecolor="none")
for v, y0, r in zip(vals_c, yi_c, slc_q):
    ax_c.text(v * 1.15, y0, f"IC50 = {v:,.0f} nM", va="center", fontsize=6.5)
ax_c.set_xscale("log")
ax_c.set_xlim(60, 40000)
ax_c.set_yticks(yi_c)
ax_c.set_yticklabels(names_c, fontsize=6.5)
ax_c.get_yticklabels()[0].set_fontweight("normal")
for tick, r in zip(ax_c.get_yticklabels(), slc_q):
    if r["npc_id"] == "NPC189588":
        tick.set_fontweight("bold")
extra = "；".join(r["compound_name"] for r in slc_assoc)
ax_c.text(70, -0.75, f"Link-only hits (no quantitative value): {extra}\nDigitoxin cross-validated with ChEMBL\n(pChEMBL 6.92 ~ 120 nM)",
          fontsize=6, va="top", color="#333333")
ax_c.set_xlabel("SLCO4C1 inhibitory activity (nM, log scale)", fontsize=7)
ax_c.set_title("SLCO4C1 hits (all cardiac glycosides / hormones)", fontsize=7.5, pad=2)
ax_c.tick_params(axis="x", labelsize=6.5)
ax_c.set_xticks([100, 1000, 10000])
ax_c.set_xticklabels(["100", "1,000", "10,000"])
ax_c.set_ylim(-2.6, len(slc_q) - 0.4)

fig.text(0.5, 0.965, "TCM active-compound screening for the core targets (NPASS 3.0)", ha="center",
         fontsize=8.5, fontweight="bold")
fig.text(0.985, 0.02, "* hit both targets (PTH1R + SLCO4C1)    * no common name annotated in NPASS, grouped by scaffold",
         ha="right", fontsize=6, color="#444444")

# ---------- alignment gate + export ----------
require_matplotlib_panel_alignment(
    fig,
    axes=[ax_a, ax_b, ax_c],
    panel_ids=["a", "b", "c"],
    json_out=OUT + ".alignment.json",
    overlay_svg=OUT + ".alignment.svg",
    tolerance_pt=1.5,
    gutter_tolerance_pt=1.5,
    strict=True,
    exemptions=[{
        "panels": ["a", "b", "c"],
        "checks": ["panel-width"],
        "reason": "intentional hero design: panel b is wider to host potency labels and log axis",
    }],
)
fig.savefig(OUT + ".svg", bbox_inches="tight")
fig.savefig(OUT + ".pdf", bbox_inches="tight")
fig.savefig(OUT + ".tiff", dpi=600, bbox_inches="tight")
print("saved:", OUT + ".{svg,pdf,tiff}")
print(f"funnel: assoc={len(uniq)} quant={len(quant)} herb={len(herb)} herbs={n_herbs}")
