# -*- coding: utf-8 -*-
"""
find_slco4c1_cavity.py
Locate the central substrate cavity of SLCO4C1 (OATP transporter) in the
AlphaFold model.

Approach for a membrane transporter:
  1. keep high-pLDDT (>=70) atoms = ordered TM-helix core
  2. PCA -> first principal axis ~ membrane normal (transmembrane bundle long axis)
  3. take the middle 50% slab along that axis (mid-membrane region)
  4. inside that slab, distance-transform -> deepest empty point = central cavity
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json
import numpy as np
import gemmi
from scipy import ndimage

REC_DIR = os.path.join(PROJECT_ROOT, "docking/receptors")
SPACING = 1.0


def main():
    st = gemmi.read_structure(f"{REC_DIR}/AF-Q6ZQN7_receptor.pdb")
    coords = []
    for ch in st[0]:
        for r in ch:
            for a in r:
                if a.element.name != "H" and a.b_iso >= 70.0:
                    coords.append((a.pos.x, a.pos.y, a.pos.z))
    coords = np.array(coords)
    print(f"high-pLDDT atoms: {len(coords)}")

    # PCA -> membrane normal
    center = coords.mean(0)
    c = coords - center
    u, s, vh = np.linalg.svd(c, full_matrices=False)
    axis = vh[0]  # first principal axis (long axis)
    print(f"principal axes stddev: {s[0]:.1f}, {s[1]:.1f}, {s[2]:.1f} A")

    proj = c @ axis
    lo_p, hi_p = np.percentile(proj, [25, 75])
    mid = coords[(proj >= lo_p) & (proj <= hi_p)]
    print(f"mid-slab atoms: {len(mid)} (proj in [{lo_p:.0f}, {hi_p:.0f}] A)")

    # bounding box of mid slab
    pad = 4.0
    lo = mid.min(0) - pad
    hi = mid.max(0) + pad
    shape = tuple(np.ceil((hi - lo) / SPACING).astype(int) + 1)

    occ = np.zeros(shape, dtype=bool)
    idx = np.clip(np.rint((mid - lo) / SPACING).astype(int), 0, np.array(shape) - 1)
    for i in idx:
        occ[tuple(i)] = True

    dt = ndimage.distance_transform_edt(~occ, sampling=SPACING)
    deepest = np.unravel_index(np.argmax(dt), dt.shape)
    dp = np.array([deepest[2], deepest[1], deepest[0]]) * SPACING + lo
    print(f"central-cavity deepest point (A): ({dp[0]:.2f}, {dp[1]:.2f}, {dp[2]:.2f}), dt={dt[deepest]:.2f} A")

    # largest enclosed cavity within slab (>=1.5 A from protein)
    cavity = (~occ) & (dt >= 1.5)
    lbl, n = ndimage.label(cavity)
    sizes = ndimage.sum(cavity, lbl, range(1, n + 1))
    biggest = int(np.argmax(sizes)) + 1
    mask = lbl == biggest
    zz, yy, xx = np.nonzero(mask)
    cg = np.array([xx.mean(), yy.mean(), zz.mean()]) * SPACING + lo
    vol = int(sizes[biggest - 1]) * SPACING**3
    print(f"largest enclosed cavity centroid: ({cg[0]:.2f}, {cg[1]:.2f}, {cg[2]:.2f}), vol={vol:.0f} A^3")

    with open(f"{REC_DIR}/AF-Q6ZQN7_site.json", "w", encoding="utf-8") as fh:
        json.dump({
            "protein": "SLCO4C1",
            "model": "AlphaFold Q6ZQN7 v6",
            "box_center": [round(float(cg[0]), 3), round(float(cg[1]), 3), round(float(cg[2]), 3)],
            "cavity_volume_A3": vol,
            "deepest_point": [round(float(dp[0]), 3), round(float(dp[1]), 3), round(float(dp[2]), 3)],
            "note": "largest enclosed cavity in mid-membrane slab (PCA axis, high-pLDDT core)",
        }, fh, indent=2, ensure_ascii=False)
    print("AF-Q6ZQN7_site.json written")


if __name__ == "__main__":
    main()
