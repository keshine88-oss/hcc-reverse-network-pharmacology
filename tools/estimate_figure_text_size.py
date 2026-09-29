# -*- coding: utf-8 -*-
"""probe_text_size.py - estimate the printed point size of body text in a finished figure (approximate)

Approach and assumptions:
  1. Text-like components from connected-component analysis: height 8-40 px, area >= 20, fill ratio < 0.80
     (high fill ratio with a nearly square shape indicates scatter dots; they must be excluded, otherwise the 5th percentile is dominated by point noise).
  2. Glyph height mainly reflects x-height -> divide by 0.55 to recover em (Arial x-height is about 0.52-0.55 em).
  3. Two conversion routes:
     - declared-dpi route: em_px / declared dpi x 72
     - print route: design width = px / declared dpi, i.e. the physical width the author designed at;
       and if the final printed width is W mm, the point size = declared point size x (W / design width)
  4. The result is an order-of-magnitude estimate for flagging text that is clearly too small; it is not a substitute for visual inspection.
"""
import os
import sys
import numpy as np
from PIL import Image
from scipy import ndimage

W_PRINT = [183.0, 90.0]     # double column / single column
XHEIGHT = 0.55


def probe(path):
    im = Image.open(path)
    g = np.asarray(im.convert("L"))
    lab, n = ndimage.label(g < 128)
    objs = ndimage.find_objects(lab)
    hs = []
    for i, sl in enumerate(objs):
        if sl is None:
            continue
        ys, xs = sl
        h, w = ys.stop - ys.start, xs.stop - xs.start
        area = int((lab[sl] == i + 1).sum())
        if not (8 <= h <= 40 and area >= 20 and w <= 100):
            continue
        fill = area / float(h * w)
        if fill > 0.80 and abs(h - w) <= 2:      # dots / scatter points
            continue
        hs.append((h, w, fill))
    return im, np.array([x[0] for x in hs], dtype=float)


def main(src):
    files = [os.path.join(src, f) for f in sorted(os.listdir(src))
             if f.lower().endswith((".png", ".tif", ".tiff"))]
    print(f"{'file':16s}{'glyph height(px)':>18s}{'declared dpi':>13s}{'design width(mm)':>18s}"
          f"{'declared pt':>11s}{'print 183mm':>10s}{'print 90mm':>9s}")
    print("-" * 84)
    for p in files:
        im, hs = probe(p)
        if hs.size < 10:
            print(f"{os.path.basename(p):16s}  too few text components, skipping")
            continue
        h = float(np.percentile(hs, 15))          # take the smaller size band
        px = im.size[0]
        d = im.info.get("dpi")
        try:
            d = float(d[0] if isinstance(d, (tuple, list)) else d)
        except (TypeError, ValueError):
            d = None
        if not d:
            print(f"{os.path.basename(p):16s}  no DPI metadata, cannot estimate")
            continue
        em_px = h / XHEIGHT
        pt_decl = em_px / d * 72
        design_mm = px / d * 25.4
        outs = []
        for W in W_PRINT:
            outs.append(pt_decl * W / design_mm)
        flag = "  WARNING: below 5 pt when printed" if min(outs) < 5 else ""
        print(f"{os.path.basename(p):16s}{h:12.1f}{d:9.0f}{design_mm:11.1f}"
              f"{pt_decl:11.1f}{outs[0]:10.1f}{outs[1]:9.1f}{flag}")
    print("\nNote: design width = pixels / declared dpi, i.e. the physical width the author designed at; "
          "printed point size = declared point size x (final width / design width).")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
