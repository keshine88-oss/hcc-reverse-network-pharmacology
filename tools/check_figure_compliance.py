# -*- coding: utf-8 -*-
"""
check_figures.py - submission figure compliance check (read-only, does not modify any file)

Usage:
    python check_figures.py <figure-directory> [--min-dpi 300] [--col-mm 183]

Judgement basis (effective print resolution = pixels / printed inches, not the DPI declared in the file):
  - halftone images (photographs / heatmaps / stains)   >= 300 dpi
  - combination figures (images + lines/text)      >= 500 dpi
  - pure line art / schematics                     >= 900 dpi (some journals require 1200)
The script cannot classify the three categories automatically, so it applies --min-dpi (300) as the ERROR threshold,
Below 500 dpi a WARN is raised saying "if this figure contains line art, the combination-art threshold applies", to avoid false negatives.

Other checks: file format, colour mode (alpha / palette), bit depth (>8 bit per channel),
physical size (exceeding the page text width), missing DPI metadata, naming and figure-number continuity, duplicate content, TIFF/EPS companions.
"""
import os
import re
import sys
import glob
import hashlib
import argparse

try:
    from PIL import Image
except ImportError:
    print("Pillow is required: pip install Pillow")
    sys.exit(2)

VECTOR_EXT = {".pdf", ".eps", ".ai", ".svg"}
RASTER_EXT = {".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp", ".gif"}
STRONG_RASTER = {".tif", ".tiff"}
WEAK_RASTER = {".png", ".jpg", ".jpeg", ".bmp", ".gif"}

HALFTONE_DPI = 300
COMBO_DPI = 500
LINEART_DPI = 900
LOWBIT_MODES = {"I", "I;16", "I;16B", "I;16L", "I;16N", "F"}   # >8 bit per channel


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def analyse(path):
    ext = os.path.splitext(path)[1].lower()
    info = {"path": path, "name": os.path.basename(path), "ext": ext,
            "size_kb": round(os.path.getsize(path) / 1024, 1),
            "kind": "vector" if ext in VECTOR_EXT else "raster"}
    if info["kind"] == "vector":
        return info
    try:
        im = Image.open(path)
    except Exception as e:
        info["error"] = f"cannot open: {e}"
        return info

    info["format"] = im.format
    info["mode"] = im.mode
    info["px"] = im.size
    info["n_frames"] = getattr(im, "n_frames", 1)
    info["compression"] = im.info.get("compression")

    dpi = im.info.get("dpi")
    d = None
    if dpi is not None:
        try:
            d = float(dpi[0] if isinstance(dpi, (tuple, list)) else dpi)
        except (TypeError, ValueError):
            d = None
    info["declared_dpi"] = d
    if d and d > 1:
        info["phys_mm"] = round(im.size[0] / d * 25.4, 1)
        info["phys_mm_h"] = round(im.size[1] / d * 25.4, 1)

    info["has_alpha"] = im.mode in ("RGBA", "LA") or (
        im.mode == "P" and "transparency" in im.info)
    info["deep_bit"] = im.mode in LOWBIT_MODES
    return info


def figure_no(name):
    """Extract the figure number from a filename: Figure1 / Fig2 / F3 are all accepted."""
    m = re.match(r"^(?:figure|fig|f)[\s_\-]*(\d+)", name, re.I)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--min-dpi", type=int, default=HALFTONE_DPI)
    ap.add_argument("--col-mm", type=float, default=183.0, help="page text width in mm (double column)")
    args = ap.parse_args()

    folder = args.folder
    if not os.path.isdir(folder):
        print(f"[error] directory does not exist: {folder}")
        return 2

    files = sorted(p for p in glob.glob(os.path.join(folder, "*"))
                   if os.path.splitext(p)[1].lower() in (VECTOR_EXT | RASTER_EXT))
    if not files:
        print(f"[error] no image files found in the directory")
        return 2

    rows, errors, warns, infos = [], [], [], []
    hashes, figmap = {}, {}

    for p in files:
        a = analyse(p)
        rows.append(a)
        n = a["name"]
        if a.get("error"):
            errors.append(f"{n}：{a['error']}")
            continue

        if a["kind"] == "vector":
            infos.append(f"{n}: vector format ({a['ext'][1:].upper()}), no dpi limit")
        else:
            if a["ext"] in WEAK_RASTER:
                warns.append(f"{n}: {a['ext'][1:].upper()} is often not accepted for main figures; use TIFF/EPS")
            eff = a["px"][0] / (args.col_mm / 25.4)

            if eff < args.min_dpi - 0.5:
                errors.append(
                    f"{n}: when printed at {args.col_mm:.0f} mm width the effective resolution is only {eff:.0f} dpi "
                    f"< {args.min_dpi} dpi (pixel width {a['px'][0]} is insufficient)")
            elif eff < COMBO_DPI:
                warns.append(
                    f"{n}: effective resolution {eff:.0f} dpi; if this figure is line art or a combination figure, "
                    f"most journals require >= {COMBO_DPI} (line art {LINEART_DPI}-1200)")
            if not a["declared_dpi"]:
                warns.append(f"{n}: no DPI metadata found (reads back as None); "
                             f"the submission system may read it as 72 dpi; please embed the DPI value")
            if a["has_alpha"]:
                errors.append(f"{n}: colour mode {a['mode']} contains a transparency channel; convert to RGB")
            elif a["mode"] == "P":
                warns.append(f"{n}: palette mode (P); convert to RGB")
            elif a["mode"] not in ("RGB", "CMYK", "L", "1"):
                warns.append(f"{n}: colour mode {a['mode']}, please confirm")
            if a["deep_bit"]:
                warns.append(f"{n}: mode {a['mode']} is >8 bit per channel; most journals accept 8-bit only")
            if a.get("phys_mm") and a["phys_mm"] > args.col_mm * 1.02:
                warns.append(f"{n}: declared physical width {a['phys_mm']} mm exceeds the page text width "
                             f"({args.col_mm:.0f} mm); the figure must be re-laid out to a column width")
            if a["n_frames"] and a["n_frames"] > 1:
                warns.append(f"{n}: multi-page file ({a['n_frames']} pages); only single-page files are accepted")

        hashes.setdefault(md5(p), []).append(n)
        fn = figure_no(a["name"])
        if fn:
            figmap.setdefault(fn, set()).add(a["ext"])

    # Duplicate content
    for h, names in hashes.items():
        if len(names) > 1:
            warns.append(f"files with identical content: {', '.join(names)}")

    # Figure-number continuity
    if figmap:
        got = sorted(figmap)
        missing = [i for i in range(1, max(got) + 1) if i not in figmap]
        if missing:
            warns.append(f"figure numbers are not contiguous; missing: {', '.join(map(str, missing))}")
        infos.append(f"figure numbers detected: {got} ({len(got)} figures)")

    # TIFF/EPS companion files
    stems = {}
    for p in files:
        s = os.path.splitext(os.path.basename(p))[0].lower()
        stems.setdefault(s, set()).add(os.path.splitext(p)[1].lower())
    for s, exts in sorted(stems.items()):
        if not (exts & (STRONG_RASTER | VECTOR_EXT)):
            infos.append(f"{s}: only {'/'.join(sorted(exts))}; no TIFF/EPS version")

    # ---- Output ----
    print(f"\nchecking directory: {folder}")
    print(f"figures: {len(files)}    page text width: {args.col_mm:.0f} mm    "
          f"minimum effective dpi: {args.min_dpi}\n")
    hdr = f"{'file':40s}{'format':7s}{'mode':6s}{'pixels':13s}{'declared dpi':13s}{'effective dpi':15s}{'width mm':9s}"
    print(hdr)
    print("-" * len(hdr))
    for a in rows:
        if a.get("error") or a["kind"] == "vector":
            tail = "vector" if a["kind"] == "vector" else "-"
            print(f"{a['name'][:39]:40s}{a['ext'][1:].upper():7s}{tail:6s}"
                  f"{'—':13s}{'—':9s}{'—':9s}{'—':8s}")
            continue
        eff = a["px"][0] / (args.col_mm / 25.4)
        d = a["declared_dpi"]
        print(f"{a['name'][:39]:40s}{str(a['format']):7s}{a['mode']:6s}"
              f"{a['px'][0]}x{a['px'][1]:<7d}{'—' if not d else f'{d:.0f}':9s}"
              f"{eff:.0f}{'':5s}{a.get('phys_mm', '—')}")

    for title, items in (("ERROR (must fix)", errors), ("WARN (recommended)", warns),
                         ("INFO (reference)", infos)):
        print(f"\n[{title}] {len(items)} items")
        for i in items:
            print("  - " + i)

    verdict = "must be fixed before submission" if errors else ("usable, optimisation recommended" if warns else "pass")
    print(f"\n=== Verdict: {verdict} ===")
    print("Thresholds: 300 / 500 / 900-1200 dpi correspond to halftone / combination / line art; "
          "the target journal's Author Guidelines take precedence.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
