# -*- coding: utf-8 -*-
"""
make_submission_figures.py - build the submission figures (300 dpi, PNG + TIFF, named F1-F6)

Specification (unified):
  - pixel width 2161 px (= 183 mm double-column text width @ 300 dpi).
    2161 px guarantees an effective resolution of at least 300 dpi as long as the final printed width is <= 183 mm.
  - RGB colour mode (alpha removed; transparent regions composited onto white, visually unchanged)
  - one PNG and one TIFF; the TIFF uses LZW compression; both carry 300 dpi metadata
  - proportional resampling; resample only when the width differs from the target (LANCZOS)

Figure sources:
  F1 <- Figure/F1.png                      (originally 2520 px)
  F2 <- Figure/F2.png                      (originally 2161 px, already the target width; F2.tif is a redundant 596 MB file and is not used)
  F3 <- Figure/F3.png                      (originally 3188 px)
  F4 <- panel_sources/F4.pdf re-rendered from vector   (the original PNG is only 1666 px, insufficient; the PDF contains vector text)
  F5 <- Figure/Figure5.png                 (originally 2161 px)
  F6 <- Figure/Figure6.png                 (originally 2956 px)

Source files are never modified; outputs are written to figures/submission/.
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys

from PIL import Image

try:
    import pypdfium2 as pdfium
except ImportError:
    pdfium = None

BASE = PROJECT_ROOT
FIGDIR = os.path.join(BASE, "figures")
OUTDIR = os.path.join(FIGDIR, "submission")

TARGET_W = 2161          # px, 183 mm @ 300 dpi
DPI = (300, 300)

RASTER_SRC = {
    "F1": os.path.join(FIGDIR, "Figure", "F1.png"),
    "F2": os.path.join(FIGDIR, "Figure", "F2.png"),
    "F3": os.path.join(FIGDIR, "Figure", "F3.png"),
    "F5": os.path.join(FIGDIR, "Figure", "Figure5.png"),
    "F6": os.path.join(FIGDIR, "Figure", "Figure6.png"),
}
PDF_SRC = {"F4": os.path.join(FIGDIR, "panel_sources", "F4.pdf")}

ORDER = ["F1", "F2", "F3", "F4", "F5", "F6"]


def flatten_rgb(im):
    """Remove alpha: composite transparent regions onto a white background."""
    if im.mode in ("RGBA", "LA") or (im.mode == "P" and "transparency" in im.info):
        im = im.convert("RGBA")
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[-1])
        return bg
    return im.convert("RGB")


def from_pdf(path, target_w):
    if pdfium is None:
        raise RuntimeError("pypdfium2 is required to render PDF pages")
    doc = pdfium.PdfDocument(path)
    page = doc[0]
    w_pt, h_pt = page.get_size()
    scale = target_w / w_pt          # pt -> px (scale = 1 when 1 pt = 1 px)
    bmp = page.render(scale=scale)
    im = bmp.to_pil().convert("RGB")
    doc.close()
    return im


def fit_width(im, target_w):
    if im.width == target_w:
        return im, 1.0
    h = max(1, round(im.height * target_w / im.width))
    return im.resize((target_w, h), Image.LANCZOS), target_w / im.width


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = []

    for key in ORDER:
        if key in PDF_SRC:
            src = PDF_SRC[key]
            im = from_pdf(src, TARGET_W)
        else:
            src = RASTER_SRC[key]
            im = flatten_rgb(Image.open(src))
            im, _ = fit_width(im, TARGET_W)

        if im.mode != "RGB":
            im = flatten_rgb(im)

        png = os.path.join(OUTDIR, f"{key}.png")
        tif = os.path.join(OUTDIR, f"{key}.tif")
        im.save(png, dpi=DPI)
        im.save(tif, format="TIFF", compression="tiff_lzw", dpi=DPI)

        rows.append((key, os.path.relpath(src, BASE), im.size,
                     round(os.path.getsize(png) / 1024, 1),
                     round(os.path.getsize(tif) / 1024, 1)))
        print(f"{key}: {os.path.basename(src):14s} -> {im.size[0]}x{im.size[1]} "
              f"PNG {rows[-1][3]:.0f}KB / TIFF {rows[-1][4]:.0f}KB")

    print(f"\noutput directory: {OUTDIR}")
    print(f"{len(rows)} figures x 2 formats = {len(rows)*2} files; "
          f"uniform {TARGET_W} px wide @300 dpi = 183 mm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
