# -*- coding: utf-8 -*-
"""
annotate_pub.py
Pixel-space annotation (two-stage pipeline, stage 2):
  1) fit projection from the 6-sphere calibration frame (least squares, 8 params)
  2) project residue CAs / outer atoms / H-bond endpoints to pixels
  3) residue labels on a ring (azimuth order preserved), leader lines with end
     dots, H-bond distance numbers hugging the dash, caption in Times New Roman
"""
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from PIL import Image, ImageDraw, ImageFont

PUB = os.path.join(PROJECT_ROOT, "work/drug/docking/pub")
FIN = f"{PUB}/final"
os.makedirs(FIN, exist_ok=True)

FONT_ARIAL = "C:/Windows/Fonts/arial.ttf"
FONT_TIMES = "C:/Windows/Fonts/timesbd.ttf"
LINE_GRAY = (70, 70, 70)

TITLES = {
    "PTH1R_ginkgolide": "Ginkgolide-type \u2212 PTH1R (\u22128.8 kcal/mol)",
    "PTH1R_cucurbitacinB": "Cucurbitacin B \u2212 PTH1R (\u22128.4 kcal/mol)",
    "PTH1R_digoxigenin": "Digoxigenin \u2212 PTH1R (\u22127.9 kcal/mol)",
    "SLCO4C1_digitoxin": "Digitoxin \u2212 SLCO4C1 (\u221212.3 kcal/mol)",
    "SLCO4C1_digoxigenin": "Digoxigenin \u2212 SLCO4C1 (\u22129.9 kcal/mol)",
}

MASKS = {
    "red": lambda r, g, b: (r > 150) & (g < 110) & (b < 110),
    "cyan": lambda r, g, b: (g > 150) & (b > 150) & (r < 110),
    "green": lambda r, g, b: (g > 150) & (r < 110) & (b < 110),
    "magenta": lambda r, g, b: (r > 150) & (b > 150) & (g < 110),
    "blue": lambda r, g, b: (b > 150) & (r < 110) & (g < 110),
    "yellow": lambda r, g, b: (r > 150) & (g > 150) & (b < 110),
}
COLORS = ["red", "cyan", "green", "magenta", "blue", "yellow"]
AXES = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)


def blob_center(arr, mask):
    ys, xs = np.nonzero(mask)
    if len(xs) < 30:
        return None
    # centroid (orthoscopic render + flat lighting => solid uniform disc)
    return float(xs.mean()), float(ys.mean())


def fit_projection(calib_png, W, H, SEP):
    im = np.asarray(Image.open(calib_png).convert("RGB"), dtype=np.int16)
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    pts_s, pts_w = [], []
    for color, axis in zip(COLORS, AXES):
        m = MASKS[color](r, g, b)
        c = blob_center(im, m)
        if c is None:
            raise RuntimeError(f"calib sphere {color} not found")
        pts_s.append(c)
        pts_w.append(SEP * axis)
    S = np.array(pts_s)  # (6,2)
    Pw = np.array(pts_w)  # (6,3)
    M = np.zeros((12, 8))
    y = np.zeros(12)
    for i in range(6):
        M[2 * i, 0] = 1
        M[2 * i, 2:5] = Pw[i]
        M[2 * i + 1, 1] = 1
        M[2 * i + 1, 5:8] = Pw[i]
        y[2 * i] = S[i, 0]
        y[2 * i + 1] = S[i, 1]
    x, res, *_ = np.linalg.lstsq(M, y, rcond=None)
    P0 = x[0:2]
    A = np.vstack([x[2:5], x[5:8]])
    pred = M @ x
    rms = float(np.sqrt(np.mean((pred - y) ** 2)))
    sv = np.linalg.svd(A, compute_uv=False)
    return P0, A, rms, sv


def rect_edge(p_from, center, hw, hh, pad=6):
    d = np.asarray(center) - np.asarray(p_from)
    ts = []
    for dd, half in ((d[0], hw + pad), (d[1], hh + pad)):
        ts.append((half / abs(dd)) if abs(dd) > 1e-9 else float("inf"))
    t = min(min(ts), 1.0)
    return np.asarray(center) - d * t


def densify(atoms_px, step=6.0):
    """pixels -> denser point cloud for clearance checks"""
    pts = []
    n = len(atoms_px)
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(atoms_px[i] - atoms_px[j])
            if d < 14:  # stick-connected neighbours
                k = max(int(d // step), 1)
                for t in range(1, k):
                    pts.append(atoms_px[i] + (atoms_px[j] - atoms_px[i]) * t / k)
    pts.extend(atoms_px)
    return np.array(pts)


def annotate(name):
    ann = json.load(open(f"{PUB}/{name}_annot.json", encoding="utf-8"))
    vw = json.load(open(f"{PUB}/{name}_view.json", encoding="utf-8"))
    view, W, H, SEP = vw["view"], vw["W"], vw["H"], vw["SEP"]
    ctr = np.array(view[12:15])

    P0, A, rms, sv = fit_projection(f"{PUB}/{name}_calib.png", W, H, SEP)
    assert abs(sv[0] - sv[1]) / max(sv[0], 1e-9) < 0.03, f"proj singular values differ: {sv}"
    assert rms < 2.5, f"proj rms too large: {rms}"
    print(f"{name}: proj rms={rms:.2f}px  sv=({sv[0]:.1f},{sv[1]:.1f})")

    def proj(p_world):
        p = np.asarray(p_world, dtype=float) - ctr
        s = P0 + A @ p
        return np.array([s[0], s[1]])

    im = Image.open(f"{PUB}/{name}.png").convert("RGB")
    draw = ImageDraw.Draw(im)
    f_lab = ImageFont.truetype(FONT_ARIAL, 34)
    f_num = ImageFont.truetype(FONT_ARIAL, 27)
    f_cap = ImageFont.truetype(FONT_TIMES, 46)

    # collect displayed-atom pixel cloud (annotated residues + ligand) for clearance
    disp = []
    for r in ann["residues"]:
        disp.append(proj(r["ca"]))
    disp.append(proj(ann["lig_centroid"]))
    for h in ann["hbonds"]:
        disp.append(proj(h["d_pos"]))
        disp.append(proj(h["a_pos"]))
    clear_pts = densify(disp)

    # residue label ring
    lig_s = proj(ann["lig_centroid"])
    cas = [np.array(proj(r["ca"])) for r in ann["residues"]]
    content_r = max(np.linalg.norm(c - lig_s) for c in cas) if cas else 100
    R = max(content_r + 90, 0.42 * min(W, H))
    order = np.argsort([np.arctan2(c[1] - lig_s[1], c[0] - lig_s[0]) for c in cas])
    n = len(cas)
    label_px = []
    for k, idx in enumerate(order):
        theta = -np.pi / 2 + 2 * np.pi * k / max(n, 1)  # start top, going clockwise
        pos = lig_s + R * np.array([np.cos(theta), np.sin(theta)])
        pos[0] = min(max(pos[0], 60), W - 60)
        pos[1] = min(max(pos[1], 60), H - 60)
        label_px.append((idx, pos))

    # draw leaders + labels
    label_boxes = []
    for idx, pos in label_px:
        r = ann["residues"][idx]
        ca_s = cas[idx]
        outer_s = proj(r["outer"])
        txt = r["label"]
        bbox = draw.textbbox((0, 0), txt, font=f_lab)
        hw, hh = (bbox[2] - bbox[0]) / 2, (bbox[3] - bbox[1]) / 2
        start = rect_edge(outer_s, pos, hw, hh)
        draw.line([tuple(outer_s), tuple(start)], fill=LINE_GRAY, width=3)
        draw.ellipse([ca_s[0] - 5, ca_s[1] - 5, ca_s[0] + 5, ca_s[1] + 5], fill=LINE_GRAY)
        draw.text((pos[0] - hw, pos[1] - hh), txt, fill="black", font=f_lab)
        label_boxes.append((pos[0] - hw, pos[1] - hh, pos[0] + hw, pos[1] + hh))

    # H-bond distance numbers hugging the dash
    for i, h in enumerate(ann["hbonds"]):
        d_s = proj(h["d_pos"])
        a_s = proj(h["a_pos"])
        seg = a_s - d_s
        L = np.linalg.norm(seg)
        u = seg / L
        nvec = np.array([-u[1], u[0]])
        txt = f"{h['dist']:.1f}"
        bbox = draw.textbbox((0, 0), txt, font=f_num)
        hw, hh = (bbox[2] - bbox[0]) / 2 + 3, (bbox[3] - bbox[1]) / 2 + 3
        placed = None
        for off_t in [0, -0.04, 0.04, -0.08, 0.08, -0.12, 0.12, -0.16, 0.16, -0.24, 0.24, -0.32, 0.32]:
            t = 0.5 + off_t
            base = d_s + seg * t
            for sgn in (1, -1):
                cand = base + nvec * sgn * (hh + 6)
                ok = True
                for p in clear_pts:
                    if abs(p[0] - cand[0]) < hw + 4 and abs(p[1] - cand[1]) < hh + 4:
                        ok = False
                        break
                if ok:
                    for bx0, by0, bx1, by1 in label_boxes:
                        if bx0 - hw < cand[0] < bx1 + hw and by0 - hh < cand[1] < by1 + hh:
                            ok = False
                            break
                if ok:
                    placed = cand
                    break
            if placed is not None:
                break
        if placed is not None:
            draw.text((placed[0] - hw, placed[1] - hh), txt, fill=(60, 60, 60), font=f_num)

    # white-margin canvas + caption
    m = 60
    cap_h = 110
    canvas = Image.new("RGB", (W + 2 * m, H + m + cap_h), "white")
    canvas.paste(im, (m, m))
    d2 = ImageDraw.Draw(canvas)
    cap = TITLES[name]
    cb = d2.textbbox((0, 0), cap, font=f_cap)
    d2.text(((canvas.width - (cb[2] - cb[0])) / 2, H + m + 28), cap, fill="black", font=f_cap)

    out = f"{FIN}/{name}_pub.png"
    canvas.save(out)
    print("saved", out)


if __name__ == "__main__":
    for name in TITLES:
        try:
            annotate(name)
        except Exception as e:
            print(f"ERROR {name}: {e!r}")
