# -*- coding: utf-8 -*-
"""
compose_combined.py
Assemble the reference-figure layout: overview (left, dashed box around the
site) + leader dashed lines + close-up (right, dashed border, residue labels
NAME-resi, H-bond distances offset from the dash). All annotations in English.
"""
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, os.path.join(PROJECT_ROOT, "docking"))
from annotate_pub import fit_projection  # noqa: E402

PUB = os.path.join(PROJECT_ROOT, "work/drug/docking/pub")
FIN = f"{PUB}/combined"
os.makedirs(FIN, exist_ok=True)

FONT_ARIAL = "C:/Windows/Fonts/arial.ttf"
FONT_TIMES = "C:/Windows/Fonts/timesbd.ttf"
LINE_GRAY = (70, 70, 70)
BLACK = (20, 20, 20)

TITLES = {
    "PTH1R_ginkgolide": "Ginkgolide C \u2212 PTH1R (\u22128.8 kcal/mol)",
    "PTH1R_cucurbitacinB": "Cucurbitacin B \u2212 PTH1R (\u22128.4 kcal/mol)",
    "PTH1R_digoxigenin": "Digoxigenin \u2212 PTH1R (\u22127.9 kcal/mol)",
    "PTH1R_kamebanin": "Kamebanin \u2212 PTH1R (\u22127.8 kcal/mol)",
    "PTH1R_roemerine": "Roemerine \u2212 PTH1R (\u22127.8 kcal/mol)",
    "PTH1R_chrysin": "Chrysin \u2212 PTH1R (\u22127.4 kcal/mol)",
    "PTH1R_tryptanthrin": "Tryptanthrin \u2212 PTH1R (\u22127.3 kcal/mol)",
    "PTH1R_anisaldehyde": "Anisaldehyde \u2212 PTH1R (\u22124.8 kcal/mol)",
    "SLCO4C1_digitoxin": "Digitoxin \u2212 SLCO4C1 (\u221212.3 kcal/mol)",
    "SLCO4C1_digoxigenin": "Digoxigenin \u2212 SLCO4C1 (\u22129.9 kcal/mol)",
}


def dashed_line(draw, p0, p1, color, width=3, dash=14, gap=9):
    p0 = np.array(p0, float)
    p1 = np.array(p1, float)
    d = p1 - p0
    L = float(np.linalg.norm(d))
    if L < 1e-6:
        return
    u = d / L
    t = 0.0
    while t < L:
        t2 = min(t + dash, L)
        a = p0 + u * t
        b = p0 + u * t2
        draw.line([tuple(a), tuple(b)], fill=color, width=width)
        t = t2 + gap


def proj_fn(P0, A, ctr):
    def f(p_world):
        p = np.asarray(p_world, float) - ctr
        s = P0 + A @ p
        return np.array([s[0], s[1]])
    return f


def rect_edge(p_from, center, hw, hh, pad=6):
    d = np.asarray(center) - np.asarray(p_from)
    ts = []
    for dd, half in ((d[0], hw + pad), (d[1], hh + pad)):
        ts.append((half / abs(dd)) if abs(dd) > 1e-9 else float("inf"))
    t = min(min(ts), 1.0)
    return np.asarray(center) - d * t


def densify(atoms_px, step=6.0):
    pts = []
    n = len(atoms_px)
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(atoms_px[i] - atoms_px[j])
            if d < 14:
                k = max(int(d // step), 1)
                for t in range(1, k):
                    pts.append(atoms_px[i] + (atoms_px[j] - atoms_px[i]) * t / k)
    pts.extend(atoms_px)
    return np.array(pts)


def annotate_closeup(name):
    """annotate the closeup in place (labels + distances); return the PIL image"""
    ann = json.load(open(f"{PUB}/{name}_annot.json", encoding="utf-8"))
    vw = json.load(open(f"{PUB}/{name}_views.json", encoding="utf-8"))
    W, H = vw["cl_size"]
    ctr = np.array(vw["closeup"][12:15])
    P0, A, rms, sv = fit_projection(f"{PUB}/{name}_closeup_calib.png", W, H, vw["SEP"])
    assert abs(sv[0] - sv[1]) / max(sv[0], 1e-9) < 0.03 and rms < 2.5, (rms, sv)

    def proj(p):
        p = np.asarray(p, float) - ctr
        return np.array(P0 + A @ p)

    im = Image.open(f"{PUB}/{name}_closeup.png").convert("RGB")
    draw = ImageDraw.Draw(im)
    f_lab = ImageFont.truetype(FONT_ARIAL, 58)
    f_num = ImageFont.truetype(FONT_ARIAL, 48)

    # clearance cloud: ligand heavy atoms + H-bond endpoints (NOT residues, so
    # that labels may sit close to their own residue)
    lig_s = proj(ann["lig_centroid"])
    clear_pts = densify([proj(a) for a in ann["lig_atoms"]]
                        + [proj(h["d_pos"]) for h in ann["hbonds"]]
                        + [proj(h["a_pos"]) for h in ann["hbonds"]])

    res = ann["residues"]
    n = len(res)
    cas = [proj(r["ca"]) for r in res]
    outers = [proj(r["outer"]) for r in res]
    # radial direction: centroid -> outermost atom of the residue
    import math
    base_ang = []
    label_sz = []
    for i, r in enumerate(res):
        vec = outers[i] - lig_s
        d = float(np.linalg.norm(vec))
        base_ang.append(math.atan2(vec[1], vec[0]) if d > 1e-6 else -np.pi / 2)
        bbox = draw.textbbox((0, 0), r["label"], font=f_lab)
        label_sz.append(((bbox[2] - bbox[0]) / 2, (bbox[3] - bbox[1]) / 2))

    # label anchored just outside its outermost atom (short leader, tight layout)
    offsets = [hw + 30 for hw, hh in label_sz]
    ang_off = [0.0] * n

    def pos_of(i):
        ang = base_ang[i] + ang_off[i]
        u = np.array([math.cos(ang), math.sin(ang)])
        return outers[i] + u * offsets[i]

    # relax: separate overlapping labels by rotating them apart (tangential),
    # keep labels off the ligand by pushing radially
    for _ in range(60):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                pi, pj = pos_of(i), pos_of(j)
                hwi, hhi = label_sz[i]
                hwj, hhj = label_sz[j]
                if abs(pi[0] - pj[0]) < hwi + hwj + 18 and abs(pi[1] - pj[1]) < hhi + hhj + 14:
                    ang_off[i] += 0.11
                    ang_off[j] -= 0.11
                    moved = True
        for i in range(n):
            pi = pos_of(i)
            hwi, hhi = label_sz[i]
            for p in clear_pts:
                if abs(p[0] - pi[0]) < hwi + 4 and abs(p[1] - pi[1]) < hhi + 4:
                    offsets[i] += 10
                    moved = True
                    break
        if not moved:
            break

    positions = []
    label_boxes = []
    for i in range(n):
        p = pos_of(i)
        hw, hh = label_sz[i]
        p = np.array([min(max(p[0], hw + 24), W - hw - 24),
                      min(max(p[1], hh + 24), H - hh - 24)])
        # if the label hugs the border, snap next to its CA instead
        hugging = (p[0] <= hw + 26 or p[0] >= W - hw - 26 or
                   p[1] <= hh + 26 or p[1] >= H - hh - 26)
        if hugging:
            for cand in (cas[i] + np.array([0, -80]), cas[i] + np.array([0, 80]),
                         cas[i] + np.array([-90, 0]), cas[i] + np.array([90, 0])):
                if hw + 8 <= cand[0] <= W - hw - 8 and hh + 8 <= cand[1] <= H - hh - 8:
                    ok = all(not (abs(q[0] - cand[0]) < hw + 4 and abs(q[1] - cand[1]) < hh + 4)
                             for q in clear_pts)
                    if ok:
                        p = cand
                        break
        positions.append(p)

    for i in range(n):
        r = res[i]
        outer_s = outers[i]
        ca_s = cas[i]
        pos = positions[i]
        hw, hh = label_sz[i]
        start = rect_edge(outer_s, pos, hw, hh)
        draw.line([tuple(outer_s), tuple(start)], fill=LINE_GRAY, width=3)
        draw.ellipse([ca_s[0] - 6, ca_s[1] - 6, ca_s[0] + 6, ca_s[1] + 6], fill=LINE_GRAY)
        draw.text((pos[0] - hw, pos[1] - hh), r["label"], fill=BLACK, font=f_lab)
        label_boxes.append((pos[0] - hw, pos[1] - hh, pos[0] + hw, pos[1] + hh))

    # H-bond dashes: drawn here from the SAME projection as the distance labels,
    # so the dashed line and its number are always in register.
    HB_YELLOW = (219, 204, 56)
    for h in ann["hbonds"]:
        dashed_line(draw, tuple(proj(h["d_pos"])), tuple(proj(h["a_pos"])),
                    HB_YELLOW, width=5, dash=16, gap=10)

    # H-bond distance numbers: hug the dash tightly, always tied to its own bond.
    # Priority: sit ON the dash (perp=0, white halo) at the midpoint, then small
    # slides along the dash, then a tight perpendicular offset. Any displacement
    # from the midpoint gets a thin connector line back to the bond midpoint so
    # the number-to-bond correspondence is unambiguous.
    for h in ann["hbonds"]:
        d_s = proj(h["d_pos"])
        a_s = proj(h["a_pos"])
        seg = a_s - d_s
        L = float(np.linalg.norm(seg))
        if L < 1e-6:
            continue
        u = seg / L
        nvec = np.array([-u[1], u[0]])
        mid = d_s + seg * 0.5
        txt = f"{h['dist']:.1f}"
        bbox = draw.textbbox((0, 0), txt, font=f_num)
        hw, hh = (bbox[2] - bbox[0]) / 2 + 2, (bbox[3] - bbox[1]) / 2 + 2

        def free(cand):
            for p in clear_pts:
                if abs(p[0] - cand[0]) < hw + 2 and abs(p[1] - cand[1]) < hh + 2:
                    return False
            for bx0, by0, bx1, by1 in label_boxes:
                if bx0 - hw < cand[0] < bx1 + hw and by0 - hh < cand[1] < by1 + hh:
                    return False
            return True

        # candidate spots, closest-to-bond first
        cands = []                                   # (position, on_dash)
        for off_t in [0.0, 0.07, -0.07, 0.14, -0.14]:          # ON the dash
            cands.append((d_s + seg * (0.5 + off_t), True))
        for off_t in [0.0, 0.07, -0.07]:                        # tight perp offset
            for sgn in (1, -1):
                cands.append((d_s + seg * (0.5 + off_t) + nvec * sgn * (hh + 1.5), False))

        placed, on_dash = None, False
        for cand, od in cands:
            if free(cand):
                placed, on_dash = cand, od
                break
        if placed is None:                            # last resort: midpoint on-dash
            placed, on_dash = mid, True

        if not on_dash:
            # thin connector from the number to the bond midpoint
            draw.line([tuple(placed), tuple(mid)], fill=(150, 150, 150), width=2)
        else:
            # white halo so the number reads clearly where it sits on the dash
            draw.rectangle([placed[0] - hw, placed[1] - hh,
                            placed[0] + hw, placed[1] + hh], fill="white")
        draw.text((placed[0] - hw + 2, placed[1] - hh + 2), txt, fill=BLACK, font=f_num)
    return im


def main():
    for name in TITLES:
        try:
            ann = json.load(open(f"{PUB}/{name}_annot.json", encoding="utf-8"))
            vw = json.load(open(f"{PUB}/{name}_views.json", encoding="utf-8"))
            ov_W, ov_H = vw["ov_size"]
            ctr = np.array(vw["overview"][12:15])
            P0, A, rms, sv = fit_projection(f"{PUB}/{name}_overview_calib.png", ov_W, ov_H, vw["SEP"])
            assert abs(sv[0] - sv[1]) / max(sv[0], 1e-9) < 0.05 and rms < 3.0, (rms, sv)

            def proj(p):
                p = np.asarray(p, float) - ctr
                return np.array(P0 + A @ p)

            ov = Image.open(f"{PUB}/{name}_overview.png").convert("RGB")
            d_ov = ImageDraw.Draw(ov)
            lig_s = proj(ann["lig_centroid"])
            # site box: projection of a 13 A sphere around the centroid
            probe = proj((np.array(ann["lig_centroid"]) + np.array([13.0, 13.0, 13.0])).tolist())
            r_px = float(np.max(np.abs(probe - lig_s))) * 0.85
            r_px = min(max(r_px, 90), ov_W * 0.42)
            box = [lig_s[0] - r_px, lig_s[1] - r_px, lig_s[0] + r_px, lig_s[1] + r_px]
            dashed_line(d_ov, (box[0], box[1]), (box[2], box[1]), BLACK, 4)
            dashed_line(d_ov, (box[2], box[1]), (box[2], box[3]), BLACK, 4)
            dashed_line(d_ov, (box[2], box[3]), (box[0], box[3]), BLACK, 4)
            dashed_line(d_ov, (box[0], box[3]), (box[0], box[1]), BLACK, 4)

            cl = annotate_closeup(name)

            # layout: equal height, gap, margins
            target_h = 1100
            ov2 = ov.resize((int(ov.width * target_h / ov.height), target_h), Image.LANCZOS)
            cl2 = cl.resize((int(cl.width * target_h / cl.height), target_h), Image.LANCZOS)
            m, gap, cap_h = 46, 90, 120
            canvas = Image.new("RGB", (m + ov2.width + gap + cl2.width + m, target_h + m + cap_h + m), "white")
            cx0 = m
            cy0 = m
            cx1 = m + ov2.width + gap
            canvas.paste(ov2, (cx0, cy0))
            canvas.paste(cl2, (cx1, cy0))
            d = ImageDraw.Draw(canvas)

            # dashed border around closeup
            bx0, by0 = cx1, cy0
            bx1, by1 = cx1 + cl2.width, cy0 + cl2.height
            dashed_line(d, (bx0, by0), (bx1, by0), BLACK, 5)
            dashed_line(d, (bx1, by0), (bx1, by1), BLACK, 5)
            dashed_line(d, (bx1, by1), (bx0, by1), BLACK, 5)
            dashed_line(d, (bx0, by1), (bx0, by0), BLACK, 5)

            # connector lines from overview box corners to closeup border
            sx = ov2.width / ov.width
            sy = ov2.height / ov.height
            tl = (cx0 + box[2] * sx, cy0 + box[1] * sy)
            bl = (cx0 + box[2] * sx, cy0 + box[3] * sy)
            dashed_line(d, tl, (bx0, by0 + 40), BLACK, 4)
            dashed_line(d, bl, (bx0, by1 - 40), BLACK, 4)

            cap = TITLES[name]
            f_cap = ImageFont.truetype(FONT_TIMES, 52)
            cb = d.textbbox((0, 0), cap, font=f_cap)
            d.text(((canvas.width - (cb[2] - cb[0])) / 2, target_h + m + 34), cap, fill="black", font=f_cap)

            out = f"{FIN}/{name}_combined.png"
            canvas.save(out)
            print("saved", out, canvas.size)
        except Exception as e:
            print(f"ERROR {name}: {e!r}")


if __name__ == "__main__":
    main()
