# -*- coding: utf-8 -*-
"""
parse_npass_compound_pages.py
Parse NPASS compound detail pages (npc_pages/NPC*.html) and extract:
  1) Activity rows against the two core HCC targets:
       NPT535  = PTH1R  (parathyroid hormone 1 receptor)
       NPT3812 = SLCO4C1 (solute carrier organic anion transporter 4C1)
  2) Species source records (TCM herb mapping)
  3) Common name + SMILES per compound

Outputs (written to G:/hcc_drug/):
  npass_tcm_activity.csv  - one row per (compound x target x activity record)
  npass_tcm_species.csv   - one row per (compound x organism source)
  npass_tcm_compounds.csv - one row per compound (name, SMILES, best potency)

Pure stdlib; no network access. Run after the page-fetch loop finishes.
"""
import csv
import html
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import re
import sys

PAGE_DIR = os.path.join(PROJECT_ROOT, "work/drug/herb/npc_pages")
OUT_DIR = os.path.join(PROJECT_ROOT, "work/drug")

TARGETS = {
    "NPT535": "PTH1R",
    "NPT3812": "SLCO4C1",
}

# section markers in the compound page (positions decide which table a row belongs to)
SECTION_MARKERS = [
    ("molecular", "<div id = \"MolLevel\""),
    ("in_vitro", "<div id=\"InVitro\""),
    ("in_vivo", "<div id=\"InVivo\""),
    ("species", "CompoundOrganismTable"),
]

TR_RE = re.compile(r"<tr align=center>(.*?)</tr>", re.DOTALL)
TD_RE = re.compile(r"<td[^>]*>(.*?)</td>", re.DOTALL)
TAG_RE = re.compile(r"<[^>]+>")
TARGET_ID_RE = re.compile(r"target\.php\?target_id=(NPT\d+)")
ORG_ID_RE = re.compile(r"organism\.php\?org_id=(NPO\d+)")


def clean(s):
    """strip tags and collapse whitespace"""
    s = TAG_RE.sub(" ", s)
    s = html.unescape(s)
    return " ".join(s.split())


def section_of(pos, markers):
    """which section does byte position pos fall into?"""
    current = "header"
    for name, mpos in markers:
        if mpos != -1 and pos > mpos:
            current = name
    return current


def parse_page(path):
    npc_id = os.path.basename(path).replace(".html", "")
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()

    # skip error pages / SPA shells (real pages are >100KB and contain the org table or activity tables)
    if "CompoundActivityTable" not in text and "CompoundOrganismTable" not in text:
        return npc_id, None

    # ---- common name ----
    name = ""
    m = re.search(r"<b>Common Name</b>(?:(?!</tr>).)*?<td align=left>([^<]*)</td>", text, re.DOTALL)
    if m:
        name = clean(m.group(1))

    # ---- SMILES ----
    smiles = ""
    m = re.search(r"<b>SMILES</b></td>\s*<td[^>]*>([^<]+)</td>", text)
    if m:
        smiles = clean(m.group(1))

    # ---- section marker positions ----
    markers = []
    for sname, marker in SECTION_MARKERS:
        markers.append((sname, text.find(marker)))
    markers.sort(key=lambda x: x[1] if x[1] != -1 else 10**9)

    activity_rows = []
    species_rows = []

    for tr in TR_RE.finditer(text):
        row_html = tr.group(1)
        cells = TD_RE.findall(row_html)
        if not cells:
            continue
        sec = section_of(tr.start(), markers)

        if sec in ("molecular", "in_vitro", "in_vivo"):
            tmatch = TARGET_ID_RE.search(cells[0])
            if not tmatch:
                continue
            tid = tmatch.group(1)
            if tid not in TARGETS:
                continue
            cells_c = [clean(c) for c in cells]
            while len(cells_c) < 9:
                cells_c.append("")
            activity_rows.append({
                "npc_id": npc_id,
                "compound_name": name,
                "section": sec,
                "gene": TARGETS[tid],
                "target_id": tid,
                "target_type": cells_c[1],
                "target_name": cells_c[2],
                "target_organism": cells_c[3],
                "activity_type": cells_c[4],
                "activity_relation": cells_c[5],
                "activity_value": cells_c[6],
                "activity_unit": cells_c[7],
                "reference": cells_c[8],
            })
        elif sec == "species":
            omatch = ORG_ID_RE.search(cells[0])
            if not omatch:
                continue
            cells_c = [clean(c) for c in cells]
            while len(cells_c) < 9:
                cells_c.append("")
            species_rows.append({
                "npc_id": npc_id,
                "compound_name": name,
                "organism_id": omatch.group(1),
                "organism_name": cells_c[1],
                "taxonomy_level": cells_c[2],
                "family": cells_c[3],
                "superkingdom": cells_c[4],
                "isolation_part": cells_c[5],
                "reference": cells_c[8],
            })

    return npc_id, {
        "npc_id": npc_id,
        "compound_name": name,
        "smiles": smiles,
        "activity": activity_rows,
        "species": species_rows,
    }


def potency_key(row):
    """sort key: numeric value in nM, smaller = more potent; non-numeric last"""
    try:
        v = float(row["activity_value"])
        unit = row["activity_unit"].lower()
        if unit in ("um", "µm"):
            v *= 1000.0
        elif unit == "mm":
            v *= 1e6
        return (0, v)
    except ValueError:
        return (1, 1e18)


def main():
    files = sorted(f for f in os.listdir(PAGE_DIR) if f.startswith("NPC") and f.endswith(".html"))
    print(f"{len(files)} page files found")

    all_act, all_species, compounds = [], [], []
    bad = []
    for f in files:
        npc_id, rec = parse_page(os.path.join(PAGE_DIR, f))
        if rec is None:
            bad.append(npc_id)
            continue
        all_act.extend(rec["activity"])
        all_species.extend(rec["species"])
        compounds.append({
            "npc_id": npc_id,
            "compound_name": rec["compound_name"],
            "smiles": rec["smiles"],
            "n_activity_records": len(rec["activity"]),
            "n_species": len(rec["species"]),
        })

    # ---- best potency per compound (for ranking) ----
    best = {}
    for r in all_act:
        k = (r["npc_id"], r["gene"])
        if k not in best or potency_key(r) < potency_key(best[k]):
            best[k] = r
    for c in compounds:
        for gene in ("PTH1R", "SLCO4C1"):
            b = best.get((c["npc_id"], gene))
            if b:
                c[f"best_{gene}_type"] = b["activity_type"]
                c[f"best_{gene}_value"] = b["activity_value"]
                c[f"best_{gene}_unit"] = b["activity_unit"]
            else:
                c[f"best_{gene}_type"] = c[f"best_{gene}_value"] = c[f"best_{gene}_unit"] = ""

    os.makedirs(OUT_DIR, exist_ok=True)

    with open(os.path.join(OUT_DIR, "npass_tcm_activity.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_act[0].keys()) if all_act else
                           ["npc_id", "compound_name", "section", "gene", "target_id", "target_type",
                            "target_name", "target_organism", "activity_type", "activity_relation",
                            "activity_value", "activity_unit", "reference"])
        w.writeheader()
        w.writerows(all_act)

    with open(os.path.join(OUT_DIR, "npass_tcm_species.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["npc_id", "compound_name", "organism_id", "organism_name",
                                           "taxonomy_level", "family", "superkingdom",
                                           "isolation_part", "reference"])
        w.writeheader()
        w.writerows(all_species)

    with open(os.path.join(OUT_DIR, "npass_tcm_compounds.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["npc_id", "compound_name", "smiles",
                                           "n_activity_records", "n_species",
                                           "best_PTH1R_type", "best_PTH1R_value", "best_PTH1R_unit",
                                           "best_SLCO4C1_type", "best_SLCO4C1_value", "best_SLCO4C1_unit"])
        w.writeheader()
        w.writerows(compounds)

    # ---- console summary ----
    print(f"parsed ok: {len(compounds)}; skipped/bad: {len(bad)}")
    print(f"activity rows (vs NPT535/NPT3812): {len(all_act)}")
    print(f"species rows: {len(all_species)}")
    n_pth = len({r['npc_id'] for r in all_act if r['gene'] == 'PTH1R'})
    n_slc = len({r['npc_id'] for r in all_act if r['gene'] == 'SLCO4C1'})
    print(f"compounds with PTH1R activity: {n_pth}; with SLCO4C1 activity: {n_slc}")
    if bad:
        print("bad pages:", ",".join(bad[:20]))


if __name__ == "__main__":
    sys.exit(main())
