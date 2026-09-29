# -*- coding: utf-8 -*-
"""Parse raw DGIdb + ChEMBL JSON and build the reverse drug-mining result tables and summary."""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json, csv, glob

OUT = os.path.join(PROJECT_ROOT, "work/drug")

# ---- DGIdb ----
d = json.load(open(OUT + r"/raw/dgidb.json", encoding="utf-8"))
dg_rows = []
for node in d["data"]["genes"]["nodes"]:
    gene = node["name"]
    for it in node["interactions"]:
        types = ";".join(sorted({t.get("type") or "" for t in it["interactionTypes"] if t.get("type")}))
        dirs = ";".join(sorted({t.get("directionality") or "" for t in it["interactionTypes"] if t.get("directionality")}))
        srcs = ";".join(sorted({s["sourceDbName"] for s in it["sources"]}))
        dg_rows.append({
            "gene": gene, "drug": it["drug"]["name"],
            "approved": it["drug"]["approved"],
            "concept_id": it["drug"].get("conceptId", ""),
            "interaction_type": types, "directionality": dirs,
            "sources": srcs, "n_sources": len(it["sources"]),
        })
dg_rows.sort(key=lambda r: (r["gene"], not r["approved"], -r["n_sources"]))
with open(OUT + "/dgidb_interactions.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(dg_rows[0].keys()))
    w.writeheader(); w.writerows(dg_rows)
print(f"DGIdb: {len(dg_rows)} rows")

# ---- ChEMBL molecule details ----
mol_info = {}
for fp in glob.glob(OUT + r"/raw/mol_*.json"):
    d = json.load(open(fp, encoding="utf-8"))
    for m in d["molecules"]:
        ms = m.get("molecule_structures") or {}
        mol_info[m["molecule_chembl_id"]] = {
            "pref_name": m.get("pref_name") or "",
            "max_phase": m.get("max_phase") if m.get("max_phase") is not None else "",
            "molecule_type": m.get("molecule_type") or "",
            "canonical_smiles": ms.get("canonical_smiles") or "",
            "first_approval": m.get("first_approval") or "",
        }
print(f"molecule details: {len(mol_info)}")

# ---- ChEMBL activities (both targets merged, strongest record per molecule) ----
acts = []
for gene, fp in [("PTH1R", OUT + r"/raw/act_pth1r_full.json"),
                 ("SLCO4C1", OUT + r"/raw/act_slco4c1.json")]:
    d = json.load(open(fp, encoding="utf-8"))
    for a in d["activities"]:
        if a.get("pchembl_value") in (None, ""):
            continue
        acts.append({
            "gene": gene,
            "molecule_chembl_id": a["molecule_chembl_id"],
            "standard_type": a.get("standard_type", ""),
            "standard_value": a.get("standard_value", ""),
            "standard_units": a.get("standard_units", ""),
            "pchembl_value": a.get("pchembl_value", ""),
            "assay_type": a.get("assay_type", ""),
        })
best = {}
for a in acts:
    key = (a["gene"], a["molecule_chembl_id"])
    p = float(a["pchembl_value"])
    if key not in best or p > float(best[key]["pchembl_value"]):
        best[key] = a
rows = []
for (gene, mid), a in best.items():
    info = mol_info.get(mid, {})
    rows.append({**a, **{k: info.get(k, "") for k in
                 ("pref_name", "max_phase", "molecule_type", "canonical_smiles", "first_approval")}})

def phase_key(v):
    try:
        return -float(v)
    except (TypeError, ValueError):
        return 1
rows.sort(key=lambda r: (r["gene"], phase_key(r["max_phase"]), -float(r["pchembl_value"])))
with open(OUT + "/chembl_ligands.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"ChEMBL ligands: {len(rows)} rows")

# ---- Summary ----
print("\n===== Summary =====")
for g in ("PTH1R", "SLCO4C1"):
    dg = [r for r in dg_rows if r["gene"] == g]
    dg_app = [r for r in dg if r["approved"]]
    ch = [r for r in rows if r["gene"] == g]
    ch_app = [r for r in ch if str(r["max_phase"]) in ("4.0", "4")]
    print(f"\n[{g}] DGIdb {len(dg)} records ({len(dg_app)} approved) | ChEMBL {len(ch)} active molecules ({len(ch_app)} approved)")
    print("  DGIdb approved drugs:")
    for r in dg_app:
        print(f"    - {r['drug']}  [{r['interaction_type'] or 'NA'}]  source:{r['sources']}")
    print("  ChEMBL approved/clinical molecules:")
    for r in [r for r in ch if str(r['max_phase']) not in ('', '0.0', '0')]:
        print(f"    - {r['pref_name'] or r['molecule_chembl_id']}  phase={r['max_phase']}  pChEMBL={r['pchembl_value']} ({r['standard_type']})  {r['molecule_type']}")
