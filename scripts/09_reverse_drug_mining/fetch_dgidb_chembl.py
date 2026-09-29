# -*- coding: utf-8 -*-
"""Reverse network pharmacology, step 1: query known ligands and drugs for PTH1R / SLCO4C1.
Data sources: DGIdb v5 (GraphQL) + ChEMBL REST API. HTTP requests go through a curl subprocess (the host proxy only allows curl).
Output: G:/hcc_drug/dgidb_interactions.csv, chembl_ligands.csv
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json, csv, time, subprocess, os

OUT = os.path.join(PROJECT_ROOT, "work/drug")
os.makedirs(OUT, exist_ok=True)
ENV = dict(os.environ, http_proxy="http://127.0.0.1:53448",
           https_proxy="http://127.0.0.1:53448")

def curl(args, tries=3):
    for i in range(tries):
        try:
            r = subprocess.run(["curl", "-s", "--max-time", "90"] + args,
                               capture_output=True, env=ENV, timeout=120)
            if r.returncode == 0 and r.stdout:
                return r.stdout.decode("utf-8", "replace")
        except Exception as e:
            print("  curl err:", e)
        print("  retry", i + 1)
        time.sleep(2 + i * 3)
    return None

def get_json(url):
    txt = curl([url])
    return json.loads(txt) if txt else None

def post_json(url, payload):
    txt = curl(["-X", "POST", url, "-H", "Content-Type: application/json",
                "-d", json.dumps(payload)])
    return json.loads(txt) if txt else None

# ---------- 1. DGIdb ----------
print("== DGIdb ==")
Q = ('{ genes(names: ["PTH1R", "SLCO4C1"]) { nodes { name '
     'interactions { drug { name approved conceptId } '
     'interactionTypes { type directionality } sources { sourceDbName } } } } }')
d = post_json("https://dgidb.org/api/graphql", {"query": Q})
dg_rows = []
if d and d.get("data"):
    for node in d["data"]["genes"]["nodes"]:
        gene = node["name"]
        for it in node["interactions"]:
            types = ";".join(sorted({t.get("type") or "" for t in it["interactionTypes"] if t.get("type")}))
            dirs = ";".join(sorted({t.get("directionality") or "" for t in it["interactionTypes"] if t.get("directionality")}))
            srcs = ";".join(sorted({s["sourceDbName"] for s in it["sources"]}))
            dg_rows.append({
                "gene": gene,
                "drug": it["drug"]["name"],
                "approved": it["drug"]["approved"],
                "concept_id": it["drug"].get("conceptId", ""),
                "interaction_type": types,
                "directionality": dirs,
                "sources": srcs,
                "n_sources": len(it["sources"]),
            })
else:
    print("  DGIdb returned an error:", str(d)[:300])

if dg_rows:
    dg_rows.sort(key=lambda r: (r["gene"], not r["approved"], -r["n_sources"]))
    with open(OUT + "/dgidb_interactions.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(dg_rows[0].keys()))
        w.writeheader(); w.writerows(dg_rows)
    for g in ("PTH1R", "SLCO4C1"):
        rows = [r for r in dg_rows if r["gene"] == g]
        app = [r for r in rows if r["approved"]]
        print(f"  {g}: {len(rows)} interactions, {len(app)} approved")

# ---------- 2. ChEMBL ----------
print("== ChEMBL ==")
TARGETS = {"PTH1R": "CHEMBL1793", "SLCO4C1": "CHEMBL2073690"}
acts = []
for gene, tid in TARGETS.items():
    d = get_json("https://www.ebi.ac.uk/chembl/api/data/activity.json"
                 f"?target_chembl_id={tid}&limit=1000")
    page = d["activities"] if d else []
    total = d["page_meta"]["total_count"] if d else 0
    print(f"  {gene} ({tid}): {total} activities, fetched {len(page)}")
    for a in page:
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
            "document_chembl_id": a.get("document_chembl_id", ""),
        })

mol_ids = sorted({a["molecule_chembl_id"] for a in acts})
print(f"  fetching {len(mol_ids)} molecule details...")
mol_info = {}
for i, mid in enumerate(mol_ids):
    d = get_json(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{mid}.json")
    if d:
        ms = d.get("molecule_structures") or {}
        mol_info[mid] = {
            "pref_name": d.get("pref_name") or "",
            "max_phase": d.get("max_phase") if d.get("max_phase") is not None else "",
            "molecule_type": d.get("molecule_type") or "",
            "canonical_smiles": ms.get("canonical_smiles") or "",
            "first_approval": d.get("first_approval") or "",
        }
    if (i + 1) % 20 == 0:
        print(f"    {i+1}/{len(mol_ids)}")

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
if rows:
    with open(OUT + "/chembl_ligands.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    for g in ("PTH1R", "SLCO4C1"):
        rr = [r for r in rows if r["gene"] == g]
        appr = [r for r in rr if str(r["max_phase"]) == "4.0" or str(r["max_phase"]) == "4"]
        clin = [r for r in rr if str(r["max_phase"]) in ("1.0", "2.0", "3.0", "1", "2", "3")]
        print(f"  {g}: {len(rr)} ligands (with pChEMBL), {len(appr)} approved, {len(clin)} clinical")
print("DONE")
