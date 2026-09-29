# -*- coding: utf-8 -*-
"""NPASS TCM active-compound screening: natural products targeting PTH1R / SLCO4C1 and their source species.
Input: G:/hcc_drug/herb/NPASS3.0_*.txt (all files must be downloaded in full)
Output: npass_tcm_hits.csv (hit compounds), npass_tcm_herbs.csv (source medicinal materials summary)
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv, os

D = os.path.join(PROJECT_ROOT, "work/drug/herb")
OUT = os.path.join(PROJECT_ROOT, "work/drug")
TARGETS = {"Q03431": "PTH1R", "Q6ZQN7": "SLCO4C1"}

def sniff(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.readline().rstrip("\n").split("\t")

# 1. Target table: target_id -> gene
tcols = sniff(f"{D}/NPASS3.0_target.txt")
print("target cols:", tcols)
tid2gene = {}
with open(f"{D}/NPASS3.0_target.txt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        joined = "\t".join(row.values())
        for acc, gene in TARGETS.items():
            if acc in joined or gene in joined.upper():
                tid2gene[row[tcols[0]]] = gene
print("target IDs hit:", tid2gene)

# 2. Activity table: select natural products acting on the two targets
acols = sniff(f"{D}/NPASS3.0_activities.txt")
print("activity cols:", acols)
hits = {}  # (gene, np_id) -> best activity row
with open(f"{D}/NPASS3.0_activities.txt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        gene = tid2gene.get(row.get("target_id", ""))
        if not gene:
            continue
        key = (gene, row["np_id"])
        hits.setdefault(key, row)  # collect everything first, sort later
print(f"activity hits (gene, np_id) pairs: {len(hits)}")

# 3. Compound information
gcols = sniff(f"{D}/NPASS3.0_naturalproducts_generalinfo.txt")
np_info = {}
with open(f"{D}/NPASS3.0_naturalproducts_generalinfo.txt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        np_info[row["np_id"]] = row

# 4. Structure
scols = sniff(f"{D}/NPASS3.0_naturalproducts_structure.txt")
print("structure cols:", scols)
np_smi = {}
need = {k[1] for k in hits}
with open(f"{D}/NPASS3.0_naturalproducts_structure.txt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["np_id"] in need:
            np_smi[row["np_id"]] = row.get("canonical_smiles") or row.get("smiles") or ""

# 5. Species source
pcols = sniff(f"{D}/NPASS3.0_naturalproducts_species_pair.txt")
spcols = sniff(f"{D}/NPASS3.0_species_info.txt")
print("species_pair cols:", pcols, "| species_info cols:", spcols)
sp_name = {}
with open(f"{D}/NPASS3.0_species_info.txt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        sp_name[row[spcols[0]]] = row
np2sp = {}
with open(f"{D}/NPASS3.0_naturalproducts_species_pair.txt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["np_id"] in need:
            np2sp.setdefault(row["np_id"], set()).add(row[pcols[1]])

# 6. Aggregated output
rows_out = []
for (gene, npid), a in hits.items():
    info = np_info.get(npid, {})
    species = []
    for sid in np2sp.get(npid, ()):  # species names (Latin + common name)
        sp = sp_name.get(sid, {})
        nm = sp.get("organism_name") or sp.get("species_name") or ""
        if nm:
            species.append(nm)
    rows_out.append({
        "target": gene,
        "np_id": npid,
        "compound": info.get("pref_name") or info.get("name_initial") or "",
        "pubchem_id": info.get("pubchem_id", ""),
        "chembl_id": info.get("chembl_id", ""),
        "activity_type": a.get("activity_type", ""),
        "activity_value": a.get("activity_value", ""),
        "activity_unit": a.get("activity_unit", ""),
        "activity_relation": a.get("activity_relation", ""),
        "assay": a.get("assay_description", "")[:120],
        "ref": a.get("ref_id", ""),
        "smiles": np_smi.get(npid, ""),
        "n_species": len(species),
        "species": "; ".join(sorted(species))[:400],
    })
rows_out.sort(key=lambda r: (r["target"], -(float(r["activity_value"]) if r["activity_value"].replace(".", "", 1).isdigit() else -1)))
with open(f"{OUT}/npass_tcm_hits.csv", "w", newline="", encoding="utf-8") as f:
    if rows_out:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader(); w.writerows(rows_out)
print(f"\nwrote {len(rows_out)} hit compounds")
for g in TARGETS.values():
    rr = [r for r in rows_out if r["target"] == g]
    print(f"\n[{g}] {len(rr)} natural products")
    for r in rr[:25]:
        print(f"  {r['compound'][:34]:34s} {r['activity_type']:6s} {r['activity_relation']:2s}{r['activity_value']:>10s} {r['activity_unit']:6s} species={r['n_species']}  e.g. {r['species'][:70]}")
