# -*- coding: utf-8 -*-
"""
build_tcm_hit_table.py
Merge NPASS parse outputs into the final TCM active-ingredient hit table.

Inputs (G:/hcc_drug/):
  npass_target_compounds.csv  - 246 target-page rows (243 PTH1R + 3 SLCO4C1, 245 unique NPC)
  npass_tcm_activity.csv      - quantitative activity rows parsed from compound pages
  npass_tcm_species.csv       - species-source rows parsed from compound pages

Outputs (G:/hcc_drug/):
  npass_tcm_hits.csv          - one row per compound-target pair, with best potency,
                                activity flag, species count, curated TCM herb annotation, SMILES
  npass_tcm_herb_pivot.csv    - one row per curated TCM herb: compounds contained, best potency
"""
import csv
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from collections import defaultdict

D = os.path.join(PROJECT_ROOT, "work/drug")

# curated species -> (Chinese herb name, Latin binomial as appears in NPASS may vary -> substring match)
TCM_SPECIES = {
    "Digitalis purpurea": "Digitalis purpurea (Maodihuang)", "Digitalis lanata": "Digitalis lanata",
    "Agastache rugosa": "Agastache rugosa (Huoxiang)", "Agastache rugosus": "Agastache rugosa (Huoxiang)",
    "Illicium verum": "Illicium verum (star anise)",
    "Baphicacanthus cusia": "Baphicacanthus cusia (Nanbanlangen)", "Isatis indigotica": "Isatis indigotica (Banlangen)",
    "Bufo gargarizans": "Bufo gargarizans (Chansu)", "Bufo bufo": "Bufo bufo (Chansu)",
    "Magnolia officinalis": "Magnolia officinalis (Houpo)", "Magnolia obovata": "Magnolia obovata (Japanese houpo)",
    "Nelumbo nucifera": "Nelumbo nucifera (lotus)",
    "Cucumis melo": "Cucumis melo (Guadi)", "Ecballium elaterium": "Ecballium elaterium (squirting cucumber)",
    "Helicteres isora": "Helicteres isora (Shanzhima)",
    "Glycyrrhiza uralensis": "Glycyrrhiza uralensis (Gancao)", "Glycyrrhiza glabra": "Glycyrrhiza glabra",
    "Glycyrrhiza inflata": "Glycyrrhiza inflata", "Glycyrrhiza aspera": "Glycyrrhiza aspera",
    "Ginkgo biloba": "Ginkgo biloba (Yinxing)",
    "Hydrangea macrophylla": "Hydrangea macrophylla (Gancha)", "Dichroa febrifuga": "Dichroa febrifuga (Changshan)",
    "Eupatorium": "Eupatorium genus (Zelan / Peilan)",
    "Garcinia mangostana": "Garcinia mangostana (mangosteen)", "Garcinia": "Garcinia genus",
    "Cotinus coggygria": "Cotinus coggygria (Huanglu)",
    "Larrea tridentata": "Larrea tridentata (creosote bush)",
    "Achillea millefolium": "Yarrow (Achillea millefolium)", "Achillea": "Achillea genus",
    "Linaria vulgaris": "Linaria vulgaris (common toadflax)",
    "Cirsium": "Cirsium genus (Daji / Xiaoji)",
    "Astragalus membranaceus": "Astragalus membranaceus (Huangqi)", "Astragalus": "Astragalus genus",
    "Scutellaria baicalensis": "Scutellaria baicalensis (Huangqin)", "Scutellaria": "Scutellaria genus",
    "Coptis chinensis": "Coptis chinensis (Huanglian)", "Phellodendron": "Phellodendron genus (Huangbai)",
    "Sophora flavescens": "Sophora flavescens (Kushen)", "Sophora": "Sophora genus",
    "Lonicera japonica": "Lonicera japonica (Jinyinhua)", "Forsythia suspensa": "Forsythia suspensa (Lianqiao)",
    "Houttuynia cordata": "Houttuynia cordata (Yuxingcao)", "Taraxacum": "Taraxacum genus (Dandelion)",
    "Angelica sinensis": "Angelica sinensis (Danggui)", "Angelica dahurica": "Angelica dahurica (Baizhi)", "Angelica gigas": "Angelica gigas (Korean angelica)",
    "Angelica acutiloba": "Angelica acutiloba (Dongdanggui)",
    "Ligusticum chuanxiong": "Ligusticum chuanxiong (Chuanxiong)", "Salvia miltiorrhiza": "Salvia miltiorrhiza (Danshen)",
    "Paeonia lactiflora": "Paeonia lactiflora (Baishao)", "Rehmannia glutinosa": "Rehmannia glutinosa (Dihuang)",
    "Rheum palmatum": "Rheum palmatum (Zhangye dahuang)", "Rheum officinale": "Rheum officinale (Yaoyong dahuang)",
    "Aconitum": "Aconitum genus (Fuzi / Chuanwu)",
    "Cinchona": "Cinchona genus",
    "Camptotheca acuminata": "Camptotheca acuminata (Xishu)", "Cephalotaxus": "Cephalotaxus genus",
    "Taxus": "Taxus genus (yew)", "Podophyllum": "Podophyllum genus",
    "Catharanthus roseus": "Catharanthus roseus (Madagascar periwinkle)",
    "Tripterygium wilfordii": "Tripterygium wilfordii (Leigongteng)",
    "Artemisia annua": "Artemisia annua (Qinghao)", "Artemisia apiacea": "Artemisia apiacea (Qinghao)",
    "Andrographis paniculata": "Andrographis paniculata (Chuanxinlian)",
    "Schisandra chinensis": "Schisandra chinensis (Wuweizi)",
    "Ganoderma lucidum": "Ganoderma lucidum (Lingzhi)", "Ganoderma": "Ganoderma genus",
    "Zingiber officinale": "Zingiber officinale (ginger)", "Curcuma longa": "Curcuma longa (turmeric)",
    "Allium sativum": "Garlic", "Piper nigrum": "Black pepper", "Piper longum": "Long pepper",
    "Mentha": "Mentha genus (mint)", "Perilla frutescens": "Perilla frutescens (Zisu)",
    "Camellia sinensis": "Tea plant (Camellia sinensis)",
    "Moringa oleifera": "Moringa oleifera", "Hemerocallis fulva": "Hemerocallis fulva (daylily)",
    "Withania somnifera": "Withania somnifera (ashwagandha)", "Tinospora cordifolia": "Tinospora cordifolia",
    "Swertia chirayita": "Swertia chirayita",
    "Aconitum heterophylloides": "Aconitum heterophylloides",
    "Gypsophila": "Gypsophila genus",
    "Arnica montana": "Arnica montana",
    "Elephantopus": "Elephantopus genus",
    "Aglaia": "Aglaia genus",
    "Streptomyces": "Streptomyces genus (microbial source)", "Micromonospora": "Micromonospora genus (microbial source)",
    "Aspergillus": "Aspergillus genus (microbial source)", "Penicillium": "Penicillium genus (microbial source)",
}


def load(path):
    with open(path, encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def main():
    tgt = load(os.path.join(D, "npass_target_compounds.csv"))
    act = load(os.path.join(D, "npass_tcm_activity.csv"))
    sp = load(os.path.join(D, "npass_tcm_species.csv"))
    comp = {r["npc_id"]: r for r in load(os.path.join(D, "npass_tcm_compounds.csv"))}

    # species per compound
    sp_map = defaultdict(set)
    for r in sp:
        if r["organism_name"] and r["organism_name"] != "n.a.":
            sp_map[r["npc_id"]].add(r["organism_name"])

    # best potency per (npc, gene)
    act_map = defaultdict(list)
    for r in act:
        act_map[(r["npc_id"], r["gene"])].append(r)

    def best_row(rows):
        def k(r):
            try:
                return (0, float(r["activity_value"]))
            except ValueError:
                return (1, 1e18)
        return sorted(rows, key=k)[0]

    def tcm_hits(npc_id):
        hits = {}
        for s in sp_map.get(npc_id, ()):  # s = organism name
            for key, herb in TCM_SPECIES.items():
                if key.lower() in s.lower():
                    hits.setdefault(herb, s)
        return hits

    rows = []
    for t in tgt:
        npc, gene = t["npc_id"], t["target"]
        rows_act = act_map.get((npc, gene), [])
        b = best_row(rows_act) if rows_act else None
        herbs = tcm_hits(npc)
        c = comp.get(npc, {})
        rows.append({
            "npc_id": npc,
            "compound_name": c.get("compound_name", "") or t.get("compound", ""),
            "gene": gene,
            "target_id": b["target_id"] if b else ("NPT535" if gene == "PTH1R" else "NPT3812"),
            "has_quantitative_activity": "yes" if b else "no (target-page association only)",
            "activity_type": b["activity_type"] if b else "",
            "activity_relation": b["activity_relation"] if b else "",
            "potency_nM": b["activity_value"] if b else "",
            "n_species": len(sp_map.get(npc, ())),
            "tcm_herbs": "; ".join(sorted(herbs)),
            "example_species": "; ".join(sorted(herbs.values())[:4]) or "; ".join(sorted(sp_map.get(npc, ()))[:4]),
            "smiles": c.get("smiles", ""),
        })

    def sort_key(r):
        try:
            return (0, float(r["potency_nM"]))
        except ValueError:
            return (1, 1e18)
    rows.sort(key=lambda r: (r["gene"] != "PTH1R", sort_key(r)))

    with open(os.path.join(D, "npass_tcm_hits.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # herb pivot: herb -> compounds (best potency across)
    pivot = defaultdict(list)
    for r in rows:
        for herb in [h.strip() for h in r["tcm_herbs"].split(";") if h.strip()]:
            pivot[herb].append(r)
    prow = []
    for herb, rs in pivot.items():
        rs2 = sorted(rs, key=sort_key)
        best = rs2[0]
        prow.append({
            "tcm_herb": herb,
            "n_compounds": len(rs),
            "compounds": "; ".join(sorted({x["compound_name"] for x in rs if x["compound_name"]})[:8]),
            "genes": "; ".join(sorted({x["gene"] for x in rs})),
            "best_compound": best["compound_name"],
            "best_gene": best["gene"],
            "best_potency_nM": best["potency_nM"],
        })
    prow.sort(key=lambda r: (sort_key({"potency_nM": r["best_potency_nM"]}), -r["n_compounds"]))
    with open(os.path.join(D, "npass_tcm_herb_pivot.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(prow[0].keys()))
        w.writeheader()
        w.writerows(prow)

    n_quant = sum(1 for r in rows if r["potency_nM"])
    n_herb = sum(1 for r in rows if r["tcm_herbs"])
    print(f"hit rows: {len(rows)} (with quantitative potency: {n_quant}; with curated TCM herb: {n_herb})")
    print(f"herb pivot rows: {len(prow)}")
    print("\n--- top 20 herbs by best compound potency ---")
    for p in prow[:20]:
        print(f"{p['tcm_herb']:22s} n={p['n_compounds']:>2d}  best={p['best_compound'][:30]:32s} {p['best_gene']:8s} {p['best_potency_nM']:>10s}")


if __name__ == "__main__":
    main()
