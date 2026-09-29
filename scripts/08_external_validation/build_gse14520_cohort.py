# -*- coding: utf-8 -*-
"""
build_gse14520_cohort.py
GSE14520 (GPL3921) external validation cohort:
  series matrix (expr) + GPL3921 annot (probe->symbol) + Extra_Supplement (survival)
Output: gse14520_allgenes.csv  (sample, OS_months, OS_event, 7 genes expression)
"""
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import gzip
import csv

BASE = os.path.join(PROJECT_ROOT, "work/val")
GENES = ["PTH1R","SLCO4C1","COL15A1","NTF3","PZP","CPEB3","COLEC10","CFHR3","MXD3","TERT","MARCO","ANGPTL6","ECM1","EDIL3","ESM1","FAM83F","GABRD","PCDH17","PKN3"]


def parse_annot():
    d = gzip.open(f"{BASE}/GPL3921.annot.gz", "rt", encoding="utf-8", errors="replace").read()
    lines = d.split("\n")
    probe2gene = {}
    hdr = None
    sym_idx = None
    for l in lines:
        if l.startswith("ID\t"):
            hdr = l.split("\t")
            sym_idx = hdr.index("Gene symbol")
            continue
        if hdr is None or not l.strip() or l.startswith(("#", "!", "^")):
            continue
        parts = l.split("\t")
        if len(parts) <= sym_idx:
            continue
        probe = parts[0].strip().strip('"')
        sym = parts[sym_idx].strip().strip('"')
        if probe and sym and sym != "---":
            probe2gene.setdefault(probe, sym)
    return probe2gene


def parse_series_matrix(probe2gene):
    """return {gene: {gsm: expr}} for the 7 genes (mean over probes)"""
    d = gzip.open(f"{BASE}/GSE14520_GPL3921_matrix.txt.gz", "rt", encoding="utf-8", errors="replace").read()
    lines = d.split("\n")
    # find table header
    tbl = None
    for i, l in enumerate(lines):
        if l.startswith("!series_matrix_table_begin"):
            tbl = i + 1
            break
    hdr = [h.strip().strip('"') for h in lines[tbl].split("\t")]
    samples = hdr[1:]  # GSM ids
    # collect probes for our genes
    gene_probes = {g: [] for g in GENES}
    probe_expr = {}  # probe -> list aligned to samples
    for l in lines[tbl + 1:]:
        if not l.strip() or l.startswith("!"):
            continue
        parts = l.split("\t")
        probe = parts[0].strip().strip('"')
        sym = probe2gene.get(probe)
        if sym in GENES:
            gene_probes[sym].append(probe)
            vals = []
            for v in parts[1:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(float("nan"))
            if len(vals) == len(samples):
                probe_expr[probe] = vals
    # mean over probes per gene
    gene_expr = {g: {} for g in GENES}
    for g in GENES:
        probes = gene_probes[g]
        if not probes:
            continue
        for j, gsm in enumerate(samples):
            vals = [probe_expr[p][j] for p in probes if p in probe_expr]
            gene_expr[g][gsm] = sum(vals) / len(vals) if vals else float("nan")
    return gene_expr, samples


def parse_survival():
    d = gzip.open(f"{BASE}/GSE14520_Extra_Supplement.txt.gz", "rt", encoding="utf-8", errors="replace").read()
    lines = d.split("\n")
    hdr = lines[0].split("\t")
    idx = {h.strip(): i for i, h in enumerate(hdr)}
    surv = {}
    for l in lines[1:]:
        if not l.strip():
            continue
        parts = l.split("\t")
        if len(parts) < len(hdr):
            continue
        gsm = parts[idx["Affy_GSM"]].strip()
        tissue = parts[idx["Tissue Type"]].strip()
        try:
            os_months = float(parts[idx["Survival months"]].strip() or "nan")
        except ValueError:
            os_months = float("nan")
        status = parts[idx["Survival status"]].strip()
        surv[gsm] = {"tissue": tissue, "os_months": os_months, "os_status": status}
    return surv


def main():
    probe2gene = parse_annot()
    print("number of probes in the annotation map:", len(probe2gene))
    gene_expr, samples = parse_series_matrix(probe2gene)
    for g in GENES:
        print(f"  {g}: {len(gene_expr[g])} samples with expression")
    surv = parse_survival()
    print("number of survival records:", len(surv))

    # match tumor samples with survival + expression
    rows = []
    for gsm in samples:
        s = surv.get(gsm)
        if not s or s["tissue"] != "Tumor":
            continue
        if s["os_months"] != s["os_months"] or s["os_status"] == "":
            continue
        rec = {"GSM": gsm, "OS_months": s["os_months"], "OS_status": s["os_status"]}
        # per-gene tolerant: keep the sample even if some genes are absent from the platform
        for g in GENES:
            v = gene_expr[g].get(gsm, float("nan"))
            rec[g] = "" if v != v else v
        rows.append(rec)

    with open(f"{BASE}/gse14520_allgenes.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    print(f"\ntumour samples (with survival + expression): {len(rows)}")
    print("OS_status values:", Counter(r["OS_status"] for r in rows))
    import statistics
    m = [r["OS_months"] for r in rows]
    print(f"OS_months: range {min(m):.1f}-{max(m):.1f}, median {statistics.median(m):.1f}")


if __name__ == "__main__":
    main()
