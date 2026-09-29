# -*- coding: utf-8 -*-
"""parse_liri.py — parse locally downloaded LIRI-JP files into liri_allgenes.csv"""
import gzip
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv
from collections import Counter

BASE = os.path.join(PROJECT_ROOT, "work/val")
LIRI = f"{BASE}/liri/release_28/data/LIRI-JP"
GENES = {"ANGPTL6", "CFHR3", "COL15A1", "CPEB3", "ECM1", "EDIL3", "ESM1", "FAM83F", "GABRD", "MARCO", "MXD3", "PCDH17", "PKN3", "PTH1R", "PZP", "TERT", "NTF3", "COLEC10", "SLCO4C1"}


def parse_donor_files():
    donors = {}
    bad = 0
    for root, _, files in os.walk(LIRI):
        if "/donor" not in root.replace("\\", "/"):
            continue
        for fn in files:
            if not fn.endswith(".gz"):
                continue
            try:
                txt = gzip.open(os.path.join(root, fn), "rt", encoding="utf-8", errors="replace").read()
            except (OSError, EOFError):
                bad += 1
                continue
            for line in txt.splitlines():
                if not line.strip():
                    continue
                p = line.split("\t")
                if len(p) < 17:
                    continue
                donors[p[0]] = (p[5], p[16])  # vital_status, survival_time
    print(f"corrupted donor files: {bad}")
    return donors


def parse_exp_files():
    exprs = {}
    bad = 0
    for root, _, files in os.walk(LIRI):
        if "/exp_seq" not in root.replace("\\", "/"):
            continue
        for fn in files:
            if not fn.endswith(".gz"):
                continue
            try:
                txt = gzip.open(os.path.join(root, fn), "rt", encoding="utf-8", errors="replace").read()
            except (OSError, EOFError):
                bad += 1
                continue
            for line in txt.splitlines():
                if not line.strip():
                    continue
                p = line.split("\t")
                if len(p) < 10:
                    continue
                did = p[0]
                gene = p[7].split("_")[0]
                if gene not in GENES:
                    continue
                try:
                    val = float(p[8])
                except ValueError:
                    continue
                exprs.setdefault(did, {})[gene] = val
    print(f"corrupted exp_seq files: {bad}")
    return exprs


def main():
    donors = parse_donor_files()
    exprs = parse_exp_files()
    print(f"donor survival records: {len(donors)}")
    print(f"donor expression records: {len(exprs)}")
    rows = []
    for did, (vital, surv) in donors.items():
        if did not in exprs:
            continue
        if not surv or not surv.strip().isdigit():
            continue
        if vital not in ("alive", "deceased"):
            continue
        rec = {"donor_id": did, "OS_status": "1" if vital == "deceased" else "0",
               "OS_days": int(surv)}
        ok = True
        for g in sorted(GENES):
            v = exprs[did].get(g)
            if v is None:
                ok = False
            rec[g] = v if v is not None else ""
        rows.append(rec)
    with open(f"{BASE}/liri_allgenes.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nmerged samples: {len(rows)}")
    print("OS_status:", Counter(r["OS_status"] for r in rows))
    if rows:
        days = [r["OS_days"] for r in rows]
        print(f"OS_days range: {min(days)}-{max(days)}")


if __name__ == "__main__":
    main()
