# -*- coding: utf-8 -*-
"""
download_liri.py — fetch ICGC LIRI-JP (release_28, open bucket) expression + survival,
extract the 7 genes and output liri_validation.csv
"""
import gzip
import json
import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import re
import sys
import urllib.request
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO

# explicit proxy (urllib does not honor http_proxy env)
_proxy = urllib.request.ProxyHandler({"http": "http://127.0.0.1:53448",
                                      "https": "http://127.0.0.1:53448"})
_opener = urllib.request.build_opener(_proxy)
urllib.request.install_opener(_opener)

BASE = os.path.join(PROJECT_ROOT, "work/val")
HOST = "https://object.genomeinformatics.org/icgc25k-open"
PREFIX = "release_28/data/LIRI-JP/"
GENES = {"PTH1R", "SLCO4C1", "COL15A1", "NTF3", "PZP", "CPEB3", "COLEC10"}


def list_objects(prefix):
    keys = []
    token = None
    while True:
        url = f"{HOST}?list-type=2&prefix={prefix}"
        if token:
            url += f"&continuation-token={urllib.parse.quote(token)}"
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read().decode("utf-8")
        for m in re.finditer(r"<Key>([^<]+)</Key>", data):
            keys.append(m.group(1))
        mt = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", data)
        if not mt:
            break
        token = mt.group(1)
    return keys


def fetch(key):
    url = f"{HOST}/{key}"
    with urllib.request.urlopen(url, timeout=90) as r:
        return gzip.open(BytesIO(r.read()), "rt", encoding="utf-8", errors="replace").read()


def parse_donor(text):
    """return dict donor_id -> (vital_status, survival_time)"""
    out = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 17:
            continue
        did = p[0]
        vital = p[5]
        surv = p[16]
        out[did] = (vital, surv)
    return out


def parse_exp(text):
    """return dict donor_id -> {gene: norm_count}; gene_id at col 8, norm count at col 9"""
    out = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 10:
            continue
        did = p[0]
        gene = p[7].split("_")[0]  # strip RefSeq suffix
        try:
            val = float(p[8])
        except ValueError:
            continue
        if gene in GENES:
            out.setdefault(did, {})[gene] = val
    return out


def main():
    print("Listing LIRI-JP files...")
    keys = list_objects(PREFIX)
    donor_keys = [k for k in keys if "/donor/" in k and k.endswith(".gz")]
    exp_keys = [k for k in keys if "/exp_seq/" in k and k.endswith(".gz")]
    print(f"donor files: {len(donor_keys)}, exp_seq files: {len(exp_keys)}")

    with ThreadPoolExecutor(max_workers=8) as ex:
        donor_texts = list(ex.map(fetch, donor_keys))
        exp_texts = list(ex.map(fetch, exp_keys))

    donors = {}
    for t in donor_texts:
        donors.update(parse_donor(t))
    exprs = {}
    for t in exp_texts:
        for did, gd in parse_exp(t).items():
            exprs.setdefault(did, {}).update(gd)
    print(f"donor survival records: {len(donors)}, donor expression records: {len(exprs)}")

    # merge
    rows = []
    for did, (vital, surv) in donors.items():
        if did not in exprs:
            continue
        if not surv or not surv.isdigit():
            continue
        if vital not in ("alive", "deceased"):
            continue
        rec = {"donor_id": did, "OS_status": "1" if vital == "deceased" else "0",
               "OS_days": int(surv)}
        ok = True
        for g in GENES:
            v = exprs[did].get(g)
            if v is None:
                ok = False
            rec[g] = v if v is not None else ""
        if ok:
            rows.append(rec)

    import csv
    with open(f"{BASE}/liri_validation.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(f"\nmerged samples: {len(rows)}")
    print("OS_status:", Counter(r["OS_status"] for r in rows))
    days = [r["OS_days"] for r in rows]
    print(f"OS_days range: {min(days)}-{max(days)}")


if __name__ == "__main__":
    main()
