import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv, os, hashlib, urllib.request, sys

manifest = os.path.join(PROJECT_ROOT, "archive/network_pharmacology_hcc/01_Data/Raw/GDCdata/manifest.csv")
base_dir = os.path.join(PROJECT_ROOT, "archive/network_pharmacology_hcc/01_Data/Raw/GDCdata/TCGA-LIHC/Transcriptome_Profiling/Gene_Expression_Quantification")

rows = list(csv.reader(open(manifest, encoding='utf-8')))[1:]
total = len(rows)
done = 0
fail = 0

def md5sum(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()

for i, r in enumerate(rows):
    fid, fname, fsize, md5, barcode, stype = r
    dest_dir = os.path.join(base_dir, fid)
    dest = os.path.join(dest_dir, fname)
    if os.path.exists(dest) and md5 and md5sum(dest) == md5:
        done += 1
        continue
    os.makedirs(dest_dir, exist_ok=True)
    url = 'https://api.gdc.cancer.gov/data/' + fid
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=300) as resp, open(dest, 'wb') as f:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        if md5 and md5sum(dest) != md5:
            os.remove(dest)
            fail += 1
            print('[%d/%d] md5 mismatch %s' % (i + 1, total, fid), flush=True)
            continue
        done += 1
        if done % 20 == 0:
            print('progress: %d/%d' % (done, total), flush=True)
    except Exception as e:
        fail += 1
        if os.path.exists(dest):
            os.remove(dest)
        print('[%d/%d] fail %s: %s' % (i + 1, total, fid, e), flush=True)

print('DONE completed=%d/%d failed=%d' % (done, total, fail), flush=True)
