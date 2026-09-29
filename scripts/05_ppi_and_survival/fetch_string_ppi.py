import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv, http.client, urllib.parse, time

genes = [r[0] for r in csv.reader(open(os.path.join(PROJECT_ROOT, "work/main/candidate_pool_final.csv"), encoding='utf-8'))][1:]
genes = [g for g in genes if g.strip()]
print('number of candidate genes:', len(genes))

params = urllib.parse.urlencode({
    'identifiers': '\n'.join(genes),
    'species': '9606',
    'required_score': '400'
}).encode()

def fetch(timeout=600):
    conn = http.client.HTTPSConnection('string-db.org', timeout=timeout)
    conn.request('POST', '/api/tsv/network', body=params,
                 headers={'User-Agent': 'Mozilla/5.0',
                          'Content-Type': 'application/x-www-form-urlencoded'})
    resp = conn.getresponse()
    chunks = []
    while True:
        chunk = resp.read(65536)
        if not chunk:
            break
        chunks.append(chunk)
    data = b''.join(chunks).decode('utf-8')
    conn.close()
    return data

t0 = time.time()
try:
    data = fetch()
except Exception as e:
    print('first attempt failed:', e, ', retrying...')
    time.sleep(3)
    data = fetch()

print('elapsed %.1f s' % (time.time() - t0))
lines = data.strip().split('\n')
open(os.path.join(PROJECT_ROOT, "work/main/ppi_string.tsv"), 'w', encoding='utf-8').write(data)
print('total PPI lines (including the header):', len(lines))
print('header:', lines[0][:150])
print('number of edges:', len(lines) - 1)
