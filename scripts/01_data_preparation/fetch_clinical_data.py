import os

PROJECT_ROOT = os.environ.get("HCC_PROJECT_ROOT",
                              os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import urllib.request, json, csv

# Fetch TCGA-LIHC survival data from the GDC cases API
url = 'https://api.gdc.cancer.gov/cases'
payload = json.dumps({
    "filters": {"project": {"project_id": ["TCGA-LIHC"]}},
    "fields": "submitter_id,demographic.gender,diagnoses.vital_status,diagnoses.days_to_death,diagnoses.days_to_last_follow_up,diagnoses.days_to_birth",
    "size": 1000
}).encode()
req = urllib.request.Request(url, data=payload, headers={'Content-Type': 'application/json'})
r = urllib.request.urlopen(req, timeout=60)
data = json.loads(r.read().decode())
hits = data['data']['hits']
print('number of cases:', len(hits))

rows = []
for h in hits:
    pid = h.get('submitter_id', '')
    diag = h.get('diagnoses', [{}])
    if diag:
        d = diag[0]
    else:
        d = {}
    vital = d.get('vital_status', '')
    dtd = d.get('days_to_death', None)
    dtlf = d.get('days_to_last_follow_up', None)
    # OS.time
    if vital == 'Dead' and dtd is not None:
        os_time = dtd
        os_status = 1
    elif vital == 'Alive' and dtlf is not None:
        os_time = dtlf
        os_status = 0
    else:
        os_time = None
        os_status = None
    rows.append({'Patient': pid, 'vital_status': vital,
                 'OS_time': os_time, 'OS': os_status})

# Write CSV
with open(os.path.join(PROJECT_ROOT, "work/main/tcga_clinical_os.csv"), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['Patient', 'vital_status', 'OS_time', 'OS'])
    w.writeheader()
    w.writerows(rows)

n_os = sum(1 for r in rows if r['OS'] is not None)
print('cases with OS information:', n_os, '/', len(rows))
print('example rows:', rows[:3])
