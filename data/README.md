# Data sources

**No data are redistributed in this repository.** All datasets used in the study are publicly
available and must be downloaded from their original providers. Placing them under the paths
below (relative to `$HCC_PROJECT_ROOT`) matches the defaults used by the scripts.

```
data/
├── raw_gdc/          # TCGA-LIHC STAR-Counts files + manifest.csv
└── intermediate/     # generated RDS objects (counts, VST matrix, batch metadata)
```

---

## 1. TCGA-LIHC — discovery cohort

| Item | Detail |
|---|---|
| Source | Genomic Data Commons (GDC) Data Portal — https://portal.gdc.cancer.gov/ |
| Project | TCGA-LIHC |
| Data type | Gene Expression Quantification (STAR-Counts) |
| Clinical | Clinical follow-up (vital status, days to death / last follow-up) |
| Samples used | 424 in total — 371 primary tumours, 3 recurrent tumours, 50 adjacent normal tissues |

Download with the GDC Data Transfer Tool or the `TCGAbiolinks` R package
(`scripts/01_data_preparation/download_tcga_expression.py` shows the file-list approach).
TCGA data are de-identified and openly accessible; no additional ethics approval or consent
is required for their use.

---

## 2. GSE14520 — external validation cohort

| Item | Detail |
|---|---|
| Source | NCBI GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520 |
| Platform | Affymetrix HT Human Genome U133A (GPL3921) |
| Samples used | 221 tumours, 85 deaths |

Files required:

- `GSE14520_GPL3921_matrix.txt.gz` — expression matrix
- `GPL3921.annot.gz` — probe-to-gene annotation
- `GSE14520_Extra_Supplement.txt.gz` — survival annotation

The U133A platform does **not** cover `CFHR3`, `ANGPTL6`, `FAM83F` or `PKN3`; these genes are
therefore reported as `NA` for this cohort.

---

## 3. ICGC-LIRI-JP — external validation cohort

| Item | Detail |
|---|---|
| Source | ICGC Data Portal — https://dcc.icgc.org/ |
| Project | LIRI-JP (liver cancer, Japan) |
| Release | release 28 |
| Expression | FKPM normalised gene expression (`exp_seq`) |
| Samples used | 222 tumours, 41 deaths |

Download the donor and `exp_seq` files for LIRI-JP and keep the directory layout unchanged;
`scripts/08_external_validation/parse_icgc_liri_expression.py` walks the release tree.

---

## 4. Structures used for docking

| Target | Structure | Source |
|---|---|---|
| PTH1R | PDB **6NBF** (cryo-EM; PTH1R–LA-PTH–Gs complex). Receptor chain R retained; G-protein, nanobody and lipid ligands removed | https://www.rcsb.org/structure/6NBF |
| SLCO4C1 | **AlphaFold** model `Q6ZQN7`, v6 (no experimental structure available). Mean pLDDT 78.7; 46.2% of residues ≥ 90 | https://alphafold.ebi.ac.uk/entry/Q6ZQN7 |

---

## 5. Small-molecule and interaction databases

| Database | Version / access | Used for |
|---|---|---|
| DGIdb | 5.0, GraphQL API — https://dgidb.org/ | drug-gene interactions |
| ChEMBL | EBI REST API — https://www.ebi.ac.uk/chembl/ | quantitative bioactivity (pChEMBL) |
| NPASS | 3.0 — http://bidd.group/NPASS/ | natural-product activity and species sources |
| PubChem | — https://pubchem.ncbi.nlm.nih.gov/ | compound identification (InChIKey), BioAssay potency |
| STRING | v12 — https://string-db.org/ | protein-protein interactions (combined score ≥ 0.4) |

Raw API responses are cached by the scripts under `work/` so that re-runs do not depend on
network availability.
