# Standard single-cohort transcriptomic workflow nominates candidate targets in hepatocellular carcinoma but fails external validation

Code repository accompanying the manuscript:

> **Standard single-cohort transcriptomic workflow nominates candidate targets in hepatocellular carcinoma but fails external validation: a reverse network-pharmacology pipeline with a methodological sensitivity analysis**

---

## What this repository contains

A complete, reproducible reverse network-pharmacology pipeline for hepatocellular carcinoma (HCC), built on the TCGA-LIHC cohort as the discovery set:

```
differential expression (DESeq2)
  -> weighted gene co-expression network analysis (WGCNA)
  -> candidate gene intersection
  -> functional enrichment (GO / KEGG)
  -> four-algorithm machine-learning consensus (random forest / SVM / LASSO / XGBoost)
  -> survival analysis (Kaplan-Meier + univariate Cox) with BH correction
  -> druggability screening
  -> reverse drug mining (DGIdb 5.0, ChEMBL)
  -> traditional Chinese medicine active-compound screening (NPASS 3.0)
  -> molecular docking (AutoDock Vina 1.1.2)
  -> external validation in GSE14520 and ICGC-LIRI-JP
  -> methodological sensitivity analysis (removal of the WGCNA grey module)
```

**Main finding (a negative/methodological result).** The prognostic signals nominated from the single TCGA-LIHC cohort do **not** replicate in either independent cohort: no gene is both significant and direction-consistent across cohorts, and SLCO4C1 and PZP show significant directional reversal. A sensitivity analysis further shows that the nomination is fragile to a *single* methodological choice — SLCO4C1 depends entirely on inclusion of the WGCNA grey (unassigned) module and drops out once it is removed. All nominated genes are therefore exploratory hypotheses, not established targets.

> The docking results are computational screening only and do not constitute measured binding or mechanistic evidence. The SLCO4C1 receptor is an AlphaFold predicted model, not an experimental structure.

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt              # Python dependencies (+ version pins)
├── R_packages.md                 # R package versions
├── .gitignore
├── data/
│   └── README.md                 # how to obtain the public datasets (no data redistributed)
├── docs/
│   └── pipeline.md               # stage-by-stage run order and inputs/outputs
├── scripts/
│   ├── 01_data_preparation/      # download and assemble TCGA-LIHC expression + clinical data
│   ├── 02_differential_expression/
│   ├── 03_coexpression_network/
│   ├── 04_machine_learning/
│   ├── 05_ppi_and_survival/
│   ├── 06_functional_enrichment/
│   ├── 07_sensitivity_analysis/  # grey module removal
│   ├── 08_external_validation/   # GSE14520 + ICGC-LIRI-JP
│   ├── 09_reverse_drug_mining/   # DGIdb, ChEMBL, NPASS
│   ├── 10_molecular_docking/     # Vina preparation, docking, H-bond analysis
│   └── 11_figure_generation/     # all main and supplementary figures
└── tools/                        # submission-figure utilities (compliance check, export)
```

---

## Data

**No data are redistributed in this repository.** All datasets are public and must be downloaded from their original sources; see [`data/README.md`](data/README.md).

| Dataset | Role | Source |
|---|---|---|
| TCGA-LIHC | discovery cohort | GDC Data Portal (STAR-Counts, 424 samples) |
| GSE14520 | external validation | GEO (Affymetrix U133A, GPL3921) |
| ICGC-LIRI-JP | external validation | ICGC Data Portal (release 28, FKPM) |
| PDB 6NBF | PTH1R receptor structure | RCSB PDB |
| AlphaFold Q6ZQN7 | SLCO4C1 predicted structure | AlphaFold DB (v6) |
| NPASS 3.0 | natural-product activity and species sources | bidd.group |
| DGIdb 5.0 | drug-gene interactions | dgidb.org (GraphQL API) |
| ChEMBL | bioactivity records | EBI REST API |

---

## Environment

**R 4.5.2** with DESeq2 1.50.2, apeglm 1.32.0, WGCNA 1.74, glmnet 5.1, xgboost 3.2.1.1, randomForest 4.7.1.2, e1071 1.7.17, clusterProfiler 4.18.4, org.Hs.eg.db 3.22.0, survival 3.8.9, survminer 0.5.2, igraph 2.3.3, pROC 1.19.1, ggplot2 4.0.3, patchwork 1.3.2.
Full list: [`R_packages.md`](R_packages.md).

**Python ≥ 3.11** with RDKit 2026.03.6, Meeko 0.8.0, gemmi 0.7.5, matplotlib 3.11.1, Pillow 12.3.0, numpy, scipy, pypdfium2.
Install: `pip install -r requirements.txt`.

**External binaries**: AutoDock Vina 1.1.2 and PyMOL 2 (used headless as `pymol -cq`) — install separately, they are not pip/R packages. Vina configuration files must be pure ASCII.

---

## Configuration

All scripts resolve paths from a single environment variable, so no absolute paths are hard-coded:

```bash
export HCC_PROJECT_ROOT=/path/to/this/repository     # Linux/macOS
set HCC_PROJECT_ROOT=D:\path\to\this\repository      # Windows
```

If it is not set, scripts fall back to the repository root (Python) or the working directory (R). Scripts write intermediate files and outputs into `work/`, `results/`, `figures/` and `docking/` under that root.

---

## Running the pipeline

See [`docs/pipeline.md`](docs/pipeline.md) for the full run order, inputs and outputs of every script. A short version:

```bash
# 1. data preparation
Rscript scripts/01_data_preparation/prepare_expression_matrix.R
python scripts/01_data_preparation/fetch_clinical_data.py

# 2. discovery analysis
Rscript scripts/02_differential_expression/run_deseq2_deg.R
Rscript scripts/03_coexpression_network/run_wgcna.R
Rscript scripts/04_machine_learning/run_ml_consensus.R
Rscript scripts/05_ppi_and_survival/survival_analysis.R

# 3. sensitivity analysis
Rscript scripts/04_machine_learning/run_ml_consensus_no_grey.R
Rscript scripts/07_sensitivity_analysis/run_survival_no_grey.R

# 4. external validation
python scripts/08_external_validation/build_gse14520_cohort.py
Rscript scripts/08_external_validation/cox_gse14520.R
Rscript scripts/08_external_validation/cox_icgc_liri.R

# 5. downstream pharmacology
python scripts/09_reverse_drug_mining/fetch_dgidb_chembl.py
python scripts/10_molecular_docking/run_vina_docking.py
```

---

## Reproducibility

- All stochastic steps (cross-validation splits, random forest and XGBoost training, LASSO cross-validation, docking) use a **fixed random seed of 42** in the main analysis.
- DESeq2, WGCNA, enrichment analysis and survival analysis are deterministic and therefore reproducible by construction.
- Docking was repeated with **4 independent seeds** (42, 1, 7, 2024) as a stability check.
- Software versions and all analysis parameters are listed in `R_packages.md` and `requirements.txt`.

---

## Notes and limitations

- This is a **single-cohort, retrospective computational study**; there are no functional experiments and no adjustment for clinical confounders.
- Prognostic associations are unadjusted univariate results and were corrected for multiple testing only within the discovery analysis.
- External validation used median-dichotomised univariate Cox models without an additional cross-cohort multiple-testing correction; it is a qualitative contrast of direction and significance, not a confirmatory test.
- Machine-learning input features were themselves selected by differential expression, so the high AUC values are an expected consequence of how the candidate pool was built and are not independent evidence of model strength.
- Molecular docking relies on a scoring function and should be regarded as preliminary screening.

---

## Citation

If you use this code, please cite the manuscript (details to be added upon publication) and this repository. See [`CITATION.cff`](CITATION.cff).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
