# R package versions

Environment used for the analyses reported in the manuscript.

- **R**: 4.5.2

| Package | Version | Used for |
|---|---|---|
| TCGAbiolinks | 2.38.0 | TCGA-LIHC download |
| DESeq2 | 1.50.2 | differential expression analysis |
| apeglm | 1.32.0 | log2 fold-change shrinkage |
| glmGamPoi | — | fast dispersion estimation (optional) |
| WGCNA | 1.74 | weighted gene co-expression network analysis |
| randomForest | 4.7.1.2 | random forest feature ranking |
| e1071 | 1.7.17 | support vector machine |
| glmnet | 5.1 | LASSO logistic regression |
| xgboost | 3.2.1.1 | gradient boosting |
| pROC | 1.19.1 | ROC curves and AUC |
| clusterProfiler | 4.18.4 | GO and KEGG enrichment |
| org.Hs.eg.db | 3.22.0 | human gene annotation |
| survival | 3.8.9 | Kaplan-Meier and Cox models |
| survminer | 0.5.2 | Kaplan-Meier plots |
| igraph | 2.3.3 | PPI network analysis |
| ggplot2 | 4.0.3 | plotting |
| patchwork | 1.3.2 | multi-panel figure assembly |
| ragg | — | high-resolution raster device (PNG/TIFF) |
| ggrepel | — | non-overlapping labels |

Install the core set with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("WGCNA", "randomForest", "e1071", "glmnet", "xgboost", "pROC",
                   "survival", "survminer", "igraph", "ggplot2", "patchwork",
                   "ragg", "ggrepel"))
```

## External tools (not installable from CRAN / PyPI)

| Tool | Version | Used for |
|---|---|---|
| AutoDock Vina | 1.1.2 | molecular docking scoring |
| PyMOL | 2 (Schrödinger) | headless rendering of binding modes (`pymol -cq`) |
