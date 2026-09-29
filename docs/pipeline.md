# Pipeline: run order, inputs and outputs

All paths below are relative to `$HCC_PROJECT_ROOT` (see the main README).
Run the stages in order; each stage consumes the outputs of the previous one.

```
work/       intermediate objects (RDS, cached API responses)
results/    analysis result tables (CSV)
figures/    generated figures
docking/    receptors, ligands, docking poses and scores
```

---

## Stage 1 — Data preparation

| Script | Input | Output |
|---|---|---|
| `01_data_preparation/download_tcga_expression.py` | GDC manifest + downloaded STAR-Counts files | file inventory for the count matrix |
| `01_data_preparation/prepare_expression_matrix.R` | raw counts | `work/count_clean.rds`, `condition_clean.rds`, `plate_clean.rds` |
| `01_data_preparation/fetch_clinical_data.py` | GDC cases API | `work/tcga_clinical.csv` |
| `01_data_preparation/trace_clinical_fields.R` | GDC clinical JSON | field-mapping diagnostics |
| `01_data_preparation/precheck_inputs.R` | count / metadata files | sanity checks before analysis |

Key preprocessing decisions: 13 mitochondrial and 1,483 ribosomal/pseudogene genes are removed
(1,496 in total, leaving 35,281 genes × 424 samples); batch metadata (plate, TSS, centre) is
extracted from the TCGA barcodes.

---

## Stage 2 — Differential expression

| Script | Input | Output |
|---|---|---|
| `02_differential_expression/run_deseq2_deg.R` | `work/count_clean.rds` | `results/deg_final.txt`, `results/deg_final_all_results.csv`, `work/dds_final.rds` |

Design formula: `~ plate + condition` (plate as a 20-level batch covariate).
Pre-filter: count ≥ 10 in ≥ 50 samples. Threshold: |log2FC| ≥ 1 and BH-adjusted P < 0.05.
log2FC shrinkage: apeglm. Result: 4,273 DEGs (2,994 up, 1,279 down).
`04_machine_learning/reproducibility_double_run.R` runs the same analysis with and without
batch correction and reports how many DEGs are gained or lost.

---

## Stage 3 — Co-expression network (WGCNA)

| Script | Input | Output |
|---|---|---|
| `03_coexpression_network/run_wgcna.R` | VST matrix | `work/wgcna_objects.rds` |
| `03_coexpression_network/diagnose_wgcna.R` | VST matrix | soft-threshold diagnostics |

Signed network on the 5,000 most variable genes (MAD); soft threshold power = 14
(scale-free R² = 0.888); minModuleSize = 30; mergeCutHeight = 0.25.
Nine modules are detected; five are tumour-related (|r| > 0.3, P < 0.05), giving 3,622 module
genes. **Note that the grey (unassigned) module is admitted by this automatic threshold** —
this is the methodological choice examined in Stage 6.

---

## Stage 4 — Candidate genes and machine-learning consensus

| Script | Input | Output |
|---|---|---|
| `04_machine_learning/run_ml_consensus.R` | DEGs ∩ module genes | `results/candidate_pool.csv`, `results/consensus_genes.csv`, `results/consensus_ranking.csv`, `results/cv_auc.csv` |
| `04_machine_learning/run_ml_consensus_no_grey.R` | candidate pool without the grey module | as above, suffixed `_no_grey` |

Candidate pool: 4,273 DEGs ∩ 3,622 module genes = 1,817 genes.
Four algorithms (random forest, linear SVM, LASSO logistic regression, XGBoost) rank the
candidate genes; a gene is a consensus hit if it appears in the top 20 of **at least two**
algorithms. This yields 16 consensus genes. Performance is measured by stratified 5-fold
out-of-fold AUC with seed 42.

---

## Stage 5 — PPI network, survival analysis and multiple-testing correction

| Script | Input | Output |
|---|---|---|
| `05_ppi_and_survival/fetch_string_ppi.py` | STRING v12 | `results/ppi_edges.tsv` |
| `05_ppi_and_survival/ppi_network_analysis.R` | PPI edges + consensus genes | hub genes, degree table |
| `05_ppi_and_survival/survival_analysis.R` | VST matrix + clinical data | `results/survival_analysis.csv` |
| `05_ppi_and_survival/survival_multiple_testing.R` | survival results | BH-adjusted P values, Schoenfeld residuals |

Survival: median-dichotomised Kaplan-Meier (log-rank) plus univariate Cox models.
Seven of the fifteen protein-coding consensus genes are significant before correction; only
five survive BH correction (all protective).

---

## Stage 6 — Sensitivity analysis (grey module removal)

| Script | Input | Output |
|---|---|---|
| `04_machine_learning/run_ml_consensus_no_grey.R` | module genes excluding grey | consensus genes without grey |
| `07_sensitivity_analysis/run_survival_no_grey.R` | expression + clinical data | `results/survival_no_grey.csv` |
| `06_functional_enrichment/run_enrichment_no_grey.R` | grey-excluded candidate pool | GO/KEGG result tables |

The grey module consists of genes not assigned to any co-expression module and has no
co-expression cohesion, yet it passes the automatic |r| > 0.3 threshold. Removing it shrinks
the module gene set from 3,622 to 2,250 and the candidate pool from 1,817 to 1,361, changes six
consensus genes, and removes **SLCO4C1** from the candidate pool entirely — the central
methodological result of the study.

---

## Stage 7 — Functional enrichment

| Script | Input | Output |
|---|---|---|
| `06_functional_enrichment/run_enrichment.R` | candidate genes | GO (BP/CC/MF) and KEGG tables |

BH-adjusted P < 0.05. Main analysis: 751 BP, 65 CC, 116 MF terms and 50 KEGG pathways.
Grey-excluded analysis: 671 BP, 62 CC, 111 MF terms and 46 KEGG pathways (Stage 6).

---

## Stage 8 — External validation

| Script | Input | Output |
|---|---|---|
| `08_external_validation/build_gse14520_cohort.py` | GSE14520 matrix + GPL3921 annotation + survival supplement | validation cohort table (221 tumours, 85 deaths) |
| `08_external_validation/download_icgc_liri.py` | ICGC release 28 | LIRI-JP donor and expression files |
| `08_external_validation/parse_icgc_liri_expression.py` | LIRI-JP files | cohort table (222 tumours, 41 deaths) |
| `08_external_validation/cox_gse14520.R` | GSE14520 cohort | per-gene HR and P |
| `08_external_validation/cox_icgc_liri.R` | LIRI-JP cohort | per-gene HR and P |
| `08_external_validation/build_cohort_comparison_table.R` | both cohorts | three-cohort comparison table |

The same median-dichotomised univariate Cox model is applied in each cohort. No gene is both
significant and direction-consistent across cohorts; SLCO4C1 and PZP reverse significantly.
No additional cross-cohort multiple-testing correction is applied — the comparison is a
qualitative contrast of direction and significance.

---

## Stage 9 — Reverse drug mining and TCM screening

| Script | Input | Output |
|---|---|---|
| `09_reverse_drug_mining/fetch_dgidb_chembl.py` | DGIdb GraphQL + ChEMBL REST | `work/drug/dgidb_interactions.csv`, `chembl_ligands.csv` |
| `09_reverse_drug_mining/parse_drug_interactions.py` | cached JSON | drug summary tables |
| `09_reverse_drug_mining/parse_npass_targets.py` | NPASS target pages | target-to-compound lists |
| `09_reverse_drug_mining/parse_npass_compound_pages.py` | NPASS compound pages | activities, species sources, SMILES |
| `09_reverse_drug_mining/build_tcm_hit_table.py` | NPASS + curated species-to-TCM map | `results/npass_tcm_hits.csv`, herb summary |

PTH1R: 47 DGIdb interactions, 142 ChEMBL ligands. SLCO4C1: 3 DGIdb interactions, 3 ChEMBL
ligands. NPASS yields 245 unique natural-product hits (PTH1R 243, SLCO4C1 3, with one shared
compound). Digitoxin was used as a single-page validation: NPASS IC50 120 nM matches the
independent ChEMBL record (pChEMBL 6.92).

---

## Stage 10 — Molecular docking

| Script | Input | Output |
|---|---|---|
| `10_molecular_docking/prepare_receptors.py` | PDB 6NBF, AlphaFold Q6ZQN7 | receptor PDBQT files |
| `10_molecular_docking/detect_binding_cavity.py` | AlphaFold model | grid-box centre for the SLCO4C1 central cavity |
| `10_molecular_docking/prepare_ligands.py` | SMILES | ligand PDBQT files (RDKit + Meeko) |
| `10_molecular_docking/run_vina_docking.py` | receptors + ligands | poses and binding energies |
| `10_molecular_docking/run_docking_multi_seed.py` | as above | stability check over 4 seeds |
| `10_molecular_docking/convert_pdbqt_to_pdb.py` | PDBQT poses | PDB poses for rendering |
| `10_molecular_docking/analyze_hydrogen_bonds.py` | poses | hydrogen bonds under strict angle criteria |

Boxes: PTH1R 22×22×22 Å centred on the LA-PTH N-terminus centroid (orthosteric pocket);
SLCO4C1 30×30×30 Å centred on the high-confidence transmembrane helix core.
Vina parameters: exhaustiveness = 16, num_modes = 9. H-bond criteria: D···A ≤ 3.5 Å,
H···A ≤ 2.6 Å, D–H···A ≥ 120°. Nine of ten ligands gave identical binding energies across the
four seeds; only cucurbitacin B varied, by 0.1 kcal/mol.

---

## Stage 11 — Figures

| Script | Figure |
|---|---|
| `11_figure_generation/figure1_differential_expression.R` | Figure 1 — DEG panel (volcano / MA / PCA / heatmap) |
| `11_figure_generation/figure2_wgcna.R`, `figure2_wgcna_enhanced.R` | Figure 2 — WGCNA panels |
| `11_figure_generation/figure3_candidate_venn.R` | Figure 3a — DEG ∩ WGCNA Venn diagram |
| `11_figure_generation/figure4_machine_learning.R` | Figure 4 — ML screening panels |
| `11_figure_generation/figure5_survival_composite.R`, `km_curve_panels.R` | Figure 5 — KM curves |
| `11_figure_generation/figure5_forest_plot.R` | Figure 5h — three-cohort forest plot |
| `11_figure_generation/figure6_binding_energy_barplot.R` | Figure 6a — binding-energy bar chart |
| `11_figure_generation/docking_figures/*` | Figure 6b–k — binding-site close-ups (PyMOL + Pillow) |
| `11_figure_generation/figure6_docking_panel_assembly.py` | Figure 6 — panel assembly |
| `11_figure_generation/tcm_screening_overview.py` | TCM screening overview (supplementary) |
| `11_figure_generation/extract_plot_data.R`, `assemble_panels.R`, `assemble_ml_panels.R` | plot-data export and composite assembly |

Docking figures use a two-stage approach: PyMOL renders the views headlessly and Pillow adds
pixel-level annotations (residue labels, hydrogen-bond dashes and distances).

---

## Tools

| Script | Purpose |
|---|---|
| `tools/make_submission_figures.py` | normalise figures to 2161 px wide (183 mm at 300 dpi), RGB, PNG + LZW TIFF |
| `tools/check_figure_compliance.py` | journal-compliance check (format, effective DPI, colour mode, size, naming) |
| `tools/estimate_figure_text_size.py` | estimate the printed point size of the smallest body text |
