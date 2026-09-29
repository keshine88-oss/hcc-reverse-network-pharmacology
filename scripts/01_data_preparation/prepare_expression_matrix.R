PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Unified data preparation (corrected): handle mitochondrial gene removal and sample alignment properly
# Fix 1: remove mitochondrial genes using an exact list (MT. dot format, 13 true mtDNA genes)
# Fix 2: align sample groups explicitly via the Barcode column of batch_info.csv
# ============================================================

cm <- readRDS(file.path(PROJECT_ROOT, "work/main/count.rds"))
bt <- read.csv(file.path(PROJECT_ROOT, "work/main/batch_info.csv"), stringsAsFactors = FALSE)

cat("raw counts:", nrow(cm), "x", ncol(cm), "\n")

# ---- Fix 1: precisely remove mitochondrial genes + ribosomal genes ----
mito_genes <- c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB",
                "MT.ND1", "MT.ND2", "MT.ND3", "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6")
is_mito <- rownames(cm) %in% mito_genes
is_ribo <- grepl("^RPL|^RPS", rownames(cm))
keep <- !(is_mito | is_ribo)
cat("mitochondrial genes removed:", sum(is_mito), " (", paste(rownames(cm)[is_mito], collapse=", "), ")\n")
cat("ribosomal genes removed:", sum(is_ribo), "\n")
cm <- cm[keep, , drop = FALSE]
cat("after removal:", nrow(cm), "x", ncol(cm), "\n")

# ---- Fix 2: align sample groups explicitly via Barcode ----
# batch_info.csv contains Barcode / Group / SampleType / Plate / TSS / Center
rownames(bt) <- bt$Barcode
bt <- bt[colnames(cm), , drop = FALSE]
condition <- factor(bt$Group, levels = c("Normal", "Tumor"))
plate <- factor(bt$Plate)

cat("group alignment: Tumour", sum(condition == "Tumor"), " Normal", sum(condition == "Normal"), "\n")
cat("plate levels:", length(levels(plate)), "\n")

# ---- Save the corrected data (reused by downstream re-runs) ----
saveRDS(cm, file.path(PROJECT_ROOT, "work/main/count_clean.rds"))
saveRDS(condition, file.path(PROJECT_ROOT, "work/main/condition_clean.rds"))
saveRDS(plate, file.path(PROJECT_ROOT, "work/main/plate_clean.rds"))
write.csv(bt, file.path(PROJECT_ROOT, "work/main/meta_clean.csv"), row.names = FALSE)

cat("\nCorrection finished; count_clean.rds / condition_clean.rds / plate_clean.rds saved\n")
cat("final input:", nrow(cm), "genes x", ncol(cm), "samples\n")
