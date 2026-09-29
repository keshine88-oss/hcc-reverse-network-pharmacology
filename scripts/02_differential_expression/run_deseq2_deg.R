PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# DEG re-run (corrected): scheme B batch correction + glmGamPoi for speed
# Input: count_clean.rds (35,281 genes, mitochondrial + ribosomal genes already removed)
# design = ~ plate + condition
# ============================================================

suppressMessages(library(DESeq2))

cm <- readRDS(file.path(PROJECT_ROOT, "work/main/count_clean.rds"))
condition <- readRDS(file.path(PROJECT_ROOT, "work/main/condition_clean.rds"))
plate <- readRDS(file.path(PROJECT_ROOT, "work/main/plate_clean.rds"))

cat("input:", nrow(cm), "genes x", ncol(cm), "samples\n")
cat("groups: Tumour", sum(condition == "Tumor"), " Normal", sum(condition == "Normal"), "\n")

coldata <- data.frame(condition = condition, plate = plate, row.names = colnames(cm))
dds <- DESeqDataSetFromMatrix(countData = round(cm), colData = coldata,
                              design = ~ plate + condition)

# Pre-filter
dds <- dds[rowSums(counts(dds) >= 10) >= 50, ]
cat("genes after pre-filtering:", nrow(dds), "\n")

# IRLS fitting (DESeq2 default parametric, Wald test)
t0 <- Sys.time()
dds <- DESeq(dds, quiet = TRUE)
t1 <- Sys.time()
cat("DESeq runtime:", round(as.numeric(difftime(t1, t0, units = "secs")), 1), "s\n")

# Results + apeglm shrinkage
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- lfcShrink(dds, coef = "condition_Tumor_vs_Normal", res = res, type = "apeglm")

# Significance
res_df <- as.data.frame(res)
deg_idx <- which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1 & !is.na(res_df$padj))
up <- which(res_df$padj < 0.05 & res_df$log2FoldChange >= 1 & !is.na(res_df$padj))
dn <- which(res_df$padj < 0.05 & res_df$log2FoldChange <= -1 & !is.na(res_df$padj))
cat("Total DEGs:", length(deg_idx), " (up", length(up), "/ down", length(dn), ")\n")

# Save
write.csv(res_df, file.path(PROJECT_ROOT, "work/main/deg_B_final_all_results.csv"))
writeLines(rownames(res_df)[deg_idx], file.path(PROJECT_ROOT, "work/main/deg_B_final.txt"))
writeLines(rownames(res_df)[up], file.path(PROJECT_ROOT, "work/main/deg_B_final_up.txt"))
writeLines(rownames(res_df)[dn], file.path(PROJECT_ROOT, "work/main/deg_B_final_dn.txt"))
saveRDS(dds, file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
saveRDS(res_df, file.path(PROJECT_ROOT, "work/main/res_B_final.rds"))

cat("Saved deg_B_final.txt / deg_B_final_all_results.csv / dds_B_final.rds\n")
cat("Done.\n")
