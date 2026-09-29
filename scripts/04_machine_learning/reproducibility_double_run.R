PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Double run: batch-corrected vs uncorrected
# Scheme A: design = ~ condition (no batch correction)
# Scheme B: design = ~ plate + condition (batch corrected)
# NOTE: for count data (raw counts) correcting via DESeq2 design covariates is the correct approach;
#       ComBat/removeBatchEffect apply only to log/VST continuous data, never to raw counts.
# ============================================================

suppressMessages(library(DESeq2))

cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")

# ---- Data preparation ----
cm <- readRDS(file.path(PROJECT_ROOT, "work/main/count.rds"))
groups <- readRDS(file.path(PROJECT_ROOT, "work/main/groups.rds"))
batch <- read.csv(file.path(PROJECT_ROOT, "work/main/batch_info.csv"), stringsAsFactors = FALSE)

cat("raw count matrix:", nrow(cm), "genes x", ncol(cm), "samples\n")

# Remove mitochondrial/ribosomal genes (plan B)
keep <- !grepl("^MT-", rownames(cm)) & !grepl("^RPL|^RPS", rownames(cm))
cat("mitochondrial/ribosomal genes removed:", sum(!keep), "\n")
cm <- cm[keep, , drop = FALSE]
cat("after removal:", nrow(cm), "genes x", ncol(cm), "samples\n")

# Sample alignment
if (!is.null(names(groups))) groups <- groups[colnames(cm)]
plate <- batch$Plate[match(colnames(cm), batch$Barcode)]
condition <- factor(groups, levels = c("Normal", "Tumor"))

cat("samples: Tumour", sum(condition == "Tumor"), " Normal", sum(condition == "Normal"), "\n")
cat("number of plate levels:", length(unique(plate)), "\n")

# ---- Scheme A: no batch correction ----
cat("\n========== Scheme A: design = ~ condition ==========\n")
gc(reset = TRUE)
tA <- system.time({
  ddsA <- DESeqDataSetFromMatrix(countData = round(cm),
                                 colData = data.frame(condition = condition),
                                 design = ~ condition)
  ddsA <- ddsA[rowSums(counts(ddsA) >= 10) >= 50, ]
  ddsA <- DESeq(ddsA, quiet = TRUE)
  resA <- results(ddsA, contrast = c("condition", "Tumor", "Normal"))
})
memA <- sum(gc()[, 6])  # max used (Mb) across Ncells+Vcells
degA <- sum(resA$padj < 0.05 & abs(resA$log2FoldChange) >= 1, na.rm = TRUE)
upA <- sum(resA$padj < 0.05 & resA$log2FoldChange >= 1, na.rm = TRUE)
dnA <- sum(resA$padj < 0.05 & resA$log2FoldChange <= -1, na.rm = TRUE)
cat("runtime:", round(tA["elapsed"], 2), "s\n")
cat("DEG count:", degA, " (up", upA, "/ down", dnA, ")\n")

# ---- Scheme B: batch correction (plate as covariate) ----
cat("\n========== Scheme B: design = ~ plate + condition ==========\n")
plate_f <- factor(plate)
gc(reset = TRUE)
tB <- system.time({
  ddsB <- DESeqDataSetFromMatrix(countData = round(cm),
                                 colData = data.frame(condition = condition, plate = plate_f),
                                 design = ~ plate + condition)
  ddsB <- ddsB[rowSums(counts(ddsB) >= 10) >= 50, ]
  ddsB <- DESeq(ddsB, quiet = TRUE)
  resB <- results(ddsB, contrast = c("condition", "Tumor", "Normal"))
})
memB <- sum(gc()[, 6])
degB <- sum(resB$padj < 0.05 & abs(resB$log2FoldChange) >= 1, na.rm = TRUE)
upB <- sum(resB$padj < 0.05 & resB$log2FoldChange >= 1, na.rm = TRUE)
dnB <- sum(resB$padj < 0.05 & resB$log2FoldChange <= -1, na.rm = TRUE)
cat("runtime:", round(tB["elapsed"], 2), "s\n")
cat("DEG count:", degB, " (up", upB, "/ down", dnB, ")\n")

# ---- Comparison ----
genesA <- rownames(resA)[resA$padj < 0.05 & abs(resA$log2FoldChange) >= 1 & !is.na(resA$padj)]
genesB <- rownames(resB)[resB$padj < 0.05 & abs(resB$log2FoldChange) >= 1 & !is.na(resB$padj)]
overlap <- length(intersect(genesA, genesB))
onlyA <- length(setdiff(genesA, genesB))
onlyB <- length(setdiff(genesB, genesA))

cat("\n========== Comparison results ==========\n")
cat("A uncorrected DEG:", length(genesA), "\n")
cat("B batch-corrected DEG:", length(genesB), "\n")
cat("overlap:", overlap, " overlap rate (relative to A):", round(overlap / length(genesA) * 100, 1), "%\n")
cat("only in A:", onlyA, " only in B:", onlyB, "\n")

# ---- Save results ----
write.csv(data.frame(
  config = c("A_no_batch", "B_batch_corrected"),
  design = c("~condition", "~plate+condition"),
  time_sec = round(c(tA["elapsed"], tB["elapsed"]), 2),
  max_mem_Mb = c(memA, memB),
  n_deg = c(degA, degB),
  n_up = c(upA, upB),
  n_down = c(dnA, dnB),
  stringsAsFactors = FALSE
), file.path(PROJECT_ROOT, "work/main/double_run_summary.csv"), row.names = FALSE)

# Save the DEG list (for later comparison)
writeLines(genesA, file.path(PROJECT_ROOT, "work/main/deg_A_no_batch.txt"))
writeLines(genesB, file.path(PROJECT_ROOT, "work/main/deg_B_batch.txt"))
# Present only in B (gained or lost after batch correction)
writeLines(setdiff(genesB, genesA), file.path(PROJECT_ROOT, "work/main/deg_only_B.txt"))
writeLines(setdiff(genesA, genesB), file.path(PROJECT_ROOT, "work/main/deg_only_A.txt"))

saveRDS(list(A = genesA, B = genesB, overlap = overlap), file.path(PROJECT_ROOT, "work/main/double_run_results.rds"))

cat("\n===== Double-run comparison finished =====\n")
