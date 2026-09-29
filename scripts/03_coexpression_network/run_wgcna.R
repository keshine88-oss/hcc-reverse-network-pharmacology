PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# WGCNA re-run (corrected): signed network + mitochondrial genes removed
# Input: dds_B_final.rds (fitted DESeq object; counts after mitochondrial removal + plate + condition)
# ============================================================

suppressMessages({
  library(DESeq2)
  library(WGCNA)
})
cor <- WGCNA::cor
options(stringsAsFactors = FALSE)

dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
cat("dds genes:", nrow(dds), "samples:", ncol(dds), "\n")

# VST transformation
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
datExpr <- t(mat)  # samples x genes
cat("VST matrix:", nrow(datExpr), "x", ncol(datExpr), "\n")

# Filter bad genes
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# top 5000 MAD
mad_vals <- apply(datExpr, 2, mad)
mad_rank <- order(mad_vals, decreasing = TRUE)
datExpr <- datExpr[, mad_rank[1:min(5000, length(mad_rank))]]
gene_names <- colnames(datExpr)
cat("WGCNA input:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# Grouping (use the dds colData; matched and aligned by barcode)
condition <- setNames(as.character(colData(dds)$condition), colnames(dds))
condition <- condition[rownames(datExpr)]
condition <- factor(condition, levels = c("Normal", "Tumor"))

# Soft-threshold selection (signed network)
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
power <- sft$powerEstimate
if (is.na(power) || power > 20) {
  r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  power <- sft$fitIndices[which(r2 > 0.85)[1], 1]
  if (is.na(power)) power <- 14
}
cat("signed network soft threshold power =", power, "(R2 =", round(sft$fitIndices[power, 2], 3), ")\n")

# Build the signed network
net <- blockwiseModules(datExpr, power = power,
                        networkType = "signed",
                        TOMType = "signed",
                        minModuleSize = 30,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        saveTOMs = FALSE,
                        verbose = 0)
moduleColors <- labels2colors(net$colors)
cat("number of modules:", length(unique(moduleColors)), "\n")
cat("module size distribution:\n")
print(table(moduleColors))

# Module-trait correlation
trait <- data.frame(Tumor = as.numeric(condition == "Tumor"),
                    Normal = as.numeric(condition == "Normal"))
rownames(trait) <- rownames(datExpr)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
moduleTraitCor <- cor(MEs0, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Significant modules (Tumour column only)
sig_idx <- which(moduleTraitPvalue[, 1] < 0.05 & abs(moduleTraitCor[, 1]) > 0.3)
sig_colors <- gsub("^ME", "", rownames(moduleTraitCor)[sig_idx])
cat("significant modules:", length(sig_colors), " (", paste(sig_colors, collapse=", "), ")\n")

# Collect module genes + hubs
module_genes <- unique(unlist(lapply(sig_colors, function(col) gene_names[moduleColors == col])))
cat("module genes:", length(module_genes), "\n")

GS <- as.numeric(cor(datExpr[, module_genes, drop=FALSE], trait$Tumor, use="p"))
names(GS) <- module_genes
hub_genes <- character(0)
for (col in sig_colors) {
  g <- gene_names[moduleColors == col]
  if (length(g) == 0) next
  me <- paste0("ME", col)
  if (!me %in% colnames(MEs0)) next
  MM <- as.numeric(cor(datExpr[, g, drop=FALSE], MEs0[, me], use="p"))
  names(MM) <- g
  sel <- names(MM)[abs(MM) > 0.8 & abs(GS[names(MM)]) > 0.2]
  if (length(sel) == 0) sel <- names(MM)[abs(MM) > 0.7 & abs(GS[names(MM)]) > 0.15]
  hub_genes <- c(hub_genes, sel)
}
hub_genes <- unique(hub_genes)
cat("hub genes:", length(hub_genes), "\n")

# Save
write.csv(data.frame(Gene = module_genes), file.path(PROJECT_ROOT, "work/main/wgcna_module_genes_final.csv"), row.names = FALSE)
write.csv(data.frame(Gene = hub_genes), file.path(PROJECT_ROOT, "work/main/wgcna_hub_genes_final.csv"), row.names = FALSE)
write.csv(data.frame(Gene = gene_names, ModuleColor = moduleColors),
          file.path(PROJECT_ROOT, "work/main/wgcna_module_colors_final.csv"), row.names = FALSE)
write.csv(data.frame(Module = rownames(moduleTraitCor),
                     Cor_Tumor = moduleTraitCor[, 1], P_Tumor = moduleTraitPvalue[, 1]),
          file.path(PROJECT_ROOT, "work/main/wgcna_module_trait_final.csv"), row.names = FALSE)
saveRDS(list(power = power, moduleColors = moduleColors, module_genes = module_genes,
             hub_genes = hub_genes, sig_colors = sig_colors, gene_names = gene_names),
        file.path(PROJECT_ROOT, "work/main/wgcna_final.rds"))

cat("\nWGCNA finished. power =", power, "| module genes", length(module_genes), "| hub", length(hub_genes), "\n")
