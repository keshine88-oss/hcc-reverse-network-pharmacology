PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Extract the full WGCNA plotting data (re-run the signed power=14 network and save all objects)
# ============================================================

suppressMessages({library(DESeq2); library(WGCNA)})
cor <- WGCNA::cor
options(stringsAsFactors = FALSE)

dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE)
datExpr <- t(assay(vsd))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# Sample clustering tree (QC)
sampleTree <- hclust(dist(datExpr), method = "average")

# top 5000 MAD
mad_vals <- apply(datExpr, 2, mad)
mad_rank <- order(mad_vals, decreasing = TRUE)
datExpr <- datExpr[, mad_rank[1:5000]]
gene_names <- colnames(datExpr)
condition <- setNames(as.character(colData(dds)$condition), colnames(dds))[rownames(datExpr)]
condition <- factor(condition, levels = c("Normal", "Tumor"))

# Soft threshold
sft <- pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "signed", verbose = 0)
r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
power <- sft$fitIndices[which(r2 > 0.88)[1], 1]
if (is.na(power)) power <- 14

# Network construction
net <- blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
                        saveTOMs = FALSE, verbose = 0)
moduleColors <- labels2colors(net$colors)

# Module eigengenes
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes

# Module-trait
trait <- data.frame(Tumor = as.numeric(condition == "Tumor"),
                    Normal = as.numeric(condition == "Normal"))
rownames(trait) <- rownames(datExpr)
moduleTraitCor <- cor(MEs0, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Significant modules
sig_idx <- which(moduleTraitPvalue[, 1] < 0.05 & abs(moduleTraitCor[, 1]) > 0.3)
sig_colors <- sub("^ME", "", rownames(moduleTraitCor)[sig_idx])

# Full GS/MM data (all modules)
gs_mm_list <- list()
for (col in unique(moduleColors)) {
  g <- gene_names[moduleColors == col]
  me <- paste0("ME", col)
  if (!me %in% colnames(MEs0)) next
  MM <- as.numeric(cor(datExpr[, g, drop = FALSE], MEs0[, me], use = "p"))
  GS <- as.numeric(cor(datExpr[, g, drop = FALSE], trait$Tumor, use = "p"))
  gs_mm_list[[col]] <- data.frame(Gene = g, Module = col, MM = MM, GS = GS, stringsAsFactors = FALSE)
}
gs_mm_all <- do.call(rbind, gs_mm_list)

# Module sizes
module_sizes <- as.data.frame(table(moduleColors), stringsAsFactors = FALSE)
names(module_sizes) <- c("Module", "Genes")

# Save all plotting data
saveRDS(list(
  sampleTree = sampleTree, sft = sft, power = power, net = net,
  moduleColors = moduleColors, gene_names = gene_names, MEs = MEs0,
  moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue,
  sig_colors = sig_colors, gs_mm_all = gs_mm_all, module_sizes = module_sizes,
  condition = condition
), file.path(PROJECT_ROOT, "work/main/wgcna_plot_data.rds"))

cat("Data extraction finished. power =", power, " modules =", length(unique(moduleColors)),
    " Significant modules =", length(sig_colors), "(", paste(sig_colors, collapse=", "), ")\n")
cat("module sizes:\n"); print(module_sizes)
