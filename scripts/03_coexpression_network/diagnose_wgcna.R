PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
suppressMessages({library(DESeq2); library(WGCNA)})
cor <- WGCNA::cor
options(stringsAsFactors = FALSE)

dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE)
datExpr <- t(assay(vsd))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
mad_vals <- apply(datExpr, 2, mad)
mad_rank <- order(mad_vals, decreasing = TRUE)
datExpr <- datExpr[, mad_rank[1:5000]]
condition <- colData(dds)$condition[rownames(datExpr)]

trait <- data.frame(Tumor = as.numeric(condition == "Tumor"),
                    Normal = as.numeric(condition == "Normal"))
rownames(trait) <- rownames(datExpr)

for (pw in c(12, 14)) {
  net <- blockwiseModules(datExpr, power = pw, networkType = "signed", TOMType = "signed",
                          minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
                          saveTOMs = FALSE, verbose = 0)
  moduleColors <- labels2colors(net$colors)
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  mtc <- cor(MEs0, trait, use = "p")
  mtp <- corPvalueStudent(mtc, nrow(datExpr))
  cat("=== power =", pw, " modules =", length(unique(moduleColors)), "===\n")
  for (i in 1:nrow(mtc)) {
    cat(sprintf("  %-10s cor(Tumor)=%+.3f p=%.2e\n", rownames(mtc)[i], mtc[i,1], mtp[i,1]))
  }
  cat("  significant:", sum(mtp[,1] < 0.05 & abs(mtc[,1]) > 0.3), "\n\n")
}
