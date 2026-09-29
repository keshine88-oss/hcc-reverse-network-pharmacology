PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Full DEG figure set (style follows the WGCNA figures: theme_pub + pal)
# Produce: D1 volcano / D2 MA / D3 top 50 heatmap / D4 PCA / D5 sample-distance heatmap
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(patchwork); library(pheatmap)
})

out_dir <- file.path(PROJECT_ROOT, "work/main/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal <- list(blue = "#3F74B5", red = "#CF4E52", grey = "#737B88", dark = "#222222", light = "#A8ADB4")

theme_pub <- function(base_size = 7.3) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.3, colour = "#333333"),
      plot.title = element_text(size = 8.7, face = "bold", hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.0, colour = "#555555", hjust = 0.5, margin = margin(b = 2)),
      panel.grid = element_blank(),
      plot.margin = margin(7, 8, 7, 8)
    )
}
theme_set(theme_pub())

# ---- Data ----
deg <- read.csv(file.path(PROJECT_ROOT, "work/main/deg_B_final_all_results.csv"), row.names = 1)
dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE)
condition <- colData(dds)$condition
mat <- assay(vsd)

deg$gene <- rownames(deg)
deg$significance <- "Not significant"
deg$significance[deg$padj < 0.05 & deg$log2FoldChange >= 1] <- "Up-regulated"
deg$significance[deg$padj < 0.05 & deg$log2FoldChange <= -1] <- "Down-regulated"
deg$significance <- factor(deg$significance, levels = c("Up-regulated", "Down-regulated", "Not significant"))
deg$logP <- -log10(deg$padj)
n_up <- sum(deg$significance == "Up-regulated")
n_dn <- sum(deg$significance == "Down-regulated")

# ============================================================
# D1 Volcano plot
# ============================================================
p_d1 <- ggplot(deg, aes(log2FoldChange, logP, colour = significance)) +
  geom_point(size = 0.5, alpha = 0.55) +
  scale_colour_manual(values = c("Up-regulated" = pal$red, "Down-regulated" = pal$blue,
                                 "Not significant" = pal$light)) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3, colour = pal$grey) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = pal$grey) +
  labs(title = "Differential expression",
       subtitle = sprintf("%d up-regulated | %d down-regulated", n_up, n_dn),
       x = "Log2 fold change", y = "-Log10(adjusted P)") +
  theme(legend.position = "none")

# ============================================================
# D2 MA plot
# ============================================================
deg$meanExpr <- log10(deg$baseMean + 1)
p_d2 <- ggplot(deg, aes(meanExpr, log2FoldChange, colour = significance)) +
  geom_point(size = 0.5, alpha = 0.55) +
  scale_colour_manual(values = c("Up-regulated" = pal$red, "Down-regulated" = pal$blue,
                                 "Not significant" = pal$light)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = pal$grey) +
  labs(title = "MA plot", x = "Log10(mean expression)", y = "Log2 fold change") +
  theme(legend.position = "none")

# ============================================================
# D3 Top 50 DEG heatmap
# ============================================================
top50 <- head(deg[order(deg$padj), ], 50)
top50_mat <- mat[top50$gene, , drop = FALSE]
top50_mat <- t(scale(t(top50_mat)))
top50_mat[top50_mat > 2.5] <- 2.5
top50_mat[top50_mat < -2.5] <- -2.5
annotation_col <- data.frame(Condition = condition, row.names = colnames(mat))
ann_colors <- list(Condition = c(Tumor = pal$red, Normal = pal$blue))
pheatmap(top50_mat, annotation_col = annotation_col, annotation_colors = ann_colors,
         cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = FALSE,
         show_rownames = TRUE, fontsize_row = 5, fontsize_col = 5,
         color = colorRampPalette(c(pal$blue, "white", pal$red))(100),
         main = "Top 50 differentially expressed genes",
         filename = file.path(out_dir, "D3_top50_heatmap.png"), width = 7, height = 8, dpi = 300)

# ============================================================
# D4 PCA
# ============================================================
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
pve <- round(100 * summary(pca)$importance[2, 1:2], 1)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = condition)
p_d4 <- ggplot(pca_df, aes(PC1, PC2, colour = Condition)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = c(Tumor = pal$red, Normal = pal$blue)) +
  labs(title = "Principal component analysis",
       subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%%", pve[1], pve[2]),
       x = sprintf("PC1 (%s%%)", pve[1]), y = sprintf("PC2 (%s%%)", pve[2])) +
  theme(legend.position = "right", legend.title = element_blank())

# ============================================================
# D5 Sample-distance heatmap
# ============================================================
sampleDists <- dist(t(mat))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(c(pal$blue, "white", pal$red))(100),
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Sample-to-sample distance",
         filename = file.path(out_dir, "D5_sample_distance.png"), width = 6, height = 6, dpi = 300)

# ---- Export ggplot figures ----
suppressMessages(library(ragg))
agg_png(file.path(out_dir, "D1_volcano.png"), width = 4.2, height = 3.6, units = "in", res = 300, background = "white")
print(p_d1); dev.off()
agg_png(file.path(out_dir, "D2_MA.png"), width = 4.2, height = 3.6, units = "in", res = 300, background = "white")
print(p_d2); dev.off()
agg_png(file.path(out_dir, "D4_PCA.png"), width = 4.2, height = 3.4, units = "in", res = 300, background = "white")
print(p_d4); dev.off()

cat("DEG figures generated. Up/down:", n_up, "/", n_dn, "\n")
