PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Full WGCNA figure set (style follows 11_wgcna_figure_enhancement.R)
# Produce: S1 sample clustering / F1 soft threshold / F2 gene tree / F3 module sizes / F4 module-trait heatmap / F5 ME clustering / F6 GS-MM
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(WGCNA)
})

dat <- readRDS(file.path(PROJECT_ROOT, "work/main/wgcna_plot_data.rds"))
sft <- dat$sft; power <- dat$power; net <- dat$net
moduleColors <- dat$moduleColors; gene_names <- dat$gene_names
MEs <- dat$MEs; moduleTraitCor <- dat$moduleTraitCor
moduleTraitPvalue <- dat$moduleTraitPvalue; sig_colors <- dat$sig_colors
gs_mm_all <- dat$gs_mm_all; module_sizes <- dat$module_sizes
sampleTree <- dat$sampleTree

out_dir <- file.path(PROJECT_ROOT, "work/main/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Colour palette (exact reproduction) ----
pal <- list(blue = "#3F74B5", red = "#CF4E52", grey = "#737B88",
            dark = "#222222", light = "#E8EDF3")
module_palette <- c(
  black = "#242424", blue = "#3F74B5", brown = "#8C5A3C", green = "#4C9A62",
  greenyellow = "#A6C94A", grey = "#A8ADB4", magenta = "#B45C9A", pink = "#E58AA8",
  purple = "#7E60A8", red = "#CF4E52", turquoise = "#31A6A1", yellow = "#E6B84A"
)

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

mods <- unique(moduleColors)
cat("modules:", paste(sort(mods), collapse=", "), "\n")

# ============================================================
# F1 soft-threshold diagnostic plot (a: scale-free fit, b: mean connectivity)
# ============================================================
fit <- as.data.frame(sft$fitIndices)
fit$power <- fit[, 1]
fit$signed_R2 <- -sign(fit[, 3]) * fit[, 2]
fit$mean_connectivity <- fit[, 5]
fit$selected <- fit$power == power

p_f1a <- ggplot(fit, aes(power, signed_R2)) +
  geom_hline(yintercept = 0.80, linetype = 2, linewidth = 0.35, colour = pal$grey) +
  geom_line(linewidth = 0.45, colour = pal$blue) +
  geom_point(aes(fill = selected), shape = 21, size = 2.0, stroke = 0.35, colour = "white") +
  geom_text(data = fit[fit$power %in% c(1:6, 10, 15, 20), ], aes(label = power), nudge_y = 0.035, size = 2.05, family = "Arial", colour = pal$dark) +
  scale_fill_manual(values = c("FALSE" = pal$blue, "TRUE" = pal$red), guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  coord_cartesian(ylim = c(min(0, min(fit$signed_R2) - 0.05), 1.02), clip = "off") +
  labs(title = "Scale-free topology fit", subtitle = sprintf("Selected power = %d; signed R\u00b2 = %.3f", power, fit$signed_R2[fit$selected]), x = "Soft-thresholding power", y = "Signed R\u00b2")

p_f1b <- ggplot(fit, aes(power, mean_connectivity)) +
  geom_line(linewidth = 0.45, colour = pal$blue) +
  geom_point(aes(fill = selected), shape = 21, size = 2.0, stroke = 0.35, colour = "white") +
  geom_text(data = fit[fit$power %in% c(1:6, 10, 15, 20), ], aes(label = power), nudge_y = max(fit$mean_connectivity) * 0.035, size = 2.05, family = "Arial", colour = pal$dark) +
  scale_fill_manual(values = c("FALSE" = pal$blue, "TRUE" = pal$red), guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  expand_limits(y = max(fit$mean_connectivity) * 1.11) +
  labs(title = "Mean connectivity", subtitle = sprintf("Mean connectivity at power %d = %.1f", power, fit$mean_connectivity[fit$selected]), x = "Soft-thresholding power", y = "Mean connectivity")

# ============================================================
# F3 module size bar plot
# ============================================================
module_sizes <- module_sizes[order(module_sizes$Genes), ]
module_sizes$Module <- factor(module_sizes$Module, levels = module_sizes$Module)
p_f3 <- ggplot(module_sizes, aes(Genes, Module, fill = Module)) +
  geom_col(width = 0.70, colour = "#555555", linewidth = 0.18) +
  geom_text(aes(label = Genes), hjust = -0.12, family = "Arial", size = 2.05, colour = pal$dark) +
  scale_fill_manual(values = module_palette, guide = "none") +
  expand_limits(x = max(module_sizes$Genes) * 1.16) +
  labs(title = "Module sizes", subtitle = "All 5,000 input genes assigned", x = "Genes", y = NULL)

# ============================================================
# F4 module-trait heatmap (positive/negative blocks, padded with blank rows for equal column height, wide cells)
# ============================================================
module_names <- sub("^ME", "", rownames(moduleTraitCor))
module_order <- module_names[order(moduleTraitCor[, 1], decreasing = TRUE)]

positive_mods <- module_names[moduleTraitCor[, 1] > 0]
negative_mods <- module_names[moduleTraitCor[, 1] <= 0]

heat <- rbind(
  data.frame(Module = module_names, Trait = "Tumor", Correlation = moduleTraitCor[, 1], Pvalue = moduleTraitPvalue[, 1]),
  data.frame(Module = module_names, Trait = "Normal", Correlation = -moduleTraitCor[, 1], Pvalue = moduleTraitPvalue[, 1])
)

# Pad with blank rows so that the Higher / Lower columns have the same number of rows
n_pos <- length(positive_mods); n_neg <- length(negative_mods)
if (n_pos < n_neg) {
  pad <- n_neg - n_pos
  pad_rows <- do.call(rbind, lapply(seq_len(pad), function(i) {
    data.frame(Module = paste0(" ", i), Trait = c("Tumor", "Normal"),
               Correlation = NA_real_, Pvalue = NA_real_)
  }))
  heat <- rbind(heat, pad_rows)
  positive_mods <- c(positive_mods, paste0(" ", seq_len(pad)))
}

heat$Block <- ifelse(heat$Module %in% positive_mods, "Higher tumor correlation", "Lower tumor correlation")
heat$Block <- factor(heat$Block, levels = c("Higher tumor correlation", "Lower tumor correlation"))
heat$Module <- factor(heat$Module, levels = c(rev(positive_mods), rev(negative_mods)))
heat$Trait <- factor(heat$Trait, levels = c("Tumor", "Normal"))

format_p <- function(x) ifelse(x < 0.001, format(x, scientific = TRUE, digits = 1), sprintf("%.3f", x))
heat$label <- ifelse(is.na(heat$Correlation), "", paste0(sprintf("%.2f", heat$Correlation), "\nP=", format_p(heat$Pvalue)))
heat$text_colour <- ifelse(abs(heat$Correlation) >= 0.38, "white", pal$dark)
heat$text_colour[is.na(heat$Correlation)] <- "white"

p_f4 <- ggplot(heat, aes(Trait, Module, fill = Correlation)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = label, colour = text_colour), family = "Arial", size = 2.15, lineheight = 0.90) +
  scale_fill_gradient2(low = pal$blue, mid = "white", high = pal$red, midpoint = 0,
                       limits = c(-0.7, 0.7), na.value = "white", name = "Pearson r") +
  scale_colour_identity() +
  scale_x_discrete(position = "top") +
  facet_wrap(~Block, nrow = 1, scales = "free_y") +
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust = 0.5,
                                barwidth = unit(48, "mm"), barheight = unit(3.2, "mm"), ticks.colour = "white")) +
  labs(title = "Complete module\u2013trait matrix", subtitle = sprintf("All %d modules; each cell reports Pearson r and nominal P value", length(module_names)), x = NULL, y = NULL) +
  theme_minimal(base_size = 7.2, base_family = "Arial") +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(face = "bold", colour = pal$dark),
    axis.text.y = element_text(colour = pal$dark, size = 7.0), legend.position = "bottom",
    strip.background = element_blank(), strip.text = element_text(size = 6.8, face = "bold", colour = "#555555"),
    panel.spacing.x = unit(7, "mm"), legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.2), plot.title = element_text(size = 8.7, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 7.0, colour = "#555555", hjust = 0.5),
    plot.margin = margin(7, 5, 7, 5)
  )

# ============================================================
# F5 module eigengene clustering tree
# ============================================================
module_eigengenes <- MEs
eigen_names <- sub("^ME", "", names(module_eigengenes))
if (!all(eigen_names %in% c("black","blue","brown","green","greenyellow","grey","magenta","pink","purple","red","turquoise","yellow"))) {
  module_id_map <- c(`0`="grey",`1`="turquoise",`2`="blue",`3`="brown",`4`="yellow",`5`="green",
                     `6`="red",`7`="black",`8`="pink",`9`="magenta",`10`="purple",`11`="greenyellow")
  names(module_eigengenes) <- paste0("ME", unname(module_id_map[eigen_names]))
}
module_tree <- hclust(as.dist(1 - cor(module_eigengenes, use = "pairwise.complete.obs")), method = "average")

dendrogram_segments <- function(hclust_object) {
  dendrogram_object <- as.dendrogram(hclust_object)
  leaf_index <- 0L; segment_index <- 0L; segment_list <- list()
  walk <- function(node) {
    node_height <- attr(node, "height"); if (is.null(node_height)) node_height <- 0
    if (is.leaf(node)) { leaf_index <<- leaf_index + 1L; return(c(x = leaf_index, y = 0)) }
    children <- lapply(node, walk)
    child_x <- vapply(children, function(c) c[["x"]], numeric(1))
    child_y <- vapply(children, function(c) c[["y"]], numeric(1))
    parent_x <- mean(range(child_x))
    for (i in seq_along(children)) {
      segment_index <<- segment_index + 1L
      segment_list[[segment_index]] <<- data.frame(x = child_x[i], y = child_y[i], xend = child_x[i], yend = node_height)
    }
    segment_index <<- segment_index + 1L
    segment_list[[segment_index]] <<- data.frame(x = min(child_x), y = node_height, xend = max(child_x), yend = node_height)
    c(x = parent_x, y = node_height)
  }
  walk(dendrogram_object); do.call(rbind, segment_list)
}
module_tree_segments <- dendrogram_segments(module_tree)
module_order_labels <- labels(as.dendrogram(module_tree))
module_leaf_labels <- sub("^ME", "", module_order_labels)
module_leaf_df <- data.frame(x = seq_along(module_leaf_labels), y = 0, Module = module_leaf_labels)

p_f5 <- ggplot() +
  geom_segment(data = module_tree_segments, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "#46515D", linewidth = 0.38, lineend = "square") +
  geom_point(data = module_leaf_df, aes(x = x, y = 0, fill = Module),
             shape = 21, size = 3.0, stroke = 0.35, colour = "#FFFFFF") +
  scale_fill_manual(values = module_palette, guide = "none") +
  scale_x_continuous(breaks = seq_along(module_leaf_labels), labels = module_leaf_labels,
                     expand = expansion(mult = c(0.035, 0.035))) +
  coord_cartesian(ylim = c(-0.055 * max(module_tree$height), max(module_tree$height) * 1.02), clip = "off") +
  labs(title = "Module eigengene clustering", subtitle = sprintf("Relationships among all %d locked WGCNA modules", length(module_leaf_labels))) +
  labs(x = "Module", y = "Eigengene dissimilarity") +
  theme_classic(base_size = 7.3, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, colour = "#333333", size = 6.8),
    axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
    axis.title = element_text(size = 7.3), axis.text.y = element_text(size = 7.0, colour = "#333333"),
    plot.title = element_text(size = 8.7, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 7.0, colour = "#555555", hjust = 0.5),
    plot.margin = margin(7, 8, 10, 8)
  )

# ============================================================
# F6 GS vs MM scatter plot (significant modules)
# ============================================================
gs_mm <- gs_mm_all[gs_mm_all$Module %in% sig_colors, ]
gs_mm$IsHub <- abs(gs_mm$MM) > 0.8 & abs(gs_mm$GS) > 0.2
p_f6 <- ggplot(gs_mm, aes(MM, GS)) +
  geom_point(aes(colour = IsHub), alpha = 0.6, size = 0.9) +
  scale_colour_manual(values = c(`FALSE` = "#A8ADB4", `TRUE` = pal$red), guide = "none") +
  geom_vline(xintercept = c(-0.8, 0.8), linetype = 2, colour = pal$grey, linewidth = 0.35) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, colour = pal$grey, linewidth = 0.35) +
  facet_wrap(~Module, nrow = 1) +
  labs(title = "Gene significance vs module membership",
       subtitle = "Hub genes: |MM| > 0.8 and |GS| > 0.2 (red)",
       x = "Module membership (MM)", y = "Gene significance (GS)")

# ============================================================
# Assemble and export
# ============================================================
tag_panel <- function(plot, label) {
  plot + labs(tag = label) + theme(plot.tag = element_text(size = 10, face = "bold"),
                                   plot.tag.position = c(0, 1))
}

fig_main <- (tag_panel(p_f1a, "a") | tag_panel(p_f1b, "b") | tag_panel(p_f3, "c")) /
  tag_panel(p_f5, "d") / tag_panel(p_f4, "e") +
  plot_layout(heights = c(0.86, 0.76, 1.10))

W <- 183/25.4; H <- 170/25.4
suppressMessages(library(ragg))
agg_png(file.path(out_dir, "Figure2_WGCNA.png"), width = W, height = H, units = "in", res = 300, background = "white")
print(fig_main); dev.off()

# Export a single figure
agg_png(file.path(out_dir, "F1_soft_threshold.png"), width = 7, height = 3.2, units = "in", res = 300, background = "white")
print(tag_panel(p_f1a, "a") | tag_panel(p_f1b, "b")); dev.off()
agg_png(file.path(out_dir, "F3_module_sizes.png"), width = 4.5, height = 3, units = "in", res = 300, background = "white")
print(p_f3); dev.off()
agg_png(file.path(out_dir, "F4_module_trait.png"), width = 7, height = 4, units = "in", res = 300, background = "white")
print(p_f4); dev.off()
agg_png(file.path(out_dir, "F5_eigengene_tree.png"), width = 7, height = 3, units = "in", res = 300, background = "white")
print(p_f5); dev.off()
agg_png(file.path(out_dir, "F6_GSMM_scatter.png"), width = 8, height = 2.6, units = "in", res = 300, background = "white")
print(p_f6); dev.off()

# ============================================================
# S1 sample clustering + F2 gene tree (base R, pdf)
# ============================================================
pdf(file.path(out_dir, "S1_sample_clustering.pdf"), width = 7, height = 4)
plot(sampleTree, main = "Sample clustering", xlab = "", sub = "", labels = FALSE)
abline(h = 150, col = "red", lwd = 1.5, lty = 2)
dev.off()

pdf(file.path(out_dir, "F2_gene_dendrogram.pdf"), width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

cat("\nAll figures generated.\n")
