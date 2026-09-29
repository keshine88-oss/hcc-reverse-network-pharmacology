options(stringsAsFactors = FALSE, scipen = 999)

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", command_args[grep("^--file=", command_args)][1])
script_file <- normalizePath(file_arg, winslash = "/", mustWork = FALSE)
source(file.path(dirname(script_file), "00_config.R"))

extra_library <- "C:/Users/ZKX/Documents/Codex/2026-08-12/qin/outputs/VPS_revision_final/01_environment/R_library/4.5"
if (dir.exists(extra_library)) .libPaths(c(extra_library, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ggplotify)
  library(ragg)
})

project_root <- config$project_root
supplementary_root <- file.path(config$dirs$figures, "supplementary")
source_root <- file.path(config$dirs$figures, "source_data", "Supplementary_Figure_S2")
dir.create(supplementary_root, recursive = TRUE, showWarnings = FALSE)
dir.create(source_root, recursive = TRUE, showWarnings = FALSE)

network_candidates <- c(
  file.path(project_root, "03_WGCNA", "Results", "WGCNA_network.RData"),
  file.path(project_root, "Archive_To_Delete_After_Check", "Old_Manuscripts", "WGCNA_network.RData")
)
network_file <- network_candidates[file.exists(network_candidates)][1]
if (is.na(network_file)) stop("Locked WGCNA_network.RData was not found.", call. = FALSE)
required_objects <- c("net", "moduleColors", "moduleTraitCor", "moduleTraitPvalue", "sft", "power", "gs_mm_all")
wgcna_env <- new.env(parent = emptyenv())
load(network_file, envir = wgcna_env)
missing_objects <- setdiff(required_objects, ls(wgcna_env))
if (length(missing_objects)) stop("Locked WGCNA object is missing: ", paste(missing_objects, collapse = ", "), call. = FALSE)
for (object_name in required_objects) assign(object_name, get(object_name, envir = wgcna_env), envir = environment())

pal <- list(
  blue = "#3F74B5", blue_light = "#A9C7E7", red = "#CF4E52",
  grey = "#737B88", light = "#E8EDF3", dark = "#222222"
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
      plot.tag = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(7, 8, 7, 8)
    )
}
theme_set(theme_pub())

# Panels a-b: the selected power must balance scale-free fit and connectivity.
fit <- as.data.frame(sft$fitIndices)
fit$power <- fit[, 1]
fit$signed_R2 <- -sign(fit[, 3]) * fit[, 2]
fit$mean_connectivity <- fit[, 5]
fit$selected <- fit$power == power

p_a <- ggplot(fit, aes(power, signed_R2)) +
  geom_hline(yintercept = 0.80, linetype = 2, linewidth = 0.35, colour = pal$grey) +
  geom_line(linewidth = 0.45, colour = pal$blue) +
  geom_point(aes(fill = selected), shape = 21, size = 2.0, stroke = 0.35, colour = "white") +
  geom_text(data = fit[fit$power %in% c(1:6, 10, 15, 20), ], aes(label = power), nudge_y = 0.035, size = 2.05, family = "Arial", colour = pal$dark) +
  scale_fill_manual(values = c(`FALSE` = pal$blue, `TRUE` = pal$red), guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  coord_cartesian(ylim = c(min(0, min(fit$signed_R2) - 0.05), 1.02), clip = "off") +
  labs(title = "Scale-free topology fit", subtitle = sprintf("Selected power = %d; signed R² = %.3f", power, fit$signed_R2[fit$selected]), x = "Soft-thresholding power", y = "Signed R²")

p_b <- ggplot(fit, aes(power, mean_connectivity)) +
  geom_line(linewidth = 0.45, colour = pal$blue) +
  geom_point(aes(fill = selected), shape = 21, size = 2.0, stroke = 0.35, colour = "white") +
  geom_text(data = fit[fit$power %in% c(1:6, 10, 15, 20), ], aes(label = power), nudge_y = max(fit$mean_connectivity) * 0.035, size = 2.05, family = "Arial", colour = pal$dark) +
  scale_fill_manual(values = c(`FALSE` = pal$blue, `TRUE` = pal$red), guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  expand_limits(y = max(fit$mean_connectivity) * 1.11) +
  labs(title = "Mean connectivity", subtitle = sprintf("Mean connectivity at power %d = %.1f", power, fit$mean_connectivity[fit$selected]), x = "Soft-thresholding power", y = "Mean connectivity")

# Dendrogram panel: preserve the locked tree and module assignment without rerunning WGCNA.
gene_tree <- net$dendrograms[[1]]
block_genes <- net$blockGenes[[1]]
tree_colors <- moduleColors[block_genes]
tree_height <- max(gene_tree$height)
module_palette <- c(
  black = "#242424", blue = "#3F74B5", brown = "#8C5A3C", green = "#4C9A62",
  greenyellow = "#A6C94A", grey = "#A8ADB4", magenta = "#B45C9A", pink = "#E58AA8",
  purple = "#7E60A8", red = "#CF4E52", turquoise = "#31A6A1", yellow = "#E6B84A"
)
dendrogram_segments <- function(hclust_object) {
  dendrogram_object <- as.dendrogram(hclust_object)
  leaf_index <- 0L
  segment_index <- 0L
  segment_list <- list()
  walk <- function(node) {
    node_height <- attr(node, "height")
    if (is.null(node_height)) node_height <- 0
    if (is.leaf(node)) {
      leaf_index <<- leaf_index + 1L
      return(c(x = leaf_index, y = 0))
    }
    children <- lapply(node, walk)
    child_x <- vapply(children, function(child) child[["x"]], numeric(1))
    child_y <- vapply(children, function(child) child[["y"]], numeric(1))
    parent_x <- mean(range(child_x))
    for (i in seq_along(children)) {
      segment_index <<- segment_index + 1L
      segment_list[[segment_index]] <<- data.frame(
        x = child_x[i], y = child_y[i], xend = child_x[i], yend = node_height
      )
    }
    segment_index <<- segment_index + 1L
    segment_list[[segment_index]] <<- data.frame(
      x = min(child_x), y = node_height, xend = max(child_x), yend = node_height
    )
    c(x = parent_x, y = node_height)
  }
  walk(dendrogram_object)
  do.call(rbind, segment_list)
}

tree_segments <- dendrogram_segments(gene_tree)
# Split long vertical branches into height bands. This preserves every branch
# while preventing thousands of opaque low-level segments from merging into a
# black rectangle at the final 183-mm page width.
split_vertical_segments <- function(segments, cuts = c(0.25, 0.55, 0.78)) {
  pieces <- vector("list", nrow(segments))
  for (i in seq_len(nrow(segments))) {
    row <- segments[i, ]
    if (isTRUE(all.equal(row$x, row$xend)) && row$yend > row$y) {
      boundaries <- sort(unique(c(row$y, cuts[cuts > row$y & cuts < row$yend], row$yend)))
      pieces[[i]] <- data.frame(
        x = row$x,
        y = head(boundaries, -1),
        xend = row$xend,
        yend = tail(boundaries, -1)
      )
    } else {
      pieces[[i]] <- row
    }
  }
  result <- do.call(rbind, pieces)
  result$height_band <- cut(
    (result$y + result$yend) / 2,
    breaks = c(-Inf, cuts, Inf),
    labels = c("terminal", "lower", "middle", "backbone")
  )
  result
}
tree_segments <- split_vertical_segments(tree_segments)
# The full 5,000-gene dendrogram is retained as source data. For the publication
# panel, show the biologically interpretable 12-module eigengene hierarchy.
module_eigengenes <- net$MEs
eigengene_module_ids <- sub("^ME", "", names(module_eigengenes))
module_id_map <- c(
  `0` = "grey", `1` = "turquoise", `2` = "blue", `3` = "brown",
  `4` = "yellow", `5` = "green", `6` = "red", `7` = "black",
  `8` = "pink", `9` = "magenta", `10` = "purple", `11` = "greenyellow"
)
names(module_eigengenes) <- paste0("ME", unname(module_id_map[eigengene_module_ids]))
module_tree <- stats::hclust(stats::as.dist(1 - stats::cor(module_eigengenes, use = "pairwise.complete.obs")), method = "average")
module_tree_segments <- dendrogram_segments(module_tree)
module_order_labels <- labels(as.dendrogram(module_tree))
module_leaf_labels <- sub("^ME", "", module_order_labels)
module_leaf_df <- data.frame(
  x = seq_along(module_leaf_labels),
  y = -0.035 * max(module_tree$height),
  Module = module_leaf_labels
)
p_c <- ggplot() +
  geom_segment(
    data = module_tree_segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "#46515D", linewidth = 0.38, lineend = "square"
  ) +
  geom_point(
    data = module_leaf_df,
    aes(x = x, y = 0, fill = Module),
    shape = 21, size = 3.0, stroke = 0.35, colour = "#FFFFFF"
  ) +
  scale_fill_manual(values = module_palette, guide = "none") +
  scale_x_continuous(
    breaks = seq_along(module_leaf_labels), labels = module_leaf_labels,
    expand = expansion(mult = c(0.035, 0.035))
  ) +
  coord_cartesian(ylim = c(-0.055 * max(module_tree$height), max(module_tree$height) * 1.02), clip = "off") +
  labs(title = "Module eigengene clustering", subtitle = "Relationships among all 12 locked WGCNA modules") +
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

# Panel d: show the complete two-column module-trait matrix with both r and P.
module_names <- sub("^ME", "", rownames(moduleTraitCor))
module_order <- module_names[order(moduleTraitCor[, 1], decreasing = TRUE)]
heat <- rbind(
  data.frame(Module = module_names, Trait = "Tumor", Correlation = moduleTraitCor[, 1], Pvalue = moduleTraitPvalue[, 1]),
  data.frame(Module = module_names, Trait = "Normal", Correlation = -moduleTraitCor[, 1], Pvalue = moduleTraitPvalue[, 1])
)
heat$Module <- factor(heat$Module, levels = rev(module_order))
heat$Trait <- factor(heat$Trait, levels = c("Tumor", "Normal"))
heat$Block <- ifelse(
  as.character(heat$Module) %in% module_order[seq_len(6)],
  "Higher tumor correlation", "Lower tumor correlation"
)
heat$Block <- factor(heat$Block, levels = c("Higher tumor correlation", "Lower tumor correlation"))
format_p <- function(x) ifelse(x < 0.001, format(x, scientific = TRUE, digits = 1), sprintf("%.3f", x))
heat$label <- paste0(sprintf("%.2f", heat$Correlation), "\nP=", format_p(heat$Pvalue))
heat$text_colour <- ifelse(abs(heat$Correlation) >= 0.38, "white", pal$dark)

p_d <- ggplot(heat, aes(Trait, Module, fill = Correlation)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = label, colour = text_colour), family = "Arial", size = 2.15, lineheight = 0.90) +
  scale_fill_gradient2(low = pal$blue, mid = "white", high = pal$red, midpoint = 0, limits = c(-0.7, 0.7), name = "Pearson r") +
  scale_colour_identity() +
  scale_x_discrete(position = "top") +
  facet_wrap(~Block, nrow = 1, scales = "free_y") +
  guides(fill = guide_colourbar(
    direction = "horizontal", title.position = "top", title.hjust = 0.5,
    barwidth = unit(48, "mm"), barheight = unit(3.2, "mm"), ticks.colour = "white"
  )) +
  labs(title = "Complete module–trait matrix", subtitle = "All 12 modules; each cell reports Pearson r and nominal P value", x = NULL, y = NULL) +
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

# Panel e: module sizes confirm the 5,000-gene partition and retain module-color semantics.
module_sizes <- as.data.frame(table(moduleColors), stringsAsFactors = FALSE)
names(module_sizes) <- c("Module", "Genes")
module_sizes <- module_sizes[order(module_sizes$Genes), ]
module_sizes$Module <- factor(module_sizes$Module, levels = module_sizes$Module)
p_e <- ggplot(module_sizes, aes(Genes, Module, fill = Module)) +
  geom_col(width = 0.70, colour = "#555555", linewidth = 0.18) +
  geom_text(aes(label = Genes), hjust = -0.12, family = "Arial", size = 2.05, colour = pal$dark) +
  scale_fill_manual(values = module_palette, guide = "none") +
  expand_limits(x = max(module_sizes$Genes) * 1.16) +
  labs(title = "Module sizes", subtitle = "All 5,000 input genes assigned", x = "Genes", y = NULL)

tag_panel <- function(plot, label) {
  plot + labs(tag = label) +
    theme(plot.tag = element_text(size = 10, face = "bold"),
          plot.tag.position = c(0, 1))
}
top_row <- tag_panel(p_a, "a") | tag_panel(p_b, "b") | tag_panel(p_e, "c")
tree_row <- tag_panel(p_c, "d")
heatmap_row <- tag_panel(p_d, "e")
figure_s2 <- (top_row / tree_row / heatmap_row) +
  plot_layout(heights = c(0.86, 0.76, 1.10)) &
  theme(plot.margin = margin(7, 8, 7, 8))

name <- "Supplementary_Figure_S2_WGCNA_diagnostics"
width_mm = 183
height_mm <- 170
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file.path(supplementary_root, paste0(name, ".svg")), width = width_in, height = height_in, system_fonts = list(sans = "Arial"))
} else {
  grDevices::svg(file.path(supplementary_root, paste0(name, ".svg")), width = width_in, height = height_in, family = "Arial", onefile = TRUE)
}
print(figure_s2)
dev.off()

grDevices::cairo_pdf(file.path(supplementary_root, paste0(name, ".pdf")), width = width_in, height = height_in, family = "Arial")
print(figure_s2)
dev.off()

ragg::agg_tiff(file.path(supplementary_root, paste0(name, ".tiff")), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw", background = "white")
print(figure_s2)
dev.off()

ragg::agg_png(file.path(supplementary_root, paste0(name, "_preview.png")), width = width_in, height = height_in, units = "in", res = 300, background = "white")
print(figure_s2)
dev.off()

# Panel-level source data and traceability.
write.csv(fit[, c("power", "signed_R2", "mean_connectivity", "selected")], file.path(source_root, "soft_threshold_diagnostics.csv"), row.names = FALSE)
write.csv(heat[, c("Module", "Trait", "Correlation", "Pvalue")], file.path(source_root, "complete_module_trait_matrix.csv"), row.names = FALSE)
write.csv(module_sizes, file.path(source_root, "module_sizes.csv"), row.names = FALSE)
tree_source <- data.frame(
  leaf_order = seq_along(gene_tree$order),
  block_index = block_genes[gene_tree$order],
  module = tree_colors[gene_tree$order]
)
write.csv(tree_source, file.path(source_root, "dendrogram_leaf_module_assignment.csv"), row.names = FALSE)
write.csv(tree_segments, file.path(source_root, "dendrogram_full_branch_segments.csv"), row.names = FALSE)
write.csv(module_tree_segments, file.path(source_root, "module_eigengene_dendrogram_segments.csv"), row.names = FALSE)
file.copy(network_file, file.path(source_root, "WGCNA_network_locked.RData"), overwrite = TRUE)

md5_value <- unname(tools::md5sum(network_file))
qa_rows <- data.frame(
  check = c(
    "network_object_md5", "input_genes", "selected_power", "signed_scale_free_R2",
    "mean_connectivity", "modules", "significant_modules", "historical_hubs",
    "module_eigengenes_in_display_tree", "full_gene_dendrogram_segments_in_source_data",
    "figure_width_mm", "figure_height_mm", "vector_exports", "raster_export"
  ),
  value = c(
    md5_value, length(moduleColors), power, sprintf("%.6f", fit$signed_R2[fit$selected]),
    sprintf("%.6f", fit$mean_connectivity[fit$selected]), length(unique(moduleColors)),
    sum(abs(moduleTraitCor[, 1]) > 0.30 & moduleTraitPvalue[, 1] < 0.05),
    sum(gs_mm_all$IsHub), ncol(module_eigengenes), nrow(tree_segments),
    width_mm, height_mm, "SVG; PDF", "TIFF 600 dpi"
  )
)
write.csv(qa_rows, file.path(config$dirs$logs, "wgcna_figure_enhancement_numeric_QA.csv"), row.names = FALSE)

contract <- c(
  "# Supplementary Figure S2 contract",
  "",
  "Core conclusion: the locked WGCNA candidate source is technically auditable through soft-threshold diagnostics, the gene dendrogram and module assignments, the complete module–trait matrix, and module sizes.",
  "Figure archetype: asymmetric quantitative grid with a dendrogram hero panel.",
  "Backend: R only.",
  "Final size: 183 × 170 mm; SVG, PDF, 600-dpi TIFF and 300-dpi preview PNG.",
  "Panel a: signed scale-free topology fit across powers 1–20.",
  "Panel b: mean connectivity across powers 1–20.",
  "Panel c: locked module sizes.",
  "Panel d: the eigengene dendrogram for all 12 locked modules; the complete 5,000-gene dendrogram, all leaf assignments and all branch segments remain supplied as source data.",
  "Panel e: complete tumor/normal module–trait matrix with Pearson r and nominal P.",
  "Reviewer boundary: this figure documents upstream discovery provenance; it is not resampled feature-selection validation or independent prognostic evidence.",
  sprintf("Locked selected power: %d; signed R²: %.3f; mean connectivity: %.1f.", power, fit$signed_R2[fit$selected], fit$mean_connectivity[fit$selected])
)
writeLines(contract, file.path(config$dirs$logs, "wgcna_figure_enhancement_contract.md"), useBytes = TRUE)

cat("Enhanced Supplementary Figure S2 written to", supplementary_root, "\n")
