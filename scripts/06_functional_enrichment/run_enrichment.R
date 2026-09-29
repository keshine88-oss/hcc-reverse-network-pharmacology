PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# GO/KEGG enrichment analysis (on the 1,817 candidate genes; ClusterProfiler style)
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ggplot2); library(ragg)
})

out_dir <- file.path(PROJECT_ROOT, "work/main/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate <- read.csv(file.path(PROJECT_ROOT, "work/main/candidate_pool_final.csv"), stringsAsFactors = FALSE)$Gene
cat("number of candidate genes:", length(candidate), "\n")

# symbol -> Entrez
gene_map <- bitr(candidate, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez <- gene_map$ENTREZID
cat("mapped to Entrez:", length(entrez), "\n")

# ---- GO enrichment (BP/CC/MF) ----
ego_bp <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ego_cc <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "CC",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ego_mf <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "MF",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)

# ---- KEGG enrichment ----
ekegg <- enrichKEGG(entrez, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

cat("GO-BP enriched terms:", if (is.null(ego_bp)) 0 else nrow(ego_bp), "\n")
cat("GO-CC enriched terms:", if (is.null(ego_cc)) 0 else nrow(ego_cc), "\n")
cat("GO-MF enriched terms:", if (is.null(ego_mf)) 0 else nrow(ego_mf), "\n")
cat("KEGG enriched terms:", if (is.null(ekegg)) 0 else nrow(ekegg), "\n")

# ---- Save results ----
if (!is.null(ego_bp)) write.csv(as.data.frame(ego_bp), file.path(out_dir, "GO_BP_enrichment.csv"))
if (!is.null(ego_cc)) write.csv(as.data.frame(ego_cc), file.path(out_dir, "GO_CC_enrichment.csv"))
if (!is.null(ego_mf)) write.csv(as.data.frame(ego_mf), file.path(out_dir, "GO_MF_enrichment.csv"))
if (!is.null(ekegg)) write.csv(as.data.frame(ekegg), file.path(out_dir, "KEGG_enrichment.csv"))

# ---- GO bar plot (top 10 BP/CC/MF merged) ----
go_list <- list(BP = ego_bp, CC = ego_cc, MF = ego_mf)
go_list <- go_list[!sapply(go_list, is.null)]
if (length(go_list) > 0) {
  go_plot_data <- do.call(rbind, lapply(names(go_list), function(ont) {
    d <- as.data.frame(go_list[[ont]])
    d <- d[order(d$p.adjust), ][1:min(8, nrow(d)), ]
    d$Category <- ont
    d
  }))
  # Truncate over-long descriptions to avoid crowded labels (and deduplicate to avoid factor errors)
  go_plot_data$Description <- ifelse(nchar(go_plot_data$Description) > 55,
                                     paste0(substr(go_plot_data$Description, 1, 55), "..."),
                                     go_plot_data$Description)
  go_plot_data$Description <- make.unique(go_plot_data$Description)
  go_plot_data$Description <- factor(go_plot_data$Description,
                                     levels = rev(go_plot_data$Description[order(go_plot_data$Category, go_plot_data$p.adjust)]))
  p_go <- ggplot(go_plot_data, aes(Description, -log10(p.adjust), fill = Category)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(BP = "#3F74B5", CC = "#CF4E52", MF = "#4C9A62")) +
    labs(x = NULL, y = "-Log10(adjusted P)", title = "GO enrichment analysis") +
    theme_classic(base_size = 7.3, base_family = "Arial") +
    theme(axis.text.y = element_text(size = 6.5), legend.position = "right",
          plot.title = element_text(size = 8.7, face = "bold", hjust = 0.5))
  agg_png(file.path(out_dir, "GO_enrichment.png"), width = 8, height = 7, units = "in", res = 300, background = "white")
  print(p_go); dev.off()
}

# ---- KEGG bubble plot ----
if (!is.null(ekegg) && nrow(ekegg) > 0) {
  p_kegg <- dotplot(ekegg, showCategory = 20, font.size = 8) +
    labs(title = "KEGG pathway enrichment") +
    theme_classic(base_size = 7.3, base_family = "Arial") +
    theme(axis.text.y = element_text(size = 6.8),
          plot.title = element_text(size = 8.7, face = "bold", hjust = 0.5))
  agg_png(file.path(out_dir, "KEGG_enrichment.png"), width = 7, height = 6, units = "in", res = 300, background = "white")
  print(p_kegg); dev.off()
}

cat("Enrichment analysis and figures finished.\n")
