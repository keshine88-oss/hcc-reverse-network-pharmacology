PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Candidate gene Venn diagram (DEG intersect WGCNA)
# ============================================================

suppressMessages({library(VennDiagram); library(grid)})

out_dir <- file.path(PROJECT_ROOT, "work/main/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(out_dir, "DEG_WGCNA_venn.png"), width = 6, height = 6, units = "in", res = 300)
grid.newpage()
draw.pairwise.venn(
  area1 = 4273, area2 = 3622, cross.area = 1817,
  category = c("DEG\n(4,273)", "WGCNA module genes\n(3,622)"),
  fill = c("#CF4E52", "#3F74B5"),
  alpha = 0.45,
  lty = "blank",
  cex = 1.8, cat.cex = 1.3, cat.fontface = "bold",
  fontfamily = "Arial", cat.fontfamily = "Arial",
  cat.pos = c(-20, 20), cat.dist = 0.05,
  cat.col = c("#222222", "#222222")
)
dev.off()

cat("Venn diagram generated. DEG 4273 | WGCNA 3622 | intersection 1817\n")
