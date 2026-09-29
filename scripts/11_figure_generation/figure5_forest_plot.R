PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# make_forest_validation.R - three-cohort external-validation forest plot
# Style aligned with the existing figures: Arial + theme_classic(base_size=7) + colours #CF4E52/#3F74B5
suppressMessages({library(ggplot2); library(ragg)})

out_dir <- file.path(PROJECT_ROOT, "work/val/forest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data: 7 genes x 3 cohorts (HR, low, high, P) ----
df <- data.frame(
  gene = rep(c("COL15A1","NTF3","PZP","CPEB3","COLEC10","SLCO4C1","PTH1R"), each = 3),
  cohort = rep(c("TCGA","GSE14520","LIRI-JP"), times = 7),
  HR  = c(0.58,0.853,1.342,  0.59,0.896,1.606,  0.61,0.727,2.353,
          0.63,0.795,1.043,  0.64,0.699,1.313,  1.42,1.392,0.383,
          0.70,1.076,0.814),
  low = c(0.41,0.557,0.724,  0.42,0.585,0.862,  0.43,0.473,1.217,
          0.44,0.519,0.565,  0.45,0.456,0.709,  1.01,0.906,0.196,
          0.497,0.703,0.439),
  high= c(0.82,1.307,2.488,  0.83,1.371,2.991,  0.87,1.117,4.547,
          0.89,1.217,1.927,  0.90,1.073,2.429,  2.01,2.137,0.752,
          0.998,1.647,1.510),
  P   = c(0.0021,0.4657,0.3499, 0.0026,0.6126,0.1356, 0.0054,0.1457,0.0109,
          0.0080,0.2910,0.8924, 0.0110,0.1016,0.3865, 0.0450,0.1311,0.0053,
          0.0480,0.7365,0.5150),
  stringsAsFactors = FALSE
)
df$sig <- ifelse(df$P < 0.05, "P < 0.05", "n.s.")
df$gene <- factor(df$gene, levels = rev(c("COL15A1","NTF3","PZP","CPEB3","COLEC10","SLCO4C1","PTH1R")))
df$cohort <- factor(df$cohort, levels = c("TCGA","GSE14520","LIRI-JP"))

cohort_col <- c(TCGA = "#CF4E52", GSE14520 = "#3F74B5", "LIRI-JP" = "#4C9F70")
dodge <- position_dodge(width = 0.62)

p <- ggplot(df, aes(y = gene, x = HR, colour = cohort, shape = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y",
                width = 0.22, linewidth = 0.38, position = dodge) +
  geom_point(size = 1.9, stroke = 0.5, position = dodge) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                labels = c("0.25","0.5","1","2","4")) +
  scale_colour_manual(values = cohort_col, name = NULL) +
  scale_shape_manual(values = c("P < 0.05" = 16, "n.s." = 1), name = NULL) +
  labs(x = "Hazard ratio (95% CI, log scale)", y = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(face = "italic"),
        legend.position = "top", legend.direction = "horizontal",
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 6.2),
        legend.margin = margin(0, 0, 2, 0),
        panel.grid = element_blank())

w <- 180 / 25.4; h <- 108 / 25.4
ragg::agg_png(file.path(out_dir, "Figure6_forest_validation.png"),
              width = w, height = h, units = "in", res = 300, background = "white")
print(p); dev.off()
grDevices::cairo_pdf(file.path(out_dir, "Figure6_forest_validation.pdf"),
                     width = w, height = h, family = "Arial")
print(p); dev.off()
ragg::agg_tiff(file.path(out_dir, "Figure6_forest_validation.tiff"),
               width = w, height = h, units = "in", res = 600, background = "white")
print(p); dev.off()
cat("saved forest plot\n")
