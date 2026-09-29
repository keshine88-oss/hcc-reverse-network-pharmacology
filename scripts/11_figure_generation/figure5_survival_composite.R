PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# figure5_survival_composite.R - Figure 5: survival composite (a-g KM + h forest plot)
# Style aligned with the existing figures: Arial + theme_classic(base_size=7) + colours #CF4E52/#3F74B5
suppressMessages({
  library(DESeq2); library(survival); library(survminer)
  library(ggplot2); library(patchwork); library(ragg)
})

out_dir <- file.path(PROJECT_ROOT, "work/val/fig5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ================= a-g: TCGA KM curves =================
dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE); mat <- assay(vsd)
condition <- colData(dds)$condition
tumor_samps <- colnames(mat)[condition == "Tumor"]
pid <- substr(tumor_samps, 1, 12)

clin <- read.csv(file.path(PROJECT_ROOT, "work/main/tcga_clinical.csv"), stringsAsFactors = FALSE)
names(clin)[1] <- sub("\ufeff", "", names(clin)[1])
clin <- clin[!duplicated(clin$submitter_id), ]
clin$OS.time <- ifelse(clin$vital_status == "Dead", clin$days_to_death,
                       ifelse(clin$vital_status == "Alive", clin$days_to_last_follow_up, NA))
clin$OS <- ifelse(clin$vital_status == "Dead", 1, ifelse(clin$vital_status == "Alive", 0, NA))
clin_os <- clin[!is.na(clin$OS.time) & !is.na(clin$OS) & clin$OS.time > 0, ]

genes <- c("COL15A1","NTF3","PZP","CPEB3","COLEC10","SLCO4C1","PTH1R")
pal <- c(High = "#CF4E52", Low = "#3F74B5")
tags <- letters[1:7]

km_list <- list()
for (i in seq_along(genes)) {
  g <- genes[i]
  expr <- mat[g, tumor_samps]
  df <- data.frame(patient = pid, expr = expr, stringsAsFactors = FALSE)
  df <- aggregate(expr ~ patient, data = df, FUN = mean)
  df <- merge(df, clin_os[, c("submitter_id","OS.time","OS")], by.x = "patient", by.y = "submitter_id")
  df$group <- factor(ifelse(df$expr >= median(df$expr), "High", "Low"), levels = c("Low","High"))
  fit <- survfit(Surv(OS.time, OS) ~ group, data = df)
  lr <- survdiff(Surv(OS.time, OS) ~ group, data = df)
  pval <- 1 - pchisq(lr$chisq, 1)
  p_label <- ifelse(pval < 0.001, "P < 0.001", sprintf("P = %.3f", pval))

  p <- ggsurvplot(fit, data = df, palette = pal, legend.labs = c("Low","High"),
                  pval = FALSE, risk.table = FALSE, conf.int = FALSE,
                  xlab = "Time (days)", ylab = "Overall survival",
                  title = g, font.family = "Arial", fontsize = 3.2,
                  ggtheme = theme_classic(base_size = 6.5, base_family = "Arial") +
                    theme(plot.title = element_text(face = "bold.italic", size = 7.5, hjust = 0.5),
                          axis.title = element_text(size = 6.2),
                          axis.text = element_text(size = 5.6, colour = "black"),
                          legend.position = c(0.82, 0.86),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 5.4),
                          legend.key.size = unit(0.26, "cm"),
                          legend.background = element_blank(),
                          plot.margin = margin(3, 3, 2, 2)))
  p$plot <- p$plot + annotate("text", x = max(df$OS.time) * 0.55, y = 0.12,
                              label = p_label, size = 1.9, family = "Arial") +
    labs(tag = tags[i])
  km_list[[i]] <- p$plot
}

km_grid <- wrap_plots(km_list, ncol = 4)

# ================= h: three-cohort forest plot =================
fdf <- data.frame(
  gene = rep(c("COL15A1","NTF3","PZP","CPEB3","COLEC10","SLCO4C1","PTH1R"), each = 3),
  cohort = rep(c("TCGA","GSE14520","LIRI-JP"), times = 7),
  HR  = c(0.58,0.853,1.342, 0.59,0.896,1.606, 0.61,0.727,2.353, 0.63,0.795,1.043,
          0.64,0.699,1.313, 1.42,1.392,0.383, 0.70,1.076,0.814),
  low = c(0.41,0.557,0.724, 0.42,0.585,0.862, 0.43,0.473,1.217, 0.44,0.519,0.565,
          0.45,0.456,0.709, 1.01,0.906,0.196, 0.497,0.703,0.439),
  high= c(0.82,1.307,2.488, 0.83,1.371,2.991, 0.87,1.117,4.547, 0.89,1.217,1.927,
          0.90,1.073,2.429, 2.01,2.137,0.752, 0.998,1.647,1.510),
  P   = c(0.0021,0.4657,0.3499, 0.0026,0.6126,0.1356, 0.0054,0.1457,0.0109,
          0.0080,0.2910,0.8924, 0.0110,0.1016,0.3865, 0.0450,0.1311,0.0053,
          0.0480,0.7365,0.5150),
  stringsAsFactors = FALSE)
fdf$sig <- ifelse(fdf$P < 0.05, "P < 0.05", "n.s.")
fdf$gene <- factor(fdf$gene, levels = rev(c("COL15A1","NTF3","PZP","CPEB3","COLEC10","SLCO4C1","PTH1R")))
fdf$cohort <- factor(fdf$cohort, levels = c("TCGA","GSE14520","LIRI-JP"))
cohort_col <- c(TCGA = "#CF4E52", GSE14520 = "#3F74B5", "LIRI-JP" = "#4C9F70")
dodge <- position_dodge(width = 0.62)

forest <- ggplot(fdf, aes(y = gene, x = HR, colour = cohort, shape = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y",
                width = 0.22, linewidth = 0.38, position = dodge) +
  geom_point(size = 1.9, stroke = 0.5, position = dodge) +
  scale_x_log10(breaks = c(0.25,0.5,1,2,4), labels = c("0.25","0.5","1","2","4")) +
  scale_colour_manual(values = cohort_col, name = NULL) +
  scale_shape_manual(values = c("P < 0.05" = 16, "n.s." = 1), name = NULL) +
  labs(x = "Hazard ratio (95% CI, log scale)", y = NULL, tag = "h") +
  theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        axis.text = element_text(colour = "black", size = 6),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(size = 6.6),
        legend.position = "top", legend.direction = "horizontal",
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6),
        legend.margin = margin(0,0,1,0), panel.grid = element_blank(),
        plot.margin = margin(3, 8, 3, 3))

# ================= Panel assembly =================
final <- km_grid / forest +
  plot_layout(heights = c(2, 1.25)) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))

w <- 183 / 25.4; h <- 168 / 25.4
ragg::agg_png(file.path(out_dir, "Figure5_survival.png"),
              width = w, height = h, units = "in", res = 300, background = "white")
print(final); dev.off()
grDevices::cairo_pdf(file.path(out_dir, "Figure5_survival_composite.pdf"),
                     width = w, height = h, family = "Arial")
print(final); dev.off()
cat("saved composite\n")
