PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
suppressMessages({library(DESeq2); library(survival); library(survminer)})

out_dir <- file.path(PROJECT_ROOT, "work/main/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data ----
dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
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

# the 7 significant genes
genes <- c("COL15A1", "NTF3", "PZP", "CPEB3", "COLEC10", "SLCO4C1", "PTH1R")

pal <- c(High = "#CF4E52", Low = "#3F74B5")

plot_list <- list()
for (g in genes) {
  expr <- mat[g, tumor_samps]
  df <- data.frame(patient = pid, expr = expr, stringsAsFactors = FALSE)
  df <- aggregate(expr ~ patient, data = df, FUN = mean)
  df <- merge(df, clin_os[, c("submitter_id", "OS.time", "OS")], by.x = "patient", by.y = "submitter_id")
  df$group <- ifelse(df$expr >= median(df$expr), "High", "Low")
  df$group <- factor(df$group, levels = c("Low", "High"))
  fit <- survfit(Surv(OS.time, OS) ~ group, data = df)
  lr <- survdiff(Surv(OS.time, OS) ~ group, data = df)
  pval <- 1 - pchisq(lr$chisq, 1)
  p_label <- ifelse(pval < 0.001, "P < 0.001", sprintf("P = %.3f", pval))

  p <- ggsurvplot(fit, data = df, palette = pal, legend.labs = c("Low", "High"),
                  pval = FALSE, risk.table = FALSE, conf.int = FALSE,
                  xlab = "Time (days)", ylab = "Overall survival",
                  title = g, font.family = "Arial", fontsize = 4,
                  ggtheme = theme_classic(base_size = 7, base_family = "Arial") +
                    theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
                          legend.position = c(0.8, 0.85),
                          legend.title = element_blank()))
  # Manually add P-value annotations
  p$plot <- p$plot + annotate("text", x = max(df$OS.time) * 0.5, y = 0.15,
                              label = p_label, size = 2.5, family = "Arial")
  plot_list[[g]] <- p$plot
}

# Assemble a 2x4 grid (7 panels + 1 empty)
suppressMessages(library(patchwork))
combined <- wrap_plots(plot_list, ncol = 4, nrow = 2)
suppressMessages(library(ragg))
agg_png(file.path(out_dir, "KM_survival_7genes.png"), width = 11, height = 5.5, units = "in", res = 300, background = "white")
print(combined); dev.off()

cat("KM survival curves generated (7 significant genes).\n")
