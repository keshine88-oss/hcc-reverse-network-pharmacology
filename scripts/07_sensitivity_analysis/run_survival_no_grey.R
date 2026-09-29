PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
suppressMessages({library(DESeq2); library(survival); library(survminer)})

# ---- Clinical data ----
clin <- read.csv(file.path(PROJECT_ROOT, "work/main/tcga_clinical.csv"), stringsAsFactors = FALSE)
# Strip the BOM
names(clin)[1] <- sub("\ufeff", "", names(clin)[1])
# Take the first diagnosis per patient
clin <- clin[!duplicated(clin$submitter_id), ]
# Compute OS.time / OS
clin$OS.time <- ifelse(clin$vital_status == "Dead", clin$days_to_death,
                       ifelse(clin$vital_status == "Alive", clin$days_to_last_follow_up, NA))
clin$OS <- ifelse(clin$vital_status == "Dead", 1,
                  ifelse(clin$vital_status == "Alive", 0, NA))
clin <- clin[!is.na(clin$OS.time) & !is.na(clin$OS) & clin$OS.time > 0, ]

# ---- Expression data (VST, tumour samples) ----
dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
condition <- colData(dds)$condition
tumor_samps <- colnames(mat)[condition == "Tumor"]
patient_id <- substr(tumor_samps, 1, 12)  # first 12 characters = patient ID

# ---- 16 consensus genes (grey module excluded; lncRNA LINC03082 removed) ----
genes <- c("ANGPTL6","CFHR3","COL15A1","CPEB3","ECM1","EDIL3","ESM1","FAM83F",
           "GABRD","MARCO","MXD3","PCDH17","PKN3","PTH1R","PZP","TERT")

# ---- Survival analysis ----
result <- data.frame(Gene = character(), logrank_P = numeric(),
                     HR = numeric(), HR_low = numeric(), HR_high = numeric(),
                     Cox_P = numeric(), stringsAsFactors = FALSE)

for (g in genes) {
  if (!g %in% rownames(mat)) next
  expr <- mat[g, tumor_samps]
  df <- data.frame(patient = patient_id, expr = expr, stringsAsFactors = FALSE)
  # Average multiple samples from the same patient
  df <- aggregate(expr ~ patient, data = df, FUN = mean)
  df <- merge(df, clin[, c("submitter_id", "OS.time", "OS")], by.x = "patient", by.y = "submitter_id")
  if (nrow(df) < 20) next
  df$group <- ifelse(df$expr >= median(df$expr), "High", "Low")
  df$group <- factor(df$group, levels = c("Low", "High"))

  # KM + log-rank
  fit <- survfit(Surv(OS.time, OS) ~ group, data = df)
  lr <- survdiff(Surv(OS.time, OS) ~ group, data = df)
  logrank_p <- 1 - pchisq(lr$chisq, 1)

  # Cox
  cox <- coxph(Surv(OS.time, OS) ~ group, data = df)
  cox_sum <- summary(cox)
  hr <- cox_sum$conf.int[1, 1]
  hr_low <- cox_sum$conf.int[1, 3]
  hr_high <- cox_sum$conf.int[1, 4]
  cox_p <- cox_sum$coefficients[1, 5]

  result <- rbind(result, data.frame(Gene = g, logrank_P = logrank_p,
                                     HR = hr, HR_low = hr_low, HR_high = hr_high,
                                     Cox_P = cox_p, stringsAsFactors = FALSE))
}

result <- result[order(result$logrank_P), ]
result$signif <- ifelse(result$logrank_P < 0.05, "*", "")
cat("=== Survival analysis of the 16 consensus genes without grey (KM + Cox) ===\n")
print(result, row.names = FALSE)

write.csv(result, file.path(PROJECT_ROOT, "work/main/survival_noGrey.csv"), row.names = FALSE)
cat("\nSignificant genes (log-rank P < 0.05):", sum(result$logrank_P < 0.05), "\n")
cat("significant genes:", paste(result$Gene[result$logrank_P < 0.05], collapse = ", "), "\n")
