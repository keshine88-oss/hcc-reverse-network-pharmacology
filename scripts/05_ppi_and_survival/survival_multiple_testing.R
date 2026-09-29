PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
suppressMessages({library(DESeq2); library(survival)})

# ---- Clinical data (reusing the logic of survival_analysis.R) ----
clin <- read.csv(file.path(PROJECT_ROOT, "work/main/tcga_clinical.csv"), stringsAsFactors = FALSE)
names(clin)[1] <- sub("﻿", "", names(clin)[1])
clin <- clin[!duplicated(clin$submitter_id), ]
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
patient_id <- substr(tumor_samps, 1, 12)

genes <- c("TPM2","ACTG2","COL15A1","COLEC10","CPEB3","ECM1","EDIL3","ESM1",
           "GABRD","NTF3","PCDH17","PKN3","PTH1R","PZP","SLCO4C1")

result <- data.frame(Gene = character(), logrank_P = numeric(), HR = numeric(),
                     HR_low = numeric(), HR_high = numeric(), Cox_P = numeric(),
                     zph_group_P = numeric(), zph_global_P = numeric(),
                     zph_continuous_P = numeric(), stringsAsFactors = FALSE)

for (g in genes) {
  if (!g %in% rownames(mat)) next
  expr <- mat[g, tumor_samps]
  df <- data.frame(patient = patient_id, expr = expr, stringsAsFactors = FALSE)
  df <- aggregate(expr ~ patient, data = df, FUN = mean)
  df <- merge(df, clin[, c("submitter_id", "OS.time", "OS")], by.x = "patient", by.y = "submitter_id")
  if (nrow(df) < 20) next
  df$group <- ifelse(df$expr >= median(df$expr), "High", "Low")
  df$group <- factor(df$group, levels = c("Low", "High"))

  lr <- survdiff(Surv(OS.time, OS) ~ group, data = df)
  logrank_p <- 1 - pchisq(lr$chisq, 1)

  cox <- coxph(Surv(OS.time, OS) ~ group, data = df)
  s <- summary(cox)
  hr <- s$conf.int[1,1]; hr_low <- s$conf.int[1,3]; hr_high <- s$conf.int[1,4]
  cox_p <- s$coefficients[1,5]

  # PH assumption test (grouped variable; first row is the covariate, last row is GLOBAL)
  zph <- cox.zph(cox)
  zph_tab <- zph$table
  zph_group <- zph_tab[1, "p"]
  zph_global <- zph_tab[nrow(zph_tab), "p"]

  # PH assumption test (continuous expression, robustness check)
  cox_c <- coxph(Surv(OS.time, OS) ~ expr, data = df)
  zph_c <- cox.zph(cox_c)
  zph_cont <- zph_c$table[1, "p"]

  result <- rbind(result, data.frame(Gene = g, logrank_P = logrank_p, HR = hr,
                                     HR_low = hr_low, HR_high = hr_high, Cox_P = cox_p,
                                     zph_group_P = zph_group, zph_global_P = zph_global,
                                     zph_continuous_P = zph_cont))
}

result <- result[order(result$logrank_P), ]

# ---- BH multiple-testing correction (15 tests) ----
result$logrank_BH <- p.adjust(result$logrank_P, method = "BH")
result$Cox_BH <- p.adjust(result$Cox_P, method = "BH")

cat("=== Survival multiple-testing correction + PH assumption test (n =", nrow(result), "genes) ===\n\n")
print(result, row.names = FALSE, digits = 4)

write.csv(result, file.path(PROJECT_ROOT, "work/main/survival_stats_correction.csv"), row.names = FALSE)

cat("\n--- Summary ---\n")
cat("unadjusted log-rank P < 0.05: ", sum(result$logrank_P < 0.05), "\n")
cat("log-rank P < 0.05 after BH correction: ", sum(result$logrank_BH < 0.05), "\n")
cat("unadjusted Cox P < 0.05: ", sum(result$Cox_P < 0.05), "\n")
cat("Cox P < 0.05 after BH correction: ", sum(result$Cox_BH < 0.05), "\n")
cat("PH test P < 0.05 (proportional hazards violated, grouped): ", sum(result$zph_group_P < 0.05, na.rm = TRUE), "\n")
cat("PH global test P < 0.05: ", sum(result$zph_global_P < 0.05, na.rm = TRUE), "\n")
