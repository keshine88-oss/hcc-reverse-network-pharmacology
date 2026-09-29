PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
library(survival)
d <- read.csv(file.path(PROJECT_ROOT, "work/val/liri_allgenes.csv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))
genes <- c("PTH1R","SLCO4C1","COL15A1","NTF3","PZP","CPEB3","COLEC10","CFHR3","MXD3","TERT","MARCO","ANGPTL6","ECM1","EDIL3","ESM1","FAM83F","GABRD","PCDH17","PKN3")
d$OS_event <- as.numeric(as.character(d$OS_status)); d$OS_days <- as.numeric(d$OS_days)
cat("=== LIRI-JP per-gene external validation (n =", nrow(d), ", events =", sum(d$OS_event), ") ===\n\n")
for (g in genes) {
  if (!g %in% names(d)) { cat(sprintf("%-9s (not covered by cohort)\n", g)); next }
  sub <- d[!is.na(d[[g]]), c("OS_days","OS_event", g)]
  if (nrow(sub) < 20) { cat(sprintf("%-9s (insufficient samples)\n", g)); next }
  med <- median(sub[[g]])
  sub$grp <- factor(ifelse(sub[[g]] >= med, "High", "Low"), levels = c("Low","High"))
  s <- summary(coxph(Surv(OS_days, OS_event) ~ grp, data = sub))
  cat(sprintf("%-9s n=%3d  HR=%.3f (%.3f-%.3f)  P=%.4f  %s\n", g, nrow(sub),
              s$conf.int[1,1], s$conf.int[1,3], s$conf.int[1,4], s$coefficients[1,5],
              ifelse(s$conf.int[1,1] < 1, "protective", "risk")))
}
