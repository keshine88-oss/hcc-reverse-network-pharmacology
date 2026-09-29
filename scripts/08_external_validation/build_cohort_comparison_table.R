PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
library(survival)

cox_one <- function(d, g) {
  sub <- d[!is.na(d[[g]]), c("T", "E", g)]
  if (nrow(sub) < 20) return(c(HR = NA, P = NA))
  med <- median(sub[[g]])
  sub$grp <- factor(ifelse(sub[[g]] >= med, "High", "Low"), levels = c("Low", "High"))
  s <- summary(coxph(Surv(T, E) ~ grp, data = sub))
  c(HR = s$conf.int[1, 1], P = s$coefficients[1, 5])
}

genes <- c("PTH1R","SLCO4C1","COL15A1","NTF3","PZP","CPEB3","COLEC10",
           "CFHR3","MXD3","TERT","MARCO","ANGPTL6","ECM1","EDIL3","ESM1","FAM83F","GABRD","PCDH17","PKN3")

gse <- read.csv(file.path(PROJECT_ROOT, "work/val/gse14520_allgenes.csv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))
gse$E <- as.numeric(as.character(gse$OS_status)); gse$T <- as.numeric(gse$OS_months)

liri <- read.csv(file.path(PROJECT_ROOT, "work/val/liri_allgenes.csv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))
liri$E <- as.numeric(as.character(liri$OS_status)); liri$T <- as.numeric(liri$OS_days)

out <- data.frame()
for (g in genes) {
  a <- cox_one(gse, g); b <- cox_one(liri, g)
  out <- rbind(out, data.frame(Gene = g,
    GSE14520_HR = round(a[1], 3), GSE14520_P = signif(a[2], 3),
    LIRIJP_HR = round(b[1], 3), LIRIJP_P = signif(b[2], 3)))
}
print(out, row.names = FALSE)
write.csv(out, file.path(PROJECT_ROOT, "work/val/two_cohort_HR_all19genes.csv"), row.names = FALSE)
cat("\nsaved: two_cohort_HR_all19genes.csv\n")
