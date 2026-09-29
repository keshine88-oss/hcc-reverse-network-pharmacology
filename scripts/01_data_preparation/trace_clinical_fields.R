PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
suppressMessages(library(DESeq2))

dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
condition <- colData(dds)$condition
tumor_samps <- colnames(dds)[condition == "Tumor"]
cat("number of tumour samples:", length(tumor_samps), "\n")
pid <- substr(tumor_samps, 1, 12)
cat("unique patient IDs (first 12 characters) among tumour samples:", length(unique(pid)), "\n")

clin <- read.csv(file.path(PROJECT_ROOT, "work/main/tcga_clinical.csv"), stringsAsFactors = FALSE)
names(clin)[1] <- sub("\ufeff", "", names(clin)[1])
clin <- clin[!duplicated(clin$submitter_id), ]
cat("unique patients in the clinical data:", nrow(clin), "\n")
cat("tumour patient IDs matched to clinical patient IDs:", sum(unique(pid) %in% clin$submitter_id), "\n")

clin$OS.time <- ifelse(clin$vital_status == "Dead", clin$days_to_death,
                       ifelse(clin$vital_status == "Alive", clin$days_to_last_follow_up, NA))
clin$OS <- ifelse(clin$vital_status == "Dead", 1,
                  ifelse(clin$vital_status == "Alive", 0, NA))
clin_os <- clin[!is.na(clin$OS.time) & !is.na(clin$OS) & clin$OS.time > 0, ]
cat("patients with valid OS in the clinical data:", nrow(clin_os), "\n")
cat("tumour patients matched to a valid OS record:", sum(unique(pid) %in% clin_os$submitter_id), "\n")

# Sample size of the survival analysis for a specific gene
mat <- assay(vst(dds, blind = FALSE))
g <- "PTH1R"
expr <- mat[g, tumor_samps]
df <- data.frame(patient = pid, expr = expr, stringsAsFactors = FALSE)
df <- aggregate(expr ~ patient, data = df, FUN = mean)
df <- merge(df, clin_os[, c("submitter_id", "OS.time", "OS")], by.x = "patient", by.y = "submitter_id")
cat("Actual sample size for the PTH1R survival analysis:", nrow(df), "\n")
cat("of which Dead:", sum(df$OS == 1), " Alive:", sum(df$OS == 0), "\n")
