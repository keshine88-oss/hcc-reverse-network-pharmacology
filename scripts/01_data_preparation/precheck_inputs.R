PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
cm <- readRDS(file.path(PROJECT_ROOT, "work/main/count.rds"))
grp <- readRDS(file.path(PROJECT_ROOT, "work/main/groups.rds"))
bt <- read.csv(file.path(PROJECT_ROOT, "work/main/batch_info.csv"), stringsAsFactors = FALSE)

cat("=== 1. data integrity ===\n")
cat("count dim:", nrow(cm), "x", ncol(cm), "\n")
cat("NA count:", sum(is.na(cm)), "\n")
cat("negative count:", sum(cm < 0), "\n")
cat("all integer:", all(cm == round(cm)), "\n")

cat("=== 2. sample alignment ===\n")
cat("colnames == names(grp):", identical(colnames(cm), names(grp)), "\n")
cat("groups length:", length(grp), " count cols:", ncol(cm), "\n")
cat("batch rows:", nrow(bt), "\n")
cat("batch covers count colnames:", all(colnames(cm) %in% bt$Barcode), "\n")

cat("=== 3. group distribution ===\n")
print(table(grp))

cat("=== 4. gene name pattern ===\n")
cat("MT- prefix:", sum(grepl("^MT-", rownames(cm))), "\n")
cat("RPL/RPS prefix:", sum(grepl("^RPL|^RPS", rownames(cm))), "\n")

cat("=== 5. plate distribution (top) ===\n")
print(head(sort(table(bt$Plate), decreasing = TRUE), 25))
