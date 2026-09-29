PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# enrich_noGrey_export.R - export the complete enrichment results after grey removal (for the supplement)
# NOTE: R cannot write to paths containing non-ASCII characters, so results are written to an ASCII directory and copied back into the project
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
options(timeout = 600)

OUT <- file.path(PROJECT_ROOT, "work/val/enrich_noGrey")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cand <- read.csv(file.path(PROJECT_ROOT, "work/main/candidate_pool_noGrey.csv"), stringsAsFactors = FALSE)$Gene
cat("number of candidate genes:", length(cand), "\n")

gm <- bitr(cand, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("mapped to Entrez:", nrow(gm), "\n")

summ <- data.frame()
for (ont in c("BP", "CC", "MF")) {
  e <- enrichGO(gm$ENTREZID, OrgDb = org.Hs.eg.db, ont = ont,
                pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  df <- as.data.frame(e)
  cat("GO", ont, "significant terms:", nrow(df), "\n")
  if (nrow(df) > 0) {
    df <- df[order(df$p.adjust), ]
    write.csv(df, file.path(OUT, sprintf("enrich_noGrey_GO_%s.csv", ont)), row.names = FALSE)
  }
  summ <- rbind(summ, data.frame(Ontology = paste0("GO_", ont), Terms = nrow(df)))
}

k <- enrichKEGG(gm$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dk <- as.data.frame(k)
cat("KEGG significant pathways:", nrow(dk), "\n")
if (nrow(dk) > 0) {
  dk <- dk[order(dk$p.adjust), ]
  write.csv(dk, file.path(OUT, "enrich_noGrey_KEGG.csv"), row.names = FALSE)
}
summ <- rbind(summ, data.frame(Ontology = "KEGG", Terms = nrow(dk)))

print(summ, row.names = FALSE)
write.csv(summ, file.path(OUT, "enrich_noGrey_summary.csv"), row.names = FALSE)

cat("\n--- GO BP top3 ---\n")
print(head(read.csv(file.path(OUT, "enrich_noGrey_GO_BP.csv"))[, c("Description", "p.adjust", "Count")], 3), row.names = FALSE)
cat("\n--- KEGG top5 ---\n")
print(head(read.csv(file.path(OUT, "enrich_noGrey_KEGG.csv"))[, c("Description", "p.adjust", "Count")], 5), row.names = FALSE)
cat("\nExport finished ->", OUT, "\n")
