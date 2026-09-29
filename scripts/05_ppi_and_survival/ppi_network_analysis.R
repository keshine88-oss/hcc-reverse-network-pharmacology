PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
suppressMessages(library(igraph))

ppi <- read.delim(file.path(PROJECT_ROOT, "work/main/ppi_string.tsv"), stringsAsFactors = FALSE)
edges <- ppi[, c("preferredName_A", "preferredName_B", "score")]
g <- graph_from_data_frame(edges[, 1:2], directed = FALSE)
E(g)$weight <- edges$score

cat("=== PPI network overview ===\n")
cat("nodes:", vcount(g), " edges:", ecount(g), "\n")

deg <- degree(g)
btw <- betweenness(g, normalized = TRUE)
cat("Degree range:", min(deg), "-", max(deg), " mean", round(mean(deg), 1), " SD", round(sd(deg), 1), "\n")

# Topological filtering (Degree mean + 1 SD)
deg_thr <- mean(deg) + sd(deg)
hub_deg <- names(deg)[deg >= deg_thr]
cat("Degree threshold (mean + 1 SD):", round(deg_thr, 1), " -> selected", length(hub_deg), "genes\n")

# Louvain community detection
set.seed(42)
comm <- cluster_louvain(g, weights = E(g)$weight)
comm_sizes <- sort(table(membership(comm)), decreasing = TRUE)
cat("number of communities:", length(comm_sizes), "\n")
cat("largest community size:", max(comm_sizes), "nodes\n")

# Degree top 20 hub genes
top20 <- names(sort(deg, decreasing = TRUE))[1:20]
cat("\n=== Degree top 20 hub genes ===\n")
print(data.frame(Gene = top20, Degree = sort(deg, decreasing = TRUE)[1:20]))

# Degree of the 16 consensus genes in the network
consensus <- c("TPM2","ACTG2","COL15A1","COLEC10","CPEB3","ECM1","EDIL3","ESM1",
               "GABRD","LINC03082","NTF3","PCDH17","PKN3","PTH1R","PZP","SLCO4C1")
cat("\n=== Degree of the 16 consensus genes in the PPI network ===\n")
cons_deg <- data.frame(Gene = consensus, Degree = NA_integer_, InNetwork = FALSE, stringsAsFactors = FALSE)
for (i in seq_along(consensus)) {
  if (consensus[i] %in% V(g)$name) {
    cons_deg$Degree[i] <- deg[consensus[i]]
    cons_deg$InNetwork[i] <- TRUE
  }
}
print(cons_deg)

# Intersection of hub genes (Degree top 20) with the consensus genes
overlap <- intersect(top20, consensus)
cat("\n=== Degree top 20 intersect the 16 consensus genes ===\n")
cat("overlap:", if (length(overlap) == 0) "none" else paste(overlap, collapse = ", "), "\n")

# Save results
saveRDS(list(g = g, deg = deg, btw = btw, comm = comm, top20 = top20,
             deg_thr = deg_thr, hub_deg = hub_deg, consensus_deg = cons_deg),
        file.path(PROJECT_ROOT, "work/main/ppi_network.rds"))
write.csv(data.frame(Gene = top20, Degree = sort(deg, decreasing = TRUE)[1:20]),
          file.path(PROJECT_ROOT, "work/main/ppi_hub_top20.csv"), row.names = FALSE)
cat("\nNetwork analysis finished; results saved.\n")
