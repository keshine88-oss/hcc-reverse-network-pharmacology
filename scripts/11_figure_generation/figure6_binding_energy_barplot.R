PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# figure6_binding_energy_barplot.R - Figure 6a: binding-energy bar chart (all 10 ligands x 2 targets)
suppressMessages({library(ggplot2); library(ragg)})

out_dir <- file.path(PROJECT_ROOT, "work/val/fig6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- data.frame(
  ligand = c("Digitoxin","Ginkgolide C","Digoxigenin","Cucurbitacin B","Digoxigenin",
             "Kamebanin","Roemerine","Chrysin","Tryptanthrin","Anisaldehyde"),
  target = c("SLCO4C1","PTH1R","SLCO4C1","PTH1R","PTH1R",
             "PTH1R","PTH1R","PTH1R","PTH1R","PTH1R"),
  score  = c(-12.3,-8.8,-9.9,-8.4,-7.9,-7.8,-7.8,-7.4,-7.3,-4.8),
  stringsAsFactors = FALSE
)
df$label <- ifelse(df$target == "PTH1R", paste0(df$ligand, "  "), df$ligand)
df <- df[order(df$score), ]
df$label <- factor(df$label, levels = df$label)
tcol <- c(PTH1R = "#CF4E52", SLCO4C1 = "#3F74B5")

p <- ggplot(df, aes(x = label, y = score, fill = target)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f", score)), hjust = 1.2, size = 2.1,
            family = "Arial", colour = "white") +
  coord_flip() +
  scale_fill_manual(values = tcol, name = NULL) +
  labs(x = NULL, y = "AutoDock Vina binding energy (kcal/mol)", tag = "a") +
  theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        axis.text = element_text(colour = "black", size = 6),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(size = 6.6),
        legend.position = "top", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6), legend.margin = margin(0,0,1,0),
        panel.grid = element_blank(), plot.margin = margin(3, 6, 3, 3),
        plot.tag = element_text(size = 8, face = "bold"))

w <- 180 / 25.4; h <- 68 / 25.4
ragg::agg_png(file.path(out_dir, "Figure6a_binding_energy.png"),
              width = w, height = h, units = "in", res = 300, background = "white")
print(p); dev.off()
cat("saved bar chart\n")
