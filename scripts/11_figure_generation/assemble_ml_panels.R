PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Assemble the ML composite (Figure 3): six panels A-F
# ============================================================

suppressMessages(library(magick))

fig_dir <- file.path(PROJECT_ROOT, "work/main/figures")

img_a <- image_read(file.path(fig_dir, "M1_ROC.png"))
img_b <- image_read(file.path(fig_dir, "M2_upset_venn.png"))
img_c <- image_read(file.path(fig_dir, "M3_LASSO.png"))
img_d <- image_read(file.path(fig_dir, "M4_consensus_heatmap.png"))
img_e <- image_read(file.path(fig_dir, "M5_importance.png"))
img_f <- image_read(file.path(fig_dir, "M6_boxplot.png"))

add_label <- function(img, label) {
  image_annotate(img, label, size = 60, font = "Arial", color = "black",
                 weight = 700, location = "+30+20")
}

# Top row (three panels): ROC / Venn / LASSO (unified to equal height)
h_top <- 1400
img_a <- image_resize(img_a, paste0("x", h_top))
img_b <- image_resize(img_b, paste0("x", h_top))
img_c <- image_resize(img_c, paste0("x", h_top))
img_a <- add_label(img_a, "A")
img_b <- add_label(img_b, "B")
img_c <- add_label(img_c, "C")

# Middle row (two panels): heatmap / importance (wide, unified height)
h_mid <- 1000
img_d <- image_resize(img_d, paste0("x", h_mid))
img_e <- image_resize(img_e, paste0("x", h_mid))
img_d <- add_label(img_d, "D")
img_e <- add_label(img_e, "E")

# Bottom row (one panel): box plots
h_bot <- 1000
img_f <- image_resize(img_f, paste0("x", h_bot))
img_f <- add_label(img_f, "F")

row1 <- image_append(c(img_a, img_b, img_c), stack = FALSE)
row2 <- image_append(c(img_d, img_e), stack = FALSE)
fig3 <- image_append(c(row1, row2, img_f), stack = TRUE)
fig3 <- image_border(fig3, "white", "30x30")
image_write(fig3, file.path(fig_dir, "Figure4_ML.png"))

cat("ML composite assembled.\n")
