PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Assemble the DEG composite (Figure 1): volcano + MA + PCA + heatmap
# ============================================================

suppressMessages(library(magick))

fig_dir <- file.path(PROJECT_ROOT, "work/main/figures")
out_dir <- file.path(PROJECT_ROOT, "work/main/figures")

img_a <- image_read(file.path(fig_dir, "D1_volcano.png"))
img_b <- image_read(file.path(fig_dir, "D2_MA.png"))
img_c <- image_read(file.path(fig_dir, "D4_PCA.png"))
img_d <- image_read(file.path(fig_dir, "D3_top50_heatmap.png"))

target_h <- 1600
img_a <- image_resize(img_a, paste0("x", target_h))
img_b <- image_resize(img_b, paste0("x", target_h))
img_c <- image_resize(img_c, paste0("x", target_h))
img_d <- image_resize(img_d, paste0("x", target_h))

add_label <- function(img, label) {
  image_annotate(img, label, size = 60, font = "Arial", color = "black",
                 weight = 700, location = "+30+20")
}
img_a <- add_label(img_a, "A")
img_b <- add_label(img_b, "B")
img_c <- add_label(img_c, "C")
img_d <- add_label(img_d, "D")

row1 <- image_append(c(img_a, img_b), stack = FALSE)
row2 <- image_append(c(img_c, img_d), stack = FALSE)
fig1 <- image_append(c(row1, row2), stack = TRUE)
fig1 <- image_border(fig1, "white", "30x30")
image_write(fig1, file.path(out_dir, "Figure1_DEG.png"))

cat("DEG composite assembled.\n")
