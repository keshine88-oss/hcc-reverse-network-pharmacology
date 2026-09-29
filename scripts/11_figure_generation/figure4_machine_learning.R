PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Full ML figure set (M1-M6, unified style)
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2); library(randomForest); library(e1071)
  library(glmnet); library(xgboost); library(pROC)
  library(ggplot2); library(patchwork); library(pheatmap); library(VennDiagram)
})

out_dir <- file.path(PROJECT_ROOT, "work/main/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal <- list(blue = "#3F74B5", red = "#CF4E52", grey = "#737B88", dark = "#222222", light = "#A8ADB4")
theme_pub <- function(base_size = 7.3) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 0.3, colour = "#333333"),
          plot.title = element_text(size = 8.7, face = "bold", hjust = 0.5, margin = margin(b = 3)),
          plot.subtitle = element_text(size = 7.0, colour = "#555555", hjust = 0.5),
          panel.grid = element_blank(), plot.margin = margin(7, 8, 7, 8))
}
theme_set(theme_pub())

SEED <- 42
set.seed(SEED)

# ---- Data preparation ----
candidates <- read.csv(file.path(PROJECT_ROOT, "work/main/candidate_pool_final.csv"), stringsAsFactors = FALSE)$Gene
consensus <- read.csv(file.path(PROJECT_ROOT, "work/main/consensus_genes_final.csv"), stringsAsFactors = FALSE)$Gene
dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
X <- t(mat[candidates, , drop = FALSE])
condition <- setNames(as.character(colData(dds)$condition), colnames(dds))
y <- as.integer(condition[rownames(X)] == "Tumor")

# ---- Stratified 5-fold ----
make_folds <- function(y, k, seed) {
  set.seed(seed); folds <- vector("list", k)
  for (cls in sort(unique(y))) {
    idx <- sample(which(y == cls)); n <- length(idx)
    sizes <- rep(floor(n/k), k); if (n %% k > 0) sizes[1:(n%%k)] <- sizes[1:(n%%k)] + 1
    st <- 1
    for (i in 1:k) { folds[[i]] <- c(folds[[i]], idx[st:(st+sizes[i]-1)]); st <- st + sizes[i] }
  }
  folds
}
folds <- make_folds(y, 5, SEED)

# Collect out-of-fold predicted probabilities (used for the ROC curves)
oof_predict <- function(train_fn, predict_fn) {
  prob <- numeric(nrow(X))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(nrow(X)), te)
    m <- train_fn(tr); prob[te] <- predict_fn(m, te)
  }
  prob
}

# RF
rf_train <- function(tr) { set.seed(SEED); randomForest(X[tr,,drop=FALSE], as.factor(y[tr]), ntree=500, mtry=floor(sqrt(ncol(X))), nodesize=1) }
rf_pred <- function(m, te) predict(m, X[te,,drop=FALSE], type="prob")[,"1"]
rf_oof <- oof_predict(rf_train, rf_pred)
rf_fit <- rf_train(seq_len(nrow(X)))
rf_imp <- importance(rf_fit)[, "MeanDecreaseGini"]

# SVM
svm_train <- function(tr) { set.seed(SEED); svm(X[tr,,drop=FALSE], as.factor(y[tr]), kernel="linear", cost=1, scale=TRUE, probability=TRUE) }
svm_pred <- function(m, te) attr(predict(m, X[te,,drop=FALSE], probability=TRUE), "probabilities")[,"1"]
svm_oof <- oof_predict(svm_train, svm_pred)
svm_fit <- svm_train(seq_len(nrow(X)))
w <- t(svm_fit$SV) %*% svm_fit$coefs; svm_imp <- abs(as.numeric(w)); names(svm_imp) <- colnames(X)

# LASSO (save the coefficient path)
set.seed(SEED)
lasso_cvfit <- cv.glmnet(X, y, family="binomial", alpha=1, nfolds=10)
lambda_min <- lasso_cvfit$lambda.min
lasso_fit <- glmnet(X, y, family="binomial", alpha=1, lambda=lasso_cvfit$lambda)
lasso_pred <- function(m, te) as.numeric(predict(m, X[te,,drop=FALSE], type="response"))
lasso_train <- function(tr) { set.seed(SEED); glmnet(X[tr,,drop=FALSE], y[tr], family="binomial", alpha=1, lambda=lambda_min) }
lasso_oof <- oof_predict(lasso_train, lasso_pred)
lasso_coef <- as.numeric(coef(lasso_cvfit, s="lambda.min")[-1]); names(lasso_coef) <- colnames(X)

# XGBoost
xgb_params <- function() list(objective="binary:logistic", eval_metric="auc", max_depth=6, eta=0.3, subsample=1, colsample_bytree=1, scale_pos_weight=sum(y==0)/sum(y==1))
xgb_train <- function(tr) { set.seed(SEED); xgb.train(params=xgb_params(), data=xgb.DMatrix(X[tr,,drop=FALSE], label=y[tr]), nrounds=100, verbose=0) }
xgb_pred <- function(m, te) predict(m, X[te,,drop=FALSE])
xgb_oof <- oof_predict(xgb_train, xgb_pred)
xgb_fit <- xgb_train(seq_len(nrow(X)))
xgb_imp_mat <- xgb.importance(model=xgb_fit); xgb_imp <- setNames(xgb_imp_mat$Gain, xgb_imp_mat$Feature)

# ============================================================
# M1 ROC curves
# ============================================================
roc_list <- list(
  RF = roc(y, rf_oof, quiet=TRUE), SVM = roc(y, svm_oof, quiet=TRUE),
  LASSO = roc(y, lasso_oof, quiet=TRUE), XGBoost = roc(y, xgb_oof, quiet=TRUE)
)
roc_df <- do.call(rbind, lapply(names(roc_list), function(nm) {
  r <- roc_list[[nm]]
  data.frame(Model=nm, FPR=1-r$specificities, TPR=r$sensitivities,
             AUC=sprintf("%s\nAUC=%.4f", nm, as.numeric(r$auc)))
}))
p_m1 <- ggplot(roc_df, aes(FPR, TPR, colour=Model)) +
  geom_line(linewidth=0.5) +
  geom_abline(slope=1, intercept=0, linetype=2, colour=pal$grey, linewidth=0.3) +
  scale_colour_manual(values=c(RF="#3F74B5", SVM="#4C9A62", LASSO="#CF4E52", XGBoost="#9467BD")) +
  labs(title="ROC curves", x="False positive rate", y="True positive rate") +
  theme(legend.position="right")

# ============================================================
# M2 Overlap of the four algorithms' top 20 (four-set Venn diagram)
# ============================================================
rank <- read.csv(file.path(PROJECT_ROOT, "work/main/consensus_ranking_final.csv"), stringsAsFactors = FALSE)
top20 <- function(col) rank$Gene[!is.na(rank[[col]]) & rank[[col]] <= 20]
sets <- list(RF=top20("RF"), SVM=top20("SVM"), LASSO=top20("LASSO"), XGBoost=top20("XGB"))
png(file.path(out_dir, "M2_upset_venn.png"), width=7, height=7, units="in", res=300)
grid.newpage()
draw.quad.venn(area1=length(sets$RF), area2=length(sets$SVM), area3=length(sets$LASSO), area4=length(sets$XGBoost),
  n12=length(intersect(sets$RF,sets$SVM)), n13=length(intersect(sets$RF,sets$LASSO)),
  n14=length(intersect(sets$RF,sets$XGBoost)), n23=length(intersect(sets$SVM,sets$LASSO)),
  n24=length(intersect(sets$SVM,sets$XGBoost)), n34=length(intersect(sets$LASSO,sets$XGBoost)),
  n123=length(Reduce(intersect, sets[c("RF","SVM","LASSO")])),
  n124=length(Reduce(intersect, sets[c("RF","SVM","XGBoost")])),
  n134=length(Reduce(intersect, sets[c("RF","LASSO","XGBoost")])),
  n234=length(Reduce(intersect, sets[c("SVM","LASSO","XGBoost")])),
  n1234=length(Reduce(intersect, sets)),
  category=c("RF","SVM","LASSO","XGBoost"),
  fill=c("#3F74B5","#4C9A62","#CF4E52","#9467BD"), alpha=0.45, lty="blank",
  fontfamily="sans", cat.fontfamily="sans", cex=1.4, cat.cex=1.1)
dev.off()

# ============================================================
# M3 LASSO coefficient path + CV
# ============================================================
coef_mat <- as.matrix(lasso_fit$beta)
# Keep genes with non-zero coefficients (or the top variables)
nonzero_rows <- rownames(coef_mat)[rowSums(coef_mat != 0) > 0]
if (length(nonzero_rows) > 30) {
  top_var <- names(sort(abs(lasso_coef), decreasing=TRUE))[1:30]
  coef_mat <- coef_mat[top_var, , drop=FALSE]
}
path_df <- data.frame()
for (i in seq_len(nrow(coef_mat))) {
  path_df <- rbind(path_df, data.frame(lambda=rep(lasso_fit$lambda, each=1),
                     coef=coef_mat[i,], gene=rep(rownames(coef_mat)[i], length(lasso_fit$lambda))))
}
p_m3a <- ggplot(path_df, aes(log(lambda), coef, group=gene)) +
  geom_line(linewidth=0.3, colour=pal$blue, alpha=0.6) +
  geom_vline(xintercept=log(lambda_min), linetype=2, colour=pal$red, linewidth=0.4) +
  labs(title="LASSO coefficient path", x="Log(lambda)", y="Coefficient")
cv_df <- data.frame(lambda=lasso_cvfit$lambda, cvm=lasso_cvfit$cvm,
                    cvup=lasso_cvfit$cvup, cvlo=lasso_cvfit$cvlo)
p_m3b <- ggplot(cv_df, aes(log(lambda), cvm)) +
  geom_errorbar(aes(ymin=cvlo, ymax=cvup), width=0, colour=pal$light, linewidth=0.4) +
  geom_point(size=0.8, colour=pal$blue) +
  geom_vline(xintercept=log(lasso_cvfit$lambda.min), linetype=2, colour=pal$red, linewidth=0.4) +
  geom_vline(xintercept=log(lasso_cvfit$lambda.1se), linetype=2, colour=pal$grey, linewidth=0.4) +
  labs(title="Cross-validation", x="Log(lambda)", y="Binomial deviance")
p_m3 <- p_m3a / p_m3b

# ============================================================
# M4 Consensus gene expression heatmap
# ============================================================
cons_mat <- mat[consensus, , drop=FALSE]
cons_mat <- t(scale(t(cons_mat)))
cons_mat[cons_mat > 2.5] <- 2.5; cons_mat[cons_mat < -2.5] <- -2.5
ann_col <- data.frame(Condition=condition, row.names=colnames(mat))
ann_colors <- list(Condition=c(Tumor=pal$red, Normal=pal$blue))
pheatmap(cons_mat, annotation_col=ann_col, annotation_colors=ann_colors,
         cluster_rows=TRUE, cluster_cols=TRUE, show_colnames=FALSE,
         fontsize_row=6, fontsize_col=5,
         color=colorRampPalette(c(pal$blue,"white",pal$red))(100),
         main="16 consensus key genes",
         filename=file.path(out_dir, "M4_consensus_heatmap.png"), width=7, height=5, dpi=300)

# ============================================================
# M5 Feature importance (top 15 per algorithm)
# ============================================================
top15_rf <- names(sort(rf_imp, decreasing=TRUE))[1:15]
imp_df <- rbind(
  data.frame(Gene=top15_rf, Importance=sort(rf_imp[top15_rf], decreasing=TRUE), Model="RF"),
  data.frame(Gene=names(sort(svm_imp, decreasing=TRUE))[1:15], Importance=sort(svm_imp, decreasing=TRUE)[1:15], Model="SVM"),
  data.frame(Gene=names(sort(abs(lasso_coef), decreasing=TRUE))[1:15], Importance=sort(abs(lasso_coef), decreasing=TRUE)[1:15], Model="LASSO"),
  data.frame(Gene=names(sort(xgb_imp, decreasing=TRUE))[1:15], Importance=sort(xgb_imp, decreasing=TRUE)[1:15], Model="XGBoost")
)
imp_df$Gene <- factor(imp_df$Gene, levels=rev(unique(imp_df$Gene)))
p_m5 <- ggplot(imp_df, aes(Gene, Importance, fill=Model)) +
  geom_col(width=0.7) + coord_flip() +
  facet_wrap(~Model, nrow=1, scales="free_x") +
  scale_fill_manual(values=c(RF="#3F74B5", SVM="#4C9A62", LASSO="#CF4E52", XGBoost="#9467BD")) +
  labs(title="Feature importance (top 15)", x=NULL, y="Importance") +
  theme(legend.position="none", axis.text.y=element_text(size=5.5))

# ============================================================
# M6 Consensus gene expression box plots
# ============================================================
box_df <- data.frame()
for (g in consensus) {
  box_df <- rbind(box_df, data.frame(Gene=g, Expression=mat[g,], Condition=condition))
}
box_df$Gene <- factor(box_df$Gene, levels=consensus)
p_m6 <- ggplot(box_df, aes(Condition, Expression, fill=Condition)) +
  geom_boxplot(outlier.size=0.3, linewidth=0.3, alpha=0.8) +
  scale_fill_manual(values=c(Tumor=pal$red, Normal=pal$blue)) +
  facet_wrap(~Gene, nrow=2, scales="free_y") +
  labs(title="Consensus gene expression", x=NULL, y="VST expression") +
  theme(legend.position="none", axis.text.x=element_text(size=6))

# ---- Export ggplot figures ----
suppressMessages(library(ragg))
agg_png(file.path(out_dir, "M1_ROC.png"), width=4.2, height=3.4, units="in", res=300, background="white"); print(p_m1); dev.off()
agg_png(file.path(out_dir, "M3_LASSO.png"), width=5, height=4.5, units="in", res=300, background="white"); print(p_m3); dev.off()
agg_png(file.path(out_dir, "M5_importance.png"), width=10, height=4, units="in", res=300, background="white"); print(p_m5); dev.off()
agg_png(file.path(out_dir, "M6_boxplot.png"), width=10, height=5, units="in", res=300, background="white"); print(p_m6); dev.off()

cat("ML figure set generated. Consensus genes:", length(consensus), "\n")
