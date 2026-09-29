PROJECT_ROOT <- Sys.getenv("HCC_PROJECT_ROOT", unset = normalizePath("."))
# ============================================================
# Step 3: candidate pool + four-algorithm ML consensus screening (final corrected version)
# Candidate pool = new DEG (scheme B) intersect new WGCNA module genes
# Data: VST of dds_B_final (mitochondrial genes removed + batch corrected)
# ============================================================

suppressMessages({
  library(DESeq2)
  library(randomForest)
  library(e1071)
  library(glmnet)
  library(xgboost)
  library(pROC)
})

SEED <- 42
set.seed(SEED)

# ---- Candidate pool ----
deg <- readLines(file.path(PROJECT_ROOT, "work/main/deg_B_final.txt"))
mc <- read.csv(file.path(PROJECT_ROOT, "work/main/wgcna_module_colors_final.csv"), stringsAsFactors = FALSE)
wg_all <- read.csv(file.path(PROJECT_ROOT, "work/main/wgcna_module_genes_final.csv"), stringsAsFactors = FALSE)$Gene
grey <- mc$Gene[mc$ModuleColor == "grey"]
wg <- setdiff(wg_all, grey)   # sensitivity analysis: remove the grey (unassigned) module
cat("grey removed:", length(grey), " remaining module genes:", length(wg), "
")
candidates <- intersect(deg, wg)
cat("DEGs:", length(deg), " WGCNA module genes:", length(wg), " candidate pool:", length(candidates), "\n")

# ---- VST expression matrix ----
dds <- readRDS(file.path(PROJECT_ROOT, "work/main/dds_B_final.rds"))
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
candidates <- candidates[candidates %in% rownames(mat)]
X <- t(mat[candidates, , drop = FALSE])  # samples x genes
condition <- setNames(as.character(colData(dds)$condition), colnames(dds))
y <- as.integer(condition[rownames(X)] == "Tumor")  # Normal=0, Tumor=1
cat("candidate genes:", length(candidates), " samples:", nrow(X), " Tumour:", sum(y==1), " Normal:", sum(y==0), "\n")

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
cv_auc <- function(pf) {
  sapply(seq_along(folds), function(i) {
    te <- folds[[i]]; tr <- setdiff(seq_len(nrow(X)), te)
    prob <- pf(tr, te)
    as.numeric(pROC::auc(pROC::roc(y[te], prob, quiet = TRUE)))
  })
}

# ---- RF ----
set.seed(SEED)
rf_fit <- randomForest(x = X, y = as.factor(y), ntree = 500, mtry = floor(sqrt(ncol(X))), nodesize = 1, importance = TRUE)
rf_imp <- importance(rf_fit)[, "MeanDecreaseGini"]; rf_rank <- sort(rf_imp, decreasing = TRUE)
rf_cv <- cv_auc(function(tr, te) { set.seed(SEED); m <- randomForest(X[tr,,drop=FALSE], as.factor(y[tr]), ntree=500, mtry=floor(sqrt(ncol(X))), nodesize=1); predict(m, X[te,,drop=FALSE], type="prob")[,"1"] })
cat("RF AUC:", round(mean(rf_cv),4), "+/-", round(sd(rf_cv),4), "\n")

# ---- SVM ----
set.seed(SEED)
svm_fit <- svm(x = X, y = as.factor(y), kernel = "linear", cost = 1, scale = TRUE, probability = TRUE)
w <- t(svm_fit$SV) %*% svm_fit$coefs; svm_imp <- abs(as.numeric(w)); names(svm_imp) <- colnames(X); svm_rank <- sort(svm_imp, decreasing = TRUE)
svm_cv <- cv_auc(function(tr, te) { set.seed(SEED); m <- svm(X[tr,,drop=FALSE], as.factor(y[tr]), kernel="linear", cost=1, scale=TRUE, probability=TRUE); attr(predict(m, X[te,,drop=FALSE], probability=TRUE), "probabilities")[,"1"] })
cat("SVM AUC:", round(mean(svm_cv),4), "+/-", round(sd(svm_cv),4), "\n")

# ---- LASSO ----
set.seed(SEED)
lasso_cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)
lambda_min <- lasso_cvfit$lambda.min
coef_min <- coef(lasso_cvfit, s = "lambda.min")[-1, ]
lasso_rank <- data.frame(gene = names(coef_min), coef = coef_min, stringsAsFactors = FALSE)
lasso_rank <- lasso_rank[order(-abs(lasso_rank$coef)), ]
lasso_genes <- lasso_rank$gene[lasso_rank$coef != 0]
lasso_cv <- cv_auc(function(tr, te) { set.seed(SEED); m <- glmnet(X[tr,,drop=FALSE], y[tr], family="binomial", alpha=1, lambda=lambda_min); as.numeric(predict(m, X[te,,drop=FALSE], type="response")) })
cat("LASSO AUC:", round(mean(lasso_cv),4), "+/-", round(sd(lasso_cv),4), " selected", length(lasso_genes), "genes lambda.min=", round(lambda_min,6), "\n")

# ---- XGBoost ----
xgb_params <- list(objective="binary:logistic", eval_metric="auc", max_depth=6, eta=0.3, nrounds=100, subsample=1, colsample_bytree=1, scale_pos_weight=sum(y==0)/sum(y==1))
set.seed(SEED)
xgb_fit <- xgb.train(params=list(objective=xgb_params$objective, eval_metric=xgb_params$eval_metric, max_depth=xgb_params$max_depth, eta=xgb_params$eta, subsample=xgb_params$subsample, colsample_bytree=xgb_params$colsample_bytree, scale_pos_weight=xgb_params$scale_pos_weight), data=xgb.DMatrix(X, label=y), nrounds=xgb_params$nrounds, verbose=0)
xgb_imp_mat <- xgb.importance(model = xgb_fit); xgb_imp <- setNames(xgb_imp_mat$Gain, xgb_imp_mat$Feature); xgb_rank <- sort(xgb_imp, decreasing = TRUE)
xgb_cv <- cv_auc(function(tr, te) { set.seed(SEED); m <- xgb.train(params=list(objective="binary:logistic", eval_metric="auc", max_depth=6, eta=0.3, subsample=1, colsample_bytree=1, scale_pos_weight=sum(y==0)/sum(y==1)), data=xgb.DMatrix(X[tr,,drop=FALSE], label=y[tr]), nrounds=100, verbose=0); predict(m, X[te,,drop=FALSE]) })
cat("XGBoost AUC:", round(mean(xgb_cv),4), "+/-", round(sd(xgb_cv),4), "\n")

# ---- Consensus ----
TOP_N <- 20
top_rf <- names(rf_rank)[1:min(TOP_N, length(rf_rank))]
top_svm <- names(svm_rank)[1:min(TOP_N, length(svm_rank))]
top_lasso <- lasso_rank$gene[1:min(TOP_N, nrow(lasso_rank))]
top_xgb <- names(xgb_rank)[1:min(TOP_N, length(xgb_rank))]
hit <- table(c(top_rf, top_svm, top_lasso, top_xgb))
consensus <- names(hit)[hit >= 2]
cat("\nConsensus genes (hit by >= 2 algorithms):", length(consensus), "\n")
print(sort(hit[hit >= 2], decreasing = TRUE))

# ---- Save ----
rank_df <- data.frame(Gene = candidates,
  RF = match(candidates, names(rf_rank)),
  SVM = match(candidates, names(svm_rank)),
  LASSO = match(candidates, lasso_rank$gene),
  XGB = match(candidates, names(xgb_rank)), stringsAsFactors = FALSE)
rank_df$hit <- as.integer(hit[rank_df$Gene]); rank_df$hit[is.na(rank_df$hit)] <- 0
rank_df <- rank_df[order(-rank_df$hit, rank_df$RF, rank_df$SVM, rank_df$XGB), ]
write.csv(rank_df, file.path(PROJECT_ROOT, "work/main/consensus_ranking_noGrey.csv"), row.names = FALSE)
write.csv(data.frame(Gene = consensus), file.path(PROJECT_ROOT, "work/main/consensus_genes_noGrey.csv"), row.names = FALSE)
write.csv(data.frame(Algorithm=c("RF","SVM","LASSO","XGBoost"),
  AUC_mean=c(mean(rf_cv),mean(svm_cv),mean(lasso_cv),mean(xgb_cv)),
  AUC_sd=c(sd(rf_cv),sd(svm_cv),sd(lasso_cv),sd(xgb_cv))), file.path(PROJECT_ROOT, "work/main/cv_auc_noGrey.csv"), row.names = FALSE)
write.csv(data.frame(Gene = candidates), file.path(PROJECT_ROOT, "work/main/candidate_pool_noGrey.csv"), row.names = FALSE)

cat("\nDone. Candidate pool", length(candidates), "consensus genes", length(consensus), "\n")
