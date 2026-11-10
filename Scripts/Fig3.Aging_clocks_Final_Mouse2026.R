########
########

######## for Mouse ########
########

##### 校准去掉 ######
# 校准：训练集预测按 sample 聚合，拟合 lm(actual ~ predicted) 得 intercept/slope

######## 
######## Mouse 预处理 ########
########

conda activate ArchR
R

#########
library(Seurat)
library(FNN)
library(Matrix)

setwd("/projects/hmz-aging/Mouse_New_Datasets")
Mouse_seurat_clcl <- readRDS("Mouse_seurat_merge_clcl_2025")

# 去除 sample 名称中的 _T1 和 _T2 后缀
Mouse_seurat_clcl$sample <- gsub("_T[12]$", "", Mouse_seurat_clcl$sample)

######### 打印总体有哪些 sample / celltype #########
cat("\n========== 全部 sample（去掉 _T1/_T2 后缀后）==========\n")
print(sort(unique(Mouse_seurat_clcl$sample)))
cat("\n每个 sample 的细胞数:\n")
print(table(Mouse_seurat_clcl$sample))
cat("\n每个 celltype 的细胞数:\n")
print(table(Mouse_seurat_clcl$celltype))

######### 按 celltype 拆分 #########
Mouse_seurat_clcl_list <- SplitObject(Mouse_seurat_clcl, split.by = "celltype")

######### 打印每个 celltype 下有哪些 sample #########
cat("\n========== 每个 celltype 下的 sample ==========\n")
for (ct in names(Mouse_seurat_clcl_list)) {
    cat(sprintf("\n--- %s ---\n", ct))
    print(table(Mouse_seurat_clcl_list[[ct]]$sample))
}


######################################################################
######## 辅助函数 1：标准化 + PCA ########
######################################################################
NormalizeAndPCA <- function(
  seu,
  norm.method = "LogNormalize",
  scale.factor = 1e4,
  nfeatures = 2000,
  npcs = 30
) {
  # 1. 标准化
  seu <- NormalizeData(
    object = seu,
    normalization.method = norm.method,
    scale.factor = scale.factor
  )

  # 2. 识别高变基因
  seu <- FindVariableFeatures(
    object = seu,
    selection.method = "vst",
    nfeatures = nfeatures
  )

  # 3. 数据缩放（所有基因）
  seu <- ScaleData(
    object = seu,
    features = rownames(seu)
  )

  # 4. PCA 降维（动态调整 npcs）
  # SVD 要求: npcs 必须**严格小于** min(n_genes, n_cells)
  n_cells <- ncol(seu)
  n_genes <- nrow(seu)
  max_npcs <- min(n_genes, n_cells) - 1

  if (npcs > max_npcs) {
    npcs_adjusted <- max(1, max_npcs)  # 至少 1 个 PC
    message(sprintf("    调整 npcs: %d → %d (n_cells=%d, n_genes=%d)",
                    npcs, npcs_adjusted, n_cells, n_genes))
    npcs <- npcs_adjusted
  }

  if (npcs < 1) {
    stop(sprintf("细胞数太少 (n_cells=%d)，无法进行 PCA", n_cells))
  }

  seu <- RunPCA(
    object = seu,
    features = VariableFeatures(seu),
    npcs = npcs
  )

  return(seu)
}


######################################################################
######## 辅助函数 2：KNN smoothing -> pseudo-cell ########
######################################################################
PseudoBulkByKNN_fast <- function(
  seu,
  k = 15,
  dims = 1:30,
  assay = "RNA",
  reduction = "pca"
) {
  library(Seurat)
  library(FNN)
  library(Matrix)
  # 1. PCA embedding
  emb <- Embeddings(seu, reduction)[, dims, drop = FALSE]

  # 2. kNN（不含自身）
  knn <- get.knn(emb, k = k)
  nn_idx <- knn$nn.index
  ncells <- nrow(nn_idx)

  # 3. 构建邻接稀疏矩阵 A：对每个目标细胞 j，把自身 j 和它的 k 个邻居 i 置 1
  rows <- as.integer(c(t(nn_idx), seq_len(ncells)))
  cols <- rep(seq_len(ncells), each = k + 1)
  A <- sparseMatrix(i = rows, j = cols, x = 1,
                    dims = c(ncells, ncells))

  # 4. 原始 counts（genes × cells）
  counts_mat <- GetAssayData(seu, assay = assay, layer = "counts")

  # 5. 一步乘法得到伪-bulk counts
  pseudo_counts <- counts_mat %*% A

  # 6. 新 Seurat 对象
  new_seu <- CreateSeuratObject(
    counts = pseudo_counts,
    meta.data = seu@meta.data,
    assay = assay
  )
  return(new_seu)
}


######################################################################
######## 主循环：每个 celltype 按 sample smooth 成 pseudo-cell ########
######################################################################
Mouse_seurat_knn_pseudobulk <- list()

for(i in 1:length(Mouse_seurat_clcl_list)){
    print(i)
    print(names(Mouse_seurat_clcl_list)[i])
    #####
    tmp = Mouse_seurat_clcl_list[[i]]
    #####
    tmp_sp = SplitObject(tmp, split.by = "sample")

    ##### 打印该 celltype 拆分出的 sample #####
    cat(sprintf("  [%s] samples: %s\n",
                names(Mouse_seurat_clcl_list)[i],
                paste(names(tmp_sp), collapse = ", ")))

    ##### 过滤细胞数<2的样本 #######
    sample_ncells <- sapply(tmp_sp, ncol)
    removed_samples <- names(sample_ncells)[sample_ncells < 2]
    if(length(removed_samples) > 0) {
        message(sprintf("  %s: 去掉细胞数<2的样本: %s", names(Mouse_seurat_clcl_list[i]), paste(removed_samples, collapse=", ")))
        tmp_sp <- tmp_sp[sample_ncells >= 2]
    }
    if(length(tmp_sp) == 0) {
        message(sprintf("  %s: 所有样本都被过滤，跳过", names(Mouse_seurat_clcl_list[i])))
        next
    }
    ##### first normalize and get PCA #######
    if(names(Mouse_seurat_clcl_list[i]) %in% c("Rod","AC","HC","BC","Cone","RGC")){
        k = 15
        npcs = 25
    }
    if(names(Mouse_seurat_clcl_list[i]) %in% c("MG")){
        k = 5  # MG: K = 5
        npcs = 5
    }
    if(names(Mouse_seurat_clcl_list[i]) %in% c("RPE","Microglia")){
        k = 2  # RPE 和 Microglia: K = 2
        npcs = 5
    }
    print(paste("Cell type:", names(Mouse_seurat_clcl_list[i])))
    print(paste("K =", k))
    print(paste("npcs (initial) =", npcs))

    ##### 动态调整 npcs 和 K：检查每个样本的细胞数 #######
    sample_ncells_check <- sapply(tmp_sp, ncol)
    min_ncells <- min(sample_ncells_check)
    print(paste("Min cells per sample:", min_ncells))

    # 特殊情况：如果细胞数太少（< 5），跳过 KNN smooth
    if (min_ncells < 5) {
        message(sprintf("  警告: 最小样本细胞数 (%d) < 5，跳过 KNN smooth，直接 merge", min_ncells))
        tmp_sp_res_pseudo_merge = merge(tmp_sp[[1]], tmp_sp[-1], add.cell.ids = names(tmp_sp))
        Mouse_seurat_knn_pseudobulk <- c(Mouse_seurat_knn_pseudobulk, list(tmp_sp_res_pseudo_merge))
        next
    }

    # npcs 必须严格小于最小的样本细胞数
    if (npcs >= min_ncells) {
        npcs_adjusted <- max(1, min_ncells - 1)
        message(sprintf("  警告: npcs (%d) >= 最小样本细胞数 (%d)，调整为 %d",
                        npcs, min_ncells, npcs_adjusted))
        npcs <- npcs_adjusted
    }

    # K 必须小于细胞数（KNN 需要邻居）
    if (k >= min_ncells) {
        k_adjusted <- max(1, min_ncells - 1)
        message(sprintf("  警告: K (%d) >= 最小样本细胞数 (%d)，调整为 %d",
                        k, min_ncells, k_adjusted))
        k <- k_adjusted
    }

    print(paste("npcs (adjusted) =", npcs))
    print(paste("K (adjusted) =", k))

    ########
    tmp_sp_res = sapply(tmp_sp, NormalizeAndPCA, npcs=npcs)
    ##### next smooth cells to get pseudo-bulk !!!! #######
    tmp_sp_res_pseudo = sapply(tmp_sp_res, PseudoBulkByKNN_fast, k=k, dims=1:npcs)
    #####
    tmp_sp_res_pseudo_merge = merge(tmp_sp_res_pseudo[[1]], tmp_sp_res_pseudo[-1], add.cell.ids = names(tmp_sp_res_pseudo))
    #####
    Mouse_seurat_knn_pseudobulk <- c(Mouse_seurat_knn_pseudobulk, list(tmp_sp_res_pseudo_merge))
}

names(Mouse_seurat_knn_pseudobulk) <- names(Mouse_seurat_clcl_list)
saveRDS(Mouse_seurat_knn_pseudobulk, file="Mouse_seurat_knn_pseudobulk_2025")

##### 打印每个 celltype pseudobulk 后的 sample #####
cat("\n========== KNN pseudobulk 后每个 celltype 的 sample ==========\n")
for (ct in names(Mouse_seurat_knn_pseudobulk)) {
    cat(sprintf("\n--- %s ---\n", ct))
    print(table(Mouse_seurat_knn_pseudobulk[[ct]]$sample))
}

####### counts matrix ########
mat <- LayerData(Mouse_seurat_knn_pseudobulk$MG, layer = "counts")
range(mat)                     # 最小 0
all(mat@x == round(mat@x))     # TRUE 说明是整数（raw counts 求和）
head(Matrix::colSums(mat))     # 每个 pseudo-cell 的总 counts



#######
####### mouse 接下来 我要 做 nexsted leave one out ####
####### 不要 最后的那个标准化 ###
#######   去掉两处标准化：(1) 每折内的 Z-score 标准化（glmnet 也设 standardize=FALSE）
#######                   (2) 最后的样本级校准（只保留原始预测，无 predicted_calibrated）

library(Seurat)
library(glmnet)
library(Matrix)

setwd("/projects/hmz-aging/Mouse_New_Datasets")

######### 载入数据、ambient / whitelist 基因 #########
Mouse_seurat_knn_pseudobulk <- readRDS("Mouse_seurat_knn_pseudobulk_2025")

M_ambient_genes <- readRDS("Mouse_ambient_RNA_Jul22")
M_ambient_genes <- M_ambient_genes[order(M_ambient_genes$Avg, decreasing = TRUE), ]
M_ambient_genes_Top50 <- M_ambient_genes$Gene[which(M_ambient_genes$Avg > 0.0025)]
M_ambient_genes_Top50 <- sapply(strsplit(M_ambient_genes_Top50, split = "~~"), function(x) x[[2]])

load("Mouse_AgingMap_Genes_Jul9")
load("SID_Mouse_Jul9")
whitelist_genes <- c(Mouse_AgingMap_Genes, SID_Mouse)
whitelist_genes <- whitelist_genes[!duplicated(whitelist_genes)]
cat(sprintf("ambient genes: %d, whitelist genes: %d\n", length(M_ambient_genes_Top50), length(whitelist_genes)))

######### 给每个 celltype 加 age 列（从 sample: Wk<age>_Rep<n> 提取；106 -> 108）#########
for (ct in names(Mouse_seurat_knn_pseudobulk)) {
    sample_names <- Mouse_seurat_knn_pseudobulk[[ct]]$sample
    ages <- as.numeric(gsub("^Wk([0-9]+)_.*", "\\1", sample_names))
    ages[ages == 106] <- 108
    Mouse_seurat_knn_pseudobulk[[ct]]$age <- ages
    cat(sprintf("  %s: 年龄范围 = %s 周\n", ct, paste(sort(unique(ages)), collapse = ", ")))
}

######### 复用 pseudobulk / 相关性计算等辅助函数 #########
source("Mouse_AgingClock_Functions.R")


######################################################################
######## 精简版 Nested LOOCV：保留 Z-score 标准化、无校准 ########
######## 覆盖 Mouse_AgingClock_Functions.R 里的同名函数 ########
######################################################################
Nested_FeatureSelection_LOOCV <- function(seurat_obj,
                                          whitelist_genes = NULL,
                                          ambient_genes = NULL,
                                          cor_pval_cutoff = 0.05,
                                          cor_abs_cutoff = 0.2,
                                          cor_top_n = 500,
                                          min_mean_expr = 0.1,
                                          alphas = seq(0.1, 0.9, by = 0.2),
                                          nfolds_inner = 5,
                                          standardize = TRUE,
                                          exclude_ages = NULL,
                                          exclude_samples = NULL,
                                          n_cores = NULL) {
    require(glmnet)
    require(Matrix)
    require(Seurat)

    use_parallel <- TRUE
    if (is.null(n_cores)) {
        if (!requireNamespace("pbmcapply", quietly = TRUE)) { use_parallel <- FALSE; n_cores <- 1 }
        else { n_cores <- parallel::detectCores() - 1 }
    }
    if (n_cores < 1) { n_cores <- 1; use_parallel <- FALSE }
    else if (n_cores == 1) { use_parallel <- FALSE }
    else if (!requireNamespace("pbmcapply", quietly = TRUE)) { use_parallel <- FALSE; n_cores <- 1 }

    message("========================================")
    message("Nested Feature Selection + LOOCV（保留 Z-score 标准化 / 无校准）")
    message(sprintf("%s", if (use_parallel) sprintf("并行: %d 核", n_cores) else "串行"))
    message("========================================")

    # 1. 准备数据
    if (length(Layers(seurat_obj)) > 1) seurat_obj <- JoinLayers(seurat_obj)
    counts <- LayerData(seurat_obj, layer = "counts")
    cell_samples <- seurat_obj$sample
    cell_ages <- seurat_obj$age
    if (is.null(cell_samples) || is.null(cell_ages)) stop("Seurat object 缺少 'sample' 或 'age' 列")

    # 过滤指定年龄
    if (!is.null(exclude_ages)) {
        keep <- !cell_ages %in% exclude_ages
        if (sum(!keep) > 0) {
            counts <- counts[, keep, drop = FALSE]; cell_samples <- cell_samples[keep]; cell_ages <- cell_ages[keep]
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
            message(sprintf("  排除年龄 %s: 去掉 %d 个细胞", paste(exclude_ages, collapse=","), sum(!keep)))
        }
    }
    # 过滤指定样本
    if (!is.null(exclude_samples)) {
        keep <- !cell_samples %in% exclude_samples
        if (sum(!keep) > 0) {
            counts <- counts[, keep, drop = FALSE]; cell_samples <- cell_samples[keep]; cell_ages <- cell_ages[keep]
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
            message(sprintf("  排除样本 %s: 去掉 %d 个细胞", paste(exclude_samples, collapse=","), sum(!keep)))
        }
    }
    # 去 ambient genes
    if (!is.null(ambient_genes)) {
        ambient_genes <- intersect(ambient_genes, rownames(counts))
        if (length(ambient_genes) > 0) counts <- counts[!rownames(counts) %in% ambient_genes, ]
    }
    # 单细胞 depth normalization + log2（建模用；注意这不是我们要去掉的“标准化”）
    depth <- Matrix::colSums(counts)
    counts_norm <- sweep(counts, 2, depth / 1e4, FUN = "/")
    counts_log2 <- log2(1 + counts_norm)
    message(sprintf("  单细胞数: %d, 基因数: %d", ncol(counts_log2), nrow(counts_log2)))

    # 2. pseudobulk（特征选择用）
    pseudobulk_all <- Calculate_Pseudobulk_Normalized(seurat_obj, ambient_genes)
    unique_samples <- unique(cell_samples)
    sample_ages <- sapply(unique_samples, function(s) unique(cell_ages[cell_samples == s])[1])
    n_samples <- length(unique_samples)
    message(sprintf("  样本数: %d, 年龄: %s", n_samples, paste(sort(unique(sample_ages)), collapse=", ")))

    # 3. LOOCV
    process_one_loocv <- function(i) {
        test_sample <- unique_samples[i]
        train_samples <- unique_samples[-i]
        train_pseudobulk <- pseudobulk_all[, train_samples, drop = FALSE]

        # 3.1 特征选择（仅训练集）
        cor_results <- suppressMessages(Calculate_Gene_Age_Correlation(train_pseudobulk, method = "spearman"))
        cor_results_sorted <- cor_results[order(abs(cor_results$correlation), decreasing = TRUE), ]
        # 强相关基因：|cor| > cor_abs_cutoff 且 p < cor_pval_cutoff；若不足 cor_top_n，则退回用 |cor| 前 cor_top_n
        genes_strong <- cor_results_sorted$gene[abs(cor_results_sorted$correlation) > cor_abs_cutoff &
                                                cor_results_sorted$p_value < cor_pval_cutoff]
        if (length(genes_strong) >= cor_top_n) {
            sig_genes <- genes_strong
        } else {
            sig_genes <- head(cor_results_sorted$gene, cor_top_n)
        }
        if (!is.null(whitelist_genes)) {
            whitelist_available <- intersect(whitelist_genes, rownames(counts_log2))
            candidate_genes_all <- unique(c(sig_genes, whitelist_available))
        } else {
            candidate_genes_all <- sig_genes
        }
        if (length(candidate_genes_all) < 10)
            candidate_genes_all <- head(cor_results_sorted$gene, min(100, nrow(cor_results)))
        gene_mean_expr <- rowMeans(train_pseudobulk)
        candidate_genes_expr <- gene_mean_expr[candidate_genes_all]
        candidate_genes <- candidate_genes_all[!is.na(candidate_genes_expr) & candidate_genes_expr >= min_mean_expr]
        n_filtered_out_lowexpr <- sum(candidate_genes_expr < min_mean_expr, na.rm = TRUE)
        if (length(candidate_genes) == 0) stop("过滤后没有候选基因！")

        # 3.2 训练/测试单细胞矩阵（仅候选基因，log2 normalized，不做 z-score）
        train_cell_idx <- which(cell_samples %in% train_samples)
        test_cell_idx <- which(cell_samples == test_sample)
        if (length(train_cell_idx) == 0) stop("训练集没有细胞")
        if (length(test_cell_idx) == 0) stop(sprintf("测试样本 %s 没有细胞", test_sample))
        candidate_genes <- intersect(candidate_genes, rownames(counts_log2))
        X_train <- t(as.matrix(counts_log2[candidate_genes, train_cell_idx, drop = FALSE]))
        y_train <- cell_ages[train_cell_idx]
        train_samples_vec <- cell_samples[train_cell_idx]
        ## Z-score 标准化（基于训练集的 mean/sd，测试集沿用同一套，防泄漏）
        if (standardize) {
            train_means <- colMeans(X_train)
            train_sds <- apply(X_train, 2, sd)
            train_sds[train_sds == 0] <- 1
            X_train <- scale(X_train, center = train_means, scale = train_sds)
        }

        # 3.3 IPF 权重：平衡 sample 与 age
        weights <- rep(1, length(train_cell_idx))
        for (iter in 1:10) {
            sample_weight_sum <- tapply(weights, train_samples_vec, sum)
            weights <- weights * (mean(sample_weight_sum) / sample_weight_sum[train_samples_vec])
            age_weight_sum <- tapply(weights, y_train, sum)
            weights <- weights * (mean(age_weight_sum) / age_weight_sum[as.character(y_train)])
        }
        weights <- weights * length(weights) / sum(weights)

        # 3.4 内层 CV 调 alpha（glmnet standardize=FALSE，真正不做标准化）
        best_alpha <- NULL; best_cv_error <- Inf; best_lambda <- NULL
        for (alpha in alphas) {
            cv_fit <- tryCatch(
                cv.glmnet(X_train, y_train, weights = weights, alpha = alpha,
                          nfolds = nfolds_inner, family = "gaussian", standardize = FALSE),
                error = function(e) NULL)
            if (!is.null(cv_fit)) {
                min_error <- min(cv_fit$cvm)
                if (min_error < best_cv_error) { best_cv_error <- min_error; best_alpha <- alpha; best_lambda <- cv_fit$lambda.1se }
            }
        }
        if (is.null(best_alpha)) { best_alpha <- 0.5; best_lambda <- NULL }

        # 3.5 训练最终模型 + 预测测试样本（原始预测，无校准）
        final_model <- glmnet(X_train, y_train, weights = weights, alpha = best_alpha,
                              family = "gaussian", standardize = FALSE)
        X_test <- t(as.matrix(counts_log2[candidate_genes, test_cell_idx, drop = FALSE]))
        if (standardize) X_test <- scale(X_test, center = train_means, scale = train_sds)
        s_use <- if (!is.null(best_lambda)) best_lambda else min(final_model$lambda)
        pred <- predict(final_model, newx = X_test, s = s_use)

        list(test_sample = test_sample, test_cell_idx = test_cell_idx,
             predictions = as.numeric(pred),
             selected_genes = candidate_genes, best_alpha = best_alpha, best_cv_error = best_cv_error,
             n_genes = length(candidate_genes), n_genes_filtered_lowexpr = n_filtered_out_lowexpr,
             n_train_cells = length(train_cell_idx), n_test_cells = length(test_cell_idx),
             actual_age = sample_ages[i])
    }

    if (use_parallel) {
        loocv_results <- pbmcapply::pbmclapply(1:n_samples, process_one_loocv,
                                               mc.cores = n_cores, mc.set.seed = TRUE, mc.preschedule = FALSE,
                                               ignore.interactive = TRUE)
    } else {
        loocv_results <- lapply(1:n_samples, function(i) {
            message(sprintf(">>> LOOCV %d/%d <<<", i, n_samples))
            res <- process_one_loocv(i)
            message(sprintf("    %s | genes=%d | pred=%.2f, actual=%.2f",
                            res$test_sample, res$n_genes, mean(res$predictions), res$actual_age))
            res
        })
    }

    # 整合结果
    valid <- sapply(loocv_results, function(x) is.list(x) && !is.null(x$test_sample))
    if (!all(valid)) stop(sprintf("LOOCV Round %s 无效", paste(which(!valid), collapse=", ")))

    cell_predictions <- data.frame(
        cell_id = colnames(counts_log2), sample = cell_samples, actual = cell_ages,
        predicted = NA, stringsAsFactors = FALSE)
    selected_genes_list <- list(); best_alphas <- numeric(n_samples)
    for (i in 1:n_samples) {
        r <- loocv_results[[i]]
        cell_predictions$predicted[r$test_cell_idx] <- r$predictions
        selected_genes_list[[r$test_sample]] <- r$selected_genes
        best_alphas[i] <- r$best_alpha
    }

    cell_pred_complete <- cell_predictions[complete.cases(cell_predictions), ]
    sample_predictions <- aggregate(predicted ~ sample + actual, data = cell_pred_complete, FUN = mean)

    # 性能（仅原始预测）
    rmse_s <- sqrt(mean((sample_predictions$actual - sample_predictions$predicted)^2))
    mae_s  <- mean(abs(sample_predictions$actual - sample_predictions$predicted))
    pcc_s  <- cor(sample_predictions$actual, sample_predictions$predicted)
    r2_s   <- 1 - sum((sample_predictions$actual - sample_predictions$predicted)^2) /
                  sum((sample_predictions$actual - mean(sample_predictions$actual))^2)
    rmse_c <- sqrt(mean((cell_pred_complete$actual - cell_pred_complete$predicted)^2))
    mae_c  <- mean(abs(cell_pred_complete$actual - cell_pred_complete$predicted))
    pcc_c  <- cor(cell_pred_complete$actual, cell_pred_complete$predicted)
    r2_c   <- 1 - sum((cell_pred_complete$actual - cell_pred_complete$predicted)^2) /
                  sum((cell_pred_complete$actual - mean(cell_pred_complete$actual))^2)
    message(sprintf("  样本级: RMSE=%.3f MAE=%.3f PCC=%.3f R2=%.3f", rmse_s, mae_s, pcc_s, r2_s))
    message(sprintf("  细胞级: RMSE=%.3f MAE=%.3f PCC=%.3f R2=%.3f", rmse_c, mae_c, pcc_c, r2_c))

    all_selected <- unlist(selected_genes_list)
    gene_freq <- sort(table(all_selected), decreasing = TRUE)

    list(
        cell_predictions = cell_pred_complete,
        sample_predictions = sample_predictions,
        selected_genes = selected_genes_list,
        gene_selection_freq = gene_freq,
        best_alphas = best_alphas,
        metrics = data.frame(
            level = c("sample", "cell"),
            RMSE = c(rmse_s, rmse_c), MAE = c(mae_s, mae_c),
            PCC = c(pcc_s, pcc_c), R2 = c(r2_s, r2_c)),
        parameters = list(cor_pval_cutoff = cor_pval_cutoff, cor_abs_cutoff = cor_abs_cutoff, cor_top_n = cor_top_n,
                          min_mean_expr = min_mean_expr, alphas = alphas, nfolds_inner = nfolds_inner,
                          standardize = standardize, calibration = "none",
                          exclude_ages = exclude_ages, exclude_samples = exclude_samples,
                          n_cores = n_cores)
    )
}


######################################################################
######## 主循环：对每个 celltype 跑 Nested LOOCV 并保存 ########
######################################################################
all_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod", "Microglia", "RPE")
available_celltypes <- intersect(all_celltypes, names(Mouse_seurat_knn_pseudobulk))
cat(sprintf("将运行: %s\n", paste(available_celltypes, collapse = ", ")))

all_results <- list()
for (celltype in available_celltypes) {
    message("\n========================================")
    message(sprintf(">>> %s <<<", celltype))
    message("========================================")
    seurat_obj <- Mouse_seurat_knn_pseudobulk[[celltype]]
    n_samples <- length(unique(seurat_obj$sample))
    message(sprintf("  细胞数=%d, 样本数=%d, 基因数=%d", ncol(seurat_obj), n_samples, nrow(seurat_obj)))
    if (n_samples < 3) { message(sprintf("  跳过 %s: 样本数不足", celltype)); next }

    tryCatch({
        results <- Nested_FeatureSelection_LOOCV(
            seurat_obj = seurat_obj,
            whitelist_genes = whitelist_genes,
            ambient_genes = M_ambient_genes_Top50,
            cor_pval_cutoff = 0.05,
            cor_abs_cutoff = 0.2,
            cor_top_n = 750,
            min_mean_expr = 0.05,
            alphas = seq(0.1, 0.9, by = 0.2),
            nfolds_inner = 5,
            standardize = TRUE,
            exclude_ages = NULL,
            exclude_samples = NULL,
            n_cores = 20
        )
        output_dir <- file.path(getwd(), celltype)
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        saveRDS(results, file.path(output_dir, sprintf("%s_aging_clock_results_2026.rds", celltype)))
        write.csv(results$cell_predictions, file.path(output_dir, sprintf("%s_cell_predictions_2026.csv", celltype)), row.names = FALSE)
        write.csv(results$sample_predictions, file.path(output_dir, sprintf("%s_sample_predictions_2026.csv", celltype)), row.names = FALSE)
        message("\n性能指标:"); print(results$metrics)
        all_results[[celltype]] <- results
        message(sprintf("%s 完成\n", celltype))
    }, error = function(e) {
        message(sprintf("%s 出错: %s\n", celltype, e$message))
    })
}
message("\n=== 所有 celltype Nested LOOCV（保留标准化/无校准）完成 ===")





#####
##### 按照之前的 Mouse 逻辑，用 所有 sample 训练一遍 ，生成最后的 aging clocks ######
##### 每个 cell type 生成 一个 R script，在 HPC 运行  ####
##### 除了 RPE 和 Microglia ####
#####



######################################################################
######## 检查 Mouse aging clocks 是否跑完（各 celltype 文件夹里有没有输出）########
######################################################################
base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod", "Microglia", "RPE")

# 每个 celltype 期望的输出文件（LOOCV + Final）
expected_files <- function(ct) {
    c(
        LOOCV_rds      = sprintf("%s_aging_clock_results_2026.rds", ct),
        LOOCV_cellpred = sprintf("%s_cell_predictions_2026.csv", ct),
        LOOCV_samppred = sprintf("%s_sample_predictions_2026.csv", ct),
        FINAL_rds      = sprintf("%s_final_aging_clock_model_2026.rds", ct),
        FINAL_coefs    = sprintf("%s_final_model_coefs_2026.csv", ct)
    )
}

check_rows <- list()
cat("\n==================== Mouse Aging Clocks 输出检查 ====================\n")
for (ct in celltypes) {
    ct_dir <- file.path(base_dir, ct)
    exp <- expected_files(ct)
    cat(sprintf("\n--- %s (%s) ---\n", ct, if (dir.exists(ct_dir)) ct_dir else paste0(ct_dir, "  [目录不存在]")))
    for (tag in names(exp)) {
        fp <- file.path(ct_dir, exp[[tag]])
        ok <- file.exists(fp)
        sz <- if (ok) file.info(fp)$size else NA
        mt <- if (ok) format(file.info(fp)$mtime, "%Y-%m-%d %H:%M") else NA
        cat(sprintf("  [%s] %-12s %s%s\n",
                    if (ok) "OK " else "MISS", tag, exp[[tag]],
                    if (ok) sprintf("  (%.1f KB, %s)", sz/1024, mt) else ""))
        check_rows[[length(check_rows) + 1]] <- data.frame(
            celltype = ct, item = tag, file = exp[[tag]], exists = ok,
            size_KB = if (ok) round(sz/1024, 1) else NA, mtime = mt,
            stringsAsFactors = FALSE)
    }
}

# 汇总表
check_df <- do.call(rbind, check_rows)

cat("\n==================== 汇总 ====================\n")
# 按 celltype 汇总 LOOCV / Final 是否齐全
for (ct in celltypes) {
    sub <- check_df[check_df$celltype == ct, ]
    loocv_ok <- all(sub$exists[grepl("^LOOCV", sub$item)])
    final_ok <- all(sub$exists[grepl("^FINAL", sub$item)])
    cat(sprintf("  %-10s LOOCV: %-6s  Final: %-6s\n", ct,
                if (loocv_ok) "完成" else "缺失",
                if (final_ok) "完成" else "缺失"))
}

n_total <- nrow(check_df)
n_ok <- sum(check_df$exists)
cat(sprintf("\n文件齐全度: %d / %d 存在\n", n_ok, n_total))
missing <- check_df[!check_df$exists, c("celltype", "file")]
if (nrow(missing) > 0) {
    cat("\n缺失文件:\n")
    print(missing, row.names = FALSE)
} else {
    cat("\n全部输出文件都在！\n")
}
# 可选：写出检查结果
# write.csv(check_df, file.path(base_dir, "AgingClock_output_check_2026.csv"), row.names = FALSE)


######################################################################
######## 只看 Final 模型有没有保存（major cell type，去掉 RPE/Microglia）########
######################################################################
base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
major_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")

cat("\n==================== Mouse Final 模型检查（major cell type）====================\n")
n_done <- 0
for (ct in major_celltypes) {
    fp <- file.path(base_dir, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    ok <- file.exists(fp)
    if (ok) {
        info <- file.info(fp)
        cat(sprintf("  [OK ] %-6s %.1f KB  %s\n", ct, info$size/1024, format(info$mtime, "%Y-%m-%d %H:%M")))
        n_done <- n_done + 1
    } else {
        cat(sprintf("  [MISS] %-6s 未保存 -> %s\n", ct, fp))
    }
}
cat(sprintf("\nFinal 模型完成: %d / %d\n", n_done, length(major_celltypes)))


######################################################################
######## 只看 LOOCV 结果有没有保存（major cell type，去掉 RPE/Microglia）########
######################################################################
base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
major_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")

cat("\n==================== Mouse LOOCV 结果检查（major cell type）====================\n")
n_done <- 0
for (ct in major_celltypes) {
    fp <- file.path(base_dir, ct, sprintf("%s_aging_clock_results_2026.rds", ct))
    ok <- file.exists(fp)
    if (ok) {
        info <- file.info(fp)
        cat(sprintf("  [OK ] %-6s %.1f KB  %s\n", ct, info$size/1024, format(info$mtime, "%Y-%m-%d %H:%M")))
        n_done <- n_done + 1
    } else {
        cat(sprintf("  [MISS] %-6s 未保存 -> %s\n", ct, fp))
    }
}
cat(sprintf("\nLOOCV 结果完成: %d / %d\n", n_done, length(major_celltypes)))


##########  然后我们 copy metformin Dox to the project ####
##########


######################################################################
######## 读入 LOOCV 结果，画 violin plot（每个 celltype 一个分面）########
######## x = 真实年龄, y = 预测年龄, 按年龄 viridis 配色, 内嵌白色箱线图 ########
######################################################################
library(ggplot2)
library(viridis)

base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
violin_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")   # 单行顺序

# 读入每个 celltype 的 LOOCV cell 级预测
pred_all <- list()
for (ct in violin_celltypes) {
    fp <- file.path(base_dir, ct, sprintf("%s_aging_clock_results_2026.rds", ct))
    if (!file.exists(fp)) { message(ct, ": LOOCV 结果不存在，跳过 -> ", fp); next }
    res <- readRDS(fp)
    df <- res$cell_predictions          # cell_id, sample, actual, predicted
    df$celltype <- ct
    pred_all[[ct]] <- df[, c("celltype", "actual", "predicted")]
}
pred_df <- do.call(rbind, pred_all)

# celltype 与 age 顺序
pred_df$celltype  <- factor(pred_df$celltype, levels = violin_celltypes)
pred_df$actual_f  <- factor(pred_df$actual, levels = sort(unique(pred_df$actual)))

p_violin <- ggplot(pred_df, aes(x = actual_f, y = predicted, fill = actual_f)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.9) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", lwd = 0.3) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "D") +
    facet_wrap(~ celltype, nrow = 1, scales = "free_x") +
    coord_cartesian(ylim = c(-10, 125)) +
    labs(x = "Ages", y = "Predicted age") +
    theme_classic(base_size = 12) +
    theme(
        strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text       = element_text(face = "bold", size = 15),
        axis.text.x      = element_text(size = 16, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y      = element_text(size = 13, color = "black"),
        axis.title       = element_text(size = 18, face = "bold"),
        legend.position  = "none",
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

ggsave("Mouse_LOOCV_violin_by_celltype_2026.png", p_violin, width = 18, height = 3.2, dpi = 300)
message("Saved: Mouse_LOOCV_violin_by_celltype_2026.png")


######################################################################
######## 按画图顺序，逐个打印每个 celltype 的 RMSE / MAE / R2 / Pearson ########
######## （cell 级，和 violin 一致；基于 LOOCV cell_predictions 现算）      ########
######################################################################
metric_tbl <- data.frame()
for (ct in violin_celltypes) {
    d <- pred_df[pred_df$celltype == ct & is.finite(pred_df$predicted), ]
    if (nrow(d) == 0) { message(ct, ": 无预测数据，跳过"); next }
    a <- d$actual; p <- d$predicted
    rmse    <- sqrt(mean((a - p)^2))
    mae     <- mean(abs(a - p))
    r2      <- 1 - sum((a - p)^2) / sum((a - mean(a))^2)
    pearson <- cor(a, p)
    cat(sprintf("%-5s  RMSE=%.3f  MAE=%.3f  R2=%.3f  Pearson=%.3f\n",
                ct, rmse, mae, r2, pearson))
    metric_tbl <- rbind(metric_tbl,
        data.frame(celltype = ct, RMSE = rmse, MAE = mae, R2 = r2, Pearson = pearson))
}
write.csv(metric_tbl, "Mouse_LOOCV_metrics_by_celltype_2026.csv", row.names = FALSE)
message("Saved: Mouse_LOOCV_metrics_by_celltype_2026.csv")



