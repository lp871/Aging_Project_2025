########
########

######## for Mouse ########
########

##### Calibration removed ######
# Calibration: aggregate training-set predictions by sample and fit lm(actual ~ predicted) to obtain intercept/slope

######## 
######## Mouse preprocessing ########
########

conda activate ArchR
R

#########
library(Seurat)
library(FNN)
library(Matrix)

setwd("/projects/hmz-aging/Mouse_New_Datasets")
Mouse_seurat_clcl <- readRDS("Mouse_seurat_merge_clcl_2025")

# Remove the _T1 and _T2 suffixes from sample names
Mouse_seurat_clcl$sample <- gsub("_T[12]$", "", Mouse_seurat_clcl$sample)

######### Print all samples / cell types #########
cat("\n========== All samples (after removing _T1/_T2 suffixes)==========\n")
print(sort(unique(Mouse_seurat_clcl$sample)))
cat("\nCells per sample:\n")
print(table(Mouse_seurat_clcl$sample))
cat("\nCells per celltype:\n")
print(table(Mouse_seurat_clcl$celltype))

######### Split by celltype #########
Mouse_seurat_clcl_list <- SplitObject(Mouse_seurat_clcl, split.by = "celltype")

######### Print samples within each celltype #########
cat("\n========== Samples per celltype ==========\n")
for (ct in names(Mouse_seurat_clcl_list)) {
    cat(sprintf("\n--- %s ---\n", ct))
    print(table(Mouse_seurat_clcl_list[[ct]]$sample))
}


######################################################################
######## Helper function 1: normalization + PCA ########
######################################################################
NormalizeAndPCA <- function(
  seu,
  norm.method = "LogNormalize",
  scale.factor = 1e4,
  nfeatures = 2000,
  npcs = 30
) {
  # 1. Normalization
  seu <- NormalizeData(
    object = seu,
    normalization.method = norm.method,
    scale.factor = scale.factor
  )

  # 2. Identify highly variable genes
  seu <- FindVariableFeatures(
    object = seu,
    selection.method = "vst",
    nfeatures = nfeatures
  )

  # 3. Scale data (all genes)
  seu <- ScaleData(
    object = seu,
    features = rownames(seu)
  )

  # 4. PCA dimensionality reduction (npcs adjusted dynamically)
  # SVD requirement: npcs must be **strictly less than** min(n_genes, n_cells)
  n_cells <- ncol(seu)
  n_genes <- nrow(seu)
  max_npcs <- min(n_genes, n_cells) - 1

  if (npcs > max_npcs) {
    npcs_adjusted <- max(1, max_npcs)  # at least 1 PC
    message(sprintf("    Adjusting npcs: %d → %d (n_cells=%d, n_genes=%d)",
                    npcs, npcs_adjusted, n_cells, n_genes))
    npcs <- npcs_adjusted
  }

  if (npcs < 1) {
    stop(sprintf("Too few cells (n_cells=%d); cannot run PCA", n_cells))
  }

  seu <- RunPCA(
    object = seu,
    features = VariableFeatures(seu),
    npcs = npcs
  )

  return(seu)
}


######################################################################
######## Helper function 2: KNN smoothing -> pseudo-cell ########
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

  # 2. kNN (excluding self)
  knn <- get.knn(emb, k = k)
  nn_idx <- knn$nn.index
  ncells <- nrow(nn_idx)

  # 3. Build sparse adjacency matrix A: for each target cell j, set j itself and its k neighbors i to 1
  rows <- as.integer(c(t(nn_idx), seq_len(ncells)))
  cols <- rep(seq_len(ncells), each = k + 1)
  A <- sparseMatrix(i = rows, j = cols, x = 1,
                    dims = c(ncells, ncells))

  # 4. Raw counts (genes × cells)
  counts_mat <- GetAssayData(seu, assay = assay, layer = "counts")

  # 5. A single matrix multiplication yields pseudo-bulk counts
  pseudo_counts <- counts_mat %*% A

  # 6. New Seurat object
  new_seu <- CreateSeuratObject(
    counts = pseudo_counts,
    meta.data = seu@meta.data,
    assay = assay
  )
  return(new_seu)
}


######################################################################
######## Main loop: smooth each celltype into pseudo-cells by sample ########
######################################################################
Mouse_seurat_knn_pseudobulk <- list()

for(i in 1:length(Mouse_seurat_clcl_list)){
    print(i)
    print(names(Mouse_seurat_clcl_list)[i])
    #####
    tmp = Mouse_seurat_clcl_list[[i]]
    #####
    tmp_sp = SplitObject(tmp, split.by = "sample")

    ##### Print samples split from this celltype #####
    cat(sprintf("  [%s] samples: %s\n",
                names(Mouse_seurat_clcl_list)[i],
                paste(names(tmp_sp), collapse = ", ")))

    ##### Filter out samples with <2 cells #######
    sample_ncells <- sapply(tmp_sp, ncol)
    removed_samples <- names(sample_ncells)[sample_ncells < 2]
    if(length(removed_samples) > 0) {
        message(sprintf("  %s: removed samples with <2 cells: %s", names(Mouse_seurat_clcl_list[i]), paste(removed_samples, collapse=", ")))
        tmp_sp <- tmp_sp[sample_ncells >= 2]
    }
    if(length(tmp_sp) == 0) {
        message(sprintf("  %s: all samples filtered out, skipping", names(Mouse_seurat_clcl_list[i])))
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
        k = 2  # RPE and Microglia: K = 2
        npcs = 5
    }
    print(paste("Cell type:", names(Mouse_seurat_clcl_list[i])))
    print(paste("K =", k))
    print(paste("npcs (initial) =", npcs))

    ##### Dynamically adjust npcs and K: check cell count per sample #######
    sample_ncells_check <- sapply(tmp_sp, ncol)
    min_ncells <- min(sample_ncells_check)
    print(paste("Min cells per sample:", min_ncells))

    # Special case: if too few cells (< 5), skip KNN smoothing
    if (min_ncells < 5) {
        message(sprintf("  Warning: minimum cells per sample (%d) < 5; skipping KNN smoothing and merging directly", min_ncells))
        tmp_sp_res_pseudo_merge = merge(tmp_sp[[1]], tmp_sp[-1], add.cell.ids = names(tmp_sp))
        Mouse_seurat_knn_pseudobulk <- c(Mouse_seurat_knn_pseudobulk, list(tmp_sp_res_pseudo_merge))
        next
    }

    # npcs must be strictly less than the minimum cells per sample
    if (npcs >= min_ncells) {
        npcs_adjusted <- max(1, min_ncells - 1)
        message(sprintf("  Warning: npcs (%d) >= minimum cells per sample (%d); adjusting to %d",
                        npcs, min_ncells, npcs_adjusted))
        npcs <- npcs_adjusted
    }

    # K must be less than the cell count (KNN needs neighbors)
    if (k >= min_ncells) {
        k_adjusted <- max(1, min_ncells - 1)
        message(sprintf("  Warning: K (%d) >= minimum cells per sample (%d); adjusting to %d",
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

##### Print samples per celltype after pseudobulking #####
cat("\n========== Samples per celltype after KNN pseudobulk ==========\n")
for (ct in names(Mouse_seurat_knn_pseudobulk)) {
    cat(sprintf("\n--- %s ---\n", ct))
    print(table(Mouse_seurat_knn_pseudobulk[[ct]]$sample))
}

####### counts matrix ########
mat <- LayerData(Mouse_seurat_knn_pseudobulk$MG, layer = "counts")
range(mat)                     # minimum is 0
all(mat@x == round(mat@x))     # TRUE means integers (summed raw counts)
head(Matrix::colSums(mat))     # total counts per pseudo-cell



#######
####### Mouse: next, nested leave-one-out ####
####### Without the final normalization ###
#######   Keep per-fold Z-score standardization (glmnet itself uses standardize=FALSE)
#######   Remove the final sample-level calibration (raw predictions only; no predicted_calibrated)

library(Seurat)
library(glmnet)
library(Matrix)

setwd("/projects/hmz-aging/Mouse_New_Datasets")

######### Load data and ambient / whitelist genes #########
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

######### Add age column to each celltype (parsed from sample: Wk<age>_Rep<n>; 106 -> 108)#########
for (ct in names(Mouse_seurat_knn_pseudobulk)) {
    sample_names <- Mouse_seurat_knn_pseudobulk[[ct]]$sample
    ages <- as.numeric(gsub("^Wk([0-9]+)_.*", "\\1", sample_names))
    ages[ages == 106] <- 108
    Mouse_seurat_knn_pseudobulk[[ct]]$age <- ages
    cat(sprintf("  %s: age range = %s weeks\n", ct, paste(sort(unique(ages)), collapse = ", ")))
}

######### Reuse helper functions (pseudobulk / correlation, etc.) #########
source("Mouse_AgingClock_Functions.R")


######################################################################
######## Streamlined Nested LOOCV: Z-score standardization kept, no calibration ########
######## Overrides the same-named function in Mouse_AgingClock_Functions.R ########
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
    message("Nested Feature Selection + LOOCV (Z-score standardization kept / no calibration)")
    message(sprintf("%s", if (use_parallel) sprintf("Parallel: %d cores", n_cores) else "Sequential"))
    message("========================================")

    # 1. Prepare data
    if (length(Layers(seurat_obj)) > 1) seurat_obj <- JoinLayers(seurat_obj)
    counts <- LayerData(seurat_obj, layer = "counts")
    cell_samples <- seurat_obj$sample
    cell_ages <- seurat_obj$age
    if (is.null(cell_samples) || is.null(cell_ages)) stop("Seurat object is missing 'sample' or 'age' column")

    # Filter out specified ages
    if (!is.null(exclude_ages)) {
        keep <- !cell_ages %in% exclude_ages
        if (sum(!keep) > 0) {
            counts <- counts[, keep, drop = FALSE]; cell_samples <- cell_samples[keep]; cell_ages <- cell_ages[keep]
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
            message(sprintf("  Excluding ages %s: removed %d cells", paste(exclude_ages, collapse=","), sum(!keep)))
        }
    }
    # Filter out specified samples
    if (!is.null(exclude_samples)) {
        keep <- !cell_samples %in% exclude_samples
        if (sum(!keep) > 0) {
            counts <- counts[, keep, drop = FALSE]; cell_samples <- cell_samples[keep]; cell_ages <- cell_ages[keep]
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
            message(sprintf("  Excluding samples %s: removed %d cells", paste(exclude_samples, collapse=","), sum(!keep)))
        }
    }
    # Remove ambient genes
    if (!is.null(ambient_genes)) {
        ambient_genes <- intersect(ambient_genes, rownames(counts))
        if (length(ambient_genes) > 0) counts <- counts[!rownames(counts) %in% ambient_genes, ]
    }
    # Single-cell depth normalization + log2 (for modeling; note this is not the "standardization" being removed)
    depth <- Matrix::colSums(counts)
    counts_norm <- sweep(counts, 2, depth / 1e4, FUN = "/")
    counts_log2 <- log2(1 + counts_norm)
    message(sprintf("  Single cells: %d, genes: %d", ncol(counts_log2), nrow(counts_log2)))

    # 2. Pseudobulk (for feature selection)
    pseudobulk_all <- Calculate_Pseudobulk_Normalized(seurat_obj, ambient_genes)
    unique_samples <- unique(cell_samples)
    sample_ages <- sapply(unique_samples, function(s) unique(cell_ages[cell_samples == s])[1])
    n_samples <- length(unique_samples)
    message(sprintf("  Samples: %d, ages: %s", n_samples, paste(sort(unique(sample_ages)), collapse=", ")))

    # 3. LOOCV
    process_one_loocv <- function(i) {
        test_sample <- unique_samples[i]
        train_samples <- unique_samples[-i]
        train_pseudobulk <- pseudobulk_all[, train_samples, drop = FALSE]

        # 3.1 Feature selection (training set only)
        cor_results <- suppressMessages(Calculate_Gene_Age_Correlation(train_pseudobulk, method = "spearman"))
        cor_results_sorted <- cor_results[order(abs(cor_results$correlation), decreasing = TRUE), ]
        # Strongly correlated genes: |cor| > cor_abs_cutoff and p < cor_pval_cutoff; if fewer than cor_top_n, fall back to the top cor_top_n by |cor|
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
        if (length(candidate_genes) == 0) stop("No candidate genes left after filtering!")

        # 3.2 Training/test single-cell matrices (candidate genes only, log2 normalized, no z-score)
        train_cell_idx <- which(cell_samples %in% train_samples)
        test_cell_idx <- which(cell_samples == test_sample)
        if (length(train_cell_idx) == 0) stop("Training set has no cells")
        if (length(test_cell_idx) == 0) stop(sprintf("Test sample %s has no cells", test_sample))
        candidate_genes <- intersect(candidate_genes, rownames(counts_log2))
        X_train <- t(as.matrix(counts_log2[candidate_genes, train_cell_idx, drop = FALSE]))
        y_train <- cell_ages[train_cell_idx]
        train_samples_vec <- cell_samples[train_cell_idx]
        ## Z-score standardization (training-set mean/sd, applied unchanged to the test set to prevent leakage)
        if (standardize) {
            train_means <- colMeans(X_train)
            train_sds <- apply(X_train, 2, sd)
            train_sds[train_sds == 0] <- 1
            X_train <- scale(X_train, center = train_means, scale = train_sds)
        }

        # 3.3 IPF weights: balance sample and age
        weights <- rep(1, length(train_cell_idx))
        for (iter in 1:10) {
            sample_weight_sum <- tapply(weights, train_samples_vec, sum)
            weights <- weights * (mean(sample_weight_sum) / sample_weight_sum[train_samples_vec])
            age_weight_sum <- tapply(weights, y_train, sum)
            weights <- weights * (mean(age_weight_sum) / age_weight_sum[as.character(y_train)])
        }
        weights <- weights * length(weights) / sum(weights)

        # 3.4 Inner CV to tune alpha (glmnet standardize=FALSE, i.e. truly no standardization)
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

        # 3.5 Train final model + predict test sample (raw predictions, no calibration)
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

    # Combine results
    valid <- sapply(loocv_results, function(x) is.list(x) && !is.null(x$test_sample))
    if (!all(valid)) stop(sprintf("LOOCV Round %s invalid", paste(which(!valid), collapse=", ")))

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

    # Performance (raw predictions only)
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
    message(sprintf("  Sample-level: RMSE=%.3f MAE=%.3f PCC=%.3f R2=%.3f", rmse_s, mae_s, pcc_s, r2_s))
    message(sprintf("  Cell-level: RMSE=%.3f MAE=%.3f PCC=%.3f R2=%.3f", rmse_c, mae_c, pcc_c, r2_c))

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
######## Main loop: run Nested LOOCV for each celltype and save ########
######################################################################
all_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod", "Microglia", "RPE")
available_celltypes <- intersect(all_celltypes, names(Mouse_seurat_knn_pseudobulk))
cat(sprintf("Will run: %s\n", paste(available_celltypes, collapse = ", ")))

all_results <- list()
for (celltype in available_celltypes) {
    message("\n========================================")
    message(sprintf(">>> %s <<<", celltype))
    message("========================================")
    seurat_obj <- Mouse_seurat_knn_pseudobulk[[celltype]]
    n_samples <- length(unique(seurat_obj$sample))
    message(sprintf("  cells=%d, samples=%d, genes=%d", ncol(seurat_obj), n_samples, nrow(seurat_obj)))
    if (n_samples < 3) { message(sprintf("  Skipping %s: not enough samples", celltype)); next }

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
        message("\nPerformance metrics:"); print(results$metrics)
        all_results[[celltype]] <- results
        message(sprintf("%s done\n", celltype))
    }, error = function(e) {
        message(sprintf("%s failed: %s\n", celltype, e$message))
    })
}
message("\n=== Nested LOOCV (standardization kept / no calibration) done for all celltypes ===")





#####
##### Following the previous Mouse logic, train once on all samples to generate the final aging clocks ######
##### Generate one R script per cell type to run on HPC  ####
##### Except RPE and Microglia ####
#####



######################################################################
######## Check whether Mouse aging clocks have finished (outputs present in each celltype folder)########
######################################################################
base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod", "Microglia", "RPE")

# Expected output files per celltype (LOOCV + Final)
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
cat("\n==================== Mouse Aging Clocks output check ====================\n")
for (ct in celltypes) {
    ct_dir <- file.path(base_dir, ct)
    exp <- expected_files(ct)
    cat(sprintf("\n--- %s (%s) ---\n", ct, if (dir.exists(ct_dir)) ct_dir else paste0(ct_dir, "  [directory not found]")))
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

# Summary table
check_df <- do.call(rbind, check_rows)

cat("\n==================== Summary ====================\n")
# Summarize LOOCV / Final completeness by celltype
for (ct in celltypes) {
    sub <- check_df[check_df$celltype == ct, ]
    loocv_ok <- all(sub$exists[grepl("^LOOCV", sub$item)])
    final_ok <- all(sub$exists[grepl("^FINAL", sub$item)])
    cat(sprintf("  %-10s LOOCV: %-6s  Final: %-6s\n", ct,
                if (loocv_ok) "done" else "missing",
                if (final_ok) "done" else "missing"))
}

n_total <- nrow(check_df)
n_ok <- sum(check_df$exists)
cat(sprintf("\nFile completeness: %d / %d present\n", n_ok, n_total))
missing <- check_df[!check_df$exists, c("celltype", "file")]
if (nrow(missing) > 0) {
    cat("\nMissing files:\n")
    print(missing, row.names = FALSE)
} else {
    cat("\nAll output files are present!\n")
}
# Optional: write out check results
# write.csv(check_df, file.path(base_dir, "AgingClock_output_check_2026.csv"), row.names = FALSE)


######################################################################
######## Check only whether Final models are saved (major cell types, excluding RPE/Microglia)########
######################################################################
base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
major_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")

cat("\n==================== Mouse Final model check (major cell types)====================\n")
n_done <- 0
for (ct in major_celltypes) {
    fp <- file.path(base_dir, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    ok <- file.exists(fp)
    if (ok) {
        info <- file.info(fp)
        cat(sprintf("  [OK ] %-6s %.1f KB  %s\n", ct, info$size/1024, format(info$mtime, "%Y-%m-%d %H:%M")))
        n_done <- n_done + 1
    } else {
        cat(sprintf("  [MISS] %-6s not saved -> %s\n", ct, fp))
    }
}
cat(sprintf("\nFinal models done: %d / %d\n", n_done, length(major_celltypes)))


######################################################################
######## Check only whether LOOCV results are saved (major cell types, excluding RPE/Microglia)########
######################################################################
base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
major_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")

cat("\n==================== Mouse LOOCV results check (major cell types)====================\n")
n_done <- 0
for (ct in major_celltypes) {
    fp <- file.path(base_dir, ct, sprintf("%s_aging_clock_results_2026.rds", ct))
    ok <- file.exists(fp)
    if (ok) {
        info <- file.info(fp)
        cat(sprintf("  [OK ] %-6s %.1f KB  %s\n", ct, info$size/1024, format(info$mtime, "%Y-%m-%d %H:%M")))
        n_done <- n_done + 1
    } else {
        cat(sprintf("  [MISS] %-6s not saved -> %s\n", ct, fp))
    }
}
cat(sprintf("\nLOOCV results done: %d / %d\n", n_done, length(major_celltypes)))


##########  Then copy metformin Dox to the project ####
##########


######################################################################
######## Load LOOCV results and draw violin plots (one facet per celltype)########
######## x = actual age, y = predicted age, viridis colors by age, embedded white boxplots ########
######################################################################
library(ggplot2)
library(viridis)

base_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
violin_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")   # single-row order

# Load cell-level LOOCV predictions for each celltype
pred_all <- list()
for (ct in violin_celltypes) {
    fp <- file.path(base_dir, ct, sprintf("%s_aging_clock_results_2026.rds", ct))
    if (!file.exists(fp)) { message(ct, ": LOOCV results not found, skipping -> ", fp); next }
    res <- readRDS(fp)
    df <- res$cell_predictions          # cell_id, sample, actual, predicted
    df$celltype <- ct
    pred_all[[ct]] <- df[, c("celltype", "actual", "predicted")]
}
pred_df <- do.call(rbind, pred_all)

# celltype and age order
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
######## In plotting order, print RMSE / MAE / R2 / Pearson for each celltype ########
######## (cell-level, consistent with violin; computed from LOOCV cell_predictions) ########
######################################################################
metric_tbl <- data.frame()
for (ct in violin_celltypes) {
    d <- pred_df[pred_df$celltype == ct & is.finite(pred_df$predicted), ]
    if (nrow(d) == 0) { message(ct, ": no prediction data, skipping"); next }
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



