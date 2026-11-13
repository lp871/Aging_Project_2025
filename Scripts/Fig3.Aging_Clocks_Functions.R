########################################################################
######## Aging Clocks Functions for Human Data
######## All R functions used for Elastic Net LODO
########################################################################

# Load required packages
library(Matrix)
library(Seurat)
library(glmnet)
library(parallel)
library(ggplot2)
library(viridis)

# 1. Read data from an mtx folder and create a Seurat object
load_mtx_to_seurat <- function(mtx_dir) {
  mat <- readMM(file.path(mtx_dir, "matrix.mtx"))
  mat <- as(mat, "CsparseMatrix")
  genes <- read.table(file.path(mtx_dir, "features.tsv"), header = FALSE, stringsAsFactors = FALSE)$V1
  barcodes <- read.table(file.path(mtx_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)$V1
  rownames(mat) <- genes
  colnames(mat) <- barcodes
  meta <- read.table(file.path(mtx_dir, "metadata.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(meta) <- meta[[1]]
  CreateSeuratObject(counts = mat, meta.data = meta)
}

# 2. Batch-process a single celltype (read mtx -> Seurat -> save rds)
process_one_celltype <- function(ct, mtx_base, seurat_out_dir) {
  mtx_dir <- file.path(mtx_base, paste0(ct, "_clean4_smoothed"))
  if (!dir.exists(mtx_dir)) { message(ct, ": mtx dir not found"); return(invisible(NULL)) }
  message("Reading: ", ct)
  seu <- load_mtx_to_seurat(mtx_dir)
  message("  ", ncol(seu), " cells, ", nrow(seu), " genes")
  rds_out <- file.path(seurat_out_dir, paste0(ct, "_clean4_smoothed.rds"))
  saveRDS(seu, rds_out)
  message("  Saved: ", rds_out)
  invisible(NULL)
}

# 3. Add Age_numeric and age_interval columns (10-year and 30-year)
seurat_add_age_numeric_and_interval <- function(seu, age_col = "age") {
  meta <- seu@meta.data
  if (!age_col %in% colnames(meta)) stop(paste("age column not found:", age_col))
  message("First 5 values of original age column: ", paste(head(meta[[age_col]], 5), collapse = ", "))
  age_str <- as.character(meta[[age_col]])
  age_numeric <- as.numeric(gsub("[^0-9.]", "", age_str))
  meta$Age_numeric <- age_numeric
  message("First 5 values of converted Age_numeric: ", paste(head(age_numeric, 5), collapse = ", "))
  age_interval_10yr <- cut(age_numeric, breaks = seq(0, 90, by = 10), labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90"), right = FALSE, include.lowest = TRUE)
  meta$age_interval_10yr <- as.character(age_interval_10yr)
  age_interval_30yr <- cut(age_numeric, breaks = c(0, 30, 60, 90), labels = c("0-30", "30-60", "60-90"), right = FALSE, include.lowest = TRUE)
  meta$age_interval_30yr <- as.character(age_interval_30yr)
  seu@meta.data <- meta
  seu
}

# 4. Add cell-level balance weights (IPF iterative balancing)
seurat_add_balance_weights <- function(x, donor_col = "donor_id", sex_col = "sex", age_interval_col = "age_interval_10yr", target_per_donor = 1, target_per_interval = 1, max_iter = 100, tol = 1e-6) {
  is_seurat <- inherits(x, "Seurat")
  meta <- if (is_seurat) x@meta.data else x
  if (!all(c(donor_col, sex_col, age_interval_col) %in% names(meta))) stop("Required columns not found")
  donor <- as.character(meta[[donor_col]])
  sex <- as.character(meta[[sex_col]])
  aint <- as.character(meta[[age_interval_col]])
  keep <- !is.na(donor) & donor != "" & !is.na(sex) & sex != "" & !is.na(aint) & aint != ""
  message("Adding balance weights (IPF): balancing age_interval (", age_interval_col, "), sex, donor")
  message("Valid cells: ", sum(keep), "/", length(donor))
  meta$weight <- 1.0
  meta$weight[!keep] <- 0
  unique_donors <- unique(donor[keep])
  unique_aint <- unique(aint[keep])
  n_donors <- length(unique_donors)
  n_aint <- length(unique_aint)
  message("unique donors: ", n_donors, ", unique age_interval: ", n_aint)
  for (iter in 1:max_iter) {
    s_donor <- tapply(meta$weight, donor, sum, na.rm = TRUE)
    for (d in names(s_donor)) {
      if (s_donor[d] > 0) meta$weight[donor == d] <- meta$weight[donor == d] * (target_per_donor / s_donor[d])
    }
    s_int <- tapply(meta$weight, aint, sum, na.rm = TRUE)
    for (iv in names(s_int)) {
      idx <- aint == iv
      if (s_int[iv] == 0) next
      meta$weight[idx] <- meta$weight[idx] * (target_per_interval / s_int[iv])
    }
    for (iv in unique_aint) {
      idx <- aint == iv
      s_sex <- tapply(meta$weight[idx], sex[idx], sum, na.rm = TRUE)
      half <- target_per_interval / 2
      for (s in names(s_sex)) {
        if (s_sex[s] > 0) meta$weight[idx & sex == s] <- meta$weight[idx & sex == s] * (half / s_sex[s])
      }
    }
    s_donor_new <- tapply(meta$weight, donor, sum, na.rm = TRUE)
    s_int_new <- tapply(meta$weight, aint, sum, na.rm = TRUE)
    s_int_new <- s_int_new[!is.na(names(s_int_new)) & names(s_int_new) != ""]
    if (length(s_int_new) > 0 && max(abs(s_donor_new - target_per_donor)) < tol && max(abs(s_int_new - target_per_interval)) < tol) break
  }
  message("===== balance weights convergence: marginal sums (check) =====")
  msg_donor <- paste(round(head(sort(unique(tapply(meta$weight, donor, sum))), 5), 3), collapse = ", ")
  message("Per-donor weight sum (first 5): ", msg_donor)
  message("Per age_interval weight sum:")
  print(tapply(meta$weight, aint, sum, na.rm = TRUE))
  message("Per (age_interval, sex) weight sum:")
  print(tapply(meta$weight, list(aint, sex), sum, na.rm = TRUE))
  if (is_seurat) {
    x@meta.data <- meta
    return(x)
  }
  meta
}

# 5. Downsample donors in old intervals (60-70, 70-80, 80-90), stratified by sex
seurat_downsample_old_donors_10yr <- function(seu, donor_col = "donor_id", age_interval_col = "age_interval_10yr", sex_col = "sex", seed = 123) {
  meta <- seu@meta.data
  if (!donor_col %in% names(meta) || !age_interval_col %in% names(meta)) stop("columns not found")
  aint <- as.character(meta[[age_interval_col]])
  donor <- as.character(meta[[donor_col]])
  if (sex_col %in% names(meta)) {
    sex <- as.character(meta[[sex_col]])
    keep <- !is.na(aint) & aint != "" & !is.na(donor) & donor != "" & !is.na(sex) & sex != ""
    aint <- aint[keep]
    donor <- donor[keep]
    sex <- sex[keep]
    cells_keep <- rownames(meta)[keep]
    group <- paste(aint, sex, sep = "_")
    uniq_donors_by_group <- lapply(split(donor, group), unique)
    n_donors_by_group <- sapply(uniq_donors_by_group, length)
    message("Number of donors per (10yr interval, sex) (before downsampling):")
    print(n_donors_by_group)
    young_intervals <- c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60")
    old_intervals <- c("60-70", "70-80", "80-90")
    target_by_sex <- list()
    for (s in unique(sex)) {
      young_groups <- paste(young_intervals, s, sep = "_")
      young_counts <- n_donors_by_group[names(n_donors_by_group) %in% young_groups]
      if (length(young_counts) > 0) {
        target_by_sex[[s]] <- round(mean(young_counts))
      }
    }
    message("Mean number of donors per sex at young time points (0-60 years):")
    print(unlist(target_by_sex))
    set.seed(seed)
    selected_donors <- character(0)
    message("\nDownsampling strategy:")
    for (grp_name in names(uniq_donors_by_group)) {
      donors_in_grp <- uniq_donors_by_group[[grp_name]]
      parts <- strsplit(grp_name, "_")[[1]]
      int_name <- parts[1]
      sex_name <- parts[2]
      if (int_name %in% old_intervals && sex_name %in% names(target_by_sex)) {
        target_n <- target_by_sex[[sex_name]]
        n_select <- min(target_n, length(donors_in_grp))
        selected_donors <- c(selected_donors, sample(donors_in_grp, n_select))
        message("  [Downsampled] ", grp_name, ": ", length(donors_in_grp), " -> ", n_select, " donors")
      } else {
        selected_donors <- c(selected_donors, donors_in_grp)
        if (int_name %in% old_intervals) {
          message("  [Kept all] ", grp_name, ": ", length(donors_in_grp), " donors (no target for this sex)")
        } else {
          message("  [Kept all] ", grp_name, ": ", length(donors_in_grp), " donors (young time point)")
        }
      }
    }
  } else {
    warning("sex column not found, downsampling without sex stratification")
    keep <- !is.na(aint) & aint != "" & !is.na(donor) & donor != ""
    aint <- aint[keep]
    donor <- donor[keep]
    cells_keep <- rownames(meta)[keep]
    uniq_donors_by_aint <- lapply(split(donor, aint), unique)
    n_donors <- sapply(uniq_donors_by_aint, length)
    message("Number of donors per 10yr interval (before downsampling):")
    print(n_donors)
    young_intervals <- c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60")
    old_intervals <- c("60-70", "70-80", "80-90")
    young_counts <- n_donors[names(n_donors) %in% young_intervals]
    if (length(young_counts) == 0) {
      warning("No young intervals (0-60 years) found, cannot compute target")
      return(seu)
    }
    target_n <- round(mean(young_counts))
    message("Mean number of donors at young time points (0-60 years): ", target_n)
    set.seed(seed)
    selected_donors <- character(0)
    for (int_name in names(uniq_donors_by_aint)) {
      donors_in_int <- uniq_donors_by_aint[[int_name]]
      if (int_name %in% old_intervals) {
        n_select <- min(target_n, length(donors_in_int))
        selected_donors <- c(selected_donors, sample(donors_in_int, n_select))
      } else {
        selected_donors <- c(selected_donors, donors_in_int)
      }
    }
  }
  cells_selected <- cells_keep[donor %in% selected_donors]
  seu_sub <- subset(seu, cells = cells_selected)
  message("Kept ", ncol(seu_sub), " cells after downsampling (originally ", ncol(seu), ")")
  return(seu_sub)
}

# 6. Downsample donors by (age_interval_30yr, sex) group (fallback function, currently unused)
seurat_downsample_donors_by_interval <- function(seu, donor_col = "donor_id", age_interval_col = "age_interval_30yr", sex_col = "sex", seed = 123) {
  meta <- seu@meta.data
  if (!donor_col %in% names(meta) || !age_interval_col %in% names(meta)) stop("columns not found")
  if (!sex_col %in% names(meta)) {
    warning("sex column not found, downsampling by age_interval only (ignoring sex)")
    aint <- as.character(meta[[age_interval_col]])
    donor <- as.character(meta[[donor_col]])
    keep <- !is.na(aint) & aint != "" & !is.na(donor) & donor != ""
    aint <- aint[keep]
    donor <- donor[keep]
    cells_keep <- rownames(meta)[keep]
    donors_by_aint <- split(donor, aint)
    uniq_donors_by_aint <- lapply(donors_by_aint, unique)
    n_donors <- sapply(uniq_donors_by_aint, length)
    message("Number of donors per age_interval (before downsampling):")
    print(n_donors)
    min_n <- min(n_donors)
    message("Downsampling each interval to ", min_n, " donors")
    set.seed(seed)
    selected_donors <- unlist(lapply(uniq_donors_by_aint, function(d) sample(d, min(min_n, length(d)))))
    cells_selected <- cells_keep[donor %in% selected_donors]
    seu_sub <- subset(seu, cells = cells_selected)
    message("Kept ", ncol(seu_sub), " cells after downsampling (originally ", ncol(seu), ")")
    return(seu_sub)
  }
  aint <- as.character(meta[[age_interval_col]])
  sex <- as.character(meta[[sex_col]])
  donor <- as.character(meta[[donor_col]])
  keep <- !is.na(aint) & aint != "" & !is.na(sex) & sex != "" & !is.na(donor) & donor != ""
  aint <- aint[keep]
  sex <- sex[keep]
  donor <- donor[keep]
  cells_keep <- rownames(meta)[keep]
  group <- paste(aint, sex, sep = "_")
  donors_by_group <- split(donor, group)
  uniq_donors_by_group <- lapply(donors_by_group, unique)
  n_donors <- sapply(uniq_donors_by_group, length)
  message("Number of donors per (age_interval, sex) group (before downsampling):")
  print(n_donors)
  min_n <- min(n_donors)
  message("Downsampling each (age_interval, sex) group to ", min_n, " donors")
  set.seed(seed)
  selected_donors <- unlist(lapply(uniq_donors_by_group, function(d) sample(d, min(min_n, length(d)))))
  cells_selected <- cells_keep[donor %in% selected_donors]
  seu_sub <- subset(seu, cells = cells_selected)
  message("Kept ", ncol(seu_sub), " cells after downsampling (originally ", ncol(seu), ")")
  return(seu_sub)
}

# 7. Elastic Net LODO (Leave-One-Donor-Out) cross-validation
elastic_net_LODO <- function(seu, donor_col = "donor_id", age_col = "age", alphas = c(0.1, 0.5, 1), ncores = 5, scale.factor = 1e4) {
  seu <- seurat_add_age_numeric_and_interval(seu, age_col = age_col)
  meta <- seu@meta.data
  if (any(is.na(meta$Age_numeric))) warning("Some ages could not be parsed to numeric")
  counts <- GetAssayData(seu, layer = "counts")
  if (is.null(counts)) counts <- GetAssayData(seu, layer = "counts", assay = "RNA")
  mt_genes <- grep("^MT-|^mt-", rownames(counts), value = TRUE)
  counts <- counts[!rownames(counts) %in% mt_genes, , drop = FALSE]
  gene_pct <- Matrix::rowSums(counts > 0) / ncol(counts)
  counts <- counts[gene_pct >= 0.01, , drop = FALSE]
  counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
  col_sums <- Matrix::colSums(counts)
  counts_norm <- sweep(counts, 2, col_sums, FUN = "/") * scale.factor
  X <- Matrix::t(log1p(counts_norm))
  y <- meta$Age_numeric
  donors <- as.character(meta[[donor_col]])
  cell_names <- rownames(meta)
  weights_all <- if ("weight" %in% colnames(meta)) meta$weight else rep(1, nrow(meta))
  if ("weight" %in% colnames(meta)) {
    message("Using cell weights (from seurat_add_balance_weights)")
  } else {
    message("No weight column detected, all cells weighted 1")
  }
  all_donors <- unique(donors)
  nfolds_inner <- 5
  process_one_donor <- function(test_donor) {
    idx_test <- which(donors == test_donor)
    idx_train <- which(donors != test_donor)
    train_donors <- unique(donors[idx_train])
    n_folds <- min(nfolds_inner, length(train_donors))
    donor_to_fold <- setNames(sample(rep(seq_len(n_folds), length.out = length(train_donors))), train_donors)
    foldid <- donor_to_fold[donors[idx_train]]
    X_train <- X[idx_train, , drop = FALSE]
    y_train <- y[idx_train]
    weights_train <- weights_all[idx_train]
    X_test <- X[idx_test, , drop = FALSE]
    X_train_scaled <- scale(X_train)
    center <- attr(X_train_scaled, "scaled:center")
    scale_val <- attr(X_train_scaled, "scaled:scale")
    scale_val[scale_val == 0] <- 1
    X_test_scaled <- sweep(X_test, 2, center, "-")
    X_test_scaled <- sweep(X_test_scaled, 2, scale_val, "/")
    best_cvm <- Inf
    best_fit <- NULL
    best_alpha <- NA_real_
    for (a in alphas) {
      fit <- tryCatch(cv.glmnet(X_train_scaled, y_train, alpha = a, foldid = foldid, weights = weights_train, standardize = FALSE), error = function(e) NULL)
      if (is.null(fit)) next
      cvm_min <- min(fit$cvm, na.rm = TRUE)
      if (cvm_min < best_cvm) {
        best_cvm <- cvm_min
        best_fit <- fit
        best_alpha <- a
      }
    }
    if (is.null(best_fit)) return(list(pred = NULL, alpha = NA_real_, lambda = NA_real_, ncoef = NA_integer_, coef = NULL))
    y_pred <- as.numeric(predict(best_fit, newx = X_test_scaled, s = "lambda.min"))
    coef_vec <- as.numeric(coef(best_fit, s = "lambda.min"))
    ncoef <- sum(coef_vec != 0)
    coef_df <- data.frame(gene = rownames(coef(best_fit, s = "lambda.min")), coef = coef_vec, stringsAsFactors = FALSE)
    coef_df <- coef_df[coef_df$coef != 0 & coef_df$gene != "(Intercept)", ]
    list(
      pred = data.frame(cell = cell_names[idx_test], actual = y[idx_test], predicted = y_pred, donor = test_donor, alpha_chosen = best_alpha, stringsAsFactors = FALSE),
      alpha = best_alpha,
      lambda = best_fit$lambda.min,
      ncoef = ncoef,
      coef = coef_df
    )
  }
  n_donors <- length(all_donors)
  if (ncores <= 1) {
    res_list <- vector("list", n_donors)
    pb <- utils::txtProgressBar(0, n_donors, style = 3)
    for (i in seq_len(n_donors)) {
      res_list[[i]] <- process_one_donor(all_donors[i])
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
  } else if (requireNamespace("pbmcapply", quietly = TRUE)) {
    res_list <- pbmcapply::pbmclapply(all_donors, process_one_donor, mc.cores = ncores)
  } else {
    message("LODO: running ", n_donors, " donors in parallel (install pbmcapply for progress bar)")
    res_list <- mclapply(all_donors, process_one_donor, mc.cores = ncores)
  }
  pred_df <- do.call(rbind, lapply(res_list, function(r) r$pred))
  if (is.null(pred_df) || nrow(pred_df) == 0) return(list(predictions = NULL, metrics = NULL, model_params = NULL, coefs = NULL))
  alpha_vec <- sapply(res_list, function(r) if (is.null(r$pred)) NA else r$alpha)
  lambda_vec <- sapply(res_list, function(r) if (is.null(r$pred)) NA else r$lambda)
  ncoef_vec <- sapply(res_list, function(r) if (is.null(r$pred)) NA else r$ncoef)
  coef_list <- lapply(seq_along(res_list), function(i){
    if (!is.null(res_list[[i]]$coef) && nrow(res_list[[i]]$coef) > 0) {
      df <- res_list[[i]]$coef
      df$donor <- all_donors[i]
      return(df)
    }
    NULL
  })
  coef_df <- do.call(rbind, coef_list[!sapply(coef_list, is.null)])
  rmse <- sqrt(mean((pred_df$actual - pred_df$predicted)^2))
  mae <- mean(abs(pred_df$actual - pred_df$predicted))
  pcc <- cor(pred_df$actual, pred_df$predicted, use = "pairwise.complete.obs")
  r2 <- 1 - sum((pred_df$actual - pred_df$predicted)^2) / sum((pred_df$actual - mean(pred_df$actual))^2)
  list(
    predictions = pred_df,
    metrics = c(RMSE = rmse, MAE = mae, PCC = pcc, R2 = r2),
    model_params = data.frame(donor_left_out = all_donors, alpha_chosen = alpha_vec, lambda_min = lambda_vec, n_nonzero_coef = ncoef_vec, stringsAsFactors = FALSE),
    coefs = coef_df
  )
}

message("========================================")
message("Aging Clocks Functions loaded successfully!")
message("Available functions:")
message("  - load_mtx_to_seurat()")
message("  - process_one_celltype()")
message("  - seurat_add_age_numeric_and_interval()")
message("  - seurat_add_balance_weights()")
message("  - seurat_downsample_old_donors_10yr()")
message("  - seurat_downsample_donors_by_interval()")
message("  - elastic_net_LODO()")
message("========================================")
