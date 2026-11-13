################################################################################
# Mouse Aging Clock Functions
# Full implementation of Nested Feature Selection + LOOCV
################################################################################

# Function: compute normalized pseudo-bulk data for each sample
# Input:
#   - seurat_obj: Seurat object (Seurat5)
#   - ambient_genes: list of ambient RNA genes (to be removed)
# Output:
#   - matrix (genes x samples), each value is log2(1 + depth_normalized_count)
Calculate_Pseudobulk_Normalized <- function(seurat_obj, ambient_genes = NULL) {
    require(Matrix)
    require(Seurat)
    
    # Seurat5: Join layers if necessary
    if(length(Layers(seurat_obj)) > 1) {
        message("  Multiple layers detected, joining...")
        seurat_obj <- JoinLayers(seurat_obj)
    }
    
    # Get count matrix and sample metadata
    counts <- LayerData(seurat_obj, layer = "counts")
    samples <- seurat_obj$sample
    
    if(is.null(samples)) {
        stop("Seurat object is missing the 'sample' column")
    }
    
    # Remove ambient genes
    if(!is.null(ambient_genes)) {
        ambient_genes <- intersect(ambient_genes, rownames(counts))
        if(length(ambient_genes) > 0) {
            message(sprintf("  Removing %d ambient genes", length(ambient_genes)))
            counts <- counts[!rownames(counts) %in% ambient_genes, ]
        }
    }
    
    # Compute pseudobulk for each sample (sum counts by sample)
    sample_names <- unique(samples)
    pseudobulk <- sapply(sample_names, function(s) {
        idx <- which(samples == s)
        Matrix::rowSums(counts[, idx, drop = FALSE])
    })
    
    # Depth normalization (scale to 1e4)
    depth <- colSums(pseudobulk)
    pseudobulk_norm <- sweep(pseudobulk, 2, depth / 1e4, FUN = "/")
    
    # Log2 transformation: log2(1 + x)
    pseudobulk_log2 <- log2(1 + pseudobulk_norm)
    
    message(sprintf("  Generated pseudobulk: %d genes x %d samples", nrow(pseudobulk_log2), ncol(pseudobulk_log2)))
    
    return(pseudobulk_log2)
}


# Function: compute the correlation between each gene and age
# Input:
#   - pseudobulk_matrix: pseudobulk matrix (genes x samples), from Calculate_Pseudobulk_Normalized
#   - method: correlation method, "pearson" or "spearman" (default "spearman")
# Output:
#   - data frame with gene, correlation, p_value, q_value, sorted by absolute correlation
Calculate_Gene_Age_Correlation <- function(pseudobulk_matrix, method = "spearman") {
    
    # Extract age from colnames (format: Wk106_Rep1 or Wk5_Rep2)
    sample_names <- colnames(pseudobulk_matrix)
    print(sample_names)
    # Extract the number after Wk
    ages <- as.numeric(gsub("^Wk([0-9]+)_.*", "\\1", sample_names))
    
    if(any(is.na(ages))) {
        stop("Cannot extract age from colnames; check that the format is 'Wk<number>_Rep'")
    }
    
    message(sprintf("  Detected %d samples, age range: %s", 
                    length(ages), paste(unique(sort(ages)), collapse = ", ")))
    
    # Compute gene-age correlation for each gene (vectorized to avoid memory issues)
    gene_names <- rownames(pseudobulk_matrix)
    n_genes <- nrow(pseudobulk_matrix)
    
    message(sprintf("  Computing age correlation for %d genes (%s)...", n_genes, method))
    
    # Vectorized correlation computation
    if (method == "spearman") {
        # Spearman: rank ages and the expression of each gene first
        ages_rank <- rank(ages)
        # Compute Spearman correlation row-wise on the matrix
        cor_vec <- apply(pseudobulk_matrix, 1, function(expr) {
            cor(ages_rank, rank(expr), method = "pearson")  # Spearman = Pearson of ranks
        })
    } else {
        # Pearson: compute directly
        cor_vec <- apply(pseudobulk_matrix, 1, function(expr) {
            cor(ages, expr, method = "pearson")
        })
    }
    
    # Compute p-value (using the cor.test formula)
    n <- length(ages)
    t_stat <- cor_vec * sqrt(n - 2) / sqrt(1 - cor_vec^2)
    p_values <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
    
    # Build result table
    cor_table <- data.frame(
        gene = gene_names,
        correlation = cor_vec,
        p_value = p_values,
        stringsAsFactors = FALSE
    )
    
    # Add FDR-adjusted q-value
    cor_table$q_value <- p.adjust(cor_table$p_value, method = "BH")
    
    # Sort by absolute correlation (descending)
    cor_table <- cor_table[order(abs(cor_table$correlation), decreasing = TRUE), ]
    
    message(sprintf("  Done! Number of significantly correlated genes (p < 0.01): %d", 
                    sum(cor_table$p_value < 0.01, na.rm = TRUE)))
    
    return(cor_table)
}


# Function: Nested Feature Selection + LOOCV for Aging Clock (Single-Cell Level, Optimized)
# Input:
#   - seurat_obj: Seurat object (containing single-cell data)
#   - whitelist_genes: whitelist genes (always included)
#   - ambient_genes: list of ambient RNA genes (to be removed)
#   - cor_pval_cutoff: p-value cutoff for gene-age correlation (default 0.05)
#   - cor_top_n: additionally keep the top N most correlated genes (default 500)
#   - min_mean_expr: mean pseudobulk expression cutoff (default 0.1, log2(CPM+1) scale)
#   - alphas: range of Elastic Net alpha values (default seq(0.1, 0.9, by=0.2))
#   - nfolds_inner: number of inner CV folds (for alpha tuning, default 5)
#   - standardize: whether to apply Z-score standardization (default TRUE)
#   - exclude_ages: vector of ages to exclude (e.g. c(24, 28); default NULL, no exclusion)
#   - exclude_samples: vector of sample names to exclude (e.g. c("6mo_rep3", "24mo_rep1"); default NULL, no exclusion)
#   - n_cores: number of parallel cores (NULL=auto-detect, 1=serial, >1=parallel)
# Output:
#   - list containing:
#     * cell_predictions: cell-level predictions (with predicted and predicted_calibrated)
#     * sample_predictions: sample-level summary (with predicted and predicted_calibrated)
#     * selected_genes: genes selected in each LOOCV round
#     * best_alphas: best alpha in each LOOCV round
#     * calibration_params: calibration parameters in each LOOCV round (intercept, slope)
#     * metrics: raw and calibrated performance metrics
# Features:
#   - Three-step feature filtering: p-value filter → Top N filter → low-expression filter
#   - Low-expression filtering is applied last; whitelist genes must also pass the expression cutoff
#   - Uses weights to balance the contribution of each sample and each age
#   - Avoids information leakage in feature selection (features selected independently in each LOOCV round)
#   - Uses lambda.1se (more conservative, prevents overfitting) instead of lambda.min
#   - Supports parallel processing (requires the pbmcapply package)
#   - Calibration performed independently in each LOOCV round (fit on training set, applied to test set)
Nested_FeatureSelection_LOOCV <- function(seurat_obj, 
                                          whitelist_genes = NULL,
                                          ambient_genes = NULL,
                                          cor_pval_cutoff = 0.05,
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
    
    # Check parallel package and set number of cores
    use_parallel <- TRUE
    
    # Auto-detect if the user did not specify the number of cores
    if (is.null(n_cores)) {
        if (!requireNamespace("pbmcapply", quietly = TRUE)) {
            message("Note: install the pbmcapply package to enable parallel processing and progress bars")
            message("      install.packages('pbmcapply')")
            use_parallel <- FALSE
            n_cores <- 1
        } else {
            n_cores <- parallel::detectCores() - 1
            message(sprintf("Auto-detected %d available cores", parallel::detectCores()))
        }
    }
    
    # Validate number of cores
    if (n_cores < 1) {
        n_cores <- 1
        use_parallel <- FALSE
    } else if (n_cores == 1) {
        use_parallel <- FALSE
    } else {
        # Check pbmcapply package
        if (!requireNamespace("pbmcapply", quietly = TRUE)) {
            message("Warning: pbmcapply is not installed, falling back to serial mode")
            message("      install.packages('pbmcapply')")
            use_parallel <- FALSE
            n_cores <- 1
        }
    }
    
    message("========================================")
    message("Starting Nested Feature Selection + LOOCV (Optimized)")
    if(use_parallel) {
        message(sprintf("Parallel mode: using %d cores", n_cores))
    } else {
        message("Serial mode")
    }
    message("========================================")
    
    # 1. Prepare data
    message("\n[Step 1/5] Preparing data...")
    
    # Join layers (Seurat5)
    if(length(Layers(seurat_obj)) > 1) {
        message("  Multiple layers detected, joining...")
        seurat_obj <- JoinLayers(seurat_obj)
    }
    
    # Get single-cell data
    counts <- LayerData(seurat_obj, layer = "counts")
    cell_samples <- seurat_obj$sample
    cell_ages <- seurat_obj$age
    
    if(is.null(cell_samples) || is.null(cell_ages)) {
        stop("Seurat object is missing the 'sample' or 'age' column")
    }
    
    # Filter samples of the specified ages (if specified)
    if (!is.null(exclude_ages)) {
        cells_to_keep <- !cell_ages %in% exclude_ages
        n_cells_removed <- sum(!cells_to_keep)
        
        if (n_cells_removed > 0) {
            # Count samples for each excluded age
            excluded_summary <- sapply(exclude_ages, function(age) {
                n_samples <- length(unique(cell_samples[cell_ages == age]))
                n_cells <- sum(cell_ages == age)
                sprintf("%dmo(%d samples,%d cells)", age, n_samples, n_cells)
            })
            
            message(sprintf("  Excluding specified ages: %s", paste(excluded_summary, collapse=", ")))
            message(sprintf("  Removed %d cells in total", n_cells_removed))
            
            # Filter data
            counts <- counts[, cells_to_keep, drop = FALSE]
            cell_samples <- cell_samples[cells_to_keep]
            cell_ages <- cell_ages[cells_to_keep]
            
            # Also update seurat_obj metadata (for subsequent pseudobulk computation)
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
        } else {
            message(sprintf("  No samples of the specified ages (%s) found, no filtering needed", 
                            paste(exclude_ages, collapse=", ")))
        }
    }
    
    # Filter specified samples (if specified)
    if (!is.null(exclude_samples)) {
        cells_to_keep <- !cell_samples %in% exclude_samples
        n_cells_removed <- sum(!cells_to_keep)
        
        if (n_cells_removed > 0) {
            # Summarize each excluded sample
            excluded_summary <- sapply(exclude_samples, function(sample) {
                if (sample %in% cell_samples) {
                    n_cells <- sum(cell_samples == sample)
                    age <- unique(cell_ages[cell_samples == sample])[1]
                    sprintf("%s(%dmo,%d cells)", sample, age, n_cells)
                } else {
                    sprintf("%s(not found)", sample)
                }
            })
            
            message(sprintf("  Excluding specified samples: %s", paste(excluded_summary, collapse=", ")))
            message(sprintf("  Removed %d cells in total", n_cells_removed))
            
            # Filter data
            counts <- counts[, cells_to_keep, drop = FALSE]
            cell_samples <- cell_samples[cells_to_keep]
            cell_ages <- cell_ages[cells_to_keep]
            
            # Also update seurat_obj metadata (for subsequent pseudobulk computation)
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
        } else {
            message(sprintf("  Specified samples (%s) not found, no filtering needed", 
                            paste(exclude_samples, collapse=", ")))
        }
    }
    
    # Remove ambient genes
    if(!is.null(ambient_genes)) {
        ambient_genes <- intersect(ambient_genes, rownames(counts))
        if(length(ambient_genes) > 0) {
            message(sprintf("  Removing %d ambient genes", length(ambient_genes)))
            counts <- counts[!rownames(counts) %in% ambient_genes, ]
        }
    }
    
    # Depth normalization + Log2 transformation (single-cell level)
    message("  Applying depth normalization + log2 transformation to single-cell data...")
    depth <- Matrix::colSums(counts)
    counts_norm <- sweep(counts, 2, depth / 1e4, FUN = "/")
    counts_log2 <- log2(1 + counts_norm)
    
    message(sprintf("  Cells: %d, genes: %d", ncol(counts_log2), nrow(counts_log2)))
    
    # 2. Compute pseudobulk (for feature selection)
    message("\n[Step 2/5] Computing pseudobulk (for feature selection)...")
    pseudobulk_all <- Calculate_Pseudobulk_Normalized(seurat_obj, ambient_genes)
    
    # Get the list of unique samples
    unique_samples <- unique(cell_samples)
    sample_ages <- sapply(unique_samples, function(s) unique(cell_ages[cell_samples == s])[1])
    n_samples <- length(unique_samples)
    
    message(sprintf("  Samples: %d, age range: %s", n_samples, paste(unique(sort(sample_ages)), collapse=", ")))
    
    # 3. LOOCV loop (by sample, parallelized)
    message("\n[Step 3/5] Starting Leave-One-Sample-Out Cross-Validation...")
    
    # Define the function for a single LOOCV round
    process_one_loocv <- function(i) {
        test_sample <- unique_samples[i]
        train_samples <- unique_samples[-i]
        
        # 3.1 Feature Selection: compute gene-age correlation on training-set pseudobulk
        train_pseudobulk <- pseudobulk_all[, train_samples, drop = FALSE]
        
        # Compute gene-age correlation (on all genes, output suppressed)
        cor_results <- suppressMessages(Calculate_Gene_Age_Correlation(train_pseudobulk, method = "spearman"))
        
        # Strategy 1: significantly correlated genes (p < cutoff)
        sig_genes_pval <- cor_results$gene[cor_results$p_value < cor_pval_cutoff]
        
        # Strategy 2: top N most correlated genes (by absolute value)
        cor_results_sorted <- cor_results[order(abs(cor_results$correlation), decreasing = TRUE), ]
        sig_genes_top <- head(cor_results_sorted$gene, cor_top_n)
        
        # Combine both strategies
        sig_genes <- unique(c(sig_genes_pval, sig_genes_top))
        
        # Add whitelist genes
        if (!is.null(whitelist_genes)) {
            whitelist_available <- intersect(whitelist_genes, rownames(counts_log2))
            candidate_genes_all <- unique(c(sig_genes, whitelist_available))
        } else {
            candidate_genes_all <- sig_genes
        }
        
        # Check the number of candidate genes
        if (length(candidate_genes_all) < 10) {
            candidate_genes_all <- head(cor_results_sorted$gene, min(100, nrow(cor_results)))
        }
        
        # Finally filter low-expression genes (applied to all genes, including whitelist)
        gene_mean_expr <- rowMeans(train_pseudobulk)
        candidate_genes_expr <- gene_mean_expr[candidate_genes_all]
        
        # Filter low-expression genes
        candidate_genes <- candidate_genes_all[!is.na(candidate_genes_expr) & candidate_genes_expr >= min_mean_expr]
        
        # Record filtering statistics
        n_filtered_out_lowexpr <- sum(candidate_genes_expr < min_mean_expr, na.rm = TRUE)
        
        # Check whether candidate genes are empty
        if (length(candidate_genes) == 0) {
            stop("No candidate genes left after filtering! Lower min_mean_expr or cor_pval_cutoff")
        }
        
        # 3.2 Prepare training data (single-cell level, candidate genes only)
        train_cell_idx <- which(cell_samples %in% train_samples)
        test_cell_idx <- which(cell_samples == test_sample)
        
        # Check cell indices
        if (length(train_cell_idx) == 0) {
            stop(sprintf("Training set has no cells! train_samples: %s", paste(train_samples, collapse=", ")))
        }
        if (length(test_cell_idx) == 0) {
            stop(sprintf("Test sample %s has no cells!", test_sample))
        }
        
        # Ensure candidate genes exist in counts_log2
        candidate_genes <- intersect(candidate_genes, rownames(counts_log2))
        
        if (length(candidate_genes) == 0) {
            stop("Candidate genes not found in single-cell data!")
        }
        
        # Safely extract data
        X_train <- t(as.matrix(counts_log2[candidate_genes, train_cell_idx, drop = FALSE]))
        y_train <- cell_ages[train_cell_idx]
        train_samples_vec <- cell_samples[train_cell_idx]
        
        # Z-score standardization (based on training set)
        if (standardize) {
            train_means <- colMeans(X_train)
            train_sds <- apply(X_train, 2, sd)
            train_sds[train_sds == 0] <- 1  # avoid division by 0
            X_train <- scale(X_train, center = train_means, scale = train_sds)
        }
        
        # Compute weights: balance sample and age
        weights <- rep(1, length(train_cell_idx))
        
        # Iterative balancing (IPF-like)
        n_iter <- 10
        for (iter in 1:n_iter) {
            # Balance sample: equal total weight per sample
            sample_weight_sum <- tapply(weights, train_samples_vec, sum)
            target_per_sample <- mean(sample_weight_sum)
            sample_scale <- target_per_sample / sample_weight_sum[train_samples_vec]
            weights <- weights * sample_scale
            
            # Balance age: equal total weight per age
            age_weight_sum <- tapply(weights, y_train, sum)
            target_per_age <- mean(age_weight_sum)
            age_scale <- target_per_age / age_weight_sum[as.character(y_train)]
            weights <- weights * age_scale
        }
        
        # Normalize weights (sum = number of cells)
        weights <- weights * length(weights) / sum(weights)
        
        # 3.3 Nested CV: tune alpha (inner CV, on single cells)
        best_alpha <- NULL
        best_cv_error <- Inf
        best_lambda <- NULL
        
        for (alpha in alphas) {
            cv_fit <- tryCatch({
                cv.glmnet(X_train, y_train, 
                         weights = weights,
                         alpha = alpha, 
                         nfolds = nfolds_inner,
                         family = "gaussian",
                         standardize = FALSE)  # already standardized manually
            }, error = function(e) {
                return(NULL)
            })
            
            if (!is.null(cv_fit)) {
                min_error <- min(cv_fit$cvm)
                if (min_error < best_cv_error) {
                    best_cv_error <- min_error
                    best_alpha <- alpha
                    best_lambda <- cv_fit$lambda.1se
                }
            }
        }
        
        if (is.null(best_alpha)) {
            best_alpha <- 0.5
            best_lambda <- NULL
        }
        
        # 3.4 Train the final model with the best alpha (on all training cells)
        final_model <- glmnet(X_train, y_train, 
                             weights = weights,
                             alpha = best_alpha, 
                             family = "gaussian",
                             standardize = FALSE)
        
        # 3.5 Predict on training set (for calibration)
        if (!is.null(best_lambda)) {
            train_pred <- predict(final_model, newx = X_train, s = best_lambda)
        } else {
            train_pred <- predict(final_model, newx = X_train, s = min(final_model$lambda))
        }
        
        # 3.6 Fit calibration model on training set
        train_pred_df <- data.frame(
            actual = y_train,
            predicted = as.numeric(train_pred)
        )
        # Aggregate to sample level by sample
        train_pred_sample <- aggregate(cbind(actual, predicted) ~ train_samples_vec, 
                                       data = cbind(train_pred_df, train_samples_vec), 
                                       FUN = mean)
        colnames(train_pred_sample)[1] <- "sample"
        
        # Fit calibration model
        calib_model_i <- lm(actual ~ predicted, data = train_pred_sample)
        calib_intercept <- coef(calib_model_i)[1]
        calib_slope <- coef(calib_model_i)[2]
        
        # 3.7 Predict all cells of the test sample
        X_test <- t(as.matrix(counts_log2[candidate_genes, test_cell_idx, drop = FALSE]))
        
        # Check test data dimensions
        if (ncol(X_test) != ncol(X_train)) {
            stop(sprintf("Number of test features (%d) does not match training set (%d)!", ncol(X_test), ncol(X_train)))
        }
        
        # Use the same standardization parameters
        if (standardize) {
            X_test <- scale(X_test, center = train_means, scale = train_sds)
        }
        
        # Predict with the best lambda
        if (!is.null(best_lambda)) {
            pred <- predict(final_model, newx = X_test, s = best_lambda)
        } else {
            pred <- predict(final_model, newx = X_test, s = min(final_model$lambda))
        }
        
        # 3.8 Apply calibration to the test sample
        pred_calibrated <- calib_intercept + calib_slope * as.numeric(pred)
        
        # Return results
        return(list(
            test_sample = test_sample,
            test_cell_idx = test_cell_idx,
            predictions = as.numeric(pred),
            predictions_calibrated = as.numeric(pred_calibrated),
            selected_genes = candidate_genes,
            best_alpha = best_alpha,
            best_cv_error = best_cv_error,
            calib_intercept = calib_intercept,
            calib_slope = calib_slope,
            n_genes = length(candidate_genes),
            n_genes_filtered_lowexpr = n_filtered_out_lowexpr,
            n_train_cells = length(train_cell_idx),
            n_test_cells = length(test_cell_idx),
            actual_age = sample_ages[i]
        ))
    }
    
           # Run LOOCV in parallel or serially
           if (use_parallel) {
               message(sprintf("Parallel mode: using %d cores, initializing...\n", n_cores))
               
               # Use pbmcapply::pbmclapply (with progress bar)
               loocv_results <- tryCatch({
                   pbmcapply::pbmclapply(1:n_samples, process_one_loocv, 
                                        mc.cores = n_cores, 
                                        mc.set.seed = TRUE,
                                        mc.preschedule = FALSE,
                                        ignore.interactive = TRUE)
               }, error = function(e) {
                   message("\nWarning: pbmcapply failed, trying parallel::mclapply...")
                   message(sprintf("Error message: %s", e$message))
                   
                   # Fallback: use parallel::mclapply (no progress bar)
                   tryCatch({
                       message("Using parallel::mclapply (no progress bar)...\n")
                       parallel::mclapply(1:n_samples, process_one_loocv, 
                                         mc.cores = n_cores, 
                                         mc.set.seed = TRUE,
                                         mc.preschedule = FALSE)
                   }, error = function(e2) {
                       message("\nWarning: parallel execution failed, falling back to serial mode")
                       message(sprintf("Error message: %s", e2$message))
                       message("Re-running in serial mode...\n")
                       
                       # Final fallback to serial
                       lapply(1:n_samples, function(i) {
                           message(sprintf(">>> LOOCV Round %d/%d <<<", i, n_samples))
                           result <- process_one_loocv(i)
                           message(sprintf("    Done: %s | genes=%d | pred=%.2f, actual=%.2f", 
                                           result$test_sample, result$n_genes, 
                                           mean(result$predictions), result$actual_age))
                           return(result)
                       })
                   })
               })
           } else {
               message("Serial mode: showing detailed progress for each round...\n")
               loocv_results <- lapply(1:n_samples, function(i) {
                   message(sprintf(">>> LOOCV Round %d/%d <<<", i, n_samples))
                   result <- process_one_loocv(i)
                   message(sprintf("    Done: %s | genes=%d | pred=%.2f→%.2f, actual=%.2f", 
                                   result$test_sample, result$n_genes, 
                                   mean(result$predictions), mean(result$predictions_calibrated),
                                   result$actual_age))
                   return(result)
               })
           }
    
    # Combine results
    message("\n[Step 3.5/5] Combining parallel results...")
    
    # Check result validity
    valid_results <- sapply(loocv_results, function(x) is.list(x) && !is.null(x$test_sample))
    if (!all(valid_results)) {
        invalid_idx <- which(!valid_results)
        stop(sprintf("LOOCV Round %s returned invalid results; check the error messages", 
                    paste(invalid_idx, collapse=", ")))
    }
    
    cell_predictions <- data.frame(
        cell_id = colnames(counts_log2),
        sample = cell_samples,
        actual = cell_ages,
        predicted = NA,
        predicted_calibrated = NA,
        stringsAsFactors = FALSE
    )
    selected_genes_list <- list()
    best_alphas <- numeric(n_samples)
    calibration_params <- data.frame(
        sample = character(n_samples),
        intercept = numeric(n_samples),
        slope = numeric(n_samples),
        stringsAsFactors = FALSE
    )
    
    for (i in 1:n_samples) {
        result <- loocv_results[[i]]
        cell_predictions$predicted[result$test_cell_idx] <- result$predictions
        cell_predictions$predicted_calibrated[result$test_cell_idx] <- result$predictions_calibrated
        selected_genes_list[[result$test_sample]] <- result$selected_genes
        best_alphas[i] <- result$best_alpha
        calibration_params$sample[i] <- result$test_sample
        calibration_params$intercept[i] <- result$calib_intercept
        calibration_params$slope[i] <- result$calib_slope
        
        # Print sample-level summary
        sample_pred_mean <- mean(result$predictions, na.rm = TRUE)
        sample_pred_calib_mean <- mean(result$predictions_calibrated, na.rm = TRUE)
        message(sprintf("  %s: pred=%.2f→%.2f, actual=%.2f, err=%.2f→%.2f | genes=%d(-%d), calib=(%.2f,%.2f)", 
                       result$test_sample, 
                       sample_pred_mean,
                       sample_pred_calib_mean,
                       result$actual_age, 
                       abs(sample_pred_mean - result$actual_age),
                       abs(sample_pred_calib_mean - result$actual_age),
                       result$n_genes,
                       result$n_genes_filtered_lowexpr,
                       result$calib_intercept,
                       result$calib_slope))
    }
    
    # 4. Summarize sample-level results
    message("\n[Step 4/5] Summarizing sample-level predictions...")
    cell_pred_complete <- cell_predictions[complete.cases(cell_predictions), ]
    
    # Summarize raw and calibrated predictions
    sample_predictions <- aggregate(cbind(predicted, predicted_calibrated) ~ sample + actual, 
                                   data = cell_pred_complete, 
                                   FUN = mean)
    
    # 5. Compute overall performance
    message("\n[Step 5/5] Computing overall prediction performance...")
    
    # Sample-level performance - raw predictions
    rmse_sample <- sqrt(mean((sample_predictions$actual - sample_predictions$predicted)^2))
    mae_sample <- mean(abs(sample_predictions$actual - sample_predictions$predicted))
    pcc_sample <- cor(sample_predictions$actual, sample_predictions$predicted, method = "pearson")
    r2_sample <- 1 - sum((sample_predictions$actual - sample_predictions$predicted)^2) / 
                     sum((sample_predictions$actual - mean(sample_predictions$actual))^2)
    
    # Sample-level performance - calibrated predictions
    rmse_sample_calib <- sqrt(mean((sample_predictions$actual - sample_predictions$predicted_calibrated)^2))
    mae_sample_calib <- mean(abs(sample_predictions$actual - sample_predictions$predicted_calibrated))
    pcc_sample_calib <- cor(sample_predictions$actual, sample_predictions$predicted_calibrated, method = "pearson")
    r2_sample_calib <- 1 - sum((sample_predictions$actual - sample_predictions$predicted_calibrated)^2) / 
                           sum((sample_predictions$actual - mean(sample_predictions$actual))^2)
    
    message("  Sample-level performance (raw):")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_sample, mae_sample, pcc_sample, r2_sample))
    
    message("  Sample-level performance (calibrated):")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_sample_calib, mae_sample_calib, pcc_sample_calib, r2_sample_calib))
    
    # Cell-level performance - raw predictions
    rmse_cell <- sqrt(mean((cell_pred_complete$actual - cell_pred_complete$predicted)^2))
    mae_cell <- mean(abs(cell_pred_complete$actual - cell_pred_complete$predicted))
    pcc_cell <- cor(cell_pred_complete$actual, cell_pred_complete$predicted, method = "pearson")
    r2_cell <- 1 - sum((cell_pred_complete$actual - cell_pred_complete$predicted)^2) / 
                   sum((cell_pred_complete$actual - mean(cell_pred_complete$actual))^2)
    
    # Cell-level performance - calibrated predictions
    rmse_cell_calib <- sqrt(mean((cell_pred_complete$actual - cell_pred_complete$predicted_calibrated)^2))
    mae_cell_calib <- mean(abs(cell_pred_complete$actual - cell_pred_complete$predicted_calibrated))
    pcc_cell_calib <- cor(cell_pred_complete$actual, cell_pred_complete$predicted_calibrated, method = "pearson")
    r2_cell_calib <- 1 - sum((cell_pred_complete$actual - cell_pred_complete$predicted_calibrated)^2) / 
                         sum((cell_pred_complete$actual - mean(cell_pred_complete$actual))^2)
    
    message("  Cell-level performance (raw):")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_cell, mae_cell, pcc_cell, r2_cell))
    
    message("  Cell-level performance (calibrated):")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_cell_calib, mae_cell_calib, pcc_cell_calib, r2_cell_calib))
    
    message("\n========================================")
    message("Done!")
    message("========================================\n")
    
    # 6. Additional statistics
    message("\n[Additional info] Feature selection statistics:")
    n_genes_selected <- sapply(selected_genes_list, length)
    message(sprintf("  Genes selected per round: mean=%.1f, median=%.0f, range=[%d, %d]",
                   mean(n_genes_selected), median(n_genes_selected),
                   min(n_genes_selected), max(n_genes_selected)))
    
    n_genes_filtered <- sapply(loocv_results, function(x) x$n_genes_filtered_lowexpr)
    message(sprintf("  Genes removed by low-expression filter (mean_expr < %.2f): mean=%.1f, range=[%d, %d]",
                   min_mean_expr,
                   mean(n_genes_filtered), 
                   min(n_genes_filtered), max(n_genes_filtered)))
    
    # Find the most frequently selected genes
    all_selected <- unlist(selected_genes_list)
    gene_freq <- sort(table(all_selected), decreasing = TRUE)
    message(sprintf("  Top 5 most frequently selected genes: %s", 
                   paste(paste0(names(head(gene_freq, 5)), "(", head(gene_freq, 5), ")"), collapse=", ")))
    
    return(list(
        cell_predictions = cell_pred_complete,
        sample_predictions = sample_predictions,
        selected_genes = selected_genes_list,
        gene_selection_freq = gene_freq,
        best_alphas = best_alphas,
        calibration_params = calibration_params,
        metrics = data.frame(
            level = c("sample", "sample_calibrated", "cell", "cell_calibrated"),
            RMSE = c(rmse_sample, rmse_sample_calib, rmse_cell, rmse_cell_calib),
            MAE = c(mae_sample, mae_sample_calib, mae_cell, mae_cell_calib),
            PCC = c(pcc_sample, pcc_sample_calib, pcc_cell, pcc_cell_calib),
            R2 = c(r2_sample, r2_sample_calib, r2_cell, r2_cell_calib)
        ),
        parameters = list(
            cor_pval_cutoff = cor_pval_cutoff,
            cor_top_n = cor_top_n,
            min_mean_expr = min_mean_expr,
            alphas = alphas,
            nfolds_inner = nfolds_inner,
            standardize = standardize,
            exclude_ages = exclude_ages,
            exclude_samples = exclude_samples,
            use_weights = TRUE,
            parallel = use_parallel,
            n_cores = n_cores,
            calibration = "per_loocv"
        )
    ))
}

################################################################################
# Usage examples
################################################################################
# 
# # Load this file
# source("Mouse_AgingClock_Functions.R")
# 
# # Example 1: auto-detect number of cores (recommended)
# results <- Nested_FeatureSelection_LOOCV(
#     seurat_obj = your_seurat_object,
#     whitelist_genes = your_whitelist,      # optional
#     ambient_genes = ambient_gene_list,     # optional
#     cor_pval_cutoff = 0.05,                # p-value cutoff
#     cor_top_n = 500,                       # top N correlated genes
#     min_mean_expr = 0.1,                   # minimum mean expression
#     alphas = seq(0.1, 0.9, by = 0.2),      # Elastic Net alpha range
#     nfolds_inner = 5,                      # number of inner CV folds
#     standardize = TRUE,                    # Z-score standardization
#     n_cores = NULL                         # NULL=auto-detect
# )
# 
# # Example 2: manually specify number of cores
# results <- Nested_FeatureSelection_LOOCV(
#     seurat_obj = your_seurat_object,
#     n_cores = 8                            # use 8 cores
# )
# 
# # Example 3: serial mode (single core)
# results <- Nested_FeatureSelection_LOOCV(
#     seurat_obj = your_seurat_object,
#     n_cores = 1                            # serial mode
# )
# 
# # View results
# print(results$metrics)
# print(results$sample_predictions)
# print(results$calibration_params)
# 
# # Save results
# write.csv(results$cell_predictions, "cell_predictions.csv", row.names = FALSE)
# write.csv(results$sample_predictions, "sample_predictions.csv", row.names = FALSE)
# write.csv(results$calibration_params, "calibration_params.csv", row.names = FALSE)
# saveRDS(results, "aging_clock_results.rds")
# 
# # Notes on setting the number of cores:
# # - n_cores = NULL: auto-detect, uses (total cores - 1)
# # - n_cores = 1: serial mode, lowest memory usage
# # - n_cores = 4: use 4 cores
# # - n_cores = 8: use 8 cores (recommended for servers)
# # - n_cores = 16: use 16 cores (high-performance servers)
# 
# # Parameter tuning tips:
# # 1. If performance is poor (R2 < 0.7), try:
# #    - increasing cor_top_n (e.g. 1000)
# #    - lowering min_mean_expr (e.g. 0.05)
# #    - using a denser alphas grid (e.g. seq(0.1, 0.9, by=0.1))
# # 2. If overfitting (good on training, poor on test), try:
# #    - reducing cor_top_n (e.g. 300)
# #    - increasing min_mean_expr (e.g. 0.5)
# # 3. If too slow, try:
# #    - increasing n_cores
# #    - reducing cor_top_n
# #    - reducing the number of alphas
# 
################################################################################
