################################################################################
# Mouse Aging Clock Functions
# 包含 Nested Feature Selection + LOOCV 的完整实现
################################################################################

# Function: 计算每个 sample 的 pseudo-bulk data（normalized）
# 输入:
#   - seurat_obj: Seurat 对象（Seurat5）
#   - ambient_genes: ambient RNA 基因列表（要去掉的基因）
# 输出:
#   - 矩阵（genes x samples），每个值是 log2(1 + depth_normalized_count)
Calculate_Pseudobulk_Normalized <- function(seurat_obj, ambient_genes = NULL) {
    require(Matrix)
    require(Seurat)
    
    # Seurat5: Join layers if necessary
    if(length(Layers(seurat_obj)) > 1) {
        message("  检测到多个 layers，正在合并...")
        seurat_obj <- JoinLayers(seurat_obj)
    }
    
    # 获取 count matrix 和 sample metadata
    counts <- LayerData(seurat_obj, layer = "counts")
    samples <- seurat_obj$sample
    
    if(is.null(samples)) {
        stop("Seurat object 缺少 'sample' 列")
    }
    
    # 去掉 ambient genes
    if(!is.null(ambient_genes)) {
        ambient_genes <- intersect(ambient_genes, rownames(counts))
        if(length(ambient_genes) > 0) {
            message(sprintf("  去掉 %d 个 ambient genes", length(ambient_genes)))
            counts <- counts[!rownames(counts) %in% ambient_genes, ]
        }
    }
    
    # 计算每个 sample 的 pseudobulk（按 sample 汇总 counts）
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
    
    message(sprintf("  生成 pseudobulk: %d genes x %d samples", nrow(pseudobulk_log2), ncol(pseudobulk_log2)))
    
    return(pseudobulk_log2)
}


# Function: 计算每个基因与年龄的相关系数
# 输入:
#   - pseudobulk_matrix: pseudobulk 矩阵（genes x samples），来自 Calculate_Pseudobulk_Normalized
#   - method: 相关系数方法，"pearson" 或 "spearman"（默认 "spearman"）
# 输出:
#   - 数据框，包含 gene, correlation, p_value, q_value, 按相关系数绝对值排序
Calculate_Gene_Age_Correlation <- function(pseudobulk_matrix, method = "spearman") {
    
    # 从 colnames 提取年龄信息（格式：Wk106_Rep1 或 Wk5_Rep2）
    sample_names <- colnames(pseudobulk_matrix)
    print(sample_names)
    # 提取 Wk 后面的数字
    ages <- as.numeric(gsub("^Wk([0-9]+)_.*", "\\1", sample_names))
    
    if(any(is.na(ages))) {
        stop("无法从 colnames 提取年龄信息，请检查格式是否为 'Wk数字_Rep'")
    }
    
    message(sprintf("  检测到 %d 个样本，年龄范围: %s", 
                    length(ages), paste(unique(sort(ages)), collapse = ", ")))
    
    # 计算每个基因与年龄的相关系数（向量化，避免内存问题）
    gene_names <- rownames(pseudobulk_matrix)
    n_genes <- nrow(pseudobulk_matrix)
    
    message(sprintf("  正在计算 %d 个基因的年龄相关系数 (%s)...", n_genes, method))
    
    # 向量化计算相关系数
    if (method == "spearman") {
        # Spearman: 先对 ages 和每个基因的表达量排序
        ages_rank <- rank(ages)
        # 对矩阵按行计算 Spearman 相关
        cor_vec <- apply(pseudobulk_matrix, 1, function(expr) {
            cor(ages_rank, rank(expr), method = "pearson")  # Spearman = Pearson of ranks
        })
    } else {
        # Pearson: 直接计算
        cor_vec <- apply(pseudobulk_matrix, 1, function(expr) {
            cor(ages, expr, method = "pearson")
        })
    }
    
    # 计算 p-value（使用 cor.test 的公式）
    n <- length(ages)
    t_stat <- cor_vec * sqrt(n - 2) / sqrt(1 - cor_vec^2)
    p_values <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
    
    # 构建结果表
    cor_table <- data.frame(
        gene = gene_names,
        correlation = cor_vec,
        p_value = p_values,
        stringsAsFactors = FALSE
    )
    
    # 添加 FDR 校正的 q-value
    cor_table$q_value <- p.adjust(cor_table$p_value, method = "BH")
    
    # 按相关系数绝对值排序（从大到小）
    cor_table <- cor_table[order(abs(cor_table$correlation), decreasing = TRUE), ]
    
    message(sprintf("  完成！显著相关基因数 (p < 0.01): %d", 
                    sum(cor_table$p_value < 0.01, na.rm = TRUE)))
    
    return(cor_table)
}


# Function: Nested Feature Selection + LOOCV for Aging Clock (Single-Cell Level, Optimized)
# 输入:
#   - seurat_obj: Seurat 对象（包含单细胞数据）
#   - whitelist_genes: 白名单基因（必须包含的基因）
#   - ambient_genes: ambient RNA 基因列表（要去掉的基因）
#   - cor_pval_cutoff: 基因-年龄相关性的 p-value 阈值（默认 0.05）
#   - cor_top_n: 额外保留相关性最强的 top N 基因（默认 500）
#   - min_mean_expr: pseudobulk 平均表达量阈值（默认 0.1，log2(CPM+1) scale）
#   - alphas: Elastic Net 的 alpha 参数范围（默认 seq(0.1, 0.9, by=0.2)）
#   - nfolds_inner: 内层 CV 折数（用于 alpha 调优，默认 5）
#   - standardize: 是否 Z-score 标准化（默认 TRUE）
#   - exclude_ages: 要排除的年龄向量（如 c(24, 28)，默认 NULL 不排除）
#   - exclude_samples: 要排除的样本名称向量（如 c("6mo_rep3", "24mo_rep1")，默认 NULL 不排除）
#   - n_cores: 并行核心数（NULL=自动检测，1=串行，>1=并行）
# 输出:
#   - 列表，包含:
#     * cell_predictions: 细胞级预测（含 predicted 和 predicted_calibrated）
#     * sample_predictions: 样本级汇总（含 predicted 和 predicted_calibrated）
#     * selected_genes: 每次 LOOCV 选择的基因
#     * best_alphas: 每次 LOOCV 的最佳 alpha
#     * calibration_params: 每次 LOOCV 的校准参数（intercept, slope）
#     * metrics: 原始和校准后的性能指标
# 特点:
#   - 三层特征过滤：p-value 筛选 → Top N 筛选 → 低表达过滤
#   - 低表达过滤在最后，包括 whitelist 基因也要通过表达量阈值
#   - 使用 weights 平衡每个 sample 和每个 age 的贡献
#   - 避免特征选择的信息泄露（每次 LOOCV 独立选择特征）
#   - 使用 lambda.1se（更保守，防止过拟合）而不是 lambda.min
#   - 支持并行处理（需要 pbmcapply 包）
#   - 每次 LOOCV 独立执行校准（在训练集上拟合，应用到测试集）
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
    
    # 检查并行包和设置核心数
    use_parallel <- TRUE
    
    # 如果用户没有指定核心数，自动检测
    if (is.null(n_cores)) {
        if (!requireNamespace("pbmcapply", quietly = TRUE)) {
            message("提示: 安装 pbmcapply 包可以启用并行处理和进度条")
            message("      install.packages('pbmcapply')")
            use_parallel <- FALSE
            n_cores <- 1
        } else {
            n_cores <- parallel::detectCores() - 1
            message(sprintf("自动检测到 %d 个可用核心", parallel::detectCores()))
        }
    }
    
    # 验证核心数
    if (n_cores < 1) {
        n_cores <- 1
        use_parallel <- FALSE
    } else if (n_cores == 1) {
        use_parallel <- FALSE
    } else {
        # 检查 pbmcapply 包
        if (!requireNamespace("pbmcapply", quietly = TRUE)) {
            message("警告: pbmcapply 包未安装，降级为串行模式")
            message("      install.packages('pbmcapply')")
            use_parallel <- FALSE
            n_cores <- 1
        }
    }
    
    message("========================================")
    message("开始 Nested Feature Selection + LOOCV (Optimized)")
    if(use_parallel) {
        message(sprintf("并行模式: 使用 %d 个核心", n_cores))
    } else {
        message("串行模式")
    }
    message("========================================")
    
    # 1. 准备数据
    message("\n[步骤 1/5] 准备数据...")
    
    # 合并 layers（Seurat5）
    if(length(Layers(seurat_obj)) > 1) {
        message("  检测到多个 layers，正在合并...")
        seurat_obj <- JoinLayers(seurat_obj)
    }
    
    # 获取单细胞数据
    counts <- LayerData(seurat_obj, layer = "counts")
    cell_samples <- seurat_obj$sample
    cell_ages <- seurat_obj$age
    
    if(is.null(cell_samples) || is.null(cell_ages)) {
        stop("Seurat object 缺少 'sample' 或 'age' 列")
    }
    
    # 过滤指定年龄的样本（如果指定）
    if (!is.null(exclude_ages)) {
        cells_to_keep <- !cell_ages %in% exclude_ages
        n_cells_removed <- sum(!cells_to_keep)
        
        if (n_cells_removed > 0) {
            # 统计每个被排除年龄的样本数
            excluded_summary <- sapply(exclude_ages, function(age) {
                n_samples <- length(unique(cell_samples[cell_ages == age]))
                n_cells <- sum(cell_ages == age)
                sprintf("%dmo(%d个样本,%d个细胞)", age, n_samples, n_cells)
            })
            
            message(sprintf("  排除指定年龄: %s", paste(excluded_summary, collapse=", ")))
            message(sprintf("  总共去掉 %d 个细胞", n_cells_removed))
            
            # 过滤数据
            counts <- counts[, cells_to_keep, drop = FALSE]
            cell_samples <- cell_samples[cells_to_keep]
            cell_ages <- cell_ages[cells_to_keep]
            
            # 同时更新 seurat_obj 的元数据（用于后续 pseudobulk 计算）
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
        } else {
            message(sprintf("  未检测到指定年龄 (%s) 的样本，无需过滤", 
                            paste(exclude_ages, collapse=", ")))
        }
    }
    
    # 过滤指定样本（如果指定）
    if (!is.null(exclude_samples)) {
        cells_to_keep <- !cell_samples %in% exclude_samples
        n_cells_removed <- sum(!cells_to_keep)
        
        if (n_cells_removed > 0) {
            # 统计每个被排除样本的信息
            excluded_summary <- sapply(exclude_samples, function(sample) {
                if (sample %in% cell_samples) {
                    n_cells <- sum(cell_samples == sample)
                    age <- unique(cell_ages[cell_samples == sample])[1]
                    sprintf("%s(%dmo,%d个细胞)", sample, age, n_cells)
                } else {
                    sprintf("%s(未找到)", sample)
                }
            })
            
            message(sprintf("  排除指定样本: %s", paste(excluded_summary, collapse=", ")))
            message(sprintf("  总共去掉 %d 个细胞", n_cells_removed))
            
            # 过滤数据
            counts <- counts[, cells_to_keep, drop = FALSE]
            cell_samples <- cell_samples[cells_to_keep]
            cell_ages <- cell_ages[cells_to_keep]
            
            # 同时更新 seurat_obj 的元数据（用于后续 pseudobulk 计算）
            seurat_obj <- subset(seurat_obj, cells = colnames(counts))
        } else {
            message(sprintf("  未检测到指定样本 (%s)，无需过滤", 
                            paste(exclude_samples, collapse=", ")))
        }
    }
    
    # 去掉 ambient genes
    if(!is.null(ambient_genes)) {
        ambient_genes <- intersect(ambient_genes, rownames(counts))
        if(length(ambient_genes) > 0) {
            message(sprintf("  去掉 %d 个 ambient genes", length(ambient_genes)))
            counts <- counts[!rownames(counts) %in% ambient_genes, ]
        }
    }
    
    # Depth normalization + Log2 transformation (单细胞级别)
    message("  对单细胞数据进行 depth normalization + log2 transformation...")
    depth <- Matrix::colSums(counts)
    counts_norm <- sweep(counts, 2, depth / 1e4, FUN = "/")
    counts_log2 <- log2(1 + counts_norm)
    
    message(sprintf("  单细胞数: %d, 基因数: %d", ncol(counts_log2), nrow(counts_log2)))
    
    # 2. 计算 pseudobulk（用于特征选择）
    message("\n[步骤 2/5] 计算 pseudobulk（用于特征选择）...")
    pseudobulk_all <- Calculate_Pseudobulk_Normalized(seurat_obj, ambient_genes)
    
    # 获取唯一的 sample 列表
    unique_samples <- unique(cell_samples)
    sample_ages <- sapply(unique_samples, function(s) unique(cell_ages[cell_samples == s])[1])
    n_samples <- length(unique_samples)
    
    message(sprintf("  样本数: %d, 年龄范围: %s", n_samples, paste(unique(sort(sample_ages)), collapse=", ")))
    
    # 3. LOOCV 循环（按 sample，并行处理）
    message("\n[步骤 3/5] 开始 Leave-One-Sample-Out Cross-Validation...")
    
    # 定义单次 LOOCV 处理函数
    process_one_loocv <- function(i) {
        test_sample <- unique_samples[i]
        train_samples <- unique_samples[-i]
        
        # 3.1 Feature Selection: 基于训练集的 pseudobulk 计算基因-年龄相关性
        train_pseudobulk <- pseudobulk_all[, train_samples, drop = FALSE]
        
        # 计算基因-年龄相关性（在所有基因上，抑制输出）
        cor_results <- suppressMessages(Calculate_Gene_Age_Correlation(train_pseudobulk, method = "spearman"))
        
        # 策略1: 显著相关基因 (p < cutoff)
        sig_genes_pval <- cor_results$gene[cor_results$p_value < cor_pval_cutoff]
        
        # 策略2: Top N 相关性最强的基因（按绝对值）
        cor_results_sorted <- cor_results[order(abs(cor_results$correlation), decreasing = TRUE), ]
        sig_genes_top <- head(cor_results_sorted$gene, cor_top_n)
        
        # 合并两种策略
        sig_genes <- unique(c(sig_genes_pval, sig_genes_top))
        
        # 合并 whitelist 基因
        if (!is.null(whitelist_genes)) {
            whitelist_available <- intersect(whitelist_genes, rownames(counts_log2))
            candidate_genes_all <- unique(c(sig_genes, whitelist_available))
        } else {
            candidate_genes_all <- sig_genes
        }
        
        # 检查候选基因数
        if (length(candidate_genes_all) < 10) {
            candidate_genes_all <- head(cor_results_sorted$gene, min(100, nrow(cor_results)))
        }
        
        # 最后过滤低表达基因（所有基因都过滤，包括 whitelist）
        gene_mean_expr <- rowMeans(train_pseudobulk)
        candidate_genes_expr <- gene_mean_expr[candidate_genes_all]
        
        # 过滤低表达基因
        candidate_genes <- candidate_genes_all[!is.na(candidate_genes_expr) & candidate_genes_expr >= min_mean_expr]
        
        # 统计过滤情况
        n_filtered_out_lowexpr <- sum(candidate_genes_expr < min_mean_expr, na.rm = TRUE)
        
        # 检查候选基因是否为空
        if (length(candidate_genes) == 0) {
            stop("过滤后没有候选基因！请降低 min_mean_expr 或 cor_pval_cutoff")
        }
        
        # 3.2 准备训练数据（单细胞级别，仅使用候选基因）
        train_cell_idx <- which(cell_samples %in% train_samples)
        test_cell_idx <- which(cell_samples == test_sample)
        
        # 检查细胞索引
        if (length(train_cell_idx) == 0) {
            stop(sprintf("训练集没有细胞！train_samples: %s", paste(train_samples, collapse=", ")))
        }
        if (length(test_cell_idx) == 0) {
            stop(sprintf("测试样本 %s 没有细胞！", test_sample))
        }
        
        # 确保候选基因在 counts_log2 中存在
        candidate_genes <- intersect(candidate_genes, rownames(counts_log2))
        
        if (length(candidate_genes) == 0) {
            stop("候选基因在单细胞数据中不存在！")
        }
        
        # 安全提取数据
        X_train <- t(as.matrix(counts_log2[candidate_genes, train_cell_idx, drop = FALSE]))
        y_train <- cell_ages[train_cell_idx]
        train_samples_vec <- cell_samples[train_cell_idx]
        
        # Z-score 标准化（基于训练集）
        if (standardize) {
            train_means <- colMeans(X_train)
            train_sds <- apply(X_train, 2, sd)
            train_sds[train_sds == 0] <- 1  # 避免除以 0
            X_train <- scale(X_train, center = train_means, scale = train_sds)
        }
        
        # 计算 weights：平衡 sample 和 age
        weights <- rep(1, length(train_cell_idx))
        
        # 迭代平衡（类似 IPF）
        n_iter <- 10
        for (iter in 1:n_iter) {
            # 平衡 sample：每个 sample 的总权重相等
            sample_weight_sum <- tapply(weights, train_samples_vec, sum)
            target_per_sample <- mean(sample_weight_sum)
            sample_scale <- target_per_sample / sample_weight_sum[train_samples_vec]
            weights <- weights * sample_scale
            
            # 平衡 age：每个 age 的总权重相等
            age_weight_sum <- tapply(weights, y_train, sum)
            target_per_age <- mean(age_weight_sum)
            age_scale <- target_per_age / age_weight_sum[as.character(y_train)]
            weights <- weights * age_scale
        }
        
        # 归一化 weights（总和 = 样本数）
        weights <- weights * length(weights) / sum(weights)
        
        # 3.3 Nested CV: 调优 alpha（内层 CV，在单细胞上）
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
                         standardize = FALSE)  # 已经手动标准化
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
        
        # 3.4 用最佳 alpha 训练最终模型（在全部训练细胞上）
        final_model <- glmnet(X_train, y_train, 
                             weights = weights,
                             alpha = best_alpha, 
                             family = "gaussian",
                             standardize = FALSE)
        
        # 3.5 预测训练集（用于校准）
        if (!is.null(best_lambda)) {
            train_pred <- predict(final_model, newx = X_train, s = best_lambda)
        } else {
            train_pred <- predict(final_model, newx = X_train, s = min(final_model$lambda))
        }
        
        # 3.6 在训练集上拟合校准模型
        train_pred_df <- data.frame(
            actual = y_train,
            predicted = as.numeric(train_pred)
        )
        # 按 sample 聚合到样本级别
        train_pred_sample <- aggregate(cbind(actual, predicted) ~ train_samples_vec, 
                                       data = cbind(train_pred_df, train_samples_vec), 
                                       FUN = mean)
        colnames(train_pred_sample)[1] <- "sample"
        
        # 拟合校准模型
        calib_model_i <- lm(actual ~ predicted, data = train_pred_sample)
        calib_intercept <- coef(calib_model_i)[1]
        calib_slope <- coef(calib_model_i)[2]
        
        # 3.7 预测测试样本的所有细胞
        X_test <- t(as.matrix(counts_log2[candidate_genes, test_cell_idx, drop = FALSE]))
        
        # 检查测试数据维度
        if (ncol(X_test) != ncol(X_train)) {
            stop(sprintf("测试集特征数 (%d) 与训练集 (%d) 不匹配！", ncol(X_test), ncol(X_train)))
        }
        
        # 使用相同的标准化参数
        if (standardize) {
            X_test <- scale(X_test, center = train_means, scale = train_sds)
        }
        
        # 使用最佳 lambda 预测
        if (!is.null(best_lambda)) {
            pred <- predict(final_model, newx = X_test, s = best_lambda)
        } else {
            pred <- predict(final_model, newx = X_test, s = min(final_model$lambda))
        }
        
        # 3.8 应用校准到测试样本
        pred_calibrated <- calib_intercept + calib_slope * as.numeric(pred)
        
        # 返回结果
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
    
           # 并行或串行执行 LOOCV
           if (use_parallel) {
               message(sprintf("并行模式: 使用 %d 个核心，正在初始化...\n", n_cores))
               
               # 使用 pbmcapply::pbmclapply（带进度条）
               loocv_results <- tryCatch({
                   pbmcapply::pbmclapply(1:n_samples, process_one_loocv, 
                                        mc.cores = n_cores, 
                                        mc.set.seed = TRUE,
                                        mc.preschedule = FALSE,
                                        ignore.interactive = TRUE)
               }, error = function(e) {
                   message("\n警告: pbmcapply 执行失败，尝试使用 parallel::mclapply...")
                   message(sprintf("错误信息: %s", e$message))
                   
                   # 备选方案：使用 parallel::mclapply（无进度条）
                   tryCatch({
                       message("使用 parallel::mclapply（无进度条）...\n")
                       parallel::mclapply(1:n_samples, process_one_loocv, 
                                         mc.cores = n_cores, 
                                         mc.set.seed = TRUE,
                                         mc.preschedule = FALSE)
                   }, error = function(e2) {
                       message("\n警告: 并行执行失败，降级为串行模式")
                       message(sprintf("错误信息: %s", e2$message))
                       message("正在重新以串行模式运行...\n")
                       
                       # 最终降级为串行
                       lapply(1:n_samples, function(i) {
                           message(sprintf(">>> LOOCV Round %d/%d <<<", i, n_samples))
                           result <- process_one_loocv(i)
                           message(sprintf("    完成: %s | genes=%d | pred=%.2f, actual=%.2f", 
                                           result$test_sample, result$n_genes, 
                                           mean(result$predictions), result$actual_age))
                           return(result)
                       })
                   })
               })
           } else {
               message("串行模式: 将显示每轮详细进度...\n")
               loocv_results <- lapply(1:n_samples, function(i) {
                   message(sprintf(">>> LOOCV Round %d/%d <<<", i, n_samples))
                   result <- process_one_loocv(i)
                   message(sprintf("    完成: %s | genes=%d | pred=%.2f→%.2f, actual=%.2f", 
                                   result$test_sample, result$n_genes, 
                                   mean(result$predictions), mean(result$predictions_calibrated),
                                   result$actual_age))
                   return(result)
               })
           }
    
    # 整合结果
    message("\n[步骤 3.5/5] 整合并行结果...")
    
    # 检查结果有效性
    valid_results <- sapply(loocv_results, function(x) is.list(x) && !is.null(x$test_sample))
    if (!all(valid_results)) {
        invalid_idx <- which(!valid_results)
        stop(sprintf("LOOCV Round %s 返回了无效结果，请检查错误信息", 
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
        
        # 打印样本级汇总
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
    
    # 4. 汇总样本级结果
    message("\n[步骤 4/5] 汇总样本级预测...")
    cell_pred_complete <- cell_predictions[complete.cases(cell_predictions), ]
    
    # 汇总原始预测和校准后预测
    sample_predictions <- aggregate(cbind(predicted, predicted_calibrated) ~ sample + actual, 
                                   data = cell_pred_complete, 
                                   FUN = mean)
    
    # 5. 计算整体性能
    message("\n[步骤 5/5] 计算整体预测性能...")
    
    # 样本级性能 - 原始预测
    rmse_sample <- sqrt(mean((sample_predictions$actual - sample_predictions$predicted)^2))
    mae_sample <- mean(abs(sample_predictions$actual - sample_predictions$predicted))
    pcc_sample <- cor(sample_predictions$actual, sample_predictions$predicted, method = "pearson")
    r2_sample <- 1 - sum((sample_predictions$actual - sample_predictions$predicted)^2) / 
                     sum((sample_predictions$actual - mean(sample_predictions$actual))^2)
    
    # 样本级性能 - 校准后预测
    rmse_sample_calib <- sqrt(mean((sample_predictions$actual - sample_predictions$predicted_calibrated)^2))
    mae_sample_calib <- mean(abs(sample_predictions$actual - sample_predictions$predicted_calibrated))
    pcc_sample_calib <- cor(sample_predictions$actual, sample_predictions$predicted_calibrated, method = "pearson")
    r2_sample_calib <- 1 - sum((sample_predictions$actual - sample_predictions$predicted_calibrated)^2) / 
                           sum((sample_predictions$actual - mean(sample_predictions$actual))^2)
    
    message("  样本级性能（原始）:")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_sample, mae_sample, pcc_sample, r2_sample))
    
    message("  样本级性能（校准后）:")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_sample_calib, mae_sample_calib, pcc_sample_calib, r2_sample_calib))
    
    # 细胞级性能 - 原始预测
    rmse_cell <- sqrt(mean((cell_pred_complete$actual - cell_pred_complete$predicted)^2))
    mae_cell <- mean(abs(cell_pred_complete$actual - cell_pred_complete$predicted))
    pcc_cell <- cor(cell_pred_complete$actual, cell_pred_complete$predicted, method = "pearson")
    r2_cell <- 1 - sum((cell_pred_complete$actual - cell_pred_complete$predicted)^2) / 
                   sum((cell_pred_complete$actual - mean(cell_pred_complete$actual))^2)
    
    # 细胞级性能 - 校准后预测
    rmse_cell_calib <- sqrt(mean((cell_pred_complete$actual - cell_pred_complete$predicted_calibrated)^2))
    mae_cell_calib <- mean(abs(cell_pred_complete$actual - cell_pred_complete$predicted_calibrated))
    pcc_cell_calib <- cor(cell_pred_complete$actual, cell_pred_complete$predicted_calibrated, method = "pearson")
    r2_cell_calib <- 1 - sum((cell_pred_complete$actual - cell_pred_complete$predicted_calibrated)^2) / 
                         sum((cell_pred_complete$actual - mean(cell_pred_complete$actual))^2)
    
    message("  细胞级性能（原始）:")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_cell, mae_cell, pcc_cell, r2_cell))
    
    message("  细胞级性能（校准后）:")
    message(sprintf("    RMSE: %.4f, MAE: %.4f, PCC: %.4f, R2: %.4f", 
                   rmse_cell_calib, mae_cell_calib, pcc_cell_calib, r2_cell_calib))
    
    message("\n========================================")
    message("完成！")
    message("========================================\n")
    
    # 6. 额外统计信息
    message("\n[额外信息] 特征选择统计:")
    n_genes_selected <- sapply(selected_genes_list, length)
    message(sprintf("  每次选择的基因数: 平均=%.1f, 中位数=%.0f, 范围=[%d, %d]",
                   mean(n_genes_selected), median(n_genes_selected),
                   min(n_genes_selected), max(n_genes_selected)))
    
    n_genes_filtered <- sapply(loocv_results, function(x) x$n_genes_filtered_lowexpr)
    message(sprintf("  低表达过滤基因数 (mean_expr < %.2f): 平均=%.1f, 范围=[%d, %d]",
                   min_mean_expr,
                   mean(n_genes_filtered), 
                   min(n_genes_filtered), max(n_genes_filtered)))
    
    # 找出最常被选择的基因
    all_selected <- unlist(selected_genes_list)
    gene_freq <- sort(table(all_selected), decreasing = TRUE)
    message(sprintf("  最常被选择的前5个基因: %s", 
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
# 使用示例
################################################################################
# 
# # 加载此文件
# source("Mouse_AgingClock_Functions.R")
# 
# # 示例 1: 自动检测核心数（推荐）
# results <- Nested_FeatureSelection_LOOCV(
#     seurat_obj = your_seurat_object,
#     whitelist_genes = your_whitelist,      # 可选
#     ambient_genes = ambient_gene_list,     # 可选
#     cor_pval_cutoff = 0.05,                # p-value 阈值
#     cor_top_n = 500,                       # Top N 相关基因
#     min_mean_expr = 0.1,                   # 最低平均表达量
#     alphas = seq(0.1, 0.9, by = 0.2),      # Elastic Net alpha 范围
#     nfolds_inner = 5,                      # 内层 CV 折数
#     standardize = TRUE,                    # Z-score 标准化
#     n_cores = NULL                         # NULL=自动检测
# )
# 
# # 示例 2: 手动指定核心数
# results <- Nested_FeatureSelection_LOOCV(
#     seurat_obj = your_seurat_object,
#     n_cores = 8                            # 使用 8 个核心
# )
# 
# # 示例 3: 串行模式（单核）
# results <- Nested_FeatureSelection_LOOCV(
#     seurat_obj = your_seurat_object,
#     n_cores = 1                            # 串行模式
# )
# 
# # 查看结果
# print(results$metrics)
# print(results$sample_predictions)
# print(results$calibration_params)
# 
# # 保存结果
# write.csv(results$cell_predictions, "cell_predictions.csv", row.names = FALSE)
# write.csv(results$sample_predictions, "sample_predictions.csv", row.names = FALSE)
# write.csv(results$calibration_params, "calibration_params.csv", row.names = FALSE)
# saveRDS(results, "aging_clock_results.rds")
# 
# # 核心数设置说明：
# # - n_cores = NULL: 自动检测，使用 (总核心数 - 1)
# # - n_cores = 1: 串行模式，内存占用最小
# # - n_cores = 4: 使用 4 个核心
# # - n_cores = 8: 使用 8 个核心（服务器推荐）
# # - n_cores = 16: 使用 16 个核心（高性能服务器）
# 
# # 参数调优建议：
# # 1. 如果性能不佳（R2 < 0.7），尝试：
# #    - 增加 cor_top_n (如 1000)
# #    - 降低 min_mean_expr (如 0.05)
# #    - 增加 alphas 密度 (如 seq(0.1, 0.9, by=0.1))
# # 2. 如果过拟合（训练好但测试差），尝试：
# #    - 减少 cor_top_n (如 300)
# #    - 增加 min_mean_expr (如 0.5)
# # 3. 如果运行太慢，尝试：
# #    - 增加 n_cores
# #    - 减少 cor_top_n
# #    - 减少 alphas 数量
# 
################################################################################
