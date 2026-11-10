#!/usr/bin/env Rscript
########################################################################
# Zebrafish Injury — Aging Clock Prediction (2026)
# 用 2026 final aging clock 模型预测 zebrafish NMDA injury 数据的衰老评分
# celltypes: MG, Rod, Cone
# 流程（与 final model 训练完全一致，无 calibration）:
#   raw counts -> 按 sample KNN smoothing -> 去 ambient -> depth norm(1e4)+log2
#   -> 按 genes_used 对齐 -> 用训练存的 means/sds 做 Z-score -> predict
# HPC 用法: Rscript Run_Zebrafish_Injury_AgingClock_2026.R
########################################################################

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(FNN)
    library(glmnet)     # 必须: 否则 predict 找不到 predict.glmnet
    library(ggplot2)
})

## 路径配置 ------------------------------------------------------------
injury_dir <- "/projects/hmz-aging/Injury_Datasets"
model_dir  <- "/projects/hmz-aging/Zebrafish_New_Datasets"   # 模型: <CT>/<CT>_final_aging_clock_model_2026.rds
clock_celltypes <- c("MG", "Rod", "Cone")
control_sample  <- "Control_2"


########################################################################
# 1) 读入 zebrafish injury 数据，提取 NMDA 样本，整理 celltype (MG/Rod/Cone)
########################################################################
setwd(injury_dir)
Injury_merge_seurat_all_2023Dec <- readRDS("Injury_merge_seurat_all_2023Dec(1)")

all_samples <- names(table(Injury_merge_seurat_all_2023Dec$sample))
# NMDA 相关样本（第 1 与 8:13 个；含 Control_2）
NMDA_samples_seurat <- Injury_merge_seurat_all_2023Dec[
    , which(Injury_merge_seurat_all_2023Dec$sample %in% all_samples[c(1, 8:13)])]
message("NMDA samples:"); print(table(NMDA_samples_seurat$sample))

# RestMG / ActivateMG 合并为 MG
NMDA_samples_seurat$celltype2 <- NMDA_samples_seurat$new_celltype
NMDA_samples_seurat$celltype2[NMDA_samples_seurat$celltype2 %in% c("RestMG")] <- "MG"

NMDA_samples_seurat_cl <- NMDA_samples_seurat[, NMDA_samples_seurat$celltype2 %in% clock_celltypes]
NMDA_samples_seurat_cl$celltype <- NMDA_samples_seurat_cl$celltype2
NMDA_samples_seurat_cl_sp <- SplitObject(NMDA_samples_seurat_cl, split.by = "celltype")


########################################################################
# 2) 按 nFeature_RNA 的 10%–90% 分位过滤细胞
########################################################################
Clean_the_cells_by_counts <- function(obj_list) {
    out <- list()
    for (i in seq_along(obj_list)) {
        tmp <- obj_list[[i]]
        lo <- quantile(tmp$nFeature_RNA, 0.1); hi <- quantile(tmp$nFeature_RNA, 0.9)
        tmp <- tmp[, which(tmp$nFeature_RNA > lo & tmp$nFeature_RNA < hi)]
        message(sprintf("  %s: 保留 %d cells (nFeature %g–%g)", names(obj_list)[i], ncol(tmp), lo, hi))
        out[[names(obj_list)[i]]] <- tmp
    }
    out
}
NMDA_sel <- Clean_the_cells_by_counts(NMDA_samples_seurat_cl_sp)


########################################################################
# 3) 基因名整理：取 "~~" 后的 gene symbol，合并重复 symbol 的 raw counts
########################################################################
for (ct in names(NMDA_sel)) {
    seu <- NMDA_sel[[ct]]
    counts <- GetAssayData(seu, layer = "counts")
    old_names <- rownames(counts)
    new_names <- sapply(old_names, function(x) {
        parts <- strsplit(x, "~~")[[1]]; if (length(parts) >= 2) parts[2] else x
    }, USE.NAMES = FALSE)
    dup_symbols <- unique(new_names[duplicated(new_names)])
    message(sprintf("  %s: 总基因 %d, 重复 symbol %d", ct, length(new_names), length(dup_symbols)))
    if (length(dup_symbols) > 0) {
        uniq_idx <- which(!new_names %in% dup_symbols)
        uniq_counts <- counts[uniq_idx, , drop = FALSE]; rownames(uniq_counts) <- new_names[uniq_idx]
        dup_list <- lapply(dup_symbols, function(sym) Matrix::colSums(counts[which(new_names == sym), , drop = FALSE]))
        dup_mat <- as(do.call(rbind, dup_list), "dgCMatrix"); rownames(dup_mat) <- dup_symbols
        final_counts <- rbind(uniq_counts, dup_mat)
    } else {
        rownames(counts) <- new_names; final_counts <- counts
    }
    NMDA_sel[[ct]] <- CreateSeuratObject(counts = final_counts, meta.data = seu@meta.data)
}


########################################################################
# 4) KNN smoothing（按 sample；与 aging clock 训练预处理一致）
########################################################################
NormalizeAndPCA <- function(seu, scale.factor = 1e4, nfeatures = 2000, npcs = 30) {
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
    max_npcs <- min(nrow(seu), ncol(seu)) - 1
    if (npcs > max_npcs) npcs <- max(1, max_npcs)
    if (npcs < 1) stop(sprintf("细胞数太少 (n_cells=%d)，无法 PCA", ncol(seu)))
    RunPCA(seu, features = VariableFeatures(seu), npcs = npcs, verbose = FALSE)
}

PseudoBulkByKNN_fast <- function(seu, k = 15, dims = 1:30, assay = "RNA", reduction = "pca") {
    emb <- Embeddings(seu, reduction)[, dims, drop = FALSE]
    nn_idx <- get.knn(emb, k = k)$nn.index
    ncells <- nrow(nn_idx)
    rows <- as.integer(c(t(nn_idx), seq_len(ncells)))
    cols <- rep(seq_len(ncells), each = k + 1)
    A <- sparseMatrix(i = rows, j = cols, x = 1, dims = c(ncells, ncells))
    pseudo_counts <- GetAssayData(seu, assay = assay, layer = "counts") %*% A
    CreateSeuratObject(counts = pseudo_counts, meta.data = seu@meta.data, assay = assay)
}

# 对 celltype list 中每个 celltype 按 sample 列 KNN smooth（含动态调参）
smooth_injury_list <- function(seurat_list, sample_col = "sample", species_label = "Zebrafish") {
    smoothed <- list()
    for (ct in names(seurat_list)) {
        message(sprintf("\n========== %s - %s ==========", species_label, ct))
        tmp_sp <- SplitObject(seurat_list[[ct]], split.by = sample_col)
        nc <- sapply(tmp_sp, ncol)
        if (any(nc < 2)) { message("  去掉 <2 细胞样本: ", paste(names(nc)[nc < 2], collapse = ", ")); tmp_sp <- tmp_sp[nc >= 2] }
        if (length(tmp_sp) == 0) { message("  ", ct, ": 所有样本被过滤，跳过"); next }
        # 每个 celltype 的 K / npcs（与训练一致）
        if (ct %in% c("Rod", "Cone", "RGC", "AC", "HC", "BC")) { k <- 15; npcs <- 25 }
        else if (ct == "MG") { k <- 5; npcs <- 5 }
        else { k <- 2; npcs <- 5 }
        min_ncells <- min(sapply(tmp_sp, ncol))
        message(sprintf("  min cells/sample=%d, K=%d, npcs=%d", min_ncells, k, npcs))
        if (min_ncells < 5) {
            message("  最小样本细胞数 <5，跳过 smooth 直接 merge")
            smoothed[[ct]] <- if (length(tmp_sp) == 1) tmp_sp[[1]] else merge(tmp_sp[[1]], tmp_sp[-1], add.cell.ids = names(tmp_sp))
            next
        }
        if (npcs >= min_ncells) npcs <- max(1, min_ncells - 1)
        if (k    >= min_ncells) k    <- max(1, min_ncells - 1)
        tmp_pca    <- sapply(tmp_sp, NormalizeAndPCA, npcs = npcs)
        tmp_smooth <- sapply(tmp_pca, PseudoBulkByKNN_fast, k = k, dims = 1:npcs)
        smoothed[[ct]] <- if (length(tmp_smooth) == 1) tmp_smooth[[1]] else merge(tmp_smooth[[1]], tmp_smooth[-1], add.cell.ids = names(tmp_smooth))
        message(sprintf("  完成: %d pseudo-cells", ncol(smoothed[[ct]])))
    }
    smoothed
}

Zebrafish_injury_smoothed <- smooth_injury_list(NMDA_sel, sample_col = "sample", species_label = "Zebrafish")
setwd(injury_dir)
saveRDS(Zebrafish_injury_smoothed, file = "Zebrafish_injury_smoothed.rds")
message("已保存: Zebrafish_injury_smoothed.rds")


########################################################################
# 5) 载入 ambient 基因 + 2026 final 模型，逐 celltype 预测（无 calibration）
########################################################################
setwd(model_dir)
load("Zebrafish_ambient_genes_Top50.rda")   # -> Z_ambient_genes_Top50

# 预测流程 = 训练流程；2026 无校准
predict_final_clock <- function(final_model_obj, seu_ct, ambient_genes = NULL) {
    counts <- GetAssayData(seu_ct, assay = "RNA", layer = "counts")
    if (!is.null(ambient_genes)) {
        ag <- intersect(ambient_genes, rownames(counts))
        if (length(ag) > 0) counts <- counts[!rownames(counts) %in% ag, , drop = FALSE]
    }
    depth <- Matrix::colSums(counts); depth[depth == 0] <- 1
    counts_norm <- sweep(counts, 2, depth / 1e4, FUN = "/")
    counts_log2 <- log2(1 + counts_norm)
    genes   <- final_model_obj$genes_used
    present <- intersect(genes, rownames(counts_log2))
    missing <- setdiff(genes, rownames(counts_log2))
    mat <- counts_log2[present, , drop = FALSE]
    if (length(missing) > 0) {
        zero_mat <- Matrix::Matrix(0, nrow = length(missing), ncol = ncol(counts_log2), sparse = TRUE)
        rownames(zero_mat) <- missing
        mat <- rbind(mat, zero_mat)
    }
    mat <- mat[genes, , drop = FALSE]
    X <- t(as.matrix(mat))
    if (isTRUE(final_model_obj$standardization$enabled)) {
        X <- scale(X, center = final_model_obj$standardization$means,
                   scale = final_model_obj$standardization$sds)
    }
    pred <- as.numeric(predict(final_model_obj$model, newx = X, s = final_model_obj$lambda))
    names(pred) <- colnames(mat)
    list(pred = pred, n_missing_genes = length(missing), n_genes = length(genes))
}

Zebrafish_injury_predictions <- list()
for (ct in clock_celltypes) {
    if (!ct %in% names(Zebrafish_injury_smoothed)) { message(ct, ": 无 smoothed 数据，跳过"); next }
    model_file <- file.path(model_dir, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    if (!file.exists(model_file)) { message(ct, ": 2026 模型不存在，跳过 -> ", model_file); next }
    final_model <- readRDS(model_file)
    seu <- Zebrafish_injury_smoothed[[ct]]
    if (length(Layers(seu)) > 1) seu <- JoinLayers(seu)
    res <- predict_final_clock(final_model, seu, ambient_genes = Z_ambient_genes_Top50)
    cell_pred <- data.frame(
        cell                 = names(res$pred),
        sample               = as.character(seu@meta.data[["sample"]]),
        predicted            = res$pred,
        predicted_calibrated = res$pred,   # 2026 无校准: 别名
        celltype             = ct,
        stringsAsFactors = FALSE
    )
    sample_pred <- aggregate(cbind(predicted, predicted_calibrated) ~ sample, data = cell_pred, FUN = mean)
    sample_pred$celltype <- ct
    Zebrafish_injury_predictions[[ct]] <- list(cell_predictions = cell_pred, sample_predictions = sample_pred)
    message(sprintf("%s: %d cells, %d samples (缺失基因 %d/%d)",
                    ct, nrow(cell_pred), nrow(sample_pred), res$n_missing_genes, res$n_genes))
    print(sample_pred)
}

setwd(injury_dir)
saveRDS(Zebrafish_injury_predictions, file = "Zebrafish_injury_predictions_2026.rds")
message("已保存: Zebrafish_injury_predictions_2026.rds")


########################################################################
# 6) Delta aging：以 Control_2 的均值归一化（injury vs control）
########################################################################
compute_delta_aging <- function(predictions_list, control_sample, species_label = "Zebrafish") {
    all_delta <- list()
    for (ct in names(predictions_list)) {
        cp <- predictions_list[[ct]]$cell_predictions
        is_ctrl <- cp$sample == control_sample
        ctrl_mean <- mean(cp$predicted[is_ctrl])
        message(sprintf("  %s %s: control mean=%.2f (%d cells), injury %d cells",
                        species_label, ct, ctrl_mean, sum(is_ctrl), sum(!is_ctrl)))
        cp$delta_aging <- cp$predicted - ctrl_mean
        cp$group <- ifelse(is_ctrl, "Control", "Injury")
        all_delta[[ct]] <- cp
    }
    do.call(rbind, all_delta)
}

message("\nZebrafish sample names:")
for (ct in names(Zebrafish_injury_predictions))
    message(sprintf("  %s: %s", ct, paste(unique(Zebrafish_injury_predictions[[ct]]$cell_predictions$sample), collapse = ", ")))

Zebrafish_delta <- compute_delta_aging(Zebrafish_injury_predictions, control_sample = control_sample)
# Control_2 即 0 小时对照，展示时重命名为 NMDA_0hr（delta 已相对 Control_2 计算）
Zebrafish_delta$sample[Zebrafish_delta$sample == "Control_2"] <- "NMDA_0hr"
saveRDS(Zebrafish_delta, file = "Zebrafish_injury_delta_aging_2026.rds")

# sample 级汇总
zf_sample_delta <- aggregate(cbind(predicted, delta_aging) ~ sample + group + celltype,
                             data = Zebrafish_delta, FUN = mean)
zf_sample_delta <- zf_sample_delta[order(zf_sample_delta$celltype, zf_sample_delta$group), ]
message("\nZebrafish - Sample Level Delta Aging Score"); print(zf_sample_delta, row.names = FALSE)
write.csv(zf_sample_delta, "Zebrafish_injury_sample_delta_aging_2026.csv", row.names = FALSE)


########################################################################
# 7) 画图：delta by celltype + 时间点分面
########################################################################
# (a) 每个 celltype: Control vs Injury 的 delta boxplot
plot_delta_boxplot <- function(delta_df, species_label = "Zebrafish") {
    ct_order <- clock_celltypes[clock_celltypes %in% delta_df$celltype]
    delta_df$celltype <- factor(delta_df$celltype, levels = ct_order)
    ggplot(delta_df, aes(x = celltype, y = delta_aging, fill = group)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        scale_fill_manual(values = c("Control" = "#4DBBD5", "Injury" = "#E64B35")) +
        labs(title = sprintf("%s: Injury-Induced Aging Acceleration", species_label),
             subtitle = "Delta = predicted age - control mean (per cell type)",
             x = "Cell Type", y = "Delta Aging Score", fill = "Group") +
        theme_bw() +
        theme(axis.text.x = element_text(size = 12),
              plot.title = element_text(size = 13, face = "bold"),
              plot.subtitle = element_text(size = 10, color = "grey40"))
}
ggsave("Zebrafish_Injury_Delta_Aging_Boxplot_2026.png",
       plot_delta_boxplot(Zebrafish_delta), width = 5, height = 4, dpi = 300)
message("已保存: Zebrafish_Injury_Delta_Aging_Boxplot_2026.png")

# (b) 时间点分面：每个 celltype 一行，x = injury 时间点（Nature/NPG 配色）
# Nature(NPG) 调色板
npg_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
              "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

plot_timepoint_facet <- function(delta_df, sample_order, ct_order, species_label = "Zebrafish") {
    delta_df <- delta_df[delta_df$sample %in% sample_order & delta_df$celltype %in% ct_order, ]
    delta_df$sample   <- factor(delta_df$sample, levels = sample_order)
    delta_df$celltype <- factor(delta_df$celltype, levels = ct_order)
    fill_cols <- setNames(rep(npg_cols, length.out = length(sample_order)), sample_order)
    ggplot(delta_df, aes(x = sample, y = delta_aging, fill = sample)) +
        geom_boxplot(outlier.size = 0.3, alpha = 0.85) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        facet_wrap(~ celltype, ncol = 1, scales = "free_y") +
        scale_fill_manual(values = fill_cols) +
        labs(x = NULL, y = "Relative aging changes to control (months)") +
        theme_classic(base_size = 15) +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
              axis.text.x  = element_text(angle = 45, hjust = 1, size = 15, color = "black"),
              axis.text.y  = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 18, face = "bold"),
              strip.text   = element_text(size = 17, face = "bold"),
              strip.background = element_rect(color = "black", fill = "white"),
              legend.position = "none")
}
# x 轴顺序（按图片给定）
zf_sample_order <- c("NMDA_0hr", "NMDA_36hr", "NMDA_54hr", "NMDA_72hr", "NMDA_96hr", "NMDA_7D", "NMDA_14D")
zf_ct_order     <- c("MG", "Rod", "Cone")
ggsave("Zebrafish_Injury_Timepoint_Boxplot_2026.png",
       plot_timepoint_facet(Zebrafish_delta, zf_sample_order, zf_ct_order),
       width = 4.5, height = length(zf_ct_order) * 2.6 + 1, dpi = 300)
message("已保存: Zebrafish_Injury_Timepoint_Boxplot_2026.png")

message("\n=== Zebrafish injury aging clock (2026) 全部完成 ===")


#######
####### 接下来参考 Aging_Clocks_injury2026 ######
####### 计算和 aging 相关的 injury genes #########
#######
########################################################################
# 8) aging 相关的 injury genes:
#    (a) 每个时间点 vs Control_2 找 injury DEGs (MG/Rod/Cone)  —— 用未 smooth 的原始 counts
#    (b) 合并 aging clock 基因权重
#    (c) injury DEG 且是 clock gene = aging 相关 injury genes
########################################################################
library(Seurat)
library(Matrix)
setwd(injury_dir)
# DEGs 在未 smooth 的原始对象 NMDA_sel 上计算（细胞相互独立，统计更严谨）

# (a) injury DEGs：每个时间点 vs Control_2（p_val_adj<0.01, |log2FC|>0.2）
find_injury_degs <- function(seurat_list, control_sample) {
    all_degs <- list()
    for (ct in names(seurat_list)) {
        message(sprintf("\n========== Zebrafish - %s ==========", ct))
        seu <- seurat_list[[ct]]
        # 仅 v5 Assay5 且多 layer 时才 JoinLayers（v4 Assay 不支持）
        if (inherits(seu[[DefaultAssay(seu)]], "Assay5") && length(Layers(seu)) > 1)
            seu <- JoinLayers(seu)
        seu <- NormalizeData(seu, verbose = FALSE)
        Idents(seu) <- seu$sample
        injury_samples <- setdiff(unique(as.character(seu$sample)), control_sample)
        ct_degs <- list()
        for (tp in injury_samples) {
            markers <- tryCatch(
                FindMarkers(seu, ident.1 = tp, ident.2 = control_sample,
                            min.pct = 0.05, logfc.threshold = 0.2, test.use = "wilcox"),
                error = function(e) { message("    FindMarkers 出错: ", e$message); NULL })
            if (is.null(markers) || nrow(markers) == 0) next
            sig <- markers[markers$p_val_adj < 0.01, , drop = FALSE]
            if (nrow(sig) == 0) next
            # 归一基因名到纯 symbol（去掉 "ENSDARG...~~" 前缀），以便和 clock 权重对齐
            sig$gene <- sub(".*~~", "", rownames(sig)); sig$timepoint <- tp; sig$celltype <- ct
            sig$direction <- ifelse(sig$avg_log2FC > 0, "up", "down")
            ct_degs[[tp]] <- sig
            message(sprintf("  %s vs %s: %d DEGs", tp, control_sample, nrow(sig)))
        }
        if (length(ct_degs) > 0) { all_degs[[ct]] <- do.call(rbind, ct_degs); rownames(all_degs[[ct]]) <- NULL }
    }
    all_degs
}

Zebrafish_injury_degs <- find_injury_degs(NMDA_samples_seurat_cl_sp, control_sample = control_sample)
# 把 DEG 基因名归一到纯 symbol（去掉 "ENSDARG...~~" 前缀），以便和 clock 权重对齐
for (ct in names(Zebrafish_injury_degs))
    Zebrafish_injury_degs[[ct]]$gene <- sub(".*~~", "", Zebrafish_injury_degs[[ct]]$gene)
saveRDS(Zebrafish_injury_degs, file = "Zebrafish_injury_degs_list_2026.rds")

# (b) 只做 MG：载入之前预测用的 MG final aging clock 权重（从 final 模型 $coefs 提取）
zf_gene_weights <- list()
{
    ct <- "MG"
    mf <- file.path(model_dir, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    if (!file.exists(mf)) stop("缺少 MG final 模型: ", mf)
    cf <- readRDS(mf)$coefs                                   # gene, coef
    cf <- cf[cf$gene != "(Intercept)" & cf$coef != 0, , drop = FALSE]
    df <- data.frame(gene              = cf$gene,
                     clock_coefficient = cf$coef,
                     clock_abs_coef    = abs(cf$coef),
                     stringsAsFactors  = FALSE)
    df <- df[order(-df$clock_abs_coef), , drop = FALSE]
    df$clock_rank <- seq_len(nrow(df))
    zf_gene_weights[[ct]] <- df
    message(sprintf("  %s: %d clock genes", ct, nrow(df)))
}

# (c) DEGs 合并 clock weights，标注 is_clock_gene
merge_degs_with_weights <- function(degs_list, weights_list) {
    out <- list()
    for (ct in names(degs_list)) {
        deg_df <- degs_list[[ct]]; wt <- weights_list[[ct]]
        if (is.null(wt)) {
            deg_df$clock_coefficient <- NA; deg_df$clock_abs_coef <- NA
            deg_df$clock_rank <- NA; deg_df$is_clock_gene <- FALSE
        } else {
            wt_tmp <- wt
            if ("coefficient" %in% colnames(wt_tmp) && !"clock_coefficient" %in% colnames(wt_tmp)) {
                colnames(wt_tmp)[colnames(wt_tmp) == "coefficient"]     <- "clock_coefficient"
                colnames(wt_tmp)[colnames(wt_tmp) == "abs_coefficient"] <- "clock_abs_coef"
                colnames(wt_tmp)[colnames(wt_tmp) == "rank"]            <- "clock_rank"
            }
            mc <- intersect(c("gene", "clock_coefficient", "clock_abs_coef", "clock_rank"), colnames(wt_tmp))
            m <- merge(deg_df, wt_tmp[, mc, drop = FALSE], by = "gene", all.x = TRUE)
            m$is_clock_gene <- !is.na(m$clock_coefficient)
            deg_df <- m
        }
        out[[ct]] <- deg_df
        message(sprintf("  %s: %d DEGs, clock overlap=%d", ct, nrow(deg_df), sum(deg_df$is_clock_gene, na.rm = TRUE)))
    }
    out
}



Zebrafish_injury_degs_table <- merge_degs_with_weights(Zebrafish_injury_degs, zf_gene_weights)
zf_deg_table <- do.call(rbind, Zebrafish_injury_degs_table); rownames(zf_deg_table) <- NULL
write.csv(zf_deg_table, "Zebrafish_injury_DEGs_with_clock_weights_2026.csv", row.names = FALSE)
saveRDS(Zebrafish_injury_degs_table, file = "Zebrafish_injury_degs_table_2026.rds")

# aging 相关 injury genes = injury DEG 且是 clock gene
zf_aging_injury_genes <- zf_deg_table[zf_deg_table$is_clock_gene %in% TRUE, , drop = FALSE]
write.csv(zf_aging_injury_genes, "Zebrafish_aging_related_injury_genes_2026.csv", row.names = FALSE)
message(sprintf("\naging 相关 injury genes: %d 行, %d unique genes",
                nrow(zf_aging_injury_genes), length(unique(zf_aging_injury_genes$gene))))




##############
############## 接下来，我们要把 这些 zf_aging_injury_genes 按照 weights 分为 两类 #####
############## 一类是 大于 0 的，一类是 小于 0 的 ###
##############
############## 用 kmeans 聚类，聚类矩阵 是个 拼接矩阵，左边是 z-score 后的 每个 时间点的 pseudobulk，右边是 aging 时间点的 pseudobulk ###
############## 先参考下 之前的 代码 构建 矩阵 ###
########################################################################
# 9) 构建拼接矩阵 (injury z-score | aging/training z-score) + 按权重正负分两类 + kmeans
########################################################################
library(Matrix)
ct_use <- "MG"      # 聚类/heatmap 聚焦的 celltype

# pseudobulk: 每个 sample 求 raw counts 之和 -> depth norm 1e4 -> log2(1+x)
Make_Pseudobulk_Matrix <- function(seurat_obj, ambient_genes = NULL, sample_col = "sample") {
    if (inherits(seurat_obj[[DefaultAssay(seurat_obj)]], "Assay5") && length(Layers(seurat_obj)) > 1)
        seurat_obj <- JoinLayers(seurat_obj)
    counts  <- GetAssayData(seurat_obj, layer = "counts")
    samples <- as.character(seurat_obj@meta.data[[sample_col]])
    if (!is.null(ambient_genes)) {
        amb <- intersect(ambient_genes, rownames(counts))
        if (length(amb) > 0) counts <- counts[!rownames(counts) %in% amb, ]
    }
    sn <- unique(samples)
    pb <- sapply(sn, function(s) Matrix::rowSums(counts[, samples == s, drop = FALSE]))
    log2(1 + sweep(pb, 2, colSums(pb) / 1e4, FUN = "/"))
}
zscore_rows <- function(m) { r <- t(scale(t(m))); r[is.na(r)] <- 0; r }

# --- (a) injury pseudobulk（smoothed；列按时间点排序）---
setwd(model_dir)
load("Zebrafish_ambient_genes_Top50.rda")   # -> Z_ambient_genes_Top50
inj_pb <- Make_Pseudobulk_Matrix(Zebrafish_injury_smoothed[[ct_use]], Z_ambient_genes_Top50, "sample")
rownames(inj_pb) <- sub(".*~~", "", rownames(inj_pb))
zf_injury_col_order <- c("Control_2", "NMDA_36hr", "NMDA_54hr", "NMDA_72hr", "NMDA_96hr", "NMDA_7D")
inj_pb <- inj_pb[, zf_injury_col_order[zf_injury_col_order %in% colnames(inj_pb)], drop = FALSE]

# --- (b) aging/training pseudobulk（clock 训练数据；合并 replicates 与 28/30，按 age 排序）---
Zebrafish_seurat_knn_pseudobulk <- readRDS("Zebrafish_seurat_knn_pseudobulk_2025")
merge_zf_training_replicates <- function(pb) {
    cn <- colnames(pb)
    pb <- pb[, !cn %in% c("3_Rep2", "6_Rep3", "24_Rep1"), drop = FALSE]   # 去掉异常 replicate
    grp <- gsub("_Rep[0-9]+$", "", colnames(pb)); grp[grp %in% c("28", "30")] <- "28-30"
    ug <- unique(grp)
    merged <- sapply(ug, function(g) { idx <- which(grp == g); if (length(idx) == 1) pb[, idx] else rowMeans(pb[, idx, drop = FALSE]) })
    colnames(merged) <- ug; merged
}
trn_pb <- Make_Pseudobulk_Matrix(Zebrafish_seurat_knn_pseudobulk[[ct_use]], Z_ambient_genes_Top50, "sample")
rownames(trn_pb) <- sub(".*~~", "", rownames(trn_pb))
trn_pb <- merge_zf_training_replicates(trn_pb)
col_nums <- as.numeric(gsub("-.*", "", gsub("[^0-9-]", "", colnames(trn_pb))))
trn_pb <- trn_pb[, order(col_nums), drop = FALSE]

# --- (c) aging 相关 injury genes 按 clock 权重正负分两类 ---
genes_ct  <- zf_aging_injury_genes[zf_aging_injury_genes$celltype == ct_use, ]
genes_ct  <- genes_ct[!duplicated(genes_ct$gene), ]
genes_pos <- genes_ct$gene[genes_ct$clock_coefficient > 0]
genes_neg <- genes_ct$gene[genes_ct$clock_coefficient < 0]
message(sprintf("%s: aging injury genes = %d (pos %d / neg %d)",
                ct_use, nrow(genes_ct), length(genes_pos), length(genes_neg)))

######
###### 我们 筛选下这些基因，在 正的 weights 里面的 基因 ，我们计算下它的 z-score 和 每个 时间点 predicted aging 的 Spearman 相关系数 ##
###### 然后 设置一个 cutoff，选出 相关基因 ##
###### 在 负的 weights 里面，就选出 相关性 负的最低的 基因 ###
###### 选完之后只用这些基因 画 heatmaps #####
######
########################################################################
# 10) 基于 Spearman 相关性筛选基因 -> ComplexHeatmap（pos / neg 各一个 PNG）
########################################################################
library(ComplexHeatmap)
library(circlize)

coef_lookup <- setNames(genes_ct$clock_coefficient, genes_ct$gene)

# --- (d) 每个 injury 时间点的 predicted aging 向量（MG），与 inj_pb 列对齐 ---
mg_delta <- zf_sample_delta[zf_sample_delta$celltype == ct_use, ]
pred_map <- setNames(mg_delta$predicted, mg_delta$sample)
inj_key  <- ifelse(colnames(inj_pb) == "Control_2", "NMDA_0hr", colnames(inj_pb))  # Control_2 = NMDA_0hr
pred_vec <- pred_map[inj_key]
message("predicted aging（按 injury 时间点）: ",
        paste(sprintf("%s=%.1f", colnames(inj_pb), pred_vec), collapse = ", "))

# --- (e) 每个基因 injury z-score 与 predicted aging 的 Spearman 相关系数 ---
inj_z_all <- zscore_rows(inj_pb)
spearman_vs_pred <- function(genes) {
    g <- intersect(genes, rownames(inj_z_all))
    setNames(sapply(g, function(gene)
        suppressWarnings(cor(inj_z_all[gene, ], pred_vec, method = "spearman"))), g)
}
cor_pos <- spearman_vs_pred(genes_pos)   # 正权重基因
cor_neg <- spearman_vs_pred(genes_neg)   # 负权重基因

# --- (f) 设 cutoff：正权重取相关性最高(>cutoff)，负权重取相关性最低(<-cutoff) ---
cor_cutoff <- 0.25
sel_pos <- names(cor_pos)[!is.na(cor_pos) & cor_pos >  cor_cutoff]
sel_neg <- names(cor_neg)[!is.na(cor_neg) & cor_neg < -cor_cutoff]
sel_pos <- sel_pos[order(-cor_pos[sel_pos])]   # 相关性从高到低
sel_neg <- sel_neg[order( cor_neg[sel_neg])]   # 相关性从低到高（最负在上）
message(sprintf("正权重: %d/%d 通过 (cor > %.2f); 负权重: %d/%d 通过 (cor < %.2f)",
                length(sel_pos), length(cor_pos), cor_cutoff,
                length(sel_neg), length(cor_neg), -cor_cutoff))

# --- (g) 合并 heatmap：上 pos / 下 neg；左 Injury / 右 Aging（配色不同）；最右 clock weight ---
sel_genes <- c(sel_pos, sel_neg)
if (length(sel_genes) < 2) {
    message("筛选基因太少，跳过 heatmap")
} else {
    # 仅保留两边矩阵都有的基因；组内按 clock weight 从大到小
    sel_genes <- Reduce(intersect, list(sel_genes, rownames(inj_pb), rownames(trn_pb)))
    sel_pos   <- intersect(sel_pos, sel_genes)
    sel_neg   <- intersect(sel_neg, sel_genes)
    sel_pos   <- sel_pos[order(-coef_lookup[sel_pos])]   # 权重大 → 小
    sel_neg   <- sel_neg[order(-coef_lookup[sel_neg])]   # 权重大 → 小（负值里越接近 0 越靠上）
    sel_genes <- c(sel_pos, sel_neg)                     # 上 pos，下 neg

    inj_z <- zscore_rows(inj_pb[sel_genes, , drop = FALSE])
    trn_z <- zscore_rows(trn_pb[sel_genes, , drop = FALSE])

    # 左边列名与 boxplot 一致：Control_2 -> NMDA_0hr
    colnames(inj_z) <- ifelse(colnames(inj_z) == "Control_2", "NMDA_0hr", colnames(inj_z))

    row_split <- factor(
        ifelse(sel_genes %in% sel_pos, "Weight > 0", "Weight < 0"),
        levels = c("Weight > 0", "Weight < 0")
    )

    coef_val  <- coef_lookup[sel_genes]; coef_val[is.na(coef_val)] <- 0
    coef_disp <- pmax(pmin(coef_val, 0.5), -0.5)
    weight_fill <- ifelse(coef_disp > 0, "#E64B35", "#3C5488")
    weight_fill[coef_disp == 0] <- "grey70"

    gap_mm <- 6   # 上下间隔 = 左右间隔
    row_ha <- rowAnnotation(
        spacer = anno_empty(border = FALSE, width = unit(5, "mm")),  # 与 heatmap 拉开距离
        Clock_Weight = anno_barplot(
            coef_disp, baseline = 0,
            gp = gpar(fill = weight_fill, col = NA),
            axis_param = list(at = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5"),
                              labels_rot = 0, gp = gpar(fontsize = 14)),
            width = unit(2.0, "cm")
        ),
        gap = unit(2, "mm"),
        annotation_name_gp = gpar(fontsize = 15),
        show_annotation_name = c(spacer = FALSE, Clock_Weight = TRUE)
    )

    # Nature/NPG：左 Injury 蓝-白-红；右 Aging 青绿-白-橙（左右颜色不同）
    inj_col_fun <- colorRamp2(c(-2, 0, 2), c("#3C5488", "#F7F7F7", "#E64B35"))
    age_col_fun <- colorRamp2(c(-2, 0, 2), c("#00A087", "#F7F7F7", "#F39B7F"))

    show_rn <- length(sel_genes) <= 100
    ht_inj <- Heatmap(
        inj_z, name = "Injury Z", col = inj_col_fun,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_split = row_split, row_gap = unit(gap_mm, "mm"),
        row_title_gp = gpar(fontsize = 16, fontface = "bold"),
        show_row_names = show_rn, row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_rot = 45,
        column_title = "Injury", column_title_gp = gpar(fontsize = 18, fontface = "bold"),
        width = unit(4.5, "cm"), height = unit(min(12, 0.22 * length(sel_genes) + 2), "cm"),
        heatmap_legend_param = list(title = "Injury Z", title_gp = gpar(fontsize = 14),
                                    labels_gp = gpar(fontsize = 13))
    )
    ht_age <- Heatmap(
        trn_z, name = "Aging Z", col = age_col_fun,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_split = row_split, row_gap = unit(gap_mm, "mm"), row_title = NULL,
        right_annotation = row_ha,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_rot = 45,
        column_title = "Aging", column_title_gp = gpar(fontsize = 18, fontface = "bold"),
        width = unit(5.5, "cm"), height = unit(min(12, 0.22 * length(sel_genes) + 2), "cm"),
        heatmap_legend_param = list(title = "Aging Z", title_gp = gpar(fontsize = 14),
                                    labels_gp = gpar(fontsize = 13))
    )

    setwd(injury_dir)
    out_file <- "Zebrafish_MG_aging_injury_corr_heatmap_2026.png"
    png(out_file, width = 1400, height = 1100, res = 150)
    draw(ht_inj + ht_age,
         column_title = "Zebrafish MG: aging-related injury genes",
         column_title_gp = gpar(fontsize = 18, fontface = "bold"),
         ht_gap = unit(gap_mm, "mm"))
    dev.off()
    message(sprintf("已保存: %s (pos=%d, neg=%d)", file.path(injury_dir, out_file), length(sel_pos), length(sel_neg)))
}

# 保存筛选结果（统一写到 injury_dir）
setwd(injury_dir)
sel_df <- rbind(
    data.frame(gene = sel_pos, spearman = cor_pos[sel_pos],
               clock_coefficient = coef_lookup[sel_pos], group = "posWeight_posCorr"),
    data.frame(gene = sel_neg, spearman = cor_neg[sel_neg],
               clock_coefficient = coef_lookup[sel_neg], group = "negWeight_negCorr")
)
write.csv(sel_df, "Zebrafish_MG_aging_injury_selected_genes_2026.csv", row.names = FALSE)
message("已保存: ", file.path(injury_dir, "Zebrafish_MG_aging_injury_selected_genes_2026.csv"))

########################################################################
# 11) Zebrafish MG → Human ortholog → enrichGO BP（无 KEGG；不 simplify；取 top 10）
# 新终端: conda activate clusterProfiler ，从本段开始跑
########################################################################
if (!exists("injury_dir")) injury_dir <- "/projects/hmz-aging/Injury_Datasets"
if (!exists("model_dir"))  model_dir  <- "/projects/hmz-aging/Zebrafish_New_Datasets"

suppressMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})

setwd(injury_dir)

# --- 始终从 Zebrafish CSV 读入（勿复用同会话里其他物种的 sel_pos/sel_neg）---
sel_csv_candidates <- c(
    file.path(injury_dir, "Zebrafish_MG_aging_injury_selected_genes_2026.csv"),
    file.path(model_dir,  "Zebrafish_MG_aging_injury_selected_genes_2026.csv"),
    "Zebrafish_MG_aging_injury_selected_genes_2026.csv"
)
sel_csv <- sel_csv_candidates[file.exists(sel_csv_candidates)][1]
if (is.na(sel_csv) || !nzchar(sel_csv))
    stop("找不到 Zebrafish_MG_aging_injury_selected_genes_2026.csv，请先跑完筛选段")
zf_sel_df <- read.csv(sel_csv, stringsAsFactors = FALSE)
zf_sel_pos <- unique(zf_sel_df$gene[zf_sel_df$group == "posWeight_posCorr"])
zf_sel_neg <- unique(zf_sel_df$gene[zf_sel_df$group == "negWeight_negCorr"])
zf_sel_pos <- zf_sel_pos[!is.na(zf_sel_pos) & zf_sel_pos != ""]
zf_sel_neg <- zf_sel_neg[!is.na(zf_sel_neg) & zf_sel_neg != ""]
message(sprintf("从 CSV 读入 Zebrafish 基因: %s (pos=%d, neg=%d)",
                sel_csv, length(zf_sel_pos), length(zf_sel_neg)))
message("  示例 pos: ", paste(head(zf_sel_pos, 5), collapse = ", "))

zf_gene_groups <- list(
    posWeight = zf_sel_pos,
    negWeight = zf_sel_neg
)

# --- Zebrafish → Human ortholog ---
load_hz_ortholog <- function() {
    cands <- c(
        file.path(injury_dir, "HMZ_ortholog_combined_2025"),
        file.path(model_dir,  "HMZ_ortholog_combined_2025"),
        "/projects/hmz-aging/HMZ_ortholog_combined_2025",
        "HMZ_ortholog_combined_2025"
    )
    fp <- cands[file.exists(cands)][1]
    if (is.na(fp) || !nzchar(fp)) return(NULL)
    obj <- readRDS(fp)
    hz <- if (is.list(obj) && !is.null(obj$HZ)) obj$HZ else if (is.list(obj) && !is.null(obj$human_zebrafish)) obj$human_zebrafish else NULL
    if (is.null(hz)) return(NULL)
    cn <- tolower(colnames(hz)); colnames(hz) <- cn
    if (!all(c("human", "zebrafish") %in% colnames(hz))) {
        if (all(c("human_gene", "zebrafish_gene") %in% colnames(hz))) {
            hz <- data.frame(human = hz$human_gene, zebrafish = hz$zebrafish_gene, stringsAsFactors = FALSE)
        } else return(NULL)
    }
    hz <- hz[!is.na(hz$human) & !is.na(hz$zebrafish) & hz$human != "" & hz$zebrafish != "", ]
    hz <- hz[!duplicated(hz$zebrafish), ]
    message("已读入 ortholog: ", fp, " (HZ pairs=", nrow(hz), ")")
    hz
}

hz_map <- load_hz_ortholog()

zebrafish_to_human <- function(zf_genes) {
    zf_genes <- unique(as.character(zf_genes))
    zf_genes <- zf_genes[!is.na(zf_genes) & zf_genes != ""]
    zf_genes_clean <- sub(".*~~", "", zf_genes)
    out <- data.frame(zebrafish = zf_genes_clean, human = NA_character_, stringsAsFactors = FALSE)

    if (!is.null(hz_map)) {
        m <- match(out$zebrafish, hz_map$zebrafish)
        out$human <- hz_map$human[m]
        miss <- which(is.na(out$human))
        if (length(miss) > 0) {
            m2 <- match(tolower(out$zebrafish[miss]), tolower(hz_map$zebrafish))
            out$human[miss] <- hz_map$human[m2]
        }
    }

    still <- which(is.na(out$human))
    if (length(still) > 0 && requireNamespace("babelgene", quietly = TRUE)) {
        bg <- tryCatch(
            babelgene::orthologs(genes = out$zebrafish[still], species = "zebrafish", human = TRUE),
            error = function(e) NULL
        )
        if (!is.null(bg) && nrow(bg) > 0) {
            sc <- intersect(c("symbol", "zebrafish_symbol", "input"), colnames(bg))[1]
            hc <- intersect(c("human_symbol", "human", "symbol_human"), colnames(bg))[1]
            if (!is.na(sc) && !is.na(hc)) {
                bg <- bg[!duplicated(bg[[sc]]), ]
                m3 <- match(out$zebrafish[still], bg[[sc]])
                out$human[still] <- as.character(bg[[hc]][m3])
            }
        }
    }

    still2 <- which(is.na(out$human))
    if (length(still2) > 0) out$human[still2] <- toupper(out$zebrafish[still2])

    out <- out[!is.na(out$human) & out$human != "", ]
    out <- out[!duplicated(out$human), ]
    out
}

parse_ratio <- function(r) {
    parts <- as.numeric(strsplit(as.character(r), "/")[[1]])
    parts[1] / parts[2]
}

# enrichGO BP；按 p.adjust 取 top 10（不 simplify）
run_go_human <- function(zf_genes, group_name, p_cut = 0.05, q_cut = 0.2, top_n = 10) {
    message(sprintf("\n========== %s: zebrafish genes=%d ==========", group_name, length(zf_genes)))
    map_df <- zebrafish_to_human(zf_genes)
    human_genes <- unique(map_df$human)
    message(sprintf("  mapped to human: %d / %d", length(human_genes),
                    length(unique(sub(".*~~", "", zf_genes)))))
    if (length(human_genes) < 5) { message("  human 基因太少，跳过"); return(NULL) }

    entrez <- na.omit(AnnotationDbi::mapIds(
        org.Hs.eg.db, keys = human_genes, column = "ENTREZID",
        keytype = "SYMBOL", multiVals = "first"
    ))
    entrez <- unique(as.character(entrez))
    message(sprintf("  human EntrezID: %d", length(entrez)))
    if (length(entrez) < 5) { message("  Entrez 太少，跳过"); return(NULL) }

    ego <- enrichGO(
        gene          = entrez,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        keyType       = "ENTREZID",
        pAdjustMethod = "BH",
        pvalueCutoff  = p_cut,
        qvalueCutoff  = q_cut,
        readable      = TRUE
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        message("  无 GO BP 结果"); return(NULL)
    }
    d <- as.data.frame(ego)
    d <- d[d$pvalue < p_cut, , drop = FALSE]
    if (nrow(d) == 0) { message("  无显著 GO"); return(NULL) }

    d$GeneRatio_val  <- vapply(d$GeneRatio, parse_ratio, numeric(1))
    d$BgRatio_val    <- vapply(d$BgRatio, parse_ratio, numeric(1))
    d$FoldEnrichment <- d$GeneRatio_val / d$BgRatio_val
    d$neglog10_padj  <- -log10(pmax(d$p.adjust, 1e-300))
    d$Database <- "GO:BP"
    d$group <- group_name
    d <- d[order(d$p.adjust, d$pvalue), ]
    message(sprintf("  GO BP sig: %d (画图取 top %d)", nrow(d), top_n))

    list(
        ego = ego,
        GO = d,
        go_top = head(d, top_n),
        ortholog_map = map_df,
        human_genes = human_genes,
        entrez = entrez
    )
}

zf_enrich_hs <- list(
    posWeight = run_go_human(zf_gene_groups$posWeight, "posWeight", top_n = 10),
    negWeight = run_go_human(zf_gene_groups$negWeight, "negWeight", top_n = 10)
)
saveRDS(zf_enrich_hs, "Zebrafish_MG_aging_injury_GO_humanOrtholog_2026.rds")

all_tabs <- list(); map_tabs <- list()
for (nm in names(zf_enrich_hs)) {
    x <- zf_enrich_hs[[nm]]
    if (is.null(x)) next
    if (!is.null(x$GO) && nrow(x$GO) > 0) {
        keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
                            "qvalue", "geneID", "Count", "FoldEnrichment", "neglog10_padj",
                            "Database", "group"), colnames(x$GO))
        all_tabs[[nm]] <- x$GO[, keep, drop = FALSE]
    }
    if (!is.null(x$ortholog_map)) {
        md <- x$ortholog_map; md$group <- nm
        map_tabs[[nm]] <- md
    }
}
if (length(all_tabs) > 0) {
    write.csv(do.call(rbind, all_tabs),
              "Zebrafish_MG_aging_injury_GO_humanOrtholog_2026.csv", row.names = FALSE)
    message("已保存: Zebrafish_MG_aging_injury_GO_humanOrtholog_2026.csv")
}
if (length(map_tabs) > 0) {
    write.csv(do.call(rbind, map_tabs),
              "Zebrafish_MG_aging_injury_zf2human_map_2026.csv", row.names = FALSE)
    message("已保存: Zebrafish_MG_aging_injury_zf2human_map_2026.csv")
}

message("\n=== 第 11 段完成（Human ortholog GO BP, top 10）===")
message("下一步: 跑第 12 段画 Dotplot")


########################################################################
# 12) 发表级 Dotplot：前 10 个 GO BP（pos / neg 各一张；Nature/NPG 配色）
#     x = GeneRatio, size = Count, color = p.adjust
########################################################################
if (!exists("injury_dir")) injury_dir <- "/projects/hmz-aging/Injury_Datasets"
setwd(injury_dir)

suppressMessages(library(ggplot2))

enrich_rds <- "Zebrafish_MG_aging_injury_GO_humanOrtholog_2026.rds"
if (!exists("zf_enrich_hs") || is.null(zf_enrich_hs)) {
    if (!file.exists(enrich_rds)) stop("找不到 ", enrich_rds, "，请先跑第 11 段")
    zf_enrich_hs <- readRDS(enrich_rds)
    message("已读入: ", enrich_rds)
}

top_n_plot <- 10
wrap_term <- function(x, width = 40) paste(strwrap(x, width = width), collapse = "\n")
# Nature/NPG：显著(小 padj)=红，不显著(大 padj)=蓝
npg_lo <- "#E64B35"
npg_hi <- "#3C5488"
grp_title <- c(
    posWeight = "Positive clock weight",
    negWeight = "Negative clock weight"
)

for (nm in names(zf_enrich_hs)) {
    x <- zf_enrich_hs[[nm]]
    if (is.null(x) || is.null(x$GO) || nrow(x$GO) == 0) {
        message(nm, ": 无 GO 可画"); next
    }
    df <- head(x$GO[order(x$GO$p.adjust, x$GO$pvalue), ], top_n_plot)
    if (!"GeneRatio_val" %in% colnames(df)) {
        df$GeneRatio_val <- vapply(df$GeneRatio, function(r) {
            parts <- as.numeric(strsplit(as.character(r), "/")[[1]]); parts[1] / parts[2]
        }, numeric(1))
    }
    df$Description_wrapped <- vapply(df$Description, wrap_term, character(1))
    df <- df[order(df$GeneRatio_val), ]
    df$Description_wrapped <- factor(df$Description_wrapped, levels = unique(df$Description_wrapped))

    p <- ggplot(df, aes(x = GeneRatio_val, y = Description_wrapped,
                        size = Count, color = p.adjust)) +
        geom_point(alpha = 0.92) +
        scale_color_gradient(
            low = npg_lo, high = npg_hi, name = "p.adjust",
            guide = guide_colorbar(reverse = FALSE, barwidth = 0.45, barheight = 4)
        ) +
        scale_size_continuous(range = c(2.5, 7), name = "Count",
                              breaks = pretty(df$Count, n = 4)) +
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.10))) +
        labs(
            title = sprintf("Zebrafish MG — %s", grp_title[[nm]]),
            subtitle = "GO BP top 10 (human orthologs)",
            x = "GeneRatio",
            y = NULL
        ) +
        theme_classic(base_size = 12) +
        theme(
            plot.title    = element_text(size = 13, face = "bold", color = "black", hjust = 0),
            plot.subtitle = element_text(size = 10, color = "grey35", margin = margin(b = 6)),
            axis.text.y   = element_text(size = 13, color = "black", lineheight = 0.95),
            axis.text.x   = element_text(size = 9, color = "black"),   # GeneRatio 刻度缩小，避免重合
            axis.title.x  = element_text(size = 12, face = "bold", margin = margin(t = 6)),
            axis.line     = element_blank(),
            axis.ticks    = element_line(linewidth = 0.5, color = "black"),
            panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.9),
            legend.title  = element_text(size = 10, face = "bold"),
            legend.text   = element_text(size = 9),
            legend.position = "right",
            legend.key.size = unit(0.35, "cm"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor   = element_blank(),
            plot.margin   = margin(8, 10, 8, 8)
        )

    out_png <- sprintf("Zebrafish_MG_GO_humanOrtholog_%s_dotplot_top10_2026.png", nm)
    ggsave(out_png, p,
           width = 5.8, height = max(3.5, 0.35 * nrow(df) + 1.3),
           dpi = 600, bg = "white")
    message("已保存: ", out_png, " (600 dpi)")
}

message("\n=== Zebrafish GO Dotplot (human ortholog, top 10) 完成 ===")
