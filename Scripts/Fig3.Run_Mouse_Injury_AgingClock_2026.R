#!/usr/bin/env Rscript
########################################################################
# Mouse Injury — Aging Clock Prediction (2026)
# 用 2026 final aging clock 模型预测 mouse NMDA injury 数据的衰老评分
# celltypes: MG (仅 Resting MG), Rod, Cone
# 流程（与 final model 训练完全一致，无 calibration）:
#   ENSMUSG->symbol -> nFeature 过滤 -> 按 Condition KNN smoothing
#   -> 去 ambient -> depth norm(1e4)+log2 -> 按 genes_used 对齐
#   -> 用训练存的 means/sds 做 Z-score -> predict
# HPC 用法: Rscript Run_Mouse_Injury_AgingClock_2026.R
########################################################################

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(FNN)
    library(glmnet)     # 必须: 否则 predict 找不到 predict.glmnet
    library(ggplot2)
})

## 路径配置 ------------------------------------------------------------
mouse_dir  <- "/projects/hmz-aging/Mouse_New_Datasets"     # 原始 injury 数据 + 模型 + ambient
injury_dir <- "/projects/hmz-aging/Injury_Datasets"        # gene features + 输出
model_dir  <- mouse_dir                                     # 模型: <CT>/<CT>_final_aging_clock_model_2026.rds
clock_celltypes <- c("MG", "Rod", "Cone")
control_sample  <- "mmP60"
sample_col      <- "Condition"


########################################################################
# 1) 读入 mouse NMDA injury 数据，加 celltype 并映射为 MG/Rod/Cone 等
########################################################################
setwd(mouse_dir)
load("Mouse_NMDA_pbmc_seurat.RData")     # -> pbmc
Mouse_NMDA <- UpdateSeuratObject(pbmc)

cells_meta <- read.table("Mouse_NMDA_cell_features.tsv", sep = "\t")
m <- match(colnames(Mouse_NMDA), cells_meta$V1)
Mouse_NMDA$celltype <- cells_meta$V3[m]

# 细粒度 celltype -> clock celltype
map_list <- list(
    MG   = c("Resting MG"),
    Cone = c("Cones"),
    AC   = c("GABAergic AC", "Glycinergic AC"),
    HC   = c("HC"),
    RGC  = c("RGC"),
    BC   = c("OFF cone BC", "ON cone BC", "Rod BC"),
    Rod  = c("Rods")
)
Mouse_NMDA$celltype2 <- Mouse_NMDA$celltype
for (cl in names(map_list)) Mouse_NMDA$celltype2[Mouse_NMDA$celltype %in% map_list[[cl]]] <- cl

Mouse_NMDA_cl <- Mouse_NMDA[, Mouse_NMDA$celltype2 %in% names(map_list)]
Mouse_NMDA_list <- SplitObject(Mouse_NMDA_cl, split.by = "celltype2")


########################################################################
# 2) ENSMUSG ID -> gene symbol，合并重复 symbol 的 raw counts
########################################################################
Genes_meta <- read.table(file.path(injury_dir, "Mouse_gene_features.tsv"), sep = "\t")

convert_ensmusg_to_symbol <- function(seurat_list, genes_meta) {
    id_to_symbol <- setNames(genes_meta$V2, genes_meta$V1)
    new_list <- list()
    for (ct in names(seurat_list)) {
        seu <- seurat_list[[ct]]
        counts <- GetAssayData(seu, layer = "counts")
        old_genes <- rownames(counts)
        matched <- id_to_symbol[old_genes]
        new_names <- ifelse(!is.na(matched), matched, old_genes)
        dup_symbols <- unique(new_names[duplicated(new_names)])
        message(sprintf("  %s: 基因 %d, 匹配 symbol %d, 重复 %d",
                        ct, length(old_genes), sum(!is.na(matched)), length(dup_symbols)))
        if (length(dup_symbols) > 0) {
            uniq_idx <- which(!new_names %in% dup_symbols)
            uniq_counts <- counts[uniq_idx, , drop = FALSE]; rownames(uniq_counts) <- new_names[uniq_idx]
            dup_list <- lapply(dup_symbols, function(sym) Matrix::colSums(counts[which(new_names == sym), , drop = FALSE]))
            dup_mat <- as(do.call(rbind, dup_list), "dgCMatrix"); rownames(dup_mat) <- dup_symbols
            final_counts <- rbind(uniq_counts, dup_mat)
        } else {
            rownames(counts) <- new_names; final_counts <- counts
        }
        new_list[[ct]] <- CreateSeuratObject(counts = final_counts, meta.data = seu@meta.data)
    }
    new_list
}
Mouse_injury_symbol <- convert_ensmusg_to_symbol(Mouse_NMDA_list, Genes_meta)


########################################################################
# 3) nFeature_RNA 过滤 + 只保留 MG/Rod/Cone
########################################################################
Mouse_injury_cl <- list()
for (ct in names(Mouse_injury_symbol)) {
    tmp <- Mouse_injury_symbol[[ct]]
    tmp <- tmp[, which(tmp$nFeature_RNA > 1000)]
    message(sprintf("  %s: -> %d cells (nFeature_RNA > 1000)", ct, ncol(tmp)))
    Mouse_injury_cl[[ct]] <- tmp
}
keep_cts <- clock_celltypes[clock_celltypes %in% names(Mouse_injury_cl)]
Mouse_injury_sel <- Mouse_injury_cl[keep_cts]


########################################################################
# 4) KNN smoothing（按 Condition；与 aging clock 训练预处理一致）
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

smooth_injury_list <- function(seurat_list, sample_col = "Condition", species_label = "Mouse") {
    smoothed <- list()
    for (ct in names(seurat_list)) {
        message(sprintf("\n========== %s - %s ==========", species_label, ct))
        tmp_sp <- SplitObject(seurat_list[[ct]], split.by = sample_col)
        nc <- sapply(tmp_sp, ncol)
        if (any(nc < 2)) { message("  去掉 <2 细胞样本: ", paste(names(nc)[nc < 2], collapse = ", ")); tmp_sp <- tmp_sp[nc >= 2] }
        if (length(tmp_sp) == 0) { message("  ", ct, ": 所有样本被过滤，跳过"); next }
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

Mouse_injury_smoothed <- smooth_injury_list(Mouse_injury_sel, sample_col = sample_col, species_label = "Mouse")

# 统一 sample 列 = Condition；MG 只保留 Resting MG（与之前一致：smooth 后再取中心为 Resting MG 的 pseudo-cell）
for (ct in names(Mouse_injury_smoothed)) Mouse_injury_smoothed[[ct]]$sample <- as.character(Mouse_injury_smoothed[[ct]][[sample_col]][, 1])
if ("MG" %in% names(Mouse_injury_smoothed)) {
    Mouse_injury_smoothed$MG <- Mouse_injury_smoothed$MG[, which(Mouse_injury_smoothed$MG$celltype == "Resting MG")]
    message(sprintf("MG 仅保留 Resting MG: %d pseudo-cells", ncol(Mouse_injury_smoothed$MG)))
}

setwd(injury_dir)
saveRDS(Mouse_injury_smoothed, file = "Mouse_injury_smoothed.rds")
message("已保存: Mouse_injury_smoothed.rds")


########################################################################
# 5) 载入 ambient 基因 + 2026 final 模型，逐 celltype 预测（无 calibration）
########################################################################
setwd(mouse_dir)
load("Mouse_ambient_genes_Top50.rda")    # -> M_ambient_genes_Top50

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

Mouse_injury_predictions <- list()
for (ct in clock_celltypes) {
    if (!ct %in% names(Mouse_injury_smoothed)) { message(ct, ": 无 smoothed 数据，跳过"); next }
    model_file <- file.path(model_dir, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    if (!file.exists(model_file)) { message(ct, ": 2026 模型不存在，跳过 -> ", model_file); next }
    final_model <- readRDS(model_file)
    seu <- Mouse_injury_smoothed[[ct]]
    if (length(Layers(seu)) > 1) seu <- JoinLayers(seu)
    res <- predict_final_clock(final_model, seu, ambient_genes = M_ambient_genes_Top50)
    cell_pred <- data.frame(
        cell                 = names(res$pred),
        sample               = as.character(seu@meta.data[["sample"]]),
        predicted            = res$pred,
        predicted_calibrated = res$pred,
        celltype             = ct,
        stringsAsFactors = FALSE
    )
    sample_pred <- aggregate(cbind(predicted, predicted_calibrated) ~ sample, data = cell_pred, FUN = mean)
    sample_pred$celltype <- ct
    Mouse_injury_predictions[[ct]] <- list(cell_predictions = cell_pred, sample_predictions = sample_pred)
    message(sprintf("%s: %d cells, %d samples (缺失基因 %d/%d)",
                    ct, nrow(cell_pred), nrow(sample_pred), res$n_missing_genes, res$n_genes))
    print(sample_pred)
}

setwd(injury_dir)
saveRDS(Mouse_injury_predictions, file = "Mouse_injury_predictions_2026.rds")
message("已保存: Mouse_injury_predictions_2026.rds")



#

########################################################################
# 6) Delta aging：以 mmP60 的均值归一化（injury vs control）
########################################################################
compute_delta_aging <- function(predictions_list, control_sample, species_label = "Mouse") {
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

message("\nMouse sample names:")
for (ct in names(Mouse_injury_predictions))
    message(sprintf("  %s: %s", ct, paste(unique(Mouse_injury_predictions[[ct]]$cell_predictions$sample), collapse = ", ")))

Mouse_delta <- compute_delta_aging(Mouse_injury_predictions, control_sample = control_sample)
# 展示重命名: mmP60 -> NMDA_0h; mmNMDA03 -> NMDA_03 ...（delta 已相对 mmP60 计算）
Mouse_delta$sample[Mouse_delta$sample == "mmP60"] <- "NMDA_0h"
Mouse_delta$sample <- sub("^mmNMDA(.*)$", "NMDA_\\1h", Mouse_delta$sample)
saveRDS(Mouse_delta, file = "Mouse_injury_delta_aging_2026.rds")

mouse_sample_delta <- aggregate(cbind(predicted, delta_aging) ~ sample + group + celltype,
                                data = Mouse_delta, FUN = mean)
mouse_sample_delta <- mouse_sample_delta[order(mouse_sample_delta$celltype, mouse_sample_delta$group), ]
message("\nMouse - Sample Level Delta Aging Score"); print(mouse_sample_delta, row.names = FALSE)
write.csv(mouse_sample_delta, "Mouse_injury_sample_delta_aging_2026.csv", row.names = FALSE)


########################################################################
# 7) 画图：delta by celltype + 时间点分面（Nature/NPG 配色，大字）
########################################################################
npg_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
              "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# (a) 每个 celltype: Control vs Injury 的 delta boxplot
plot_delta_boxplot <- function(delta_df, species_label = "Mouse") {
    ct_order <- clock_celltypes[clock_celltypes %in% delta_df$celltype]
    delta_df$celltype <- factor(delta_df$celltype, levels = ct_order)
    ggplot(delta_df, aes(x = celltype, y = delta_aging, fill = group)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        scale_fill_manual(values = c("Control" = "#4DBBD5", "Injury" = "#E64B35")) +
        labs(title = sprintf("%s: Injury-Induced Aging Acceleration", species_label),
             x = "Cell Type", y = "Delta Aging Score", fill = "Group") +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(size = 13, color = "black"),
              plot.title = element_text(size = 14, face = "bold"))
}
ggsave("Mouse_Injury_Delta_Aging_Boxplot_2026.png",
       plot_delta_boxplot(Mouse_delta), width = 5, height = 4, dpi = 300)
message("已保存: Mouse_Injury_Delta_Aging_Boxplot_2026.png")

# (b) 时间点分面：每个 celltype 一行，x = injury 时间点
plot_timepoint_facet <- function(delta_df, sample_order, ct_order, species_label = "Mouse") {
    delta_df <- delta_df[delta_df$sample %in% sample_order & delta_df$celltype %in% ct_order, ]
    delta_df$sample   <- factor(delta_df$sample, levels = sample_order)
    delta_df$celltype <- factor(delta_df$celltype, levels = ct_order)
    fill_cols <- setNames(rep(npg_cols, length.out = length(sample_order)), sample_order)
    ggplot(delta_df, aes(x = sample, y = delta_aging, fill = sample)) +
        geom_boxplot(outlier.size = 0.3, alpha = 0.85) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        facet_wrap(~ celltype, ncol = 1, scales = "free_y") +
        scale_fill_manual(values = fill_cols) +
        labs(x = NULL, y = "Relative aging changes to control (weeks)") +
        theme_classic(base_size = 15) +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
              axis.text.x  = element_text(angle = 45, hjust = 1, size = 15, color = "black"),
              axis.text.y  = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 18, face = "bold"),
              strip.text   = element_text(size = 17, face = "bold"),
              strip.background = element_rect(color = "black", fill = "white"),
              legend.position = "none")
}
# 画图时把 6h–36h（06/12/24/36h）的细胞合并成一个组 NMDA_6-36h
Mouse_delta_plot <- Mouse_delta
Mouse_delta_plot$sample[Mouse_delta_plot$sample %in%
    c("NMDA_06h", "NMDA_12h", "NMDA_24h", "NMDA_36h")] <- "NMDA_6-36h"

# x 轴顺序（NMDA_0h 为对照；6-36h 合并；mmNMDA48FI 不画）
mouse_sample_order <- c("NMDA_0h", "NMDA_03h", "NMDA_6-36h", "NMDA_48h", "NMDA_72h")
mouse_ct_order     <- c("MG", "Rod", "Cone")
ggsave("Mouse_Injury_Timepoint_Boxplot_2026.png",
       plot_timepoint_facet(Mouse_delta_plot, mouse_sample_order, mouse_ct_order),
       width = 4.5, height = length(mouse_ct_order) * 3 + 1, dpi = 300)
message("已保存: Mouse_Injury_Timepoint_Boxplot_2026.png")

message("\n=== Mouse injury aging clock (2026) 全部完成 ===")



####### 
setwd(injury_dir)
saveRDS(Mouse_injury_predictions, file = "Mouse_injury_predictions_2026.rds")
message("已保存: Mouse_injury_predictions_2026.rds")


########################################################################
# 8) 比较 48h (mmNMDA48) vs 48FIh (mmNMDA48FI) 的 predicted age
#    boxplot + Nature/NPG 配色 + Wilcoxon p 值
########################################################################
library(ggplot2)
setwd(injury_dir)
Mouse_injury_predictions <- readRDS("Mouse_injury_predictions_2026.rds")

cell_all <- do.call(rbind, lapply(Mouse_injury_predictions, function(x) x$cell_predictions))
cmp_samples <- c("mmNMDA48", "mmNMDA48FI")
d <- cell_all[cell_all$sample %in% cmp_samples, , drop = FALSE]
d$grp <- factor(ifelse(d$sample == "mmNMDA48", "NMDA_48h", "NMDA_48FIh"),
                levels = c("NMDA_48h", "NMDA_48FIh"))
ct_order <- c("MG", "Rod", "Cone")
d$celltype <- factor(d$celltype, levels = ct_order[ct_order %in% unique(d$celltype)])

# 每个 celltype 的 Wilcoxon p（48h vs 48FIh）
fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.2g", p)
pdat <- do.call(rbind, lapply(levels(d$celltype), function(ct) {
    dd <- d[d$celltype == ct, ]
    pv <- tryCatch(wilcox.test(predicted ~ grp, data = droplevels(dd))$p.value, error = function(e) NA)
    rng <- range(dd$predicted, na.rm = TRUE)
    data.frame(celltype = ct, label = fmt_p(pv), y = rng[2] + diff(rng) * 0.06)
}))
pdat$celltype <- factor(pdat$celltype, levels = levels(d$celltype))

box_cols <- c("NMDA_48h" = "#4DBBD5", "NMDA_48FIh" = "#E64B35")

p_cmp <- ggplot(d, aes(x = grp, y = predicted, fill = grp)) +
    geom_boxplot(width = 0.6, outlier.size = 0.4, alpha = 0.85) +
    scale_fill_manual(values = box_cols, guide = "none") +
    facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
    geom_text(data = pdat, aes(x = 1.5, y = y, label = label),
              inherit.aes = FALSE, vjust = -0.2, size = 4.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = "Predicted age (weeks)") +
    theme_classic(base_size = 15) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.text.x  = element_text(angle = 30, hjust = 1, size = 14, color = "black"),
          axis.text.y  = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 16, face = "bold"),
          strip.text   = element_text(size = 16, face = "bold"),
          strip.background = element_rect(color = "black", fill = "white"))

ggsave("Mouse_Injury_48h_vs_48FIh_boxplot_2026.png", p_cmp,
       width = 2.6 * nlevels(d$celltype) + 1, height = 4.2, dpi = 300)
message("已保存: Mouse_Injury_48h_vs_48FIh_boxplot_2026.png")


#####
##### 和 zebraFish 一样，先 计算 injury DEGs，NMDA_6-36h 合在一起计算 DEGs #####
#####
########################################################################
# 9) aging 相关的 injury genes (Mouse):
#    (a) 在未 smooth 的原始 counts (Mouse_injury_sel) 上找 injury DEGs
#        - MG 仅 Resting MG；mmNMDA06/12/24/36 合并为 mmNMDA6-36h；排除 mmNMDA48FI
#    (b) 从 2026 final model 提取 clock 基因权重
#    (c) injury DEG 且是 clock gene = aging 相关 injury genes
########################################################################
library(Seurat)
library(Matrix)

# 原始未 smooth 对象: Mouse_injury_sel（第 3 段生成）；统一 sample=Condition
mouse_deg_input <- Mouse_injury_sel
for (ct in names(mouse_deg_input)) mouse_deg_input[[ct]]$sample <- as.character(mouse_deg_input[[ct]]$Condition)
# MG 仅保留 Resting MG（与 clock/heatmap 口径一致）
if ("MG" %in% names(mouse_deg_input)) {
    mg <- mouse_deg_input$MG
    mouse_deg_input$MG <- mg[, which(mg$celltype == "Resting MG")]
}
# 合并 6-36h（mmNMDA06/12/24/36 -> mmNMDA6-36h），排除 mmNMDA48FI
merge6_36 <- c("mmNMDA06", "mmNMDA12", "mmNMDA24", "mmNMDA36")
for (ct in names(mouse_deg_input)) {
    s <- mouse_deg_input[[ct]]$sample
    s[s %in% merge6_36] <- "mmNMDA6-36h"
    mouse_deg_input[[ct]]$sample <- s
    mouse_deg_input[[ct]] <- mouse_deg_input[[ct]][, s != "mmNMDA48FI"]
}

# (a) injury DEGs：每个时间点 vs mmP60（p_val_adj<0.01, |log2FC|>0.2）
find_injury_degs <- function(seurat_list, control_sample) {
    all_degs <- list()
    for (ct in names(seurat_list)) {
        message(sprintf("\n========== Mouse - %s ==========", ct))
        seu <- seurat_list[[ct]]
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
            sig$gene <- rownames(sig); sig$timepoint <- tp; sig$celltype <- ct
            sig$direction <- ifelse(sig$avg_log2FC > 0, "up", "down")
            ct_degs[[tp]] <- sig
            message(sprintf("  %s vs %s: %d DEGs", tp, control_sample, nrow(sig)))
        }
        if (length(ct_degs) > 0) { all_degs[[ct]] <- do.call(rbind, ct_degs); rownames(all_degs[[ct]]) <- NULL }
    }
    all_degs
}

Mouse_injury_degs <- find_injury_degs(mouse_deg_input, control_sample = control_sample)  # control_sample = "mmP60"
setwd(injury_dir)
saveRDS(Mouse_injury_degs, file = "Mouse_injury_degs_list_2026.rds")

# (b) 从 2026 final model 提取 clock 基因权重（MG/Rod/Cone）
mouse_gene_weights <- list()
for (ct in clock_celltypes) {
    mf <- file.path(model_dir, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    if (!file.exists(mf)) next
    fm <- readRDS(mf)
    lambda_idx <- which.min(abs(fm$model$lambda - fm$lambda))
    beta_vec <- as.numeric(fm$model$beta[, lambda_idx])
    wt_df <- data.frame(gene = rownames(fm$model$beta), clock_coefficient = beta_vec, stringsAsFactors = FALSE)
    wt_df <- wt_df[wt_df$clock_coefficient != 0, ]
    wt_df$clock_abs_coef <- abs(wt_df$clock_coefficient)
    wt_df <- wt_df[order(wt_df$clock_abs_coef, decreasing = TRUE), ]
    wt_df$clock_rank <- seq_len(nrow(wt_df))
    mouse_gene_weights[[ct]] <- wt_df
    message(sprintf("  Mouse %s: %d clock genes", ct, nrow(wt_df)))
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
            mc <- intersect(c("gene", "clock_coefficient", "clock_abs_coef", "clock_rank"), colnames(wt))
            m <- merge(deg_df, wt[, mc, drop = FALSE], by = "gene", all.x = TRUE)
            m$is_clock_gene <- !is.na(m$clock_coefficient)
            deg_df <- m
        }
        out[[ct]] <- deg_df
        message(sprintf("  %s: %d DEGs, clock overlap=%d", ct, nrow(deg_df), sum(deg_df$is_clock_gene, na.rm = TRUE)))
    }
    out
}

Mouse_injury_degs_table <- merge_degs_with_weights(Mouse_injury_degs, mouse_gene_weights)
mouse_deg_table <- do.call(rbind, Mouse_injury_degs_table); rownames(mouse_deg_table) <- NULL
write.csv(mouse_deg_table, "Mouse_injury_DEGs_with_clock_weights_2026.csv", row.names = FALSE)
saveRDS(Mouse_injury_degs_table, file = "Mouse_injury_degs_table_2026.rds")

# aging 相关 injury genes = injury DEG 且是 clock gene
mouse_aging_injury_genes <- mouse_deg_table[mouse_deg_table$is_clock_gene %in% TRUE, , drop = FALSE]
write.csv(mouse_aging_injury_genes, "Mouse_aging_related_injury_genes_2026.csv", row.names = FALSE)
message(sprintf("\naging 相关 injury genes: %d 行, %d unique genes",
                nrow(mouse_aging_injury_genes), length(unique(mouse_aging_injury_genes$gene))))


##############
############## 接下来，把 mouse_aging_injury_genes 按 weights 分为正/负两类 #####
############## 用 injury z-score 与各时间点 predicted aging 的 Spearman 相关筛选 #####
############## 只做 MG；heatmap 左 Injury | 右 Aging；上 pos / 下 neg #####
########################################################################
# 10) MG: 构建 injury/aging pseudobulk + Spearman 筛选 + ComplexHeatmap
########################################################################
library(Matrix)
library(ComplexHeatmap)
library(circlize)

ct_use <- "MG"

# pseudobulk: raw counts 之和 -> depth norm 1e4 -> log2(1+x)
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

# --- (a) injury pseudobulk（smoothed MG；合并 6-36h；排除 48FI；列名与 boxplot 一致）---
# 用 KNN smoothed 后的 injury 数据做 pb（与预测一致）；若不在内存则读 RDS
setwd(injury_dir)
if (!exists("Mouse_injury_smoothed")) {
    Mouse_injury_smoothed <- readRDS("Mouse_injury_smoothed.rds")
    message("已读入: Mouse_injury_smoothed.rds")
}
setwd(mouse_dir)
load("Mouse_ambient_genes_Top50.rda")   # -> M_ambient_genes_Top50
inj_pb <- Make_Pseudobulk_Matrix(Mouse_injury_smoothed[[ct_use]], M_ambient_genes_Top50, "sample")
# 去掉 48FI
inj_pb <- inj_pb[, colnames(inj_pb) != "mmNMDA48FI", drop = FALSE]
# 合并 6-36h
merge6_36 <- c("mmNMDA06", "mmNMDA12", "mmNMDA24", "mmNMDA36")
present6  <- intersect(merge6_36, colnames(inj_pb))
if (length(present6) > 0) {
    merged_col <- rowMeans(inj_pb[, present6, drop = FALSE])
    keep_cols  <- setdiff(colnames(inj_pb), merge6_36)
    inj_pb     <- cbind(inj_pb[, keep_cols, drop = FALSE], `mmNMDA6-36h` = merged_col)
}
# 重命名为 boxplot 标注
rename_inj <- c(mmP60 = "NMDA_0h", mmNMDA03 = "NMDA_03h", `mmNMDA6-36h` = "NMDA_6-36h",
                mmNMDA48 = "NMDA_48h", mmNMDA72 = "NMDA_72h")
cn <- colnames(inj_pb)
cn <- ifelse(cn %in% names(rename_inj), rename_inj[cn], cn)
colnames(inj_pb) <- cn
mouse_injury_col_order <- c("NMDA_0h", "NMDA_03h", "NMDA_6-36h", "NMDA_48h", "NMDA_72h")
inj_pb <- inj_pb[, mouse_injury_col_order[mouse_injury_col_order %in% colnames(inj_pb)], drop = FALSE]

# --- (b) aging/training pseudobulk（clock 训练数据；合并 replicates / Wk106-108）---
Mouse_seurat_knn_pseudobulk <- readRDS("Mouse_seurat_knn_pseudobulk_2025")
merge_mouse_training_replicates <- function(pb) {
    grp <- gsub("_Rep[0-9]+$", "", colnames(pb))
    grp[grp %in% c("Wk106", "Wk108")] <- "Wk106-108"
    ug <- unique(grp)
    merged <- sapply(ug, function(g) {
        idx <- which(grp == g)
        if (length(idx) == 1) pb[, idx] else rowMeans(pb[, idx, drop = FALSE])
    })
    colnames(merged) <- ug
    merged
}
trn_pb <- Make_Pseudobulk_Matrix(Mouse_seurat_knn_pseudobulk[[ct_use]], M_ambient_genes_Top50, "sample")
trn_pb <- merge_mouse_training_replicates(trn_pb)
mouse_training_col_order <- c("Wk5", "Wk12", "Wk17", "Wk32", "Wk49", "Wk68", "Wk91", "Wk106-108", "Wk120")
trn_pb <- trn_pb[, mouse_training_col_order[mouse_training_col_order %in% colnames(trn_pb)], drop = FALSE]

# --- (c) aging 相关 injury genes 按 clock 权重正负分两类（仅 MG）---
if (!exists("mouse_aging_injury_genes")) {
    setwd(injury_dir)
    mouse_aging_injury_genes <- read.csv("Mouse_aging_related_injury_genes_2026.csv",
                                         stringsAsFactors = FALSE)
    message("已读入: Mouse_aging_related_injury_genes_2026.csv")
}
genes_ct  <- mouse_aging_injury_genes[mouse_aging_injury_genes$celltype == ct_use, ]
genes_ct  <- genes_ct[!duplicated(genes_ct$gene), ]
genes_pos <- genes_ct$gene[genes_ct$clock_coefficient > 0]
genes_neg <- genes_ct$gene[genes_ct$clock_coefficient < 0]
message(sprintf("%s: aging injury genes = %d (pos %d / neg %d)",
                ct_use, nrow(genes_ct), length(genes_pos), length(genes_neg)))

coef_lookup <- setNames(genes_ct$clock_coefficient, genes_ct$gene)

# --- (d) 每个 injury 时间点的 predicted aging（MG），与 inj_pb 列对齐（含 6-36h 合并）---
if (!exists("Mouse_delta")) {
    setwd(injury_dir)
    Mouse_delta <- readRDS("Mouse_injury_delta_aging_2026.rds")
    message("已读入: Mouse_injury_delta_aging_2026.rds")
}
md <- Mouse_delta[Mouse_delta$celltype == ct_use, ]
md$sample[md$sample %in% c("NMDA_06h", "NMDA_12h", "NMDA_24h", "NMDA_36h")] <- "NMDA_6-36h"
md <- md[!md$sample %in% "NMDA_48FIh", ]
pred_map <- setNames(aggregate(predicted ~ sample, md, mean)$predicted,
                     aggregate(predicted ~ sample, md, mean)$sample)
pred_vec <- pred_map[colnames(inj_pb)]
message("predicted aging（按 injury 时间点）: ",
        paste(sprintf("%s=%.1f", colnames(inj_pb), pred_vec), collapse = ", "))

# --- (e) Spearman：基因 injury z-score vs predicted aging ---
inj_z_all <- zscore_rows(inj_pb)
spearman_vs_pred <- function(genes) {
    g <- intersect(genes, rownames(inj_z_all))
    setNames(sapply(g, function(gene)
        suppressWarnings(cor(inj_z_all[gene, ], pred_vec, method = "spearman"))), g)
}
cor_pos <- spearman_vs_pred(genes_pos)
cor_neg <- spearman_vs_pred(genes_neg)

# --- (f) cutoff 筛选 ---
cor_cutoff <- 0.25
sel_pos <- names(cor_pos)[!is.na(cor_pos) & cor_pos >  cor_cutoff]
sel_neg <- names(cor_neg)[!is.na(cor_neg) & cor_neg < -cor_cutoff]
sel_pos <- sel_pos[order(-cor_pos[sel_pos])]
sel_neg <- sel_neg[order( cor_neg[sel_neg])]
message(sprintf("正权重: %d/%d 通过 (cor > %.2f); 负权重: %d/%d 通过 (cor < %.2f)",
                length(sel_pos), length(cor_pos), cor_cutoff,
                length(sel_neg), length(cor_neg), -cor_cutoff))

# --- (g) 合并 heatmap：上 pos / 下 neg；左 Injury / 右 Aging；最右 clock weight ---
setwd(injury_dir)
sel_genes <- c(sel_pos, sel_neg)
if (length(sel_genes) < 2) {
    message("筛选基因太少，跳过 heatmap")
} else {
    sel_genes <- Reduce(intersect, list(sel_genes, rownames(inj_pb), rownames(trn_pb)))
    sel_pos   <- intersect(sel_pos, sel_genes)
    sel_neg   <- intersect(sel_neg, sel_genes)
    sel_pos   <- sel_pos[order(-coef_lookup[sel_pos])]
    sel_neg   <- sel_neg[order(-coef_lookup[sel_neg])]
    sel_genes <- c(sel_pos, sel_neg)

    inj_z <- zscore_rows(inj_pb[sel_genes, , drop = FALSE])
    trn_z <- zscore_rows(trn_pb[sel_genes, , drop = FALSE])
    # 右边 Aging 列名：只显示周龄数字；Wk106-108 -> 108
    cn_age <- colnames(trn_z)
    cn_age <- sub("^Wk", "", cn_age)
    cn_age[cn_age == "106-108"] <- "108"
    colnames(trn_z) <- cn_age

    row_split <- factor(
        ifelse(sel_genes %in% sel_pos, "Weight > 0", "Weight < 0"),
        levels = c("Weight > 0", "Weight < 0")
    )

    coef_val  <- coef_lookup[sel_genes]; coef_val[is.na(coef_val)] <- 0
    coef_disp <- pmax(pmin(coef_val, 0.5), -0.5)
    weight_fill <- ifelse(coef_disp > 0, "#E64B35", "#3C5488")
    weight_fill[coef_disp == 0] <- "grey70"

    gap_mm <- 6
    # 与斑马鱼一致：anno_barplot
    row_ha <- rowAnnotation(
        spacer = anno_empty(border = FALSE, width = unit(5, "mm")),
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

    inj_col_fun <- colorRamp2(c(-2, 0, 2), c("#3C5488", "#F7F7F7", "#E64B35"))
    age_col_fun <- colorRamp2(c(-2, 0, 2), c("#00A087", "#F7F7F7", "#F39B7F"))

    # 尺寸与斑马鱼一致
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
        width = unit(5.2, "cm"), height = unit(min(14, 0.26 * length(sel_genes) + 2.5), "cm"),
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
        width = unit(6.2, "cm"), height = unit(min(14, 0.26 * length(sel_genes) + 2.5), "cm"),
        heatmap_legend_param = list(title = "Aging Z", title_gp = gpar(fontsize = 14),
                                    labels_gp = gpar(fontsize = 13))
    )

    out_file <- "Mouse_MG_aging_injury_corr_heatmap_2026.png"
    png(out_file, width = 1600, height = 1300, res = 150)
    draw(ht_inj + ht_age,
         column_title = "Mouse MG: aging-related injury genes",
         column_title_gp = gpar(fontsize = 18, fontface = "bold"),
         ht_gap = unit(gap_mm, "mm"))
    dev.off()
    message(sprintf("已保存: %s (pos=%d, neg=%d)", out_file, length(sel_pos), length(sel_neg)))
}

setwd(injury_dir)
sel_df <- rbind(
    data.frame(gene = sel_pos, spearman = cor_pos[sel_pos],
               clock_coefficient = coef_lookup[sel_pos], group = "posWeight_posCorr"),
    data.frame(gene = sel_neg, spearman = cor_neg[sel_neg],
               clock_coefficient = coef_lookup[sel_neg], group = "negWeight_negCorr")
)
write.csv(sel_df, "Mouse_MG_aging_injury_selected_genes_2026.csv", row.names = FALSE)
message("已保存: ", file.path(injury_dir, "Mouse_MG_aging_injury_selected_genes_2026.csv"))

########################################################################
# 11) Mouse MG → Human ortholog → enrichGO BP（无 KEGG；取 top 10）
# 流程对齐 Zebrafish；新终端: conda activate clusterProfiler
########################################################################
if (!exists("injury_dir")) injury_dir <- "/projects/hmz-aging/Injury_Datasets"
if (!exists("mouse_dir"))  mouse_dir  <- "/projects/hmz-aging/Mouse_New_Datasets"

suppressMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})

setwd(injury_dir)

# --- 始终从 Mouse CSV 读入（勿复用同会话里 Zebrafish 的 sel_pos/sel_neg）---
sel_csv_candidates <- c(
    file.path(injury_dir, "Mouse_MG_aging_injury_selected_genes_2026.csv"),
    file.path(mouse_dir,  "Mouse_MG_aging_injury_selected_genes_2026.csv"),
    "Mouse_MG_aging_injury_selected_genes_2026.csv"
)
sel_csv <- sel_csv_candidates[file.exists(sel_csv_candidates)][1]
if (is.na(sel_csv) || !nzchar(sel_csv))
    stop("找不到 Mouse_MG_aging_injury_selected_genes_2026.csv，请先跑完筛选段")
mouse_sel_df <- read.csv(sel_csv, stringsAsFactors = FALSE)
mouse_sel_pos <- unique(mouse_sel_df$gene[mouse_sel_df$group == "posWeight_posCorr"])
mouse_sel_neg <- unique(mouse_sel_df$gene[mouse_sel_df$group == "negWeight_negCorr"])
mouse_sel_pos <- mouse_sel_pos[!is.na(mouse_sel_pos) & mouse_sel_pos != ""]
mouse_sel_neg <- mouse_sel_neg[!is.na(mouse_sel_neg) & mouse_sel_neg != ""]
message(sprintf("从 CSV 读入 Mouse 基因: %s (pos=%d, neg=%d)",
                sel_csv, length(mouse_sel_pos), length(mouse_sel_neg)))
message("  示例 pos: ", paste(head(mouse_sel_pos, 5), collapse = ", "))

mouse_gene_groups <- list(
    posWeight = mouse_sel_pos,
    negWeight = mouse_sel_neg
)

# --- Mouse → Human ortholog ---
load_hm_ortholog <- function() {
    cands <- c(
        file.path(injury_dir, "HMZ_ortholog_combined_2025"),
        file.path(mouse_dir,  "HMZ_ortholog_combined_2025"),
        "/projects/hmz-aging/HMZ_ortholog_combined_2025",
        "HMZ_ortholog_combined_2025"
    )
    fp <- cands[file.exists(cands)][1]
    if (is.na(fp) || !nzchar(fp)) return(NULL)
    obj <- readRDS(fp)
    hm <- if (is.list(obj) && !is.null(obj$HM)) obj$HM else if (is.list(obj) && !is.null(obj$human_mouse)) obj$human_mouse else NULL
    if (is.null(hm)) return(NULL)
    cn <- tolower(colnames(hm)); colnames(hm) <- cn
    if (!all(c("human", "mouse") %in% colnames(hm))) {
        if (all(c("human_gene", "mouse_gene") %in% colnames(hm))) {
            hm <- data.frame(human = hm$human_gene, mouse = hm$mouse_gene, stringsAsFactors = FALSE)
        } else return(NULL)
    }
    hm <- hm[!is.na(hm$human) & !is.na(hm$mouse) & hm$human != "" & hm$mouse != "", ]
    hm <- hm[!duplicated(hm$mouse), ]
    message("已读入 ortholog: ", fp, " (HM pairs=", nrow(hm), ")")
    hm
}

hm_map <- load_hm_ortholog()

mouse_to_human <- function(mouse_genes) {
    mouse_genes <- unique(as.character(mouse_genes))
    mouse_genes <- mouse_genes[!is.na(mouse_genes) & mouse_genes != ""]
    out <- data.frame(mouse = mouse_genes, human = NA_character_, stringsAsFactors = FALSE)

    if (!is.null(hm_map)) {
        m <- match(out$mouse, hm_map$mouse)
        out$human <- hm_map$human[m]
        miss <- which(is.na(out$human))
        if (length(miss) > 0) {
            m2 <- match(tolower(out$mouse[miss]), tolower(hm_map$mouse))
            out$human[miss] <- hm_map$human[m2]
        }
    }

    still <- which(is.na(out$human))
    if (length(still) > 0 && requireNamespace("babelgene", quietly = TRUE)) {
        bg <- tryCatch(
            babelgene::orthologs(genes = out$mouse[still], species = "mouse", human = TRUE),
            error = function(e) NULL
        )
        if (!is.null(bg) && nrow(bg) > 0) {
            sc <- intersect(c("symbol", "mouse_symbol", "input"), colnames(bg))[1]
            hc <- intersect(c("human_symbol", "human", "symbol_human"), colnames(bg))[1]
            if (!is.na(sc) && !is.na(hc)) {
                bg <- bg[!duplicated(bg[[sc]]), ]
                m3 <- match(out$mouse[still], bg[[sc]])
                out$human[still] <- as.character(bg[[hc]][m3])
            }
        }
    }

    still2 <- which(is.na(out$human))
    if (length(still2) > 0) out$human[still2] <- toupper(out$mouse[still2])

    out <- out[!is.na(out$human) & out$human != "", ]
    out <- out[!duplicated(out$human), ]
    out
}

parse_ratio <- function(r) {
    parts <- as.numeric(strsplit(as.character(r), "/")[[1]])
    parts[1] / parts[2]
}

# enrichGO BP；按 p.adjust 取 top 10
run_go_human <- function(mouse_genes, group_name, p_cut = 0.05, q_cut = 0.2, top_n = 10) {
    message(sprintf("\n========== %s: mouse genes=%d ==========", group_name, length(mouse_genes)))
    map_df <- mouse_to_human(mouse_genes)
    human_genes <- unique(map_df$human)
    message(sprintf("  mapped to human: %d / %d", length(human_genes), length(unique(mouse_genes))))
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

mouse_enrich_hs <- list(
    posWeight = run_go_human(mouse_gene_groups$posWeight, "posWeight", top_n = 10),
    negWeight = run_go_human(mouse_gene_groups$negWeight, "negWeight", top_n = 10)
)
saveRDS(mouse_enrich_hs, "Mouse_MG_aging_injury_GO_humanOrtholog_2026.rds")

all_tabs <- list(); map_tabs <- list()
for (nm in names(mouse_enrich_hs)) {
    x <- mouse_enrich_hs[[nm]]
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
              "Mouse_MG_aging_injury_GO_humanOrtholog_2026.csv", row.names = FALSE)
    message("已保存: Mouse_MG_aging_injury_GO_humanOrtholog_2026.csv")
}
if (length(map_tabs) > 0) {
    write.csv(do.call(rbind, map_tabs),
              "Mouse_MG_aging_injury_mouse2human_map_2026.csv", row.names = FALSE)
    message("已保存: Mouse_MG_aging_injury_mouse2human_map_2026.csv")
}

message("\n=== 第 11 段完成（Mouse→Human GO BP, top 10）===")
message("下一步: 跑第 12 段画 Dotplot")


########################################################################
# 12) 发表级 Dotplot：前 10 个 GO BP（pos / neg 各一张；Nature/NPG 配色）
#     x = GeneRatio, size = Count, color = p.adjust（与 Zebrafish 一致）
########################################################################
if (!exists("injury_dir")) injury_dir <- "/projects/hmz-aging/Injury_Datasets"
setwd(injury_dir)

suppressMessages(library(ggplot2))

enrich_rds <- "Mouse_MG_aging_injury_GO_humanOrtholog_2026.rds"
if (!exists("mouse_enrich_hs") || is.null(mouse_enrich_hs)) {
    if (!file.exists(enrich_rds)) stop("找不到 ", enrich_rds, "，请先跑第 11 段")
    mouse_enrich_hs <- readRDS(enrich_rds)
    message("已读入: ", enrich_rds)
}

top_n_plot <- 10
wrap_term <- function(x, width = 40) paste(strwrap(x, width = width), collapse = "\n")
npg_lo <- "#E64B35"
npg_hi <- "#3C5488"
grp_title <- c(
    posWeight = "Positive clock weight",
    negWeight = "Negative clock weight"
)

for (nm in names(mouse_enrich_hs)) {
    x <- mouse_enrich_hs[[nm]]
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
            title = sprintf("Mouse MG — %s", grp_title[[nm]]),
            subtitle = "GO BP top 10 (human orthologs)",
            x = "GeneRatio",
            y = NULL
        ) +
        theme_classic(base_size = 12) +
        theme(
            plot.title    = element_text(size = 13, face = "bold", color = "black", hjust = 0),
            plot.subtitle = element_text(size = 10, color = "grey35", margin = margin(b = 6)),
            axis.text.y   = element_text(size = 13, color = "black", lineheight = 0.95),
            axis.text.x   = element_text(size = 9, color = "black"),
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

    out_png <- sprintf("Mouse_MG_GO_humanOrtholog_%s_dotplot_top10_2026.png", nm)
    ggsave(out_png, p,
           width = 5.8, height = max(3.5, 0.35 * nrow(df) + 1.3),
           dpi = 600, bg = "white")
    message("已保存: ", out_png, " (600 dpi)")
}

message("\n=== Mouse GO Dotplot (human ortholog, top 10) 完成 ===")
