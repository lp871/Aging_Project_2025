#!/usr/bin/env Rscript
########################################################################
# Mouse Injury — Aging Clock Prediction (2026)
# Predict aging scores for mouse NMDA injury data using the 2026 final aging clock models
# celltypes: MG (Resting MG only), Rod, Cone
# Pipeline (identical to final model training, no calibration):
#   ENSMUSG->symbol -> nFeature filtering -> per-Condition KNN smoothing
#   -> remove ambient -> depth norm(1e4)+log2 -> align to genes_used
#   -> Z-score with stored training means/sds -> predict
# HPC usage: Rscript Run_Mouse_Injury_AgingClock_2026.R
########################################################################

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(FNN)
    library(glmnet)     # required: otherwise predict cannot find predict.glmnet
    library(ggplot2)
})

## Path configuration ------------------------------------------------------------
mouse_dir  <- "/projects/hmz-aging/Mouse_New_Datasets"     # raw injury data + models + ambient
injury_dir <- "/projects/hmz-aging/Injury_Datasets"        # gene features + outputs
model_dir  <- mouse_dir                                     # models: <CT>/<CT>_final_aging_clock_model_2026.rds
clock_celltypes <- c("MG", "Rod", "Cone")
control_sample  <- "mmP60"
sample_col      <- "Condition"


########################################################################
# 1) Load mouse NMDA injury data, add celltypes and map them to MG/Rod/Cone etc.
########################################################################
setwd(mouse_dir)
load("Mouse_NMDA_pbmc_seurat.RData")     # -> pbmc
Mouse_NMDA <- UpdateSeuratObject(pbmc)

cells_meta <- read.table("Mouse_NMDA_cell_features.tsv", sep = "\t")
m <- match(colnames(Mouse_NMDA), cells_meta$V1)
Mouse_NMDA$celltype <- cells_meta$V3[m]

# Fine-grained celltype -> clock celltype
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
# 2) ENSMUSG ID -> gene symbol; sum raw counts of duplicated symbols
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
        message(sprintf("  %s: genes %d, matched symbols %d, duplicated %d",
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
# 3) nFeature_RNA filtering + keep only MG/Rod/Cone
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
# 4) KNN smoothing (per Condition; same preprocessing as aging clock training)
########################################################################
NormalizeAndPCA <- function(seu, scale.factor = 1e4, nfeatures = 2000, npcs = 30) {
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
    max_npcs <- min(nrow(seu), ncol(seu)) - 1
    if (npcs > max_npcs) npcs <- max(1, max_npcs)
    if (npcs < 1) stop(sprintf("Too few cells (n_cells=%d) to run PCA", ncol(seu)))
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
        if (any(nc < 2)) { message("  Removing samples with <2 cells: ", paste(names(nc)[nc < 2], collapse = ", ")); tmp_sp <- tmp_sp[nc >= 2] }
        if (length(tmp_sp) == 0) { message("  ", ct, ": all samples filtered out, skipping"); next }
        if (ct %in% c("Rod", "Cone", "RGC", "AC", "HC", "BC")) { k <- 15; npcs <- 25 }
        else if (ct == "MG") { k <- 5; npcs <- 5 }
        else { k <- 2; npcs <- 5 }
        min_ncells <- min(sapply(tmp_sp, ncol))
        message(sprintf("  min cells/sample=%d, K=%d, npcs=%d", min_ncells, k, npcs))
        if (min_ncells < 5) {
            message("  Smallest sample has <5 cells; skipping smoothing and merging directly")
            smoothed[[ct]] <- if (length(tmp_sp) == 1) tmp_sp[[1]] else merge(tmp_sp[[1]], tmp_sp[-1], add.cell.ids = names(tmp_sp))
            next
        }
        if (npcs >= min_ncells) npcs <- max(1, min_ncells - 1)
        if (k    >= min_ncells) k    <- max(1, min_ncells - 1)
        tmp_pca    <- sapply(tmp_sp, NormalizeAndPCA, npcs = npcs)
        tmp_smooth <- sapply(tmp_pca, PseudoBulkByKNN_fast, k = k, dims = 1:npcs)
        smoothed[[ct]] <- if (length(tmp_smooth) == 1) tmp_smooth[[1]] else merge(tmp_smooth[[1]], tmp_smooth[-1], add.cell.ids = names(tmp_smooth))
        message(sprintf("  Done: %d pseudo-cells", ncol(smoothed[[ct]])))
    }
    smoothed
}

Mouse_injury_smoothed <- smooth_injury_list(Mouse_injury_sel, sample_col = sample_col, species_label = "Mouse")

# Unify sample column = Condition; keep only Resting MG for MG (as before: after smoothing, keep pseudo-cells centered on Resting MG)
for (ct in names(Mouse_injury_smoothed)) Mouse_injury_smoothed[[ct]]$sample <- as.character(Mouse_injury_smoothed[[ct]][[sample_col]][, 1])
if ("MG" %in% names(Mouse_injury_smoothed)) {
    Mouse_injury_smoothed$MG <- Mouse_injury_smoothed$MG[, which(Mouse_injury_smoothed$MG$celltype == "Resting MG")]
    message(sprintf("MG restricted to Resting MG: %d pseudo-cells", ncol(Mouse_injury_smoothed$MG)))
}

setwd(injury_dir)
saveRDS(Mouse_injury_smoothed, file = "Mouse_injury_smoothed.rds")
message("Saved: Mouse_injury_smoothed.rds")


########################################################################
# 5) Load ambient genes + 2026 final models and predict per celltype (no calibration)
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
    if (!ct %in% names(Mouse_injury_smoothed)) { message(ct, ": no smoothed data, skipping"); next }
    model_file <- file.path(model_dir, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    if (!file.exists(model_file)) { message(ct, ": 2026 model not found, skipping -> ", model_file); next }
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
    message(sprintf("%s: %d cells, %d samples (missing genes %d/%d)",
                    ct, nrow(cell_pred), nrow(sample_pred), res$n_missing_genes, res$n_genes))
    print(sample_pred)
}

setwd(injury_dir)
saveRDS(Mouse_injury_predictions, file = "Mouse_injury_predictions_2026.rds")
message("Saved: Mouse_injury_predictions_2026.rds")



#

########################################################################
# 6) Delta aging: normalize to the mmP60 mean (injury vs control)
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
# Rename for display: mmP60 -> NMDA_0h; mmNMDA03 -> NMDA_03 ... (delta already computed relative to mmP60)
Mouse_delta$sample[Mouse_delta$sample == "mmP60"] <- "NMDA_0h"
Mouse_delta$sample <- sub("^mmNMDA(.*)$", "NMDA_\\1h", Mouse_delta$sample)
saveRDS(Mouse_delta, file = "Mouse_injury_delta_aging_2026.rds")

mouse_sample_delta <- aggregate(cbind(predicted, delta_aging) ~ sample + group + celltype,
                                data = Mouse_delta, FUN = mean)
mouse_sample_delta <- mouse_sample_delta[order(mouse_sample_delta$celltype, mouse_sample_delta$group), ]
message("\nMouse - Sample Level Delta Aging Score"); print(mouse_sample_delta, row.names = FALSE)
write.csv(mouse_sample_delta, "Mouse_injury_sample_delta_aging_2026.csv", row.names = FALSE)


########################################################################
# 7) Plots: delta by celltype + faceted by timepoint (Nature/NPG palette, large fonts)
########################################################################
npg_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
              "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# (a) Per celltype: delta boxplot of Control vs Injury
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
message("Saved: Mouse_Injury_Delta_Aging_Boxplot_2026.png")

# (b) Timepoint facets: one row per celltype, x = injury timepoint
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
# For plotting, pool cells from 6h–36h (06/12/24/36h) into a single group NMDA_6-36h
Mouse_delta_plot <- Mouse_delta
Mouse_delta_plot$sample[Mouse_delta_plot$sample %in%
    c("NMDA_06h", "NMDA_12h", "NMDA_24h", "NMDA_36h")] <- "NMDA_6-36h"

# x-axis order (NMDA_0h is the control; 6-36h pooled; mmNMDA48FI not plotted)
mouse_sample_order <- c("NMDA_0h", "NMDA_03h", "NMDA_6-36h", "NMDA_48h", "NMDA_72h")
mouse_ct_order     <- c("MG", "Rod", "Cone")
ggsave("Mouse_Injury_Timepoint_Boxplot_2026.png",
       plot_timepoint_facet(Mouse_delta_plot, mouse_sample_order, mouse_ct_order),
       width = 4.5, height = length(mouse_ct_order) * 3 + 1, dpi = 300)
message("Saved: Mouse_Injury_Timepoint_Boxplot_2026.png")

message("\n=== Mouse injury aging clock (2026) all done ===")



####### 
setwd(injury_dir)
saveRDS(Mouse_injury_predictions, file = "Mouse_injury_predictions_2026.rds")
message("Saved: Mouse_injury_predictions_2026.rds")


########################################################################
# 8) Compare predicted age between 48h (mmNMDA48) and 48FIh (mmNMDA48FI)
#    boxplot + Nature/NPG palette + Wilcoxon p-value
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

# Wilcoxon p per celltype (48h vs 48FIh)
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
message("Saved: Mouse_Injury_48h_vs_48FIh_boxplot_2026.png")


#####
##### As in zebrafish, first compute injury DEGs, with NMDA_6-36h pooled for DEG calculation #####
#####
########################################################################
# 9) Aging-related injury genes (Mouse):
#    (a) Find injury DEGs on unsmoothed raw counts (Mouse_injury_sel)
#        - MG restricted to Resting MG; mmNMDA06/12/24/36 pooled as mmNMDA6-36h; mmNMDA48FI excluded
#    (b) Extract clock gene weights from the 2026 final models
#    (c) Injury DEG that is also a clock gene = aging-related injury gene
########################################################################
library(Seurat)
library(Matrix)

# Unsmoothed raw object: Mouse_injury_sel (generated in Section 3); unify sample=Condition
mouse_deg_input <- Mouse_injury_sel
for (ct in names(mouse_deg_input)) mouse_deg_input[[ct]]$sample <- as.character(mouse_deg_input[[ct]]$Condition)
# Keep only Resting MG for MG (consistent with clock/heatmap definition)
if ("MG" %in% names(mouse_deg_input)) {
    mg <- mouse_deg_input$MG
    mouse_deg_input$MG <- mg[, which(mg$celltype == "Resting MG")]
}
# Pool 6-36h (mmNMDA06/12/24/36 -> mmNMDA6-36h) and exclude mmNMDA48FI
merge6_36 <- c("mmNMDA06", "mmNMDA12", "mmNMDA24", "mmNMDA36")
for (ct in names(mouse_deg_input)) {
    s <- mouse_deg_input[[ct]]$sample
    s[s %in% merge6_36] <- "mmNMDA6-36h"
    mouse_deg_input[[ct]]$sample <- s
    mouse_deg_input[[ct]] <- mouse_deg_input[[ct]][, s != "mmNMDA48FI"]
}

# (a) Injury DEGs: each timepoint vs mmP60 (p_val_adj<0.01, |log2FC|>0.2)
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
                error = function(e) { message("    FindMarkers error: ", e$message); NULL })
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

# (b) Extract clock gene weights from the 2026 final models (MG/Rod/Cone)
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

# (c) Merge DEGs with clock weights and flag is_clock_gene
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

# Aging-related injury genes = injury DEGs that are also clock genes
mouse_aging_injury_genes <- mouse_deg_table[mouse_deg_table$is_clock_gene %in% TRUE, , drop = FALSE]
write.csv(mouse_aging_injury_genes, "Mouse_aging_related_injury_genes_2026.csv", row.names = FALSE)
message(sprintf("\nAging-related injury genes: %d rows, %d unique genes",
                nrow(mouse_aging_injury_genes), length(unique(mouse_aging_injury_genes$gene))))


##############
############## Next, split mouse_aging_injury_genes into positive/negative groups by weight #####
############## Filter by Spearman correlation between injury z-score and predicted aging at each timepoint #####
############## MG only; heatmap Injury left | Aging right; pos on top / neg below #####
########################################################################
# 10) MG: build injury/aging pseudobulk + Spearman filtering + ComplexHeatmap
########################################################################
library(Matrix)
library(ComplexHeatmap)
library(circlize)

ct_use <- "MG"

# pseudobulk: sum of raw counts -> depth norm 1e4 -> log2(1+x)
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

# --- (a) Injury pseudobulk (smoothed MG; 6-36h pooled; 48FI excluded; column names match the boxplot) ---
# Build pb from KNN-smoothed injury data (same as prediction); read RDS if not in memory
setwd(injury_dir)
if (!exists("Mouse_injury_smoothed")) {
    Mouse_injury_smoothed <- readRDS("Mouse_injury_smoothed.rds")
    message("Loaded: Mouse_injury_smoothed.rds")
}
setwd(mouse_dir)
load("Mouse_ambient_genes_Top50.rda")   # -> M_ambient_genes_Top50
inj_pb <- Make_Pseudobulk_Matrix(Mouse_injury_smoothed[[ct_use]], M_ambient_genes_Top50, "sample")
# Remove 48FI
inj_pb <- inj_pb[, colnames(inj_pb) != "mmNMDA48FI", drop = FALSE]
# Pool 6-36h
merge6_36 <- c("mmNMDA06", "mmNMDA12", "mmNMDA24", "mmNMDA36")
present6  <- intersect(merge6_36, colnames(inj_pb))
if (length(present6) > 0) {
    merged_col <- rowMeans(inj_pb[, present6, drop = FALSE])
    keep_cols  <- setdiff(colnames(inj_pb), merge6_36)
    inj_pb     <- cbind(inj_pb[, keep_cols, drop = FALSE], `mmNMDA6-36h` = merged_col)
}
# Rename to the boxplot labels
rename_inj <- c(mmP60 = "NMDA_0h", mmNMDA03 = "NMDA_03h", `mmNMDA6-36h` = "NMDA_6-36h",
                mmNMDA48 = "NMDA_48h", mmNMDA72 = "NMDA_72h")
cn <- colnames(inj_pb)
cn <- ifelse(cn %in% names(rename_inj), rename_inj[cn], cn)
colnames(inj_pb) <- cn
mouse_injury_col_order <- c("NMDA_0h", "NMDA_03h", "NMDA_6-36h", "NMDA_48h", "NMDA_72h")
inj_pb <- inj_pb[, mouse_injury_col_order[mouse_injury_col_order %in% colnames(inj_pb)], drop = FALSE]

# --- (b) Aging/training pseudobulk (clock training data; merge replicates / Wk106-108) ---
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

# --- (c) Split aging-related injury genes into two groups by clock weight sign (MG only) ---
if (!exists("mouse_aging_injury_genes")) {
    setwd(injury_dir)
    mouse_aging_injury_genes <- read.csv("Mouse_aging_related_injury_genes_2026.csv",
                                         stringsAsFactors = FALSE)
    message("Loaded: Mouse_aging_related_injury_genes_2026.csv")
}
genes_ct  <- mouse_aging_injury_genes[mouse_aging_injury_genes$celltype == ct_use, ]
genes_ct  <- genes_ct[!duplicated(genes_ct$gene), ]
genes_pos <- genes_ct$gene[genes_ct$clock_coefficient > 0]
genes_neg <- genes_ct$gene[genes_ct$clock_coefficient < 0]
message(sprintf("%s: aging injury genes = %d (pos %d / neg %d)",
                ct_use, nrow(genes_ct), length(genes_pos), length(genes_neg)))

coef_lookup <- setNames(genes_ct$clock_coefficient, genes_ct$gene)

# --- (d) Predicted aging for each injury timepoint (MG), aligned to inj_pb columns (with 6-36h pooled) ---
if (!exists("Mouse_delta")) {
    setwd(injury_dir)
    Mouse_delta <- readRDS("Mouse_injury_delta_aging_2026.rds")
    message("Loaded: Mouse_injury_delta_aging_2026.rds")
}
md <- Mouse_delta[Mouse_delta$celltype == ct_use, ]
md$sample[md$sample %in% c("NMDA_06h", "NMDA_12h", "NMDA_24h", "NMDA_36h")] <- "NMDA_6-36h"
md <- md[!md$sample %in% "NMDA_48FIh", ]
pred_map <- setNames(aggregate(predicted ~ sample, md, mean)$predicted,
                     aggregate(predicted ~ sample, md, mean)$sample)
pred_vec <- pred_map[colnames(inj_pb)]
message("Predicted aging (by injury timepoint): ",
        paste(sprintf("%s=%.1f", colnames(inj_pb), pred_vec), collapse = ", "))

# --- (e) Spearman: gene injury z-score vs predicted aging ---
inj_z_all <- zscore_rows(inj_pb)
spearman_vs_pred <- function(genes) {
    g <- intersect(genes, rownames(inj_z_all))
    setNames(sapply(g, function(gene)
        suppressWarnings(cor(inj_z_all[gene, ], pred_vec, method = "spearman"))), g)
}
cor_pos <- spearman_vs_pred(genes_pos)
cor_neg <- spearman_vs_pred(genes_neg)

# --- (f) Cutoff filtering ---
cor_cutoff <- 0.25
sel_pos <- names(cor_pos)[!is.na(cor_pos) & cor_pos >  cor_cutoff]
sel_neg <- names(cor_neg)[!is.na(cor_neg) & cor_neg < -cor_cutoff]
sel_pos <- sel_pos[order(-cor_pos[sel_pos])]
sel_neg <- sel_neg[order( cor_neg[sel_neg])]
message(sprintf("Positive weights: %d/%d passed (cor > %.2f); negative weights: %d/%d passed (cor < %.2f)",
                length(sel_pos), length(cor_pos), cor_cutoff,
                length(sel_neg), length(cor_neg), -cor_cutoff))

# --- (g) Combined heatmap: pos on top / neg below; Injury left / Aging right; clock weight at far right ---
setwd(injury_dir)
sel_genes <- c(sel_pos, sel_neg)
if (length(sel_genes) < 2) {
    message("Too few selected genes, skipping heatmap")
} else {
    sel_genes <- Reduce(intersect, list(sel_genes, rownames(inj_pb), rownames(trn_pb)))
    sel_pos   <- intersect(sel_pos, sel_genes)
    sel_neg   <- intersect(sel_neg, sel_genes)
    sel_pos   <- sel_pos[order(-coef_lookup[sel_pos])]
    sel_neg   <- sel_neg[order(-coef_lookup[sel_neg])]
    sel_genes <- c(sel_pos, sel_neg)

    inj_z <- zscore_rows(inj_pb[sel_genes, , drop = FALSE])
    trn_z <- zscore_rows(trn_pb[sel_genes, , drop = FALSE])
    # Right-side Aging column names: show age in weeks only; Wk106-108 -> 108
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
    # Same as zebrafish: anno_barplot
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

    # Sizes matched to zebrafish
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
    message(sprintf("Saved: %s (pos=%d, neg=%d)", out_file, length(sel_pos), length(sel_neg)))
}

setwd(injury_dir)
sel_df <- rbind(
    data.frame(gene = sel_pos, spearman = cor_pos[sel_pos],
               clock_coefficient = coef_lookup[sel_pos], group = "posWeight_posCorr"),
    data.frame(gene = sel_neg, spearman = cor_neg[sel_neg],
               clock_coefficient = coef_lookup[sel_neg], group = "negWeight_negCorr")
)
write.csv(sel_df, "Mouse_MG_aging_injury_selected_genes_2026.csv", row.names = FALSE)
message("Saved: ", file.path(injury_dir, "Mouse_MG_aging_injury_selected_genes_2026.csv"))

########################################################################
# 11) Mouse MG → Human ortholog → enrichGO BP (no KEGG; top 10)
# Pipeline matches Zebrafish; new terminal: conda activate clusterProfiler
########################################################################
if (!exists("injury_dir")) injury_dir <- "/projects/hmz-aging/Injury_Datasets"
if (!exists("mouse_dir"))  mouse_dir  <- "/projects/hmz-aging/Mouse_New_Datasets"

suppressMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})

setwd(injury_dir)

# --- Always read from the Mouse CSV (do not reuse Zebrafish sel_pos/sel_neg from the same session) ---
sel_csv_candidates <- c(
    file.path(injury_dir, "Mouse_MG_aging_injury_selected_genes_2026.csv"),
    file.path(mouse_dir,  "Mouse_MG_aging_injury_selected_genes_2026.csv"),
    "Mouse_MG_aging_injury_selected_genes_2026.csv"
)
sel_csv <- sel_csv_candidates[file.exists(sel_csv_candidates)][1]
if (is.na(sel_csv) || !nzchar(sel_csv))
    stop("Mouse_MG_aging_injury_selected_genes_2026.csv not found; run the gene selection section first")
mouse_sel_df <- read.csv(sel_csv, stringsAsFactors = FALSE)
mouse_sel_pos <- unique(mouse_sel_df$gene[mouse_sel_df$group == "posWeight_posCorr"])
mouse_sel_neg <- unique(mouse_sel_df$gene[mouse_sel_df$group == "negWeight_negCorr"])
mouse_sel_pos <- mouse_sel_pos[!is.na(mouse_sel_pos) & mouse_sel_pos != ""]
mouse_sel_neg <- mouse_sel_neg[!is.na(mouse_sel_neg) & mouse_sel_neg != ""]
message(sprintf("Loaded Mouse genes from CSV: %s (pos=%d, neg=%d)",
                sel_csv, length(mouse_sel_pos), length(mouse_sel_neg)))
message("  Example pos: ", paste(head(mouse_sel_pos, 5), collapse = ", "))

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
    message("Loaded ortholog: ", fp, " (HM pairs=", nrow(hm), ")")
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

# enrichGO BP; take top 10 by p.adjust
run_go_human <- function(mouse_genes, group_name, p_cut = 0.05, q_cut = 0.2, top_n = 10) {
    message(sprintf("\n========== %s: mouse genes=%d ==========", group_name, length(mouse_genes)))
    map_df <- mouse_to_human(mouse_genes)
    human_genes <- unique(map_df$human)
    message(sprintf("  mapped to human: %d / %d", length(human_genes), length(unique(mouse_genes))))
    if (length(human_genes) < 5) { message("  Too few human genes, skipping"); return(NULL) }

    entrez <- na.omit(AnnotationDbi::mapIds(
        org.Hs.eg.db, keys = human_genes, column = "ENTREZID",
        keytype = "SYMBOL", multiVals = "first"
    ))
    entrez <- unique(as.character(entrez))
    message(sprintf("  human EntrezID: %d", length(entrez)))
    if (length(entrez) < 5) { message("  Too few Entrez IDs, skipping"); return(NULL) }

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
        message("  No GO BP results"); return(NULL)
    }
    d <- as.data.frame(ego)
    d <- d[d$pvalue < p_cut, , drop = FALSE]
    if (nrow(d) == 0) { message("  No significant GO terms"); return(NULL) }

    d$GeneRatio_val  <- vapply(d$GeneRatio, parse_ratio, numeric(1))
    d$BgRatio_val    <- vapply(d$BgRatio, parse_ratio, numeric(1))
    d$FoldEnrichment <- d$GeneRatio_val / d$BgRatio_val
    d$neglog10_padj  <- -log10(pmax(d$p.adjust, 1e-300))
    d$Database <- "GO:BP"
    d$group <- group_name
    d <- d[order(d$p.adjust, d$pvalue), ]
    message(sprintf("  GO BP sig: %d (plotting top %d)", nrow(d), top_n))

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
    message("Saved: Mouse_MG_aging_injury_GO_humanOrtholog_2026.csv")
}
if (length(map_tabs) > 0) {
    write.csv(do.call(rbind, map_tabs),
              "Mouse_MG_aging_injury_mouse2human_map_2026.csv", row.names = FALSE)
    message("Saved: Mouse_MG_aging_injury_mouse2human_map_2026.csv")
}

message("\n=== Section 11 done (Mouse→Human GO BP, top 10) ===")
message("Next step: run Section 12 to draw the dotplot")


########################################################################
# 12) Publication-quality dotplot: top 10 GO BP terms (one each for pos / neg; Nature/NPG palette)
#     x = GeneRatio, size = Count, color = p.adjust (same as Zebrafish)
########################################################################
if (!exists("injury_dir")) injury_dir <- "/projects/hmz-aging/Injury_Datasets"
setwd(injury_dir)

suppressMessages(library(ggplot2))

enrich_rds <- "Mouse_MG_aging_injury_GO_humanOrtholog_2026.rds"
if (!exists("mouse_enrich_hs") || is.null(mouse_enrich_hs)) {
    if (!file.exists(enrich_rds)) stop(enrich_rds, " not found; run Section 11 first")
    mouse_enrich_hs <- readRDS(enrich_rds)
    message("Loaded: ", enrich_rds)
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
        message(nm, ": no GO terms to plot"); next
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
    message("Saved: ", out_png, " (600 dpi)")
}

message("\n=== Mouse GO Dotplot (human ortholog, top 10) done ===")
