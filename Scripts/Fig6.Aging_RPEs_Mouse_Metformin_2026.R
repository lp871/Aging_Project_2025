########
######## metformin bulk mouse RPE ##########
########

conda activate seurat_env
R

#########
#########

setwd("/projects/hmz-aging/Metformin/")

#########
#########
library(SummarizedExperiment)


#########
#########

bharti_human <- readRDS("mitra_metformin_bulkRNA_iRPE.rds")
bharti_mouse <- readRDS("mitra_metformin_bulkRNA_mouseRPE.rds")

#########
#########
#########

colData(bharti_mouse)
meta = colData(bharti_mouse)
table(meta$Genotype_Treatment_Tissue)
assay(bharti_mouse)[1:10,1:10]

#########
#########

meta = colData(bharti_human)
table(meta$treatment)
assay(bharti_human)[1:10,1:10]


######## PCA - bharti_mouse ########

library(DESeq2)
library(ggplot2)

counts_mouse <- as.matrix(assay(bharti_mouse))
meta_mouse <- as.data.frame(colData(bharti_mouse))
cat("Dim:", dim(counts_mouse), "\n")

# # DESeq2 归一化 + variance stabilizing transformation
# dds <- DESeqDataSetFromMatrix(countData = counts_mouse,
#                                colData = meta_mouse,
#                                design = ~ 1)
# vsd <- vst(dds, blind = TRUE)

# # 去掉 0 方差基因后做 PCA
# vsd_mat <- assay(vsd)
# gene_var <- apply(vsd_mat, 1, var)
# vsd_mat <- vsd_mat[gene_var > 0, ]
# pca_res <- prcomp(t(vsd_mat), scale. = TRUE)
# pca_df <- data.frame(
#   PC1 = pca_res$x[, 1],
#   PC2 = pca_res$x[, 2],
#   group = meta_mouse$Genotype_Treatment_Tissue
# )

# var_pct <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

# p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#     aspect.ratio = 1,
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 11, face = "bold")
#   ) +
#   labs(x = paste0("PC1 (", var_pct[1], "%)"),
#        y = paste0("PC2 (", var_pct[2], "%)"),
#        color = "Group",
#        title = "PCA - Mouse Bulk RNA-seq")

# ggsave("bharti_mouse_PCA.png", plot = p_pca, width = 8, height = 6, dpi = 300)

# # 分 Retina 和 RPE 单独做 PCA
# meta_mouse$tissue <- ifelse(grepl("Retina", meta_mouse$Genotype_Treatment_Tissue), "Retina", "RPE_Choroid")

# # Retina samples
# retina_idx <- which(meta_mouse$tissue == "Retina")
# counts_retina <- counts_mouse[, retina_idx]
# meta_retina <- meta_mouse[retina_idx, ]
# dds_retina <- DESeqDataSetFromMatrix(countData = counts_retina, colData = meta_retina, design = ~ 1)
# vsd_retina <- vst(dds_retina, blind = TRUE)
# vsd_retina_mat <- assay(vsd_retina)
# vsd_retina_mat <- vsd_retina_mat[apply(vsd_retina_mat, 1, var) > 0, ]
# pca_retina <- prcomp(t(vsd_retina_mat), scale. = TRUE)
# var_retina <- round(100 * summary(pca_retina)$importance[2, 1:2], 1)

# pca_retina_df <- data.frame(
#   PC1 = pca_retina$x[, 1], PC2 = pca_retina$x[, 2],
#   group = meta_retina$Genotype_Treatment_Tissue
# )

# p_retina <- ggplot(pca_retina_df, aes(x = PC1, y = PC2, color = group)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#         aspect.ratio = 1,
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 11, face = "bold")) +
#   labs(x = paste0("PC1 (", var_retina[1], "%)"),
#        y = paste0("PC2 (", var_retina[2], "%)"),
#        color = "Group", title = "PCA - Retina")

# ggsave("bharti_mouse_PCA_Retina.png", plot = p_retina, width = 8, height = 6, dpi = 300)

# RPE/Choroid samples
rpe_idx <- which(meta_mouse$Tissue == "RPE/Choroid")
all.equal(colnames(counts_mouse),rownames(meta_mouse))
counts_rpe <- counts_mouse[, rpe_idx]
meta_rpe <- meta_mouse[rpe_idx, ]
dds_rpe <- DESeqDataSetFromMatrix(countData = counts_rpe, colData = meta_rpe, design = ~ 1)
vsd_rpe <- vst(dds_rpe, blind = TRUE)
vsd_rpe_mat <- assay(vsd_rpe)
vsd_rpe_mat <- vsd_rpe_mat[apply(vsd_rpe_mat, 1, var) > 0, ]


# pca_rpe <- prcomp(t(vsd_rpe_mat), scale. = TRUE)
# var_rpe <- round(100 * summary(pca_rpe)$importance[2, 1:2], 1)

# pca_rpe_df <- data.frame(
#   PC1 = pca_rpe$x[, 1], PC2 = pca_rpe$x[, 2],
#   group = meta_rpe$Genotype_Treatment_Tissue
# )

# p_rpe <- ggplot(pca_rpe_df, aes(x = PC1, y = PC2, color = group)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#         aspect.ratio = 1,
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 11, face = "bold")) +
#   labs(x = paste0("PC1 (", var_rpe[1], "%)"),
#        y = paste0("PC2 (", var_rpe[2], "%)"),
#        color = "Group", title = "PCA - RPE/Choroid")

# ggsave("bharti_mouse_PCA_RPE.png", plot = p_rpe, width = 8, height = 6, dpi = 300)

# # UMAP based on PCA (99% variance)
# library(uwot)

# cum_var <- cumsum(summary(pca_res)$importance[2, ])
# n_pcs <- which(cum_var >= 0.99)[1]
# cat("Using", n_pcs, "PCs for 99% variance\n")

# umap_res <- umap(pca_res$x[, 1:n_pcs], n_neighbors = 3, min_dist = 0.1)

# umap_df <- data.frame(
#   UMAP1 = umap_res[, 1],
#   UMAP2 = umap_res[, 2],
#   group = meta_mouse$Genotype_Treatment_Tissue
# )

# p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#     aspect.ratio = 1,
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 11, face = "bold")
#   ) +
#   labs(color = "Group", title = paste0("UMAP (", n_pcs, " PCs, 99% variance)"))

# ggsave("bharti_mouse_UMAP.png", plot = p_umap, width = 8, height = 6, dpi = 300)


######## 用 Aging DEGs 计算 bulk RPE 的 aging score ########

# # 读入 aging DEGs
# aging_degs <- read.csv("/projects/hmz-aging/RPE_smoking/RPE_aging_consistent_DEGs.csv")
# up_aging <- aging_degs$gene_symbol[aging_degs$direction == "Up"]
# down_aging <- aging_degs$gene_symbol[aging_degs$direction == "Down"]

# # 只用 RPE/Choroid 样本的 VST 归一化数据
# # vsd_rpe_mat 已经在前面 PCA 时计算好了
# rpe_expr <- vsd_rpe_mat

# # 提取基因名（如果 rownames 有 ENSEMBL 前缀，需要转换）
# gene_names <- rownames(rpe_expr)
# cat("Example gene names:", head(gene_names, 5), "\n")

# # 匹配 aging DEGs
# up_found <- intersect(up_aging, gene_names)
# down_found <- intersect(down_aging, gene_names)
# cat("Up genes found:", length(up_found), "/", length(up_aging), "\n")
# cat("Down genes found:", length(down_found), "/", length(down_aging), "\n")

# # 计算 aging score: 每个样本的 z-score 平均值
# # 先对每个基因做 z-score（跨样本标准化）
# rpe_z <- t(scale(t(rpe_expr)))

# # Aging score = mean(z-score of up genes) - mean(z-score of down genes)
# up_score <- colMeans(rpe_z[up_found, , drop = FALSE], na.rm = TRUE)
# down_score <- colMeans(rpe_z[down_found, , drop = FALSE], na.rm = TRUE)
# aging_score_bulk <- up_score - down_score

# # 打印每个样本的 aging score
# score_df <- data.frame(
#   sample = names(aging_score_bulk),
#   group = meta_rpe$Genotype_Treatment_Tissue,
#   aging_score = aging_score_bulk
# )
# print(score_df[order(score_df$group), ])

# # 画 barplot
# library(ggplot2)
# p_score <- ggplot(score_df, aes(x = reorder(sample, aging_score), y = aging_score, fill = group)) +
#   geom_col() +
#   coord_flip() +
#   theme_classic() +
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#         legend.text = element_text(size = 10)) +
#   labs(x = "", y = "Aging Score", fill = "Group", title = "Bulk RPE Aging Score")

# ggsave("bharti_mouse_RPE_aging_score.png", plot = p_score, width = 8, height = 5, dpi = 300)

######## ssGSEA aging score ########

rpe_expr <- vsd_rpe_mat

# 读入 aging DEGs
aging_degs <- read.csv("/projects/hmz-aging/RPE_smoking/RPE_aging_consistent_DEGs.csv")
up_aging <- aging_degs$gene_symbol[aging_degs$direction == "Up"]
down_aging <- aging_degs$gene_symbol[aging_degs$direction == "Down"]

# 只用 RPE/Choroid 样本的 VST 归一化数据
# vsd_rpe_mat 已经在前面 PCA 时计算好了
rpe_expr <- vsd_rpe_mat

# 提取基因名（如果 rownames 有 ENSEMBL 前缀，需要转换）
gene_names <- rownames(rpe_expr)
cat("Example gene names:", head(gene_names, 5), "\n")

# 匹配 aging DEGs
up_found <- intersect(up_aging, gene_names)
down_found <- intersect(down_aging, gene_names)
cat("Up genes found:", length(up_found), "/", length(up_aging), "\n")
cat("Down genes found:", length(down_found), "/", length(down_aging), "\n")


library(GSVA)

gs <- list(Up_aging = up_found, Down_aging = down_found)
param <- ssgseaParam(exprData = rpe_expr, geneSets = gs, alpha = 0)
ssgsea_res <- gsva(param, verbose = FALSE)

aging_ssgsea <- ssgsea_res["Up_aging", ] - ssgsea_res["Down_aging", ]

score_ssgsea <- data.frame(
  sample = names(aging_ssgsea),
  group = meta_rpe$Genotype_Treatment_Tissue,
  aging_score = aging_ssgsea
)
print(score_ssgsea[order(score_ssgsea$group), ])

# 400mg vs 其他样品的 aging score 比较
score_ssgsea$is_400mg <- grepl("400 mg", score_ssgsea$group)

met_scores <- score_ssgsea$aging_score[score_ssgsea$is_400mg]
other_scores <- score_ssgsea$aging_score[!score_ssgsea$is_400mg]

cat("400mg scores:", met_scores, "\n")
cat("Other scores:", other_scores, "\n")
cat("400mg mean:", mean(met_scores), "| Other mean:", mean(other_scores), "\n")

# Wilcoxon test (单侧: 400mg 是否显著低于其他)
wt <- wilcox.test(met_scores, other_scores, alternative = "less")
cat("Wilcox p-value (400mg < others):", wt$p.value, "\n")

# t-test
tt <- t.test(met_scores, other_scores, alternative = "less")
cat("t-test p-value (400mg < others):", tt$p.value, "\n")
cat("t-test 95% CI upper:", tt$conf.int[2], "\n")



####
#### sample 1，2 #####
####
####

#### 散点图: Metformin- vs Metformin+ ####

score_ssgsea$condition <- ifelse(grepl("400 mg", score_ssgsea$group), "Metformin+", "Metformin-")
score_ssgsea$condition <- factor(score_ssgsea$condition, levels = c("Metformin-", "Metformin+"))

# t-test: Metformin+ < Metformin-
pval <- t.test(aging_score ~ condition, data = score_ssgsea, alternative = "greater")$p.value
pval_label <- paste0("p = ", format(pval, digits = 3))

y_max <- max(score_ssgsea$aging_score) + 0.05 * diff(range(score_ssgsea$aging_score))

group_colors_met <- c(
  "ABCA4KO_400 mg_RPE/Choroid" = "#E64B35",
  "ABCA4KO_untreated_RPE/Choroid" = "#4DBBD5",
  "WT_untreated_RPE/Choroid" = "#00A087"
)

p_dot <- ggplot(score_ssgsea, aes(x = condition, y = aging_score, color = group)) +
  geom_jitter(size = 4, width = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, linewidth = 0.6, color = "black") +
  scale_color_manual(values = group_colors_met) +
  geom_segment(aes(x = 1, xend = 2, y = y_max, yend = y_max), color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  geom_text(aes(x = 1.5, y = y_max + 0.07 * diff(range(score_ssgsea$aging_score)), label = pval_label),
            inherit.aes = FALSE, size = 3.5) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1.2,
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  ) +
  labs(x = "", y = "Aging Score (ssGSEA)", title = "RPE Aging Score")

ggsave("metformin_aging_score_dotplotXXX.png", plot = p_dot, width = 6, height = 7, dpi = 300)



########----------------------------------------------------------########
########----------------------------------------------------------########
########----------------------------------------------------------########



#########
######### Next is the DEGs heatmaps ##########----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########

#### Metformin DEGs: treated (2) vs untreated (4) ####

library(limma)

# 分组: treated vs untreated
group_met <- ifelse(grepl("400 mg", meta_rpe$Genotype_Treatment_Tissue), "treated", "untreated")
group_met <- factor(group_met, levels = c("untreated", "treated"))

library(edgeR)
dge <- DGEList(counts = counts_rpe)
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ group_met)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# 最宽松条件: 不设 logFC 和 p 值阈值
met_degs <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
met_degs$gene <- rownames(met_degs)

cat("Total genes tested:", nrow(met_degs), "\n")
cat("p < 0.05:", sum(met_degs$P.Value < 0.05), "\n")
cat("p < 0.1:", sum(met_degs$P.Value < 0.1), "\n")
cat("adj.P < 0.1:", sum(met_degs$adj.P.Val < 0.1), "\n")

# 用 p < 0.1 作为宽松阈值
met_degs_sig <- met_degs[met_degs$P.Value < 0.1, ]
cat("Sig DEGs (p<0.1):", nrow(met_degs_sig), "\n")
cat("  Up:", sum(met_degs_sig$logFC > 0), "\n")
cat("  Down:", sum(met_degs_sig$logFC < 0), "\n")

#### 与 aging DEGs 比较: 找方向相反的 ####
####

aging_degs <- read.csv("/projects/hmz-aging/RPE_smoking/RPE_aging_consistent_DEGs.csv")

# merge
overlap_met <- merge(met_degs_sig, aging_degs, by.x = "gene", by.y = "gene_symbol")

# 方向相反: metformin 下调 + aging 上调, 或 metformin 上调 + aging 下调
# => metformin 逆转了衰老
overlap_met$met_direction <- ifelse(overlap_met$logFC > 0, "Up", "Down")
opposite_met <- overlap_met[overlap_met$met_direction != overlap_met$direction, ]
same_met <- overlap_met[overlap_met$met_direction == overlap_met$direction, ]

cat("\nOverlap total:", nrow(overlap_met), "\n")
cat("Opposite direction (metformin reverses aging):", nrow(opposite_met), "\n")
cat("  Aging Up + Met Down:", sum(opposite_met$direction == "Up"), "\n")
cat("  Aging Down + Met Up:", sum(opposite_met$direction == "Down"), "\n")
cat("Same direction:", nrow(same_met), "\n")

print(opposite_met[, c("gene", "logFC", "P.Value", "direction", "met_direction")])

write.csv(opposite_met, "Metformin_reverses_aging_genes.csv", row.names = FALSE)



#######
####### Heatmap: Metformin consistent genes + Aging ####
#######

library(ComplexHeatmap)
library(circlize)
library(Seurat)

reversed_up = opposite_met$gene[opposite_met$direction == "Up"]
reversed_down = opposite_met$gene[opposite_met$direction == "Down"]

# 只取方向和 aging 相反的基因
# reversed_up: met up + aging down
# reversed_down: met down + aging up
all_met_genes <- c(reversed_up, reversed_down)
all_met_genes <- all_met_genes[all_met_genes %in% rownames(expr)]
cat("Reversed genes to plot:", length(all_met_genes), "\n")

# 左边: Metformin bulk samples, 按 WT → KO_untreated → KO_400mg 排序
col_order <- c(wt_samples, ko_samples, met_samples)
met_mat <- expr[all_met_genes, col_order]
met_z <- t(scale(t(met_mat)))

# 行分组: Aging Up (met down) / Aging Down (met up)
met_row_dir <- ifelse(all_met_genes %in% reversed_down, "Aging Up\n(Met Down)", "Aging Down\n(Met Up)")
met_row_split <- factor(met_row_dir, levels = c("Aging Up\n(Met Down)", "Aging Down\n(Met Up)"))

# 列注释
col_group <- c(rep("WT_untreated", length(wt_samples)),
               rep("ABCA4KO_untreated", length(ko_samples)),
               rep("ABCA4KO_400mg", length(met_samples)))
col_anno_met <- HeatmapAnnotation(
  Group = col_group,
  col = list(Group = c("WT_untreated" = "#4DBBD5", "ABCA4KO_untreated" = "#00A087", "ABCA4KO_400mg" = "#E64B35")),
  annotation_name_gp = gpar(fontsize = 10)
)

col_fun_left <- colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))

# 右边: Aging RPE pseudobulk (重新计算)
Mouse_RPE_rpe <- Mouse_RPE_aging[, Mouse_RPE_aging$celltype %in% c("RPE")]
Mouse_RPE_rpe$sample <- factor(Mouse_RPE_rpe$sample, levels = c("5W_RPE", "49W_RPE", "91W_RPE"))
counts_ag <- GetAssayData(Mouse_RPE_rpe, assay = "RNA", layer = "counts")
sample_ag <- Mouse_RPE_rpe$sample
bulk_ag <- sapply(levels(sample_ag), function(s) {
  Matrix::rowSums(counts_ag[, sample_ag == s, drop = FALSE])
})
lib_ag <- colSums(bulk_ag)
mat_aging <- log2(t(t(bulk_ag) / lib_ag * 1e6) + 1)

aging_gene_map2 <- setNames(sub(".*~~", "", rownames(mat_aging)), rownames(mat_aging))
aging_rows2 <- names(aging_gene_map2)[aging_gene_map2 %in% all_met_genes]
aging_mat2 <- mat_aging[aging_rows2, c("5W_RPE", "49W_RPE", "91W_RPE"), drop = FALSE]
rownames(aging_mat2) <- aging_gene_map2[aging_rows2]
aging_z2 <- t(scale(t(aging_mat2)))

common_final <- intersect(all_met_genes, rownames(aging_z2))
show_names <- length(common_final) <= 100
met_row_dir2 <- ifelse(common_final %in% reversed_down, "Aging Up\n(Met Down)", "Aging Down\n(Met Up)")
met_row_split2 <- factor(met_row_dir2, levels = c("Aging Up\n(Met Down)", "Aging Down\n(Met Up)"))

ht_left <- Heatmap(met_z[common_final, ],
                   name = "Metformin\nZ-score", col = colorRamp2(c(-2, 0, 2), c("#00A087", "white", "#DC0000")),
                   rect_gp = gpar(col = "white", lwd = 0.5),
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   row_split = met_row_split2, row_gap = unit(3, "mm"),
                   row_title = c("Aging Up\n(Met Down)", "Aging Down\n(Met Up)"),
                   row_title_gp = gpar(fontsize = 11),
                   column_title = "Metformin RPE", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                   top_annotation = col_anno_met,
                   show_row_names = show_names, row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                   row_names_side = "left",
                   column_names_rot = 45, column_names_gp = gpar(fontsize = 9))

ht_right <- Heatmap(aging_z2[common_final, ],
                    name = "Aging\nZ-score", col = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")),
                    rect_gp = gpar(col = "white", lwd = 0.5),
                    cluster_columns = FALSE, cluster_rows = FALSE,
                    row_split = met_row_split2, row_gap = unit(3, "mm"),
                    column_title = "Aging RPE (5W→49W→91W)", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                    show_row_names = FALSE,
                    column_names_rot = 45, column_names_gp = gpar(fontsize = 9))

ht_height <- min(40, max(6, length(common_final) * 0.15 + 3))

cat("Heatmap genes:", length(common_final), "| Height:", ht_height, "in\n")
show_names <- length(common_final) <= 100

png("Metformin_consistent_heatmapXXX.png", width = 8, height = ht_height, units = "in", res = 300)
draw(ht_left + ht_right, ht_gap = unit(8, "mm"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     merge_legend = TRUE)
dev.off()



#######
#######
#######


load(file.path("/projects/hmz-aging/RPE_smoking", "Handa_snRNA_Final 1.RData"))
Mouse_RPE_aging <- readRDS(file.path("/projects/hmz-aging/RPE_smoking", "Mouse_RPE_RNA_Final_Jul16"))
Mouse_RPE_rpe <- Mouse_RPE_aging[, Mouse_RPE_aging$celltype %in% c("RPE")]
Mouse_RPE_rpe$sample <- factor(Mouse_RPE_rpe$sample, levels = c("5W_RPE", "49W_RPE", "91W_RPE"))




#######