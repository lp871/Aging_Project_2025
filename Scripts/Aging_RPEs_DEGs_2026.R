########
########
########

conda activate seurat_env
R


library(Seurat)
setwd("/projects/hmz-aging/RPE_smoking")
load(file.path("Handa_snRNA_Final 1.RData"))
Mouse_RPE <- readRDS(file.path("Mouse_RPE_RNA_Final_Jul16"))



#####
#####----------------------------------------------------------#####
#####

table(Mouse_RPE$sample)
table(Mouse_RPE$celltype)

Mouse_RPE$celltype <- dplyr::recode(Mouse_RPE$celltype,
  "RPE" = "Retinal pigment epithelial cells",
  "Stromal" = "Stromal cells",
  "VE" = "Vascular endothelial cells",
  "Microglia" = "Macrophage",
  "ciliary_epithelial" = "Ciliary epithelial cells"
)
table(Mouse_RPE$celltype)




#####
#####----------------------------------------------------------#####
#####




library(ggplot2)

umap_df <- data.frame(
  UMAP1 = Embeddings(Mouse_RPE, "umap")[, 1],
  UMAP2 = Embeddings(Mouse_RPE, "umap")[, 2],
  celltype = Mouse_RPE$celltype,
  sample = factor(Mouse_RPE$sample, levels = c("5W_RPE", "49W_RPE", "91W_RPE"))
)

# Nature 风格配色
celltype_colors <- c(
  "Retinal pigment epithelial cells" = "#E64B35",
  "Stromal cells" = "#4DBBD5",
  "Vascular endothelial cells" = "#00A087",
  "Macrophage" = "#3C5488",
  "Ciliary epithelial cells" = "#F39B7F"
)

sample_colors <- c(
  "5W_RPE" = "#4DBBD5",
  "49W_RPE" = "#E64B35",
  "91W_RPE" = "#3C5488"
)

square_theme <- theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1,
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm"),
    plot.title = element_text(size = 13)
  )

p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = celltype_colors) +
  square_theme +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell Type", title = "UMAP by Cell Type")

p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sample)) +
  geom_point(size = 0.2, alpha = 0.6) +
  scale_color_manual(values = sample_colors) +
  square_theme +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Sample", title = "UMAP by Sample")

ggsave("Mouse_RPE_UMAP_celltype.png", plot = p1, width = 6, height = 5, dpi = 300)
ggsave("Mouse_RPE_UMAP_sample.png", plot = p2, width = 4, height = 3, dpi = 300)

##### Dot plot - cell type markers #####

marker_symbols <- c(
  "Pax6", "Otx1",            # Ciliary epithelial cells
  "C1qa", "Csf1r",           # Macrophage
  "Rpe65", "Lrat",           # RPE
  "Col1a1", "Dcn",           # Stromal cells
  "Pecam1", "Cdh5"           # Vascular endothelial cells
)

all_genes <- rownames(Mouse_RPE)
gene_symbols <- sub(".*~~", "", all_genes)
markers <- sapply(marker_symbols, function(s) {
  idx <- which(gene_symbols == s)
  if (length(idx) > 0) all_genes[idx[1]] else NA
})
markers <- markers[!is.na(markers)]

Mouse_RPE$celltype <- factor(Mouse_RPE$celltype, levels = rev(c(
  "Ciliary epithelial cells",
  "Macrophage",
  "Retinal pigment epithelial cells",
  "Stromal cells",
  "Vascular endothelial cells"
)))
Idents(Mouse_RPE) <- "celltype"

# 去名、去重，避免 DotPlot 因命名向量触发分面 / 行名重复报错
markers <- unique(unname(markers))

library(ggplot2)

p_dot <- DotPlot(Mouse_RPE, features = markers, assay = "RNA") +
  facet_null() +
  scale_x_discrete(labels = function(x) sub(".*~~", "", x)) +
  scale_color_gradientn(colors = c("grey90", "#FCBBA1", "#CB181D", "#67000D")) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  ) +
  labs(x = "", y = "")

ggsave("Mouse_RPE_dotplot_markers.png", plot = p_dot, width = 8, height = 4, dpi = 300)

#####



##### 下面找差异基因 ####
##### 画Heatmaps ######
#####

# ============================== METHODS ==============================
# Identification of RPE aging-associated genes with a consistent
# age-dependent direction
#
# Retinal pigment epithelial (RPE) cells were extracted from the annotated
# mouse dataset (celltype == "Retinal pigment epithelial cells") and split
# by age into three groups (5W, 49W and 91W). Expression values were
# re-normalized on the RNA assay using Seurat's global-scaling normalization
# (NormalizeData, "LogNormalize", scale factor 1e4). Differentially expressed
# genes (DEGs) were identified for the three pairwise age comparisons
# (49W vs 5W, 91W vs 49W and 91W vs 5W) with FindMarkers (Wilcoxon rank-sum
# test) using a minimum detection rate of 5% (min.pct = 0.05) and a
# log-fold-change threshold of 0.1. For each comparison, genes with a nominal
# p-value < 0.05 were retained as significant.
#
# To obtain robust, monotonically age-regulated genes, we first took the
# intersection of the significant genes across all three comparisons, and then
# kept only those whose log2 fold changes had the same sign in every
# comparison (i.e. consistently up- or down-regulated with age).
#
# These candidates were further validated at the pseudobulk level. Raw UMI
# counts were summed across all cells within each age group to form pseudobulk
# profiles, normalized to a common library size (counts per 10,000) and
# log2-transformed (log2(CPM + 1)). Only genes for which the direction of the
# pseudobulk difference agreed with the single-cell log2 fold change in all
# three comparisons were retained. Genes were classified as up- or
# down-regulated with aging according to the sign of the 91W-vs-5W log2 fold
# change.
#
# For visualization, the pseudobulk log2-expression matrix of the final gene
# set was Z-score scaled per gene across the three age groups and displayed as
# a heatmap (ComplexHeatmap), with genes split into up- and down-regulated
# blocks and age groups ordered 5W -> 49W -> 91W.
# =====================================================================


table(Mouse_RPE$celltype)

Mouse_RPE_cl = Mouse_RPE[,which(Mouse_RPE$celltype %in% c("Retinal pigment epithelial cells"))]
table(Mouse_RPE_cl$celltype)
table(Mouse_RPE_cl$sample)

##### RPE Aging Genes - Consistent Direction across 5W, 49W, 91W #####

library(ComplexHeatmap)
library(circlize)
library(dplyr)

Mouse_RPE_cl$sample <- factor(Mouse_RPE_cl$sample, levels = c("5W_RPE", "49W_RPE", "91W_RPE"))
Idents(Mouse_RPE_cl) <- "sample"

DefaultAssay(Mouse_RPE_cl) <- "RNA"
Mouse_RPE_cl <- NormalizeData(Mouse_RPE_cl)


# 5W vs 49W
markers_5v49 <- FindMarkers(Mouse_RPE_cl, ident.1 = "49W_RPE", ident.2 = "5W_RPE",
                             logfc.threshold = 0.1, min.pct = 0.05)
markers_5v49$gene <- rownames(markers_5v49)

# 49W vs 91W
markers_49v91 <- FindMarkers(Mouse_RPE_cl, ident.1 = "91W_RPE", ident.2 = "49W_RPE",
                              logfc.threshold = 0.1, min.pct = 0.05)
markers_49v91$gene <- rownames(markers_49v91)

# 5W vs 91W
markers_5v91 <- FindMarkers(Mouse_RPE_cl, ident.1 = "91W_RPE", ident.2 = "5W_RPE",
                             logfc.threshold = 0.1, min.pct = 0.05)
markers_5v91$gene <- rownames(markers_5v91)

sig_5v49 <- markers_5v49 %>% filter(p_val < 0.05)
sig_49v91 <- markers_49v91 %>% filter(p_val < 0.05)
sig_5v91 <- markers_5v91 %>% filter(p_val < 0.05)

common_genes <- Reduce(intersect, list(sig_5v49$gene, sig_49v91$gene, sig_5v91$gene))

df_merge <- data.frame(
  gene = common_genes,
  logFC_5v49 = sig_5v49[common_genes, "avg_log2FC"],
  logFC_49v91 = sig_49v91[common_genes, "avg_log2FC"],
  logFC_5v91 = sig_5v91[common_genes, "avg_log2FC"]
)


########### ############# #########
########### ############# #########

# 方向一致：5->49 和 49->91 和 5->91 同号
consistent_genes <- df_merge %>%
  filter(sign(logFC_5v49) == sign(logFC_49v91) & sign(logFC_5v49) == sign(logFC_5v91))

cat("Consistent logFC genes:", nrow(consistent_genes), "\n")

####### 428 genes #######

all_consistent <- consistent_genes$gene

# Pseudobulk: 全基因聚合 → 用全基因文库大小归一化 → log2 → 再提取差异基因
counts_mat <- GetAssayData(Mouse_RPE_cl, assay = "RNA", layer = "counts")

sample_ids <- Mouse_RPE_cl$sample
# 对所有基因按年龄组聚合成 pseudobulk
bulk_counts_all <- sapply(levels(sample_ids), function(s) {
  rowSums(counts_mat[, sample_ids == s, drop = FALSE])
})

# 文库大小按全部基因计算（每组所有基因的 counts 之和）
lib_size <- colSums(bulk_counts_all)
cpm_all <- t(t(bulk_counts_all) / lib_size * 1e4)
mat_all <- log2(cpm_all + 1)

# 归一化后再提取差异基因，并固定年龄顺序
mat <- mat_all[all_consistent, c("5W_RPE", "49W_RPE", "91W_RPE")]

# pseudobulk 差值方向也要和 logFC 一致
consistent_genes$pb_5v49 <- mat[consistent_genes$gene, "49W_RPE"] - mat[consistent_genes$gene, "5W_RPE"]
consistent_genes$pb_49v91 <- mat[consistent_genes$gene, "91W_RPE"] - mat[consistent_genes$gene, "49W_RPE"]
consistent_genes$pb_5v91 <- mat[consistent_genes$gene, "91W_RPE"] - mat[consistent_genes$gene, "5W_RPE"]

consistent_genes <- consistent_genes %>%
  filter(sign(pb_5v49) == sign(logFC_5v49) &
         sign(pb_49v91) == sign(logFC_49v91) &
         sign(pb_5v91) == sign(logFC_5v91))


cutoff = 0.75
consistent_genes = consistent_genes[which(abs(consistent_genes$logFC_5v91) > cutoff),]


up_genes <- consistent_genes %>% filter(logFC_5v91 > 0) %>% pull(gene)
down_genes <- consistent_genes %>% filter(logFC_5v91 < 0) %>% pull(gene)

cat("After cor filter - Up:", length(up_genes), "| Down:", length(down_genes), "\n")

# 更新 mat 只保留过滤后的基因
mat <- mat[c(up_genes, down_genes), ]

# Z-score per gene (row)
mat_z <- t(scale(t(mat)))

# 上调/下调分组注释
mat_z_ordered <- rbind(mat_z[up_genes, ], mat_z[down_genes, ])
gene_group_ordered <- c(rep("Up", length(up_genes)), rep("Down", length(down_genes)))

row_split <- factor(gene_group_ordered, levels = c("Up", "Down"))

col_fun <- colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))

png("RPE_aging_consistent_genes_heatmapXXX.png", width = 5, height = 6, units = "in", res = 300)
ht <- Heatmap(mat_z_ordered,
              name = "Z-score",
              col = col_fun,
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              row_split = row_split,
              row_gap = unit(3, "mm"),
              row_title = c("Up with aging", "Down with aging"),
              column_title = "RPE Aging Genes (Consistent Direction)",
              column_names_rot = 45,
              show_row_names = ifelse(nrow(mat_z_ordered) <= 80, TRUE, FALSE),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 12),
              heatmap_legend_param = list(title = "Z-score"))
draw(ht)
dev.off()

# 保存 heatmap 里的 DEGs
heatmap_degs <- data.frame(
  gene = c(up_genes, down_genes),
  direction = c(rep("Up", length(up_genes)), rep("Down", length(down_genes))),
  gene_symbol = sub(".*~~", "", c(up_genes, down_genes))
)


write.csv(heatmap_degs, "RPE_aging_consistent_DEGs.csv", row.names = FALSE)
cat("Saved", nrow(heatmap_degs), "DEGs to RPE_aging_consistent_DEGs.csv\n")




# ##### GO 分析 #####

# setwd("/projects/hmz-aging/RPE_smoking")
# heatmap_degs <- read.csv("RPE_aging_consistent_DEGs.csv")

# library(clusterProfiler)
# library(org.Mm.eg.db)

# # gene symbol → Entrez ID
# gene_symbols <- heatmap_degs$gene_symbol
# entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# up_symbols <- heatmap_degs$gene_symbol[heatmap_degs$direction == "Up"]
# down_symbols <- heatmap_degs$gene_symbol[heatmap_degs$direction == "Down"]

# up_entrez <- entrez$ENTREZID[entrez$SYMBOL %in% up_symbols]
# down_entrez <- entrez$ENTREZID[entrez$SYMBOL %in% down_symbols]

# # GO BP - Up genes
# go_up <- enrichGO(gene = up_entrez,
#                   OrgDb = org.Mm.eg.db,
#                   ont = "BP",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff = 0.05,
#                   readable = TRUE)

# # GO BP - Down genes
# go_down <- enrichGO(gene = down_entrez,
#                     OrgDb = org.Mm.eg.db,
#                     ont = "BP",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff = 0.05,
#                     readable = TRUE)

# library(ggplot2)

# # 过滤 GeneRatio > 0.03
# filter_go <- function(go_res, ratio_cutoff = 0.03) {
#   df <- as.data.frame(go_res)
#   ratio_val <- sapply(df$GeneRatio, function(x) {
#     parts <- as.numeric(strsplit(x, "/")[[1]])
#     parts[1] / parts[2]
#   })
#   df <- df[ratio_val > ratio_cutoff, ]
#   df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
#     parts <- as.numeric(strsplit(x, "/")[[1]])
#     parts[1] / parts[2]
#   })
#   df
# }

# go_up_df <- filter_go(go_up)
# go_down_df <- filter_go(go_down)

# nature_go_plot <- function(df, title, top_n = 15) {
#   df <- head(df[order(-df$GeneRatioNum), ], top_n)
#   df$Description <- factor(df$Description, levels = rev(df$Description))
#   ggplot(df, aes(x = GeneRatioNum, y = Description, size = Count, color = -log10(p.adjust))) +
#     geom_point() +
#     scale_color_gradientn(colors = c("#4DBBD5", "#00A087", "#E64B35", "#B2182B"),
#                           name = "-log10(adj.P)") +
#     scale_size_continuous(range = c(2, 6), name = "Count") +
#     theme_classic() +
#     theme(
#       panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#       axis.text.y = element_text(size = 12, color = "black"),
#       axis.text.x = element_text(size = 10, color = "black"),
#       axis.title = element_text(size = 11),
#       plot.title = element_text(size = 13, face = "bold"),
#       legend.text = element_text(size = 9),
#       legend.title = element_text(size = 10)
#     ) +
#     labs(x = "Gene Ratio", y = "", title = title)
# }

# p_up <- nature_go_plot(go_up_df, "GO BP - Up with Aging")
# p_down <- nature_go_plot(go_down_df, "GO BP - Down with Aging")

# ggsave("GO_BP_up_aging.png", plot = p_up, width = 7, height = 5, dpi = 300)
# ggsave("GO_BP_down_aging.png", plot = p_down, width = 7, height = 4, dpi = 300)

# cat("Up genes:", length(up_genes), "| Down genes:", length(down_genes), "\n")
# cat("Up GO terms:", nrow(go_up_df), "| Down GO terms:", nrow(go_down_df), "\n")



#####
#####

library(Seurat)
setwd("/projects/hmz-aging/RPE_smoking")
load(file.path("Handa_snRNA_Final 1.RData"))

smoking_RPE = combined 

table(smoking_RPE$group)
table(smoking_RPE$group_1)

#### 1 RPE
#### 9 Sick RPE

cluster_num <- as.numeric(sub(".*_", "", smoking_RPE$group_1))
smoking_RPE$celltype <- NA
smoking_RPE$celltype[cluster_num == 1] <- "RPE"
smoking_RPE$celltype[cluster_num == 9] <- "SickRPE"
table(smoking_RPE$celltype)

smoking_RPE$celltype <- paste(smoking_RPE$group, smoking_RPE$celltype, sep = "_")
table(smoking_RPE$celltype)

smoking_RPE_cl <- smoking_RPE[, !grepl("NA", smoking_RPE$celltype)]
table(smoking_RPE_cl$celltype)

DefaultAssay(smoking_RPE_cl) <- "RNA"
smoking_RPE_cl <- NormalizeData(smoking_RPE_cl)

Idents(smoking_RPE_cl) <- "celltype"
avg <- AverageExpression(smoking_RPE_cl, features = "Apoe", assays = "RNA", layer = "data")
print(avg$RNA)


#####
#####
#####




##### 用 heatmap DEGs 计算 aging score #####
heatmap_degs <- read.csv("RPE_aging_consistent_DEGs.csv")

up_hm <- heatmap_degs$gene_symbol[heatmap_degs$direction == "Up"]
down_hm <- heatmap_degs$gene_symbol[heatmap_degs$direction == "Down"]

smoking_RPE_cl <- AddModuleScore(smoking_RPE_cl,
                                  features = list(Up = up_hm, Down = down_hm),
                                  name = "aging_score")
smoking_RPE_cl$aging_score <- smoking_RPE_cl$aging_score1 - smoking_RPE_cl$aging_score2

tapply(smoking_RPE_cl$aging_score, smoking_RPE_cl$celltype, mean)


########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
################### UMAP ###############################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################

conda activate seurat_env
R


library(Seurat)
setwd("/projects/hmz-aging/RPE_smoking")
load(file.path("Handa_snRNA_Final 1.RData"))

smoking_RPE = combined 

table(smoking_RPE$group)
table(smoking_RPE$group_1)

#### 1 RPE ##########
#### 9 Sick RPE #####

cluster_num <- as.numeric(sub(".*_", "", smoking_RPE$group_1))
smoking_RPE$celltype <- NA
smoking_RPE$celltype[cluster_num == 1] <- "RPE"
smoking_RPE$celltype[cluster_num == 9] <- "SickRPE"
table(smoking_RPE$celltype)

smoking_RPE$celltype <- paste(smoking_RPE$group, smoking_RPE$celltype, sep = "_")
table(smoking_RPE$celltype)

smoking_RPE_cl <- smoking_RPE[, !grepl("NA", smoking_RPE$celltype)]
table(smoking_RPE_cl$celltype)

DefaultAssay(smoking_RPE_cl) <- "RNA"
smoking_RPE_cl <- NormalizeData(smoking_RPE_cl)






# 排列顺序：con/cse pair, 3d→6d→10d→10d_old, RPE先 SickRPE后
level_order <- c(
  "con_3d_RPE",        "cse_3d_RPE",
  "con_3d_SickRPE",    "cse_3d_SickRPE",
  "con_6d_RPE",        "cse_6d_RPE",
  "con_6d_SickRPE",    "cse_6d_SickRPE",
  "con_10d_RPE",       "cse_10d_RPE",
  "con_10d_SickRPE",   "cse_10d_SickRPE",
  "con_10d_old_RPE",   "cse_10d_old_RPE",
  "con_10d_old_SickRPE", "cse_10d_old_SickRPE"
)
# 只保留实际存在的 level
level_order <- level_order[level_order %in% unique(smoking_RPE_cl$celltype)]

smoking_RPE_cl$celltype <- factor(smoking_RPE_cl$celltype, levels = level_order)





# 只保留指定样品
keep_ct <- c("con_3d_RPE", "cse_3d_RPE",
             "con_6d_RPE", "cse_6d_RPE", "cse_6d_SickRPE",
             "con_10d_old_RPE", "cse_10d_old_RPE", "con_10d_old_SickRPE")
smoking_RPE_clcl <- smoking_RPE_cl[, smoking_RPE_cl$celltype %in% keep_ct]
smoking_RPE_clcl$celltype <- droplevels(factor(smoking_RPE_clcl$celltype, levels = keep_ct))
table(smoking_RPE_clcl$celltype)

# 重新跑 UMAP
DefaultAssay(smoking_RPE_clcl) <- "RNA"
smoking_RPE_clcl <- NormalizeData(smoking_RPE_clcl)
smoking_RPE_clcl <- FindVariableFeatures(smoking_RPE_clcl, nfeatures = 2000)
smoking_RPE_clcl <- ScaleData(smoking_RPE_clcl)
smoking_RPE_clcl <- RunPCA(smoking_RPE_clcl, npcs = 30)
smoking_RPE_clcl <- RunUMAP(smoking_RPE_clcl, dims = 1:10)

# Plot UMAP by celltype
#
library(ggplot2)

umap_df2 <- data.frame(
  UMAP1 = Embeddings(smoking_RPE_clcl, "umap")[, 1],
  UMAP2 = Embeddings(smoking_RPE_clcl, "umap")[, 2],
  celltype = smoking_RPE_clcl$celltype
)

p_umap <- ggplot(umap_df2, aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = c(
    "con_3d_RPE" = "#4DBBD5", "cse_3d_RPE" = "#E64B35",
    "con_6d_RPE" = "#00A087", "cse_6d_RPE" = "#3C5488", "cse_6d_SickRPE" = "#F39B7F",
    "con_10d_old_RPE" = "#7E6148", "cse_10d_old_RPE" = "#B09C85", "con_10d_old_SickRPE" = "#DC0000"
  )) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1,
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL)

ggsave("smoking_RPE_clcl_UMAPXXX.png", plot = p_umap, width = 6, height = 4, dpi = 300)

########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
################### UMAP ###############################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################



##### 分组 violin plot + 显著性检验 #####

library(ggplot2)
library(dplyr)

group_colors <- c("con" = "#4DBBD5", "cse" = "#E64B35", "SickRPE" = "#7E6148")

smoking_RPE_clcl <- AddModuleScore(smoking_RPE_clcl,
                                    features = list(Up = up_hm, Down = down_hm),
                                    name = "aging_sc")
smoking_RPE_clcl$aging_score <- smoking_RPE_clcl$aging_sc1 - smoking_RPE_clcl$aging_sc2

df_all <- data.frame(
  celltype = as.character(smoking_RPE_clcl$celltype),
  aging_score = smoking_RPE_clcl$aging_score
)

group1 <- c("con_3d_RPE", "cse_3d_RPE")
group2 <- c("con_6d_RPE", "cse_6d_RPE", "cse_6d_SickRPE")
group3 <- c("con_10d_old_RPE", "cse_10d_old_RPE", "con_10d_old_SickRPE")

assign_color_group <- function(ct) {
  ifelse(grepl("SickRPE", ct), "SickRPE",
         ifelse(grepl("^con_", ct), "con", "cse"))
}

plot_violin_group <- function(df_all, celltypes, title, filename) {
  df <- df_all %>% filter(celltype %in% celltypes)
  df$celltype <- factor(df$celltype, levels = celltypes)
  df$color_group <- assign_color_group(df$celltype)

  con_ct <- celltypes[grepl("^con_.*RPE$", celltypes) & !grepl("Sick", celltypes)][1]
  other_cts <- setdiff(celltypes, con_ct)

  # Wilcoxon test vs control
  pval_df <- do.call(rbind, lapply(other_cts, function(x) {
    con_vals <- df$aging_score[df$celltype == con_ct]
    test_vals <- df$aging_score[df$celltype == x]
    pv <- wilcox.test(test_vals, con_vals, alternative = "greater")$p.value
    label <- ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", "ns")))
    data.frame(group1 = con_ct, group2 = x, p = pv, label = label)
  }))

  y_max <- max(df$aging_score, na.rm = TRUE)
  y_step <- (y_max - min(df$aging_score, na.rm = TRUE)) * 0.1
  pval_df$x1 <- match(pval_df$group1, celltypes)
  pval_df$x2 <- match(pval_df$group2, celltypes)
  # 从 violin 顶端再多留出一段间距，避免横线和图形重叠
  pval_df$y <- y_max + y_step * (seq_len(nrow(pval_df)) + 1.5)

  p <- ggplot(df, aes(x = celltype, y = aging_score, fill = color_group)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
    scale_fill_manual(values = group_colors) +
    geom_segment(data = pval_df, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = pval_df, aes(x = (x1 + x2) / 2, y = y + y_step * 0.3, label = label),
              inherit.aes = FALSE, size = 4) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    ) +
    labs(x = "", y = "Aging Score", title = title)

  ggsave(filename, plot = p, width = 3, height = 4.5, dpi = 300)
  return(p)
}

p_3d <- plot_violin_group(df_all, group1, "3d", "aging_violin_3dXXX.png")
p_6d <- plot_violin_group(df_all, group2, "6d", "aging_violin_6dXXX.png")
p_10d_old <- plot_violin_group(df_all, group3, "10d_old", "aging_violin_10d_oldXXX.png")









#####################################################
############ Next is the aging DEGs Heatmaps #########
#####################################################




########
##### Smoking DEGs - 变化趋势一致 across 3d, 6d, 10d_old #####
########

Idents(smoking_RPE_clcl) <- "celltype"

# cse vs con: 3d
markers_3d <- FindMarkers(smoking_RPE_clcl, ident.1 = "cse_3d_RPE", ident.2 = "con_3d_RPE",
                           logfc.threshold = 0.1, min.pct = 0.05)
markers_3d$gene <- rownames(markers_3d)

# cse vs con: 6d
markers_6d <- FindMarkers(smoking_RPE_clcl, ident.1 = "cse_6d_RPE", ident.2 = "con_6d_RPE",
                           logfc.threshold = 0.1, min.pct = 0.05)
markers_6d$gene <- rownames(markers_6d)

# cse vs con: 10d_old
markers_10d_old <- FindMarkers(smoking_RPE_clcl, ident.1 = "cse_10d_old_RPE", ident.2 = "con_10d_old_RPE",
                                logfc.threshold = 0.1, min.pct = 0.05)
markers_10d_old$gene <- rownames(markers_10d_old)

sig_3d <- markers_3d %>% filter(p_val < 0.05)
sig_6d <- markers_6d %>% filter(p_val < 0.05)
sig_10d_old <- markers_10d_old %>% filter(p_val < 0.05)

common_smoke <- Reduce(intersect, list(sig_3d$gene, sig_6d$gene, sig_10d_old$gene))
cat("Common significant genes:", length(common_smoke), "\n")

df_smoke <- data.frame(
  gene = common_smoke,
  logFC_3d = sig_3d[common_smoke, "avg_log2FC"],
  logFC_6d = sig_6d[common_smoke, "avg_log2FC"],
  logFC_10d_old = sig_10d_old[common_smoke, "avg_log2FC"]
)

# 方向一致：3d, 6d, 10d_old 同号
consistent_smoke <- df_smoke %>%
  filter(sign(logFC_3d) == sign(logFC_6d) & sign(logFC_3d) == sign(logFC_10d_old))

cat("Consistent logFC genes:", nrow(consistent_smoke), "\n")

# Pseudobulk 验证方向
counts_clcl <- GetAssayData(smoking_RPE_clcl, assay = "RNA", layer = "counts")
pb_celltypes <- c("con_3d_RPE", "cse_3d_RPE", "con_6d_RPE", "cse_6d_RPE", "con_10d_old_RPE", "cse_10d_old_RPE")
ct_ids <- smoking_RPE_clcl$celltype

bulk_smoke <- sapply(pb_celltypes, function(ct) {
  idx <- which(ct_ids == ct)
  if (length(idx) == 1) return(counts_clcl[, idx])
  Matrix::rowSums(counts_clcl[, idx, drop = FALSE])
})

lib_smoke <- colSums(bulk_smoke)
cpm_smoke <- t(t(bulk_smoke) / lib_smoke * 1e4)
mat_smoke <- log2(cpm_smoke + 1)

# pseudobulk 差值方向也要和 logFC 一致
consistent_smoke$pb_3d <- mat_smoke[consistent_smoke$gene, "cse_3d_RPE"] - mat_smoke[consistent_smoke$gene, "con_3d_RPE"]
consistent_smoke$pb_6d <- mat_smoke[consistent_smoke$gene, "cse_6d_RPE"] - mat_smoke[consistent_smoke$gene, "con_6d_RPE"]
consistent_smoke$pb_10d_old <- mat_smoke[consistent_smoke$gene, "cse_10d_old_RPE"] - mat_smoke[consistent_smoke$gene, "con_10d_old_RPE"]

consistent_smoke <- consistent_smoke %>%
  filter(sign(pb_3d) == sign(logFC_3d) &
         sign(pb_6d) == sign(logFC_6d) &
         sign(pb_10d_old) == sign(logFC_10d_old))

up_smoke <- consistent_smoke %>% filter(logFC_3d > 0) %>% pull(gene)
down_smoke <- consistent_smoke %>% filter(logFC_3d < 0) %>% pull(gene)

cat("After pb filter - Up:", length(up_smoke), "| Down:", length(down_smoke), "\n")

# 保存
smoke_degs <- data.frame(
  gene = c(up_smoke, down_smoke),
  direction = c(rep("Up", length(up_smoke)), rep("Down", length(down_smoke)))
)
write.csv(smoke_degs, "Smoking_consistent_DEGs.csv", row.names = FALSE)



##### Smoking ∩ Aging: 方向一致的基因 #####

heatmap_degs <- read.csv("RPE_aging_consistent_DEGs.csv")
smoke_degs <- read.csv("Smoking_consistent_DEGs.csv")

# 检查列名
cat("Aging DEGs columns:", colnames(heatmap_degs), "\n")
cat("Smoking DEGs columns:", colnames(smoke_degs), "\n")

# 用 gene_symbol 匹配 smoking 的 gene
overlap <- merge(heatmap_degs, smoke_degs, by.x = "gene_symbol", by.y = "gene", suffixes = c("_aging", "_smoking"))

# 方向一致：aging 上调 + smoking 上调，或 aging 下调 + smoking 下调
same_direction <- overlap %>% filter(direction_aging == direction_smoking)

# 方向相反的也看看
opposite_direction <- overlap %>% filter(direction_aging != direction_smoking)

cat("Overlap total:", nrow(overlap), "\n")
cat("Same direction (smoking accelerates aging):", nrow(same_direction), "\n")
cat("  Both Up:", sum(same_direction$direction_aging == "Up"), "\n")
cat("  Both Down:", sum(same_direction$direction_aging == "Down"), "\n")
cat("Opposite direction:", nrow(opposite_direction), "\n")

print(same_direction)

write.csv(same_direction, "Smoking_Aging_same_direction_genes.csv", row.names = FALSE)
write.csv(opposite_direction, "Smoking_Aging_opposite_direction_genes.csv", row.names = FALSE)

# Heatmap: smoking (左) + aging (右)，方向一致的基因

# 重新读入数据并计算 pseudobulk
library(Seurat)
library(dplyr)
setwd("/projects/hmz-aging/RPE_smoking")

load(file.path("Handa_snRNA_Final 1.RData"))
Mouse_RPE <- readRDS(file.path("Mouse_RPE_RNA_Final_Jul16"))

smoking_RPE <- combined
cluster_num <- as.numeric(sub(".*_", "", smoking_RPE$group_1))
smoking_RPE$celltype <- NA
smoking_RPE$celltype[cluster_num == 1] <- "RPE"
smoking_RPE$celltype[cluster_num == 9] <- "SickRPE"
smoking_RPE$celltype <- paste(smoking_RPE$group, smoking_RPE$celltype, sep = "_")
smoking_RPE_cl <- smoking_RPE[, !grepl("NA", smoking_RPE$celltype)]
DefaultAssay(smoking_RPE_cl) <- "RNA"
smoking_RPE_cl <- NormalizeData(smoking_RPE_cl)

keep_ct <- c("con_3d_RPE", "cse_3d_RPE",
             "con_6d_RPE", "cse_6d_RPE", "cse_6d_SickRPE",
             "con_10d_old_RPE", "cse_10d_old_RPE", "con_10d_old_SickRPE")
smoking_RPE_clcl <- smoking_RPE_cl[, smoking_RPE_cl$celltype %in% keep_ct]
smoking_RPE_clcl$celltype <- factor(smoking_RPE_clcl$celltype, levels = keep_ct)

# Mouse RPE pseudobulk (aging)
Mouse_RPE_cl <- Mouse_RPE[, Mouse_RPE$celltype %in% c("RPE")]
Mouse_RPE_cl$sample <- factor(Mouse_RPE_cl$sample, levels = c("5W_RPE", "49W_RPE", "91W_RPE"))

counts_aging <- GetAssayData(Mouse_RPE_cl, assay = "RNA", layer = "counts")
sample_ids_aging <- Mouse_RPE_cl$sample
bulk_aging <- sapply(levels(sample_ids_aging), function(s) {
  Matrix::rowSums(counts_aging[, sample_ids_aging == s, drop = FALSE])
})
lib_aging <- colSums(bulk_aging)
cpm_aging <- t(t(bulk_aging) / lib_aging * 1e4)
mat <- log2(cpm_aging + 1)

# Smoking pseudobulk
smoke_cols <- c("con_3d_RPE", "con_6d_RPE", "con_10d_old_RPE",
                "cse_3d_RPE", "cse_6d_RPE", "cse_10d_old_RPE")
counts_smoke <- GetAssayData(smoking_RPE_clcl, assay = "RNA", layer = "counts")
ct_ids <- smoking_RPE_clcl$celltype
bulk_smoke_all <- sapply(smoke_cols, function(ct) {
  idx <- which(ct_ids == ct)
  Matrix::rowSums(counts_smoke[, idx, drop = FALSE])
})
lib_smoke <- colSums(bulk_smoke_all)
cpm_smoke <- t(t(bulk_smoke_all) / lib_smoke * 1e4)
mat_smoke <- log2(cpm_smoke + 1)

# 读入 DEGs
heatmap_degs <- read.csv("RPE_aging_consistent_DEGs.csv")
smoke_degs <- read.csv("Smoking_consistent_DEGs.csv")

overlap <- merge(heatmap_degs, smoke_degs, by.x = "gene_symbol", by.y = "gene", suffixes = c("_aging", "_smoking"))
same_direction <- overlap %>% filter(direction_aging == direction_smoking)
opposite_direction <- overlap %>% filter(direction_aging != direction_smoking)

cat("Same direction:", nrow(same_direction), "\n")
cat("Opposite direction:", nrow(opposite_direction), "\n")



if (nrow(same_direction) > 0) {
  library(ComplexHeatmap)
  library(circlize)

  # 用 gene_symbol 匹配
  same_symbols <- same_direction$gene_symbol

  # Aging pseudobulk: 需要从 ENSMUSG~~Symbol 的 rownames 中找
  aging_gene_map <- setNames(sub(".*~~", "", rownames(mat)), rownames(mat))
  aging_rows <- names(aging_gene_map)[aging_gene_map %in% same_symbols]
  aging_mat_sub <- mat[aging_rows, c("5W_RPE", "49W_RPE", "91W_RPE"), drop = FALSE]
  rownames(aging_mat_sub) <- aging_gene_map[aging_rows]
  aging_z <- t(scale(t(aging_mat_sub)))

  # Smoking pseudobulk: 分三块，每块 con + cse，各自 z-score
  smoke_genes <- same_symbols[same_symbols %in% rownames(mat_smoke)]
  common_sym <- intersect(rownames(aging_z), smoke_genes)

  # 上调/下调分组
  dir_map <- setNames(same_direction$direction_aging, same_direction$gene_symbol)
  up_sym_ht <- common_sym[dir_map[common_sym] == "Up"]
  down_sym_ht <- common_sym[dir_map[common_sym] == "Down"]
  ordered_sym <- c(up_sym_ht, down_sym_ht)

  row_split_ht <- factor(c(rep("Up", length(up_sym_ht)), rep("Down", length(down_sym_ht))),
                          levels = c("Up", "Down"))

  # 三个时间点合并，列顺序：前三 con_（3d/6d/10d），后三 cse_（3d/6d/10d），一起算 z-score
  smoke_col_order <- c("con_3d_RPE", "con_6d_RPE", "con_10d_old_RPE",
                       "cse_3d_RPE", "cse_6d_RPE", "cse_10d_old_RPE")
  mat_smoke_sub <- mat_smoke[ordered_sym, smoke_col_order]
  z_smoke <- t(scale(t(mat_smoke_sub)))

  # 颜色
  col_smoke <- colorRamp2(c(-2, 0, 2), c("#00A087", "white", "#E64B35"))
  col_aging <- colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#DC0000"))

  # smoking 列分组注释：前三 Control，后三 Smoking
  smoke_group <- c(rep("Control", 3), rep("Smoking", 3))
  col_anno_smoke <- HeatmapAnnotation(
    Group = smoke_group,
    col = list(Group = c("Control" = "#4DBBD5", "Smoking" = "#E64B35")),
    annotation_name_gp = gpar(fontsize = 10)
  )

  # 数值标注函数（aging heatmap 仍用）
  cell_label <- function(m) {
    function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.1f", m[i, j]), x, y, gp = gpar(fontsize = 7, col = "black"))
    }
  }

  ht_smoke <- Heatmap(z_smoke, name = "Smoking\nZ-score", col = col_smoke,
                      rect_gp = gpar(col = "white", lwd = 0.5),
                      top_annotation = col_anno_smoke,
                      cluster_columns = FALSE, cluster_rows = FALSE,
                      row_split = row_split_ht, row_gap = unit(3, "mm"),
                      row_title = c("Up", "Down"),
                      row_title_gp = gpar(fontsize = 12),
                      column_title = "Smoking", column_title_gp = gpar(fontsize = 12),
                      show_row_names = TRUE, row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                      row_names_side = "left",
                      column_names_rot = 45, column_names_gp = gpar(fontsize = 11))

  ht_aging <- Heatmap(aging_z[ordered_sym, ], name = "Aging\nZ-score", col = col_aging,
                      rect_gp = gpar(col = "white", lwd = 0.5),
                      cluster_columns = FALSE, cluster_rows = FALSE,
                      row_split = row_split_ht, row_gap = unit(3, "mm"),
                      column_title = "Aging (5W→49W→91W)", column_title_gp = gpar(fontsize = 12),
                      show_row_names = FALSE,
                      column_names_rot = 45, column_names_gp = gpar(fontsize = 11))

  png("Smoking_Aging_overlap_heatmapXXX.png", width = 12,
      height = max(6, length(ordered_sym) * 0.15 + 3), units = "in", res = 300)
  draw(ht_smoke + ht_aging,
       column_title = "Smoking & Aging: Consistent DEGs",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       ht_gap = unit(8, "mm"))
  dev.off()
}





