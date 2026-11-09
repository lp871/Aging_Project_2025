#######
####### 我们先找一下 Human RPE 的 scRNAseq ###### 
####### metformin human bulk iRPE #############
####### 然后构建 RPE Aging clocks ##############
#######


# conda activate seurat4
# R


#######
####### 分别 把 男女的 sample 分开 ##########
#######


#######
#######

# H_ambient_genes <- readRDS("/zp1/data/plyu3/Aging_Clocks_Final/Human_ambient_RNA_Jul29")
# H_ambient_genes <- H_ambient_genes[order(H_ambient_genes$Avg,decreasing=T),]
# H_ambient_genes_Top50 = H_ambient_genes$Gene[which(H_ambient_genes$Avg > 0.005)]
# H_ambient_genes_Top50 = sapply(strsplit(H_ambient_genes_Top50,split="~~"),function(x) x[[2]])

#######
####### 在 python 里面读入 然后变成 seurat #######
#######


##### /projects/hmz-aging/Metformin/
##### Human_RPE_atlas.h5ad
##### python 里面读入 #######

conda activate scvi-env
python

import os
os.chdir("/projects/hmz-aging/Metformin/")

import scanpy as sc
import anndata as ad

adata = sc.read_h5ad("/projects/hmz-aging/Metformin/Human_RPE_atlas.h5ad")

print(adata)
print(adata.obs.columns.tolist())
print(adata.obs.head())

# 查看细胞类型和样本信息
if 'cell_type' in adata.obs.columns:
    print(adata.obs['cell_type'].value_counts())

if 'subclass' in adata.obs.columns:
    print(adata.obs['subclass'].value_counts())

# 检查 X 是否为 raw counts
import numpy as np
x_sample = adata.X[:10, :10]
if hasattr(x_sample, 'toarray'):
    x_sample = x_sample.toarray()
print("X sample values:\n", x_sample)
print("X dtype:", adata.X.dtype)
print("X min:", adata.X.min(), "| max:", adata.X.max())
x_dense = x_sample if not hasattr(x_sample, 'toarray') else x_sample.toarray()
print("Has integers?", np.allclose(x_dense, np.round(x_dense)))

# 检查是否有 raw layer
if adata.raw is not None:
    print("\nadata.raw exists:", adata.raw.X.shape)

if 'counts' in adata.layers:
    print("adata.layers['counts'] exists")

print("Available layers:", list(adata.layers.keys()))

# 导出 mtx + tsv 给 Seurat 读入
import scipy.io
import pandas as pd

out_dir = "/projects/hmz-aging/Metformin/Human_RPE_atlas_mtx"
os.makedirs(out_dir, exist_ok=True)

# 导出 raw count matrix (genes x cells)
from scipy.sparse import csc_matrix
mat = csc_matrix(adata.X.T)
scipy.io.mmwrite(os.path.join(out_dir, "matrix.mtx"), mat)

# 导出 gene names
genes = pd.DataFrame({"gene": adata.var_names})
genes.to_csv(os.path.join(out_dir, "features.tsv"), sep="\t", header=False, index=False)

# 导出 cell barcodes
barcodes = pd.DataFrame({"barcode": adata.obs_names})
barcodes.to_csv(os.path.join(out_dir, "barcodes.tsv"), sep="\t", header=False, index=False)

# 导出 metadata
adata.obs.to_csv(os.path.join(out_dir, "metadata.tsv"), sep="\t")

print("Exported to:", out_dir)
print("Files:", os.listdir(out_dir))


##### 回到 R 读入 #####
# conda activate seurat_env
# R

library(Seurat)
library(Matrix)

mtx_dir <- "/projects/hmz-aging/Metformin/Human_RPE_atlas_mtx"

mat <- readMM(file.path(mtx_dir, "matrix.mtx"))
genes <- read.table(file.path(mtx_dir, "features.tsv"), stringsAsFactors = FALSE)$V1
barcodes <- read.table(file.path(mtx_dir, "barcodes.tsv"), stringsAsFactors = FALSE)$V1
rownames(mat) <- genes
colnames(mat) <- barcodes

meta <- read.table(file.path(mtx_dir, "metadata.tsv"), sep = "\t", header = TRUE, row.names = 1)

human_rpe <- CreateSeuratObject(counts = mat, meta.data = meta)
human_rpe

saveRDS(human_rpe, file.path(mtx_dir, "Human_RPE_atlas_seurat.rds"))
cat("Saved to:", file.path(mtx_dir, "Human_RPE_atlas_seurat.rds"), "\n")


#########
######### 新终端重新读入 #######
######### we remove RPE donors which cells < 100 ######
#########


conda activate seurat_env
R

library(Seurat)
library(ggplot2)
library(dplyr)

setwd("/projects/hmz-aging/Metformin/Human_RPE_atlas_mtx")

human_rpe <- readRDS("Human_RPE_atlas_seurat.rds")
human_rpe

table(human_rpe$Donor)

# 去掉小于 100 个 cell 的 donor
donor_counts <- table(human_rpe$Donor)
keep_donors <- names(donor_counts)[donor_counts >= 100]
cat("Donors before:", length(donor_counts), "| After:", length(keep_donors), "\n")
cat("Removed:", setdiff(names(donor_counts), keep_donors), "\n")

human_rpe_cl <- human_rpe[, human_rpe$Donor %in% keep_donors]
human_rpe_cl$Donor <- droplevels(factor(human_rpe_cl$Donor))
table(human_rpe_cl$Donor)

#####
#####
#####
table(human_rpe_cl$Age)

##### Pseudobulk by Donor #####

# Pseudobulk: 按 Donor 汇总 raw counts
Idents(human_rpe_cl) <- "Donor"
donors <- levels(Idents(human_rpe_cl))

counts_h <- GetAssayData(human_rpe_cl, assay = "RNA", layer = "counts")
pb <- sapply(donors, function(d) {
  idx <- which(Idents(human_rpe_cl) == d)
  if (length(idx) == 1) return(counts_h[, idx])
  Matrix::rowSums(counts_h[, idx, drop = FALSE])
})
colnames(pb) <- donors

# Donor metadata table
meta_donor <- human_rpe_cl@meta.data %>%
  group_by(Donor) %>%
  summarise(
    n_cells = n(),
    Age = first(Age),
    Gender = first(Gender),
    Ethnicity = first(Ethnicity)
  ) %>%
  as.data.frame()
rownames(meta_donor) <- meta_donor$Donor

print(meta_donor[order(meta_donor$Age), ])

# 创建 pseudobulk Seurat 对象
pb_seurat <- CreateSeuratObject(counts = pb, meta.data = meta_donor[colnames(pb), ])
pb_seurat


#####
##### Human RPE DEGs: Young ≤ 30 vs Old ≥ 60 #####
#####


library(limma)
library(edgeR)



##### 连续回归: Age + Gender #####

pb_counts_raw = pb_seurat[['RNA']]$counts

# 用全部 donor
dge_all <- DGEList(counts = pb_counts_raw)
dge_all <- calcNormFactors(dge_all, method = "TMM")

keep_all <- filterByExpr(dge_all)
dge_all <- dge_all[keep_all, , keep.lib.sizes = FALSE]
cat("Genes after filtering:", nrow(dge_all), "\n")

age_cont <- as.numeric(pb_seurat$Age)
gender <- factor(pb_seurat$Gender)

design_cont <- model.matrix(~ age_cont + gender)
cat("Design matrix:\n")
print(head(design_cont))

v_cont <- voom(dge_all, design_cont, plot = FALSE)
fit_cont <- lmFit(v_cont, design_cont)
fit_cont <- eBayes(fit_cont)

# coef=2 是 age 的效应（控制了 gender）
degs_cont <- topTable(fit_cont, coef = 2, number = Inf, sort.by = "P")
degs_cont$gene <- rownames(degs_cont)

cat("\n=== Continuous Age Regression (controlling Gender) ===\n")
cat("Total tested:", nrow(degs_cont), "\n")
cat("p < 0.05:", sum(degs_cont$P.Value < 0.05), "\n")
cat("adj.P < 0.05:", sum(degs_cont$adj.P.Val < 0.05), "\n")
cat("adj.P < 0.1:", sum(degs_cont$adj.P.Val < 0.1), "\n")
cat("  Up with age (adj.P<0.05):", sum(degs_cont$adj.P.Val < 0.05 & degs_cont$logFC > 0), "\n")
cat("  Down with age (adj.P<0.05):", sum(degs_cont$adj.P.Val < 0.05 & degs_cont$logFC < 0), "\n")

head(degs_cont[order(degs_cont$P.Value), ], 20)

# 用 P < 0.01 筛选
degs_sig <- degs_cont[degs_cont$P.Value < 0.005, ]
cat("\n=== P < 0.005 ===\n")
cat("Total sig:", nrow(degs_sig), "\n")
cat("Up with age:", sum(degs_sig$logFC > 0), "\n")
cat("Down with age:", sum(degs_sig$logFC < 0), "\n")

human_aging_up <- degs_sig$gene[degs_sig$logFC > 0]
human_aging_down <- degs_sig$gene[degs_sig$logFC < 0]

write.csv(degs_sig, "Human_RPE_aging_DEGs_continuous.csv", row.names = FALSE)

# 打印 top 30 up with age (按 logFC)
up_sorted <- degs_sig[degs_sig$logFC > 0, ]
up_sorted <- up_sorted[order(-up_sorted$logFC), ]
print(head(up_sorted[, c("gene", "logFC", "P.Value", "adj.P.Val")], 30))

##### Heatmap: Human RPE aging DEGs #####

# Normalize: CPM (1e4) + log2(x+1)
pb_counts <- as.matrix(pb_counts_raw)
lib_size <- colSums(pb_counts)
pb_norm <- sweep(pb_counts, 2, lib_size, FUN = "/") * 1e4
pb_norm <- log2(pb_norm + 1)

library(ComplexHeatmap)
library(circlize)

male_idx <- which(pb_seurat$Gender == "Male")
female_idx <- which(pb_seurat$Gender == "Female")
male_order <- male_idx[order(as.numeric(pb_seurat$Age[male_idx]))]
female_order <- female_idx[order(as.numeric(pb_seurat$Age[female_idx]))]

human_aging_up <- human_aging_up[human_aging_up %in% rownames(pb_norm)]
human_aging_down <- human_aging_down[human_aging_down %in% rownames(pb_norm)]
cat("Up genes in pb_norm:", length(human_aging_up), "| Down:", length(human_aging_down), "\n")

male_mat <- pb_norm[c(human_aging_up, human_aging_down), male_order]
male_z <- t(scale(t(male_mat)))

female_mat <- pb_norm[c(human_aging_up, human_aging_down), female_order]
female_z <- t(scale(t(female_mat)))

row_split_h <- factor(c(rep("Up with aging", length(human_aging_up)),
                         rep("Down with aging", length(human_aging_down))),
                       levels = c("Up with aging", "Down with aging"))

anno_male <- HeatmapAnnotation(
  Age = as.numeric(pb_seurat$Age[male_order]),
  col = list(Age = colorRamp2(c(0, 50, 100), c("#4DBBD5", "white", "#E64B35"))),
  annotation_name_gp = gpar(fontsize = 11)
)
anno_female <- HeatmapAnnotation(
  Age = as.numeric(pb_seurat$Age[female_order]),
  col = list(Age = colorRamp2(c(0, 50, 100), c("#4DBBD5", "white", "#E64B35"))),
  annotation_name_gp = gpar(fontsize = 11)
)

col_male <- colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
col_female <- colorRamp2(c(-2, 0, 2), c("#00A087", "white", "#DC0000"))

show_rn <- nrow(male_z) <= 80

ht_male <- Heatmap(male_z, name = "Male\nZ-score", col = col_male,
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   row_split = row_split_h, row_gap = unit(3, "mm"),
                   row_title = c("Up", "Down"), row_title_gp = gpar(fontsize = 11),
                   column_title = "Male", column_title_gp = gpar(fontsize = 13, fontface = "bold"),
                   top_annotation = anno_male,
                   show_row_names = show_rn, row_names_gp = gpar(fontsize = 5),
                   row_names_side = "left",
                   show_column_names = FALSE)

ht_female <- Heatmap(female_z, name = "Female\nZ-score", col = col_female,
                     cluster_columns = FALSE, cluster_rows = FALSE,
                     row_split = row_split_h, row_gap = unit(3, "mm"),
                     column_title = "Female", column_title_gp = gpar(fontsize = 13, fontface = "bold"),
                     top_annotation = anno_female,
                     show_row_names = FALSE,
                     show_column_names = FALSE)

ht_height <- min(30, max(6, nrow(male_z) * 0.08 + 3))


setwd("/projects/hmz-aging/Metformin/Human_RPE_atlas_mtx")

png("Human_RPE_aging_DEGs_heatmapXXX.png", width = 6, height = 6, units = "in", res = 300)
draw(ht_male + ht_female, ht_gap = unit(8, "mm"))
dev.off()

# 保存 heatmap 中的基因
heatmap_genes_h <- data.frame(
  gene = c(human_aging_up, human_aging_down),
  direction = c(rep("Up", length(human_aging_up)), rep("Down", length(human_aging_down)))
)
write.csv(heatmap_genes_h, "Human_RPE_aging_heatmap_genes.csv", row.names = FALSE)
cat("Saved:", nrow(heatmap_genes_h), "genes (Up:", length(human_aging_up), "| Down:", length(human_aging_down), ")\n")




#######
#######
#######
library(SummarizedExperiment)

setwd("/projects/hmz-aging/Metformin/")
bharti_human <- readRDS("mitra_metformin_bulkRNA_iRPE.rds")

meta = colData(bharti_human)
table(meta$treatment)

assay(bharti_human)[1:10,1:10]

##### 用 Human RPE aging DEGs 预测 bharti_human 的 aging score #####

library(GSVA)
library(ggplot2)

# 读入 human aging 基因 list
heatmap_genes_h <- read.csv("/projects/hmz-aging/Metformin/Human_RPE_atlas_mtx/Human_RPE_aging_heatmap_genes.csv")
up_h <- heatmap_genes_h$gene[heatmap_genes_h$direction == "Up"]
down_h <- heatmap_genes_h$gene[heatmap_genes_h$direction == "Down"]

# 表达矩阵是非整数的 counts，round 后用 DESeq2 VST 归一化
library(DESeq2)
counts_human <- round(as.matrix(assay(bharti_human)))
meta_human <- as.data.frame(colData(bharti_human))

dds_human <- DESeqDataSetFromMatrix(countData = counts_human, colData = meta_human, design = ~ 1)
vsd_human <- vst(dds_human, blind = TRUE)
expr_human <- assay(vsd_human)
cat("VST expression range:", round(range(expr_human), 2), "\n")

# 匹配基因
up_found_h <- intersect(up_h, rownames(expr_human))
down_found_h <- intersect(down_h, rownames(expr_human))
cat("Up found:", length(up_found_h), "/", length(up_h), "\n")
cat("Down found:", length(down_found_h), "/", length(down_h), "\n")

# ssGSEA
gs_human <- list(Up_aging = up_found_h, Down_aging = down_found_h)
param_h <- ssgseaParam(exprData = expr_human, geneSets = gs_human, alpha = 0)
ssgsea_h <- gsva(param_h, verbose = FALSE)

aging_score_h <- ssgsea_h["Up_aging", ] - ssgsea_h["Down_aging", ]

score_h <- data.frame(
  sample = colnames(expr_human),
  treatment = meta_human$treatment,
  aging_score = aging_score_h
)
print(score_h[order(score_h$treatment), ])

# 散点图（去掉 POS 和 POS+MET）
score_h <- score_h[!score_h$treatment %in% c("POS", "POS+MET"), ]
score_h$treatment <- factor(score_h$treatment,
                            levels = c("untreated", "MET", "CCHS", "CCHS+MET"))
# 计算配对 p 值
comp_pairs <- list(c("untreated", "MET"), c("CCHS", "CCHS+MET"))
pval_annot <- data.frame(x1 = numeric(), x2 = numeric(), y = numeric(), label = character())
y_range <- diff(range(score_h$aging_score))
y_base <- max(score_h$aging_score) + 0.05 * y_range

for (i in seq_along(comp_pairs)) {
  g1 <- score_h$aging_score[score_h$treatment == comp_pairs[[i]][1]]
  g2 <- score_h$aging_score[score_h$treatment == comp_pairs[[i]][2]]
  pv <- t.test(g1, g2)$p.value
  x1_pos <- which(levels(score_h$treatment) == comp_pairs[[i]][1])
  x2_pos <- which(levels(score_h$treatment) == comp_pairs[[i]][2])
  pval_annot <- rbind(pval_annot, data.frame(
    x1 = x1_pos, x2 = x2_pos,
    y = y_base + (i - 1) * 0.06 * y_range,
    label = paste0("p = ", signif(pv, 2))
  ))
}

npg_colors_h <- c("untreated" = "#E64B35", "MET" = "#4DBBD5",
                   "CCHS" = "#F39B7F", "CCHS+MET" = "#8491B4")

p_human <- ggplot(score_h, aes(x = treatment, y = aging_score, color = treatment)) +
  geom_jitter(size = 2, width = 0.15) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, linewidth = 0.6, color = "black") +
  scale_color_manual(values = npg_colors_h) +
  geom_segment(data = pval_annot, aes(x = x1, xend = x2, y = y, yend = y),
               color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  geom_text(data = pval_annot, aes(x = (x1 + x2) / 2, y = y + 0.05 * y_range, label = label),
            size = 3, inherit.aes = FALSE) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(x = "", y = "Aging Score (ssGSEA)", title = "iRPE Aging Score (Human DEGs)")

ggsave("bharti_human_aging_scoreXXX.png", plot = p_human, width = 4, height = 5, dpi = 300)

##### 配对比较：MET 处理 vs 对照 #####
comparisons <- list(
  c("CCHS", "CCHS+MET"),
  c("untreated", "MET")
)

for (comp in comparisons) {
  g1 <- score_h$aging_score[score_h$treatment == comp[1]]
  g2 <- score_h$aging_score[score_h$treatment == comp[2]]
  
  cat("\n===", comp[1], "vs", comp[2], "===\n")
  cat(comp[1], "mean:", round(mean(g1), 4), "| ", comp[2], "mean:", round(mean(g2), 4), "\n")
  cat("Difference:", round(mean(g2) - mean(g1), 4), "\n")
  
  wt <- wilcox.test(g1, g2)
  cat("Wilcoxon p =", signif(wt$p.value, 3), "\n")
  
  tt <- t.test(g1, g2)
  cat("t-test p =", signif(tt$p.value, 3), "\n")
}




#######
#######
#######

# ##### Human RPE aging DEGs GO 分析 #####

# setwd("/projects/hmz-aging/Metformin/Human_RPE_atlas_mtx")

# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ggplot2)

# heatmap_genes_h <- read.csv("Human_RPE_aging_heatmap_genes.csv")
# up_symbols_h <- heatmap_genes_h$gene[heatmap_genes_h$direction == "Up"]
# down_symbols_h <- heatmap_genes_h$gene[heatmap_genes_h$direction == "Down"]

# # 使用 MSigDB Hallmark + GO BP + KEGG
# library(msigdbr)

# hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
# gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")

# msig_all <- rbind(hallmark, gobp)
# msig_t2g <- msig_all[, c("gs_name", "gene_symbol")]

# # enricher ORA - Up genes
# enr_up_h <- enricher(gene = up_symbols_h,
#                      TERM2GENE = msig_t2g,
#                      pAdjustMethod = "BH",
#                      pvalueCutoff = 0.05,
#                      qvalueCutoff = 0.2,
#                      minGSSize = 5)

# # enricher ORA - Down genes
# enr_down_h <- enricher(gene = down_symbols_h,
#                        TERM2GENE = msig_t2g,
#                        pAdjustMethod = "BH",
#                        pvalueCutoff = 0.05,
#                        qvalueCutoff = 0.2,
#                        minGSSize = 5)

# cat("Up enriched terms:", nrow(as.data.frame(enr_up_h)), "\n")
# cat("Down enriched terms:", nrow(as.data.frame(enr_down_h)), "\n")

# # 过滤和画图函数
# filter_enr <- function(enr_res, ratio_cutoff = 0.03) {
#   df <- as.data.frame(enr_res)
#   if (nrow(df) == 0) return(df)
#   df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
#     parts <- as.numeric(strsplit(x, "/")[[1]])
#     parts[1] / parts[2]
#   })
#   df <- df[df$GeneRatioNum > ratio_cutoff, ]
#   df$Description <- gsub("^HALLMARK_|^GOBP_|^KEGG_", "", df$ID)
#   df$Description <- gsub("_", " ", df$Description)
#   df$Description <- tolower(df$Description)
#   substr(df$Description, 1, 1) <- toupper(substr(df$Description, 1, 1))
#   df
# }

# nature_msig_plot <- function(df, title, top_n = 15) {
#   if (nrow(df) == 0) { cat("No terms for:", title, "\n"); return(NULL) }
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

# enr_up_df <- filter_enr(enr_up_h)
# enr_down_df <- filter_enr(enr_down_h)

# p_up_h <- nature_msig_plot(enr_up_df, "Up with Aging (Human RPE)")
# p_down_h <- nature_msig_plot(enr_down_df, "Down with Aging (Human RPE)")

# if (!is.null(p_up_h)) ggsave("Human_MSigDB_up_aging.png", plot = p_up_h, width = 8, height = 5, dpi = 300)
# if (!is.null(p_down_h)) ggsave("Human_MSigDB_down_aging.png", plot = p_down_h, width = 8, height = 5, dpi = 300)




##### CCHS vs CCHS+MET 差异基因，找 MET 逆转衰老的基因 #####

library(limma)
library(edgeR)

counts_bh <- round(as.matrix(assay(bharti_human)))
meta_bh <- as.data.frame(colData(bharti_human))

# 只取 CCHS 和 CCHS+MET
idx_cchs <- meta_bh$treatment %in% c("CCHS", "CCHS+MET")
counts_cchs <- counts_bh[, idx_cchs]
meta_cchs <- meta_bh[idx_cchs, ]
group_cchs <- factor(meta_cchs$treatment, levels = c("CCHS", "CCHS+MET"))

dge_cchs <- DGEList(counts = counts_cchs)
dge_cchs <- calcNormFactors(dge_cchs, method = "TMM")
design_cchs <- model.matrix(~ group_cchs)
v_cchs <- voom(dge_cchs, design_cchs, plot = FALSE)
fit_cchs <- lmFit(v_cchs, design_cchs)
fit_cchs <- eBayes(fit_cchs)

cchs_degs <- topTable(fit_cchs, coef = 2, number = Inf, sort.by = "P")
cchs_degs$gene <- rownames(cchs_degs)

cat("=== CCHS vs CCHS+MET ===\n")
cat("Total genes:", nrow(cchs_degs), "\n")
cat("p < 0.05:", sum(cchs_degs$P.Value < 0.05), "\n")
cat("p < 0.1:", sum(cchs_degs$P.Value < 0.1), "\n")

cchs_degs_sig <- cchs_degs[cchs_degs$P.Value < 0.05, ]
cat("Sig DEGs (p<0.05):", nrow(cchs_degs_sig), "| Up:", sum(cchs_degs_sig$logFC > 0), "| Down:", sum(cchs_degs_sig$logFC < 0), "\n")

#### 与 human aging DEGs 比较: 找方向相反的 ####

aging_degs_h <- read.csv("/projects/hmz-aging/Metformin/Human_RPE_atlas_mtx/Human_RPE_aging_heatmap_genes.csv")

overlap_cchs <- merge(cchs_degs_sig, aging_degs_h, by.x = "gene", by.y = "gene")
overlap_cchs$met_direction <- ifelse(overlap_cchs$logFC > 0, "Up", "Down")
opposite_cchs <- overlap_cchs[overlap_cchs$met_direction != overlap_cchs$direction, ]
same_cchs <- overlap_cchs[overlap_cchs$met_direction == overlap_cchs$direction, ]

cat("\nOverlap total:", nrow(overlap_cchs), "\n")
cat("Opposite direction (MET reverses aging):", nrow(opposite_cchs), "\n")
cat("  Aging Up + Met Down:", sum(opposite_cchs$direction == "Up"), "\n")
cat("  Aging Down + Met Up:", sum(opposite_cchs$direction == "Down"), "\n")
cat("Same direction:", nrow(same_cchs), "\n")

if (nrow(opposite_cchs) > 0) print(opposite_cchs[, c("gene", "logFC", "P.Value", "direction", "met_direction")])

write.csv(opposite_cchs, "CCHS_MET_reverses_aging_genes.csv", row.names = FALSE)


##### Heatmap: CCHS+MET 逆转衰老基因 #####

library(ComplexHeatmap)
library(circlize)

# opposite 基因按 aging 方向分组
opp_genes <- opposite_cchs$gene
opp_aging_up <- opposite_cchs$gene[opposite_cchs$direction == "Up"]
opp_aging_down <- opposite_cchs$gene[opposite_cchs$direction == "Down"]
all_opp <- c(opp_aging_up, opp_aging_down)

# 左边: Metformin bulk (CCHS / CCHS+MET)
cchs_samples <- rownames(meta_bh)[meta_bh$treatment == "CCHS"]
cchs_met_samples <- rownames(meta_bh)[meta_bh$treatment == "CCHS+MET"]
met_col_order <- c(cchs_samples, cchs_met_samples)

all_opp_in_expr <- all_opp[all_opp %in% rownames(expr_human)]
met_mat_h <- expr_human[all_opp_in_expr, met_col_order]
met_z_h <- t(scale(t(met_mat_h)))

# 中间/右边: Aging pseudobulk (male / female)
all_opp_in_pb <- all_opp_in_expr[all_opp_in_expr %in% rownames(pb_norm)]

male_idx <- which(pb_seurat$Gender == "Male")
female_idx <- which(pb_seurat$Gender == "Female")
male_order <- male_idx[order(as.numeric(pb_seurat$Age[male_idx]))]
female_order <- female_idx[order(as.numeric(pb_seurat$Age[female_idx]))]

common_genes <- intersect(all_opp_in_expr, all_opp_in_pb)
common_up <- common_genes[common_genes %in% opp_aging_up]
common_down <- common_genes[common_genes %in% opp_aging_down]
common_ordered <- c(common_up, common_down)
cat("Heatmap genes:", length(common_ordered), "(Up:", length(common_up), "| Down:", length(common_down), ")\n")

# Z-score matrices
met_z_final <- t(scale(t(expr_human[common_ordered, met_col_order])))
male_z_final <- t(scale(t(pb_norm[common_ordered, male_order])))
female_z_final <- t(scale(t(pb_norm[common_ordered, female_order])))

# 行分组
row_split_opp <- factor(c(rep("Aging Up\n(Met Down)", length(common_up)),
                           rep("Aging Down\n(Met Up)", length(common_down))),
                         levels = c("Aging Up\n(Met Down)", "Aging Down\n(Met Up)"))

# 左边列注释
met_group <- c(rep("CCHS", length(cchs_samples)), rep("CCHS+MET", length(cchs_met_samples)))
anno_met <- HeatmapAnnotation(
  Group = met_group,
  col = list(Group = c("CCHS" = "#E64B35", "CCHS+MET" = "#4DBBD5")),
  show_annotation_name = FALSE
)

# 中间/右边 Age 注释
anno_male <- HeatmapAnnotation(
  Age = as.numeric(pb_seurat$Age[male_order]),
  col = list(Age = colorRamp2(c(0, 50, 100), c("#4DBBD5", "white", "#E64B35"))),
  show_annotation_name = FALSE
)
anno_female <- HeatmapAnnotation(
  Age = as.numeric(pb_seurat$Age[female_order]),
  col = list(Age = colorRamp2(c(0, 50, 100), c("#4DBBD5", "white", "#E64B35"))),
  show_annotation_name = FALSE
)

show_rn <- length(common_ordered) <= 80

ht_met <- Heatmap(met_z_final, name = "Met\nZ-score",
                  col = colorRamp2(c(-2, 0, 2), c("#00A087", "white", "#DC0000")),
                  cluster_columns = FALSE, cluster_rows = FALSE,
                  row_split = row_split_opp, row_gap = unit(3, "mm"),
                  row_title_gp = gpar(fontsize = 14),
                  column_title = "Metformin iRPE", column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                  top_annotation = anno_met,
                  show_row_names = show_rn, row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                  row_names_side = "left",
                  show_column_names = FALSE)

ht_male <- Heatmap(male_z_final, name = "Male\nZ-score",
                   col = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")),
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   row_split = row_split_opp, row_gap = unit(3, "mm"),
                   column_title = "Aging Male", column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                   top_annotation = anno_male,
                   show_row_names = FALSE,
                   show_column_names = FALSE)

# 右侧标注关键基因
highlight_genes <- list(
  "Inflammatory / NF-kB" = c("OSMR","IL6R","RELB","NFKBIA","CEBPD","LGALS3","SERPINA3","SERPING1","CHI3L1","ICAM1","TNFAIP1"),
  "Epithelial remodeling" = c("ITGA3","ITGB4","PXN","ACTN4","PKP2","PVR","FOSL2","MAFK","TGIF1","AHNAK"),
  "Mito / OXPHOS / FAO" = c("NDUFB9","NDUFB2","UQCRC2","MRPL1","MRPL15","MRPL16","MRPL19","MRPL42","MCUR1","ACADM","ACADVL","GOT2","FH"),
  "Vesicle / lysosome" = c("GULP1","VPS26C","AP3M1","STXBP4","RAB4A","SCAMP1","SCAMP3","PPT1"),
  "RPE homeostasis" = c("CLDN10","COL4A5","NECTIN3","LRP8","AMOTL1","ESRRG","TGFBR1","SMAD2","SMAD5","NOTCH2","ARRB1","ERBB4","PDGFD")
)

all_highlight <- unlist(highlight_genes)
gene_to_category <- rep(names(highlight_genes), sapply(highlight_genes, length))
names(gene_to_category) <- all_highlight

# 构建右侧 row annotation
mark_idx <- which(common_ordered %in% all_highlight)
mark_labels <- common_ordered[mark_idx]

anno_right <- rowAnnotation(
  Genes = anno_mark(
    at = mark_idx,
    labels = mark_labels,
    labels_gp = gpar(fontsize = 11, fontface = "italic"),
    link_width = unit(3, "mm")
  )
)

ht_female <- Heatmap(female_z_final, name = "Female\nZ-score",
                     col = colorRamp2(c(-2, 0, 2), c("#00A087", "white", "#DC0000")),
                     cluster_columns = FALSE, cluster_rows = FALSE,
                     row_split = row_split_opp, row_gap = unit(3, "mm"),
                     column_title = "Aging Female", column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                     top_annotation = anno_female,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     right_annotation = anno_right)

ht_height <- min(14, max(4, length(common_ordered) * 0.08 + 1.5))
png("CCHS_MET_reverses_aging_heatmapXXX.png", width = 8, height = ht_height, units = "in", res = 300)
draw(ht_met + ht_male + ht_female, ht_gap = unit(3, "mm"))
dev.off()

######
######
######
