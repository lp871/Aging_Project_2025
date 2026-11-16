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

# Nature-style color palette
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

# Unname and deduplicate to avoid DotPlot faceting / duplicate row-name errors caused by named vectors
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



##### Find DEGs below ####
##### Draw heatmaps ######
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

# Consistent direction: 5->49, 49->91 and 5->91 have the same sign
consistent_genes <- df_merge %>%
  filter(sign(logFC_5v49) == sign(logFC_49v91) & sign(logFC_5v49) == sign(logFC_5v91))

cat("Consistent logFC genes:", nrow(consistent_genes), "\n")

####### 428 genes #######

all_consistent <- consistent_genes$gene

# Pseudobulk: aggregate all genes -> normalize by whole-transcriptome library size -> log2 -> then extract DEGs
counts_mat <- GetAssayData(Mouse_RPE_cl, assay = "RNA", layer = "counts")

sample_ids <- Mouse_RPE_cl$sample
# Aggregate all genes into pseudobulk by age group
bulk_counts_all <- sapply(levels(sample_ids), function(s) {
  rowSums(counts_mat[, sample_ids == s, drop = FALSE])
})

# Library size computed from all genes (sum of counts across all genes per group)
lib_size <- colSums(bulk_counts_all)
cpm_all <- t(t(bulk_counts_all) / lib_size * 1e4)
mat_all <- log2(cpm_all + 1)

# Extract DEGs after normalization and fix the age order
mat <- mat_all[all_consistent, c("5W_RPE", "49W_RPE", "91W_RPE")]

# Pseudobulk difference direction must also agree with logFC
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

# Update mat to keep only the filtered genes
mat <- mat[c(up_genes, down_genes), ]

# Z-score per gene (row)
mat_z <- t(scale(t(mat)))

# Up/down group annotation
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

# Save DEGs shown in the heatmap
heatmap_degs <- data.frame(
  gene = c(up_genes, down_genes),
  direction = c(rep("Up", length(up_genes)), rep("Down", length(down_genes))),
  gene_symbol = sub(".*~~", "", c(up_genes, down_genes))
)


write.csv(heatmap_degs, "RPE_aging_consistent_DEGs.csv", row.names = FALSE)
cat("Saved", nrow(heatmap_degs), "DEGs to RPE_aging_consistent_DEGs.csv\n")




# ##### GO analysis #####

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

# # Filter GeneRatio > 0.03
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




##### Compute aging score using heatmap DEGs (excluding mt genes) #####
heatmap_degs <- read.csv("RPE_aging_consistent_DEGs.csv")
mt_drop <- grepl("^mt-", heatmap_degs$gene_symbol, ignore.case = TRUE)
cat("dropped mt aging features:", paste(heatmap_degs$gene_symbol[mt_drop], collapse = ", "),
    "(", sum(mt_drop), ")\n")
heatmap_degs <- heatmap_degs[!mt_drop, ]

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






# Ordering: con/cse pairs, 3d->6d->10d->10d_old, RPE first, then SickRPE
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
# Keep only levels that are present
level_order <- level_order[level_order %in% unique(smoking_RPE_cl$celltype)]

smoking_RPE_cl$celltype <- factor(smoking_RPE_cl$celltype, levels = level_order)





# Keep only the specified samples
keep_ct <- c("con_3d_RPE", "cse_3d_RPE",
             "con_6d_RPE", "cse_6d_RPE", "cse_6d_SickRPE",
             "con_10d_old_RPE", "cse_10d_old_RPE", "con_10d_old_SickRPE")
smoking_RPE_clcl <- smoking_RPE_cl[, smoking_RPE_cl$celltype %in% keep_ct]
smoking_RPE_clcl$celltype <- droplevels(factor(smoking_RPE_clcl$celltype, levels = keep_ct))
table(smoking_RPE_clcl$celltype)

# Re-run UMAP
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



##### Grouped violin plot + significance tests #####

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
    pv <- wilcox.test(test_vals, con_vals, alternative = "two.sided")$p.value
    label <- ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", "ns")))
    data.frame(group1 = con_ct, group2 = x, p = pv, label = label)
  }))

  y_max <- max(df$aging_score, na.rm = TRUE)
  y_step <- (y_max - min(df$aging_score, na.rm = TRUE)) * 0.1
  pval_df$x1 <- match(pval_df$group1, celltypes)
  pval_df$x2 <- match(pval_df$group2, celltypes)
  # Leave extra space above the violin top to avoid overlap between brackets and the plot
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
##### Smoking & Aging: pseudobulk direction only (no DEG test, mt genes removed) #####
########
# Aging features: RPE_aging_consistent_DEGs.csv without ^mt- genes.
# Pseudobulk log2(CP10k+1) differences:
#   3d:  cse_3d_RPE          - con_3d_RPE
#   6d:  cse_6d_SickRPE      - con_6d_RPE
#   10d: con_10d_old_SickRPE - con_10d_old_RPE
# Keep genes whose three differences share the same non-zero sign and match the
# aging direction; Table 10 genes additionally require |mean difference| > 0.2.

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
setwd("/projects/hmz-aging/RPE_smoking")

heatmap_degs <- read.csv("RPE_aging_consistent_DEGs.csv", stringsAsFactors = FALSE)
heatmap_degs <- heatmap_degs[!grepl("^mt-", heatmap_degs$gene_symbol, ignore.case = TRUE), ]
heatmap_degs <- heatmap_degs[!duplicated(heatmap_degs$gene_symbol), ]
cat("Aging features (no mt):", nrow(heatmap_degs), "\n")

# Smoking pseudobulk
load(file.path("Handa_snRNA_Final 1.RData"))
smoking_RPE <- combined
rm(combined); invisible(gc())
cluster_num <- as.numeric(sub(".*_", "", smoking_RPE$group_1))
smoking_RPE$celltype <- NA
smoking_RPE$celltype[cluster_num == 1] <- "RPE"
smoking_RPE$celltype[cluster_num == 9] <- "SickRPE"
smoking_RPE$celltype <- paste(smoking_RPE$group, smoking_RPE$celltype, sep = "_")

pb_cols <- c("con_3d_RPE", "cse_3d_RPE",
             "con_6d_RPE", "cse_6d_SickRPE",
             "con_10d_old_RPE", "con_10d_old_SickRPE")
smoking_pb <- smoking_RPE[, smoking_RPE$celltype %in% pb_cols]
DefaultAssay(smoking_pb) <- "RNA"
if (inherits(smoking_pb[["RNA"]], "Assay5") && length(Layers(smoking_pb[["RNA"]], search = "counts")) > 1)
  smoking_pb[["RNA"]] <- JoinLayers(smoking_pb[["RNA"]])
counts_smoke <- GetAssayData(smoking_pb, assay = "RNA", layer = "counts")
bulk_smoke <- sapply(pb_cols, function(ct) {
  Matrix::rowSums(counts_smoke[, smoking_pb$celltype == ct, drop = FALSE])
})
mat_smoke <- log2(t(t(bulk_smoke) / colSums(bulk_smoke) * 1e4) + 1)
rm(smoking_RPE, smoking_pb, counts_smoke, bulk_smoke); invisible(gc())

genes <- heatmap_degs$gene_symbol[heatmap_degs$gene_symbol %in% rownames(mat_smoke)]
cat("Aging features present in smoking matrix:", length(genes),
    "| missing:", sum(!heatmap_degs$gene_symbol %in% rownames(mat_smoke)), "\n")

d3 <- mat_smoke[genes, "cse_3d_RPE"] - mat_smoke[genes, "con_3d_RPE"]
d6 <- mat_smoke[genes, "cse_6d_SickRPE"] - mat_smoke[genes, "con_6d_RPE"]
d10 <- mat_smoke[genes, "con_10d_old_SickRPE"] - mat_smoke[genes, "con_10d_old_RPE"]
same <- sign(d3) == sign(d6) & sign(d6) == sign(d10) & sign(d3) != 0

pb_dir <- data.frame(
  gene = genes[same],
  direction = ifelse(d3[same] > 0, "Up", "Down"),
  direction_aging = heatmap_degs$direction[match(genes[same], heatmap_degs$gene_symbol)],
  d_3d = d3[same],
  d_6d = d6[same],
  d_10d = d10[same],
  stringsAsFactors = FALSE
)
pb_dir$mean_d <- rowMeans(pb_dir[, c("d_3d", "d_6d", "d_10d")])
pb_dir <- pb_dir[order(pb_dir$direction != "Up", -abs(pb_dir$mean_d)), ]
write.csv(pb_dir, "pseudobulk_same_direction_genes.csv", row.names = FALSE)

same_direction <- pb_dir %>% filter(direction == direction_aging)
opposite_direction <- pb_dir %>% filter(direction != direction_aging)
cat("Same pseudobulk direction:", nrow(pb_dir), "\n")
cat("  Match aging direction:", nrow(same_direction),
    "(Up", sum(same_direction$direction == "Up"),
    "| Down", sum(same_direction$direction == "Down"), ")\n")
cat("  Opposite aging direction:", nrow(opposite_direction), "\n")
write.csv(same_direction, "Smoking_Aging_same_direction_genes.csv", row.names = FALSE)
write.csv(opposite_direction, "Smoking_Aging_opposite_direction_genes.csv", row.names = FALSE)

mean_d_cutoff <- 0.2
table10_genes <- same_direction %>% filter(abs(mean_d) > mean_d_cutoff)
cat("Same as aging & |mean log2FC| >", mean_d_cutoff, ":", nrow(table10_genes), "\n")
print(table10_genes$gene)
write.csv(table10_genes, "pseudobulk_meanLog2FC0.2_genes.csv", row.names = FALSE)

# Aging pseudobulk
Mouse_RPE <- readRDS(file.path("Mouse_RPE_RNA_Final_Jul16"))
Mouse_RPE_cl <- Mouse_RPE[, Mouse_RPE$celltype %in% c("RPE")]
Mouse_RPE_cl$sample <- factor(Mouse_RPE_cl$sample, levels = c("5W_RPE", "49W_RPE", "91W_RPE"))
counts_aging <- GetAssayData(Mouse_RPE_cl, assay = "RNA", layer = "counts")
bulk_aging <- sapply(levels(Mouse_RPE_cl$sample), function(s) {
  Matrix::rowSums(counts_aging[, Mouse_RPE_cl$sample == s, drop = FALSE])
})
mat <- log2(t(t(bulk_aging) / colSums(bulk_aging) * 1e4) + 1)
rm(Mouse_RPE, Mouse_RPE_cl, counts_aging, bulk_aging); invisible(gc())
aging_gene_map <- setNames(sub(".*~~", "", rownames(mat)), rownames(mat))

plot_smoking_aging_heatmap <- function(gene_df, file, title) {
  gene_df <- gene_df[gene_df$gene %in% aging_gene_map, ]
  up_sym_ht <- gene_df$gene[gene_df$direction == "Up"]
  down_sym_ht <- gene_df$gene[gene_df$direction == "Down"]
  ordered_sym <- c(up_sym_ht, down_sym_ht)
  row_split_ht <- factor(c(rep("Up", length(up_sym_ht)), rep("Down", length(down_sym_ht))),
                         levels = c("Up", "Down"))

  z_smoke <- t(scale(t(mat_smoke[ordered_sym, pb_cols, drop = FALSE])))
  aging_rows <- names(aging_gene_map)[match(ordered_sym, aging_gene_map)]
  aging_sub <- mat[aging_rows, c("5W_RPE", "49W_RPE", "91W_RPE"), drop = FALSE]
  rownames(aging_sub) <- ordered_sym
  aging_z <- t(scale(t(aging_sub)))

  fs <- if (length(ordered_sym) > 80) 5 else if (length(ordered_sym) > 40) 7 else 10
  smoke_group <- c("Control", "Smoking", "Control", "Sick", "Control", "Sick")
  col_anno_smoke <- HeatmapAnnotation(
    Group = smoke_group,
    col = list(Group = c(Control = "#4DBBD5", Smoking = "#E64B35", Sick = "#7E6148")),
    annotation_name_gp = gpar(fontsize = 10)
  )

  ht_smoke <- Heatmap(z_smoke, name = "Smoking\nZ-score",
                      col = colorRamp2(c(-2, 0, 2), c("#00A087", "white", "#E64B35")),
                      rect_gp = gpar(col = "white", lwd = 0.5),
                      top_annotation = col_anno_smoke,
                      cluster_columns = FALSE, cluster_rows = FALSE,
                      row_split = row_split_ht, row_gap = unit(3, "mm"),
                      row_title = c("Up", "Down"), row_title_gp = gpar(fontsize = 12),
                      column_title = "Smoking (pseudobulk)", column_title_gp = gpar(fontsize = 12),
                      show_row_names = TRUE, row_names_gp = gpar(fontsize = fs, fontface = "italic"),
                      row_names_side = "left",
                      column_names_rot = 45, column_names_gp = gpar(fontsize = 11))
  ht_aging <- Heatmap(aging_z, name = "Aging\nZ-score",
                      col = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#DC0000")),
                      rect_gp = gpar(col = "white", lwd = 0.5),
                      cluster_columns = FALSE, cluster_rows = FALSE,
                      row_split = row_split_ht, row_gap = unit(3, "mm"),
                      column_title = "Aging (5W -> 49W -> 91W)", column_title_gp = gpar(fontsize = 12),
                      show_row_names = FALSE,
                      column_names_rot = 45, column_names_gp = gpar(fontsize = 11))

  png(file, width = 10, height = max(5, length(ordered_sym) * 0.16 + 3), units = "in", res = 300)
  draw(ht_smoke + ht_aging, column_title = title,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       ht_gap = unit(8, "mm"))
  dev.off()

  data.frame(gene = ordered_sym,
             category = c(rep("Up with aging & smoking", length(up_sym_ht)),
                          rep("Down with aging & smoking", length(down_sym_ht))))
}

smoke_hm_genes <- plot_smoking_aging_heatmap(
  same_direction, "Smoking_Aging_pseudobulk_direction_heatmap.png",
  "Aging features, same pseudobulk direction as aging, no mt")
write.csv(smoke_hm_genes, "Smoking_heatmap_genes_category.csv", row.names = FALSE)

table10_hm_genes <- plot_smoking_aging_heatmap(
  table10_genes, "Smoking_Aging_pseudobulk_meanLog2FC0.2_heatmap.png",
  "Same direction as aging, |mean log2FC| > 0.2, no mt")
cat("Heatmap genes:", nrow(smoke_hm_genes), "| |mean log2FC| > 0.2:", nrow(table10_hm_genes), "\n")





