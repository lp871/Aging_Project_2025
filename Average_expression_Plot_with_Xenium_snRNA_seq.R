############
######
############

######
###### load the smoothed RNA expression ######
######

setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Xenium5k_Mouse_Obj_pseudobulk_list <- readRDS("Xenium5k_Mouse_Obj_pseudobulk_list_May2025")


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat_knn_pseudobulk <- readRDS("Mouse_seurat_knn_pseudobulk_2025")

####
####

setwd("")


####
####

Genes = c("Grik3","C4b","Cd36","A2m","Wnt6","Ifi27")

#####
#####

Xenium5k_MG = Xenium5k_Mouse_Obj_pseudobulk_list$MG
snRNAseq_MG = Mouse_seurat_knn_pseudobulk$MG

######
###### seurat_obj = Xenium5k_MG
######

Xenium5k_MG_Gexp = Get_Gene_expressions(Genes,Xenium5k_MG)
Xenium5k_MG_Gexp$age = sapply(strsplit(as.character(Xenium5k_MG_Gexp$cell),split="_"),function(x) x[[1]])

snRNAseq_MG_Gexp = Get_Gene_expressions(Genes,snRNAseq_MG)
snRNAseq_MG_Gexp$age = sapply(strsplit(as.character(snRNAseq_MG_Gexp$cell),split="_"),function(x) x[[1]])

#######
library(ggplot2)

Xenium5k_MG_Gexp$age = factor(Xenium5k_MG_Gexp$age,levels=c('5','12','52','112'))
snRNAseq_MG_Gexp$age = factor(snRNAseq_MG_Gexp$age,levels=c('5','12','17','32','49','68','91','108','120'))

Xenium5k_MG_Gexp$gene = factor(Xenium5k_MG_Gexp$gene,Genes)
snRNAseq_MG_Gexp$gene = factor(snRNAseq_MG_Gexp$gene,Genes)


p <- ggplot(Xenium5k_MG_Gexp, aes(x = age, y = expression, fill=age)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  facet_wrap(~ gene, nrow = 1) +
  # 给绘图区加黑色边框
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# 显示并保存
print(p)
ggsave("Xenium5k_MG_Gexp.png", plot = p, width = 16, height = 3, dpi = 300)


#################
#################
#################

p <- ggplot(snRNAseq_MG_Gexp, aes(x = age, y = expression, fill=age)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  facet_wrap(~ gene, nrow = 1) +
  # 给绘图区加黑色边框
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, 2))

# 显示并保存
print(p)
ggsave("snRNAseq_MG_Gexp.png", plot = p, width = 16, height = 3, dpi = 300)


####


Get_Gene_expressions <- function(Genes,seurat_obj){
    ##########
    layer_list <- seurat_obj[["RNA"]]@layers
    layer_list_merge = do.call("cbind",layer_list)
    ##########
    rownames(layer_list_merge) = rownames(seurat_obj)
    colnames(layer_list_merge) = colnames(seurat_obj)
    mat = layer_list_merge
    total_counts <- colSums(mat)
    scale.factor = 1e4
    mat_norm <- sweep(mat, 2, total_counts, FUN = "/") * scale.factor
    mat_norm <- log2(mat_norm + 1)
    ##########
    genes_index = which(rownames(mat_norm) %in% Genes == T)
    ########## rownames(layer_list_merge)[genes_index]
    ########## rownames(layer_list_merge)[grep("fi27",rownames(layer_list_merge))]
    layer_list_merge_cl = mat_norm[genes_index,]
    ##########
    library(reshape2)
    long_df <- melt(
        as.matrix(layer_list_merge_cl),
        varnames   = c("gene", "cell"),
        value.name = "expression"
    )
    ##########
    return(long_df)
}













