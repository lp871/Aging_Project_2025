##########
########## process the OSK samples, UMAP and Other #####
##########



ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4

R
library(Seurat)
setwd("/zp1/data/plyu3/Xenium5k_Mouse/OSK_New")
OSK_obj <- readRDS("OSK_obj_2024Dec12")


#########
#########
#########

table(OSK_obj$celltype)

####### rm low number of cell types ######
#######

rm = c("Schwann","Adipose Cells","Corneal Epithelial")

OSK_obj_cl = OSK_obj[,which(OSK_obj$celltype %in% rm == F)]

table(OSK_obj_cl$celltype)


cell_colors <- c(
  "AC"                  = "#E41A1C",
  "Astrocyte"           = "#66C2A5",
  "BC"                  = "#1F78B4",
  "Cone"                = "#FF7F00",
  "Fibroblasts"         = "#66A61E",
  "HC"                  = "#01796F",  # updated to deep teal
  "Keratocytes"         = "#1B9E77",
  "Macrophages"         = "#D95F02",
  "Melanocyte"          = "#7570B3",
  "MG"                  = "#A04000",  # updated to warm sienna
  "Pericyte"            = "#FB9A99",
  "RGC"                 = "#FFC20A",  # updated to richer gold
  "Rod"                 = "#A6CEE3",
  "RPE"                 = "black",
  "smooth muscle cell"  = "#B2DF8A",
  "VE"                  = "#E31A1C"
)

########
######## Plot the UMAP plots ######
########


umap_df <- seuratToUMAP(OSK_obj_cl,id_col='celltype',sample_col = "sample")

########
########
########

library(ggplot2)

ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype, fill = celltype)) +
  geom_point(
    shape  = 21,
    size   = 0.3,
    stroke = 0.2,
    alpha  = 0.1
  ) +
  scale_color_manual(values = cell_colors) +
  scale_fill_manual(values  = cell_colors, guide = "none") + 
  # hide the separate fill scale
  guides(
    color = guide_legend(
      override.aes = list(
        alpha = 1,
        shape  = 21,
        size   = 3,
        stroke = 0.2,
        fill   = cell_colors  # enforce fill colors in legend
      )
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  ) +
  labs(
    x     = "UMAP 1",
    y     = "UMAP 2"
  )


ggsave("test1.png",height=5,width=8)


####
####
####


umap_df$sample <- factor(
  umap_df$sample,
  levels = c("OSKR1", "OSKR2", "OSKR3", "OSKR4"),
  labels = c("Control", "Control_Dox", "Control_Diet", "Diet")
)

#####
#####
sample_colors <- c(
  "Control"   = "#E64B35",
  "Control_Dox"  = "#4DBBD5",
  "Control_Diet"  = "#00A087",
  "Diet" = "#3C5488"
)

ggplot(umap_df, aes(x = umap_1, y = umap_2, color = sample)) +
  geom_point(
    shape = 16,        # 实心圆
    size  = 0.1,       # 点大小
    alpha = 0.25       # 点透明度
  ) +
  scale_color_manual(values = sample_colors) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 16,
        size  = 4,
        alpha = 1
      )
    )
  ) +
  facet_wrap(~ sample, nrow = 2, ncol = 2) +   # 按 sample 拆分为 2×2 面板
  theme_classic(base_size = 14) +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10),
    strip.background = element_rect(fill = "grey90", color = "black"),  # 面板标题背景
    strip.text       = element_text(size = 12)                         # 面板标题字体
  ) +
  labs(
    x     = "UMAP 1",
    y     = "UMAP 2",
    color = "Sample"
  )
ggsave("test2.png",height=5,width=8)


###########
###########

画一下 空间分布 #####
###########


source("/zp1/data/plyu3/All_Functions_2025/R/source_all.R")
source_all_r("/zp1/data/plyu3/All_Functions_2025/R/")

############ load the total of the OSK R1 R2 R3 R4 #######
############

setwd("/zp1/data/plyu3/Xenium5k_Mouse/OSK_New")
R1 = readRDS("OSK_R1_xenium.obj.total.rds")
R2 = readRDS("OSK_R2_xenium.obj.total.rds")
R3 = readRDS("OSK_R3_xenium.obj.total.rds")
R4 = readRDS("OSK_R4_xenium.obj.total.rds")

R1_region = Get_regions_dims_from_seurat(R1)
R2_region = Get_regions_dims_from_seurat(R2)
R3_region = Get_regions_dims_from_seurat(R3)
R4_region = Get_regions_dims_from_seurat(R4)

R1_region$Cell_ID <- paste0("OSKR1_",R1_region$Cell_ID)
R2_region$Cell_ID <- paste0("OSKR2_",R2_region$Cell_ID)
R3_region$Cell_ID <- paste0("OSKR3_",R3_region$Cell_ID)
R4_region$Cell_ID <- paste0("OSKR4_",R4_region$Cell_ID)

M_region = rbind(R1_region,R2_region,R3_region,R4_region)
OSK_obj_cl <- Add_regions_dims_to_seurat(M_region,OSK_obj_cl)

#######
#######
#######

embedding <- Embeddings(OSK_obj_cl, "spatial")
df <- cbind(OSK_obj_cl@meta.data,as.data.frame(embedding))
df_sp = split(df,df$sample)

names(df_sp)
df_sp_new = list()

for(i in 1:length(df_sp)){
    #####
    tmp = df_sp[[i]]
    print(dim(tmp))
    #m = match(tmp$cells,colnames(Xenium5k_Mouse_Obj))
    #print(summary(m))
    #tmp$celltype = Xenium5k_Mouse_Obj$celltype[m]
    #tmp_cl = tmp[is.na(tmp$celltype) == F,]
    #print(dim(tmp_cl))
    keep_types <- c("AC", "BC", "Rod", "Cone", "RGC", "MG", "HC", "RPE")
    #####
    tmp$celltype <- ifelse(tmp$celltype %in% keep_types, as.character(tmp$celltype), "Other")
    tmp$celltype = as.factor(tmp$celltype)
    #####
    df_sp_new = c(df_sp_new,list(tmp))
    #####
}

table(df_sp_new[[1]]$celltype)

cell_colors <- c(
  "AC"                  = "#E41A1C",
  "BC"                  = "#1F78B4",
  "Cone"                = "#FF7F00",
  "HC"                  = "#01796F",  # updated to deep teal
  "MG"                  = "#A04000",  # updated to warm sienna
  "RGC"                 = "#FFC20A",  # updated to richer gold
  "Rod"                 = "#A6CEE3",
  "RPE"                 = "black",
  "Other"               = "grey"
)



Plot = df_sp_new[[4]]

table(Plot$celltype)

ggplot(Plot, aes(x = spatial_1, y = spatial_2, color = celltype, fill = celltype)) +
  geom_point(
    shape  = 21,
    size   = 0.2,
    stroke = 0.2,
    alpha  = 1
  ) +
  scale_color_manual(values = cell_colors) +
  scale_fill_manual(values  = cell_colors) + 
  # hide the separate fill scale
  guides(
    color = guide_legend(
      override.aes = list(
        alpha = 1,
        shape  = 21,
        size   = 3,
        stroke = 0.2
      )
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  ) +
  labs(
    x     = "",
    y     = ""
  )


ggsave("R4.png",height=15,width=15)

#####
#####
#####

setwd("/zp1/data/plyu3/Xenium5k_Mouse/OSK_New")
saveRDS(OSK_obj_cl,file="OSK_obj_cl_May2025.rds")


######
######


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4

R
library(Seurat)
setwd("/zp1/data/plyu3/Xenium5k_Mouse/OSK_New")

OSK_obj_cl = readRDS("OSK_obj_cl_May2025.rds")

#####
all_genes_5k = rownames(OSK_obj_cl)

k = which(all_genes_5k %in% c("Oct4","Sox2","Klf4") == T)

all_genes_5k[k]

##### No Oct4 ####
#####

DefaultAssay(OSK_obj_cl) = "SCT"
OSK_obj_cl <- NormalizeData(OSK_obj_cl)

######
getGeneStatsBySampleCelltype <- function(
  seu,
  gene,
  sample_col    = "sample",
  celltype_col  = "celltype",
  assay         = "SCT",
  slot          = "data"
) {
  # 基因检查
  if (!gene %in% rownames(seu@assays[[assay]]@data)) {
    stop("基因 '", gene, "' 不在 ", assay, " assay 的 ", slot, " 中。")
  }
  # 列检查
  md <- seu@meta.data
  if (!sample_col   %in% colnames(md)) stop("meta.data 中没找到样本列 '", sample_col, "'.")
  if (!celltype_col %in% colnames(md)) stop("meta.data 中没找到细胞类型列 '", celltype_col, "'.")
  
  # 提取表达
  expr <- Seurat::GetAssayData(seu, assay = assay, slot = slot)[gene, ]
  df <- data.frame(
    sample     = md[[sample_col]],
    celltype   = md[[celltype_col]],
    expression = as.numeric(expr),
    stringsAsFactors = FALSE
  )
  
  # 分组汇总
  stats_df <- df %>%
    group_by(sample, celltype) %>%
    summarise(
      n_total        = n(),
      n_expressing   = sum(expression > 0),
      pct_expressing = n_expressing / n_total * 100,
      avg_expression = mean(expression),
      .groups = "drop"
    )
  
  return(as.data.frame(stats_df))
}

Sox2_Plot = getGeneStatsBySampleCelltype(OSK_obj_cl,"Sox2",celltype_col  = "celltype",assay = "SCT",slot="data",sample_col    = "sample")
Klf4_Plot = getGeneStatsBySampleCelltype(OSK_obj_cl,"Klf4",celltype_col  = "celltype",assay = "SCT",slot="data",sample_col    = "sample")

###See MG cells:
Sox2_Plot_MG = Sox2_Plot[which(Sox2_Plot$celltype=="MG"),]
Klf4_Plot_MG = Klf4_Plot[which(Klf4_Plot$celltype=="MG"),]

tapply(Sox2_Plot_MG$expression,Sox2_Plot_MG$sample,summary)
tapply(Klf4_Plot_MG$expression,Klf4_Plot_MG$sample,summary)

Sox2_Plot_MG$gene = "Sox2"
Klf4_Plot_MG$gene = "Klf4"

Total_Plot = rbind(Sox2_Plot_MG,Klf4_Plot_MG)

ggplot(Sox2_Plot_MG, aes(
  x    = sample,
  y    = gene,
  size = pct_expressing,
  color= avg_expression
)) +
  geom_point() +
  scale_size_area(
    max_size = 10,
    name  = "% Expressing"
  ) +
  scale_color_gradient(
    limits = c(0.2,0.6),
    breaks = c(0.3,0.5),
    low  = "blue",
    high = "red",
    name = "Avg Expr"
  )+
  theme_classic(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title.y     = element_text(size = 12),
    legend.position  = "bottom",            # 放到下方
    legend.direction = "horizontal",        # 横向排列
    legend.box       = "horizontal",        # 如果有多个图例，横向排列它们
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(t = 10, r = 10, b = 20, l = 10)  # 给图例留空隙
  ) + xlab("") + ylab("")

  ggsave("Plot.png",width=7,height=3)



ggplot(Klf4_Plot_MG, aes(
  x    = sample,
  y    = gene,
  size = pct_expressing,
  color= avg_expression
)) +
  geom_point() +
  scale_size_area(
    max_size = 10,
    name  = "% Expressing"
  ) +
  scale_color_gradient(
    limits = c(0.01,0.04),
    breaks = c(0.01,0.04),
    low  = "blue",
    high = "red",
    name = "Avg Expr"
  )+
  theme_classic(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title.y     = element_text(size = 12),
    legend.position  = "bottom",            # 放到下方
    legend.direction = "horizontal",        # 横向排列
    legend.box       = "horizontal",        # 如果有多个图例，横向排列它们
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(t = 10, r = 10, b = 20, l = 10)  # 给图例留空隙
  ) + xlab("") + ylab("")

  ggsave("Plot.png",width=7.5,height=3)

#####
#####
#####




#####
#####
#####


head(Sox2_Plot)

getGeneExpressionDF <- function(seu, gene,
                                sample_col   = "sample",
                                celltype_col = "celltype") {
  # 检查基因是否存在
  if (!gene %in% rownames(seu@assays$SCT@counts)) {
    stop(paste0("基因 '", gene, "' 不在 SCT assay 的表达矩阵中。"))
  }
  # 检查 sample 列
  if (!sample_col %in% colnames(seu@meta.data)) {
    stop(paste0("meta.data 中没有找到样本列 '", sample_col, "'。"))
  }
  # 检查 celltype 列
  if (!celltype_col %in% colnames(seu@meta.data)) {
    stop(paste0("meta.data 中没有找到细胞类型列 '", celltype_col, "'。"))
  }
  
  # 提取表达值
  expr_vals <- Seurat::GetAssayData(seu, assay = "SCT", slot = "data")[gene, ]
  
  # 构建输出数据框
  df <- data.frame(
    cell       = names(expr_vals),
    expression = as.numeric(expr_vals),
    sample     = seu@meta.data[[sample_col]],
    celltype   = seu@meta.data[[celltype_col]],
    stringsAsFactors = FALSE
  )
  
  return(df)
}





########
########
########
########











































