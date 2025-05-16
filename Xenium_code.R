
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

##########
conda activate seurat4
R

library(Seurat)
setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")

R14_sct_merge <- readRDS(file="R14_sct_merge_2024")

library(Seurat)
setwd("/zp1/data/plyu3/Xenium5k_Mouse/Processed_data")

setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Xenium5k_Mouse_Obj <- readRDS("obj_cl2_2024")

##########
head(Xenium5k_Mouse_Obj@meta.data)

##########
table(Xenium5k_Mouse_Obj$celltype)

##########

cell_colors <- c(
  "AC"                  = "#E41A1C",
  "Adipose Cells"       = "#E7298A",
  "Astrocyte"           = "#66C2A5",
  "BC"                  = "#1F78B4",
  "Cone"                = "#FF7F00",
  "Corneal Epithelial"  = "#A65628",
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
  "Schwann"             = "#E6AB02",
  "smooth muscle cell"  = "#B2DF8A",
  "VE"                  = "#E31A1C"
)

######
######
######
seuratToUMAP <- function(seurat_obj, reduction = "umap", id_col = NULL, sample_col = NULL) {
  # Extract embeddings
  embedding <- Embeddings(seurat_obj, reduction)
  df <- as.data.frame(embedding)
  colnames(df) <- c(paste0(reduction, "_1"), paste0(reduction, "_2"))
  
  # Add cell-type annotation
  if (!is.null(id_col) && id_col %in% colnames(seurat_obj@meta.data)) {
    df$celltype <- seurat_obj@meta.data[[id_col]]
  } else {
    df$celltype <- Idents(seurat_obj)
  }
  
  # Add sample annotation if specified
  if (!is.null(sample_col) && sample_col %in% colnames(seurat_obj@meta.data)) {
    df$sample <- seurat_obj@meta.data[[sample_col]]
  }
  
  # Add cell barcodes as a column
  df$cell_id <- rownames(df)
  
  # Reorder columns: include sample if present
  cols <- c("cell_id", paste0(reduction, "_1"), paste0(reduction, "_2"), "celltype")
  if (!is.null(sample_col) && sample_col %in% colnames(seurat_obj@meta.data)) {
    cols <- c(cols, "sample")
  }
  df <- df[, cols]
  
  return(df)
}



# Then you can call:
umap_df <- seuratToUMAP(Xenium5k_Mouse_Obj,id_col='celltype',sample_col = "sample1")

###
###
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

######
###### #######
######

###### 重新写一个 function，colored by sample ######
######

umap_df$sample <- factor(
  umap_df$sample,
  levels = c("R1", "R2", "R3", "R4"),
  labels = c("5wk", "12wk", "52wk", "112wk")
)

#####
#####
sample_colors <- c(
  "5wk"   = "#E64B35",
  "12wk"  = "#4DBBD5",
  "52wk"  = "#00A087",
  "112wk" = "#3C5488"
)
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = sample)) +
  geom_point(
    shape  = 16,        # solid circle, no separate border
    size   = 0.1,       # small plot points
    alpha  = 0.25        # adjust opacity as desired
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
  theme_classic(base_size = 14) +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  ) +
  labs(
    x     = "UMAP 1",
    y     = "UMAP 2",
    color = "Sample"
  )

ggsave("test2.png",height=5,width=8)

#######
#######
#######
#######

load("R1_4_seurat_Merged_dim_sp_2024")

#######
names(R1_4_seurat_Merged_dim_sp)

####### add the cell type to R1_4_seurat_Merged_dim_sp from Xenium5k_Mouse_Obj ######
R1_4_seurat_Merged_dim_sp_list = list()

for(i in 1:length(R1_4_seurat_Merged_dim_sp)){
    #####
    tmp = R1_4_seurat_Merged_dim_sp[[i]]
    print(dim(tmp))
    m = match(tmp$cells,colnames(Xenium5k_Mouse_Obj))
    #print(summary(m))
    tmp$celltype = Xenium5k_Mouse_Obj$celltype[m]
    tmp_cl = tmp[is.na(tmp$celltype) == F,]
    print(dim(tmp_cl))
    keep_types <- c("AC", "BC", "Rod", "Cone", "RGC", "MG", "HC", "RPE")
    #####
    tmp_cl$celltype <- ifelse(tmp_cl$celltype %in% keep_types, as.character(tmp_cl$celltype), "Other")
    tmp_cl$celltype = as.factor(tmp_cl$celltype)
    #####
    R1_4_seurat_Merged_dim_sp_list = c(R1_4_seurat_Merged_dim_sp_list,list(tmp_cl))
    #####
}

names(R1_4_seurat_Merged_dim_sp_list) = names(R1_4_seurat_Merged_dim_sp)

#########
######### 5wk ##########
#########


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


Plot5wk = R1_4_seurat_Merged_dim_sp_list[[4]]

ggplot(Plot5wk, aes(x = x, y = y, color = celltype, fill = celltype)) +
  geom_point(
    shape  = 21,
    size   = 0.1,
    stroke = 0.2,
    alpha  = 0.5
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


ggsave("112wk_spatial.png",height=5,width=20)




