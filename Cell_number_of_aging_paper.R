###########
########### we will see how many donors we have for human samples ######
###########

###########



ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

######
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_seurat <- readRDS("Zebrafish_seurat_merge_clcl_2025")

table(Zebrafish_seurat$sample)
#####
Zebrafish_seurat_cl = Zebrafish_seurat[,which(Zebrafish_seurat$sample %in% c("28_Rep1","28_Rep2") == F)]
table(Zebrafish_seurat_cl$sample)
library(Seurat)



  
#######
####### Next for the mouse #######
#######



setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat <- readRDS("Mouse_seurat_merge_clcl_2025")
Mouse_meta <- readRDS("Mouse_meta_2025.rds")
celltype = Mouse_meta$celltype[match(colnames(Mouse_seurat),Mouse_meta$X)]
all.equal(as.character(Mouse_seurat$celltype),celltype)

setwd("/zp1/data/plyu3/Mouse_RPE_Jul12/RNA")
Mouse_RPE_merge_cl3 <- readRDS("Mouse_RPE_RNA_Final_Jul16")

cell = which(Mouse_RPE_merge_cl3$celltype == "RPE")
####### 16293 RPE + 152,796
#######
#######


setwd("/zp1/data/share/Human_aging_new")
library(Seurat)

Human_AC_clcl <- readRDS("Human_AC_clcl")
Human_HC_clcl <- readRDS("Human_HC_clcl")
Human_BC_clcl <- readRDS("Human_BC_clcl")
Human_Rod_clcl <- readRDS("Human_Rod_clcl")
Human_Cone_clcl <- readRDS("Human_Cone_clcl")
Human_RPE_clcl <- readRDS("Human_RPE_clcl")
Human_Microglia_clcl <- readRDS("Human_Microglia_clcl")
Human_MG_clcl <- readRDS("Human_MG_clcl")
Human_RGC_clcl <- readRDS("Human_RGC_clcl")

######
######
######

seurat_list <- list(
  AC       = Human_AC_clcl,
  HC       = Human_HC_clcl,
  BC       = Human_BC_clcl,
  Rod      = Human_Rod_clcl,
  Cone     = Human_Cone_clcl,
  RPE      = Human_RPE_clcl,
  Microglia= Human_Microglia_clcl,
  MG       = Human_MG_clcl,
  RGC      = Human_RGC_clcl
)

# 2. 统计每个对象的细胞数
cell_counts <- sapply(seurat_list, function(obj) {
  # 方法1：直接取列数
  ncol(obj)
  # 或者：length(Cells(obj))
  # length(Cells(obj))
})

# 3. 输出各个簇的细胞数
print(cell_counts)

# 4. 计算并输出总细胞数
total_cells <- sum(cell_counts)
cat("Total cells across all objects:", total_cells, "\n")



seurat_list <- list(
  AC       = Human_AC_clcl,
  HC       = Human_HC_clcl,
  BC       = Human_BC_clcl,
  Rod      = Human_Rod_clcl,
  Cone     = Human_Cone_clcl,
  #RPE      = Human_RPE_clcl,
  Microglia= Human_Microglia_clcl,
  MG       = Human_MG_clcl,
  RGC      = Human_RGC_clcl
)


donor_list <- lapply(seurat_list, function(obj) {
  if ("donor_id" %in% colnames(obj@meta.data)) {
    unique(obj@meta.data$donor_id)
  } else {
    NA
  }
})

# 3. 打印每个簇的 donor_id
for (name in names(donor_list)) {
  cat(paste0(name, " donors: "), paste0(donor_list[[name]], collapse = ", "), "\n")
}

# 4. 所有对象的 donor_id 并集
all_donors <- sort(unique(unlist(donor_list)))
cat("All donors across objects:", paste0(all_donors, collapse = ","), "\n")


######
###### Next we will find the Mouse RPE datasets ########
######










