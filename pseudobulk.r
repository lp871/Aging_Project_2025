#######
pseudo-bulk for each time point #####
#######
library(Seurat)
library(Matrix)



#######
#######
#######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

######
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_seurat <- readRDS("Zebrafish_seurat_merge_clcl_2025")

Zebrafish_meta <- readRDS("Zebrafish_meta_2025.rds")

celltype = Zebrafish_meta$celltype[match(colnames(Zebrafish_seurat),Zebrafish_meta$cell_id)]
all.equal(as.character(Zebrafish_seurat$celltype),celltype)

######

celltype_need = c("MG","RGC","AC","HC","Cone","Rod","BC","Microglia","RPE")

###### add time points #####
######
Zebrafish_seurat$age = sapply(strsplit(Zebrafish_seurat$sample,split="_"),function(x) x[[1]])
Zebrafish_seurat$age = as.numeric(Zebrafish_seurat$age)
table(Zebrafish_seurat$age)

###### merge 106 to 108 #######
###### seurat_obj = Mouse_seurat ###
######

###### rm the mt- gene #####
k = grep("^mt-",rownames(Zebrafish_seurat))
Zebrafish_seurat_cl = Zebrafish_seurat[-k,]

###### Mouse_MG_test = readRDS("/zp1/data/plyu3/Aging_Clocks_Final_2025/Mouse/Mouse_MG_seurat2025.rds")
######

Z_ambient_genes <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/RNA_ambient/Zebrafish_ambient_RNA_Jul29")
Z_ambient_genes <- Z_ambient_genes[order(Z_ambient_genes$Avg,decreasing=T),]
Z_ambient_genes_Top50 = Z_ambient_genes$Gene[which(Z_ambient_genes$Avg > 0.005)]
Z_ambient_genes_Top50 = sapply(strsplit(Z_ambient_genes_Top50,split="~~"),function(x) x[[2]])


######
k = which(rownames(Zebrafish_seurat_cl) %in% Z_ambient_genes_Top50 == T)
Zebrafish_seurat_clcl = Zebrafish_seurat_cl[-k,]

###### seurat_obj = Mouse_seurat_clcl

getPseudoBulkByCellTime <- function(
  seurat_obj,
  assay_name   = "RNA",
  slot_name    = "counts",
  celltype_col = "celltype",
  time_col = "age"
) {
  library(Matrix)
  # 1. 提取数据
  counts <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = slot_name)
  meta   <- seurat_obj@meta.data
  # 2. 校验
  if (!all(c(celltype_col, time_col) %in% colnames(meta))) {
    stop("在 meta.data 中未找到列：",
         paste(setdiff(c(celltype_col, time_col), colnames(meta)), collapse = ", "))
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  # 3. 获取所有 celltype
  ct_levels <- c("MG","RGC","AC","HC","Cone","Rod","BC","Microglia","RPE")
  # 4. 对每个 celltype，按 time 汇总
  pb_list <- lapply(ct_levels, function(ct) {
    # 4.1 过滤出当前 celltype 的细胞
    cells_ct <- rownames(meta)[meta[[celltype_col]] == ct]
    if (length(cells_ct) == 0) {
      # 如果没有细胞，返回 NULL
      return(NULL)
    }
    sub_counts <- counts[, cells_ct, drop = FALSE]
    sub_meta   <- meta[cells_ct, , drop = FALSE]
    
    # 4.2 构造 time 分组因子
    grp_time <- factor(sub_meta[[time_col]], 
                       levels = sort(unique(sub_meta[[time_col]])))
    
    # 4.3 设计矩阵（细胞 × time）
    design <- model.matrix(~0 + grp_time)
    colnames(design) <- levels(grp_time)
    
    # 4.4 稀疏矩阵乘法汇总
    Matrix(sub_counts %*% Matrix(design, sparse = TRUE), sparse = TRUE)
  })
  
  # 5. 设置 list 名称，并剔除可能的 NULL
  names(pb_list) <- ct_levels
  ####
  ####
  ####
  pb_list_norm = normalizePBList(pb_list)
  ####
  return(pb_list_norm)
}

#######
#######
#######
normalizePBList <- function(pb_list, scale = 1e4, pseudocount = 1) {
  lapply(pb_list, function(mat) {
    # 1. 计算每个样本（列）的总计数
    lib_sizes <- Matrix::colSums(mat)
    # 2. 按列除以 library size 并乘以 scale → CPM
    #    注意使用 t() 以保持稀疏矩阵结构，转换后可能变为普通矩阵
    cpm <- t( t(mat) / lib_sizes * scale )
    # 3. log2 转换
    log2(cpm + pseudocount)
  }) -> normalized_list
  # 保留原来的 names
  names(normalized_list) <- names(pb_list)
  normalized_list
}

Zebrafish_pseudo_count = getPseudoBulkByCellTime(Zebrafish_seurat_clcl)
names(Zebrafish_pseudo_count)

scale(Zebrafish_pseudo_count[['MG']]["sparc",])
Zebrafish_pseudo_count[['MG']]["b2m",]

saveRDS(Zebrafish_pseudo_count,file="Zebrafish_pseudo_count_2025")


##
##### Next Mouse ##########---------#########
##

##
######
##

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat <- readRDS("Mouse_seurat_merge_clcl_2025")
Mouse_meta <- readRDS("Mouse_meta_2025.rds")
celltype = Mouse_meta$celltype[match(colnames(Mouse_seurat),Mouse_meta$X)]
all.equal(as.character(Mouse_seurat$celltype),celltype)

######

celltype_need = c("MG","RGC","AC","HC","Cone","Rod","BC","Microglia","RPE")

###### add time points #####
######
Mouse_seurat$age = sapply(strsplit(Mouse_seurat$sample,split="_"),function(x) x[[1]])
Mouse_seurat$age = as.numeric(gsub("Wk","",Mouse_seurat$age))
table(Mouse_seurat$age)

###### merge 106 to 108 #######
###### seurat_obj = Mouse_seurat ###
######
Mouse_seurat$age[which(Mouse_seurat$age == 106)] = 108
table(Mouse_seurat$age)

###### rm the mt- gene #####
k = grep("^mt-",rownames(Mouse_seurat))
Mouse_seurat_cl = Mouse_seurat[-k,]

###### Mouse_MG_test = readRDS("/zp1/data/plyu3/Aging_Clocks_Final_2025/Mouse/Mouse_MG_seurat2025.rds")
######

M_ambient_genes <- readRDS("/zp1/data/plyu3/Aging_Clocks_Final/Mouse_ambient_RNA_Jul22")
M_ambient_genes <- M_ambient_genes[order(M_ambient_genes$Avg,decreasing=T),]
M_ambient_genes_Top50 = M_ambient_genes$Gene[which(M_ambient_genes$Avg > 0.005)]
M_ambient_genes_Top50 = sapply(strsplit(M_ambient_genes_Top50,split="~~"),function(x) x[[2]])


######
k = which(rownames(Mouse_seurat_cl) %in% M_ambient_genes_Top50 == T)
Mouse_seurat_clcl = Mouse_seurat_cl[-k,]

######
Mouse_pseudo_count = getPseudoBulkByCellTime(Mouse_seurat_clcl)
names(Mouse_pseudo_count)

scale(Mouse_pseudo_count[['MG']]["Apoe",])
saveRDS(Mouse_pseudo_count,file="Mouse_pseudo_count_2025")


######


