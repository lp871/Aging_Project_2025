#########
######### Aging clocks smooth and sample cells ########
#########



######### for Mouse: ##########
#########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat_clcl <- readRDS("Mouse_seurat_clcl_2025")





#########
#########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

#########

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat_clcl <- readRDS("Mouse_seurat_clcl_2025")

#########
table(Mouse_seurat_clcl$sample)
table(Mouse_seurat_clcl$celltype)

######### RPE 5 cells, Microglia 5 cells HC 5 cells ######
######### Other 25 cells ####

library(Seurat)
Mouse_seurat_clcl_list = SplitObject(Mouse_seurat_clcl, split.by = "celltype")

Mouse_seurat_knn_pseudobulk <- list()

names(Mouse_seurat_knn_pseudobulk) = names(Mouse_seurat_clcl_list)

saveRDS(Mouse_seurat_knn_pseudobulk,file="Mouse_seurat_knn_pseudobulk_2025")
###########

for(i in 1:length(Mouse_seurat_clcl_list)){
    print(i)
    print(names(Mouse_seurat_clcl_list)[i])
    #####
    tmp = Mouse_seurat_clcl_list[[i]]
    #####
    tmp_sp = SplitObject(tmp, split.by = "age")
    ##### first normalize and get PCA #######
    if(names(Mouse_seurat_clcl_list[i]) %in% c("Rod","AC","MG","BC","Cone","RGC")){
        k = 15
        npcs = 25
    }
    if(names(Mouse_seurat_clcl_list[i]) %in% c("HC","Microglia","RPE")){
        k = 2
        npcs = 5
    }
    ###############
    tmp_sp_res = sapply(tmp_sp,NormalizeAndPCA,npcs=npcs)
    ##### next smooth cells to get pseudo-bulk !!!! #######
    ####
    tmp_sp_res_pseudo = sapply(tmp_sp_res,PseudoBulkByKNN_fast,k=k,dims=1:npcs)
    #####
    #####
    tmp_sp_res_pseudo_merge = merge(tmp_sp_res_pseudo[[1]],tmp_sp_res_pseudo[-1],add.cell.ids = names(tmp_sp_res_pseudo))
    #####
    Mouse_seurat_knn_pseudobulk <- c(Mouse_seurat_knn_pseudobulk,list(tmp_sp_res_pseudo_merge))
}

######
######

NormalizeAndPCA <- function(
  seu,
  norm.method = "LogNormalize",
  scale.factor = 1e4,
  nfeatures = 2000,
  npcs = 30
) {
  # 1. 标准化
  seu <- NormalizeData(
    object = seu,
    normalization.method = norm.method,
    scale.factor = scale.factor
  )
  
  # 2. 识别高变基因
  seu <- FindVariableFeatures(
    object = seu,
    selection.method = "vst",
    nfeatures = nfeatures
  )
  
  # 3. 数据缩放（所有基因）
  seu <- ScaleData(
    object = seu,
    features = rownames(seu)
  )
  
  # 4. PCA 降维
  seu <- RunPCA(
    object = seu,
    features = VariableFeatures(seu),
    npcs = npcs
  )
  
  # 返回带有 PCA 结果的 Seurat 对象
  return(seu)
}


#######
####### seu = tmp_sp_res[[1]]

PseudoBulkByKNN_fast <- function(
  seu,
  k = 15,
  dims = 1:30,
  assay = "RNA",
  reduction = "pca"
) {
  #####
  library(Seurat)
  library(FNN)
  library(Matrix)
  # 1. PCA embedding
  emb <- Embeddings(seu, reduction)[, dims, drop = FALSE]
  
  # 2. kNN（不含自身）
  knn <- get.knn(emb, k = k)
  nn_idx <- knn$nn.index
  ncells <- nrow(nn_idx)
  
  # 3. 构建邻接稀疏矩阵 A
  #    对每个目标细胞 j，把自身 j 和它的 k 个邻居 i 都在 A[i,j]=1
  rows <- as.integer(c(t(nn_idx), seq_len(ncells)))
  cols <- rep(seq_len(ncells), each = k + 1)
  A <- sparseMatrix(i = rows, j = cols, x = 1,
                    dims = c(ncells, ncells))
  
  # 4. 原始 counts（genes × cells）
  counts_mat <- GetAssayData(seu, assay = assay, layer = "counts")

  
  # 5. 一步乘法得到伪‐bulk counts
  pseudo_counts <- counts_mat %*% A
  
  # 6. 新 Seurat 对象
  new_seu <- CreateSeuratObject(
    counts = pseudo_counts,
    meta.data = seu@meta.data,
    assay = assay
  )
  return(new_seu)
}




#########
#########------------------ Next for the Zebrafish #########--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

#########

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_seurat_clcl <- readRDS("Zebrafish_seurat_clcl_2025")

#########
table(Zebrafish_seurat_clcl$sample)
table(Zebrafish_seurat_clcl$celltype)

#########
#########


library(Seurat)
Zebrafish_seurat_clcl_list = SplitObject(Zebrafish_seurat_clcl, split.by = "celltype")

Zebrafish_seurat_knn_pseudobulk <- list()


###########
########### 所有 sample 去掉 28mo #####
########### RPE 去掉 24mo 只有一个 cell ######
###########


for(i in 1:length(Zebrafish_seurat_clcl_list)){
    print(i)
    print(names(Zebrafish_seurat_clcl_list)[i])
    #####
    tmp = Zebrafish_seurat_clcl_list[[i]]
    #####
    tmp_sp = SplitObject(tmp, split.by = "age")
    tmp_sp = tmp_sp[-which(names(tmp_sp) == "28")]
    ##### first normalize and get PCA #######
    if(names(Zebrafish_seurat_clcl_list[i]) %in% c("Rod","AC","HC","BC","Cone","RGC")){
        k = 15
        npcs = 25
    }
    if(names(Zebrafish_seurat_clcl_list[i]) %in% c("MG")){
        k = 2
        npcs = 5
    }
    if(names(Zebrafish_seurat_clcl_list[i]) %in% c("RPE","Microglia")){
        ###### small number of cells #######
        ######
        if(names(Zebrafish_seurat_clcl_list[i]) %in% c("RPE")){
            tmp_sp = tmp_sp[-which(names(tmp_sp) == "24")]
        }
        tmp_sp_res_pseudo_merge = merge(tmp_sp[[1]],tmp_sp[-1],add.cell.ids = names(tmp_sp))
        Zebrafish_seurat_knn_pseudobulk <- c(Zebrafish_seurat_knn_pseudobulk,list(tmp_sp_res_pseudo_merge))
        next
    }
    print(k)
    print(npcs)
    ########
    tmp_sp_res = sapply(tmp_sp,NormalizeAndPCA,npcs=npcs)
    ##### next smooth cells to get pseudo-bulk !!!! #######
    ####
    tmp_sp_res_pseudo = sapply(tmp_sp_res,PseudoBulkByKNN_fast,k=k,dims=1:npcs)
    #####
    #####
    tmp_sp_res_pseudo_merge = merge(tmp_sp_res_pseudo[[1]],tmp_sp_res_pseudo[-1],add.cell.ids = names(tmp_sp_res_pseudo))
    #####
    Zebrafish_seurat_knn_pseudobulk <- c(Zebrafish_seurat_knn_pseudobulk,list(tmp_sp_res_pseudo_merge))
}

names(Zebrafish_seurat_knn_pseudobulk) <- names(Zebrafish_seurat_clcl_list)

saveRDS(Zebrafish_seurat_knn_pseudobulk,file="Zebrafish_seurat_knn_pseudobulk_2025")


#######
########### process each human samples !!! Next for Human !!! ##########################
#######


#######
#######


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/share/Human_aging_new")
library(Seurat)



Human_AC_clcl2 <- readRDS("Human_AC_clcl2")
Human_HC_clcl2 <- readRDS("Human_HC_clcl2")
Human_BC_clcl2 <- readRDS("Human_BC_clcl2")
Human_Rod_clcl2 <- readRDS("Human_Rod_clcl2")


Human_Cone_clcl2 <- readRDS("Human_Cone_clcl2")
Human_RPE_clcl2 <- readRDS("Human_RPE_clcl2")
Human_Microglia_clcl2 <- readRDS("Human_Microglia_clcl2")
Human_MG_clcl2 <- readRDS("Human_MG_clcl2")
Human_RGC_clcl2 <- readRDS("Human_RGC_clcl2")


########
######## tmp = Human_HC_clcl2 ########
########
setwd("/zp1/data/share/Human_aging_new")
Human_AC_pseudo_bulk = Process_human_knn_sample_fun(Human_AC_clcl2,k=15,npcs=25)
saveRDS(Human_AC_pseudo_bulk,file="Human_AC_pseudo_bulk_2025")

Human_BC_pseudo_bulk = Process_human_knn_sample_fun(Human_BC_clcl2)
saveRDS(Human_BC_pseudo_bulk,file="Human_BC_pseudo_bulk_2025")

Human_Rod_pseudo_bulk = Process_human_knn_sample_fun(Human_Rod_clcl2,k=15,npcs=25)
saveRDS(Human_Rod_pseudo_bulk,file="Human_Rod_pseudo_bulk_2025")

Human_RPE_pseudo_bulk = Process_human_knn_sample_fun(Human_RPE_clcl2)
saveRDS(Human_RPE_pseudo_bulk,file="Human_RPE_pseudo_bulk_2025")

Human_Microglia_pseudo_bulk = Process_human_knn_sample_fun(Human_Microglia_clcl2,k=2,npcs=5)
saveRDS(Human_Microglia_pseudo_bulk,file="Human_Microglia_pseudo_bulk_2025")

Human_HC_pseudo_bulk = Process_human_knn_sample_fun(Human_HC_clcl2,k=15,npcs=25)
saveRDS(Human_HC_pseudo_bulk,file="Human_HC_pseudo_bulk_2025")

Human_Cone_pseudo_bulk = Process_human_knn_sample_fun(Human_Cone_clcl2,k=15,npcs=25)
saveRDS(Human_Cone_pseudo_bulk,file="Human_Cone_pseudo_bulk_2025")

Human_MG_pseudo_bulk = Process_human_knn_sample_fun(Human_MG_clcl2,k=15,npcs=25)
saveRDS(Human_MG_pseudo_bulk,file="Human_MG_pseudo_bulk_2025")

###########

Human_RGC_pseudo_bulk = Process_human_knn_sample_fun(Human_RGC_clcl2,k=15,npcs=25)
saveRDS(Human_RGC_pseudo_bulk,file="Human_RGC_pseudo_bulk_2025")


#############
#######
#############



NormalizeAndPCA <- function(
  seu,
  norm.method = "LogNormalize",
  scale.factor = 1e4,
  nfeatures = 2000,
  npcs = 30
) {
  # 1. 标准化
  seu <- NormalizeData(
    object = seu,
    normalization.method = norm.method,
    scale.factor = scale.factor
  )
  
  # 2. 识别高变基因
  seu <- FindVariableFeatures(
    object = seu,
    selection.method = "vst",
    nfeatures = nfeatures
  )
  
  # 3. 数据缩放（所有基因）
  seu <- ScaleData(
    object = seu,
    features = rownames(seu)
  )
  
  # 4. PCA 降维
  seu <- RunPCA(
    object = seu,
    features = VariableFeatures(seu),
    npcs = npcs
  )
  
  # 返回带有 PCA 结果的 Seurat 对象
  return(seu)
}


Process_human_knn_sample_fun <- function(tmp,k,npcs){
    ####
    ####
    tmp_sp = SplitObject(tmp, split.by = "donor")
    ####
    k = k
    npcs = npcs
    #### test  
    tmp_sp_res = sapply(tmp_sp,NormalizeAndPCA,npcs=npcs)
    ####
    tmp_sp_res_pseudo = sapply(tmp_sp_res,PseudoBulkByKNN_fast_Human,k=k,dims=1:npcs)
    ####
    tmp_sp_res_pseudo_merge = merge(tmp_sp_res_pseudo[[1]],tmp_sp_res_pseudo[-1],add.cell.ids = names(tmp_sp_res_pseudo))
    ####
    return(tmp_sp_res_pseudo_merge)
}




####
PseudoBulkByKNN_fast_Human <- function(
  seu,
  k = 15,
  dims = 1:30,
  assay = "RNA",
  reduction = "pca",
  max.cells = 750,
  seed = 42
) {
  library(Seurat)
  library(FNN)
  library(Matrix)
  # 全部细胞条码
  all_cells <- colnames(seu)
  n_all    <- length(all_cells)
  
  # 1. 在全体细胞上提取 PCA embedding
  emb_all <- Embeddings(seu, reduction)[, dims, drop = FALSE]
  
  # 2. 全体 kNN（不含自身）
  knn_all  <- get.knn(emb_all, k = k)
  nn_index <- knn_all$nn.index  # 矩阵：n_all × k
  
  # 3. 确定目标细胞：如果超出 max.cells，就随机抽样
  if (n_all > max.cells) {
    set.seed(seed)
    target_cells <- sample(all_cells, max.cells, replace = FALSE)
  } else {
    target_cells <- all_cells
  }
  # 抽样后目标数
  n_target <- length(target_cells)
  # 目标在全体中的索引
  target_idx <- match(target_cells, all_cells)
  
  # 4. 构建邻接稀疏矩阵 A (rows = 全体细胞, cols = 目标 pseudobulk 样本)
  #    对第 j 个目标（col j），取其自身 idx 和最近邻 nn_index[idx, ]
  rows <- vector("integer", length = n_target * (k + 1))
  cols <- vector("integer", length = n_target * (k + 1))
  ptr  <- 1
  for (j in seq_len(n_target)) {
    idx    <- target_idx[j]
    neighs <- nn_index[idx, ]
    group  <- c(idx, neighs)
    len    <- length(group) # k + 1
    rows[ptr:(ptr + len - 1)] <- group
    cols[ptr:(ptr + len - 1)] <- j
    ptr <- ptr + len
  }
  A <- sparseMatrix(
    i    = rows,
    j    = cols,
    x    = 1,
    dims = c(n_all, n_target)
  )
  
  # 5. 原始 counts（genes × 全体细胞）
  counts_all <- GetAssayData(seu, assay = assay, layer = "counts")
  
  # 6. 乘法得到伪‐bulk counts（genes × 目标数）
  pseudo_counts <- counts_all %*% A
  
  # 7. 构造新的 Seurat 对象，只包含这 n_target 列
  meta_sub <- seu@meta.data[target_cells, , drop = FALSE]
  new_seu   <- CreateSeuratObject(
    counts    = pseudo_counts,
    meta.data = meta_sub,
    assay     = assay,
    project   = paste0(seu@project.name, "_pseudo")
  )
  
  return(new_seu)
}
####






















