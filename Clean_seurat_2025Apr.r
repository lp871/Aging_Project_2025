#########
#### clean the seurat, merge the Gene names #####
#########


#### for Mouse ####
####


ssh plyu3@10.181.57.115
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")

Mouse_meta = readRDS("Mouse_meta_2025.rds")
Mouse_seurat = readRDS("Mouse_aging_RNA_raw.rds")

####### merge gene ######
####### mat = Mat

collapse_duplicate_genes_fast <- function(mat) {
  # 输入检查：必须是 dgCMatrix
  if (!inherits(mat, "dgCMatrix")) {
    stop("输入必须是 dgCMatrix 类型的稀疏矩阵")
  }
  
  rn <- rownames(mat)
  orig_dim <- dim(mat)
  cat("合并前矩阵维度：", orig_dim[1], "行 (基因) ×", orig_dim[2], "列 (细胞)\n")
  
  # 计算每个基因出现次数
  freq <- table(rn)
  dup_genes <- names(freq[freq > 1])
  
  # 1) 保留单拷贝基因（按原始顺序）
  single_idx <- which(!rn %in% dup_genes)
  mat_singles <- mat[single_idx, , drop = FALSE]
  
  # 2) 仅对重复基因循环：按它们首次出现顺序
  first_pos <- sapply(dup_genes, function(g) match(g, rn))
  dup_genes <- dup_genes[order(first_pos)]
  ##
  collapsed_rows <- lapply(dup_genes, function(g) {
    library(Matrix)
    idxs <- which(rn == g)
    if (length(idxs) == 1) {
      # 虽然 dup_genes 中应全是 length>1，但为保险起见仍保留
      mat[idxs, , drop = FALSE]
    } else {
      # 明确调用 Matrix 包的 colSums 方法，避免降维为向量
      summed <- Matrix::colSums(mat[idxs, , drop = FALSE])
      summed
    }
  })
  
  # 合并重复基因行
  mat_dup_collapsed <- do.call(cbind, collapsed_rows)
  mat_dup_collapsed = t(mat_dup_collapsed)
  rownames(mat_dup_collapsed) <- dup_genes
  all.equal(colnames(mat_dup_collapsed),colnames(mat))
  
  # 3) 合并单拷贝行与合并后的重复行
  mat_collapsed <- rbind(mat_singles, mat_dup_collapsed)
  
  new_dim <- dim(mat_collapsed)
  cat("合并后矩阵维度：", new_dim[1], "行 (基因) ×", new_dim[2], "列 (细胞)\n")
  cat("共合并重复基因：", length(dup_genes), "个\n")
  if (length(dup_genes) > 0) {
    cat("重复基因列表：", paste(dup_genes, collapse = ", "), "\n")
  }
  
  return(mat_collapsed)
}

############
length(Mouse_seurat)
##### 32285 features ######
Mouse_seurat_merge = merge(
  x = Mouse_seurat[[1]],
  y = Mouse_seurat[-1]
)
##### clean the cells ######
print(dim(Mouse_meta))
k = which(colnames(Mouse_seurat_merge) %in% Mouse_meta$X == T)
print(length(k))
#####
Mouse_seurat_merge_cl = Mouse_seurat_merge[,k]
##### 
Mat = Mouse_seurat_merge_cl[['RNA']]@counts
rownames(Mat) = sapply(strsplit(rownames(Mat),split="~~"),function(x) x[[2]])
k = which(duplicated(rownames(Mat)) == T)
#####
Mat_rmdup = collapse_duplicate_genes_fast(Mat)
#####
#####

library(Seurat)
Mouse_seurat_merge_clcl = CreateSeuratObject(Mat_rmdup)

##### add annotations #####
m = match(colnames(Mouse_seurat_merge_clcl),Mouse_meta$X)
summary(m)
Mouse_seurat_merge_clcl$celltype = Mouse_meta$celltype
Mouse_seurat_merge_clcl$sample = Mouse_meta$sample
Mouse_seurat_merge_clcl$UMAP1 = Mouse_meta$UMAP1
Mouse_seurat_merge_clcl$UMAP2 = Mouse_meta$UMAP2

saveRDS(Mouse_seurat_merge_clcl,file="Mouse_seurat_merge_clcl_2025")

#####
#####
k1 = which(rownames(Mouse_seurat_merge_clcl) == "Zc3h11a")
k2 = which(rownames(Mat) == "Zc3h11a")
#####
k1_v = Mouse_seurat_merge_clcl[['RNA']]$counts[k1,]
k2_v = Mat[k2,]
k2_v = colSums(k2_v)
#####
all.equal(k1_v,k2_v) ### OK! ###
#####
#####
#####








#####
##### Next for Zebrafish #######
#####


setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")

Zebrafish_meta = readRDS("Zebrafish_meta_2025.rds")
Zebrafish_seurat = readRDS("Zebrafish_list_raw_counts_2025.rds")

######
######
length(Zebrafish_seurat)
##### 32520 features ######
Zebrafish_seurat_merge = merge(
  x = Zebrafish_seurat[[1]],
  y = Zebrafish_seurat[-1]
)

######
###### features #######
######

#####
##### clean the cells ######
#####

print(dim(Zebrafish_meta))
k = which(colnames(Zebrafish_seurat_merge) %in% Zebrafish_meta$cell_id == T)
print(length(k))
#####
Zebrafish_seurat_merge_cl = Zebrafish_seurat_merge[,k]
##### 
Mat = Zebrafish_seurat_merge_cl[['RNA']]@counts
rownames(Mat) = sapply(strsplit(rownames(Mat),split="~~"),function(x) x[[2]])
k = which(duplicated(rownames(Mat)) == T)
rownames(Mat)[k]
#####
Mat_rmdup = collapse_duplicate_genes_fast(Mat)

#####
#####


library(Seurat)
Zebrafish_seurat_merge_clcl = CreateSeuratObject(Mat_rmdup)

##### add annotations #####
m = match(colnames(Zebrafish_seurat_merge_clcl),Zebrafish_meta$cell_id)
summary(m)
Zebrafish_seurat_merge_clcl$celltype = Mouse_meta$celltype
Zebrafish_seurat_merge_clcl$sample = Mouse_meta$sample
Zebrafish_seurat_merge_clcl$UMAP1 = Mouse_meta$UMAP1
Zebrafish_seurat_merge_clcl$UMAP2 = Mouse_meta$UMAP2

saveRDS(Zebrafish_seurat_merge_clcl,file="Zebrafish_seurat_merge_clcl_2025")

###########
########### End #####
###########











