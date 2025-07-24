###############
###############
###############

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
library(Seurat)


######## load the mouse seurat: #####

library(Seurat)
library(Matrix)



#######
#######
#######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

###### for Zebrafish datasets ######
library(Seurat)
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_seurat <- readRDS("Zebrafish_seurat_merge_clcl_2025")
Zebrafish_meta <- readRDS("Zebrafish_meta_2025.rds")
celltype = Zebrafish_meta$celltype[match(colnames(Zebrafish_seurat),Zebrafish_meta$cell_id)]
all.equal(as.character(Zebrafish_seurat$celltype),celltype)

######
celltype_need = c("MG","RGC","AC","HC","Cone","Rod","BC","Microglia","RPE")

Zebrafish_seurat$age = sapply(strsplit(Zebrafish_seurat$sample,split="_"),function(x) x[[1]])
Zebrafish_seurat$age = as.numeric(Zebrafish_seurat$age)
table(Zebrafish_seurat$age)

###### rm the mt- gene #####
k = grep("^mt-",rownames(Zebrafish_seurat))
Zebrafish_seurat_cl = Zebrafish_seurat[-k,]

######
######

Z_ambient_genes <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/RNA_ambient/Zebrafish_ambient_RNA_Jul29")
Z_ambient_genes <- Z_ambient_genes[order(Z_ambient_genes$Avg,decreasing=T),]
Z_ambient_genes_Top50 = Z_ambient_genes$Gene[which(Z_ambient_genes$Avg > 0.005)]
Z_ambient_genes_Top50 = sapply(strsplit(Z_ambient_genes_Top50,split="~~"),function(x) x[[2]])


######
k = which(rownames(Zebrafish_seurat_cl) %in% Z_ambient_genes_Top50 == T)
Zebrafish_seurat_clcl = Zebrafish_seurat_cl[-k,]

######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### filter yang and old #######
######
table(Zebrafish_seurat_clcl$sample)
table(Zebrafish_seurat_clcl$age)

Zebrafish_YO = Zebrafish_seurat_clcl[,which(Zebrafish_seurat_clcl$age %in% c(1,3,36,48))]
table(Zebrafish_YO$age)

Zebrafish_YO_list = SplitObject(Zebrafish_YO,split.by="celltype")
names(Zebrafish_YO_list)

seurat_obj = Zebrafish_YO_list[[1]]

BC_res = downsampleCellsAndUMI(Zebrafish_YO_list[[1]])


downsampleCellsAndUMI <- function(
  seurat_obj,
  umi_target = 1000,
  min_cells   = 50,
  max_cells  = 1000,
  seed        = 123
) {
  set.seed(seed)
  ####
  counts = seurat_obj[['RNA']]$counts
  sample_id = seurat_obj$age

  # 1. 计算每个细胞的总 UMI
  umi_totals <- colSums(counts, na.rm = TRUE)
  
  # 2. 筛选 UMI > umi_target 的细胞
  keep1 <- names(umi_totals)[umi_totals > umi_target]
  if (length(keep1) == 0) stop("没有细胞的 UMI 总数超过 ", umi_target)
  
  counts_filt <- counts[, keep1, drop = FALSE]
  sample_filt <- sample_id[keep1]
  
  # 3. 计算每个 Sample 可用细胞数
  tbl <- table(sample_filt)
  n_min <- min(tbl)
  if (n_min < min_cells) {
    stop("每个 sample 至少要保留 ", min_cells, " 个细胞，目前最小只有 ", n_min)
  }
  #####
  if(n_min > max_cells){
    n_min = max_cells
  }
  #####
  # 4. 对每个 Sample 不放回抽样 n_min 个细胞
  selected <- unlist(lapply(names(tbl), function(s) {
    cells_s <- names(sample_filt)[sample_filt == s]
    sample(cells_s, size = n_min, replace = FALSE)
  }), use.names = FALSE)
  
  counts_sel <- counts_filt[, selected, drop = FALSE]
  
  # 5. 对每个细胞进行 UMI 下采样到 umi_target（multinomial）
  #    下采样后每列之和等于 umi_target
  counts_ds <- apply(counts_sel, 2, function(v) {
    # v 是一个基因 × 1 列向量
    # rmultinom 返回基因数 × 1 矩阵
    as.vector(rmultinom(1, size = umi_target, prob = v))
  })
  colnames(counts_ds) <- selected
  rownames(counts_ds) <- rownames(counts_sel)
  ####
  sample_id_new = sample_id[match(colnames(counts_ds),names(sample_id))]
  cell_metadata = data.frame(cell=names(sample_id_new),age_group=sample_id_new)
  ####
  list(
    downsampled_counts = counts_ds,  # genes × selected cells
    selected_cells     = selected,   # 被保留的细胞名
    per_sample_n       = setNames(rep(n_min, length(tbl)), names(tbl)),
    cell_metadata = cell_metadata
  )
}



######

computeNoiseByAgeOnly <- function(
  ds_out,             # downsampleCellsAndUMI() 的输出列表
                      #   ds_out$downsampled_counts: 基因 × 细胞 矩阵
                      #   ds_out$cell_metadata: data.frame，含列 cell, age_group
  n_bins = 10,
  exclude_bins = c(1, n_bins),
  pct_lowest = 0.10
) {
  library(matrixStats)
  library(dplyr)
  
  # 1. 提取下采样后矩阵和注释
  counts_ds <- ds_out$downsampled_counts
  meta_sel  <- ds_out$cell_metadata
  
  # 保证 meta 与矩阵列一致
  stopifnot(all(meta_sel$cell %in% colnames(counts_ds)))
  meta_sel <- meta_sel %>%
    arrange(match(cell, colnames(counts_ds)))
  stopifnot(all(meta_sel$cell == colnames(counts_ds)))
  
  # 2. 基因分箱与筛选
  gene_means <- rowMeans(counts_ds, na.rm = TRUE)
  bin_id     <- cut(rank(gene_means, ties.method = "first"),
                    breaks = n_bins, labels = FALSE)
  keep_bins  <- setdiff(seq_len(n_bins), exclude_bins)
  
  low_var_genes <- unlist(lapply(keep_bins, function(b) {
    idx    <- which(bin_id == b)
    submat <- counts_ds[idx, , drop = FALSE]
    cvs    <- rowSds(submat) / rowMeans(submat)
    n_sel  <- ceiling(length(cvs) * pct_lowest)
    idx[order(cvs, decreasing = FALSE)[1:n_sel]]
  }))
  low_var_genes <- unique(low_var_genes)
  
  # 3. 截取基因并 sqrt 变换
  mat_sqrt <- sqrt(counts_ds[low_var_genes, , drop = FALSE])
  
  # 4. 按 age_group 分组计算 Euclidean 距离噪声
  noise_df <- meta_sel %>%
    group_by(age_group) %>%
    do({
      cells <- .$cell
      subm  <- mat_sqrt[, cells, drop = FALSE]
      mu    <- rowMeans(subm)
      dists <- sqrt(colSums((subm - mu)^2))
      data.frame(cell = cells, noise = dists)
    }) %>%
    ungroup()
    #####
    return(data.frame(noise_df))
}

#####
##### for Zebrafish ######
#####

Zebrafish_YO_list = Zebrafish_YO_list[-6]
Zebrafish_YO_list = Zebrafish_YO_list[-8]

names(Zebrafish_YO_list)

Zebrafish_Total_res = list()
for(i in 1:length(Zebrafish_YO_list)){
    ####
    CT = names(Zebrafish_YO_list)[i]
    print(CT)
    ####
    tmp_res = downsampleCellsAndUMI(Zebrafish_YO_list[[i]])
    tmp_res_noise = computeNoiseByAgeOnly(tmp_res)
    ####
    tmp_res_noise$CT = CT
    ####
    Zebrafish_Total_res <- c(Zebrafish_Total_res,list(tmp_res_noise))
}
#####
Zebrafish_Total_res = do.call("rbind",Zebrafish_Total_res)
#####

Zebrafish_Total_res$age_group = paste0(Zebrafish_Total_res$age_group,"mo")
Zebrafish_Total_res$age_group = factor(Zebrafish_Total_res$age_group,levels=c("1mo","3mo","36mo","48mo"))

library(ggplot2)
ggplot(Zebrafish_Total_res,aes(x=age_group,y=noise)) + geom_boxplot() + theme_classic() +
  facet_wrap(~ CT, nrow = 2) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + xlab("age group") + ylab("transcriptional noise")
ggsave("Zebrafish_res.png",height=4,width=7)




#####
##### transriptional noise ######
#####



############### Next for Mouse #############



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



Mouse_YO = Mouse_seurat_clcl[,which(Mouse_seurat_clcl$age %in% c(5,12,108,120))]
table(Mouse_YO$age)

Mouse_YO_list = SplitObject(Mouse_YO,split.by="celltype")
names(Mouse_YO_list)

Mouse_YO_list = Mouse_YO_list[-8]
Mouse_YO_list = Mouse_YO_list[-8]

Mouse_Total_res = list()
for(i in 1:length(Mouse_YO_list)){
    ####
    CT = names(Mouse_YO_list)[i]
    print(CT)
    ####
    tmp_res = downsampleCellsAndUMI(Mouse_YO_list[[i]])
    tmp_res_noise = computeNoiseByAgeOnly(tmp_res)
    ####
    tmp_res_noise$CT = CT
    ####
    Mouse_Total_res <- c(Mouse_Total_res,list(tmp_res_noise))
}
#####
Mouse_Total_res = do.call("rbind",Mouse_Total_res)
#####
#####

Mouse_Total_res$age_group = paste0(Mouse_Total_res$age_group,"wk")
Mouse_Total_res$age_group = factor(Mouse_Total_res$age_group,levels=c("5wk","12wk","108wk","120wk"))

library(ggplot2)
ggplot(Mouse_Total_res,aes(x=age_group,y=noise)) + geom_boxplot() + theme_classic() +
  facet_wrap(~ CT, nrow = 2) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x    = element_text(angle = 60, hjust = 1)) + xlab("age group") + ylab("transcriptional noise")
ggsave("Mouse_res.png",height=4,width=7)



######
###### Next we will load the Human samples ########-----------------------------------------------------------------------------------------------------------------------------------------------
######

setwd("/zp1/data/share/Human_aging_new")
library(Seurat)

Human_RPE_clcl2 <- readRDS(file="Human_RPE_clcl2")
Human_Rod_clcl2 <- readRDS(file="Human_Rod_clcl2")
Human_Cone_clcl2 <- readRDS(file="Human_Cone_clcl2")
Human_MG_clcl2 <- readRDS(file="Human_MG_clcl2")
Human_AC_clcl2 <- readRDS(file="Human_AC_clcl2")
Human_BC_clcl2 <- readRDS(file="Human_BC_clcl2")
Human_HC_clcl2 <- readRDS(file="Human_HC_clcl2")
Human_Microglia_clcl2 <- readRDS(file="Human_Microglia_clcl2")
Human_RGC_clcl2 <- readRDS(file="Human_RGC_clcl2")

###########
########### subset the intervals #####
###########

table(Human_RPE_clcl2_cl$age2)
Human_RPE_clcl2$age = paste0(Human_RPE_clcl2$age2,"_",Human_RPE_clcl2$donor)
Human_RPE_clcl2_cl = Human_RPE_clcl2[,which(Human_RPE_clcl2$age2 %in% c(10,25,80,90) == T)]

Human_Rod_clcl2$age = paste0(Human_Rod_clcl2$age2,"_",Human_Rod_clcl2$donor)
Human_Rod_clcl2_cl = Human_Rod_clcl2[,which(Human_Rod_clcl2$age2 %in% c(10,25,80,90) == T)]

Human_Cone_clcl2$age = paste0(Human_Cone_clcl2$age2,"_",Human_Cone_clcl2$donor)
Human_Cone_clcl2_cl = Human_Cone_clcl2[,which(Human_Cone_clcl2$age2 %in% c(10,25,80,90) == T)]

Human_AC_clcl2$age = paste0(Human_AC_clcl2$age2,"_",Human_AC_clcl2$donor)
Human_AC_clcl2_cl = Human_AC_clcl2[,which(Human_AC_clcl2$age2 %in% c(10,25,80,90) == T)]

Human_HC_clcl2$age = paste0(Human_HC_clcl2$age2,"_",Human_HC_clcl2$donor)
Human_HC_clcl2_cl = Human_HC_clcl2[,which(Human_HC_clcl2$age2 %in% c(10,25,80,90) == T)]

Human_RGC_clcl2$age = paste0(Human_RGC_clcl2$age2,"_",Human_RGC_clcl2$donor)
Human_RGC_clcl2_cl = Human_RGC_clcl2[,which(Human_RGC_clcl2$age2 %in% c(10,25,80,90) == T)]

Human_MG_clcl2$age = paste0(Human_MG_clcl2$age2,"_",Human_MG_clcl2$donor)
Human_MG_clcl2_cl = Human_MG_clcl2[,which(Human_MG_clcl2$age2 %in% c(10,25,80,90) == T)]

Human_BC_clcl2$age = paste0(Human_BC_clcl2$age2,"_",Human_BC_clcl2$donor)
Human_BC_clcl2_cl = Human_BC_clcl2[,which(Human_BC_clcl2$age2 %in% c(10,25,80,90) == T)]

########
########

Human_List = list(RPE=Human_RPE_clcl2_cl,AC=Human_AC_clcl2_cl,HC=Human_HC_clcl2_cl,BC=Human_BC_clcl2_cl,Rod=Human_Rod_clcl2_cl,Cone=Human_Cone_clcl2_cl,MG=Human_MG_clcl2_cl,RGC=Human_RGC_clcl2_cl)

Human_Total_res = list()
for(i in 1:length(Human_List)){
    ####
    CT = names(Human_List)[i]
    print(CT)
    ####
    tmp_res = downsampleCellsAndUMI(Human_List[[i]])
    tmp_res_noise = computeNoiseByAgeOnly(tmp_res)
    ####
    tmp_res_noise$CT = CT
    ####
    Human_Total_res <- c(Human_Total_res,list(tmp_res_noise))
}
#####
Human_Total_res = do.call("rbind",Human_Total_res)

#######
Human_Total_res$age_group2 = sapply(strsplit(Human_Total_res$age_group,split="_"),function(x) x[[1]])

#######

library(ggplot2)
ggplot(Human_Total_res,aes(x=age_group2,y=noise)) + geom_boxplot() + theme_classic() +
  facet_wrap(~ CT, nrow = 2) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x    = element_text(angle = 60, hjust = 1)) + xlab("age group") + ylab("transcriptional noise")
ggsave("Human_res.png",height=4,width=7)




ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R
########
########

