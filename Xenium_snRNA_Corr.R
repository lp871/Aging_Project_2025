#############
#############
##
#############
#############


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate seurat4
R


#############
#############



setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_pseudo_count <- readRDS("Mouse_pseudo_count_2025")


Mouse_pseudo_count_input = list()
for(i in 1:length(Mouse_pseudo_count)){
    tmp_mat = Mouse_pseudo_count[[i]]
    tmp_table = data.frame(sample=colnames(tmp_mat),age=colnames(tmp_mat))
    ####
    tmp_list = list(bulk_mat=tmp_mat,bulk_meta=tmp_table)
    Mouse_pseudo_count_input <- c(Mouse_pseudo_count_input,list(tmp_list))
    ####
}
names(Mouse_pseudo_count_input) = names(Mouse_pseudo_count)



Mouse_DEGs <- list()
for(i in 1:length(Mouse_pseudo_count_input)){
    print(i)
    #####
    res = Identify_Corr_pearson_total(Mouse_pseudo_count_input[[i]])
    #####
    Mouse_DEGs <- c(Mouse_DEGs,list(res))
}

names(Mouse_DEGs) <- names(Mouse_pseudo_count_input)

#####
#####
#####














##########

Identify_Corr_pearson_total <- function(list_process2) {
  bulk_mat  <- list_process2$bulk_mat
  bulk_meta <- list_process2$bulk_meta
  
  # 检查
  stopifnot(all(colnames(bulk_mat) == bulk_meta$sample))
  
  # 过滤：基因平均表达大于阈值
  bulk_mat <- bulk_mat
  
  genes <- rownames(bulk_mat)
  ages  <- as.numeric(bulk_meta$age)
  
  # 并行计算：每个基因做 Spearman 相关
  library(parallel)
  ##
  ncores <- 30
  res_list <- mclapply(
    seq_along(genes),
    function(i) {
      x <- bulk_mat[i, ]
      ct <- cor.test(x, ages, method = "pearson", exact = FALSE)
      list(rho = as.numeric(ct$estimate),
           p   = ct$p.value)
    },
    mc.cores = ncores
  )
  
  # 汇总结果
  rho    <- vapply(res_list, `[[`, numeric(1), "rho")
  pvalue <- vapply(res_list, `[[`, numeric(1), "p")
  
  res_table <- data.frame(
    Gene   = genes,
    rho    = rho,
    pvalue = pvalue,
    TAG    = "NotDEGs",
    stringsAsFactors = FALSE
  )
  
  # 只保留显著基因
  res_table$TAG <- "DEG"
  ##
  return(res_table)
}





##############
##############------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############


library(Seurat)

setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Xenium5k_Mouse_Obj <- readRDS("obj_cl2_2024")

head(Xenium5k_Mouse_Obj@meta.data)
tail(Xenium5k_Mouse_Obj@meta.data)
####
####
Xenium5k_Mouse_Obj$age = Xenium5k_Mouse_Obj$sample1

age = c("R1","R2","R3","R4")
names(age) = c("5",'12','52','112')

index = match(Xenium5k_Mouse_Obj$sample1,age)

Xenium5k_Mouse_Obj$age = names(age)[index]


Xenium5k_Mouse_Avg = getPseudoBulkByCellTime(Xenium5k_Mouse_Obj)

#### seurat_obj = Xenium5k_Mouse_Obj

getPseudoBulkByCellTime <- function(
  seurat_obj,
  assay_name   = "SCT",
  slot_name    = "counts",
  celltype_col = "celltype",
  time_col = "age"
) {
  library(Matrix)
  # 1. 提取数据
  counts <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = slot_name)
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
  ct_levels <- c("MG","RGC","AC","HC","Cone","Rod","BC","RPE")
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



#######
#######
#######


Mouse_X_pseudo_count_input = list()
for(i in 1:length(Xenium5k_Mouse_Avg)){
    tmp_mat = Xenium5k_Mouse_Avg[[i]]
    tmp_table = data.frame(sample=colnames(tmp_mat),age=colnames(tmp_mat))
    ####
    tmp_list = list(bulk_mat=tmp_mat,bulk_meta=tmp_table)
    Mouse_X_pseudo_count_input <- c(Mouse_X_pseudo_count_input,list(tmp_list))
    ####
}
names(Mouse_X_pseudo_count_input) = names(Xenium5k_Mouse_Avg)

#####
#####



Mouse_DEGs_5k <- list()
for(i in 1:length(Mouse_X_pseudo_count_input)){
    print(i)
    #####
    res = Identify_Corr_pearson_total(Mouse_X_pseudo_count_input[[i]])
    #####
    Mouse_DEGs_5k <- c(Mouse_DEGs_5k,list(res))
}

names(Mouse_DEGs_5k) <- names(Mouse_X_pseudo_count_input)


#########
names(Mouse_DEGs)
names(Mouse_DEGs_5k)

#########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ct_levels <- c("MG","RGC","AC","HC","Cone","Rod","BC","RPE")

Reslist = list()

for(i in 1:length(ct_levels)){
    tmp_ct = ct_levels[i]
    print(tmp_ct)
    #######
    snDEG = Mouse_DEGs[[which(names(Mouse_DEGs) == tmp_ct)]]
    snXen = Mouse_DEGs_5k[[which(names(Mouse_DEGs_5k) == tmp_ct)]]
    #######
    snDEG = snDEG[,c(1,2)]
    snXen = snXen[,c(1,2)]
    #######
    colnames(snDEG) = c("Gene","snDEG_corr")
    colnames(snXen) = c("Gene","snXen_corr")
    #######
    Cor_merge = merge(snDEG,snXen)
    #######
    filter = which(is.na(Cor_merge$snDEG_corr) == T | is.na(Cor_merge$snXen_corr) == T)
    Cor_merge = Cor_merge[-filter,]
    #######
    Cor_merge$CT = tmp_ct
    #######
    Reslist <- c(Reslist,list(Cor_merge))
}

Res_list_merge = do.call('rbind',Reslist)

Res_list_merge$CT = factor(Res_list_merge$CT,levels=c("MG","RGC","AC","HC","Rod","Cone","BC","RPE"))


library(ggplot2)
library(MASS)
library(viridis)  # for scale_color_viridis_c()

# 1) 用 kde2d 在一个网格上估算二维密度
dens <- kde2d(
  Res_list_merge$snDEG_corr,
  Res_list_merge$snXen_corr,
  n = 200
)

# 2) 把每个点映射到最近的网格上，取出密度值
ix <- findInterval(Res_list_merge$snDEG_corr, dens$x)
iy <- findInterval(Res_list_merge$snXen_corr, dens$y)
Res_list_merge$density <- dens$z[cbind(ix, iy)]

# 3) 绘图
p <- ggplot(Res_list_merge,
            aes(x = snDEG_corr,
                y = snXen_corr,
                color = density)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(option = "A") +
  facet_wrap(~ CT, nrow = 1) +
  theme_classic(base_size = 14) +
  theme(
    panel.border   = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x    = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x     = "snDEG_corr",
    y     = "snXen_corr",
    color = "Point\nDensity"
  )



# 显示并保存
print(p)
ggsave("Xenium5k_SnRNA_Corr.png", plot = p, width = 16, height = 3, dpi = 300)


setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
save(Res_list_merge,file="Res_list_merge_2025")

#######
#######

MG_Res_list_merge = Res_list_merge[which(Res_list_merge$CT=="MG"),]


MG_Res_list_merge$score = MG_Res_list_merge$snDEG_corr*MG_Res_list_merge$snXen_corr
MG_Res_list_merge = MG_Res_list_merge[order(MG_Res_list_merge$score,decreasing=T),]

######
######

Genes = c("Tgfbr2","Icam1","Hgf","Abi3bp","Moxd1","Robo2","Adgrg6","Sema6a","Tead2","Cxcr4")


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

