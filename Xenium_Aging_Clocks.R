################

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate seurat4
R


################
NormalizeAndPCASCT <- function(
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
#######



setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Xenium5k_Mouse_Obj <- readRDS("obj_cl2_2024")


####### split by cell type #####
#######

library(Seurat)

table(Xenium5k_Mouse_Obj$celltype)
table(Xenium5k_Mouse_Obj$sample1)

####### 5 12 52 112 #######
library(dplyr)

Xenium5k_Mouse_Obj@meta.data <- Xenium5k_Mouse_Obj@meta.data %>%
  mutate(
    sample1 = recode(
      sample1,
      "R1" = "5",
      "R2" = "12",
      "R3" = "52",
      "R4" = "112"
    )
  )

table(Xenium5k_Mouse_Obj$sample1)
Xenium5k_Mouse_Obj$sample = Xenium5k_Mouse_Obj$sample1

#######

Xenium5k_Mouse_Obj_list = SplitObject(Xenium5k_Mouse_Obj, split.by = "celltype")
Xenium5k_Mouse_Obj_list_cl = Xenium5k_Mouse_Obj_list[c("Rod","Cone","BC","AC","RGC","MG","RPE","HC")]


#######
Xenium5k_Mouse_Obj_pseudobulk_list <- list()

#######
for(i in 1:length(Xenium5k_Mouse_Obj_list_cl)){
    print(i)
    print(names(Xenium5k_Mouse_Obj_list_cl)[i])
    #####
    tmp = Xenium5k_Mouse_Obj_list_cl[[i]]
    #####
    tmp_sp = SplitObject(tmp, split.by = "sample")
    ##### first normalize and get PCA #######
    if(names(Xenium5k_Mouse_Obj_list_cl[i]) %in% c("Rod","AC","MG","BC","Cone","RGC")){
        k = 15
        npcs = 25
    }
    if(names(Xenium5k_Mouse_Obj_list_cl[i]) %in% c("HC","RPE")){
        k = 2
        npcs = 5
    }
    ###############
    tmp_sp_res = sapply(tmp_sp,NormalizeAndPCASCT,npcs=npcs)
    ##### next smooth cells to get pseudo-bulk !!!! #######
    ####
    tmp_sp_res_pseudo = sapply(tmp_sp_res,PseudoBulkByKNN_fast_SCT,k=k,dims=1:npcs)
    #####
    #####
    tmp_sp_res_pseudo_merge = merge(tmp_sp_res_pseudo[[1]],tmp_sp_res_pseudo[-1],add.cell.ids = names(tmp_sp_res_pseudo))
    #####
    Xenium5k_Mouse_Obj_pseudobulk_list <- c(Xenium5k_Mouse_Obj_pseudobulk_list,list(tmp_sp_res_pseudo_merge))
}

#####
getwd()

#####
setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
names(Xenium5k_Mouse_Obj_pseudobulk_list) <- names(Xenium5k_Mouse_Obj_list_cl)
saveRDS(Xenium5k_Mouse_Obj_pseudobulk_list,file="Xenium5k_Mouse_Obj_pseudobulk_list_May2025")

#####
##### seu = Xenium5k_Mouse_Obj_pseudobulk_list$HC ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#####

Mouse_MG_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$MG)
saveRDS(Mouse_MG_model,file="Mouse_MG_Xmodel.rds")

Mouse_Rod_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$Rod)
saveRDS(Mouse_Rod_model,file="Mouse_Rod_Xmodel.rds")

Mouse_Cone_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$Cone)
saveRDS(Mouse_Cone_model,file="Mouse_Cone_Xmodel.rds")

Mouse_AC_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$AC)
saveRDS(Mouse_AC_model,file="Mouse_AC_Xmodel.rds")

Mouse_HC_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$HC)
saveRDS(Mouse_HC_model,file="Mouse_HC_Xmodel.rds")

Mouse_BC_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$BC)
saveRDS(Mouse_BC_model,file="Mouse_BC_Xmodel.rds")

Mouse_RGC_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$RGC)
saveRDS(Mouse_RGC_model,file="Mouse_RGC_Xmodel.rds")

Mouse_RPE_model = Elastic_net_one_round_Xenium(seu=Xenium5k_Mouse_Obj_pseudobulk_list$RPE)
saveRDS(Mouse_RPE_model,file="Mouse_RPE_Xmodel.rds")

#########
#########--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########

Elastic_net_one_round_Xenium <- function(
  seu,
  scale.factor = 1e4,
  train_frac = 0.7,
  alpha.vec = seq(0, 1, by = 0.1),
  nfolds = 5,
  ncores = 20,
  seed = 42
) {
  library(Seurat)
  raw_layers <- seu[["RNA"]]@layers
  mat <- do.call(cbind, raw_layers)
  rownames(mat) = rownames(seu)
  colnames(mat) = colnames(seu)
  # 1. 归一化 (CPM 标准化后乘以 scale.factor) 并 log2(x+1)
  total_counts <- colSums(mat)
  mat_norm <- sweep(mat, 2, total_counts, FUN = "/") * scale.factor
  mat_norm <- log2(mat_norm + 1)
  # 2. 筛选目标基因
  mat_cl <- mat_norm
  # 3. 分层抽样：按 age 层内随机拆分训练/测试集
  set.seed(seed)
  ###
  seu$age = seu$sample
  ages <- seu$age
  samples <- colnames(mat_cl)
  train_idx <- c()
  for(a in unique(ages)){
    idx_age <- which(ages == a)
    n_age <- length(idx_age)
    n_train <- floor(n_age * train_frac)
    train_idx <- c(train_idx, sample(idx_age, n_train))
  }
  test_idx <- setdiff(seq_along(samples), train_idx)
  ####
  train_matrix <- mat_cl[, train_idx, drop = FALSE]
  test_matrix  <- mat_cl[, test_idx, drop = FALSE]
  ####
  x_train <- t(train_matrix)
  y_train <- as.numeric(ages[train_idx])
  x_test  <- t(test_matrix)
  y_test  <- as.numeric(ages[test_idx])

  # 4. 计算训练集权重 (inverse frequency) 并标准化
  age_tab <- table(y_train)
  w_train <- 1 / age_tab[as.character(y_train)]
  w_train <- w_train * length(w_train) / sum(w_train)

  # 5. 并行搜索最佳 alpha 并训练模型
  library(parallel)
  library(glmnet)
  ###
  cv_results <- mclapply(alpha.vec, function(a) {
    cv <- glmnet::cv.glmnet(
      x = x_train, y = y_train,
      alpha = a,
      weights = w_train,
      type.measure = "mse",
      nfolds = nfolds,
      standardize = FALSE
    )
    list(alpha = a,
         lambda.min = cv$lambda.min,
         cvm = min(cv$cvm),
         fit = cv)
  }, mc.cores = ncores)
  ###
  best <- Reduce(function(x, y) if (x$cvm < y$cvm) x else y, cv_results)
  best_alpha <- best$alpha
  fit <- best$fit

  # 6. 测试集预测及性能评估
  preds_test <- as.numeric(predict(fit, newx = x_test, s = "lambda.min"))
  test_mse <- mean((y_test - preds_test)^2)
  test_cor <- cor(y_test, preds_test)
  print(paste("test_cor:",test_cor))
  # 7. 提取系数
  coefs <- as.matrix(coef(fit, s = "lambda.min"))
  coef_df <- data.frame(
    gene = rownames(coefs),
    coefficient = as.vector(coefs)
  )[-1, , drop = FALSE]

  # 8. 构建预测结果表
  pred_df <- data.frame(
    sample    = colnames(test_matrix),
    actual    = y_test,
    predicted = preds_test,
    set       = "test"
  )

  # 9. 绘图并保存
  library(ggplot2)
  p <- ggplot(pred_df, aes(x = actual, y = predicted)) +
    geom_jitter(size = 0.5) +
    labs(
      title = sprintf("Test Set: Actual vs Predicted Age (alpha=%.1f)", best_alpha),
      x = "Actual Age", y = "Predicted Age"
    )
  ggsave("elastic_net_performance.png", plot = p, height = 4, width = 6)

  return(list(
    model        = fit,
    best_alpha   = best_alpha,
    coefficients = coef_df,
    performance  = list(test_mse = test_mse, test_cor = test_cor),
    predictions  = pred_df,
    train_matrix = train_matrix,
    test_matrix  = test_matrix
  ))
}







#######
#######
####### seu = tmp_sp_res[[1]]
#######
#######


PseudoBulkByKNN_fast_SCT <- function(
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
  counts_mat <- GetAssayData(seu, assay = "SCT", layer = "counts")

  
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





###############
############### Aging clocks for Xenium ###########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
############### 

setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Mouse_MG_model = readRDS("Mouse_MG_Xmodel.rds")
Mouse_AC_model = readRDS("Mouse_AC_Xmodel.rds")
Mouse_HC_model = readRDS("Mouse_HC_Xmodel.rds")
Mouse_RGC_model = readRDS("Mouse_RGC_Xmodel.rds")
Mouse_BC_model = readRDS("Mouse_BC_Xmodel.rds")
Mouse_Rod_model = readRDS("Mouse_Rod_Xmodel.rds")
Mouse_Cone_model = readRDS("Mouse_Cone_Xmodel.rds")
Mouse_RPE_model = readRDS("Mouse_RPE_Xmodel.rds")

###############----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Mouse_MG_Aging_Clock_Res <- Mouse_MG_model$predictions
Mouse_MG_Aging_Clock_Res$class <- "MG"
cor(Mouse_MG_Aging_Clock_Res$predicted, Mouse_MG_Aging_Clock_Res$actual)

Mouse_Rod_Aging_Clock_Res <- Mouse_Rod_model$predictions
Mouse_Rod_Aging_Clock_Res$class <- "Rod"
cor(Mouse_Rod_Aging_Clock_Res$predicted, Mouse_Rod_Aging_Clock_Res$actual)

Mouse_Cone_Aging_Clock_Res <- Mouse_Cone_model$predictions
Mouse_Cone_Aging_Clock_Res$class <- "Cone"
cor(Mouse_Cone_Aging_Clock_Res$predicted, Mouse_Cone_Aging_Clock_Res$actual)

Mouse_RGC_Aging_Clock_Res <- Mouse_RGC_model$predictions
Mouse_RGC_Aging_Clock_Res$class <- "RGC"
cor(Mouse_RGC_Aging_Clock_Res$predicted, Mouse_RGC_Aging_Clock_Res$actual)

Mouse_BC_Aging_Clock_Res <- Mouse_BC_model$predictions
Mouse_BC_Aging_Clock_Res$class <- "BC"
cor(Mouse_BC_Aging_Clock_Res$predicted, Mouse_BC_Aging_Clock_Res$actual)

Mouse_AC_Aging_Clock_Res <- Mouse_AC_model$predictions
Mouse_AC_Aging_Clock_Res$class <- "AC"
cor(Mouse_AC_Aging_Clock_Res$predicted, Mouse_AC_Aging_Clock_Res$actual)

Mouse_HC_Aging_Clock_Res <- Mouse_HC_model$predictions
Mouse_HC_Aging_Clock_Res$class <- "HC"
cor(Mouse_HC_Aging_Clock_Res$predicted, Mouse_HC_Aging_Clock_Res$actual)

Mouse_RPE_Aging_Clock_Res <- Mouse_RPE_model$predictions
Mouse_RPE_Aging_Clock_Res$class <- "RPE"
cor(Mouse_RPE_Aging_Clock_Res$predicted, Mouse_RPE_Aging_Clock_Res$actual)

# 合并所有 Mouse 结果
All_merge_res <- rbind(
  Mouse_MG_Aging_Clock_Res,
  Mouse_Cone_Aging_Clock_Res,
  Mouse_Rod_Aging_Clock_Res,
  Mouse_RGC_Aging_Clock_Res,
  Mouse_BC_Aging_Clock_Res,
  Mouse_AC_Aging_Clock_Res,
  Mouse_HC_Aging_Clock_Res,
  Mouse_RPE_Aging_Clock_Res
)

# 查看合并后的数据
head(All_merge_res)

# 创建一个新的 Samples 列用于绘图
All_merge_res$Samples <- All_merge_res$actual
table(All_merge_res$Samples)

# 按照时间点顺序设置因子水平
All_merge_res$Samples <- factor(
  All_merge_res$Samples,
  levels = c(5,12,52,112)
)

# 按照细胞类型设置因子水平
All_merge_res$class <- factor(
  All_merge_res$class,
  levels = c("MG", "RGC", "AC", "HC", "Rod", "Cone", "BC", "RPE")
)

# 加载绘图包
library(ggplot2)
library(viridis)

# 绘制分面小提琴＋箱线图
p <- ggplot(All_merge_res, aes(x = Samples, y = predicted, color = actual, fill = actual)) +
  facet_wrap(~ class, nrow = 1) +
  geom_violin() +
  geom_boxplot(
    width = 0.25,
    fill = "grey",
    color = "black",
    outlier.shape = NA
  ) +
  scale_color_viridis() +
  scale_fill_viridis() +
  theme_classic() +
  scale_y_continuous(limits = c(-25, 150)) +
  theme(
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border      = element_rect(color = "black", fill = NA, size = 1)
  ) +
  xlab("") +
  ylab("")

ggsave("AllMouse_results.png",height=3,width=20) 














###############
############### Aging clocks predict for Xenium OSK ###########
###############


############### load the OSK samples ##########################-----------------------------------------------------------------------------
########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate seurat4
R


setwd("/zp1/data/plyu3/Xenium5k_Mouse/OSK_New")


######## load the aging clocks !! ######


setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Mouse_MG_model = readRDS("Mouse_MG_Xmodel.rds")
Mouse_AC_model = readRDS("Mouse_AC_Xmodel.rds")
Mouse_HC_model = readRDS("Mouse_HC_Xmodel.rds")
Mouse_RGC_model = readRDS("Mouse_RGC_Xmodel.rds")
Mouse_BC_model = readRDS("Mouse_BC_Xmodel.rds")
Mouse_Rod_model = readRDS("Mouse_Rod_Xmodel.rds")
Mouse_Cone_model = readRDS("Mouse_Cone_Xmodel.rds")
Mouse_RPE_model = readRDS("Mouse_RPE_Xmodel.rds")

##########
########## OSK_obj_pseudobulk_list
##########


names(OSK_obj_pseudobulk_list)

##########


Prepare_Matrix_And_Meta_OSK <- function(OSK_obj_pseudobulk_list){
  #########
  Results <- c()
  #########
  for(i in 1:length(OSK_obj_pseudobulk_list)){
    ########
    tmp = OSK_obj_pseudobulk_list[[i]]
    ########
    tmp[["RNA"]] <- JoinLayers(tmp[["RNA"]])
    tmp_matrix <- GetAssayData(tmp, assay = "RNA", layer = "counts")
    ########
    tmp_meta = tmp@meta.data
    tmp_meta = tmp_meta[,c("sample","celltype")]
    tmp_meta$cellid = rownames(tmp_meta)
    tmp_meta$age = 0
    ########
    ######## norm the matrix #########
    ########
    tmp_matrix_sum = Matrix::colSums(tmp_matrix)
    tmp_matrix_sum_factor = tmp_matrix_sum / 1e4
    ########
    tmp_matrix_norm = sweep(tmp_matrix,2,tmp_matrix_sum_factor,FUN= "/")
    tmp_matrix_sum2 = Matrix::colSums(tmp_matrix_norm)
    ######## change the log2 ##############
    ########
    tmp_matrix_norm_Log2 = log2(tmp_matrix_norm+1)
    ########
    reslist = list(Matrix = tmp_matrix_norm_Log2, Meta = tmp_meta)
    ########
    Results <- c(Results,list(reslist))
  }
  #########
  names(Results) <- names(OSK_obj_pseudobulk_list)
  #########
  return(Results)
}

######
######
######

setwd("/zp1/data/plyu3/Xenium5k_Mouse/OSK_New")
names(OSK_obj_pseudobulk_list) <- names(OSK_obj_cl_list_cl)
saveRDS(OSK_obj_pseudobulk_list,file="OSK_obj_pseudobulk_list_May2025")

########### load the OSK ######
setwd("/zp1/data/plyu3/Xenium5k_Mouse/OSK_New")
OSK_obj_pseudobulk_list = readRDS("OSK_obj_pseudobulk_list_May2025")
OSK_obj_pseudobulk_list_prepared = Prepare_Matrix_And_Meta_OSK(OSK_obj_pseudobulk_list)

###########
####
#### tapply(MG_predict_res$preds,MG_predict_res$sample,summary)#####
####
###########

# define the cell‐types you want
cell_types <- c("MG", "Rod", "Cone", "RGC", "AC", "HC", "RPE", "BC")

# run PredictAgeFromModel for each and store in a named list
predict_results <- lapply(cell_types, function(ct) {
  # fetch the pseudobulk matrix & meta by name
  mat  <- OSK_obj_pseudobulk_list_prepared[[ct]]$Matrix
  meta <- OSK_obj_pseudobulk_list_prepared[[ct]]$Meta
  
  # fetch the corresponding model object, e.g. Mouse_Rod_model$model
  model <- get(paste0("Mouse_", ct, "_model"))$model
  
  # run prediction and tag the class
  res <- PredictAgeFromModel(mat, meta, model)
  res$class <- ct
  return(res)
})

# name each element of the list by its cell‐type
names(predict_results) <- cell_types

# if you really want each one as its own variable:
list2env(predict_results, envir = .GlobalEnv)



#########
#########
#########
#########
all_res = do.call("rbind",predict_results)

all_res = all_res[which(all_res$sample %in% c("OSKR3","OSKR4") == T),]
all_res$class = factor(all_res$class,levels=c("MG","RGC","AC","HC","Rod","Cone","BC","RPE"))
all_res$age = factor(all_res$sample,levels=c("OSKR3","OSKR4"))

########
setwd("/zp1/data/share/OSK")
library(ggplot2)
ggplot(all_res,aes(x=age ,y=preds,color=age)) + facet_wrap(~ class, nrow = 2,scales = "free_y") + scale_y_continuous(expand=c(0,0)) + geom_boxplot(outlier.shape = NA,size=1)+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("Mouse_Xenium_OSK_2025.png",height=3,width=6) 

saveRDS(all_res,file="OSK_predict_by_Xeniun_clock.rds")

####### Next we will use aging clocks from snRNA-seq to predict #####


# define the cell‐types you want
cell_types <- c("MG", "Rod", "Cone", "RGC", "AC", "HC", "RPE", "BC")

# run PredictAgeFromModel for each and store in a named list
predict_results <- lapply(cell_types, function(ct) {
  # fetch the pseudobulk matrix & meta by name
  mat  <- OSK_obj_pseudobulk_list_prepared[[ct]]$Matrix
  meta <- OSK_obj_pseudobulk_list_prepared[[ct]]$Meta
  
  # fetch the corresponding model object, e.g. Mouse_Rod_model$model
  model <- get(paste0("Mouse_", ct, "_model_Xen"))$model
  
  # run prediction and tag the class
  res <- PredictAgeFromModel(mat, meta, model)
  res$class <- ct
  return(res)
})


all_res = do.call("rbind",predict_results)
all_res$class = factor(all_res$class,levels=c("MG","RGC","AC","HC","Rod","Cone","BC","RPE"))
all_res$age = factor(all_res$sample,levels=c("OSKR1","OSKR2","OSKR3","OSKR4"))
########
setwd("/zp1/data/share/OSK")
library(ggplot2)
ggplot(all_res,aes(x=age ,y=preds,color=age)) + facet_wrap(~ class, nrow = 1,scales = "free_y") + scale_y_continuous(expand=c(0,0)) + geom_boxplot(outlier.shape = NA,size=1)+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("Mouse_Xenium_OSK_snRNAseq.png",height=2,width=12) 

saveRDS(all_res,file="OSK_predict_by_snRNAseq_clock.rds")


####### see the pvalue #######

t.test(all_res$preds[which(all_res$celltype == "MG" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "MG" & all_res$sample == "OSKR4")],alternative = "greater")
###summary(all_res$preds[which(all_res$celltype == "MG" & all_res$sample == "OSKR3")]);summary(all_res$preds[which(all_res$celltype == "MG" & all_res$sample == "OSKR4")])
t.test(all_res$preds[which(all_res$celltype == "RGC" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "RGC" & all_res$sample == "OSKR4")],alternative = "greater")
###summary(all_res$preds[which(all_res$celltype == "RGC" & all_res$sample == "OSKR3")]);summary(all_res$preds[which(all_res$celltype == "RGC" & all_res$sample == "OSKR4")])

t.test(all_res$preds[which(all_res$celltype == "AC" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "AC" & all_res$sample == "OSKR4")],alternative = "greater")
t.test(all_res$preds[which(all_res$celltype == "HC" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "HC" & all_res$sample == "OSKR4")],alternative = "greater")
t.test(all_res$preds[which(all_res$celltype == "Rod" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "Rod" & all_res$sample == "OSKR4")],alternative = "greater")
t.test(all_res$preds[which(all_res$celltype == "Cone" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "Cone" & all_res$sample == "OSKR4")],alternative = "greater")
t.test(all_res$preds[which(all_res$celltype == "BC" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "BC" & all_res$sample == "OSKR4")],alternative = "greater")
t.test(all_res$preds[which(all_res$celltype == "RPE" & all_res$sample == "OSKR3")],all_res$preds[which(all_res$celltype == "RPE" & all_res$sample == "OSKR4")],alternative = "greater")

#######
########
#######

######
####### compare Aging clocks between snRNA-seq and Xenium ##########
######


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load("Mouse_MG_model_Xen")
load("Mouse_BC_model_Xen")
load("Mouse_AC_model_Xen")
load("Mouse_HC_model_Xen")
load("Mouse_RGC_model_Xen")
load("Mouse_Rod_model_Xen")
load("Mouse_Cone_model_Xen")
load("Mouse_RPE_model_Xen")


setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Mouse_MG_model = readRDS("Mouse_MG_Xmodel.rds")
Mouse_AC_model = readRDS("Mouse_AC_Xmodel.rds")
Mouse_HC_model = readRDS("Mouse_HC_Xmodel.rds")
Mouse_RGC_model = readRDS("Mouse_RGC_Xmodel.rds")
Mouse_BC_model = readRDS("Mouse_BC_Xmodel.rds")
Mouse_Rod_model = readRDS("Mouse_Rod_Xmodel.rds")
Mouse_Cone_model = readRDS("Mouse_Cone_Xmodel.rds")
Mouse_RPE_model = readRDS("Mouse_RPE_Xmodel.rds")


########
########

Compare_weights <- function(snRNA_model = Mouse_MG_model_Xen,Xenium_model = Mouse_MG_model){
    ##########
    coefficients_sn = snRNA_model$coefficients
    coefficients_X = Xenium_model$coefficients
    ##########
    colnames(coefficients_sn) = c("gene","SnRNA_coef")
    colnames(coefficients_X) = c("gene","Xenium_coef")
    ###########
    coefficients = merge(coefficients_sn,coefficients_X)
    ##########
    coefficients_cl = coefficients[-which(coefficients$SnRNA_coef==0 & coefficients$Xenium_coef==0),]
    cor(coefficients_cl$SnRNA_coef,coefficients_cl$Xenium_coef)
    ##########
    coefficients_cl$score = coefficients_cl$SnRNA_coef*coefficients_cl$Xenium_coef
    coefficients_cl = coefficients_cl[order(coefficients_cl$score,decreasing=T),]
    ##########
    return(coefficients_cl)
}

library(ggplot2)
library(dplyr)
library(ggrepel)

####### for MG cells ###########
MG_res = Compare_weights(snRNA_model = Mouse_RPE_model_Xen,Xenium_model = Mouse_RPE_model)

gene_list1 <- MG_res[which(MG_res$SnRNA_coef > 0),]$gene[1:10]
gene_list2 <- MG_res[which(MG_res$SnRNA_coef < 0),]$gene[1:10]
gene_list <- c(gene_list1,gene_list2)



cor(MG_res$SnRNA_coef,MG_res$Xenium_coef)


k1 = which(MG_res$SnRNA_coef > 10)
k2 = which(MG_res$SnRNA_coef < -10)
MG_res$SnRNA_coef[k1] = 10
MG_res$SnRNA_coef[k2] = -10


plot_df <- MG_res %>%
  mutate(is_target = gene %in% gene_list)

library(ggplot2)
# 绘图
p <- ggplot(plot_df, aes(x = SnRNA_coef, y = Xenium_coef)) +
  geom_point(data = filter(plot_df, !is_target),
             color = "grey80", size = 1) +
  geom_point(data = filter(plot_df, is_target),
             color = "red", size = 1) +
  geom_text_repel(
    data = filter(plot_df, is_target),
    aes(label = gene),
    size          = 5,
    segment.size  = 0.5,
    segment.color = "black",
    box.padding   = 0.3,
    point.padding = 0.2,
    max.overlaps  = Inf
  ) +
  labs(
    x = "snRNAseq clock weights",
    y = "Xenium clock weights"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text    = element_text(color = "black"),
    axis.title   = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.margin  = margin(5, 5, 5, 5)
  ) +
  # 加上坐标轴范围限制
  coord_cartesian(xlim = c(-5, 5), ylim = c(-3, 3))

# 显示及保存
print(p)
ggsave("MG_highlighted_norepel_border_limits.png",
       plot = p, width = 6, height = 5, dpi = 300)



####
####











PredictAgeFromModel <- function(expr_mat,test_meta,fit_res) {
  library(glmnet)
  coef_df  = coef(fit_res, s = fit_res$lambda.min)
  coef_df  <- data.frame(
  gene        = rownames(coef_df),
  coefficient = coef_df[, 1],
  row.names   = NULL,
  stringsAsFactors = FALSE
  )
  # 1. 对 expr_mat 做与训练时相同的归一化和 log2(x+1)
  mat_norm <- expr_mat

  # 2. 对齐基因顺序：提取模型中出现的基因（含截距为 '(Intercept)'）
  #    coef_df 需含截距行，gene 为 '(Intercept)'
  #    expr_mat 只需行名与 coef gene 匹配
  mat_genes <- intersect(coef_df$gene, rownames(mat_norm))
  # 构建系数向量，确保顺序与表达矩阵一致
  coef_vec <- coef_df$coefficient[match(mat_genes, coef_df$gene)]
  # 加上截距
  intercept <- coef_df$coefficient[coef_df$gene == '(Intercept)']
  # 3. 子集表达矩阵并矩阵乘法
  sub_mat <- mat_norm[mat_genes, , drop = FALSE]
  # 预测 = t(sub_mat) %*% coef_vec + intercept
  preds <- as.numeric(crossprod(as.matrix(sub_mat), coef_vec) + intercept)
  #######
  test_meta$preds = preds
  return(test_meta)
}
















