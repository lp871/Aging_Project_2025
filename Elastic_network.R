#############
#############




library(glmnet)
library(ggplot2)
library(parallel)

#' 一轮 Elastic Net 年龄预测（含分层训练/测试拆分、alpha 搜索、归一化及输出数据集矩阵，加速版）
#'
#' @param mat         原始伪-bulk counts 矩阵，行为基因，列为样本
#' @param samples_tab 样本信息表，需包含列 'age'
#' @param Need_genes  需要纳入模型的基因向量
#' @param scale.factor 归一化缩放因子（默认 1e4）
#' @param train_frac  每个年龄层训练集比例（默认 0.7）
#' @param alpha.vec   要搜索的 alpha 值向量（默认 seq(0,1,0.1)）
#' @param nfolds      交叉验证折数（默认 5，减少计算量）
#' @param ncores      并行核心数（默认 detectCores()-1）
#' @param seed        随机种子（默认 42）
#' @return 列表，包括：
#'   - model: 最佳 alpha 和 lambda 下的 cv.glmnet 模型对象
#'   - best_alpha: 最佳 alpha
#'   - coefficients: 基因系数数据框
#'   - performance: 测试集上的 MSE 和 相关性
#'   - predictions: 含样本、实际 age、预测值和所属数据集
#'   - train_matrix: 归一化后训练集表达矩阵 (genes x samples)
#'   - test_matrix: 归一化后测试集表达矩阵 (genes x samples)
#'
seu = Mouse_seurat_knn_pseudobulk$MG
Need_genes = Mouse_features$MG

Elastic_net_one_round <- function(
  seu,
  Need_genes,
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
  mat_cl <- mat_norm[rownames(mat_norm) %in% Need_genes, , drop = FALSE]
  # 3. 分层抽样：按 age 层内随机拆分训练/测试集
  set.seed(seed)
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

  train_matrix <- mat_cl[, train_idx, drop = FALSE]
  test_matrix  <- mat_cl[, test_idx, drop = FALSE]

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


######
###### for mouse ！！！！ #################-----------------------------------------------------------------------------------------------------------------------------------------------
######


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat_knn_pseudobulk <- readRDS("Mouse_seurat_knn_pseudobulk_2025")
Mouse_features <- readRDS("Mouse_features_2025")

names(Mouse_seurat_knn_pseudobulk)
#######
#######
#######
seu = Mouse_seurat_knn_pseudobulk$MG
Need_genes = Mouse_features$MG

Mouse_MG_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$MG,Need_genes = Mouse_features$MG)
save(Mouse_MG_model,file="Mouse_MG_model")
Mouse_Rod_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$Rod,Need_genes = Mouse_features$Rod)
save(Mouse_Rod_model,file="Mouse_Rod_model")
Mouse_Cone_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$Cone,Need_genes = Mouse_features$Cone)
save(Mouse_Cone_model,file="Mouse_Cone_model")
Mouse_AC_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$AC,Need_genes = Mouse_features$AC)
save(Mouse_AC_model,file="Mouse_AC_model")
Mouse_HC_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$HC,Need_genes = Mouse_features$HC)
save(Mouse_HC_model,file="Mouse_HC_model")
Mouse_BC_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$BC,Need_genes = Mouse_features$BC)
save(Mouse_BC_model,file="Mouse_BC_model")
Mouse_RGC_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$RGC,Need_genes = Mouse_features$RGC)
save(Mouse_RGC_model,file="Mouse_RGC_model")
Mouse_RPE_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$RPE,Need_genes = Mouse_features$RPE)
save(Mouse_RPE_model,file="Mouse_RPE_model")
Mouse_Microglia_model = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$Microglia,Need_genes = Mouse_features$Microglia)
save(Mouse_Microglia_model,file="Mouse_Microglia_model")



############
############
############

######
###### for zebrafish ！！！！ #################-----------------------------------------------------------------------------------------------------------------------------------------------
######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_seurat_knn_pseudobulk <- readRDS("Zebrafish_seurat_knn_pseudobulk_2025")
Zebrafish_features <- readRDS("Zebrafish_features_2025")


#########

Zebrafish_MG_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$MG,Need_genes = Zebrafish_features$MG)
save(Zebrafish_MG_model,file="Zebrafish_MG_model")

Zebrafish_Rod_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$Rod,Need_genes = Zebrafish_features$Rod)
save(Zebrafish_Rod_model,file="Zebrafish_Rod_model")

Zebrafish_RPE_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$RPE,Need_genes = Zebrafish_features$RPE)
save(Zebrafish_RPE_model,file="Zebrafish_RPE_model")

Zebrafish_Microglia_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$Microglia,Need_genes = Zebrafish_features$Microglia)
save(Zebrafish_Microglia_model,file="Zebrafish_Microglia_model")

Zebrafish_AC_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$AC,Need_genes = Zebrafish_features$AC)
save(Zebrafish_AC_model,file="Zebrafish_AC_model")

Zebrafish_BC_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$BC,Need_genes = Zebrafish_features$BC)
save(Zebrafish_BC_model,file="Zebrafish_BC_model")

Zebrafish_HC_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$HC,Need_genes = Zebrafish_features$HC)
save(Zebrafish_HC_model,file="Zebrafish_HC_model")

Zebrafish_Cone_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$Cone,Need_genes = Zebrafish_features$Cone)
save(Zebrafish_Cone_model,file="Zebrafish_Cone_model")

Zebrafish_RGC_model = Elastic_net_one_round(seu=Zebrafish_seurat_knn_pseudobulk$RGC,Need_genes = Zebrafish_features$RGC)
save(Zebrafish_RGC_model,file="Zebrafish_RGC_model")

#######
#######
#######




#######-------------------------------------------------------------------------------------------
####### Plot the performence in R ############ --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######-------------------------------------------------------------------------------------------

#######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
load(file="Zebrafish_MG_model")
load(file="Zebrafish_Rod_model")
load(file="Zebrafish_Cone_model")
load(file="Zebrafish_AC_model")
load(file="Zebrafish_HC_model")
load(file="Zebrafish_RGC_model")
load(file="Zebrafish_RPE_model")
load(file="Zebrafish_Microglia_model")
load(file="Zebrafish_BC_model")

Zebrafish_MG_Aging_Clock_Res = Zebrafish_MG_model$predictions
Zebrafish_MG_Aging_Clock_Res$class = "MG"
cor(Zebrafish_MG_Aging_Clock_Res$predicted,Zebrafish_MG_Aging_Clock_Res$actual)

Zebrafish_Rod_Aging_Clock_Res = Zebrafish_Rod_model$predictions
Zebrafish_Rod_Aging_Clock_Res$class = "Rod"
cor(Zebrafish_Rod_Aging_Clock_Res$predicted,Zebrafish_Rod_Aging_Clock_Res$actual)

Zebrafish_Cone_Aging_Clock_Res = Zebrafish_Cone_model$predictions
Zebrafish_Cone_Aging_Clock_Res$class = "Cone"
cor(Zebrafish_Cone_Aging_Clock_Res$predicted,Zebrafish_Cone_Aging_Clock_Res$actual)

Zebrafish_RGC_Aging_Clock_Res = Zebrafish_RGC_model$predictions
Zebrafish_RGC_Aging_Clock_Res$class = "RGC"
cor(Zebrafish_RGC_Aging_Clock_Res$predicted,Zebrafish_RGC_Aging_Clock_Res$actual)

Zebrafish_BC_Aging_Clock_Res = Zebrafish_BC_model$predictions
Zebrafish_BC_Aging_Clock_Res$class = "BC"
cor(Zebrafish_BC_Aging_Clock_Res$predicted,Zebrafish_BC_Aging_Clock_Res$actual)

Zebrafish_AC_Aging_Clock_Res = Zebrafish_AC_model$predictions
Zebrafish_AC_Aging_Clock_Res$class = "AC"
cor(Zebrafish_AC_Aging_Clock_Res$predicted,Zebrafish_AC_Aging_Clock_Res$actual)

Zebrafish_HC_Aging_Clock_Res = Zebrafish_HC_model$predictions
Zebrafish_HC_Aging_Clock_Res$class = "HC"
cor(Zebrafish_HC_Aging_Clock_Res$predicted,Zebrafish_HC_Aging_Clock_Res$actual)

Zebrafish_RPE_Aging_Clock_Res = Zebrafish_RPE_model$predictions
Zebrafish_RPE_Aging_Clock_Res$class = "RPE"
cor(Zebrafish_RPE_Aging_Clock_Res$predicted,Zebrafish_RPE_Aging_Clock_Res$actual)

Zebrafish_Microglia_Aging_Clock_Res = Zebrafish_Microglia_model$predictions
Zebrafish_Microglia_Aging_Clock_Res$class = "Microglia"
cor(Zebrafish_Microglia_Aging_Clock_Res$predicted,Zebrafish_Microglia_Aging_Clock_Res$actual)


All_merge_res = rbind(Zebrafish_MG_Aging_Clock_Res,Zebrafish_Cone_Aging_Clock_Res,Zebrafish_Rod_Aging_Clock_Res,Zebrafish_RGC_Aging_Clock_Res,Zebrafish_BC_Aging_Clock_Res,Zebrafish_AC_Aging_Clock_Res,Zebrafish_HC_Aging_Clock_Res,Zebrafish_RPE_Aging_Clock_Res,Zebrafish_Microglia_Aging_Clock_Res)

##########
head(All_merge_res)
All_merge_res$Samples = All_merge_res$actual
table(All_merge_res$Samples)
All_merge_res$Samples = factor(All_merge_res$Samples,levels=c(1,3,6,12,18,22,24,30,36,48))
All_merge_res$class = factor(All_merge_res$class,levels=c("MG","RGC","AC","HC","Rod","Cone","BC","Microglia","RPE"))


#########
#########
library(ggplot2)
library(viridis)

ggplot(All_merge_res, aes(x = Samples, y = predicted, color = actual, fill = actual)) +
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
  scale_y_continuous(limits = c(-5, 55)) +
  theme_classic() +
  theme(
    strip.text          = element_text(size = 16, face = "bold"),
    axis.text.x         = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y         = element_text(size = 14),
    axis.ticks.length   = unit(0.3, "cm"),    # 刻度线变长
    axis.title         = element_text(size = 16),
    panel.border        = element_rect(color = "black", fill = NA, size = 1)
  ) +
  xlab("") +
  ylab("")

ggsave("AllZebrafish_results.png", height = 3, width = 20)


#######
####### Next we will see the Mouse ########
#######



setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load(file="Mouse_MG_model")
load(file="Mouse_Rod_model")
load(file="Mouse_Cone_model")
load(file="Mouse_AC_model")
load(file="Mouse_HC_model")
load(file="Mouse_RGC_model")
load(file="Mouse_RPE_model")
load(file="Mouse_Microglia_model")
load(file="Mouse_BC_model")

Mouse_MG_Aging_Clock_Res = Mouse_MG_model$predictions
Mouse_MG_Aging_Clock_Res$class = "MG"
cor(Mouse_MG_Aging_Clock_Res$predicted,Mouse_MG_Aging_Clock_Res$actual)

Mouse_Rod_Aging_Clock_Res = Mouse_Rod_model$predictions
Mouse_Rod_Aging_Clock_Res$class = "Rod"
cor(Mouse_Rod_Aging_Clock_Res$predicted,Mouse_Rod_Aging_Clock_Res$actual)

Mouse_Cone_Aging_Clock_Res = Mouse_Cone_model$predictions
Mouse_Cone_Aging_Clock_Res$class = "Cone"
cor(Mouse_Cone_Aging_Clock_Res$predicted,Mouse_Cone_Aging_Clock_Res$actual)

Mouse_RGC_Aging_Clock_Res = Mouse_RGC_model$predictions
Mouse_RGC_Aging_Clock_Res$class = "RGC"
cor(Mouse_RGC_Aging_Clock_Res$predicted,Mouse_RGC_Aging_Clock_Res$actual)

Mouse_BC_Aging_Clock_Res = Mouse_BC_model$predictions
Mouse_BC_Aging_Clock_Res$class = "BC"
cor(Mouse_BC_Aging_Clock_Res$predicted,Mouse_BC_Aging_Clock_Res$actual)

Mouse_AC_Aging_Clock_Res = Mouse_AC_model$predictions
Mouse_AC_Aging_Clock_Res$class = "AC"
cor(Mouse_AC_Aging_Clock_Res$predicted,Mouse_AC_Aging_Clock_Res$actual)

Mouse_HC_Aging_Clock_Res = Mouse_HC_model$predictions
Mouse_HC_Aging_Clock_Res$class = "HC"
cor(Mouse_HC_Aging_Clock_Res$predicted,Mouse_HC_Aging_Clock_Res$actual)

Mouse_RPE_Aging_Clock_Res = Mouse_RPE_model$predictions
Mouse_RPE_Aging_Clock_Res$class = "RPE"
cor(Mouse_RPE_Aging_Clock_Res$predicted,Mouse_RPE_Aging_Clock_Res$actual)

Mouse_Microglia_Aging_Clock_Res = Mouse_Microglia_model$predictions
Mouse_Microglia_Aging_Clock_Res$class = "Microglia"
cor(Mouse_Microglia_Aging_Clock_Res$predicted,Mouse_Microglia_Aging_Clock_Res$actual)


All_merge_res = rbind(Mouse_MG_Aging_Clock_Res,Mouse_Cone_Aging_Clock_Res,Mouse_Rod_Aging_Clock_Res,Mouse_RGC_Aging_Clock_Res,Mouse_BC_Aging_Clock_Res,Mouse_AC_Aging_Clock_Res,Mouse_HC_Aging_Clock_Res,Mouse_RPE_Aging_Clock_Res,Mouse_Microglia_Aging_Clock_Res)

##########
head(All_merge_res)
All_merge_res$Samples = All_merge_res$actual
table(All_merge_res$Samples)
All_merge_res$Samples = factor(All_merge_res$Samples,levels=c(5,12,17,32,49,68,91,108,120))
All_merge_res$class = factor(All_merge_res$class,levels=c("MG","RGC","AC","HC","Rod","Cone","BC","Microglia","RPE"))

library(ggplot2)
library(viridis)
ggplot(All_merge_res,aes(x=Samples,y=predicted,color=actual,fill=actual)) + facet_wrap(~ class, nrow = 1) + geom_violin() + geom_boxplot(width=0.25,fill="grey",color="black",outlier.shape = NA) + scale_color_viridis() + scale_fill_viridis() + theme_classic() + scale_y_continuous(limits=c(-10,135)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("AllMouse_results.png",height=3,width=20) 

  
ggplot(All_merge_res, aes(x = Samples, y = predicted, color = actual, fill = actual)) +
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
  scale_y_continuous(limits = c(-10,135)) +
  theme_classic() +
  theme(
    strip.text          = element_text(size = 16, face = "bold"),
    axis.text.x         = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y         = element_text(size = 14),
    axis.ticks.length   = unit(0.3, "cm"),    # 刻度线变长
    axis.title         = element_text(size = 16),
    panel.border        = element_rect(color = "black", fill = NA, size = 1)
  ) +
  xlab("") +
  ylab("")

ggsave("AllMouse_results.png", height = 3, width = 20)



########
########
######## Next for Human model:##########
########
########



ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R


setwd("/zp1/data/share/Human_aging_new")
library(Seurat)


Human_RPE_seurat_knn_pseudobulk <- readRDS("Human_RPE_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

names(Human_features)
head(Human_RPE_seurat_knn_pseudobulk@meta.data)

#######
#######
Human_RPE_seurat_knn_pseudobulk_F = Human_RPE_seurat_knn_pseudobulk[,which(Human_RPE_seurat_knn_pseudobulk$gender == "Female")]
Human_RPE_seurat_knn_pseudobulk_M = Human_RPE_seurat_knn_pseudobulk[,which(Human_RPE_seurat_knn_pseudobulk$gender == "Male")]

#### 1 ####
Human_RPE_model_F = Elastic_net_one_round_Human(seu=Human_RPE_seurat_knn_pseudobulk_F,Need_genes = Human_features$RPE_F)
save(Human_RPE_model_F,file="Human_RPE_model_F")

head(Human_RPE_model_F$)

#### 2 ####
Human_RPE_model_M = Elastic_net_one_round_Human(seu=Human_RPE_seurat_knn_pseudobulk_M,Need_genes = Human_features$RPE_M)
save(Human_RPE_model_M,file="Human_RPE_model_M")


########------ Microglia ########

setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_Microglia_seurat_knn_pseudobulk <- readRDS("Human_Microglia_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_Microglia_seurat_knn_pseudobulk_F = Human_Microglia_seurat_knn_pseudobulk[,which(Human_Microglia_seurat_knn_pseudobulk$gender == "female")]
Human_Microglia_seurat_knn_pseudobulk_M = Human_Microglia_seurat_knn_pseudobulk[,which(Human_Microglia_seurat_knn_pseudobulk$gender == "male")]

#######
#######
Human_Microglia_model_F = Elastic_net_one_round_Human(seu=Human_Microglia_seurat_knn_pseudobulk_F,Need_genes = Human_features$Microglia_F)
save(Human_Microglia_model_F,file="Human_Microglia_model_F")
Human_Microglia_model_M = Elastic_net_one_round_Human(seu=Human_Microglia_seurat_knn_pseudobulk_M,Need_genes = Human_features$Microglia_M)
save(Human_Microglia_model_M,file="Human_Microglia_model_M")

#######
#######
#######
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_Rod_seurat_knn_pseudobulk <- readRDS("Human_Rod_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_Rod_seurat_knn_pseudobulk_F = Human_Rod_seurat_knn_pseudobulk[,which(Human_Rod_seurat_knn_pseudobulk$gender == "female")]
Human_Rod_seurat_knn_pseudobulk_M = Human_Rod_seurat_knn_pseudobulk[,which(Human_Rod_seurat_knn_pseudobulk$gender == "male")]

Human_Rod_model_F = Elastic_net_one_round_Human(seu=Human_Rod_seurat_knn_pseudobulk_F,Need_genes = Human_features$Rod_F)
save(Human_Rod_model_F,file="Human_Rod_model_F")
Human_Rod_model_M = Elastic_net_one_round_Human(seu=Human_Rod_seurat_knn_pseudobulk_M,Need_genes = Human_features$Rod_M)
save(Human_Rod_model_M,file="Human_Rod_model_M")

#######
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_Cone_seurat_knn_pseudobulk <- readRDS("Human_Cone_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_Cone_seurat_knn_pseudobulk_F = Human_Cone_seurat_knn_pseudobulk[,which(Human_Cone_seurat_knn_pseudobulk$gender == "female")]
Human_Cone_seurat_knn_pseudobulk_M = Human_Cone_seurat_knn_pseudobulk[,which(Human_Cone_seurat_knn_pseudobulk$gender == "male")]

Human_Cone_model_F = Elastic_net_one_round_Human(seu=Human_Cone_seurat_knn_pseudobulk_F,Need_genes = Human_features$Cone_F)
save(Human_Cone_model_F,file="Human_Cone_model_F")
Human_Cone_model_M = Elastic_net_one_round_Human(seu=Human_Cone_seurat_knn_pseudobulk_M,Need_genes = Human_features$Cone_M)
save(Human_Cone_model_M,file="Human_Cone_model_M")

#######
#######
#######
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_BC_seurat_knn_pseudobulk <- readRDS("Human_BC_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_BC_seurat_knn_pseudobulk_F = Human_BC_seurat_knn_pseudobulk[,which(Human_BC_seurat_knn_pseudobulk$gender == "female")]
Human_BC_seurat_knn_pseudobulk_M = Human_BC_seurat_knn_pseudobulk[,which(Human_BC_seurat_knn_pseudobulk$gender == "male")]

Human_BC_model_F = Elastic_net_one_round_Human(seu=Human_BC_seurat_knn_pseudobulk_F,Need_genes = Human_features$BC_F)
save(Human_BC_model_F,file="Human_BC_model_F")
Human_BC_model_M = Elastic_net_one_round_Human(seu=Human_BC_seurat_knn_pseudobulk_M,Need_genes = Human_features$BC_M)
save(Human_BC_model_M,file="Human_BC_model_M")
#######
#######
#######
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_AC_seurat_knn_pseudobulk <- readRDS("Human_AC_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_AC_seurat_knn_pseudobulk_F = Human_AC_seurat_knn_pseudobulk[,which(Human_AC_seurat_knn_pseudobulk$gender == "female")]
Human_AC_seurat_knn_pseudobulk_M = Human_AC_seurat_knn_pseudobulk[,which(Human_AC_seurat_knn_pseudobulk$gender == "male")]

Human_AC_model_F = Elastic_net_one_round_Human(seu=Human_AC_seurat_knn_pseudobulk_F,Need_genes = Human_features$AC_F)
save(Human_AC_model_F,file="Human_AC_model_F")
Human_AC_model_M = Elastic_net_one_round_Human(seu=Human_AC_seurat_knn_pseudobulk_M,Need_genes = Human_features$AC_M)
save(Human_AC_model_M,file="Human_AC_model_M")
#######
#######
#######
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_RGC_seurat_knn_pseudobulk <- readRDS("Human_RGC_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_RGC_seurat_knn_pseudobulk_F = Human_RGC_seurat_knn_pseudobulk[,which(Human_RGC_seurat_knn_pseudobulk$gender == "female")]
Human_RGC_seurat_knn_pseudobulk_M = Human_RGC_seurat_knn_pseudobulk[,which(Human_RGC_seurat_knn_pseudobulk$gender == "male")]

Human_RGC_model_F = Elastic_net_one_round_Human(seu=Human_RGC_seurat_knn_pseudobulk_F,Need_genes = Human_features$RGC_F)
save(Human_RGC_model_F,file="Human_RGC_model_F")
Human_RGC_model_M = Elastic_net_one_round_Human(seu=Human_RGC_seurat_knn_pseudobulk_M,Need_genes = Human_features$RGC_M)
save(Human_RGC_model_M,file="Human_RGC_model_M")
#######
#######
#######
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_MG_seurat_knn_pseudobulk <- readRDS("Human_MG_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_MG_seurat_knn_pseudobulk_F = Human_MG_seurat_knn_pseudobulk[,which(Human_MG_seurat_knn_pseudobulk$gender == "female")]
Human_MG_seurat_knn_pseudobulk_M = Human_MG_seurat_knn_pseudobulk[,which(Human_MG_seurat_knn_pseudobulk$gender == "male")]

Human_MG_model_F = Elastic_net_one_round_Human(seu=Human_MG_seurat_knn_pseudobulk_F,Need_genes = Human_features$MG_F)
save(Human_MG_model_F,file="Human_MG_model_F")
Human_MG_model_M = Elastic_net_one_round_Human(seu=Human_MG_seurat_knn_pseudobulk_M,Need_genes = Human_features$MG_M)
save(Human_MG_model_M,file="Human_MG_model_M")

#######
#######
#######
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)
Human_HC_seurat_knn_pseudobulk <- readRDS("Human_HC_pseudo_bulk_2025")
Human_features <- readRDS("Human_features_2025")

Human_HC_seurat_knn_pseudobulk_F = Human_HC_seurat_knn_pseudobulk[,which(Human_HC_seurat_knn_pseudobulk$gender == "female")]
Human_HC_seurat_knn_pseudobulk_M = Human_HC_seurat_knn_pseudobulk[,which(Human_HC_seurat_knn_pseudobulk$gender == "male")]

Human_HC_model_F = Elastic_net_one_round_Human(seu=Human_HC_seurat_knn_pseudobulk_F,Need_genes = Human_features$HC_F)
save(Human_HC_model_F,file="Human_HC_model_F")
Human_HC_model_M = Elastic_net_one_round_Human(seu=Human_HC_seurat_knn_pseudobulk_M,Need_genes = Human_features$HC_M)
save(Human_HC_model_M,file="Human_HC_model_M")


#######

#######
####### 检查下 function ####### 
####### 对每一个 donor 和 每一个 age 做 weights #######
####### seu = Human_RPE_seurat_knn_pseudobulk_F
####### Need_genes = Human_features$RPE_F

Elastic_net_one_round_Human <- function(
  seu,
  Need_genes,
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
  mat_cl <- mat_norm[rownames(mat_norm) %in% Need_genes, , drop = FALSE]
  # 3. 分层抽样：按 age 层内随机拆分训练/测试集
  set.seed(seed)
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

  train_matrix <- mat_cl[, train_idx, drop = FALSE]
  test_matrix  <- mat_cl[, test_idx, drop = FALSE]

  x_train <- t(train_matrix)
  y_train <- as.numeric(ages[train_idx])
  x_test  <- t(test_matrix)
  y_test  <- as.numeric(ages[test_idx])

  # 4. 训练集权重：确保每个 age2 均匀分配给每个 donor，且同一 donor 的 cell 权重之和相等
  age2_train <- seu$age2[train_idx]
  donor_train <- seu$donor[train_idx]
  # 统计每个 age2 下的独特 donor 数
  donors_per_age2 <- tapply(donor_train, age2_train, function(x) length(unique(x)))
  # 统计每个 donor 在训练集中的 cell 数
  cells_per_donor <- table(donor_train)
  # 计算每个 cell 的初始权重： donor_share / cells_per_donor
  # donor_share = 1 / number of donors in that age2
  donor_share <- 1 / donors_per_age2[as.character(age2_train)]
  w_train <- donor_share / as.numeric(cells_per_donor[donor_train])
  # 标准化：总和等于训练集样本数
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

  # 9. 保存
  return(list(
    model        = fit,
    best_alpha   = best_alpha,
    coefficients = coef_df,
    performance  = list(test_mse = test_mse, test_cor = test_cor),
    predictions  = pred_df,
    train_matrix_namecell = colnames(train_matrix),
    test_matrix_namecell  = colnames(train_matrix)
  ))
}



########------ Next Plot the Human model !!! #######

setwd("/zp1/data/share/Human_aging_new")
load(file="Human_MG_model_F")
load(file="Human_Rod_model_F")
load(file="Human_Cone_model_F")
load(file="Human_AC_model_F")
load(file="Human_HC_model_F")
load(file="Human_BC_model_F")
load(file="Human_RPE_model_F")
load(file="Human_RGC_model_F")
load(file="Human_Microglia_model_F")
load(file="Human_MG_model_M")
load(file="Human_Rod_model_M")
load(file="Human_Cone_model_M")
load(file="Human_AC_model_M")
load(file="Human_HC_model_M")
load(file="Human_BC_model_M")
load(file="Human_RPE_model_M")
load(file="Human_RGC_model_M")
load(file="Human_Microglia_model_M")

MG_Aging_Clock_Res = rbind(Human_MG_model_F$predictions,Human_MG_model_M$predictions)
MG_Aging_Clock_Res$class = "MG"
cor(MG_Aging_Clock_Res$predicted,MG_Aging_Clock_Res$actual)

Rod_Aging_Clock_Res = rbind(Human_Rod_model_F$predictions,Human_Rod_model_M$predictions)
Rod_Aging_Clock_Res$class = "Rod"
cor(Rod_Aging_Clock_Res$predicted,Rod_Aging_Clock_Res$actual)

Cone_Aging_Clock_Res = rbind(Human_Cone_model_F$predictions,Human_Cone_model_M$predictions)
Cone_Aging_Clock_Res$class = "Cone"
cor(Cone_Aging_Clock_Res$predicted,Cone_Aging_Clock_Res$actual)

BC_Aging_Clock_Res = rbind(Human_BC_model_F$predictions,Human_BC_model_M$predictions)
BC_Aging_Clock_Res$class = "BC"
cor(BC_Aging_Clock_Res$predicted,BC_Aging_Clock_Res$actual)

AC_Aging_Clock_Res = rbind(Human_AC_model_F$predictions,Human_AC_model_M$predictions)
AC_Aging_Clock_Res$class = "AC"
cor(AC_Aging_Clock_Res$predicted,AC_Aging_Clock_Res$actual)

HC_Aging_Clock_Res = rbind(Human_HC_model_F$predictions,Human_HC_model_M$predictions)
HC_Aging_Clock_Res$class = "HC"
cor(HC_Aging_Clock_Res$predicted,HC_Aging_Clock_Res$actual)

RGC_Aging_Clock_Res = rbind(Human_RGC_model_F$predictions,Human_RGC_model_M$predictions)
RGC_Aging_Clock_Res$class = "RGC"
cor(RGC_Aging_Clock_Res$predicted,RGC_Aging_Clock_Res$actual)

RPE_Aging_Clock_Res = rbind(Human_RPE_model_F$predictions,Human_RPE_model_M$predictions)
RPE_Aging_Clock_Res$class = "RPE"
cor(RPE_Aging_Clock_Res$predicted,RPE_Aging_Clock_Res$actual)

Microglia_Aging_Clock_Res = rbind(Human_Microglia_model_F$predictions,Human_Microglia_model_M$predictions)
Microglia_Aging_Clock_Res$class = "Microglia"
cor(Microglia_Aging_Clock_Res$predicted,Microglia_Aging_Clock_Res$actual)


All_merge_res = rbind(MG_Aging_Clock_Res,Cone_Aging_Clock_Res,Rod_Aging_Clock_Res,RGC_Aging_Clock_Res,BC_Aging_Clock_Res,AC_Aging_Clock_Res,HC_Aging_Clock_Res,RPE_Aging_Clock_Res,Microglia_Aging_Clock_Res)

table(All_merge_res$class)
##########
head(All_merge_res)

Add_interval_Human <- function(All_merge_res){
    ######## rm unknown gender #####
    ########
    meta = All_merge_res
    meta$age = meta$actual
    ######## summary(meta$age)
    k1 = which(meta$age >= 0 & meta$age < 20)
    k1.5 = which(meta$age >= 20 & meta$age < 30)
    k2 = which(meta$age >= 30 & meta$age < 45)
    k3 = which(meta$age >= 45 & meta$age < 55)
    k5 = which(meta$age >= 55 & meta$age < 65)
    k7 = which(meta$age >= 65 & meta$age < 75)
    k8 = which(meta$age >= 75 & meta$age < 85)
    k9 = which(meta$age >= 85 & meta$age < 95)
    #####
    meta$age2 = 0
    meta$age2[k1] = 10
    meta$age2[k1.5] = 25
    meta$age2[k2] = 37.5
    meta$age2[k3] = 50
    meta$age2[k5] = 60
    meta$age2[k7] = 70
    meta$age2[k8] = 80
    meta$age2[k9] = 90
    ########
    ########
    ########
    return(meta)
}


All_merge_res = Add_interval_Human(All_merge_res)
table(All_merge_res$age2)

##########
##########
All_merge_res$Samples = All_merge_res$age2
table(All_merge_res$Samples)

All_merge_res$Samples = factor(All_merge_res$Samples,levels=c(10,25,37.5,50,60,70,80,90))
All_merge_res$class = factor(All_merge_res$class,levels=c("MG","RGC","AC","HC","Rod","Cone","BC","Microglia","RPE"))

library(ggplot2)
library(viridis)
ggplot(All_merge_res,aes(x=Samples,y=predicted,color=age2,fill=age2)) + facet_wrap(~ class, nrow = 1) + geom_violin() + geom_boxplot(width=0.25,fill="grey",color="black",outlier.shape = NA) + scale_color_viridis() + scale_fill_viridis() + theme_classic() + scale_y_continuous(limits=c(-10,110)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("AllHuman_results.png",height=3,width=20) 


                 
ggplot(All_merge_res, aes(x = Samples, y = predicted, color = age2, fill = age2)) +
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
  scale_y_continuous(limits = c(-10,110)) +
  theme_classic() +
  theme(
    strip.text          = element_text(size = 16, face = "bold"),
    axis.text.x         = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y         = element_text(size = 14),
    axis.ticks.length   = unit(0.3, "cm"),    # 刻度线变长
    axis.title         = element_text(size = 16),
    panel.border        = element_rect(color = "black", fill = NA, size = 1)
  ) +
  xlab("") +
  ylab("")

ggsave("AllHuman_results.png", height = 3, width = 20)
######
#####
######
##### Done!!!! #####
######

























