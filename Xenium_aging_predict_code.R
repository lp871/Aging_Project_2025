##########
#########
##########
#########
##########

#########
######### load the Xenium datasets !!!!!! ##########
#########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

##### Not OSK samples !!! #####

library(Seurat)

setwd("/zp1/data/plyu3/Xenium5k_Mouse/New")
Xenium5k_Mouse_Obj <- readRDS("obj_cl2_2024")

########## see the celldenisty ########
head(Xenium5k_Mouse_Obj@meta.data)

##########

CTs <- names(table(Xenium5k_Mouse_Obj$celltype))
CTs_Need = CTs[c(1,4,5,8,12,14,15,16)] ##### ignore RPEs ###

##########
Xenium5k_Mouse_Obj_Cl = Xenium5k_Mouse_Obj[,which(Xenium5k_Mouse_Obj$celltype %in% CTs_Need == T)]
Xenium5k_Mouse_Obj_Cl_Sp = SplitObject(Xenium5k_Mouse_Obj_Cl,"celltype")



Clean_the_cells_by_counts <- function(Xenium5k_Mouse_Obj_Cl_Sp){
  #####
  Out_list <- c()
  ######
  for(i in 1:length(Xenium5k_Mouse_Obj_Cl_Sp)){
      ########
      tmp_Sp = Xenium5k_Mouse_Obj_Cl_Sp[[i]]
      ########
      lowers = quantile(tmp_Sp$nFeature_RNA,0.15)
      highers = quantile(tmp_Sp$nFeature_RNA,0.85)
      print(names(Xenium5k_Mouse_Obj_Cl_Sp)[i])
      print(lowers)
      print(highers)
      ########
      tmp_Index = which(tmp_Sp$nFeature_RNA > lowers & tmp_Sp$nFeature_RNA < highers)
      tmp_Sp = tmp_Sp[,tmp_Index]
      ########
      print(summary(tmp_Sp$nFeature_RNA))
      Out_list <- c(Out_list,list(tmp_Sp))
  }
  #####
  names(Out_list) <- names(Xenium5k_Mouse_Obj_Cl_Sp)
  #####
  return(Out_list)
}

Xenium5k_Mouse_Obj_Cl_Sp_clocks <- Clean_the_cells_by_counts(Xenium5k_Mouse_Obj_Cl_Sp)

Prepare_Matrix_And_Meta <- function(Xenium5k_Mouse_Obj_Cl_Sp_clocks){
  #########
  Results <- c()
  #########
  for(i in 1:length(Xenium5k_Mouse_Obj_Cl_Sp_clocks)){
    ########
    tmp = Xenium5k_Mouse_Obj_Cl_Sp_clocks[[i]]
    ########
    tmp_matrix = tmp[['SCT']]@counts
    tmp_meta = tmp@meta.data
    tmp_meta = tmp_meta[,c("sample_total","sample1","celltype")]
    tmp_meta$cellid = rownames(tmp_meta)
    tmp_meta$age = 0
    vector1 = c("R1","R2","R3","R4")
    vector2 = c("5wk","12wk","13mo","28mo")
    vector3 = c("5wk","12wk","52wk","112wk")
    ########
    m = match(tmp_meta$sample1,vector1)
    tmp_meta$age = vector3[m]
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
  names(Results) <- names(Xenium5k_Mouse_Obj_Cl_Sp_clocks)
  #########
  return(Results)
}

Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared = Prepare_Matrix_And_Meta(Xenium5k_Mouse_Obj_Cl_Sp_clocks)
saveRDS(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared,file="Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared_2025")



################
################ load the mouse model !!!! ########
################

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")

load("Mouse_MG_model")
load("Mouse_RGC_model")
load("Mouse_AC_model")
load("Mouse_HC_model")
load("Mouse_Rod_model")
load("Mouse_Cone_model")
load("Mouse_BC_model")
load("Mouse_RPE_model")

############

MG_test_mat = Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$MG$Matrix
MG_test_meta = Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$MG$Meta

fit_res = Mouse_MG_model$model
expr_mat = MG_test_mat
test_meta = MG_test_meta

######
###### 
######

MG_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$MG$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$MG$Meta,Mouse_MG_model_Xen$model)
MG_predict_res$class="MG"

Rod_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$Rod$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$Rod$Meta,Mouse_Rod_model_Xen$model)
Rod_predict_res$class = "Rod"
Cone_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$Cone$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$Cone$Meta,Mouse_Cone_model_Xen$model)
Cone_predict_res$class = "Cone"
AC_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$AC$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$AC$Meta,Mouse_AC_model_Xen$model)
AC_predict_res$class = "AC"
BC_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$BC$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$BC$Meta,Mouse_BC_model_Xen$model)
BC_predict_res$class = "BC"
HC_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$HC$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$HC$Meta,Mouse_HC_model_Xen$model)
HC_predict_res$class = "HC"
RPE_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$RPE$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$RPE$Meta,Mouse_RPE_model_Xen$model)
RPE_predict_res$class = "RPE"
RGC_predict_res = PredictAgeFromModel(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$RGC$Matrix,Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$RGC$Meta,Mouse_RGC_model_Xen$model)
RGC_predict_res$class = "RGC"

#########
#########
#########
#########
all_res = rbind(MG_predict_res,Rod_predict_res,Cone_predict_res,AC_predict_res,BC_predict_res,HC_predict_res,RGC_predict_res,RPE_predict_res)
all_res$class = factor(all_res$class,levels=c("MG","RGC","AC","HC","Rod","Cone","BC","RPE"))
all_res$age = factor(all_res$age,levels=c("5wk","12wk","52wk","112wk"))
########
library(ggplot2)
ggplot(all_res,aes(x=age ,y=preds,color=age)) + facet_wrap(~ class, nrow = 1,scales = "free_y") + scale_y_continuous(expand=c(0,0)) + geom_boxplot(outlier.shape = NA,size=1)+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("Mouse_Xenium.png",height=3,width=20) 

#####
all_res$age = as.numeric(gsub("wk","",all_res$age))

all_res_sp = split(all_res,all_res$celltype)

for(i in 1:length(all_res_sp)){
    ####
    print(names(all_res_sp)[i])
    print(cor(all_res_sp[[i]]$age,all_res_sp[[i]]$preds))
}


######
######
######


tapply(RGC_predict_res$preds,RGC_predict_res$age,mean)

#######
head(MG_predict_res)







########
######## Next we will find the Model for zebrafish injury NMDA and prediction results ###################
########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R
setwd("/zp1/data/plyu3/Aging_Clocks_Final")

##########
Injury_merge_seurat_all_2023Dec <- readRDS("Injury_merge_seurat_all_2023Dec(1)")

########## OK!!! let us see the different genes expression #######

head(Injury_merge_seurat_all_2023Dec@meta.data)
table(Injury_merge_seurat_all_2023Dec$sample)

######### Extract the NMDA samples ##########
all_samples = names(table(Injury_merge_seurat_all_2023Dec$sample))
NMDA_samples_seurat <- Injury_merge_seurat_all_2023Dec[,which(Injury_merge_seurat_all_2023Dec$sample %in% all_samples[c(1,8:13)])]

#########
table(NMDA_samples_seurat$sample)

CTs <- names(table(NMDA_samples_seurat$new_celltype))
CTs_Need = c("Rod","Cone","RestMG") ##### ignore RPEs ###

##########
library(Seurat)
NMDA_samples_seurat_cl = NMDA_samples_seurat[,which(NMDA_samples_seurat$new_celltype %in% CTs_Need == T)]
NMDA_samples_seurat_cl$celltype = NMDA_samples_seurat_cl$new_celltype
NMDA_samples_seurat_cl_sp = SplitObject(NMDA_samples_seurat_cl,"celltype")

NMDA_samples_seurat_cl_sp_clocks <- Clean_the_cells_by_counts(NMDA_samples_seurat_cl_sp)
#####
#####
NMDA_samples_seurat_cl_sp_clocks_input = Prepare_Matrix_And_Meta_Injury(NMDA_samples_seurat_cl_sp_clocks)


Prepare_Matrix_And_Meta_Injury <- function(NMDA_samples_seurat_cl_sp_clocks){
  #########
  Results <- c()
  #########
  for(i in 1:length(NMDA_samples_seurat_cl_sp_clocks)){
    ########
    tmp = NMDA_samples_seurat_cl_sp_clocks[[i]]
    ########
    tmp_matrix = tmp[['RNA']]$counts
    #########
    rownames(tmp_matrix) = sapply(strsplit(rownames(tmp_matrix),split="~~"),function(x) x[[2]])
    dup = which(duplicated(rownames(tmp_matrix)) == T)
    tmp_matrix = tmp_matrix[-dup,]
    #########
    tmp_meta = tmp@meta.data
    tmp_meta = tmp_meta[,c("sample","celltype")]
    tmp_meta$cellid = rownames(tmp_meta)
    ########
    ######## norm the matrix #########
    ########
    tmp_matrix_sum = Matrix::colSums(tmp_matrix)
    tmp_matrix_sum_factor = tmp_matrix_sum / 1e4
    ########
    tmp_matrix_norm = sweep(tmp_matrix,2,tmp_matrix_sum_factor,FUN= "/")
    ######## change the log2 ##############
    tmp_matrix_norm_Log2 = log2(tmp_matrix_norm+1)
    ########
    reslist = list(Matrix = tmp_matrix_norm_Log2, Meta = tmp_meta)
    ########
    Results <- c(Results,list(reslist))
  }
  #########
  names(Results) <- names(NMDA_samples_seurat_cl_sp_clocks)
  #########
  return(Results)
}

#########
######### Next we will prdict the results from the model !!! ##########
#########

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")

load("Zebrafish_MG_model")
load("Zebrafish_Rod_model")
load("Zebrafish_Cone_model")
load("Zebrafish_RGC_model")


Injury_MG_res = PredictAgeFromModel(expr_mat=NMDA_samples_seurat_cl_sp_clocks_input$RestMG$Matrix,test_meta=NMDA_samples_seurat_cl_sp_clocks_input$RestMG$Meta,fit_res=Zebrafish_MG_model$model)
Injury_Rod_res = PredictAgeFromModel(expr_mat=NMDA_samples_seurat_cl_sp_clocks_input$Rod$Matrix,test_meta=NMDA_samples_seurat_cl_sp_clocks_input$Rod$Meta,fit_res=Zebrafish_Rod_model$model)
Injury_Cone_res = PredictAgeFromModel(expr_mat=NMDA_samples_seurat_cl_sp_clocks_input$Cone$Matrix,test_meta=NMDA_samples_seurat_cl_sp_clocks_input$Cone$Meta,fit_res=Zebrafish_Cone_model$model)
Injury_RGC_res = PredictAgeFromModel(expr_mat=NMDA_samples_seurat_cl_sp_clocks_input$RGC$Matrix,test_meta=NMDA_samples_seurat_cl_sp_clocks_input$RGC$Meta,fit_res=Zebrafish_RGC_model$model)

head(Injury_MG_res)
tapply(Injury_MG_res$preds,Injury_MG_res$sample,mean)
table()

Injury_MG_res$sample = factor(Injury_MG_res$sample,levels=c("Control_2","NMDA_36hr","NMDA_54hr","NMDA_72hr","NMDA_96hr","NMDA_7D","NMDA_14D"))

library(ggplot2)
ggplot(Injury_MG_res,aes(x=sample ,y=preds,color=sample)) + scale_y_continuous(expand=c(0,0)) + geom_boxplot(outlier.shape = NA,size=1)+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("Injury_MG_res.png",height=4,width=7) 


Injury_Rod_res$sample = factor(Injury_Rod_res$sample,levels=c("Control_2","NMDA_36hr","NMDA_54hr","NMDA_72hr","NMDA_96hr","NMDA_7D","NMDA_14D"))

library(ggplot2)
ggplot(Injury_Rod_res,aes(x=sample ,y=preds,color=sample)) + scale_y_continuous(expand=c(0,0)) + geom_boxplot(outlier.shape = NA,size=1)+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("Injury_Rod_res.png",height=4,width=7) 


Injury_Cone_res$sample = factor(Injury_Cone_res$sample,levels=c("Control_2","NMDA_36hr","NMDA_54hr","NMDA_72hr","NMDA_96hr","NMDA_7D","NMDA_14D"))

library(ggplot2)
ggplot(Injury_Cone_res,aes(x=sample ,y=preds,color=sample)) + scale_y_continuous(expand=c(0,0)) + geom_boxplot(outlier.shape = NA,size=1)+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")
ggsave("Injury_Cone_res.png",height=4,width=7) 

########
########
########





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

#####
##### Xenium5k genes #######
#####

dim(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$MG$Matrix)
dim(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$Rod$Matrix)
dim(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$Cone$Matrix)

Xenium5k_genes = rownames(Xenium5k_Mouse_Obj_Cl_Sp_clocks_Prepared$MG$Matrix)

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
save(Xenium5k_genes,file="Xenium5k_genes")

######
#############
######
###### retrain the Mouse models #####
######


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat_knn_pseudobulk <- readRDS("Mouse_seurat_knn_pseudobulk_2025")

Mouse_MG_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$MG,Need_genes = Xenium5k_genes)
save(Mouse_MG_model_Xen,file="Mouse_MG_model_Xen")
Mouse_Rod_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$Rod,Need_genes = Xenium5k_genes)
save(Mouse_Rod_model_Xen,file="Mouse_Rod_model_Xen")
Mouse_Cone_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$Cone,Need_genes = Xenium5k_genes)
save(Mouse_Cone_model_Xen,file="Mouse_Cone_model_Xen")
Mouse_AC_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$AC,Need_genes = Xenium5k_genes)
save(Mouse_AC_model_Xen,file="Mouse_AC_model_Xen")
Mouse_HC_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$HC,Need_genes = Xenium5k_genes)
save(Mouse_HC_model_Xen,file="Mouse_HC_model_Xen")
Mouse_BC_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$BC,Need_genes = Xenium5k_genes)
save(Mouse_BC_model_Xen,file="Mouse_BC_model_Xen")
Mouse_RGC_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$RGC,Need_genes = Xenium5k_genes)
save(Mouse_RGC_model_Xen,file="Mouse_RGC_model_Xen")
Mouse_RPE_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$RPE,Need_genes = Xenium5k_genes)
save(Mouse_RPE_model_Xen,file="Mouse_RPE_model_Xen")
Mouse_Microglia_model_Xen = Elastic_net_one_round(seu=Mouse_seurat_knn_pseudobulk$Microglia,Need_genes = Xenium5k_genes)
save(Mouse_Microglia_model_Xen,file="Mouse_Microglia_model_Xen")











###############
############### #################
############### 
############### #################
###############





