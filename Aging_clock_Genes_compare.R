#######
####### 
####### 读入 HMZ 的 aging model #######
#######
####### 查找 overlap genes #######
#######




######## M 的 old gene #######
########

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")

load("Mouse_MG_model")
load("Mouse_RGC_model")
load("Mouse_AC_model")
load("Mouse_HC_model")
load("Mouse_Rod_model")
load("Mouse_Cone_model")
load("Mouse_BC_model")
load("Mouse_RPE_model")
load("Mouse_Microglia_model")


##########
######

MG = Get_old_G_from_model(Mouse_MG_model$model)
RGC = Get_old_G_from_model(Mouse_RGC_model$model)
AC = Get_old_G_from_model(Mouse_AC_model$model)
HC = Get_old_G_from_model(Mouse_HC_model$model)
BC = Get_old_G_from_model(Mouse_BC_model$model)
RPE = Get_old_G_from_model(Mouse_RPE_model$model)
Rod = Get_old_G_from_model(Mouse_Rod_model$model)
Cone = Get_old_G_from_model(Mouse_Cone_model$model)
Microglia = Get_old_G_from_model(Mouse_Microglia_model$model)

G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

Get_old_G_from_model <- function(fit_res){
    library(glmnet)
    coef_df  = coef(fit_res, s = fit_res$lambda.min)
    #####
    coef_df  <- data.frame(
    gene        = rownames(coef_df),
    coefficient = coef_df[, 1],
    row.names   = NULL,
    stringsAsFactors = FALSE
    )
    GetPositiveGenes <- function(coef_df) {
    # 过滤出系数大于 0 的行，并排除截距
        pos <- coef_df$coefficient > 0 & coef_df$gene != '(Intercept)'
        genes_pos <- coef_df$gene[pos]
        return(genes_pos)
    }
    ######
    res= GetPositiveGenes(coef_df)
    res
}

Get_young_G_from_model <- function(fit_res){
    library(glmnet)
    coef_df  = coef(fit_res, s = fit_res$lambda.min)
    #####
    coef_df  <- data.frame(
    gene        = rownames(coef_df),
    coefficient = coef_df[, 1],
    row.names   = NULL,
    stringsAsFactors = FALSE
    )
    GetPositiveGenes <- function(coef_df) {
    # 过滤出系数大于 0 的行，并排除截距
        pos <- coef_df$coefficient < 0 & coef_df$gene != '(Intercept)'
        genes_pos <- coef_df$gene[pos]
        return(genes_pos)
    }
    ######
    res= GetPositiveGenes(coef_df)
    res
}





kc_index = data.frame(CT=CT,G=G)

Mouse_old_res = Get_overlap_tables_UP(kc_index)
Old_Z = data.frame(Mouse_old_res[[2]])
Old_Z = Old_Z[which(as.numeric(Old_Z$Var1) > 1),]



ggplot(Old_Z, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "pink") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 1600), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Mouse_old_overlap.png", height = 2.5, width = 4)


MG = Get_young_G_from_model(Mouse_MG_model$model)
RGC = Get_young_G_from_model(Mouse_RGC_model$model)
AC = Get_young_G_from_model(Mouse_AC_model$model)
HC = Get_young_G_from_model(Mouse_HC_model$model)
BC = Get_young_G_from_model(Mouse_BC_model$model)
RPE = Get_young_G_from_model(Mouse_RPE_model$model)
Rod = Get_young_G_from_model(Mouse_Rod_model$model)
Cone = Get_young_G_from_model(Mouse_Cone_model$model)
Microglia = Get_young_G_from_model(Mouse_Microglia_model$model)

G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

kc_index = data.frame(CT=CT,G=G)
Mouse_young_res = Get_overlap_tables_UP(kc_index)
Young_Z = data.frame(Mouse_young_res[[2]])
Young_Z = Young_Z[which(as.numeric(Young_Z$Var1) > 1),]

head(Mouse_young_res[[1]],n=20)
head(Mouse_old_res[[1]],n=20)

ggplot(Young_Z, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 1600), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Mouse_young_overlap.png", height = 2.5, width = 4)


Get_overlap_tables_UP <- function(kc_index){
    ######
    all_ct = names(table(kc_index$CT))
    all_G = names(table(kc_index$G))
    ######
    tmp_mat = matrix(0,nrow=length(all_G),ncol=length(all_ct))
    rownames(tmp_mat) <- all_G
    colnames(tmp_mat) <- all_ct
    ######
    for(i in 1:dim(tmp_mat)[1]){
        for(j in 1:dim(tmp_mat)[2]){
            tmp_ct = colnames(tmp_mat)[j]
            tmp_g = rownames(tmp_mat)[i]
            ####
            k = which(kc_index$CT == tmp_ct & kc_index$G == tmp_g)
            if(length(k) > 0){
                tmp_mat[i,j] = 1
            }
        }
    }
    ######
    tmp_mat = data.frame(tmp_mat)
    tmp_mat$sum = apply(tmp_mat,1,sum)
    tmp_mat = tmp_mat[order(tmp_mat$sum,decreasing=T),]
    #######
    #######
    #######
    counts = data.frame(table(tmp_mat$sum))
    #######
    merge_table = data.frame(tmp_mat)
    merge_table$cluster = "NotFind"
    #######
    for(i in 1:length(rownames(merge_table))){
        ######
        tmp_G = rownames(merge_table)[i]
        tmp_kc_index = kc_index[which(kc_index$G == tmp_G),]
        tmp_kc_index$cluster2 = paste(tmp_kc_index$CT,tmp_kc_index$cluster,sep="_")
        res = paste(tmp_kc_index$cluster2,collapse=";")
        merge_table$cluster[i] = res
    }
    #######
    merge_table = cbind(data.frame(gene = rownames(merge_table)),merge_table)
    #######
    reslist = list(table=tmp_mat,counts = counts,merge_table = merge_table)
    #######
    return(reslist)
}







setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")

load("Zebrafish_MG_model")
load("Zebrafish_RGC_model")
load("Zebrafish_AC_model")
load("Zebrafish_HC_model")
load("Zebrafish_Rod_model")
load("Zebrafish_Cone_model")
load("Zebrafish_BC_model")
load("Zebrafish_RPE_model")
load("Zebrafish_Microglia_model")


MG = Get_old_G_from_model(Zebrafish_MG_model$model)
RGC = Get_old_G_from_model(Zebrafish_RGC_model$model)
AC = Get_old_G_from_model(Zebrafish_AC_model$model)
HC = Get_old_G_from_model(Zebrafish_HC_model$model)
BC = Get_old_G_from_model(Zebrafish_BC_model$model)
RPE = Get_old_G_from_model(Zebrafish_RPE_model$model)
Rod = Get_old_G_from_model(Zebrafish_Rod_model$model)
Cone = Get_old_G_from_model(Zebrafish_Cone_model$model)
Microglia = Get_old_G_from_model(Zebrafish_Microglia_model$model)

MG = Get_young_G_from_model(Zebrafish_MG_model$model)
RGC = Get_young_G_from_model(Zebrafish_RGC_model$model)
AC = Get_young_G_from_model(Zebrafish_AC_model$model)
HC = Get_young_G_from_model(Zebrafish_HC_model$model)
BC = Get_young_G_from_model(Zebrafish_BC_model$model)
RPE = Get_young_G_from_model(Zebrafish_RPE_model$model)
Rod = Get_young_G_from_model(Zebrafish_Rod_model$model)
Cone = Get_young_G_from_model(Zebrafish_Cone_model$model)
Microglia = Get_young_G_from_model(Zebrafish_Microglia_model$model)



G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))


kc_index = data.frame(CT=CT,G=G)

Zebrafish_old_res = Get_overlap_tables_UP(kc_index)
Old_Z = data.frame(Zebrafish_old_res[[2]])
Old_Z = Old_Z[which(as.numeric(Old_Z$Var1) > 1),]

head(Zebrafish_old_res[[1]],n=20)

Zebrafish_young_res = Get_overlap_tables_UP(kc_index)
Young_Z = data.frame(Zebrafish_young_res[[2]])
Young_Z = Young_Z[which(as.numeric(Young_Z$Var1) > 1),]



ggplot(Old_Z, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "pink") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 1400), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Zebrafish_old_overlap.png", height = 2.5, width = 4)


ggplot(Young_Z, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 1400), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Zebrafish_young_overlap.png", height = 2.5, width = 4)



########
########

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

MG = c(Get_old_G_from_model(Human_MG_model_F$model),Get_old_G_from_model(Human_MG_model_M$model))
RGC = c(Get_old_G_from_model(Human_RGC_model_F$model),Get_old_G_from_model(Human_RGC_model_M$model))
AC = c(Get_old_G_from_model(Human_AC_model_F$model),Get_old_G_from_model(Human_AC_model_M$model))
HC = c(Get_old_G_from_model(Human_HC_model_F$model),Get_old_G_from_model(Human_HC_model_M$model))
BC = c(Get_old_G_from_model(Human_BC_model_F$model),Get_old_G_from_model(Human_BC_model_M$model))
RPE = c(Get_old_G_from_model(Human_RPE_model_F$model),Get_old_G_from_model(Human_RPE_model_M$model))
Rod = c(Get_old_G_from_model(Human_Rod_model_F$model),Get_old_G_from_model(Human_Rod_model_M$model))
Cone = c(Get_old_G_from_model(Human_Cone_model_F$model),Get_old_G_from_model(Human_Cone_model_M$model))
Microglia = c(Get_old_G_from_model(Human_Microglia_model_F$model),Get_old_G_from_model(Human_Microglia_model_M$model))



MG = c(Get_young_G_from_model(Human_MG_model_F$model),Get_young_G_from_model(Human_MG_model_M$model))
RGC = c(Get_young_G_from_model(Human_RGC_model_F$model),Get_young_G_from_model(Human_RGC_model_M$model))
AC = c(Get_young_G_from_model(Human_AC_model_F$model),Get_young_G_from_model(Human_AC_model_M$model))
HC = c(Get_young_G_from_model(Human_HC_model_F$model),Get_young_G_from_model(Human_HC_model_M$model))
BC = c(Get_young_G_from_model(Human_BC_model_F$model),Get_young_G_from_model(Human_BC_model_M$model))
RPE = c(Get_young_G_from_model(Human_RPE_model_F$model),Get_young_G_from_model(Human_RPE_model_M$model))
Rod = c(Get_young_G_from_model(Human_Rod_model_F$model),Get_young_G_from_model(Human_Rod_model_M$model))
Cone = c(Get_young_G_from_model(Human_Cone_model_F$model),Get_young_G_from_model(Human_Cone_model_M$model))
Microglia = c(Get_young_G_from_model(Human_Microglia_model_F$model),Get_young_G_from_model(Human_Microglia_model_M$model))



G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

kc_index = data.frame(CT=CT,G=G)

Human_old_res = Get_overlap_tables_UP(kc_index)
Old_Z = data.frame(Human_old_res[[2]])
Old_Z = Old_Z[which(as.numeric(Old_Z$Var1) > 1),]


Human_young_res = Get_overlap_tables_UP(kc_index)
Young_Z = data.frame(Human_young_res[[2]])
Young_Z = Young_Z[which(as.numeric(Young_Z$Var1) > 1),]



ggplot(Old_Z, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "pink") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 3000), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Human_old_overlap.png", height = 2.5, width = 4)


ggplot(Young_Z, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 3000), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Human_young_overlap.png", height = 2.5, width = 4)




head(Human_old_res[[1]],n=10)
head(Human_young_res[[1]],n=10)

##########
##########
##########

##########-------- COMPARE BETWEEN HMZ ---------########
##########


setwd("/zp1/data/share/Human_aging_new")
HMZ_ortholog_combined <- readRDS("HMZ_ortholog_combined_2025")
names(HMZ_ortholog_combined)



Compare_HMZ_pairs <- function(H,M,Z,HMZ_ortholog_combined,color="pink"){
    H = H[!duplicated(H)]
    M = M[!duplicated(M)]
    Z = Z[!duplicated(Z)]
    #####
    total_genes_H = c(HMZ_ortholog_combined$HM$human,HMZ_ortholog_combined$HZ$human,HMZ_ortholog_combined$HMZ$human)
    total_genes_M = c(HMZ_ortholog_combined$HM$mouse,HMZ_ortholog_combined$MZ$mouse,HMZ_ortholog_combined$HMZ$mouse)
    total_genes_Z = c(HMZ_ortholog_combined$MZ$zebrafish,HMZ_ortholog_combined$HZ$zebrafish,HMZ_ortholog_combined$HMZ$zebrafish)
    #####
    #####
    H_cl = H[which(H %in% total_genes_H == T)]
    M_cl = M[which(M %in% total_genes_M == T)]
    Z_cl = Z[which(Z %in% total_genes_Z == T)]
    #####
    #####
    HM_overlap = HMZ_ortholog_combined$HM[which(HMZ_ortholog_combined$HM$human %in% H_cl == T & HMZ_ortholog_combined$HM$mouse %in% M_cl == T),]
    HZ_overlap = HMZ_ortholog_combined$HZ[which(HMZ_ortholog_combined$HZ$human %in% H_cl == T & HMZ_ortholog_combined$HZ$zebrafish %in% Z_cl == T),]
    MZ_overlap = HMZ_ortholog_combined$MZ[which(HMZ_ortholog_combined$MZ$mouse %in% M_cl == T & HMZ_ortholog_combined$MZ$zebrafish %in% Z_cl == T),]
    ######
    HMZ_overlap = HMZ_ortholog_combined$HMZ[which(HMZ_ortholog_combined$HMZ$human %in% H_cl == T & HMZ_ortholog_combined$HMZ$mouse %in% M_cl == T & HMZ_ortholog_combined$HMZ$zebrafish %in% Z_cl == T),]
    ######
    ###### remove duplicates both in HMZ_overlap and HM_overlap HZ_overlap and MZ_overlap ############
    ######
    HMZ_overlap_index_HM = paste(HMZ_overlap$human,HMZ_overlap$mouse,sep="__")
    HMZ_overlap_index_HZ = paste(HMZ_overlap$human,HMZ_overlap$zebrafish,sep="__")
    HMZ_overlap_index_MZ = paste(HMZ_overlap$mouse,HMZ_overlap$zebrafish,sep="__")
    ######
    k1 = which(HM_overlap$HM_index %in% HMZ_overlap_index_HM == T)
    k2 = which(HZ_overlap$HZ_index %in% HMZ_overlap_index_HZ == T)
    k3 = which(MZ_overlap$MZ_index %in% HMZ_overlap_index_MZ == T)
    ######
    HM_overlap = HM_overlap[-k1,]
    HZ_overlap = HZ_overlap[-k2,]
    MZ_overlap = MZ_overlap[-k3,]
    ######
    ######
    ######
    H_overlap_Gs = c(HM_overlap$human,HZ_overlap$human,HMZ_overlap$human)
    M_overlap_Gs = c(HZ_overlap$mouse,MZ_overlap$mouse,HMZ_overlap$mouse)
    Z_overlap_Gs = c(MZ_overlap$zebrafish,HZ_overlap$zebrafish,HMZ_overlap$zebrafish)
    #####
    H_sp = H_cl[which(H_cl %in% H_overlap_Gs == F)]
    M_sp = M_cl[which(M_cl %in% M_overlap_Gs == F)]
    Z_sp = Z_cl[which(Z_cl %in% Z_overlap_Gs == F)]
    #####
    #####
    library(UpSetR)
    counts <- c(
    `100` = length(unique(H_sp)),   # H only
    `010` = length(unique(M_sp)),   # M only
    `001` = length(unique(Z_sp)),   # Z only
    `110` = length(unique(HM_overlap$HM_index)),   # H & M
    `101` = length(unique(HZ_overlap$HZ_index)),   # H & Z
    `011` = length(unique(MZ_overlap$MZ_index)),   # M & Z
    `111` = length(unique(HMZ_overlap$HMZ_index))    # H & M & Z
    )
    ####
    # counts 向量名就是二进制 code：1=选，0=不选，三位分别对应 H, M, Z
    dummy_ids <- lapply(names(counts), function(code) {
    # 生成 code_1, code_2, … code_n
    paste0(code, "_", seq_len(counts[code]))
    })
    names(dummy_ids) <- names(counts)
    all_ids <- unlist(dummy_ids, use.names = FALSE)

    # membership逻辑：如果 code 第1位是1，属于 H；否则不属于
    H_ids <- unlist(dummy_ids[ substr(names(dummy_ids),1,1) == "1" ])
    M_ids <- unlist(dummy_ids[ substr(names(dummy_ids),2,2) == "1" ])
    Z_ids <- unlist(dummy_ids[ substr(names(dummy_ids),3,3) == "1" ])
    set_list <- list(H = H_ids, M = M_ids, Z = Z_ids)
    ####
    ####
    library(ComplexHeatmap)

    # 生成组合矩阵
    comb_mat <- make_comb_mat(set_list)
    #####
    #####
    library(ComplexHeatmap)
    pink_fill <- gpar(fill = color)

png("test_with_labels.png", width = 1600, height = 1250, res = 350)

up <- UpSet(
  comb_mat,
  set_order  = c("Z","M","H"),
  comb_order = order(-comb_size(comb_mat)),
  top_annotation = upset_top_annotation(
    comb_mat,
    gp     = gpar(fill = color),
    height = unit(5, "cm")
  ),
  right_annotation = upset_right_annotation(comb_mat,gp     = gpar(fill = color))
)

draw(up)

# 对 comb_size 先按 comb_order 排序
ord  <- order(-comb_size(comb_mat))
vals <- comb_size(comb_mat)[ord]

decorate_annotation("intersection_size", slice = 1, {
  for(i in seq_along(vals)) {
    grid.text(
      label         = vals[i],
      x             = unit(i, "native"),                        # 现在 i 对应绘图中的第 i 根 bar
      y             = unit(vals[i], "native") + unit(1, "mm"),  # 顶部再上移 1mm
      just          = "bottom",
      default.units = "native",
      gp            = gpar(fontsize = 10)
    )
  }
})

dev.off()

    #####
    res_list = list(H_sp=H_sp,M_sp=M_sp,Z_sp=Z_sp,HM=HM_overlap,HZ=HZ_overlap,MZ=MZ_overlap,HMZ=HMZ_overlap)
    #####
    return(res_list)
}


H = c(Get_old_G_from_model(Human_MG_model_F$model),Get_old_G_from_model(Human_MG_model_M$model))
M = Get_old_G_from_model(Mouse_MG_model$model)
Z = Get_old_G_from_model(Zebrafish_MG_model$model)

MG_3sp_overlap_UP = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="pink")


H = c(Get_young_G_from_model(Human_MG_model_F$model),Get_young_G_from_model(Human_MG_model_M$model))
M = Get_young_G_from_model(Mouse_MG_model$model)
Z = Get_young_G_from_model(Zebrafish_MG_model$model)

MG_3sp_overlap_DOWN = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="lightblue")











