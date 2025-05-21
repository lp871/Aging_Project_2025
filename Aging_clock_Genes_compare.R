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
