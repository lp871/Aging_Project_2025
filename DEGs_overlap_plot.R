###########
######
###########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

############
###### prepare DEGs list #######
############


######
###### ------- for Mouse --------- ######### ---------------------------------------------------------------------------------------------------------------------
######

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load("Mouse_DEGs_Plot_Kmeans_order")

kc = Mouse_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8)
Downclusters = c(1,2,3,4)

kc$CT = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])
kc$Gene = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])

kc$Class = "Unknown"
kc$Class[which(kc$cluster %in% Upclusters)] = "UP"
kc$Class[which(kc$cluster %in% Downclusters)] = "DOWN"

#####
Mouse_UP = kc[which(kc$Class=="UP"),]
Mouse_UP_list = split(Mouse_UP,Mouse_UP$CT)

Mouse_DOWN = kc[which(kc$Class=="DOWN"),]
Mouse_DOWN_list = split(Mouse_DOWN,Mouse_DOWN$CT)

##### for Zebrafish ########
#####

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
load("Zebrafish_DEGs_Plot_Kmeans_order")

kc = Zebrafish_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8)
Downclusters = c(1,2,3,4)

kc$CT = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])
kc$Gene = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])

kc$Class = "Unknown"
kc$Class[which(kc$cluster %in% Upclusters)] = "UP"
kc$Class[which(kc$cluster %in% Downclusters)] = "DOWN"

#####
Zebrafish_UP = kc[which(kc$Class=="UP"),]
Zebrafish_UP_list = split(Zebrafish_UP,Zebrafish_UP$CT)

Zebrafish_DOWN = kc[which(kc$Class=="DOWN"),]
Zebrafish_DOWN_list = split(Zebrafish_DOWN,Zebrafish_DOWN$CT)

######
###### for Human ######
######


setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot_Kmeans_order")


#######-----for Male------- ######

kc = Human_DEGs_Plot_Kmeans_order
Upclusters = c(8,9,10,11,12)
Downclusters = c(1,2,3,4)

#######

kc$CT = sapply(strsplit(kc$genes,split="_"),function(x) x[[1]])
kc$Gene = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])

kc$Class = "Unknown"
kc$Class[which(kc$cluster %in% Upclusters)] = "UP"
kc$Class[which(kc$cluster %in% Downclusters)] = "DOWN"

#####
Human_UP = kc[which(kc$Class=="UP"),]
Human_UP_list = split(Human_UP,Human_UP$CT)

Human_DOWN = kc[which(kc$Class=="DOWN"),]
Human_DOWN_list = split(Human_DOWN,Human_DOWN$CT)

#######
####### load the conseved gene pairs !!! ######
#######

setwd("/zp1/data/share/Human_aging_new")
HMZ_ortholog_combined <- readRDS("HMZ_ortholog_combined_2025")
names(HMZ_ortholog_combined)

#######

#######

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
    if(dim(HMZ_overlap)[1] > 0){
        HM_overlap = HM_overlap[-k1,]
        HZ_overlap = HZ_overlap[-k2,]
        MZ_overlap = MZ_overlap[-k3,]
    }
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

##########
##########

H = Human_UP_list$MG$Gene
M = Mouse_UP_list$MG$Gene
Z = Zebrafish_UP_list$MG$Gene

MG_3sp_overlap_UP = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="pink")


H = Human_DOWN_list$MG$Gene
M = Mouse_DOWN_list$MG$Gene
Z = Zebrafish_DOWN_list$MG$Gene

MG_3sp_overlap_DOWN = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="lightblue")


H = Human_DOWN_list$Rod$Gene
M = Mouse_DOWN_list$Rod$Gene
Z = Zebrafish_DOWN_list$Rod$Gene

Rod_3sp_overlap_DOWN = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="lightblue")

H = Human_UP_list$Rod$Gene
M = Mouse_UP_list$Rod$Gene
Z = Zebrafish_UP_list$Rod$Gene

Rod_3sp_overlap_UP = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="pink")

######

H = Human_DOWN_list$RGC$Gene
M = Mouse_DOWN_list$RGC$Gene
Z = Zebrafish_DOWN_list$RGC$Gene

RGC_3sp_overlap_DOWN = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="lightblue")

H = Human_UP_list$RGC$Gene
M = Mouse_UP_list$RGC$Gene
Z = Zebrafish_UP_list$RGC$Gene

RGC_3sp_overlap_UP = Compare_HMZ_pairs(H,M,Z,HMZ_ortholog_combined,color="pink")


######         
###### we will output the overlap results !!! ######
######
cell_types <- c("MG","Rod", "Cone", "AC", "HC", "RGC", "RPE", "BC")

# 用 lapply 批量生成 overlap 结果，结果存在一个列表里
overlap_up_list <- lapply(cell_types, function(ct) {
  H <- Human_UP_list[[ct]]$Gene
  M <- Mouse_UP_list[[ct]]$Gene
  Z <- Zebrafish_UP_list[[ct]]$Gene
  Compare_HMZ_pairs(H = H,
                    M = M,
                    Z = Z,
                    HMZ_ortholog_combined,
                    color = "pink")
})

# 给列表元素命名，方便后续调用
names(overlap_up_list) <- paste0(cell_types, "_3sp_overlap_Old")

overlap_down_list <- lapply(cell_types, function(ct) {
  H <- Human_DOWN_list[[ct]]$Gene
  M <- Mouse_DOWN_list[[ct]]$Gene
  Z <- Zebrafish_DOWN_list[[ct]]$Gene
  Compare_HMZ_pairs(H = H,
                    M = M,
                    Z = Z,
                    HMZ_ortholog_combined,
                    color = "pink")
})

names(overlap_down_list) <- paste0(cell_types, "_3sp_overlap_Young")

##### total list ####
overlap_DEGs_total = c(overlap_up_list,overlap_down_list)

overlap_DEGs_total_toExcel = list()
for(i in 1:length(overlap_DEGs_total)){
    print(i)
    tmp = overlap_DEGs_total[[i]]
    tmp_1 = tmp$HM
    tmp_2 = tmp$HZ
    tmp_3 = tmp$MZ
    tmp_4 = tmp$HMZ
    ###
    tmp_1 = data.frame(class="HM",gene=tmp_1$"HM_index")
    tmp_2 = data.frame(class="HZ",gene=tmp_2$"HZ_index")
    tmp_3 = data.frame(class="MZ",gene=tmp_3$"MZ_index")
    if(dim(tmp_4)[1] > 0){
        tmp_4 = data.frame(class="HMZ",gene=tmp_4$"HMZ_index")
        overlap_DEGs_total_toExcel <- c(overlap_DEGs_total_toExcel,list(rbind(tmp_1,tmp_2,tmp_3,tmp_4)))
    }else{
        overlap_DEGs_total_toExcel <- c(overlap_DEGs_total_toExcel,list(rbind(tmp_1,tmp_2,tmp_3)))
    }
    ####
}
names(overlap_DEGs_total_toExcel) = names(overlap_DEGs_total)
#####
library(writexl)
write_xlsx(overlap_DEGs_total_toExcel, path = "TableS2:Overlap_DEGs_between_HMZ.xlsx")

                 









                 
######
###### Plot one gene #####
######

######----- for Human ###

setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot")

library(reshape2)
H_all_Matrix_cl_s_Plot = melt(Human_DEGs_Plot)

######
Gene_index = "REST"
Gene_index = "CDKN2C"
Gene_index = "SASH1"

######
H_all_Matrix_cl_s_Plot_sub = H_all_Matrix_cl_s_Plot[grep(Gene_index,H_all_Matrix_cl_s_Plot$Var1),]
H_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(H_all_Matrix_cl_s_Plot_sub$Var1),split=":"),function(x) x[[1]])
H_all_Matrix_cl_s_Plot_sub$time = as.numeric(H_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######
H_all_Matrix_cl_s_Plot_sub = H_all_Matrix_cl_s_Plot_sub[which(H_all_Matrix_cl_s_Plot_sub$CT == "MG_F__SASH1"),]
H_all_Matrix_cl_s_Plot_sub = H_all_Matrix_cl_s_Plot_sub[which(H_all_Matrix_cl_s_Plot_sub$CT == "Rod_F__CDKN2C"),]
H_all_Matrix_cl_s_Plot_sub = H_all_Matrix_cl_s_Plot_sub[which(H_all_Matrix_cl_s_Plot_sub$CT == "RGC_F__REST"),]

library(ggplot2)
ggplot(H_all_Matrix_cl_s_Plot_sub,aes(x=time,y=value)) + geom_point(size=3,color="blue") + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_smooth(method = "lm", se = FALSE, color = "red") + scale_x_continuous(breaks=c(10,50,100))
ggsave("RGC_F__REST.png",height=2.5,width=2.5)


######
######


######----- for Mouse ###

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load("Mouse_DEGs_Plot")

library(reshape2)
M_all_Matrix_cl_s_Plot = melt(Mouse_DEGs_Plot)

######
Gene_index = "Rest"
Gene_index = "CDKN2C"
Gene_index = "SASH1"

Gene_index = "Cdkn2c"
Gene_index = "Sash1"


######
M_all_Matrix_cl_s_Plot_sub = M_all_Matrix_cl_s_Plot[grep(Gene_index,M_all_Matrix_cl_s_Plot$Var1),]
M_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(M_all_Matrix_cl_s_Plot_sub$Var1),split=":"),function(x) x[[1]])
M_all_Matrix_cl_s_Plot_sub$time = as.numeric(M_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######
M_all_Matrix_cl_s_Plot_sub = M_all_Matrix_cl_s_Plot_sub[which(M_all_Matrix_cl_s_Plot_sub$CT == "MG__Sash1"),]

library(ggplot2)
ggplot(M_all_Matrix_cl_s_Plot_sub,aes(x=time,y=value)) + geom_point(size=3,color="blue") + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_smooth(method = "lm", se = FALSE, color = "red") + scale_x_continuous(breaks=c(12,49,120))
ggsave("MG__Sash1.png",height=2.5,width=2.5)

#######

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
load("Zebrafish_DEGs_Plot")

library(reshape2)
Z_all_Matrix_cl_s_Plot = melt(Zebrafish_DEGs_Plot)

######
Gene_index = "REST"
Gene_index = "CDKN2C"
Gene_index = "SASH1"

Gene_index = "sash1"
Gene_index = "cdkn2c"
Gene_index = "sash1a"

######
Z_all_Matrix_cl_s_Plot_sub = Z_all_Matrix_cl_s_Plot[grep(Gene_index,Z_all_Matrix_cl_s_Plot$Var1),]
Z_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(Z_all_Matrix_cl_s_Plot_sub$Var1),split=":"),function(x) x[[1]])
Z_all_Matrix_cl_s_Plot_sub$time = as.numeric(Z_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######
Z_all_Matrix_cl_s_Plot_sub = Z_all_Matrix_cl_s_Plot_sub[which(Z_all_Matrix_cl_s_Plot_sub$CT == "RGC__rest"),]

library(ggplot2)
ggplot(Z_all_Matrix_cl_s_Plot_sub,aes(x=time,y=value)) + geom_point(size=3,color="blue") + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_smooth(method = "lm", se = FALSE, color = "red") + scale_x_continuous(breaks=c(6,24,48))
ggsave("RGC_rest.png",height=2.5,width=2.5)


#######
###### print("Done!!!")
#######


