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

#######-----for M ######

kc = Human_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8,9)
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
    ###
    df <- data.frame(
    H = c(1,0,0,1,1,0,1),
    M = c(0,1,0,1,0,1,1),
    Z = c(0,0,1,0,1,1,1),
    row.names = expr)
    #####
    library(ComplexHeatmap)
    png("test.png", width = 600, height = 550, res = 150)
    print(UpSet(
        comb_mat,
        # 保留原始 H, M, Z 顺序
        set_order  = c("H","M","Z"),
        # 按交集大小降序
        comb_order = order(-comb_size(comb_mat)),
        top_annotation = upset_top_annotation(
            comb_mat,
            # 把交集大小的柱子涂成 steelblue
            gp = gpar(fill = color),
            height = unit(5, "cm")
        ),
        right_annotation = upset_right_annotation(comb_mat)
    ))
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
###### Plot one gene #####
######

######----- for Human ###

setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot_cl")

library(reshape2)
H_all_Matrix_cl_s_Plot = melt(Human_DEGs_Plot_cl)

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
H_all_Matrix_cl_s_Plot_sub = H_all_Matrix_cl_s_Plot_sub[which(H_all_Matrix_cl_s_Plot_sub$CT == "RGC_F__CDKN2C"),]

library(ggplot2)
ggplot(H_all_Matrix_cl_s_Plot_sub,aes(x=time,y=value)) + geom_point(size=1,color="blue") + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_smooth(method = "lm", se = FALSE, color = "red") + scale_x_continuous(breaks=c(50,100))
ggsave("MG_F__SASH1.png",height=3,width=3)


######
######


######----- for Mouse ###

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load("Mouse_DEGs_Plot")

library(reshape2)
M_all_Matrix_cl_s_Plot = melt(Mouse_DEGs_Plot)

######
Gene_index = "REST"
Gene_index = "CDKN2C"
Gene_index = "SASH1"

Gene_index = "cdkn2c"

######
M_all_Matrix_cl_s_Plot_sub = M_all_Matrix_cl_s_Plot[grep(Gene_index,M_all_Matrix_cl_s_Plot$Var1),]
M_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(M_all_Matrix_cl_s_Plot_sub$Var1),split=":"),function(x) x[[1]])
M_all_Matrix_cl_s_Plot_sub$time = as.numeric(M_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######
M_all_Matrix_cl_s_Plot_sub = M_all_Matrix_cl_s_Plot_sub[which(M_all_Matrix_cl_s_Plot_sub$CT == "MG__Sash1"),]

library(ggplot2)
ggplot(M_all_Matrix_cl_s_Plot_sub,aes(x=time,y=value)) + geom_point(size=1,color="blue") + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_smooth(method = "lm", se = FALSE, color = "red") + scale_x_continuous(breaks=c(50,100))
ggsave("MG__Sash1.png",height=3,width=3)

#######

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
load("Zebrafish_DEGs_Plot")

library(reshape2)
Z_all_Matrix_cl_s_Plot = melt(Zebrafish_DEGs_Plot)

######
Gene_index = "REST"
Gene_index = "CDKN2C"
Gene_index = "SASH1"

Gene_index = "Sash1"
Gene_index = "cdkn2c"
Gene_index = "sash1a"

######
Z_all_Matrix_cl_s_Plot_sub = Z_all_Matrix_cl_s_Plot[grep(Gene_index,Z_all_Matrix_cl_s_Plot$Var1),]
Z_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(Z_all_Matrix_cl_s_Plot_sub$Var1),split=":"),function(x) x[[1]])
Z_all_Matrix_cl_s_Plot_sub$time = as.numeric(Z_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######
Z_all_Matrix_cl_s_Plot_sub = Z_all_Matrix_cl_s_Plot_sub[which(Z_all_Matrix_cl_s_Plot_sub$CT == "MG__sash1a"),]

library(ggplot2)
ggplot(Z_all_Matrix_cl_s_Plot_sub,aes(x=time,y=value)) + geom_point(size=1,color="blue") + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_smooth(method = "lm", se = FALSE, color = "red") + scale_x_continuous(breaks=c(50,100))
ggsave("MG__sash1a.png",height=3,width=3)


#######
###### print("Done!!!")
#######


