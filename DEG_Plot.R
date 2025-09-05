########
######## figure2B #########
########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R


DEGs_to_BarPlot_F <- function(kc,Upclusters,Downclusters,Middleclusters){
    ###########
    kc$cell_type = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])
    ###########
    index = names(table(kc$cell_type))
    cell_type = c()
    num_DEGs = c()
    direction = c()
    ###########
    for(i in 1:length(index)){
        #######
        cell_type <- c(cell_type,rep(index[i],3))
        k1 = which(kc$cell_type == index[i] &  kc$cluster %in% Downclusters == T)
        k2 = which(kc$cell_type == index[i] &  kc$cluster %in% Upclusters == T)
        k3 = which(kc$cell_type == index[i] &  kc$cluster %in% Middleclusters == T)
        ####
        num_DEGs <- c(num_DEGs,length(k1),length(k2),length(k3))
        ####
        direction <- c(direction, c("Young", "Old","Middle"))
    }
    #####
    data = data.frame(cell_type=cell_type,num_DEGs=num_DEGs,direction=direction)
    data$direction = factor(data$direction,levels=c("Young","Middle","Old"))
    ######
    res = tapply(data$num_DEGs,data$cell_type,function(x) sum(abs(x)))
    cell_type_index = rev(names(res)[order(res)])
    ######
    data$cell_type = factor(data$cell_type,levels=cell_type_index)
    return(data)
}

####### for fish ######
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
load("Zebrafish_DEGs_Plot_Kmeans_order")

kc = Zebrafish_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8)
Downclusters = c(1,2,3,4)
Middleclusters = c(9,10)
Zebrafish_data = DEGs_to_BarPlot_F(kc,Upclusters,Downclusters,Middleclusters)

library(ggplot2)
ggplot(Zebrafish_data, aes(x = cell_type, y = num_DEGs, fill = direction)) +
  geom_col(position = "stack", width = 0.5) +                # 堆叠条形
  scale_fill_manual(values = c("Young" = "lightblue",
                               "Old"    = "pink",
                               "Middle" = "darkgreen")) +          # 指定颜色
  scale_y_continuous(expand = c(0, 0),limits=c(0,5500)) +                     # Y 轴从 0 开始，无扩展
  labs(x = "", y = "Number of DEGs") +
  theme_classic(base_size = 14) +
 theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),  # 给整图加黑框
    axis.title.y     = element_text(size = 16),                             # Y 轴标签放大
    axis.text.x      = element_text(angle = 30, hjust = 1, vjust = 1),      # X 轴标签旋转 30°
    legend.position  = "right"
  )

ggsave("Zebrafish_data_barplot.png",width=7,height=3)


####### Next mouse #####
setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load("Mouse_DEGs_Plot_Kmeans_order")

kc = Mouse_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8)
Downclusters = c(1,2,3,4)
Middleclusters = c(9,10)

Mouse_data = DEGs_to_BarPlot_F(kc,Upclusters,Downclusters,Middleclusters)

ggplot(Mouse_data, aes(x = cell_type, y = num_DEGs, fill = direction)) +
  geom_col(position = "stack", width = 0.5) +                # 堆叠条形
  scale_fill_manual(values = c("Young" = "lightblue",
                               "Old"    = "pink",
                               "Middle" = "darkgreen")) +          # 指定颜色
  scale_y_continuous(expand = c(0, 0),limits=c(0,5500)) +                     # Y 轴从 0 开始，无扩展
  labs(x = "", y = "Number of DEGs") +
  theme_classic(base_size = 14) +
 theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),  # 给整图加黑框
    axis.title.y     = element_text(size = 16),                             # Y 轴标签放大
    axis.text.x      = element_text(angle = 30, hjust = 1, vjust = 1),      # X 轴标签旋转 30°
    legend.position  = "right"
  )

ggsave("Mouse_data_barplot.png",width=7,height=3)

######
###### Next for human #######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######

setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot_Kmeans_order")


DEGs_to_BarPlot_F <- function(kc,Upclusters,Downclusters,Middleclusters){
    ###########
    kc$cell_type = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])
    kc$cell_type = sapply(strsplit(kc$cell_type,split="_"),function(x) x[[1]])
    ###########
    index = names(table(kc$cell_type))
    cell_type = c()
    num_DEGs = c()
    direction = c()
    ###########
    for(i in 1:length(index)){
        #######
        cell_type <- c(cell_type,rep(index[i],3))
        k1 = which(kc$cell_type == index[i] &  kc$cluster %in% Downclusters == T)
        k2 = which(kc$cell_type == index[i] &  kc$cluster %in% Upclusters == T)
        k3 = which(kc$cell_type == index[i] &  kc$cluster %in% Middleclusters == T)
        ####
        num_DEGs <- c(num_DEGs,length(k1),length(k2),length(k3))
        ####
        direction <- c(direction, c("Young", "Old","Middle"))
    }
    #####
    data = data.frame(cell_type=cell_type,num_DEGs=num_DEGs,direction=direction)
    data$direction = factor(data$direction,levels=c("Young","Middle","Old"))
    ######
    res = tapply(data$num_DEGs,data$cell_type,function(x) sum(abs(x)))
    cell_type_index = rev(names(res)[order(res)])
    ######
    data$cell_type = factor(data$cell_type,levels=cell_type_index)
    return(data)
}


kc = Human_DEGs_Plot_Kmeans_order[grep("_M__",Human_DEGs_Plot_Kmeans_order$genes),]
Upclusters = c(8,9,10,11,12)
Downclusters = c(1,2,3,4)
Middleclusters = c(5,6,7)
Human_M_data = DEGs_to_BarPlot_F(kc,Upclusters,Downclusters,Middleclusters)

ggplot(Human_M_data, aes(x = cell_type, y = num_DEGs, fill = direction)) +
  geom_col(position = "stack", width = 0.5) +                # 堆叠条形
  scale_fill_manual(values = c("Young" = "lightblue",
                               "Old"    = "pink",
                               "Middle" = "darkgreen")) +          # 指定颜色
  scale_y_continuous(expand = c(0, 0),limits=c(0,5500)) +                     # Y 轴从 0 开始，无扩展
  labs(x = "", y = "Number of DEGs") +
  theme_classic(base_size = 14) +
 theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),  # 给整图加黑框
    axis.title.y     = element_text(size = 16),                             # Y 轴标签放大
    axis.text.x      = element_text(angle = 30, hjust = 1, vjust = 1),      # X 轴标签旋转 30°
    legend.position  = "right"
  )
ggsave("Human_M_barplot.png",width=7,height=3)


                 

kc = Human_DEGs_Plot_Kmeans_order[grep("_F__",Human_DEGs_Plot_Kmeans_order$genes),]
Upclusters = c(8,9,10,11,12)
Downclusters = c(1,2,3,4)
Middleclusters = c(5,6,7)
Human_F_data = DEGs_to_BarPlot_F(kc,Upclusters,Downclusters,Middleclusters)

                 
ggplot(Human_F_data, aes(x = cell_type, y = num_DEGs, fill = direction)) +
  geom_col(position = "stack", width = 0.5) +                # 堆叠条形
  scale_fill_manual(values = c("Young" = "lightblue",
                               "Old"    = "pink",
                               "Middle" = "darkgreen")) +          # 指定颜色
  scale_y_continuous(expand = c(0, 0),limits=c(0,5500)) +                     # Y 轴从 0 开始，无扩展
  labs(x = "", y = "Number of DEGs") +
  theme_classic(base_size = 14) +
 theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),  # 给整图加黑框
    axis.title.y     = element_text(size = 16),                             # Y 轴标签放大
    axis.text.x      = element_text(angle = 30, hjust = 1, vjust = 1),      # X 轴标签旋转 30°
    legend.position  = "right"
  )
ggsave("Human_F_barplot.png",width=7,height=3)





#######
####### Next figure F G H #######----------------------------------------------------------------------------------------------------------------------------------------------------------------
#######

#######


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
load("Zebrafish_DEGs_Plot_Kmeans_order")

kc = Zebrafish_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8)
Downclusters = c(1,2,3,4)
#####
#####
kc$G = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])
kc$CT = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])

Z_up_kc = kc[which(kc$cluster %in% Upclusters == T),]
Z_down_kc = kc[which(kc$cluster %in% Downclusters == T),]


Z_up_Res = Get_overlap_tables_UP(Z_up_kc)
Z_down_Res = Get_overlap_tables_UP(Z_down_kc)


Old_Z = data.frame(Z_up_Res[[2]])
Young_Z = data.frame(Z_down_Res[[2]])

Old_Z = Old_Z[which(as.numeric(Old_Z$Var1) > 1),]
Young_Z = Young_Z[which(as.numeric(Young_Z$Var1) > 1),]

ggplot(Old_Z, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "pink") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 2000), expand = c(0, 0)) +
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
  scale_y_continuous(limits = c(0, 1500), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Zebrafish_young_overlap.png", height = 2.5, width = 4)


head(Z_up_Res[[1]],n=50)
head(Z_down_Res[[1]],n=50)
    

#### pfkfb1 ###

library(reshape2)
load("Zebrafish_DEGs_Plot")

Z_all_Matrix_cl_s_Plot = melt(Zebrafish_DEGs_Plot)

######
Gene_index = "pfkfb1"
######
Z_all_Matrix_cl_s_Plot_sub = Z_all_Matrix_cl_s_Plot[grep(Gene_index,Z_all_Matrix_cl_s_Plot$Var1),]
Z_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(Z_all_Matrix_cl_s_Plot_sub$Var1),split="__"),function(x) x[[1]])
Z_all_Matrix_cl_s_Plot_sub$time = as.numeric(Z_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######
Z_all_Matrix_cl_s_Plot_sub = Z_all_Matrix_cl_s_Plot_sub[which(Z_all_Matrix_cl_s_Plot_sub$CT != 'Microglia'),]

library(ggplot2)

ggplot(Z_all_Matrix_cl_s_Plot_sub, aes(x = time, y = value)) +
  geom_point(size = 3, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_continuous(breaks = c(6, 24,48)) +
  ylab("RNA z-score") +        # 设置 Y 轴标签
  xlab(NULL) +                 # 去掉 X 轴标签
  theme_classic() +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, size = 1),
    strip.text        = element_text(size = 14)  # 增大 facet 标签字体
  ) +
  facet_wrap(~ CT, ncol = 4, nrow = 2)

ggsave("Zebrafish_pfkfb1.png", height = 3.5, width = 6)



##### seprate G and CT ######


##### Next count the overlaps !!! #######
#####

                                       
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
Mouse_DOWN = kc[which(kc$Class=="DOWN"),]

Mouse_up_Res = Get_overlap_tables_UP(Mouse_UP)
Mouse_down_Res = Get_overlap_tables_UP(Mouse_DOWN)

##### 

Mouse_Res = list(Mouse_increasing= Mouse_up_Res[[3]],Mouse_decreasing=Mouse_down_Res[[3]])

#####
library(writexl)
write_xlsx(Mouse_Res, "Mouse_aging_overlaps_Mar2025.xlsx")



##### kc_index = Mouse_up_kc #######

######
######
######
Old_M = data.frame(Mouse_up_Res[[2]])
Young_M = data.frame(Mouse_down_Res[[2]])

Old_M = Old_M[which(as.numeric(Old_M$Var1) > 1),]
Young_M = Young_M[which(as.numeric(Young_M$Var1) > 1),]

ggplot(Old_M, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "pink") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 1500), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Mouse_old_overlap.png", height = 2.5, width = 4)


ggplot(Young_M, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 1500), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Mouse_young_overlap.png", height = 2.5, width = 4)

##############

##############


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load("Mouse_DEGs_Plot_Kmeans_order")

kc = Mouse_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8)
Downclusters = c(1,2,3,4)
#####
#####
kc$G = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])
kc$CT = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])

M_up_kc = kc[which(kc$cluster %in% Upclusters == T),]
M_down_kc = kc[which(kc$cluster %in% Downclusters == T),]


M_up_Res = Get_overlap_tables_UP(M_up_kc)
M_down_Res = Get_overlap_tables_UP(M_down_kc)


Old_M = data.frame(M_up_Res[[2]])
Young_M = data.frame(M_down_Res[[2]])

Old_M = Old_M[which(as.numeric(Old_M$Var1) > 3),]
Young_M = Young_M[which(as.numeric(Young_M$Var1) > 3),]

head(Mouse_up_Res[[1]])
head(Mouse_down_Res[[1]])


load("Mouse_DEGs_Plot")
######
Gene_index = "Stat1"
######

library(reshape2)
M_all_Matrix_cl_s_Plot = melt(Mouse_DEGs_Plot)
M_all_Matrix_cl_s_Plot_sub = M_all_Matrix_cl_s_Plot[grep(Gene_index,M_all_Matrix_cl_s_Plot$Var1),]
M_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(M_all_Matrix_cl_s_Plot_sub$Var1),split="__"),function(x) x[[1]])

M_all_Matrix_cl_s_Plot_sub$time = as.numeric(M_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######

                                       
library(ggplot2)

ggplot(M_all_Matrix_cl_s_Plot_sub, aes(x = time, y = value)) +
  geom_point(size = 3, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_continuous(breaks = c(12, 49,120)) +
  ylab("RNA z-score") +        # 设置 Y 轴标签
  xlab(NULL) +                 # 去掉 X 轴标签
  theme_classic() +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, size = 1),
    strip.text        = element_text(size = 14)  # 增大 facet 标签字体
  ) +
  facet_wrap(~ CT, ncol = 4, nrow = 2)

ggsave("Mouse_stat1.png", height = 3.5, width = 6)

                                       

#######
####### Next for Human #######
#######


setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot_Kmeans_order")

#######-----for M ######

kc = Human_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8,9)
Downclusters = c(1,2,3,4)

######
######
kc$CT = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])
kc$CT = sapply(strsplit(kc$CT,split="_"),function(x) x[[1]])
kc$G = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])


H_up_kc = kc[which(kc$cluster %in% Upclusters == T),]
H_down_kc = kc[which(kc$cluster %in% Downclusters == T),]

####
####

H_up_Res = Get_overlap_tables_UP(H_up_kc)
H_down_Res = Get_overlap_tables_UP(H_down_kc)

#####
#####



Old_H = data.frame(H_up_Res[[2]])
Young_H = data.frame(H_down_Res[[2]])

#####
#####

Old_H = Old_H[which(as.numeric(Old_H$Var1) > 1),]
Young_H = Young_H[which(as.numeric(Young_H$Var1) > 1),]


#####
#####
library(ggplot2)
      
ggplot(Old_H, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "pink") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 2500), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Human_old_overlap.png", height = 2.5, width = 4)
library(ggplot2)
      
ggplot(Young_H, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 3) +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_y_continuous(limits = c(0, 2500), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(size = 14)  # 增大 X 轴刻度标签字体大小
  )

ggsave("Human_young_overlap.png", height = 2.5, width = 4)





#######
head(H_up_Res[[1]],n=20)
head(H_down_Res[[1]],n=20)


load("Human_DEGs_Plot_cl")

H_all_Matrix_cl_s_Plot = melt(Human_DEGs_Plot_cl)

######
Gene_index = "MT1F"
######
H_all_Matrix_cl_s_Plot_sub = H_all_Matrix_cl_s_Plot[grep(Gene_index,H_all_Matrix_cl_s_Plot$Var1),]
H_all_Matrix_cl_s_Plot_sub$CT = sapply(strsplit(as.character(H_all_Matrix_cl_s_Plot_sub$Var1),split=":"),function(x) x[[1]])
H_all_Matrix_cl_s_Plot_sub$time = as.numeric(H_all_Matrix_cl_s_Plot_sub$Var2)
######
######
######
H_all_Matrix_cl_s_Plot_sub$CT = gsub("__MT1F","",H_all_Matrix_cl_s_Plot_sub$CT)
table(H_all_Matrix_cl_s_Plot_sub$CT)

H_all_Matrix_cl_s_Plot_sub = H_all_Matrix_cl_s_Plot_sub[-which(H_all_Matrix_cl_s_Plot_sub$CT == "Microglia_M"),]

                                            
library(ggplot2)

ggplot(H_all_Matrix_cl_s_Plot_sub, aes(x = time, y = value)) +
  geom_point(size = 3, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_continuous(breaks = c(10, 50,100)) +
  ylab("RNA z-score") +        # 设置 Y 轴标签
  xlab(NULL) +                 # 去掉 X 轴标签
  theme_classic() +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, size = 1),
    strip.text        = element_text(size = 14)  # 增大 facet 标签字体
  ) +
  facet_wrap(~ CT, ncol = 4, nrow = 2)

ggsave("H__MT1F.png", height = 3.5, width = 6)


                                       

######
######
######

library(ggplot2)
ggplot(Young_H, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 5) +
  theme_classic() + ylab("") + scale_y_continuous(limits=c(0,500))

ggsave("Human_young_overlap.png",height=4,width=3)


######
###### Done!!!! ########
######
