#####
##### ######
#####


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

#######-----------
####### for Zebrafish -------

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_pseudo_count <- readRDS("Zebrafish_pseudo_count_2025")

####### load the K-means ####
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
load(file="Zebrafish_DEGs_Plot_Kmeans_order")
load(file="Zebrafish_DEGs_Plot")

####### for Heatmap, order the clusters, young-middle-old #######
new_order2 = c(1,2,3,4,9,10,5,6,7,8)
Zebrafish_DEGs_Plot_Kmeans_order2 = Zebrafish_DEGs_Plot_Kmeans_order
Zebrafish_DEGs_Plot_Kmeans_order2$CT = sapply(strsplit(Zebrafish_DEGs_Plot_Kmeans_order2$genes,split="__"),function(x) x[[1]])
Zebrafish_DEGs_Plot_Kmeans_order2$G = sapply(strsplit(Zebrafish_DEGs_Plot_Kmeans_order2$genes,split="__"),function(x) x[[2]])

######
######

df_list <- lapply(names(Zebrafish_pseudo_count), function(ct) {
  # 取出稀疏矩阵
  mat_sparse <- Zebrafish_pseudo_count[[ct]]
  
  # 先转换成普通矩阵
  mat_dense <- as.matrix(mat_sparse)
  
  # 再转 data.frame
  df <- as.data.frame(mat_dense, check.names = FALSE)
  
  # 取基因名（假设行名是基因）
  genes <- rownames(mat_dense)
  
  # 构造新的行名：celltype__gene
  rownames(df) <- paste0(ct, "__", genes)
  
  df
})

combined_df <- do.call(rbind, df_list)

m = match(rownames(combined_df),Zebrafish_DEGs_Plot_Kmeans_order2$genes)

combined_df$cluster = Zebrafish_DEGs_Plot_Kmeans_order2$cluster[m]

######
index = c('young','young','young','young','middle','middle','old','old','old','old')
combined_df$class = index[combined_df$cluster]

######
Zebrafish_table = combined_df

#######
####### Next for Mouse ########
#######


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_pseudo_count <- readRDS("Mouse_pseudo_count_2025")


df_list <- lapply(names(Mouse_pseudo_count), function(ct) {
  # 取出稀疏矩阵
  mat_sparse <- Mouse_pseudo_count[[ct]]
  
  # 先转换成普通矩阵
  mat_dense <- as.matrix(mat_sparse)
  
  # 再转 data.frame
  df <- as.data.frame(mat_dense, check.names = FALSE)
  
  # 取基因名（假设行名是基因）
  genes <- rownames(mat_dense)
  
  # 构造新的行名：celltype__gene
  rownames(df) <- paste0(ct, "__", genes)
  
  df
})

combined_df <- do.call(rbind, df_list)

#####

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load(file="Mouse_DEGs_Plot_Kmeans_order")
load(file="Mouse_DEGs_Plot")

####### for Heatmap, order the clusters, young-middle-old #######
new_order2 = c(1,2,3,4,9,10,5,6,7,8)
Mouse_DEGs_Plot_Kmeans_order2 = Mouse_DEGs_Plot_Kmeans_order
Mouse_DEGs_Plot_Kmeans_order2$cluster = match(Mouse_DEGs_Plot_Kmeans_order$cluster,new_order2)

m = match(rownames(combined_df),Mouse_DEGs_Plot_Kmeans_order2$genes)

combined_df$cluster = Mouse_DEGs_Plot_Kmeans_order2$cluster[m]

######
index = c('young','young','young','young','middle','middle','old','old','old','old')
combined_df$class = index[combined_df$cluster]

######
Mouse_table = combined_df

####### Next for Human table: ######
#######

setwd("/zp1/data/share/Human_aging_new")

files = list.files()
files_Avg = files[grep("_Avg_2025",files)]
files_Avg_Name = gsub("Human_","",files_Avg)
files_Avg_Name = gsub("_Avg_2025","",files_Avg_Name)
files_Avg_Name = gsub("_clcl2","",files_Avg_Name)

load_object <- function(file) {
  # 检查文件存在
  if (!file.exists(file)) {
    stop("文件不存在：", file)
  }
  # 新建临时环境，不污染全局
  tmp <- new.env()
  # load 返回载入的对象名向量
  objs <- load(file, envir = tmp)
  
  # 如果只有一个对象，直接返回该对象
  if (length(objs) == 1) {
    return(tmp[[objs]])
  }
  # 否则返回一个列表，名字对应各个对象
  result <- mget(objs, envir = tmp)
  return(result)
}


Human_pseudo_count = list()
for(i in 1:length(files_Avg)){
    #########
    tmp = load_object(files_Avg[i])
    ####
    Human_pseudo_count <- c(Human_pseudo_count,list(tmp))
}

names(Human_pseudo_count) = files_Avg_Name

#######
#######

df_list <- lapply(names(Human_pseudo_count), function(ct) {
  # 取出稀疏矩阵
  mat_sparse <- Human_pseudo_count[[ct]][[1]]
  
  # 先转换成普通矩阵
  mat_dense <- as.matrix(mat_sparse)
  
  # 再转 data.frame
  df <- as.data.frame(mat_dense, check.names = FALSE)
  
  # 取基因名（假设行名是基因）
  genes <- rownames(mat_dense)
  
  # 构造新的行名：celltype__gene
  rownames(df) <- paste0(ct, "__", genes)
  df = df[,c("10","25","37.5","50","60","70","80","90")]
  df
})

combined_df <- do.call(rbind, df_list)

####
####

setwd("/zp1/data/share/Human_aging_new")
load(file="Human_DEGs_Plot_Kmeans_order")
load(file="Human_DEGs_Plot_cl")

new_order2 = c(1,2,3,4,10,11,12,5,6,7,8,9)
Human_DEGs_Plot_Kmeans_order2 = Human_DEGs_Plot_Kmeans_order
Human_DEGs_Plot_Kmeans_order2$cluster = match(Human_DEGs_Plot_Kmeans_order$cluster,new_order2)


m = match(rownames(combined_df),Human_DEGs_Plot_Kmeans_order2$genes)
combined_df$cluster = Human_DEGs_Plot_Kmeans_order2$cluster[m]

######
index = c('young','young','young','young','middle','middle','middle','old','old','old','old','old')
combined_df$class = index[combined_df$cluster]

######
Human_table = combined_df

#######
Output_table_list = list(Zebrafish=Zebrafish_table,Mouse=Mouse_table,Human=Human_table)

#######

library(openxlsx)
wb <- createWorkbook()
for (ct in names(Output_table_list)) {
  addWorksheet(wb, ct)
  writeData(wb, ct, Output_table_list[[ct]], rowNames = TRUE)
}
saveWorkbook(wb, "HMZ_DEGs_results.xlsx", overwrite = TRUE)

########
########
######## Output Overlap genes list from the same species ########
########
########


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

####

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

table_list = list(Zebrafish_Old=Z_up_Res[[1]],Zebrafish_Young=Z_down_Res[[1]],Mouse_Old=Mouse_up_Res[[1]],Mouse_Young=Mouse_down_Res[[1]],Human_Old=H_up_Res[[1]],Human_Young=H_down_Res[[1]])



library(openxlsx)
wb <- createWorkbook()
for (ct in names(table_list)) {
  addWorksheet(wb, ct)
  writeData(wb, ct, table_list[[ct]], rowNames = TRUE)
}
saveWorkbook(wb, "HMZ_DEGs_within_Species.xlsx", overwrite = TRUE)


######
######



















