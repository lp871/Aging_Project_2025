########
######## for Zebrafish aging DEGs ##########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_pseudo_count <- readRDS("Zebrafish_pseudo_count_2025")

names(Zebrafish_pseudo_count)

########
Zebrafish_pseudo_count_input = list()

for(i in 1:length(Zebrafish_pseudo_count)){
    tmp_mat = Zebrafish_pseudo_count[[i]]
    tmp_table = data.frame(sample=colnames(tmp_mat),age=colnames(tmp_mat))
    ####
    tmp_list = list(bulk_mat=tmp_mat,bulk_meta=tmp_table)
    Zebrafish_pseudo_count_input <- c(Zebrafish_pseudo_count_input,list(tmp_list))
    ####
}

names(Zebrafish_pseudo_count_input) = names(Zebrafish_pseudo_count)

######## rm 28mo there!!! #######
for(i in 1:length(Zebrafish_pseudo_count_input)){
    ###
    tmp = Zebrafish_pseudo_count_input[[i]]
    ###
    tmp_1 = tmp[[1]]
    tmp_2 = tmp[[2]]
    ####
    k = which(colnames(tmp_1) == "28")
    tmp_1 =tmp_1[,-k]
    k2 = which(tmp_2$sample == "28")
    tmp_2 = tmp_2[-k2,]
    ####
    tmp[[1]] = tmp_1
    tmp[[2]] = tmp_2
    #####
    Zebrafish_pseudo_count_input[[i]] = tmp
}



########

Smooth_FUN <- function(x,age){
    ##########
    vector = as.numeric(x)
    vector = scale(vector)
    ##########
    if(length(unique(vector)) == 1){
        df = data.frame(age=age,value=vector,predict=vector,p_value=1)
        return(df)
    }
    ##########
    df = data.frame(age=age,value=vector)
    df = df[order(df$age),]
    rownames(df) = 1:dim(df)[1]
    #######
    library(mgcv)
    model_select <- gam(value ~ s(age, k = 4), data = df, select = TRUE)
    p_value <- summary(model_select)$s.table[, 4]
    #######
    df_predict = predict(model_select, newdata = df, type = "response")
    df$predict = df_predict
    #######
    df$p_value = p_value
    ####### summary(model_select）#####
    return(df)
}

#### list_process2 = Zebrafish_pseudo_count_input[[3]]

Identify_DEGs <- function(list_process2){
    #####
    ########## 首先根据 bulk 的 data, 我们计算一下 correlation ############
    #####
    bulk_mat = list_process2$bulk_mat
    bulk_meta = list_process2$bulk_meta
    all.equal(colnames(bulk_mat),bulk_meta$sample)
    ########## filter bulk_mat ####
    avg_by_gene = apply(bulk_mat,1,mean)
    k = which(avg_by_gene > 0.01)
    bulk_mat = bulk_mat[k,]
    ##########
    Cor_table = data.frame(Gene=rownames(bulk_mat),Cor=0)
    ########## Cor_table[which(Cor_table$Gene == "Apoe"),]
    ##########
    library(parallel)
    ncores <- 30
    #### Smooth <- apply(bulk_mat,1,Smooth_FUN,age=as.numeric(bulk_meta$age))
    res_list <- mclapply(
        seq_len(nrow(bulk_mat)),
        function(i) Smooth_FUN(bulk_mat[i, ], age = as.numeric(bulk_meta$age)),
        mc.cores = ncores
    )
    ########## cutoff 0.3 #########
    ########## Smooth$Aqp4 ########
    #delta = c()
    pvalue = c()
    for(i in 1:length(res_list)){
        print(i)
        tmp = res_list[[i]]
        pvalue = c(pvalue,tmp$p_value[1])
    }
    #####
    #Cor_table$delta = delta
    Cor_table$pvalue = pvalue
    res_table = Cor_table
    #####
    res_table$TAG = "NotDEGs"
    ##### res_table[which(res_table$Gene == "stat1"),]
    print(dim(res_table))
    print(which(res_table$pvalue < 0.05))
    #####
    res_table_cl = res_table[which(res_table$pvalue < 0.05),]
    #####
    print("Done!")
    #####
    return(res_table_cl)
}


#####
Zebrafish_DEGs <- list()

for(i in 1:length(Zebrafish_pseudo_count_input)){
    print(i)
    #####
    res = Identify_DEGs(Zebrafish_pseudo_count_input[[i]])
    #####
    Zebrafish_DEGs <- c(Zebrafish_DEGs,list(res))
}

######
names(Zebrafish_DEGs) <- names(Zebrafish_pseudo_count_input)

######

######
DEGs_List = Zebrafish_DEGs
Avgs_List = Zebrafish_pseudo_count_input
col_index = as.character(c(1,3,6,12,18,22,24,30,36,48))

Merge_DEGs_genes_and_Avg_expression <- function(DEGs_List,Avgs_List,col_index){
    #########
    #########
    Avgs_List2  = list()
    for(i in 1:length(Avgs_List)){
        ####
        tmp = Avgs_List[[i]][[1]]
        rownames(tmp) = paste(names(Avgs_List)[i],rownames(tmp),sep="__")
        print(dim(tmp))
        Avgs_List2[[i]] = tmp
    }
    #########
    all_Matrix = do.call("rbind",Avgs_List2)
    #########
    for(i in 1:length(DEGs_List)){
        ####
        tmp = DEGs_List[[i]]
        tmp$Gene = paste(names(DEGs_List)[i],tmp$Gene,sep="__")
        DEGs_List[[i]] = tmp
    }
    #########
    all_DEGs = do.call("rbind",DEGs_List)
    ######### all_Avg[grep("Clu",all_Avg$Gene),] #########
    ######### all_Avg[grep("Xkr4",all_Avg$Gene),]
    all_Matrix_cl = all_Matrix[which(rownames(all_Matrix) %in% all_DEGs$Gene == T),]
    #########
    all_Matrix_cl = all_Matrix_cl[,col_index]
    #########
    all_Matrix_cl_s = t(apply(all_Matrix_cl,1,scale))
    rownames(all_Matrix_cl_s) = rownames(all_Matrix_cl)
    colnames(all_Matrix_cl_s) = colnames(all_Matrix_cl)
    #########
    return(all_Matrix_cl_s)
}

#######
col_index = as.character(c(1,3,6,12,18,22,24,30,36,48))
Zebrafish_DEGs_Plot = Merge_DEGs_genes_and_Avg_expression(Zebrafish_DEGs,Avgs_List,col_index)

####### we will ignore 28mo #########


#######
#######
####### Plot = Zebrafish_DEGs_Plot
#######

Zebrafish_DEGs_Plot_Kmeans = Kmeans_function(Zebrafish_DEGs_Plot,k=10)

Kmeans_function <- function(Plot,k=8){
    ############
    #####
    k_up = which(Plot > 4)
    k_down = which(Plot < -4)
    Plot[k_up] = 4
    Plot[k_down] = -4
    ############
    set.seed(123)
    kc = kmeans(Plot, k, iter.max = 100,nstart=20)
    ####
    kc_dat = data.frame(genes = names(kc$cluster),cluster=as.numeric(kc$cluster))
    #### 
    head(kc_dat)
    ####
    return(kc_dat)
}

######
###### 先看一下 cluster ######
######

library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Zebrafish_DEGs_Plot),split="__"),function(x) x[[1]])
names(celltype) = rownames(Zebrafish_DEGs_Plot)
celltype_colors <- c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")
ra <- rowAnnotation(CellType = celltype,
                    col = list(CellType = celltype_colors),
                    annotation_name_side = "top")

row_sp = as.factor(Zebrafish_DEGs_Plot_Kmeans$cluster)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


labels = c("MG__b2m","MG__apoeb","MG__sparc")
at = match(labels,rownames(Zebrafish_DEGs_Plot))

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
png('Zebrafish_DEGs.png',height=7000,width=4000,res=72*12)
Heatmap(Zebrafish_DEGs_Plot, name = "XX",right_annotation = ra,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

########
######## reorder !!!! ####
########

new_order = c(4,3,8,1,2,6,7,10,5,9)
Zebrafish_DEGs_Plot_Kmeans_order = Zebrafish_DEGs_Plot_Kmeans
Zebrafish_DEGs_Plot_Kmeans_order$cluster = match(Zebrafish_DEGs_Plot_Kmeans_order$cluster,new_order)


row_sp = as.factor(Zebrafish_DEGs_Plot_Kmeans_order$cluster)

#################
#################
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
load(file="Zebrafish_DEGs_Plot_Kmeans_order")
load(file="Zebrafish_DEGs_Plot")

####### for Heatmap, order the clusters, young-middle-old #######
new_order2 = c(1,2,3,4,9,10,5,6,7,8)
Zebrafish_DEGs_Plot_Kmeans_order2 = Zebrafish_DEGs_Plot_Kmeans_order
Zebrafish_DEGs_Plot_Kmeans_order2$cluster = match(Zebrafish_DEGs_Plot_Kmeans_order$cluster,new_order2)

########
library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Zebrafish_DEGs_Plot),split="__"),function(x) x[[1]])
names(celltype) = rownames(Zebrafish_DEGs_Plot)
celltype_colors <- c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")
ra <- rowAnnotation(CellType = celltype,
                    col = list(CellType = celltype_colors),
                    annotation_name_side = "top")

row_sp = as.factor(Zebrafish_DEGs_Plot_Kmeans_order2$cluster)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


labels = c("MG__b2m","MG__apoeb","MG__sparc")
at = match(labels,rownames(Zebrafish_DEGs_Plot))


setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
png('Zebrafish_DEGs3.png',height=7000,width=4000,res=72*12)
Heatmap(Zebrafish_DEGs_Plot, name = "XX",right_annotation = ra,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-4,-2,0,2,4)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

#########
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
save(Zebrafish_DEGs_Plot_Kmeans_order,file="Zebrafish_DEGs_Plot_Kmeans_order")
save(Zebrafish_DEGs_Plot,file="Zebrafish_DEGs_Plot")

######### Plot = Zebrafish_DEGs_Plot
######### Kmeans_order = Zebrafish_DEGs_Plot_Kmeans_order

Plot_dot <- function(Plot,Kmeans_order){
    ######
    library(reshape2)
    all_Matrix_cl_s_Plot = melt(Plot)
    all_Matrix_cl_s_Plot$class = Kmeans_order$cluster[match(all_Matrix_cl_s_Plot$Var1,Kmeans_order$genes)]
    #####
    levels_tag = levels(as.factor(all_Matrix_cl_s_Plot$Var2))
    all_Matrix_cl_s_Plot$Var2_1 = as.numeric(match(all_Matrix_cl_s_Plot$Var2,levels_tag))
    ######
    all_Matrix_cl_s_Plot$class = factor(all_Matrix_cl_s_Plot$class,levels=c(1:10))
    ######
    return(all_Matrix_cl_s_Plot)
}

Zebrafish_DEGs_Plot2 = Plot_dot(Zebrafish_DEGs_Plot,Zebrafish_DEGs_Plot_Kmeans_order2)


library(ggplot2)

############----------------------------------------------------------------------------------------------------------------------------------------
ggplot(Zebrafish_DEGs_Plot2, aes(x = Var2_1, y = value)) + 
  geom_point(
    size = 0, 
    position = position_jitter(width = 0.2, height = 0.2),
    alpha = 0.01
  ) + 
  stat_summary(
    fun = mean, 
    geom = "line", 
    color = "red", 
    size = 1
  ) +
  facet_wrap(
    ~ class, 
    ncol = 1, 
    strip.position = "right"
  ) +
  scale_x_discrete(breaks = NULL) +
  theme_classic() + 
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    strip.text       = element_blank(),      # remove facet (strip) labels
    legend.position  = "none",               # remove any legends
    axis.title.x     = element_blank()       # ensure x-axis title is gone
  ) +
  ylab("")

ggsave("Zebrafish_trend2.png", height = 8, width = 1.2)
############----------------------------------------------------------------------------------------------------------------------------------------


########
########
########





########
########--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######## 
######## Next for Mouse ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_pseudo_count <- readRDS("Mouse_pseudo_count_2025")

names(Mouse_pseudo_count)

########
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

######
######

Mouse_DEGs <- list()

for(i in 1:length(Mouse_pseudo_count_input)){
    print(i)
    #####
    res = Identify_DEGs(Mouse_pseudo_count_input[[i]])
    #####
    Mouse_DEGs <- c(Mouse_DEGs,list(res))
}

######
names(Mouse_DEGs) <- names(Mouse_pseudo_count_input)
######

col_index = as.character(c(5,12,17,32,49,68,91,108,120))
Mouse_DEGs_Plot = Merge_DEGs_genes_and_Avg_expression(Mouse_DEGs,Mouse_pseudo_count_input,col_index)

#######
Mouse_DEGs_Plot_Kmeans = Kmeans_function(Mouse_DEGs_Plot,k=10)


library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Mouse_DEGs_Plot),split="__"),function(x) x[[1]])
names(celltype) = rownames(Mouse_DEGs_Plot)
celltype_colors <- c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")
ra <- rowAnnotation(CellType = celltype,
                    col = list(CellType = celltype_colors),
                    annotation_name_side = "top")

row_sp = as.factor(Mouse_DEGs_Plot_Kmeans$cluster)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


labels = c("MG__Apoe","MG__Clu","MG__Aqp4","MG__Notch1","MG__B2m","MG__Stat1")
at = match(labels,rownames(Mouse_DEGs_Plot))

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
png('Mouse_DEGs.png',height=7000,width=4000,res=72*12)
Heatmap(Mouse_DEGs_Plot, name = "XX",right_annotation = ra,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()


######
######
######


new_order = c(8,6,9,3,5,10,7,1,4,2)
Mouse_DEGs_Plot_Kmeans_order = Mouse_DEGs_Plot_Kmeans
Mouse_DEGs_Plot_Kmeans_order$cluster = match(Mouse_DEGs_Plot_Kmeans_order$cluster,new_order)


row_sp = as.factor(Mouse_DEGs_Plot_Kmeans_order$cluster)


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
png('Mouse_DEGs2.png',height=7000,width=4000,res=72*12)
Heatmap(Mouse_DEGs_Plot, name = "XX",right_annotation = ra,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-4,-2,0,2,4)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
save(Mouse_DEGs_Plot_Kmeans_order,file="Mouse_DEGs_Plot_Kmeans_order")
save(Mouse_DEGs_Plot,file="Mouse_DEGs_Plot")

#########
#########
#########


Mouse_DEGs_Plot2 = Plot_dot(Mouse_DEGs_Plot,Mouse_DEGs_Plot_Kmeans_order)


library(ggplot2)
ggplot(Mouse_DEGs_Plot2, aes(x = Var2_1, y = value)) + 
  geom_point(size=0, position = position_jitter(width = 0.2, height = 0.2),alpha=0.01) + facet_wrap(~ class, ncol = 1, strip.position = "right" ) +
  stat_summary(fun = mean, geom = "line", color = "red", size = 1) +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")

ggsave("Mouse_trend2.png",height=12,width=2)


######
########### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
load(file="Mouse_DEGs_Plot_Kmeans_order")
load(file="Mouse_DEGs_Plot")

####### for Heatmap, order the clusters, young-middle-old #######
new_order2 = c(1,2,3,4,9,10,5,6,7,8)
Mouse_DEGs_Plot_Kmeans_order2 = Mouse_DEGs_Plot_Kmeans_order
Mouse_DEGs_Plot_Kmeans_order2$cluster = match(Mouse_DEGs_Plot_Kmeans_order$cluster,new_order2)

#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Mouse_DEGs_Plot),split="__"),function(x) x[[1]])
names(celltype) = rownames(Mouse_DEGs_Plot)
celltype_colors <- c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")
ra <- rowAnnotation(CellType = celltype,
                    col = list(CellType = celltype_colors),
                    annotation_name_side = "top")

row_sp = as.factor(Mouse_DEGs_Plot_Kmeans_order2$cluster)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


labels = c("MG__Apoe","MG__Clu","MG__Aqp4","MG__Notch1","MG__B2m","MG__Stat1")
at = match(labels,rownames(Mouse_DEGs_Plot))

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
png('Mouse_DEGs3.png',height=7000,width=4000,res=72*12)
Heatmap(Mouse_DEGs_Plot, name = "XX",right_annotation = ra,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-4,-2,0,2,4)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Mouse_DEGs_Plot2 = Plot_dot(Mouse_DEGs_Plot,Mouse_DEGs_Plot_Kmeans_order2)
                  

############----------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)

ggplot(Mouse_DEGs_Plot2, aes(x = Var2_1, y = value)) + 
  geom_point(
    size = 0, 
    position = position_jitter(width = 0.2, height = 0.2),
    alpha = 0.01
  ) + 
  stat_summary(
    fun = mean, 
    geom = "line", 
    color = "red", 
    size = 1
  ) +
  facet_wrap(
    ~ class, 
    ncol = 1, 
    strip.position = "right"
  ) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = c(-2,0,2)) +
  theme_classic() + 
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    strip.text       = element_blank(),      # remove facet (strip) labels
    legend.position  = "none",               # remove any legends
    axis.title.x     = element_blank()       # ensure x-axis title is gone
  ) +
  ylab("")

ggsave("Mouse_trend2.png", height = 8, width = 1.2)













                  
###########
########### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###########


########### for Human #####
###########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

########## 
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

#####
#####

#####
Human_DEGs <- list()

for(i in 1:length(Human_pseudo_count)){
    print(i)
    #####
    res = Identify_DEGs(Human_pseudo_count[[i]])
    #####
    Human_DEGs <- c(Human_DEGs,list(res))
}

######
names(Human_DEGs) <- names(Human_pseudo_count)


#######
col_index = as.character(c(10,25,37.5,50,60,70,80,90))
Human_DEGs_Plot = Merge_DEGs_genes_and_Avg_expression(Human_DEGs,Human_pseudo_count,col_index)



#######

#######
Human_DEGs_Plot_Kmeans = Kmeans_function(Human_DEGs_Plot,k=18)


library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Human_DEGs_Plot),split="__"),function(x) x[[1]])
names(celltype) = rownames(Human_DEGs_Plot)
celltype_colors <- c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")


####
celltype = sapply(strsplit(rownames(Human_DEGs_Plot),split="_"),function(x) x[[1]])
sex = sapply(strsplit(rownames(Human_DEGs_Plot),split="_"),function(x) x[[2]])

anno_df = data.frame(
  celltype = celltype,
  sex = sex
)

col_list <- list(
  celltype = c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64"),
  sex = c("M"="blue","F"="red")
)

row_anno <- rowAnnotation(df = anno_df, col = col_list)

#######
#######

row_sp = as.factor(Human_DEGs_Plot_Kmeans$cluster)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


labels = c("MG_M__APOE")
at = match(labels,rownames(Human_DEGs_Plot))

setwd("/zp1/data/share/Human_aging_new")

png('Human_DEGs.png',height=7000,width=4500,res=72*12)
Heatmap(Human_DEGs_Plot, name = "XX",right_annotation = row_anno,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()


###### rm clusters: ##########
######

rm_clusters = c(9,8,13,16,18,2)

#######
Human_DEGs_Plot_Kmeans_cl = Human_DEGs_Plot_Kmeans[which(Human_DEGs_Plot_Kmeans$cluster %in% rm_clusters == F),]
Human_DEGs_Plot_cl = Human_DEGs_Plot[which(rownames(Human_DEGs_Plot) %in% Human_DEGs_Plot_Kmeans_cl$genes == T),]

#######


####
celltype = sapply(strsplit(rownames(Human_DEGs_Plot_cl),split="_"),function(x) x[[1]])
sex = sapply(strsplit(rownames(Human_DEGs_Plot_cl),split="_"),function(x) x[[2]])

anno_df = data.frame(
  celltype = celltype,
  sex = sex
)

col_list <- list(
  celltype = c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64"),
  sex = c("M"="blue","F"="red")
)

row_anno <- rowAnnotation(df = anno_df, col = col_list)

#######
#######

row_sp = as.factor(Human_DEGs_Plot_Kmeans_cl$cluster)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


labels = c("MG_M__APOE")
at = match(labels,rownames(Human_DEGs_Plot_cl))

setwd("/zp1/data/share/Human_aging_new")

png('Human_DEGs2.png',height=7000,width=4500,res=72*12)
Heatmap(Human_DEGs_Plot_cl, name = "XX",right_annotation = row_anno,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()


######
###### Next we will order ###
######

new_order = c(3,7,6,12,14,10,11,1,5,15,4,17)
Human_DEGs_Plot_Kmeans_order = Human_DEGs_Plot_Kmeans_cl
Human_DEGs_Plot_Kmeans_order$cluster = match(Human_DEGs_Plot_Kmeans_order$cluster,new_order)


row_sp = as.factor(Human_DEGs_Plot_Kmeans_order$cluster)

####
labels = c("RPE_M__APOE","Cone_F__STAT1","RPE_M__BTN3A3","MG_M__BTN3A3","AC_F__BTN3A3","RPE_M__DKK3")
at = match(labels,rownames(Human_DEGs_Plot_cl))

rownames(Human_DEGs_Plot_cl)[grep("__DKK3",rownames(Human_DEGs_Plot_cl))]

png('Human_DEGs4.png',height=7000,width=4500,res=72*12)
Heatmap(Human_DEGs_Plot_cl, name = "XX",right_annotation = row_anno,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-4,-2,0,2,4)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

#####
#####
#####


setwd("/zp1/data/share/Human_aging_new")
save(Human_DEGs_Plot_Kmeans_order,file="Human_DEGs_Plot_Kmeans_order")
save(Human_DEGs_Plot_cl,file="Human_DEGs_Plot_cl")

######
######
table(Human_DEGs_Plot_Kmeans_order$cluster)
######

Plot_dot <- function(Plot,Kmeans_order){
    ######
    library(reshape2)
    all_Matrix_cl_s_Plot = melt(Plot)
    all_Matrix_cl_s_Plot$class = Kmeans_order$cluster[match(all_Matrix_cl_s_Plot$Var1,Kmeans_order$genes)]
    #####
    levels_tag = levels(as.factor(all_Matrix_cl_s_Plot$Var2))
    all_Matrix_cl_s_Plot$Var2_1 = as.numeric(match(all_Matrix_cl_s_Plot$Var2,levels_tag))
    ######
    all_Matrix_cl_s_Plot$class = factor(all_Matrix_cl_s_Plot$class,levels=c(1:12))
    ######
    return(all_Matrix_cl_s_Plot)
}

Human_DEGs_Plot2 = Plot_dot(Human_DEGs_Plot_cl,Human_DEGs_Plot_Kmeans_order)


library(ggplot2)
ggplot(Human_DEGs_Plot2, aes(x = Var2_1, y = value)) + 
  geom_point(size=0, position = position_jitter(width = 0.2, height = 0.2),alpha=0.01) + facet_wrap(~ class, ncol = 1, strip.position = "right" ) +
  stat_summary(fun = mean, geom = "line", color = "red", size = 1) +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + xlab("") + ylab("")

ggsave("Human_trend2.png",height=12,width=2)



#########
######### Done!!!
#########



