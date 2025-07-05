#############
#######
#############
####### Female and Male #########
#############

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

#############
############# load the DEGs for Human #####
#############


##
######### for Human: ########
##

setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot_Kmeans_order")

#######-----for M ######

kc = Human_DEGs_Plot_Kmeans_order
Upclusters = c(5,6,7,8,9)
Downclusters = c(1,2,3,4)


kc$CT = sapply(strsplit(kc$genes,split="_"),function(x) x[[1]])
kc$Gene = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])

kc$Class = "Unknown"
kc$Class[which(kc$cluster %in% Upclusters)] = "UP"
kc$Class[which(kc$cluster %in% Downclusters)] = "DOWN"

#######
extract_gender <- function(x) {
  # x: 字符串向量，比如 c("sample_F__01", "ID_M__02", ...)
  # 1. 提取 "_" 和 "__" 之间的内容
  codes <- sub(".*_([^_]+)__.*", "\\1", x)
  # 2. 映射 F->Female, M->Male, 其他映射为 NA
  gender <- ifelse(codes == "F", "Female",
                   ifelse(codes == "M", "Male", NA))
  return(gender)
}
########
Human_DEGs = kc
Human_DEGs$gender = extract_gender(kc$genes)

########
######## find the overlap genes ######
########

table(Human_DEGs$Class)
Human_DEGs_cl = Human_DEGs[which(Human_DEGs$Class!="Unknown"),]

Human_DEGs_cl$Class2 = paste(Human_DEGs_cl$Gene,Human_DEGs_cl$CT,Human_DEGs_cl$Class,sep="_")

Human_DEGs_cl_sp = split(Human_DEGs_cl,Human_DEGs_cl$gender)

########
Human_DEGs_Common_overlap_index = which(Human_DEGs_cl_sp$Female$Class2 %in% Human_DEGs_cl_sp$Male$Class2 == T)

########
Human_DEGs_Common_overlap =  Human_DEGs_cl_sp$Female$Class2[Human_DEGs_Common_overlap_index]

length(Human_DEGs_Common_overlap)

########
########
########

Index = data.frame(Female=Human_DEGs_cl_sp$Female$Class2,Male="NO")

# 拷贝一份原始向量
f <- Index$Female

# 1. 把所有 "_UP" 先替换成 "__TMP__"
f <- gsub("_UP$", "__TMP__", f)

# 2. 把所有 "_DOWN" 替换成 "_UP"
f <- gsub("_DOWN$", "_UP", f)

# 3. 再把 "__TMP__" 替换成 "_DOWN"
f <- gsub("__TMP__$", "_DOWN", f)

# 赋回去
Index$Male <- f

####
Index$TAG = "NotFind"
k = which(Index$Male %in% Human_DEGs_cl_sp$Male$Class2 == T)
Index$TAG[k] = "YES"

####
Human_DEGs_Divergent_overlap = Index[which(Index$TAG=="YES"),]


#####
##### Consistent #####
#####


#####
##### load the avg expression matrix #####
#####

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

######
######
######
#### Let us see the Human_DEGs_Common_overlap ####
######
######
Human_DEGs_Common_overlap_tab = data.frame(Index=Human_DEGs_Common_overlap,CT="NO",Gene="NO",Change="NO")
Human_DEGs_Common_overlap_tab$CT = sapply(strsplit(Human_DEGs_Common_overlap_tab$Index,split="_"),function(x) x[[2]])
Human_DEGs_Common_overlap_tab$Gene = sapply(strsplit(Human_DEGs_Common_overlap_tab$Index,split="_"),function(x) x[[1]])
Human_DEGs_Common_overlap_tab$Change = sapply(strsplit(Human_DEGs_Common_overlap_tab$Index,split="_"),function(x) x[[3]])

Human_DEGs_Common_overlap_tab$Male = paste0(Human_DEGs_Common_overlap_tab$CT,"_M__",Human_DEGs_Common_overlap_tab$Gene)
Human_DEGs_Common_overlap_tab$Female = paste0(Human_DEGs_Common_overlap_tab$CT,"_F__",Human_DEGs_Common_overlap_tab$Gene)

Avgs_List = Human_pseudo_count


#######

Human_DEGs_Common_overlap_tab_out = Human_DEGs_Common_overlap_tab[,c("CT","Gene","Change","Male","Female")]
Human_DEGs_Common_overlap_tab_out$Male = Human_DEGs_Common_overlap_tab$Change
Human_DEGs_Common_overlap_tab_out$Female = Human_DEGs_Common_overlap_tab$Change

k1 = which(Human_DEGs_Common_overlap_tab_out$Change == "UP")
k2 = which(Human_DEGs_Common_overlap_tab_out$Change == "DOWN")

Human_DEGs_Common_overlap_tab_out$Female[k1] = "Old"
Human_DEGs_Common_overlap_tab_out$Male[k1] = "Old"
Human_DEGs_Common_overlap_tab_out$Female[k2] = "Young"
Human_DEGs_Common_overlap_tab_out$Male[k2] = "Young"

Human_DEGs_Common_overlap_tab_out$Change = NULL
write_xlsx(Human_DEGs_Common_overlap_tab_out, path = "Human_SEX_common_by_celltype.xlsx")



Human_DEGs_Divergent_overlap

table(Human_DEGs_Divergent_overlap$TAG)

Human_DEGs_Divergent_overlap$CT = sapply(strsplit(Human_DEGs_Divergent_overlap$Female,split="_"),function(x) x[[2]]) 
Human_DEGs_Divergent_overlap$Gene = sapply(strsplit(Human_DEGs_Divergent_overlap$Female,split="_"),function(x) x[[1]]) 

Human_DEGs_Divergent_overlap$Index = sapply(strsplit(Human_DEGs_Divergent_overlap$Female,split="_"),function(x) x[[3]]) 

k1 = which(Human_DEGs_Divergent_overlap$Index == "DOWN")
k2 = which(Human_DEGs_Divergent_overlap$Index == "UP")

Human_DEGs_Divergent_overlap$Female[k1] = "Young"
Human_DEGs_Divergent_overlap$Female[k2] = "Old"
Human_DEGs_Divergent_overlap$Male[k1] = "Old"
Human_DEGs_Divergent_overlap$Male[k2] = "Young"

Human_DEGs_Divergent_overlap = Human_DEGs_Divergent_overlap[,c("CT","Gene","Male","Female")]

write_xlsx(Human_DEGs_Divergent_overlap, path = "Human_SEX_divergent_by_celltype.xlsx")

#######
#######
#######

DEGs_tab = Human_DEGs_Common_overlap_tab
col_index = as.character(c(10,25,37.5,50,60,70,80,90))

#######
Merge_DEGs_genes_and_Avg_expression_Common <- function(DEGs_tab,Avgs_List,col_index){
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
    #########
    all_Matrix = do.call("rbind",Avgs_List2)
    ######### all_Avg[grep("Clu",all_Avg$Gene),] #########
    ######### all_Avg[grep("Xkr4",all_Avg$Gene),] ########
    all_Matrix_cl_F = all_Matrix[which(rownames(all_Matrix) %in% DEGs_tab$Female == T),]
    all_Matrix_cl_M = all_Matrix[which(rownames(all_Matrix) %in% DEGs_tab$Male == T),]
    #########
    ######### merge and order ########
    ######### to Z-score #############
    #########
    all_Matrix_cl_F = all_Matrix_cl_F[,col_index]
    all_Matrix_cl_M = all_Matrix_cl_M[,col_index]
    #########
    #########
    all_Matrix_cl_s_F = t(apply(all_Matrix_cl_F,1,scale))
    rownames(all_Matrix_cl_s_F) = rownames(all_Matrix_cl_F)
    colnames(all_Matrix_cl_s_F) = colnames(all_Matrix_cl_F)
    #####
    all_Matrix_cl_s_M = t(apply(all_Matrix_cl_M,1,scale))
    rownames(all_Matrix_cl_s_M) = rownames(all_Matrix_cl_M)
    colnames(all_Matrix_cl_s_M) = colnames(all_Matrix_cl_M)
    #########
    rownames(all_Matrix_cl_s_F) = gsub("_F__","_",rownames(all_Matrix_cl_s_F))
    rownames(all_Matrix_cl_s_M) = gsub("_M__","_",rownames(all_Matrix_cl_s_M))
    ######### merge the 2 matrix ######
    #########
    all.equal(rownames(all_Matrix_cl_s_F),rownames(all_Matrix_cl_s_M))
    all_Matrix_cl_s_M <- all_Matrix_cl_s_M[rownames(all_Matrix_cl_s_F), , drop = FALSE]
    all_Matrix_cl_s_FM <- cbind(all_Matrix_cl_s_F, all_Matrix_cl_s_M)
    #########
    #### add colsp and add rowsp #######
    #########
    colsp = c(rep('Female',dim(all_Matrix_cl_s_F)[2]),rep('Male',dim(all_Matrix_cl_s_F)[2]))
    colsp = factor(colsp)
    #########
    DEGs_tab$Index2 = paste0(DEGs_tab$CT,"_",DEGs_tab$Gene)
    rowsp = DEGs_tab$Change[match(DEGs_tab$Index2,rownames(all_Matrix_cl_s_FM))]
    rowsp = factor(rowsp)
    #########
    return(list(Mat=all_Matrix_cl_s_FM,rowsp=rowsp,colsp=colsp))
}


rownames(Human_pseudo_count[[1]]$bulk_mat)


#####
Human_DEGs_Plot_common = Merge_DEGs_genes_and_Avg_expression_Common(Human_DEGs_Common_overlap_tab,Human_pseudo_count,col_index)

#####
#####



setwd("/zp1/data/share/Human_aging_new")


library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Human_DEGs_Plot_common$Mat),split="_"),function(x) x[[1]])
names(celltype) = rownames(Human_DEGs_Plot_common$Mat)
celltype_colors <- c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")


####
anno_df = data.frame(
  celltype = celltype
)

col_list <- list(
  celltype = c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")
)

row_anno <- rowAnnotation(df = anno_df, col = col_list)

#######
#######

row_sp = Human_DEGs_Plot_common$rowsp
col_sp = Human_DEGs_Plot_common$colsp

table(Human_DEGs_Plot_common$rowsp)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


png('Human_DEGs_common.png',height=4000,width=3500,res=72*12)
Heatmap(Human_DEGs_Plot_common$Mat, name = "XX",right_annotation = row_anno,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp,column_split=col_sp)
dev.off()

##### 左边 Female 右边 Male ######
##### 然后加上 label celltypes ###
#####






##### 左边 Female 右边 Male #######
#####
#####






##### 接下来 divergent ######
#####

DEGs_tab = Human_DEGs_Divergent_overlap
col_index = as.character(c(10,25,37.5,50,60,70,80,90))

Merge_DEGs_genes_and_Avg_expression_Divergent <- function(DEGs_tab,Avgs_List,col_index){
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
    #########
    all_Matrix = do.call("rbind",Avgs_List2)
    ######### all_Avg[grep("Clu",all_Avg$Gene),] #########
    ######### all_Avg[grep("Xkr4",all_Avg$Gene),] ########
    DEGs_tab$CT = sapply(strsplit(DEGs_tab$Female,split="_"),function(x) x[[2]])
    DEGs_tab$Gene = sapply(strsplit(DEGs_tab$Female,split="_"),function(x) x[[1]])
    DEGs_tab$Female_Index = paste0(DEGs_tab$CT,"_F__",DEGs_tab$Gene)
    DEGs_tab$Male_Index = paste0(DEGs_tab$CT,"_M__",DEGs_tab$Gene)
    #########
    all_Matrix_cl_F = all_Matrix[which(rownames(all_Matrix) %in% DEGs_tab$Female_Index == T),]
    all_Matrix_cl_M = all_Matrix[which(rownames(all_Matrix) %in% DEGs_tab$Male_Index == T),]
    #########
    ######### merge and order ########
    ######### to Z-score #############
    #########
    all_Matrix_cl_F = all_Matrix_cl_F[,col_index]
    all_Matrix_cl_M = all_Matrix_cl_M[,col_index]
    #########
    #########
    all_Matrix_cl_s_F = t(apply(all_Matrix_cl_F,1,scale))
    rownames(all_Matrix_cl_s_F) = rownames(all_Matrix_cl_F)
    colnames(all_Matrix_cl_s_F) = colnames(all_Matrix_cl_F)
    #####
    all_Matrix_cl_s_M = t(apply(all_Matrix_cl_M,1,scale))
    rownames(all_Matrix_cl_s_M) = rownames(all_Matrix_cl_M)
    colnames(all_Matrix_cl_s_M) = colnames(all_Matrix_cl_M)
    #########
    rownames(all_Matrix_cl_s_F) = gsub("_F__","_",rownames(all_Matrix_cl_s_F))
    rownames(all_Matrix_cl_s_M) = gsub("_M__","_",rownames(all_Matrix_cl_s_M))
    ######### merge the 2 matrix ######
    #########
    all.equal(rownames(all_Matrix_cl_s_F),rownames(all_Matrix_cl_s_M))
    all_Matrix_cl_s_M <- all_Matrix_cl_s_M[rownames(all_Matrix_cl_s_F), , drop = FALSE]
    all_Matrix_cl_s_FM <- cbind(all_Matrix_cl_s_F, all_Matrix_cl_s_M)
    #########
    #### add colsp and add rowsp #######
    #########
    colsp = c(rep('Female',dim(all_Matrix_cl_s_F)[2]),rep('Male',dim(all_Matrix_cl_s_F)[2]))
    colsp = factor(colsp)
    #########
    DEGs_tab$Index2 = paste0(DEGs_tab$CT,"_",DEGs_tab$Gene)
    DEGs_tab$Change = sapply(strsplit(DEGs_tab$Female,split="_"),function(x) x[[3]])
    rowsp = DEGs_tab$Change[match(DEGs_tab$Index2,rownames(all_Matrix_cl_s_FM))]
    rowsp = factor(rowsp)
    #########
    return(list(Mat=all_Matrix_cl_s_FM,rowsp=rowsp,colsp=colsp))
}

#####
#####
Human_DEGs_Plot_div = Merge_DEGs_genes_and_Avg_expression_Divergent(Human_DEGs_Divergent_overlap,Human_pseudo_count,col_index)


#####
#####

setwd("/zp1/data/share/Human_aging_new")

library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Human_DEGs_Plot_div$Mat),split="_"),function(x) x[[1]])
names(celltype) = rownames(Human_DEGs_Plot_div$Mat)
celltype_colors <- c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")


####
anno_df = data.frame(
  celltype = celltype
)

col_list <- list(
  celltype = c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64")
)

row_anno <- rowAnnotation(df = anno_df, col = col_list)

#######
#######

row_sp = Human_DEGs_Plot_div$rowsp
col_sp = Human_DEGs_Plot_div$colsp

table(Human_DEGs_Plot_div$rowsp)

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))


png('Human_DEGs_div.png',height=4000,width=3500,res=72*12)
Heatmap(Human_DEGs_Plot_div$Mat, name = "XX",right_annotation = row_anno,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp,column_split=col_sp)
dev.off()


######
###### Output the gene list ########
######
Human_DEGs_common_List = split(rownames(Human_DEGs_Plot_common$Mat),Human_DEGs_Plot_common$rowsp)

for(i in 1:length(Human_DEGs_common_List)){
    tmp = sapply(strsplit(Human_DEGs_common_List[[i]],split="_"),function(x) x[[2]])
    Human_DEGs_common_List[[i]] = tmp
}

Human_DEGs_div_List = split(rownames(Human_DEGs_Plot_div$Mat),Human_DEGs_Plot_div$rowsp)

for(i in 1:length(Human_DEGs_div_List)){
    tmp = sapply(strsplit(Human_DEGs_div_List[[i]],split="_"),function(x) x[[2]])
    Human_DEGs_div_List[[i]] = tmp
}

########
######## Next we will perform GO and KEGG analysis #######
########

setwd("/zp1/data/share/Human_aging_new")

save(Human_DEGs_common_List,file="Human_DEGs_common_List")
save(Human_DEGs_div_List,file="Human_DEGs_div_List")


########

Enrich_GO_function_Human <- function(tmp_genes){
		library(clusterProfiler)
		library(org.Hs.eg.db)
		#####
		Set1 <- bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		result <- enrichGO(gene = Set1$ENTREZID,
                   #universe = background_genes,
                   OrgDb = org.Hs.eg.db,
                   ont  = "BP",
                   keyType  = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
		)
		#####
		GOtable = result@result
		#####
        for(j in 1:dim(GOtable)[1]){
			########
			tmp_G_sp = unlist(strsplit(GOtable$geneID[j],split='/',fixed=T))
			tmp_G_sp_covert = Set1$SYMBOL[match(tmp_G_sp,Set1$ENTREZID)]
			########
			tmp_G_sp_covert_merge = paste(tmp_G_sp_covert,collapse=";")
			GOtable$geneID[j] = tmp_G_sp_covert_merge

		}
		#####
		#####
        return(GOtable)
        #####
}	

Enrich_KEGG_function_Human <- function(tmp_genes){
		library(clusterProfiler)
		library(org.Hs.eg.db)
        #####
		Set1 <- bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		result <- enrichKEGG(gene  = Set1$ENTREZID,
                 organism  = 'hsa',
                 pvalueCutoff = 0.05
	    )
		#####
		GOtable = result@result
		#####
        for(j in 1:dim(GOtable)[1]){
			########
			tmp_G_sp = unlist(strsplit(GOtable$geneID[j],split='/',fixed=T))
			tmp_G_sp_covert = Set1$SYMBOL[match(tmp_G_sp,Set1$ENTREZID)]
			########
			tmp_G_sp_covert_merge = paste(tmp_G_sp_covert,collapse=";")
			GOtable$geneID[j] = tmp_G_sp_covert_merge

		}
		#####
        return(GOtable)
        #####
}	

#######
#######
#######

setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_common_List")
load("Human_DEGs_div_List")


#####
Common_Down <- data.frame(Genes = Human_DEGs_common_List$DOWN,Female="Young",MALE="Young")
Common_Up <- data.frame(Genes = Human_DEGs_common_List$UP,Female="Old",MALE="Old")

Div_Down = data.frame(Genes = Human_DEGs_common_List$DOWN,Female="Young",MALE="Old")
Div_Up = data.frame(Genes = Human_DEGs_div_List$UP,Female="Old",MALE="Young")


#####
Common = rbind(Common_Down,Common_Up)
Div = rbind(Div_Down,Div_Up)

#####
List = list(Common=Common,divergent=Div)

library(writexl)
write_xlsx(List, path = "Human_SEX_DEGs_2025.xlsx")


#####
library(writexl)
# Write to "output.xlsx" — each list element becomes a sheet named after its list name

Human_common_GOKEGG <- list()

####
for(i in 1:length(Human_DEGs_common_List)){
    print(i)
    ######
    tmp = Human_DEGs_common_List[[i]]
    GOres = Enrich_GO_function_Human(tmp)
    KEGGres = Enrich_KEGG_function_Human(tmp)
    ######
    GOres = cbind(data.frame(Class = "GOterms"),GOres)
    KEGGres = cbind(data.frame(Class = "KEGGterms"),KEGGres)
    ######
    GOres = GOres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    KEGGres = KEGGres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    ######
    Mergeres = rbind(GOres,KEGGres)
    ######
    Human_common_GOKEGG <- c(Human_common_GOKEGG,list(Mergeres))
}

######
names(Human_common_GOKEGG) <- names(Human_DEGs_common_List)
######
######

for(i in 1:length(Human_common_GOKEGG)){
    tmp = Human_common_GOKEGG[[i]]
    tmp = tmp[tmp$pvalue < 0.05,]
    tmp = tmp[order(tmp$pvalue),]
    Human_common_GOKEGG[[i]] = tmp
}
######
######

Human_SEX_GO <- c(Human_common_GOKEGG,Human_div_GOKEGG)
names(Human_SEX_GO) <- c("Common_Young","Common_Old","Divergent_F_Young_M_Old","Divergent_F_Old_M_Young")
write_xlsx(Human_SEX_GO, path = "Human_SEX_diff_GOterms_2025.xlsx")


######
######
######

Human_div_GOKEGG <- list()

####
for(i in 1:length(Human_DEGs_div_List)){
    print(i)
    ######
    tmp = Human_DEGs_div_List[[i]]
    GOres = Enrich_GO_function_Human(tmp)
    KEGGres = Enrich_KEGG_function_Human(tmp)
    ######
    GOres = cbind(data.frame(Class = "GOterms"),GOres)
    KEGGres = cbind(data.frame(Class = "KEGGterms"),KEGGres)
    ######
    GOres = GOres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    KEGGres = KEGGres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    ######
    Mergeres = rbind(GOres,KEGGres)
    ######
    Human_div_GOKEGG <- c(Human_div_GOKEGG,list(Mergeres))
}

######
names(Human_div_GOKEGG) <- names(Human_DEGs_div_List)
######

library(dplyr)
library(ggplot2)
library(stringr)

######
######
# 1. 准备数据
df <- as.data.frame(Human_common_GOKEGG[[2]])
top10 <- df %>%
  arrange(pvalue) %>%
  slice(1:10) %>%
  mutate(log10p = -log10(pvalue),Description_wrap = str_wrap(Description, width = 30))

# 2. 画图并添加黑框 steelblue
p <- ggplot(top10, aes(x = log10p, y = reorder(Description_wrap, log10p))) +
  geom_col(fill = "red",width=0.5) +
  geom_text(aes(label = Count),
            hjust = -0.1,
            size = 3) +
  labs(x = "-log10(p-value)",
       y = NULL,
       title = "Top 10 GO terms") +
  theme_classic() +
  theme(
    # Panel（绘图区域）边框
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(5, 20, 5, 5)
  )

# 3. 保存成 PNG
ggsave(
  filename = "go_barplot.png",
  plot     = p,
  width    = 5,     # 图的宽度（英寸）
  height   = 4,     # 图的高度（英寸）
  dpi      = 300    # 分辨率
)





######
###### Next we will perform GO terms analysis #####
######


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R












