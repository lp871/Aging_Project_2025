#########-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########
######### Output DEGs tables for Mouse Zebrafish and Human ##########
#########
#########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate clusterProfiler
R


##
######### for Mouse: ########----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##
##
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
Mouse_UP = kc[which(kc$Class=="UP"),]
Mouse_UP_list = split(Mouse_UP,Mouse_UP$CT)
Mouse_DOWN = kc[which(kc$Class=="DOWN"),]
Mouse_DOWN_list = split(Mouse_DOWN,Mouse_DOWN$CT)

###############
############### Output the Excel files ####
###############

head(Mouse_UP_list[[1]])
names(Mouse_UP_list) <- paste0(names(Mouse_UP_list),"_Old")
head(Mouse_DOWN_list[[1]])
names(Mouse_DOWN_list) <- paste0(names(Mouse_DOWN_list),"_Young")

Mouse_combined = c(Mouse_UP_list,Mouse_DOWN_list)
names(Mouse_combined)

######
###### Output to excel files ######
######


######
######
######

Mouse_combined_GOKEGG <- list()
for(i in 1:length(Mouse_combined)){
    print(i)
    ######
    tmp = Mouse_combined[[i]]
    GOres = Enrich_GO_function_Mouse(tmp$Gene)
    KEGGres = Enrich_KEGG_function_Mouse(tmp$Gene)
    ######
    GOres = cbind(data.frame(Class = "GOterms"),GOres)
    KEGGres = cbind(data.frame(Class = "KEGGterms"),KEGGres)
    ######
    GOres = GOres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    KEGGres = KEGGres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    ######
    Mergeres = rbind(GOres,KEGGres)
    ######
    Mouse_combined_GOKEGG <- c(Mouse_combined_GOKEGG,list(Mergeres))
}

names(Mouse_combined_GOKEGG) <- names(Mouse_combined)


for(i in 1:length(Mouse_combined_GOKEGG)){
    tmp = Mouse_combined_GOKEGG[[i]]
    tmp = tmp[tmp$pvalue < 0.05,]
    Mouse_combined_GOKEGG[[i]] = tmp
}



library(openxlsx)
wb <- createWorkbook()
# 遍历 list，把每个元素写到单独的 sheet
for (i in seq_along(Mouse_combined_GOKEGG)) {
  sheet_name <- names(Mouse_combined_GOKEGG)[i]  # 如果 list 有名字就用名字做 sheet 名
  if (is.null(sheet_name) || sheet_name == "") {
    sheet_name <- paste0("Sheet", i)  # 如果没有名字，用 Sheet1, Sheet2...
  }
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, Mouse_combined_GOKEGG[[i]])
}

# 保存 Excel 文件
saveWorkbook(wb, "Mouse_Aging_DEGs_GOKEGG_May5_2025.xlsx", overwrite = TRUE)


#######
#######













##
######### for Zebrafish: ########
##

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


head(Zebrafish_UP_list[[1]])
names(Zebrafish_UP_list) <- paste0(names(Zebrafish_UP_list),"_Old")
head(Mouse_DOWN_list[[1]])
names(Zebrafish_DOWN_list) <- paste0(names(Zebrafish_DOWN_list),"_Young")

Zebrafish_combined = c(Zebrafish_UP_list,Zebrafish_DOWN_list)
names(Zebrafish_combined)

######
###### Output to excel files ######
######


#######
####### 
#######

Zebrafish_combined_GOKEGG <- list()
for(i in 1:length(Zebrafish_combined)){
    print(i)
    ######
    tmp = Zebrafish_combined[[i]]
    GOres = Enrich_GO_function_Fish(tmp$Gene)
    KEGGres = Enrich_KEGG_function_Fish(tmp$Gene)
    ######
    GOres = cbind(data.frame(Class = "GOterms"),GOres)
    KEGGres = cbind(data.frame(Class = "KEGGterms"),KEGGres)
    ######
    GOres = GOres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    KEGGres = KEGGres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    ######
    Mergeres = rbind(GOres,KEGGres)
    ######
    Zebrafish_combined_GOKEGG <- c(Zebrafish_combined_GOKEGG,list(Mergeres))
}

names(Zebrafish_combined_GOKEGG) <- names(Zebrafish_combined)

for(i in 1:length(Zebrafish_combined_GOKEGG)){
    tmp = Zebrafish_combined_GOKEGG[[i]]
    tmp = tmp[tmp$pvalue < 0.05,]
    Zebrafish_combined_GOKEGG[[i]] = tmp
}



library(openxlsx)
wb <- createWorkbook()
# 遍历 list，把每个元素写到单独的 sheet
for (i in seq_along(Zebrafish_combined_GOKEGG)) {
  sheet_name <- names(Zebrafish_combined_GOKEGG)[i]  # 如果 list 有名字就用名字做 sheet 名
  if (is.null(sheet_name) || sheet_name == "") {
    sheet_name <- paste0("Sheet", i)  # 如果没有名字，用 Sheet1, Sheet2...
  }
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, Zebrafish_combined_GOKEGG[[i]])
}

# 保存 Excel 文件
saveWorkbook(wb, "Zebrafish_Aging_DEGs_GOKEGG_May5_2025.xlsx", overwrite = TRUE)





##
######### for Human: ########
##


				 
setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot_Kmeans_order")

#######-----for M ######

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
extract_gender <- function(x) {
  # x: 字符串向量，比如 c("sample_F__01", "ID_M__02", ...)
  # 1. 提取 "_" 和 "__" 之间的内容
  codes <- sub(".*_([^_]+)__.*", "\\1", x)
  # 2. 映射 F->Female, M->Male, 其他映射为 NA
  gender <- ifelse(codes == "F", "Female",
                   ifelse(codes == "M", "Male", NA))
  return(gender)
}

#####
Human_UP = kc[which(kc$Class=="UP"),]
Human_UP$gender = extract_gender(Human_UP$genes)
Human_UP_list = split(Human_UP,Human_UP$CT)

#####
Human_DOWN = kc[which(kc$Class=="DOWN"),]
Human_DOWN$gender = extract_gender(Human_DOWN$genes)
Human_DOWN_list = split(Human_DOWN,Human_DOWN$CT)


#####
names(Human_UP_list) <- paste0(names(Human_UP_list),"_Old")
names(Human_DOWN_list) <- paste0(names(Human_DOWN_list),"_Young")

Human_combined = c(Human_UP_list,Human_DOWN_list)
names(Human_combined)

######
###### Output to excel files ######
######
library(openxlsx)
# Write to "output.xlsx" — each list element becomes a sheet named after its list name

Human_combined_GOKEGG <- list()

####
for(i in 1:length(Human_combined)){
    print(i)
    ######
    tmp = Human_combined[[i]]
    GOres = Enrich_GO_function_Human(tmp$Gene)
    KEGGres = Enrich_KEGG_function_Human(tmp$Gene)
    ######
    GOres = cbind(data.frame(Class = "GOterms"),GOres)
    KEGGres = cbind(data.frame(Class = "KEGGterms"),KEGGres)
    ######
    GOres = GOres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    KEGGres = KEGGres[,c("Class","ID","Description","pvalue","p.adjust","qvalue","geneID","Count")]
    ######
    Mergeres = rbind(GOres,KEGGres)
    ######
    Human_combined_GOKEGG <- c(Human_combined_GOKEGG,list(Mergeres))
}

names(Human_combined_GOKEGG) <- names(Human_combined)

######
######

for(i in 1:length(Human_combined_GOKEGG)){
    tmp = Human_combined_GOKEGG[[i]]
    tmp = tmp[tmp$pvalue < 0.05,]
    Human_combined_GOKEGG[[i]] = tmp
}

library(openxlsx)
wb <- createWorkbook()
# 遍历 list，把每个元素写到单独的 sheet
for (i in seq_along(Human_combined_GOKEGG)) {
  sheet_name <- names(Human_combined_GOKEGG)[i]  # 如果 list 有名字就用名字做 sheet 名
  if (is.null(sheet_name) || sheet_name == "") {
    sheet_name <- paste0("Sheet", i)  # 如果没有名字，用 Sheet1, Sheet2...
  }
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, Human_combined_GOKEGG[[i]])
}

# 保存 Excel 文件
saveWorkbook(wb, "Human_Aging_DEGs_GOKEGG_May5_2025.xlsx", overwrite = TRUE)


######
###### Next for the GO term and KEGG terms ########
######


######
###### Next for GO terms !!!! ######
######



###### Zebrafish ######----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###
###
######
Enrich_GO_function_Fish <- function(tmp_genes){
		library(clusterProfiler)
		library(org.Dr.eg.db)
        #####
		##### for Zebrafish: #######
		#####
		Set1 <- bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dr.eg.db")
		result <- enrichGO(gene = Set1$ENTREZID,
                   OrgDb = org.Dr.eg.db,
                   ont  = "BP",
                   keyType  = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
        )
		#####
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
		##### replace the number to gene names ########
		#####
        return(GOtable)
        #####
}
######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Enrich_KEGG_function_Fish <- function(tmp_genes){
		library(clusterProfiler)
		library(org.Dr.eg.db)
        #####
		##### for Zebrafish: #######
		#####
		Set1 <- bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dr.eg.db")
		result <- enrichKEGG(gene  = Set1$ENTREZID,
                 organism  = 'dre',
                 pvalueCutoff = 0.05
	    )
		#####
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


########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




















#############
############# Mouse GO terms ########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############
###
###

Enrich_GO_function_Mouse <- function(tmp_genes){
		library(clusterProfiler)
		library(org.Mm.eg.db)
        #####
		##### for Zebrafish: #######
		#####
		Set1 <- bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
		result <- enrichGO(gene = Set1$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont  = "BP",
                   keyType  = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
		)
		#####
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



Enrich_KEGG_function_Mouse <- function(tmp_genes){
		library(clusterProfiler)
		library(org.Mm.eg.db)
		#####
		Set1 <- bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
		result <- enrichKEGG(gene  = Set1$ENTREZID,
                 organism  = 'mmu',
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
}	


#############
############# Human GO terms ########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############

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


#### tmp_genes = Human_combined[[1]]$Gene


###
###
	
Enrich_KEGG_function_Human <- function(tmp_genes){
		library(clusterProfiler)
		library(org.Hs.eg.db)
        #####
		Set1 <- bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		result <- enrichKEGG(gene  = Set1$ENTREZID,
                 organism  = 'hsa',
                 pvalueCutoff = 0.05,
				 use_internal_data = FALSE
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




###########
###########
###########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate clusterProfiler
R


names(Zebrafish_combined_GOKEGG)
names(Human_combined_GOKEGG)
names(Mouse_combined_GOKEGG)


#########
#########
#########
kc_index = Zebrafish_combined_GOKEGG[grep("_Young",names(Zebrafish_combined_GOKEGG))]
  
Get_overlap_withinGOtables <- function(kc_index){
    ######
    all_ID = c()
	all_Terms = c()
    for(i in 1:length(kc_index)){
		all_ID = c(all_ID,kc_index[[i]]$ID)
		all_Terms = c(all_Terms,kc_index[[i]]$Description)
	}
	k = !duplicated(all_ID)
	all_ID = all_ID[k]
	all_Terms = all_Terms[k]
    ######
    tmp_mat = matrix(NA,nrow=length(all_ID),ncol=length(kc_index))
    rownames(tmp_mat) <- all_ID
    colnames(tmp_mat) <- names(kc_index)
    ######
    for(i in 1:dim(tmp_mat)[1]){
        for(j in 1:dim(tmp_mat)[2]){
            tmp_ct = colnames(tmp_mat)[j]
            tmp_g = rownames(tmp_mat)[i]
            ####
			tmp_tab = kc_index[[tmp_ct]]
            k = which(tmp_tab$ID == tmp_g)
            if(length(k) > 0){
                tmp_mat[i,j] = tmp_tab$pvalue[k]
            }
        }
    }
    ######
    tmp_mat = data.frame(tmp_mat)
	#######
    tmp_mat$overlap = apply(tmp_mat,1,function(x) length(which(is.na(x) == F)))
	tmp_mat$Description = all_Terms
    tmp_mat = tmp_mat[order(tmp_mat$overlap,decreasing=T),]
    #######
    #######
    #######
    #######
    tmp_mat = tmp_mat[which(tmp_mat$overlap > 1),]
    #######
    return(tmp_mat)
}


Zebrafish_Old_overlap = Get_overlap_withinGOtables(Zebrafish_combined_GOKEGG[grep("_Old",names(Zebrafish_combined_GOKEGG))])
Zebrafish_Young_overlap = Get_overlap_withinGOtables(Zebrafish_combined_GOKEGG[grep("_Young",names(Zebrafish_combined_GOKEGG))])

Human_Old_overlap = Get_overlap_withinGOtables(Human_combined_GOKEGG[grep("_Old",names(Human_combined_GOKEGG))])
Human_Young_overlap = Get_overlap_withinGOtables(Human_combined_GOKEGG[grep("_Young",names(Human_combined_GOKEGG))])

Mouse_Old_overlap = Get_overlap_withinGOtables(Mouse_combined_GOKEGG[grep("_Old",names(Mouse_combined_GOKEGG))])
Mouse_Young_overlap = Get_overlap_withinGOtables(Mouse_combined_GOKEGG[grep("_Young",names(Mouse_combined_GOKEGG))])

#########
Total_List = c(list(Zebrafish_Old_overlap),list(Zebrafish_Young_overlap),list(Human_Old_overlap),list(Human_Young_overlap),list(Mouse_Old_overlap),list(Mouse_Young_overlap))
names(Total_List) = c("Zebrafish_Old_overlap","Zebrafish_Young_overlap","Human_Old_overlap","Human_Young_overlap","Mouse_Old_overlap","Mouse_Young_overlap")
#########


library(openxlsx)
wb <- createWorkbook()
# 遍历 list，把每个元素写到单独的 sheet
for (i in seq_along(Total_List)) {
  sheet_name <- names(Total_List)[i]  # 如果 list 有名字就用名字做 sheet 名
  if (is.null(sheet_name) || sheet_name == "") {
    sheet_name <- paste0("Sheet", i)  # 如果没有名字，用 Sheet1, Sheet2...
  }
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, Total_List[[i]])
}

# 保存 Excel 文件
setwd("/zp1/data/share/Human_aging_new")
saveWorkbook(wb, "HMZ_overlaps_within_species_GO_May5_2025.xlsx", overwrite = TRUE)


####### Next see the 3species overlaps #######
#######

index = names(Zebrafish_combined_GOKEGG)

total_res_list = list()
for(i in index){
	print(i)
	H = Human_combined_GOKEGG[[i]]
	M = Mouse_combined_GOKEGG[[i]]
	Z = Zebrafish_combined_GOKEGG[[i]]
	##
	inputlist = list(H=H,M=M,Z=Z)
	names(inputlist) = paste0(names(inputlist),"_",i)
	##
	res = Get_overlap_withinGOtables(inputlist)
	##
	total_res_list <- c(total_res_list,list(res))
}

names(total_res_list) = index
names(total_res_list) = paste0(names(total_res_list),"_overlap")

####
merge_overlap = c(Total_List,total_res_list)
####
							

names(Zebrafish_combined_GOKEGG) <- paste0("Zebrafish_",names(Zebrafish_combined_GOKEGG))
names(Human_combined_GOKEGG) <- paste0("Human_",names(Human_combined_GOKEGG))
names(Mouse_combined_GOKEGG) <- paste0("Mouse_",names(Mouse_combined_GOKEGG))

#####

							
merge_overlap = c(Zebrafish_combined_GOKEGG,Mouse_combined_GOKEGG,Human_combined_GOKEGG,merge_overlap)

#####


library(openxlsx)
wb <- createWorkbook()
# 遍历 list，把每个元素写到单独的 sheet
for (i in seq_along(merge_overlap)) {
  sheet_name <- names(merge_overlap)[i]  # 如果 list 有名字就用名字做 sheet 名
  if (is.null(sheet_name) || sheet_name == "") {
    sheet_name <- paste0("Sheet", i)  # 如果没有名字，用 Sheet1, Sheet2...
  }
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, merge_overlap[[i]])
}

# 保存 Excel 文件
setwd("/zp1/data/share/Human_aging_new")
saveWorkbook(wb, "total_merge_GOterms.xlsx", overwrite = TRUE)






  



