#########-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########
######### Output DEGs tables for Mouse Zebrafish and Human ##########
#########
#########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
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
library(writexl)
# Write to "output.xlsx" — each list element becomes a sheet named after its list name
write_xlsx(Mouse_combined, path = "Mouse_Aging_DEGs_May5_2025.xlsx")


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


write_xlsx(Mouse_combined_GOKEGG, path = "Mouse_Aging_DEGs_GOKEGG_May5_2025.xlsx")


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
library(writexl)
# Write to "output.xlsx" — each list element becomes a sheet named after its list name
write_xlsx(Zebrafish_combined, path = "Zebrafish_Aging_DEGs_May5_2025.xlsx")


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


write_xlsx(Zebrafish_combined_GOKEGG, path = "Zebrafish_Aging_DEGs_GOKEGG_May5_2025.xlsx")





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
library(writexl)
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
write_xlsx(Human_combined_GOKEGG, path = "Human_Aging_DEGs_GOKEGG_May5_2025.xlsx")




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



