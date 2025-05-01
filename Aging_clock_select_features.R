###########
########### load the Aging MAP and SenCID ###########----------------------------------------------------------------------------------------------------------------------------------------------------------
###########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R


setwd("/zp1/data/plyu3/Aging_Clocks_Final/")
load("Human_AgingMap_Genes_Jul9")
head(Human_AgingMap_Genes)
load("Mouse_AgingMap_Genes_Jul9")
head(Mouse_AgingMap_Genes)
load("Zebrafish_AgingMap_Genes_Jul9")
head(Zebrafish_AgingMap_Genes)
load("SID_Human_Jul9")
head(SID_Human)
load("SID_Mouse_Jul9")
head(SID_Mouse)
load("SID_Zebrafish_Jul9")
head(SID_Zebrafish)



###########
########### calculate the spearman correlations for each gene based on the pseduo-bulk matrix #########---------------------------------------------------------------------------------------------------------------------------------------------------------------------
########### 


###########----------function----------------------------------------------------------------------------------------------------------------------------------------------------------------
###########----------
library(parallel)

#' 根据 bulk_mat 和 bulk_meta 计算每个基因与 age 的 Spearman 相关
#'
#' @param list_process2 包含
#'   - bulk_mat: matrix，行是基因，列是样本
#'   - bulk_meta: data.frame，需要有列 sample 和 age
#' @return data.frame，只包括 pvalue < 0.05 的基因，包含 Gene, rho, pvalue, TAG
###  list_process2 = Mouse_pseudo_count_input$MG

Identify_Corr_Spearman <- function(list_process2) {
  bulk_mat  <- list_process2$bulk_mat
  bulk_meta <- list_process2$bulk_meta
  
  # 检查
  stopifnot(all(colnames(bulk_mat) == bulk_meta$sample))
  
  # 过滤：基因平均表达大于阈值
  avg_by_gene <- rowMeans(bulk_mat)
  keep <- avg_by_gene > 0.01
  bulk_mat <- bulk_mat[keep, , drop = FALSE]
  
  genes <- rownames(bulk_mat)
  ages  <- as.numeric(bulk_meta$age)
  
  # 并行计算：每个基因做 Spearman 相关
  library(parallel)
  ##
  ncores <- 30
  res_list <- mclapply(
    seq_along(genes),
    function(i) {
      x <- bulk_mat[i, ]
      ct <- cor.test(x, ages, method = "spearman", exact = FALSE)
      list(rho = as.numeric(ct$estimate),
           p   = ct$p.value)
    },
    mc.cores = ncores
  )
  
  # 汇总结果
  rho    <- vapply(res_list, `[[`, numeric(1), "rho")
  pvalue <- vapply(res_list, `[[`, numeric(1), "p")
  
  res_table <- data.frame(
    Gene   = genes,
    rho    = rho,
    pvalue = pvalue,
    TAG    = "NotDEGs",
    stringsAsFactors = FALSE
  )
  
  # 只保留显著基因
  sig <- res_table$pvalue < 0.1
  res_sig <- res_table[sig, ]
  res_sig$TAG <- "DEG"
  ##
  message("Total genes tested: ", nrow(res_table),
          "; Significant (p<0.1): ", nrow(res_sig))
  return(res_sig)
}


##### Corrlist = Mouse_DEGs
##### AgingMap = Mouse_AgingMap_Genes
##### SID = SID_Mouse
Add_Genes <- function(Corrlist,AgingMap,SID){
    ########
    Out_list <- list()
    for(i in 1:length(Corrlist)){
        print(i)
        ########
        tmp = Corrlist[[i]]
        Genes = tmp$Gene
        Genes_all = c(Genes,AgingMap,SID)
        Genes_all = Genes_all[!duplicated(Genes_all)]
        ########
        Out_list <- c(Out_list,list(Genes_all))
    }
    #######
    names(Out_list) = names(Corrlist)
    return(Out_list)
}



###########----------for Mouse--------------------###########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########### 

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
    res = Identify_Corr_Spearman(Mouse_pseudo_count_input[[i]])
    #####
    Mouse_DEGs <- c(Mouse_DEGs,list(res))
}

######
names(Mouse_DEGs) <- names(Mouse_pseudo_count_input)

######
Mouse_features = Add_Genes(Mouse_DEGs,Mouse_AgingMap_Genes,SID_Mouse)

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
saveRDS(Mouse_features,file="Mouse_features_2025")


###########----------for Human--------------------###########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_pseudo_count <- readRDS("Mouse_pseudo_count_2025")


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
    res = Identify_Corr_Spearman(Human_pseudo_count[[i]])
    #####
    Human_DEGs <- c(Human_DEGs,list(res))
}

######
names(Human_DEGs) <- names(Human_pseudo_count)


Human_features = Add_Genes(Human_DEGs,Human_AgingMap_Genes,SID_Human)

######
######
setwd("/zp1/data/share/Human_aging_new")
saveRDS(Human_features,file="Human_features_2025")

######
######


###
###########----------for Zebrafish--------------------###########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###

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



Zebrafish_DEGs <- list()

for(i in 1:length(Zebrafish_pseudo_count_input)){
    print(i)
    #####
    res = Identify_Corr_Spearman(Zebrafish_pseudo_count_input[[i]])
    #####
    Zebrafish_DEGs <- c(Zebrafish_DEGs,list(res))
}

######
names(Zebrafish_DEGs) <- names(Zebrafish_pseudo_count_input)
Zebrafish_features = Add_Genes(Zebrafish_DEGs,Zebrafish_AgingMap_Genes,SID_Zebrafish)

#####

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
saveRDS(Zebrafish_features,file="Zebrafish_features_2025")


######
###### Done !!!! ######
###### Done !!!! ######
######













###########-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





















