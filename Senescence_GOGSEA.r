#
####### THIS FILE IS FOR Senescence GENES ####
#

### https://cells.ucsc.edu/?ds=retina+hrca+atac ####


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R


###
setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")

install.packages("jsonlite")

library(jsonlite)

###
json_data <- fromJSON("Reactome+Pathways+2024.txt")
all_symbols <- json_data$associations$gene$symbol

####
unique_symbols <- unique(all_symbols)

Human_Reactome_senescence = unique_symbols

###### covert to zebrafish and mouse #####

setwd("/zp1/data/share/Human_aging_new")
saveRDS(HMZ_ortholog_combined,file="HMZ_ortholog_combined_2025")

HMZ_ortholog_combined = readRDS("HMZ_ortholog_combined_2025")

names(HMZ_ortholog_combined)

HM = HMZ_ortholog_combined$HM

Mouse_Reactome_senescence = HM$mouse[which(HM$human %in% Human_Reactome_senescence == T)]
Mouse_Reactome_senescence = Mouse_Reactome_senescence[!duplicated(Mouse_Reactome_senescence)]

#########
#########

HZ = HMZ_ortholog_combined$HZ

Zebrafish_Reactome_senescence = HZ$zebrafish[which(HZ$human %in% Human_Reactome_senescence == T)]
Zebrafish_Reactome_senescence = Zebrafish_Reactome_senescence[!duplicated(Zebrafish_Reactome_senescence)]

#########
#########
Reactome_senescence_list = list(H=Human_Reactome_senescence,M=Mouse_Reactome_senescence,Z=Zebrafish_Reactome_senescence)

##########

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Reactome_senescence_list,file="Reactome_senescence_list")

##########

library(readxl)

excel_path <- "SenMayo41467_2022_32552_MOESM4_ESM.xlsx"


# 4.1 读取人类基因列表所在的 sheet（名称为 "human"）
df_human <- read_excel(excel_path, sheet = "human")

# 4.2 读取小鼠基因列表所在的 sheet（名称为 "mouse"）
df_mouse <- read_excel(excel_path, sheet = "mouse")


# 6.1 提取人类基因名
human_genes_raw <- df_human$`Gene(human)`       # 或者用 df_human[["Gene(human)"]]

# 去掉 NA 或 空字符
human_genes_clean <- human_genes_raw[!is.na(human_genes_raw) &
                                     human_genes_raw != ""]

# 如果需要去重：
human_genes <- unique(human_genes_clean)

# 6.2 提取小鼠基因名
mouse_genes_raw <- df_mouse$`Gene(murine)`      # 或者用 df_mouse[["Gene(murine)"]]

# 同样去掉 NA 或 空字符
mouse_genes_clean <- mouse_genes_raw[!is.na(mouse_genes_raw) &
                                     mouse_genes_raw != ""]

# 如果需要去重：
mouse_genes <- unique(mouse_genes_clean)



####### Next for Zebrafish #########
#######
HZ = HMZ_ortholog_combined$HZ

Z_genes_1 = HZ$zebrafish[which(HZ$human %in% human_genes == T)]

MZ = HMZ_ortholog_combined$MZ

Z_genes_2 = MZ$zebrafish[which(MZ$mouse %in% mouse_genes == T)]


Z_genes = c(Z_genes_1,Z_genes_2)
Z_genes = Z_genes[!duplicated(Z_genes)]

########
########

SenMayo_senescence_list = list(H=human_genes,M=mouse_genes,Z=Z_genes)


setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(SenMayo_senescence_list,file="SenMayo_senescence_list")

#########
#########

#########
######### 只做 Old genes ######
#########






######### GO term enrichment #########
#########


######### 把 HMZ 合起来 #######
#########


#########


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
#######

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Mouse_UP_list,file="Mouse_UP_list")
save(Mouse_DOWN_list,file="Mouse_DOWN_list")


setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Human_UP_list,file="Human_UP_list")
save(Human_DOWN_list,file="Human_DOWN_list")


setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Zebrafish_UP_list,file="Zebrafish_UP_list")
save(Zebrafish_DOWN_list,file="Zebrafish_DOWN_list")



#######
####### perform GO terms enrichment ######
#######


###### zebrafish #######
gene_list <- Zebrafish_UP_list$MG$Gene


setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load(file="Reactome_senescence_list")
load(file="SenMayo_senescence_list")


Z_custom_term2gene <- data.frame(
  TERM = c(
    rep("Reactome_senescence", length(Reactome_senescence_list$Z)),
    rep("SenMayo_senescence", length(SenMayo_senescence_list$Z))
  ),
  GENE = c(
    Reactome_senescence_list$Z,
    SenMayo_senescence_list$Z
  ),
  stringsAsFactors = FALSE
)


library(clusterProfiler)

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load("Zebrafish_UP_list")

###
###
### see the universe ####
###

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_pseudo_count <- readRDS("Zebrafish_pseudo_count_2025")
Z_universe = rownames(Zebrafish_pseudo_count[[1]])


#####

gene_list <- Zebrafish_UP_list$Rod$Gene

enrich_res <- enricher(
    gene           = gene_list,
    universe       = Z_universe,
    TERM2GENE      = Z_custom_term2gene[, c("TERM","GENE")],
    pAdjustMethod  = "BH",
    pvalueCutoff   = 1,
    qvalueCutoff   = 1
)

result_df <- enrich_res@result
Rod_res = result_df

########
########

celltypes <- c("AC", "BC", "Cone", "HC", "MG", "Microglia", "RGC", "Rod", "RPE")

enrich_results <- vector("list", length(celltypes))
names(enrich_results) <- celltypes

# 循环跑每个细胞类型
for (ct in celltypes) {
  # 1. 从 Zebrafish_UP_list 中提取该细胞类型对应的基因向量
  #    假设每个子列表的结构是：Zebrafish_UP_list$<CellType>$Gene
  gene_list <- Zebrafish_UP_list[[ct]]$Gene
  
  # 2. 调用 enricher() 做富集分析
   enrich_res <- HyperGeoEnrichment(
    gene          = gene_list,
    background      = Z_universe,
    term2gene_df     = Z_custom_term2gene[, c("TERM", "GENE")],
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
    enrich_results[[ct]] <- enrich_res
}

########
library(dplyr)
combined_df <- bind_rows(enrich_results, .id = "celltype")


##########
Z_combined_df = combined_df

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Z_combined_df,file="Z_combined_df")



setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load("Z_combined_df")


###########------------------------------------------------------------------------------------------------------------------------
###########------------------------------------------------------------------------------------------------------------------------
###########------------------------------------------------------------------------------------------------------------------------
###########------------------------------------------------------------------------------------------------------------------------
########### 

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load("Mouse_UP_list")


M_custom_term2gene <- data.frame(
  TERM = c(
    rep("Reactome_senescence", length(Reactome_senescence_list$M)),
    rep("SenMayo_senescence", length(SenMayo_senescence_list$M))
  ),
  GENE = c(
    Reactome_senescence_list$M,
    SenMayo_senescence_list$M
  ),
  stringsAsFactors = FALSE
)


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_pseudo_count <- readRDS("Mouse_pseudo_count_2025")
M_universe = rownames(Mouse_pseudo_count[[1]])

celltypes <- c("AC", "BC", "Cone", "HC", "MG", "Microglia", "RGC", "Rod", "RPE")
enrich_results <- vector("list", length(celltypes))
names(enrich_results) <- celltypes

# 循环跑每个细胞类型
for (ct in celltypes) {
  # 1. 从 Zebrafish_UP_list 中提取该细胞类型对应的基因向量
  #    假设每个子列表的结构是：Zebrafish_UP_list$<CellType>$Gene
  gene_list <- Mouse_UP_list[[ct]]$Gene
  print(length(gene_list))
  print(length(intersect(gene_list,M_universe)))
  print(length(intersect(gene_list,M_custom_term2gene$GENE)))
  print(length(intersect(M_universe,M_custom_term2gene$GENE)))
  # 2. 调 . . 用 enricher() 做富集分析
  enrich_res <- HyperGeoEnrichment(
    gene          = gene_list,
    background      = M_universe,
    term2gene_df     = M_custom_term2gene[, c("TERM", "GENE")],
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  enrich_results[[ct]] <- enrich_res
}

########
library(dplyr)
combined_df <- bind_rows(enrich_results, .id = "celltype")

##########
M_combined_df = combined_df

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(M_combined_df,file="M_combined_df")



###########
###########



setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load("Human_UP_list")


H_custom_term2gene <- data.frame(
  TERM = c(
    rep("Reactome_senescence", length(Reactome_senescence_list$H)),
    rep("SenMayo_senescence", length(SenMayo_senescence_list$H))
  ),
  GENE = c(
    Reactome_senescence_list$H,
    SenMayo_senescence_list$H
  ),
  stringsAsFactors = FALSE
)

######
######
setwd("/zp1/data/share/Human_aging_new")

load(file="Human_AC_clcl2_F_Avg_2025")
load(file="Human_AC_clcl2_M_Avg_2025")

H_universe = rownames(Human_AC_clcl2_F_Avg[[1]])


celltypes <- c("AC", "BC", "Cone", "HC", "MG", "Microglia", "RGC", "Rod", "RPE")
enrich_results <- vector("list", length(celltypes))
names(enrich_results) <- celltypes

# 循环跑每个细胞类型
for (ct in celltypes) {
  # 1. 从 Zebrafish_UP_list 中提取该细胞类型对应的基因向量
  #    假设每个子列表的结构是：Zebrafish_UP_list$<CellType>$Gene
  gene_list <- Human_UP_list[[ct]]$Gene
  print(length(gene_list))
  print(length(intersect(gene_list,H_universe)))
  print(length(intersect(gene_list,H_custom_term2gene$GENE)))
  print(length(intersect(H_universe,H_custom_term2gene$GENE)))
  # 2. 调 . . 用 enricher() 做富集分析
  enrich_res <- HyperGeoEnrichment(
    gene          = gene_list,
    background      = H_universe,
    term2gene_df     = H_custom_term2gene[, c("TERM", "GENE")],
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  enrich_results[[ct]] <- enrich_res
}

########
library(dplyr)
combined_df <- bind_rows(enrich_results, .id = "celltype")

##########
H_combined_df = combined_df

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(H_combined_df,file="H_combined_df")


##########
########## GSEA ######
##########



setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load("Z_combined_df")
load("M_combined_df")
load("H_combined_df")

#####
#####
combined_df = Z_combined_df

Process_combined <- function(combined_df){
    ########
    combined_df$logP = -log10(combined_df$pvalue)
    ########
    return(combined_df)
}

########
library(ggplot2)

Z_combined_df = Process_combined(Z_combined_df)

library(ggplot2)
ggplot(Z_combined_df,aes(x=celltype,y=ID)) + geom_point(aes(size= Count,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='red', guide = guide_colorbar(reverse = FALSE, order = 1)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("HMZ_MG_overlap.png",height=4,width=8)





























HyperGeoEnrichment <- function(
  gene_list,        # 目标基因（字符向量）
  term2gene_df,     # 至少两列：TERM、GENE
  background,       # 背景基因（字符向量）
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1
) {
  # 1. 确保输入格式正确
  if (!is.character(gene_list) || !is.character(background) ||
      !is.data.frame(term2gene_df) || ncol(term2gene_df) < 2) {
    return(NULL)
  }
  
  # 2. 规范列名并去重
  colnames(term2gene_df)[1:2] <- c("TERM", "GENE")
  background <- unique(background)
  gene_mapped <- intersect(unique(gene_list), background)
  n <- length(gene_mapped)
  N <- length(background)
  if (n == 0 || N == 0) return(NULL)
  
  # 3. 按 TERM 分组
  term2genes <- split(term2gene_df$GENE, term2gene_df$TERM)
  
  # 4. 计算每个 TERM 的 k、M、pvalue
  out <- lapply(names(term2genes), function(term) {
    genes_u <- intersect(unique(term2genes[[term]]), background)
    M <- length(genes_u)
    k <- length(intersect(gene_mapped, genes_u))
    pval <- if (k == 0) 1 else phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(
      ID          = term,
      GeneRatio   = paste0(k, "/", n),
      BgRatio     = paste0(M, "/", N),
      pvalue      = pval,
      geneID      = if (k > 0) paste(intersect(gene_mapped, genes_u), collapse = "/") else "",
      Count       = k,
      stringsAsFactors = FALSE
    )
  })
  
  result <- do.call(rbind, out)
  if (nrow(result) == 0) return(NULL)
  
  # 5. p.adjust & qvalue
  result$p.adjust <- p.adjust(result$pvalue, method = pAdjustMethod)
  result$qvalue   <- result$p.adjust
  
  # 6. 阈值筛选并排序
  keep <- !is.na(result$pvalue) & result$pvalue <= pvalueCutoff &
          !is.na(result$p.adjust) & result$p.adjust <= qvalueCutoff
  result <- result[keep, ]
  if (nrow(result) == 0) return(NULL)
  
  result[order(result$p.adjust), ]
}


#########
#########
#########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

######## load the Human mouse and fish #######
######## aging genes #######
########

######## for Mouse ##
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


## fit_res = Mouse_MG_model$model

Get_G_from_model <- function(fit_res){
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
        pos <- coef_df$coefficient != 0 & coef_df$gene != '(Intercept)'
        genes_pos <- coef_df[pos,]
        return(genes_pos)
    }
    ######
    res= GetPositiveGenes(coef_df)
    res = res[order(res$coefficient,decreasing=T),]
    res
}



MG = Get_G_from_model(Mouse_MG_model$model)
RGC = Get_G_from_model(Mouse_RGC_model$model)
AC = Get_G_from_model(Mouse_AC_model$model)
HC = Get_G_from_model(Mouse_HC_model$model)
BC = Get_G_from_model(Mouse_BC_model$model)
RPE = Get_G_from_model(Mouse_RPE_model$model)
Rod = Get_G_from_model(Mouse_Rod_model$model)
Cone = Get_G_from_model(Mouse_Cone_model$model)
Microglia = Get_G_from_model(Mouse_Microglia_model$model)


######
######

Mouse_clock = list(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
names(Mouse_clock) <- c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia")
setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Mouse_clock,file="Mouse_clock")


######### next for zebrafish ########


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


MG = Get_G_from_model(Zebrafish_MG_model$model)
RGC = Get_G_from_model(Zebrafish_RGC_model$model)
AC = Get_G_from_model(Zebrafish_AC_model$model)
HC = Get_G_from_model(Zebrafish_HC_model$model)
BC = Get_G_from_model(Zebrafish_BC_model$model)
RPE = Get_G_from_model(Zebrafish_RPE_model$model)
Rod = Get_G_from_model(Zebrafish_Rod_model$model)
Cone = Get_G_from_model(Zebrafish_Cone_model$model)
Microglia = Get_G_from_model(Zebrafish_Microglia_model$model)


Zebrafish_clock = list(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
names(Zebrafish_clock) <- c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia")
setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Zebrafish_clock,file="Zebrafish_clock")


######Next for Human #####

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


MG        <- Get_G_from_model(Human_MG_model_F$model)
RGC       <- Get_G_from_model(Human_RGC_model_F$model)
AC        <- Get_G_from_model(Human_AC_model_F$model)
HC        <- Get_G_from_model(Human_HC_model_F$model)
BC        <- Get_G_from_model(Human_BC_model_F$model)
RPE       <- Get_G_from_model(Human_RPE_model_F$model)
Rod       <- Get_G_from_model(Human_Rod_model_F$model)
Cone      <- Get_G_from_model(Human_Cone_model_F$model)
Microglia <- Get_G_from_model(Human_Microglia_model_F$model)


Human_F_clock = list(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
names(Human_F_clock) <- c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia")
setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Human_F_clock,file="Human_F_clock")

MG        <- Get_G_from_model(Human_MG_model_M$model)
RGC       <- Get_G_from_model(Human_RGC_model_M$model)
AC        <- Get_G_from_model(Human_AC_model_M$model)
HC        <- Get_G_from_model(Human_HC_model_M$model)
BC        <- Get_G_from_model(Human_BC_model_M$model)
RPE       <- Get_G_from_model(Human_RPE_model_M$model)
Rod       <- Get_G_from_model(Human_Rod_model_M$model)
Cone      <- Get_G_from_model(Human_Cone_model_M$model)
Microglia <- Get_G_from_model(Human_Microglia_model_M$model)



Human_M_clock = list(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
names(Human_M_clock) <- c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia")
setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
save(Human_M_clock,file="Human_M_clock")


######
######
######


setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load(file="SenMayo_senescence_list")

setwd("/zp1/data/plyu3/Aging_add_figures/Senescence")
load("Reactome_senescence_list")

######
######


H_custom_term2gene <- data.frame(
  TERM = c(
    rep("Reactome_senescence", length(Reactome_senescence_list$H)),
    rep("SenMayo_senescence", length(SenMayo_senescence_list$H))
  ),
  GENE = c(
    Reactome_senescence_list$H,
    SenMayo_senescence_list$H
  ),
  stringsAsFactors = FALSE
)

M_custom_term2gene <- data.frame(
  TERM = c(
    rep("Reactome_senescence", length(Reactome_senescence_list$M)),
    rep("SenMayo_senescence", length(SenMayo_senescence_list$M))
  ),
  GENE = c(
    Reactome_senescence_list$M,
    SenMayo_senescence_list$M
  ),
  stringsAsFactors = FALSE
)

Z_custom_term2gene <- data.frame(
  TERM = c(
    rep("Reactome_senescence", length(Reactome_senescence_list$Z)),
    rep("SenMayo_senescence", length(SenMayo_senescence_list$Z))
  ),
  GENE = c(
    Reactome_senescence_list$Z,
    SenMayo_senescence_list$Z
  ),
  stringsAsFactors = FALSE
)

#######
#######
#######
library(clusterProfiler)

geneList =Mouse_clock$MG$gene
geneList <- sort(setNames(Mouse_clock$MG$coefficient, Mouse_clock$MG$gene), decreasing = TRUE)

gsea_res <- GSEA(
  geneList  = geneList,
  TERM2GENE = M_custom_term2gene[, c("TERM","GENE")],
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  verbose = FALSE
)
MG_gsea_res = gsea_res@result

#####
#####
#####

celltypes <- c("MG", "RGC", "AC", "HC", "BC", "RPE", "Rod", "Cone", "Microglia")

# 用于存储每个细胞类型的 GSEA 结果
gsea_results <- vector("list", length(celltypes))
names(gsea_results) <- celltypes

# 预先定义一个“无富集”时要返回的示例行
no_enrichment_row <- data.frame(
  ID              = "No_enrichment",
  Description     = "No enrichment detected",
  setSize         = 0L,
  enrichmentScore = NA_real_,
  NES             = NA_real_,
  pvalue          = NA_real_,
  p.adjust        = NA_real_,
  qvalue         = NA_real_,
  rank            = NA_integer_,
  leading_edge    = "",
  core_enrichment = "",
  stringsAsFactors = FALSE
)

for (ct in celltypes) {
  # 1. 构造 geneList：取 Mouse_clock[[ct]]$gene 和 coefficient
  genes    <- Mouse_clock[[ct]]$gene
  coeffs   <- Mouse_clock[[ct]]$coefficient
  geneList <- sort(setNames(coeffs, genes), decreasing = TRUE)
  
  # 2. 运行 GSEA
  gsea_res <- tryCatch(
    GSEA(
      geneList      = geneList,
      TERM2GENE     = M_custom_term2gene[, c("TERM", "GENE")],
      pAdjustMethod = "BH",
      pvalueCutoff  = 2,
      verbose       = FALSE
    ),
    error = function(e) NULL
  )
  
  # 3. 提取结果表格并存储
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    gsea_results[[ct]] <- gsea_res@result
  } else {
    # 返回一个只有“无富集”这一行的 data.frame
    gsea_results[[ct]] <- no_enrichment_row
  }
}

Mouse_gsea_results = do.call(rbind,gsea_results)
######
######
save(Mouse_gsea_results,file="Mouse_gsea_results")



#####

celltypes <- c("MG", "RGC", "AC", "HC", "BC", "RPE", "Rod", "Cone", "Microglia")

# 用于存储每个细胞类型的 GSEA 结果
gsea_results <- vector("list", length(celltypes))
names(gsea_results) <- celltypes

# 预先定义一个“无富集”时要返回的示例行
no_enrichment_row <- data.frame(
  ID              = "No_enrichment",
  Description     = "No enrichment detected",
  setSize         = 0L,
  enrichmentScore = NA_real_,
  NES             = NA_real_,
  pvalue          = NA_real_,
  p.adjust        = NA_real_,
  qvalue         = NA_real_,
  rank            = NA_integer_,
  leading_edge    = "",
  core_enrichment = "",
  stringsAsFactors = FALSE
)

for (ct in celltypes) {
  # 1. 构造 geneList：取 Mouse_clock[[ct]]$gene 和 coefficient
  genes    <- Zebrafish_clock[[ct]]$gene
  coeffs   <- Zebrafish_clock[[ct]]$coefficient
  geneList <- sort(setNames(coeffs, genes), decreasing = TRUE)
  
  # 2. 运行 GSEA
  gsea_res <- tryCatch(
    GSEA(
      geneList      = geneList,
      TERM2GENE     = Z_custom_term2gene[, c("TERM", "GENE")],
      pAdjustMethod = "BH",
      pvalueCutoff  = 2,
      verbose       = FALSE
    ),
    error = function(e) NULL
  )
  
  # 3. 提取结果表格并存储
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    gsea_results[[ct]] <- gsea_res@result
  } else {
    # 返回一个只有“无富集”这一行的 data.frame
    gsea_results[[ct]] <- no_enrichment_row
  }
}

Zebrafish_gsea_results = do.call(rbind,gsea_results)
######
######
save(Zebrafish_gsea_results,file="Zebrafish_gsea_results")


######
######
######
######


celltypes <- c("MG", "RGC", "AC", "HC", "BC", "RPE", "Rod", "Cone", "Microglia")

# 用于存储每个细胞类型的 GSEA 结果
gsea_results <- vector("list", length(celltypes))
names(gsea_results) <- celltypes

# 预先定义一个“无富集”时要返回的示例行
no_enrichment_row <- data.frame(
  ID              = "No_enrichment",
  Description     = "No enrichment detected",
  setSize         = 0L,
  enrichmentScore = NA_real_,
  NES             = NA_real_,
  pvalue          = NA_real_,
  p.adjust        = NA_real_,
  qvalue         = NA_real_,
  rank            = NA_integer_,
  leading_edge    = "",
  core_enrichment = "",
  stringsAsFactors = FALSE
)

for (ct in celltypes) {
  # 1. 构造 geneList：取 Mouse_clock[[ct]]$gene 和 coefficient
  genes    <- Human_F_clock[[ct]]$gene
  coeffs   <- Human_F_clock[[ct]]$coefficient
  geneList <- sort(setNames(coeffs, genes), decreasing = TRUE)
  
  # 2. 运行 GSEA
  gsea_res <- tryCatch(
    GSEA(
      geneList      = geneList,
      TERM2GENE     = H_custom_term2gene[, c("TERM", "GENE")],
      pAdjustMethod = "BH",
      pvalueCutoff  = 2,
      verbose       = FALSE
    ),
    error = function(e) NULL
  )
  
  # 3. 提取结果表格并存储
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    gsea_results[[ct]] <- gsea_res@result
  } else {
    # 返回一个只有“无富集”这一行的 data.frame
    gsea_results[[ct]] <- no_enrichment_row
  }
}

Human_F_gsea_results = do.call(rbind,gsea_results)
######
######
save(Human_F_gsea_results,file="Human_F_gsea_results")


#####
#####



celltypes <- c("MG", "RGC", "AC", "HC", "BC", "RPE", "Rod", "Cone", "Microglia")

# 用于存储每个细胞类型的 GSEA 结果
gsea_results <- vector("list", length(celltypes))
names(gsea_results) <- celltypes

# 预先定义一个“无富集”时要返回的示例行
no_enrichment_row <- data.frame(
  ID              = "No_enrichment",
  Description     = "No enrichment detected",
  setSize         = 0L,
  enrichmentScore = NA_real_,
  NES             = NA_real_,
  pvalue          = NA_real_,
  p.adjust        = NA_real_,
  qvalue         = NA_real_,
  rank            = NA_integer_,
  leading_edge    = "",
  core_enrichment = "",
  stringsAsFactors = FALSE
)

for (ct in celltypes) {
  # 1. 构造 geneList：取 Mouse_clock[[ct]]$gene 和 coefficient
  genes    <- Human_M_clock[[ct]]$gene
  coeffs   <- Human_M_clock[[ct]]$coefficient
  geneList <- sort(setNames(coeffs, genes), decreasing = TRUE)
  
  # 2. 运行 GSEA
  gsea_res <- tryCatch(
    GSEA(
      geneList      = geneList,
      TERM2GENE     = H_custom_term2gene[, c("TERM", "GENE")],
      pAdjustMethod = "BH",
      pvalueCutoff  = 2,
      verbose       = FALSE
    ),
    error = function(e) NULL
  )
  
  # 3. 提取结果表格并存储
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    gsea_results[[ct]] <- gsea_res@result
  } else {
    # 返回一个只有“无富集”这一行的 data.frame
    gsea_results[[ct]] <- no_enrichment_row
  }
}

Human_M_gsea_results = do.call(rbind,gsea_results)
######
######
save(Human_M_gsea_results,file="Human_M_gsea_results")




####### 最后输出 excel files ！！！！ ########
#######

GSEA_table <- list(Zebrafish=Zebrafish_gsea_results,Mouse=Mouse_gsea_results,Human_F=Human_F_gsea_results,Human_M=Human_M_gsea_results)

########

library(openxlsx)

wb <- createWorkbook()

# 4. 把每个 data.frame 写入一个 sheet，保留行名
for (sheet in names(GSEA_table)) {
  addWorksheet(wb, sheetName = sheet)
  writeData(
    wb,
    sheet     = sheet,
    x         = GSEA_table[[sheet]],
    rowNames  = TRUE     # 这个参数会把行名一起写入第一列
  )
}

# 5. 保存为 Excel 文件（例如 "output.xlsx"，如果已存在则覆盖）
saveWorkbook(wb, file = "Clock_GSEA.xlsx", overwrite = TRUE)



###
load("H_combined_df")
load("M_combined_df")
load("Z_combined_df")


GO_table <- list(Zebrafish=Z_combined_df,Mouse=M_combined_df,Human=H_combined_df)




library(openxlsx)

wb <- createWorkbook()

# 4. 把每个 data.frame 写入一个 sheet，保留行名
for (sheet in names(GO_table)) {
  addWorksheet(wb, sheetName = sheet)
  writeData(
    wb,
    sheet     = sheet,
    x         = GO_table[[sheet]],
    rowNames  = TRUE     # 这个参数会把行名一起写入第一列
  )
}

# 5. 保存为 Excel 文件（例如 "output.xlsx"，如果已存在则覆盖）
saveWorkbook(wb, file = "DEGs_GO.xlsx", overwrite = TRUE)



