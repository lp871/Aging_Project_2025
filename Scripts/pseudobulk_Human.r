####
#########
#########
#########
####

####
#### Human expression is not the raw counts ###########
#### 

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/share/Human_aging_new")
library(Seurat)

Human_AC_clcl <- readRDS("Human_AC_clcl")
Human_HC_clcl <- readRDS("Human_HC_clcl")
Human_BC_clcl <- readRDS("Human_BC_clcl")
Human_Rod_clcl <- readRDS("Human_Rod_clcl")
Human_Cone_clcl <- readRDS("Human_Cone_clcl")
Human_RPE_clcl <- readRDS("Human_RPE_clcl")
Human_Microglia_clcl <- readRDS("Human_Microglia_clcl")
Human_MG_clcl <- readRDS("Human_MG_clcl")
Human_RGC_clcl <- readRDS("Human_RGC_clcl")

#######
####### 分别 把 男女的 sample 分开 ##########
#######


#######
#######

H_ambient_genes <- readRDS("/zp1/data/plyu3/Aging_Clocks_Final/Human_ambient_RNA_Jul29")
H_ambient_genes <- H_ambient_genes[order(H_ambient_genes$Avg,decreasing=T),]
H_ambient_genes_Top50 = H_ambient_genes$Gene[which(H_ambient_genes$Avg > 0.005)]
H_ambient_genes_Top50 = sapply(strsplit(H_ambient_genes_Top50,split="~~"),function(x) x[[2]])

#######
####### clean ambient genes ##############
####### and remove mt genes ##############
#######

Clean_By_Gene_Names_Human <- function(Seurat_Obj,ambient_genes_Top50){
    ######## remove ambient genes ######
    k2 = which(rownames(Seurat_Obj) %in% ambient_genes_Top50 == T)
    mat_cl1 = Seurat_Obj[-k2,]
    ####### remove mt genes #####
    k3 = grep("^MT-",rownames(mat_cl1))
    print(length(k3))
    mat_cl2 = mat_cl1[-k3,]
    ######
    output_seurat = mat_cl2
    ###### clean metadata #########
    output_seurat@meta.data = Seurat_Obj@meta.data[,c("celltype","sample","age","gender","donor")]
    ######
    return(output_seurat)
}


####### remove samples with low cell counts for that celltype ########
####### and add aging interval for human #############################
####### convert > >89 years to 91 years #######

# colnames(Human_AC_clcl@meta.data)
# table(Human_AC_clcl$donor_id)
# table(Human_AC_clcl$donor_age)

# tmp = Human_AC_clcl[,which(Human_AC_clcl$donor_age == ">89 years")]
# table(tmp$donor_id)
# tmp = Human_AC_clcl[,which(Human_AC_clcl$donor_age == "71 years")]
# table(tmp$donor_id)

####### process_not RPE Human_HC_clcl with metadata #########
####### seurat_obj = Human_HC_clcl
sample = c(">89 years","pooled [>89 years,65 years]","pooled [>89 years,81 years]","pooled [43 years,16 years]","pooled [47 years,30 years]","pooled [67 years,71 years]","pooled [76 years,77 years]","pooled [77 years,74 years]","pooled [80 years,83 years,76 years]","pooled [81 years,80 years,79 years]")
sample2 = c("92 years","78 years","87 years","30 years","39 years","69 years","77 years","76 years","80 years","80 years")
tab = data.frame(sample=sample,sample2=sample2)


Process_meta_NORPE <- function(seurat_obj,tab){
    #########
    ### colnames(seurat_obj)
    seurat_obj$celltype = seurat_obj$majorclass
    print(as.character(seurat_obj$majorclass[1]))
    #########
    seurat_obj$sample = seurat_obj$sampleid
    #########
    seurat_obj$age = seurat_obj$donor_age
    #########
    for(i in 1:length(tab$sample)){
        k = which(seurat_obj$age == tab$sample[i])
        if(length(k) > 1){
            seurat_obj$age[k] = tab$sample2[i]
        }
    }
    ######
    print(table(seurat_obj$age))
    ####
    seurat_obj$gender = seurat_obj$sex
    ####
    seurat_obj$donor = seurat_obj$donor_id
    ###
    return(seurat_obj)
}



Human_AC_clcl = Process_meta_NORPE(Human_AC_clcl,tab)
Human_AC_clcl1 = Clean_By_Gene_Names_Human(Human_AC_clcl,H_ambient_genes_Top50)

#

Human_BC_clcl = Process_meta_NORPE(Human_BC_clcl,tab)
Human_BC_clcl1 = Clean_By_Gene_Names_Human(Human_BC_clcl,H_ambient_genes_Top50)

#
Human_HC_clcl = Process_meta_NORPE(Human_HC_clcl,tab)
Human_HC_clcl1 = Clean_By_Gene_Names_Human(Human_HC_clcl,H_ambient_genes_Top50)

#
Human_Rod_clcl = Process_meta_NORPE(Human_Rod_clcl,tab)
Human_Rod_clcl1 = Clean_By_Gene_Names_Human(Human_Rod_clcl,H_ambient_genes_Top50)


Human_Cone_clcl = Process_meta_NORPE(Human_Cone_clcl,tab)
Human_Cone_clcl1 = Clean_By_Gene_Names_Human(Human_Cone_clcl,H_ambient_genes_Top50)


#######---------
#######---------

Process_meta_RPE <- function(seurat_obj,tab){
    #########
    ### colnames(seurat_obj)
    seurat_obj$celltype = "RPE"
    #########
    seurat_obj$sample = seurat_obj$sampleid
    #########
    seurat_obj$age = seurat_obj$Age
    #########
    print(table(seurat_obj$age))
    ####
    seurat_obj$gender = seurat_obj$Gender
    ####
    seurat_obj$donor = seurat_obj$Donor
    ###
    return(seurat_obj)
}


Human_RPE_clcl = Process_meta_RPE(Human_RPE_clcl)
Human_RPE_clcl1 = Clean_By_Gene_Names_Human(Human_RPE_clcl,H_ambient_genes_Top50)


#######---------

Human_Microglia_clcl = Process_meta_NORPE(Human_Microglia_clcl,tab)
Human_Microglia_clcl1 = Clean_By_Gene_Names_Human(Human_Microglia_clcl,H_ambient_genes_Top50)


#######---------
Human_MG_clcl = Process_meta_NORPE(Human_MG_clcl,tab)
Human_MG_clcl1 = Clean_By_Gene_Names_Human(Human_MG_clcl,H_ambient_genes_Top50)


#######---------
Human_RGC_clcl = Process_meta_NORPE(Human_RGC_clcl,tab)
Human_RGC_clcl1 = Clean_By_Gene_Names_Human(Human_RGC_clcl,H_ambient_genes_Top50)

#######
#######
#######

filterSeuratByDonor <- function(seu, donor.col = "donor", min.cells = 100) {
  # 确保元数据中有指定的 donor 列
  # 1. 统计每个 donor 的细胞数
  donor_counts <- table(seu@meta.data[[donor.col]])
  # 2. 找出满足阈值的 donor
  donors_to_keep <- names(donor_counts)[donor_counts >= min.cells]
  # 3. 子集化 Seurat 对象
  cells_to_keep <- rownames(seu@meta.data)[seu@meta.data[[donor.col]] %in% donors_to_keep]
  #
  seu_filtered <- subset(seu, cells = cells_to_keep)
  print(table(seu_filtered@meta.data[[donor.col]]))
  return(seu_filtered)
}

######
###### convert ages to numeric ###########
###### process1 = Human_Rod_clcl2

Add_interval_Human <- function(process1){
    ######## rm unknown gender #####
    k = which(process1$gender == "unknown")
    if(length(k) > 0){
        process1 = process1[,-k]
    }
    ########
    meta = process1@meta.data
    meta$age = sapply(strsplit(as.character(meta$age),split=" ",fixed=T),function(x) x[[1]])
    print(table(meta$age))
    ######## summary(meta$age)
    k1 = which(meta$age >= 0 & meta$age < 20)
    k1.5 = which(meta$age >= 20 & meta$age < 30)
    k2 = which(meta$age >= 30 & meta$age < 45)
    k3 = which(meta$age >= 45 & meta$age < 55)
    k5 = which(meta$age >= 55 & meta$age < 65)
    k7 = which(meta$age >= 65 & meta$age < 75)
    k8 = which(meta$age >= 75 & meta$age < 85)
    k9 = which(meta$age >= 85 & meta$age < 95)
    #####
    meta$age2 = 0
    meta$age2[k1] = 10
    meta$age2[k1.5] = 25
    meta$age2[k2] = 37.5
    meta$age2[k3] = 50
    meta$age2[k5] = 60
    meta$age2[k7] = 70
    meta$age2[k8] = 80
    meta$age2[k9] = 90
    ########
    meta$age3 = paste(meta$gender,meta$age2)
    ########
    process1@meta.data = meta
    print(table(meta$age3))
    ########
    return(process1)
}

####### 
Human_RPE_clcl2 = filterSeuratByDonor(Human_RPE_clcl1,donor.col = "donor",min.cells = 100)
Human_Rod_clcl2 = filterSeuratByDonor(Human_Rod_clcl1,donor.col = "donor",min.cells = 100)
Human_Cone_clcl2 = filterSeuratByDonor(Human_Cone_clcl1,donor.col = "donor",min.cells = 100)
Human_BC_clcl2 = filterSeuratByDonor(Human_BC_clcl1,donor.col = "donor",min.cells = 100)
Human_AC_clcl2 = filterSeuratByDonor(Human_AC_clcl1,donor.col = "donor",min.cells = 100)
Human_RGC_clcl2 = filterSeuratByDonor(Human_RGC_clcl1,donor.col = "donor",min.cells = 100)
Human_MG_clcl2 = filterSeuratByDonor(Human_MG_clcl1,donor.col = "donor",min.cells = 100)
Human_HC_clcl2 = filterSeuratByDonor(Human_HC_clcl1,donor.col = "donor",min.cells = 100)
Human_Microglia_clcl2 = filterSeuratByDonor(Human_Microglia_clcl1,donor.col = "donor",min.cells = 10)

######

Human_RPE_clcl2 = Add_interval_Human(Human_RPE_clcl2)
Human_Rod_clcl2 = Add_interval_Human(Human_Rod_clcl2)
Human_Cone_clcl2 = Add_interval_Human(Human_Cone_clcl2)
Human_BC_clcl2 = Add_interval_Human(Human_BC_clcl2)
Human_AC_clcl2 = Add_interval_Human(Human_AC_clcl2)
Human_RGC_clcl2 = Add_interval_Human(Human_RGC_clcl2)
Human_MG_clcl2 = Add_interval_Human(Human_MG_clcl2)
Human_HC_clcl2 = Add_interval_Human(Human_HC_clcl2)
Human_Microglia_clcl2 = Add_interval_Human(Human_Microglia_clcl2)

########
########

saveRDS(Human_RPE_clcl2,file="Human_RPE_clcl2")

########
########

saveRDS(Human_Rod_clcl2,file="Human_Rod_clcl2")
saveRDS(Human_Cone_clcl2,file="Human_Cone_clcl2")
saveRDS(Human_MG_clcl2,file="Human_MG_clcl2")
saveRDS(Human_AC_clcl2,file="Human_AC_clcl2")
saveRDS(Human_BC_clcl2,file="Human_BC_clcl2")
saveRDS(Human_HC_clcl2,file="Human_HC_clcl2")
saveRDS(Human_Microglia_clcl2,file="Human_Microglia_clcl2")
saveRDS(Human_RGC_clcl2,file="Human_RGC_clcl2")

#######
#######------------- calculate the average expression for each sample ------------##########
####### for each donor, averaged their sampleid, for each age, averaged the donor expression ###########
#######
#######

####### Pseudobulk aggregation and normalization:
####### We extracted raw counts from the Seurat RNA assay (slot “counts”) as a sparse DGCM­atrix. Counts were library‐size normalized to CP10K (10 000 / cell total counts) via sparse diagonal scaling. Cells were averaged by sample to produce sample‐level pseudobulk, then sample profiles were averaged by donor and, finally, donor profiles were averaged by age to yield time‐point matrices. All steps use sparse‐matrix operations and vectorized row‐means (Matrix R package) to minimize memory and maximize speed.
#######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

setwd("/zp1/data/share/Human_aging_new")
library(Seurat)

####
#### seurat_obj = Human_HC 
#### 就合并每一个 donor ########
####


Prepare_Avg_Meta_ByDonor_Fast2 <- function(Seurat_Obj){
  library(Matrix)
  # 1. 原始稀疏 counts 和 meta
  raw_counts <- Seurat_Obj[["RNA"]]$counts
  meta       <- Seurat_Obj@meta.data
  # 2. size factors
  sf         <- colSums(raw_counts) / 1e4
  # 3. 列归一化（稀疏乘 Diagonal）
  raw_norm   <- raw_counts %*% Diagonal(x = 1/sf)
  # 4. cell→donor 矩阵 D
  donors         <- meta$donor
  uniq_donors    <- unique(donors)
  donor_idx      <- match(donors, uniq_donors)
  cells_per_donor<- tabulate(donor_idx)
  D <- sparseMatrix(
    i = seq_along(donor_idx),
    j = donor_idx,
    x = 1 / cells_per_donor[donor_idx],
    dims = c(ncol(raw_norm), length(uniq_donors))
  )
  # 5. donor→age2 矩阵 A
  donor_ages     <- meta$age2[match(uniq_donors, donors)]
  uniq_ages      <- unique(donor_ages)
  age_idx        <- match(donor_ages, uniq_ages)
  donors_per_age <- tabulate(age_idx)
  A <- sparseMatrix(
    i = seq_along(age_idx),
    j = age_idx,
    x = 1 / donors_per_age[age_idx],
    dims = c(length(uniq_donors), length(uniq_ages))
  )
  # 6. 三步稀疏乘法
  mat_avg_sparse <- raw_norm %*% D %*% A
  # 7. 加上 dimnames，然后转 dense 并 log2(x+1)
  dimnames(mat_avg_sparse) <- list(
    rownames(raw_counts),
    uniq_ages
  )
  mat_log2 <- log2(as.matrix(mat_avg_sparse) + 1)
  # 8. 构造 metadata
  bulk_meta <- data.frame(
    sample = uniq_ages,
    age    = uniq_ages,
    stringsAsFactors = FALSE
  )
  #####
  list(
    bulk_mat  = mat_log2,   # 现在有 colnames 了
    bulk_meta = bulk_meta
  )
}

#######
#######
#######
#######

Human_AC_clcl2_F = Human_AC_clcl2[,which(Human_AC_clcl2$gender == 'female')]
Human_AC_clcl2_M = Human_AC_clcl2[,which(Human_AC_clcl2$gender == 'male')]

Human_AC_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_AC_clcl2_F)
Human_AC_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_AC_clcl2_M)

save(Human_AC_clcl2_F_Avg,file="Human_AC_clcl2_F_Avg_2025")
save(Human_AC_clcl2_M_Avg,file="Human_AC_clcl2_M_Avg_2025")

####
####

Human_HC_clcl2_F = Human_HC_clcl2[,which(Human_HC_clcl2$gender == 'female')]
Human_HC_clcl2_M = Human_HC_clcl2[,which(Human_HC_clcl2$gender == 'male')]

Human_HC_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_HC_clcl2_F)
Human_HC_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_HC_clcl2_M)

save(Human_HC_clcl2_F_Avg,file="Human_HC_clcl2_F_Avg_2025")
save(Human_HC_clcl2_M_Avg,file="Human_HC_clcl2_M_Avg_2025")


####
####

Human_BC_clcl2_F = Human_BC_clcl2[,which(Human_BC_clcl2$gender == 'female')]
Human_BC_clcl2_M = Human_BC_clcl2[,which(Human_BC_clcl2$gender == 'male')]

Human_BC_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_BC_clcl2_F)
Human_BC_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_BC_clcl2_M)

save(Human_BC_clcl2_F_Avg,file="Human_BC_clcl2_F_Avg_2025")
save(Human_BC_clcl2_M_Avg,file="Human_BC_clcl2_M_Avg_2025")
######


Human_Rod_clcl2_F = Human_Rod_clcl2[,which(Human_Rod_clcl2$gender == 'female')]
Human_Rod_clcl2_M = Human_Rod_clcl2[,which(Human_Rod_clcl2$gender == 'male')]

Human_Rod_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_Rod_clcl2_F)
Human_Rod_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_Rod_clcl2_M)

save(Human_Rod_clcl2_F_Avg,file="Human_Rod_clcl2_F_Avg_2025")
save(Human_Rod_clcl2_M_Avg,file="Human_Rod_clcl2_M_Avg_2025")

#######


Human_Cone_clcl2_F = Human_Cone_clcl2[,which(Human_Cone_clcl2$gender == 'female')]
Human_Cone_clcl2_M = Human_Cone_clcl2[,which(Human_Cone_clcl2$gender == 'male')]

Human_Cone_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_Cone_clcl2_F)
Human_Cone_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_Cone_clcl2_M)

save(Human_Cone_clcl2_F_Avg,file="Human_Cone_clcl2_F_Avg_2025")
save(Human_Cone_clcl2_M_Avg,file="Human_Cone_clcl2_M_Avg_2025")


#####

Human_RPE_clcl2_F = Human_RPE_clcl2[,which(Human_RPE_clcl2$gender == 'Female')]
Human_RPE_clcl2_M = Human_RPE_clcl2[,which(Human_RPE_clcl2$gender == 'Male')]

Human_RPE_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_RPE_clcl2_F)
Human_RPE_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_RPE_clcl2_M)

save(Human_RPE_clcl2_F_Avg,file="Human_RPE_clcl2_F_Avg_2025")
save(Human_RPE_clcl2_M_Avg,file="Human_RPE_clcl2_M_Avg_2025")


####
####

Human_Microglia_clcl2_F = Human_Microglia_clcl2[,which(Human_Microglia_clcl2$gender == 'female')]
Human_Microglia_clcl2_M = Human_Microglia_clcl2[,which(Human_Microglia_clcl2$gender == 'male')]

Human_Microglia_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_Microglia_clcl2_F)
Human_Microglia_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_Microglia_clcl2_M)

save(Human_Microglia_clcl2_F_Avg,file="Human_Microglia_clcl2_F_Avg_2025")
save(Human_Microglia_clcl2_M_Avg,file="Human_Microglia_clcl2_M_Avg_2025")
######



Human_MG_clcl2_F = Human_MG_clcl2[,which(Human_MG_clcl2$gender == 'female')]
Human_MG_clcl2_M = Human_MG_clcl2[,which(Human_MG_clcl2$gender == 'male')]

Human_MG_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_MG_clcl2_F)
Human_MG_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_MG_clcl2_M)

save(Human_MG_clcl2_F_Avg,file="Human_MG_clcl2_F_Avg_2025")
save(Human_MG_clcl2_M_Avg,file="Human_MG_clcl2_M_Avg_2025")


#####

Human_RGC_clcl2_F = Human_RGC_clcl2[,which(Human_RGC_clcl2$gender == 'female')]
Human_RGC_clcl2_M = Human_RGC_clcl2[,which(Human_RGC_clcl2$gender == 'male')]

Human_RGC_clcl2_F_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_RGC_clcl2_F)
Human_RGC_clcl2_M_Avg = Prepare_Avg_Meta_ByDonor_Fast2(Human_RGC_clcl2_M)

save(Human_RGC_clcl2_F_Avg,file="Human_RGC_clcl2_F_Avg_2025")
save(Human_RGC_clcl2_M_Avg,file="Human_RGC_clcl2_M_Avg_2025")


####
#### seurat_obj = Human_Cone_clcl2
####
#### looks library size around 4200 #####
####

















