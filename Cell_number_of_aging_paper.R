###########
########### we will see how many donors we have for human samples ######
###########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R

######
setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_seurat <- readRDS("Zebrafish_seurat_merge_clcl_2025")

table(Zebrafish_seurat$sample)
#####
Zebrafish_seurat_cl = Zebrafish_seurat[,which(Zebrafish_seurat$sample %in% c("28_Rep1","28_Rep2") == F)]
table(Zebrafish_seurat_cl$sample)
library(Seurat)

#####
##### 计算每个时间点，每个cell type 各有多少细胞 #######
#####

head(Zebrafish_seurat_cl@meta.data)

#####
Zebrafish_seurat_cl$age = sapply(strsplit(Zebrafish_seurat_cl$sample,split="_"),function(x) x[[1]])

#####
count_cells_by_celltype_age <- function(
  seurat_obj,
  celltype_column = NULL,
  age_column = "age"
) {
  # Extract metadata
  meta <- seurat_obj@meta.data

  # Determine cell type vector
  if (is.null(celltype_column)) {
    cell_types <- as.character(Idents(seurat_obj))
  } else {
    if (!celltype_column %in% colnames(meta)) {
      stop(sprintf("Column '%s' not found in meta.data", celltype_column))
    }
    cell_types <- as.character(meta[[celltype_column]])
  }

  # Check age column
  if (!age_column %in% colnames(meta)) {
    stop(sprintf("Column '%s' not found in meta.data", age_column))
  }
  ages <- as.character(meta[[age_column]])

  # Create contingency table and convert to matrix
  tbl <- table(cell_types, ages)
  mat <- as.matrix(tbl)

  # Sort columns by numeric age values
  age_vals <- as.numeric(colnames(mat))
  if (any(is.na(age_vals))) {
    warning("Some age column names are not numeric. Columns will remain in original order.")
  } else {
    ord <- order(age_vals)
    mat <- mat[, ord, drop = FALSE]
  }

  return(mat)
} 

                                 
mat <- count_cells_by_celltype_age(
  Zebrafish_seurat_cl,
  celltype_column = "celltype",
  age_column = "age"
)

openxlsx::write.xlsx(mat, file = "Fish_Cell_Counts.xlsx", rowNames = TRUE)

length(which(Zebrafish_seurat_cl$age == "6" & Zebrafish_seurat_cl$celltype=="RPE"))



  
#######
####### Next for the mouse #######
#######



setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_seurat <- readRDS("Mouse_seurat_merge_clcl_2025")
Mouse_meta <- readRDS("Mouse_meta_2025.rds")
celltype = Mouse_meta$celltype[match(colnames(Mouse_seurat),Mouse_meta$X)]
all.equal(as.character(Mouse_seurat$celltype),celltype)

head(Mouse_seurat@meta.data)

Mouse_seurat$age = sapply(strsplit(Mouse_seurat$sample,split="_"),function(x) x[[1]])
Mouse_seurat$age = gsub("106","108",  Mouse_seurat$age)
Mouse_seurat$age = gsub("Wk","",  Mouse_seurat$age)

                          
mat <- count_cells_by_celltype_age(
  Mouse_seurat,
  celltype_column = "celltype",
  age_column = "age"
)

openxlsx::write.xlsx(mat, file = "Mouse_Cell_Counts.xlsx", rowNames = TRUE)

length(which(Zebrafish_seurat_cl$age == "6" & Zebrafish_seurat_cl$celltype=="RPE"))




                                 
setwd("/zp1/data/plyu3/Mouse_RPE_Jul12/RNA")
Mouse_RPE_merge_cl3 <- readRDS("Mouse_RPE_RNA_Final_Jul16")

cell = which(Mouse_RPE_merge_cl3$celltype == "RPE")

                          
####### 16293 RPE + 152,796
#######
#######


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

######
######


H_ambient_genes <- readRDS("/zp1/data/plyu3/Aging_Clocks_Final/Human_ambient_RNA_Jul29")
H_ambient_genes <- H_ambient_genes[order(H_ambient_genes$Avg,decreasing=T),]
H_ambient_genes_Top50 = H_ambient_genes$Gene[which(H_ambient_genes$Avg > 0.005)]
H_ambient_genes_Top50 = sapply(strsplit(H_ambient_genes_Top50,split="~~"),function(x) x[[2]])


######
######                        

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

#####
                               
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

Human_RPE_clcl1 = Add_interval_Human(Human_RPE_clcl1)
Human_Rod_clcl1 = Add_interval_Human(Human_Rod_clcl1)
Human_Cone_clcl1 = Add_interval_Human(Human_Cone_clcl1)
Human_BC_clcl1 = Add_interval_Human(Human_BC_clcl1)
Human_AC_clcl1 = Add_interval_Human(Human_AC_clcl1)
Human_RGC_clcl1 = Add_interval_Human(Human_RGC_clcl1)
Human_MG_clcl1 = Add_interval_Human(Human_MG_clcl1)
Human_HC_clcl1 = Add_interval_Human(Human_HC_clcl1)
Human_Microglia_clcl1 = Add_interval_Human(Human_Microglia_clcl1)


                      
######
######
######

seurat_list <- list(
  AC       = Human_AC_clcl1,
  HC       = Human_HC_clcl1,
  BC       = Human_BC_clcl1,
  Rod      = Human_Rod_clcl1,
  Cone     = Human_Cone_clcl1,
  RPE      = Human_RPE_clcl1,
  Microglia= Human_Microglia_clcl1,
  MG       = Human_MG_clcl1,
  RGC      = Human_RGC_clcl1
)


count_cells_by_celltype_age_list <- function(
  seurat_list,
  age_column = "age",
  output_path = NULL
) {
  # Validate input
  if (!is.list(seurat_list) || length(seurat_list) == 0) {
    stop("seurat_list must be a non-empty named list of Seurat objects")
  }
  if (is.null(names(seurat_list)) || any(names(seurat_list) == "")) {
    stop("seurat_list must be a named list: each element name will be used as the row name (cell type)")
  }

  # Collect all age values across the list
  ages_all <- unlist(lapply(seurat_list, function(obj) {
    meta <- obj@meta.data
    if (!age_column %in% colnames(meta)) {
      stop(sprintf("Column '%s' not found in meta.data of one Seurat object", age_column))
    }
    as.character(meta[[age_column]])
  }))
  unique_ages <- unique(ages_all)

  # Sort ages by numeric if possible, else lexicographically
  age_nums <- suppressWarnings(as.numeric(unique_ages))
  if (all(!is.na(age_nums))) {
    sorted_ages <- unique_ages[order(age_nums)]
  } else {
    sorted_ages <- sort(unique_ages)
    warning("Age values not all numeric; sorting lexicographically")
  }

  # Initialize matrix
  mat <- matrix(
    0,
    nrow = length(seurat_list),
    ncol = length(sorted_ages),
    dimnames = list(names(seurat_list), sorted_ages)
  )

  # Fill matrix with counts
  for (ct in names(seurat_list)) {
    obj <- seurat_list[[ct]]
    print(dim(obj))
    meta <- obj@meta.data
    ages_vec <- as.character(meta[[age_column]])
    tbl <- table(ages_vec)
    for (a in names(tbl)) {
      if (a %in% colnames(mat)) {
        mat[ct, a] <- tbl[[a]]
      }
    }
  }
  ####
  return(mat)
}


                      
                          
mat <- count_cells_by_celltype_age_list(
  seurat_list,
  age_column = "age2"
)

openxlsx::write.xlsx(mat, file = "Human_Cell_Counts.xlsx", rowNames = TRUE)


                      
#####



                      
# 2. 统计每个对象的细胞数
cell_counts <- sapply(seurat_list, function(obj) {
  # 方法1：直接取列数
  ncol(obj)
  # 或者：length(Cells(obj))
  # length(Cells(obj))
})

# 3. 输出各个簇的细胞数
print(cell_counts)

# 4. 计算并输出总细胞数
total_cells <- sum(cell_counts)
cat("Total cells across all objects:", total_cells, "\n")



seurat_list <- list(
  AC       = Human_AC_clcl,
  HC       = Human_HC_clcl,
  BC       = Human_BC_clcl,
  Rod      = Human_Rod_clcl,
  Cone     = Human_Cone_clcl,
  #RPE      = Human_RPE_clcl,
  Microglia= Human_Microglia_clcl,
  MG       = Human_MG_clcl,
  RGC      = Human_RGC_clcl
)


donor_list <- lapply(seurat_list, function(obj) {
  if ("donor_id" %in% colnames(obj@meta.data)) {
    unique(obj@meta.data$donor_id)
  } else {
    NA
  }
})

# 3. 打印每个簇的 donor_id
for (name in names(donor_list)) {
  cat(paste0(name, " donors: "), paste0(donor_list[[name]], collapse = ", "), "\n")
}

# 4. 所有对象的 donor_id 并集
all_donors <- sort(unique(unlist(donor_list)))
cat("All donors across objects:", paste0(all_donors, collapse = ","), "\n")


######
###### Next we will find the Mouse RPE datasets ########
######










