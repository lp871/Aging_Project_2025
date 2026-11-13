######

library(biomaRt)

#' Retrieve and split ortholog gene pairs across three species
#'
#' @return A list of 4 data.frames
getOrthologListsHuman <- function() {
  # 1. Connect to the human Mart
  hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # 2. Define attributes to retrieve: human symbol, mouse ortholog symbol, zebrafish ortholog symbol
  attrs <- c(
    "external_gene_name",                   # human symbol
    "mmusculus_homolog_associated_gene_name", # mouse ortholog symbol
    "drerio_homolog_associated_gene_name"   # zebrafish ortholog symbol
  )
  df <- getBM(attributes = attrs, mart = hs)
  colnames(df) <- c("human", "mouse", "zebrafish")
  
  # 3. Remove duplicates
  df <- unique(df)
  
  # 4. Extract each set of ortholog pairs
  human_mouse     <- subset(df, mouse != "", select = c("human", "mouse"))
  human_zebrafish <- subset(df, zebrafish != "", select = c("human", "zebrafish"))
  mouse_zebrafish <- subset(df, mouse != "" & zebrafish != "", select = c("mouse", "zebrafish"))
  triple          <- subset(df, mouse != "" & zebrafish != "", select = c("human", "mouse", "zebrafish"))
  
  # 5. Return a list
  list(
    human_mouse     = human_mouse,
    human_zebrafish = human_zebrafish,
    mouse_zebrafish = mouse_zebrafish,
    triple          = triple
  )
}

########
########



getOrthologListsMouse <- function() {
  # 1. Connect to the mouse Mart
  mm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  # 2. Define attributes to retrieve: mouse symbol, human ortholog symbol, zebrafish ortholog symbol
  attrs <- c(
    "external_gene_name",                    # mouse symbol
    "hsapiens_homolog_associated_gene_name",# human ortholog symbol
    "drerio_homolog_associated_gene_name"   # zebrafish ortholog symbol
  )
  df <- getBM(attributes = attrs, mart = mm)
  colnames(df) <- c("mouse", "human", "zebrafish")
  
  # 3. Remove duplicates
  df <- unique(df)
  
  # 4. Extract each set of ortholog pairs
  human_mouse     <- subset(df, human != "" & mouse != "", select = c("human", "mouse"))
  human_zebrafish <- subset(df, human != "" & zebrafish != "", select = c("human", "zebrafish"))
  mouse_zebrafish <- subset(df, mouse != "" & zebrafish != "", select = c("mouse", "zebrafish"))
  triple          <- subset(df, mouse != "" & human != "" & zebrafish != "",
                             select = c("human", "mouse", "zebrafish"))
  
  # 5. Return a list
  list(
    human_mouse     = human_mouse,
    human_zebrafish = human_zebrafish,
    mouse_zebrafish = mouse_zebrafish,
    triple          = triple
  )
}




getOrthologListsZebrafish <- function() {
  # 1. Connect to the zebrafish Mart
  dr <- useMart("ensembl", dataset = "drerio_gene_ensembl")
  
  # 2. Define attributes to retrieve: zebrafish symbol, human ortholog symbol, mouse ortholog symbol
  attrs <- c(
    "external_gene_name",                      # zebrafish symbol
    "hsapiens_homolog_associated_gene_name",  # human ortholog symbol
    "mmusculus_homolog_associated_gene_name"  # mouse ortholog symbol
  )
  df <- getBM(attributes = attrs, mart = dr)
  colnames(df) <- c("zebrafish", "human", "mouse")
  
  # 3. Remove duplicates
  df <- unique(df)
  
  # 4. Extract each set of ortholog pairs
  human_mouse     <- subset(df, human != "" & mouse != "", select = c("human", "mouse"))
  human_zebrafish <- subset(df, human != "" & zebrafish != "", select = c("human", "zebrafish"))
  mouse_zebrafish <- subset(df, mouse != "" & zebrafish != "", select = c("mouse", "zebrafish"))
  triple          <- subset(df, mouse != "" & human != "" & zebrafish != "",
                             select = c("human", "mouse", "zebrafish"))
  
  # 5. Return a list
  list(
    human_mouse     = human_mouse,
    human_zebrafish = human_zebrafish,
    mouse_zebrafish = mouse_zebrafish,
    triple          = triple
  )
}

#########
#########
#########

#########

Hres = getOrthologListsHuman()
Mres = getOrthologListsMouse()
Zres = getOrthologListsZebrafish()

#########
names(Hres)
names(Mres)
names(Zres)

##########

human_mouse = rbind(Hres[[1]],Mres[[1]],Zres[[1]])
human_mouse$HM_index = paste(human_mouse$human,human_mouse$mouse,sep="__")
human_mouse = human_mouse[!duplicated(human_mouse$HM_index),]

dim(Hres[[1]])
dim(Mres[[1]])
dim(Zres[[1]])
dim(human_mouse)


human_zebrafish = rbind(Hres[[2]],Mres[[2]],Zres[[2]])
human_zebrafish$HZ_index = paste(human_zebrafish$human,human_zebrafish$zebrafish,sep="__")
human_zebrafish = human_zebrafish[!duplicated(human_zebrafish$HZ_index),]

dim(Hres[[2]])
dim(Mres[[2]])
dim(Zres[[2]])
dim(human_zebrafish)

#########

mouse_zebrafish = rbind(Hres[[3]],Mres[[3]],Zres[[3]])
mouse_zebrafish$MZ_index = paste(mouse_zebrafish$mouse,mouse_zebrafish$zebrafish,sep="__")
mouse_zebrafish = mouse_zebrafish[!duplicated(mouse_zebrafish$MZ_index),]


dim(Hres[[3]])
dim(Mres[[3]])
dim(Zres[[3]])
dim(mouse_zebrafish)

#######
#######

triple = rbind(Hres[[4]],Mres[[4]],Zres[[4]])
triple$HMZ_index = paste(triple$human,triple$mouse,triple$zebrafish,sep="__")
triple = triple[!duplicated(triple$HMZ_index),]

#######

dim(Hres[[4]])
dim(Mres[[4]])
dim(Zres[[4]])
dim(triple)

######
######
HMZ_ortholog_combined = list(HM=human_mouse,HZ=human_zebrafish,MZ=mouse_zebrafish,HMZ=triple)

#####
head(HMZ_ortholog_combined[[1]])
head(HMZ_ortholog_combined[[2]])
head(HMZ_ortholog_combined[[3]])

setwd("/zp1/data/share/Human_aging_new")
saveRDS(HMZ_ortholog_combined,file="HMZ_ortholog_combined_2025")