#######
####### This function is used to plot the GO enrichment for figure2 and figure3 ##########
#######

#######
####### figure3s selected GO terms for Seth Excel files ######
#######



ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R


library(ArchR)

/zp1/data/plyu3/GO_PLot

setwd("/zp1/data/plyu3/GO_PLot")


#####
##### 我们首先读入三个物种的 GO Plot #########
#####


##### for Mouse ######--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(openxlsx)

file   <- "TableS3_Mouse_Aging_DEGs_GOKEGG_May5_2025_Plot.xlsx"
sheets <- getSheetNames(file)

# 2. 用 lapply 按 sheet 读入
Mouse_GO <- lapply(sheets, function(sh) {
  read.xlsx(file, sheet = sh)
})

# 3. 给 list 元素命名
names(Mouse_GO) <- sheets
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##### for Zebrafish ######--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(openxlsx)

file   <- "TableS3_Zebrafish_Aging_DEGs_GOKEGG_May5_2025_Plot.xlsx"
sheets <- getSheetNames(file)

# 2. 用 lapply 按 sheet 读入
Zebrafish_GO <- lapply(sheets, function(sh) {
  read.xlsx(file, sheet = sh)
})

# 3. 给 list 元素命名
names(Zebrafish_GO) <- sheets
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##### for Human ######--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(openxlsx)

file   <- "TableS3_Human_Aging_DEGs_GOKEGG_May5_2025.xlsx"
sheets <- getSheetNames(file)

# 2. 用 lapply 按 sheet 读入
Human_GO <- lapply(sheets, function(sh) {
  read.xlsx(file, sheet = sh)
})

# 3. 给 list 元素命名
names(Human_GO) <- sheets


######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######---------
######---------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Next for Mouse Young samples #######
######

Mouse_GO_Y = Mouse_GO[grep("Young",names(Mouse_GO))]
Mouse_GO_Y = do.call('rbind',Mouse_GO_Y)
Mouse_GO_Y$celltype = sapply(strsplit(rownames(Mouse_GO_Y),split="_"), function(x) x[[1]])

Need_GO_terms = c("GO:0002181","GO:0008380","GO:0061621","GO:0043161","GO:0099111","GO:0044782","GO:0006119")

head(Mouse_GO_Y[grep("cytoplasmic translation",Mouse_GO_Y$Description),c(1,2,3,4)])
head(Mouse_GO_Y[grep("oxphos",Mouse_GO_Y$Description),c(1,2,3,4)])
head(Mouse_GO_Y[grep("oxidative phosphorylation",Mouse_GO_Y$Description),c(1,2,3,4)],n=50)
head(Mouse_GO_Y[grep("GO:0044782",Mouse_GO_Y$ID),c(1,2,3,4)],n=50)

######
######
Mouse_GO_Y_Cl = Mouse_GO_Y[which(Mouse_GO_Y$ID %in% Need_GO_terms == T),]
Mouse_GO_Y_Cl$logP = -log10(Mouse_GO_Y_Cl$pvalue)
######
######

library(ggplot2)

######
library(stringr)
Mouse_GO_Y_Cl$logP[which(Mouse_GO_Y_Cl$logP > 5)] = 5
Mouse_GO_Y_Cl$Description <- str_wrap(Mouse_GO_Y_Cl$Description, width = 30)



library(ggplot2)
ggplot(Mouse_GO_Y_Cl,aes(x=celltype,y=Description)) + geom_point(aes(size= Count,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='blue', guide = guide_colorbar(reverse = FALSE, order = 1),limits=c(0,5)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("Mouse_GO_Y.png",height=4,width=6)


##### Next for the Old samples in Mouse ######
#####


Mouse_GO_O = Mouse_GO[grep("Old",names(Mouse_GO))]
Mouse_GO_O = do.call('rbind',Mouse_GO_O)
Mouse_GO_O$celltype = sapply(strsplit(rownames(Mouse_GO_O),split="_"), function(x) x[[1]])

Need_GO_terms = c("GO:0050727","GO:0000910","GO:0048041","GO:0032868","GO:0010506","GO:0006302","GO:0042063","GO:0043114","mmu04068","GO:0016055")

head(Mouse_GO_O[grep("Wnt",Mouse_GO_O$Description),c(1,2,3,4)])
head(Mouse_GO_O[grep("double-strand break repair",Mouse_GO_O$Description),c(1,2,3,4)],n=40)
head(Mouse_GO_O[grep("regulation of autophagy",Mouse_GO_O$Description),c(1,2,3,4)],n=50)
head(Mouse_GO_O[grep("GO:0043114",Mouse_GO_O$ID),c(1,2,3,4)],n=50)

######
######
Mouse_GO_O_Cl = Mouse_GO_O[which(Mouse_GO_O$ID %in% Need_GO_terms == T),]
Mouse_GO_O_Cl$logP = -log10(Mouse_GO_O_Cl$pvalue)
######
######

library(ggplot2)

######
library(stringr)
Mouse_GO_O_Cl$logP[which(Mouse_GO_O_Cl$logP > 5)] = 5
Mouse_GO_O_Cl$Description <- str_wrap(Mouse_GO_O_Cl$Description, width = 30)


library(ggplot2)
ggplot(Mouse_GO_O_Cl,aes(x=celltype,y=Description)) + geom_point(aes(size= Count,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='red', guide = guide_colorbar(reverse = FALSE, order = 1),limits=c(0,5)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("Mouse_GO_O.png",height=4,width=6)


######
###### ----------- Next Zebrafish Young genes ----------------------------------------------------------------------------------------------------------------------------------------------------------------
######

Zebrafish_GO_Y = Zebrafish_GO[grep("Young",names(Zebrafish_GO))]
Zebrafish_GO_Y = do.call('rbind',Zebrafish_GO_Y)
Zebrafish_GO_Y$celltype = sapply(strsplit(rownames(Zebrafish_GO_Y),split="_"), function(x) x[[1]])

Need_GO_terms = c("dre00010","dre00020","GO:0001666","GO:0043484","GO:0051493")

head(Zebrafish_GO_Y[grep("cytoskeleton",Zebrafish_GO_Y$Description),c(1,2,3,4)])
head(Zebrafish_GO_Y[grep("TCA ",Zebrafish_GO_Y$Description),c(1,2,3,4)])
head(Zebrafish_GO_Y[grep("oxidative phosphorylation",Zebrafish_GO_Y$Description),c(1,2,3,4)],n=50)
head(Zebrafish_GO_Y[grep("GO:0051493",Zebrafish_GO_Y$ID),c(1,2,3,4)],n=50)

######
######
Zebrafish_GO_Y_Cl = Zebrafish_GO_Y[which(Zebrafish_GO_Y$ID %in% Need_GO_terms == T),]
Zebrafish_GO_Y_Cl$logP = -log10(Zebrafish_GO_Y_Cl$pvalue)
######
######

library(ggplot2)

######
library(stringr)
Zebrafish_GO_Y_Cl$logP[which(Zebrafish_GO_Y_Cl$logP > 5)] = 5
Zebrafish_GO_Y_Cl$Description <- str_wrap(Zebrafish_GO_Y_Cl$Description, width = 30)



library(ggplot2)
ggplot(Zebrafish_GO_Y_Cl,aes(x=celltype,y=Description)) + geom_point(aes(size= Count,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='blue', guide = guide_colorbar(reverse = FALSE, order = 1),limits=c(0,5)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("Zebrafish_GO_Y.png",height=4,width=6)

######## ----------- Next Zebrafish Old genes -----------------------------------------------------------------


Zebrafish_GO_O = Zebrafish_GO[grep("Old",names(Zebrafish_GO))]
Zebrafish_GO_O = do.call('rbind',Zebrafish_GO_O)
Zebrafish_GO_O$celltype = sapply(strsplit(rownames(Zebrafish_GO_O),split="_"), function(x) x[[1]])

Need_GO_terms = c("GO:0098609","GO:0061462","dre04520","GO:0006914","GO:0044782","GO:0061572")

head(Zebrafish_GO_O[grep("junction",Zebrafish_GO_O$Description),c(1,2,3,4)])
head(Zebrafish_GO_O[grep("filament",Zebrafish_GO_O$Description),c(1,2,3,4)],n=50)
head(Zebrafish_GO_O[grep("oxidative phosphorylation",Zebrafish_GO_O$Description),c(1,2,3,4)],n=50)
head(Zebrafish_GO_O[grep("GO:0061572",Zebrafish_GO_O$ID),c(1,2,3,4)],n=50)

######
######
Zebrafish_GO_O_Cl = Zebrafish_GO_O[which(Zebrafish_GO_O$ID %in% Need_GO_terms == T),]
Zebrafish_GO_O_Cl$logP = -log10(Zebrafish_GO_O_Cl$pvalue)
######
######

library(ggplot2)

######
library(stringr)
Zebrafish_GO_O_Cl$logP[which(Zebrafish_GO_O_Cl$logP > 5)] = 5
Zebrafish_GO_O_Cl$Description <- str_wrap(Zebrafish_GO_O_Cl$Description, width = 30)



library(ggplot2)
ggplot(Zebrafish_GO_O_Cl,aes(x=celltype,y=Description)) + geom_point(aes(size= Count,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='red', guide = guide_colorbar(reverse = FALSE, order = 1),limits=c(0,5)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("Zebrafish_GO_O.png",height=4,width=6)






###### ----------- -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### Next for the human Young samples #####
###### ----------- -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



Human_GO_Y = Human_GO[grep("Young",names(Human_GO))]
Human_GO_Y = do.call('rbind',Human_GO_Y)
Human_GO_Y$celltype = sapply(strsplit(rownames(Human_GO_Y),split="_"), function(x) x[[1]])

Need_GO_terms = c("GO:0006302","GO:0042254","GO:0061640","GO:0001774","GO:0050808","GO:0000723","GO:0007416")

head(Human_GO_Y[grep("telomere",Human_GO_Y$Description),c(1,2,3,4)])
head(Human_GO_Y[grep("TCA ",Human_GO_Y$Description),c(1,2,3,4)])
head(Human_GO_Y[grep("synapse",Human_GO_Y$Description),c(1,2,3,4)],n=50)
head(Human_GO_Y[grep("GO:0050808",Human_GO_Y$ID),c(1,2,3,4)],n=50)

######
######
Human_GO_Y_Cl = Human_GO_Y[which(Human_GO_Y$ID %in% Need_GO_terms == T),]
Human_GO_Y_Cl$logP = -log10(Human_GO_Y_Cl$pvalue)
######
######

library(ggplot2)

######
library(stringr)
Human_GO_Y_Cl$logP[which(Human_GO_Y_Cl$logP > 5)] = 5
Human_GO_Y_Cl$Description <- str_wrap(Human_GO_Y_Cl$Description, width = 30)

library(ggplot2)
ggplot(Human_GO_Y_Cl,aes(x=celltype,y=Description)) + geom_point(aes(size= Count,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='blue', guide = guide_colorbar(reverse = FALSE, order = 1),limits=c(0,5)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("Human_GO_Y.png",height=4,width=6)



###### ----------- -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### Next for the human Old samples #####
###### ----------- -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



Human_GO_O = Human_GO[grep("Old",names(Human_GO))]
Human_GO_O = do.call('rbind',Human_GO_O)
Human_GO_O$celltype = sapply(strsplit(rownames(Human_GO_O),split="_"), function(x) x[[1]])

Need_GO_terms = c("GO:0099003","GO:0090382","GO:0016236","GO:0007040","GO:0060271")

head(Human_GO_O[grep("vesicle-mediated ",Human_GO_O$Description),c(1,2,3,4)])
head(Human_GO_O[grep("lysosome",Human_GO_O$Description),c(1,2,3,4)])
head(Human_GO_O[grep("cilium",Human_GO_O$Description),c(1,2,3,4)],n=50)
head(Human_GO_O[grep("GO:0150146",Human_GO_O$ID),c(1,2,3,4)],n=50)

######
######
Human_GO_O_Cl = Human_GO_O[which(Human_GO_O$ID %in% Need_GO_terms == T),]
Human_GO_O_Cl$logP = -log10(Human_GO_O_Cl$pvalue)
######
######

library(ggplot2)

######
library(stringr)
Human_GO_O_Cl$logP[which(Human_GO_O_Cl$logP > 5)] = 5
Human_GO_O_Cl$Description <- str_wrap(Human_GO_O_Cl$Description, width = 30)

library(ggplot2)
ggplot(Human_GO_O_Cl,aes(x=celltype,y=Description)) + geom_point(aes(size= Count,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='red', guide = guide_colorbar(reverse = FALSE, order = 1),limits=c(0,5)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("Human_GO_O.png",height=4,width=6)












