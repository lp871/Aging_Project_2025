######## load the datasets: #######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R


########
######## read the cell counts excel file #########
########
setwd("/zp1/data/plyu3/Human_aging")

library(readxl)

library(dplyr)
library(tidyr)
library(ggplot2)

# 文件路径
file_path <- "Supplemental Table 1- Cell Counts Across Aging in Zebrafish, Mouse, and Human (1).xlsx"

# 只读第 2-4 个 sheet
zebrafish <- read_excel(file_path, sheet = "Zebrafish")
mouse     <- read_excel(file_path, sheet = "Mouse")
human     <- read_excel(file_path, sheet = "Human")

########
######## plot the barplot for each species #######
########

zebrafish_long <- zebrafish %>%
  rename(CellType = 1) %>%              # 第一列改名
  pivot_longer(
    cols = -CellType, 
    names_to = "Time", 
    values_to = "Count"
  )

zebrafish_long$CellType = factor(zebrafish_long$CellType,levels=c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia"))
zebrafish_long$Time = factor(zebrafish_long$Time,levels=c(1,3,6,12,18,22,24,30,36,48))

########
########

library(ggplot2)

cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

ggplot(zebrafish_long,aes(x=Time,y=Count,fill=CellType)) + geom_bar(stat="identity",width=0.75,alpha=0.5) + theme_classic() + scale_fill_manual(values=cols) +
ylab("Number of cells") + xlab("Age") + scale_y_continuous(expand=c(0,0),limits=c(0,25000))

ggsave("Fish_aging_RNA_cell_bar.png",height=3,width=6)

#######
####### Next for the Mouse #######
#######


mouse_long <- mouse %>%
  rename(CellType = 1) %>%              # 第一列改名
  pivot_longer(
    cols = -CellType, 
    names_to = "Time", 
    values_to = "Count"
  )


########


mouse_long$CellType = factor(mouse_long$CellType,levels=c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia"))
mouse_long$Time = factor(mouse_long$Time,levels=c(5,12,17,32,49,68,91,108,120))

library(ggplot2)

cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

ggplot(mouse_long,aes(x=Time,y=Count,fill=CellType)) + geom_bar(stat="identity",width=0.75,alpha=0.5) + theme_classic() + scale_fill_manual(values=cols) +
ylab("Number of cells") + xlab("Age") + scale_y_continuous(expand=c(0,0),limits=c(0,25000))

ggsave("Mouse_aging_RNA_cell_bar.png",height=3,width=6)


########
########



Human_long <- human %>%
  rename(CellType = 1) %>%              # 第一列改名
  pivot_longer(
    cols = -CellType, 
    names_to = "Time", 
    values_to = "Count"
  )


########


Human_long$CellType = factor(Human_long$CellType,levels=c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia"))
Human_long$Time = factor(Human_long$Time,levels=c(10,25,37.5,50,60,70,80,90))


library(ggplot2)

cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

ggplot(Human_long,aes(x=Time,y=Count,fill=CellType)) + geom_bar(stat="identity",width=0.75,alpha=0.5) + theme_classic() + scale_fill_manual(values=cols) +
ylab("Number of cells") + xlab("Age") + scale_y_continuous(expand=c(0,0),limits=c(0,1000000))

ggsave("Human_aging_RNA_cell_bar.png",height=3,width=6)



