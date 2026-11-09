#######
#######
#######


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R


######## load the mouse GRNs results #######
setwd("/zp1/data/plyu3/Aging_Clocks_Final")
load("Mouse_Corr_res_Sep23")
head(Mouse_Corr_res)

######## load the zebrafish GRNs results #######
setwd("/zp1/data/plyu3/Aging_Clocks_Final")
load("Fish_Corr_res_Sep23")
head(Fish_Corr_res)


######## load the human GRNs results #######
setwd("/zp1/data/plyu3/Aging_Clocks_Final")
load("Human_Corr_res_Sep23")
head(Human_Corr_res)


######## load the motifs ############---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

######## Mouse: #######

setwd("/zp1/data/plyu3/Aging_Clocks_Final/GRNs")
Mouse_Trans <- readRDS("Mouse_Trans")
Mouse_CisBP <- readRDS("Mouse_CisBP")


######## Zebrafish: #########

setwd("/zp1/data/plyu3/Aging_Clocks_Final/GRNs")
Fish_Trans <- readRDS("Fish_Trans")
Fish_CisBP <- readRDS("Fish_CisBP")

######## Human: #########

setwd("/zp1/data/plyu3/Aging_Clocks_Final/GRNs")
Human_Trans <- readRDS("Human_Trans")
Human_CisBP <- readRDS("Human_CisBP")



########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######--------- Next get bw for the mouse for each cell type ####################
########------

#########--------------------------------- for Mouse fragments !!!! --------------------------------------------------------------------------------------------------------------------------------------------
conda activate ArchR2 
R
library(ArchR)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
Mouse_ATAC_project <- readRDS("Project_all_clcl_Jul2")

#####
head(Mouse_ATAC_project@cellColData)
table(Mouse_ATAC_project$celltype)

#####
Celltypes = c("AC","BC","Cone","HC","MG","RGC","Rod")

Mouse_Fragments_Need_List = list()

for(i in 1:length(Celltypes)){
    print(i)
    ####
    cell_id = rownames(Mouse_ATAC_project@cellColData)[which(Mouse_ATAC_project$celltype == Celltypes[i])]
    ####
    Fragments_Need = getFragmentsFromProject(
        ArchRProj = Mouse_ATAC_project,
        cellNames = cell_id,
        verbose = FALSE,
        logFile = createLogFile("getFragmentsFromProject")
    )
    #####
    Mouse_Fragments_Need_List <- c(Mouse_Fragments_Need_List,list(Fragments_Need))
}

names(Mouse_Fragments_Need_List) <- Celltypes
saveRDS(Mouse_Fragments_Need_List,file="Mouse_Fragments_Need_List_2025")


source("/zp1/data/plyu3/All_Functions_2025/R/source_all.R")
source_all_r("/zp1/data/plyu3/All_Functions_2025/R/")

Celltypes = c("AC","BC","Cone","HC","MG","RGC","Rod")

Convert_frag_list_to_bedpe(Mouse_Fragments_Need_List$AC,TAG="Mouse_AC")
Convert_frag_list_to_bedpe(Mouse_Fragments_Need_List$BC,TAG="Mouse_BC")
Convert_frag_list_to_bedpe(Mouse_Fragments_Need_List$HC,TAG="Mouse_HC")
Convert_frag_list_to_bedpe(Mouse_Fragments_Need_List$RGC,TAG="Mouse_RGC")
Convert_frag_list_to_bedpe(Mouse_Fragments_Need_List$MG,TAG="Mouse_MG")
Convert_frag_list_to_bedpe(Mouse_Fragments_Need_List$Rod,TAG="Mouse_Rod")
Convert_frag_list_to_bedpe(Mouse_Fragments_Need_List$Cone,TAG="Mouse_Cone")



#########--------------------------------- for Zebrafish fragments !!!! --------------------------------------------------------------------------------------------------------------------------------------------
#########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2 
R
library(ArchR)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
##### save(Project_all_clclclcl,file="Zebrafish_ATAC_FinalFinal_2024_Jun26")

load("Zebrafish_ATAC_FinalFinal_2024_Jun26")
Zebrafish_ATAC_project = Project_all_clclclcl
table(Zebrafish_ATAC_project$celltype)
Celltypes = c("AC","BC","Cone","HC","MG","RGC","Rod")
Zebrafish_Fragments_Need_List = list()

for(i in 1:length(Celltypes)){
    print(i)
    ####
    cell_id = rownames(Zebrafish_ATAC_project@cellColData)[which(Zebrafish_ATAC_project$celltype == Celltypes[i])]
    ####
    Fragments_Need = getFragmentsFromProject(
        ArchRProj = Zebrafish_ATAC_project,
        cellNames = cell_id,
        verbose = FALSE,
        logFile = createLogFile("getFragmentsFromProject")
    )
    #####
    Zebrafish_Fragments_Need_List <- c(Zebrafish_Fragments_Need_List,list(Fragments_Need))
}

names(Zebrafish_Fragments_Need_List) <- Celltypes
saveRDS(Zebrafish_Fragments_Need_List,file="Zebrafish_Fragments_Need_List_2025")


Celltypes = c("AC","BC","Cone","HC","MG","RGC","Rod")

Convert_frag_list_to_bedpe(Zebrafish_Fragments_Need_List$AC,TAG="Zebrafish_AC")
Convert_frag_list_to_bedpe(Zebrafish_Fragments_Need_List$BC,TAG="Zebrafish_BC")
Convert_frag_list_to_bedpe(Zebrafish_Fragments_Need_List$HC,TAG="Zebrafish_HC")
Convert_frag_list_to_bedpe(Zebrafish_Fragments_Need_List$RGC,TAG="Zebrafish_RGC")
Convert_frag_list_to_bedpe(Zebrafish_Fragments_Need_List$MG,TAG="Zebrafish_MG")
Convert_frag_list_to_bedpe(Zebrafish_Fragments_Need_List$Rod,TAG="Zebrafish_Rod")
Convert_frag_list_to_bedpe(Zebrafish_Fragments_Need_List$Cone,TAG="Zebrafish_Cone")




#########--------------------------------- for Human fragments !!!! --------------------------------------------------------------------------------------------------------------------------------------------
#########

conda activate ArchR2 
R
library(ArchR)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
Human_ATAC_Project <- readRDS(file="Human_ATAC_Final_Jun2024")

#########


#####
head(Human_ATAC_Project@cellColData)
table(Human_ATAC_Project$celltype)

#####
Celltypes = c("AC","BC","Cone","HC","MG","RGC","Rod")

Human_Fragments_Need_List = list()

for(i in 1:length(Celltypes)){
    print(i)
    ####
    cell_id = rownames(Human_ATAC_Project@cellColData)[which(Human_ATAC_Project$celltype == Celltypes[i])]
    ####
    Fragments_Need = getFragmentsFromProject(
        ArchRProj = Human_ATAC_Project,
        cellNames = cell_id,
        verbose = FALSE,
        logFile = createLogFile("getFragmentsFromProject")
    )
    #####
    Human_Fragments_Need_List <- c(Human_Fragments_Need_List,list(Fragments_Need))
}

names(Human_Fragments_Need_List) <- Celltypes
saveRDS(Human_Fragments_Need_List,file="Human_Fragments_Need_List_2025")


source("/zp1/data/plyu3/All_Functions_2025/R/source_all.R")
source_all_r("/zp1/data/plyu3/All_Functions_2025/R/")

Celltypes = c("AC","BC","Cone","HC","MG","RGC","Rod")

Convert_frag_list_to_bedpe(Human_Fragments_Need_List$AC,TAG="Human_AC")
Convert_frag_list_to_bedpe(Human_Fragments_Need_List$BC,TAG="Human_BC")
Convert_frag_list_to_bedpe(Human_Fragments_Need_List$HC,TAG="Human_HC")
Convert_frag_list_to_bedpe(Human_Fragments_Need_List$RGC,TAG="Human_RGC")
Convert_frag_list_to_bedpe(Human_Fragments_Need_List$MG,TAG="Human_MG")
Convert_frag_list_to_bedpe(Human_Fragments_Need_List$Rod,TAG="Human_Rod")
Convert_frag_list_to_bedpe(Human_Fragments_Need_List$Cone,TAG="Human_Cone")


########
########
######## Output the genome size and peaks ！！ ###### 
########
########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")

genome_size = Mouse_ATAC_project@genomeAnnotation$chromSizes
genome_size_tab = data.frame(chr=seqnames(genome_size),length=end(genome_size))
write.table(genome_size_tab,file="mm10.sizes",quote=F,row.names=F,col.names=F,sep="\t")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")

genome_size = Zebrafish_ATAC_project@genomeAnnotation$chromSizes
genome_size_tab = data.frame(chr=seqnames(genome_size),length=end(genome_size))
write.table(genome_size_tab,file="danRer11.sizes",quote=F,row.names=F,col.names=F,sep="\t")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")

genome_size = Human_ATAC_Project@genomeAnnotation$chromSizes
genome_size_tab = data.frame(chr=seqnames(genome_size),length=end(genome_size))
write.table(genome_size_tab,file="hg38.sizes",quote=F,row.names=F,col.names=F,sep="\t")

#######
####### next find the peaks ！！！ #######
#######

#######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
Fish_PeakMat = readRDS("Fish_aging_PeakMat_Aug28")
Fish_PeakMat_Peaks = data.frame(GRanges(rownames(Fish_PeakMat)))
write.table(Fish_PeakMat_Peaks[,c(1:3)],file="Fish_PeakMat_Peaks.bed",sep="\t",quote=F,col.names=F,row.names=F)

######## on HPC #########
cp /projects/hmz-aging/Aging_Clocks_Final/All_Peaks_GRNs_Mouse_Sep10 /home/plyu3/TRANSFER

######## not working !!##
U[9C20&&
########


setwd("/projects/hmz-aging/Aging_Clocks_Final")
All_Peaks_GRNs_Mouse <- readRDS("All_Peaks_GRNs_Mouse_Sep10")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
All_Peaks_GRNs_Mouse <- readRDS("All_Peaks_GRNs_Mouse_Sep10")
All_Peaks_GRNs_Mouse <- data.frame(All_Peaks_GRNs_Mouse)
write.table(All_Peaks_GRNs_Mouse[,c(1:3)],file="Mouse_PeakMat_Peaks.bed",sep="\t",quote=F,col.names=F,row.names=F)

#######
####### Next Human Peaks！！！ #######
#######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
All_Peaks_GRNs_Human <- readRDS("All_Peaks_GRNs_Human_Sep10")
All_Peaks_GRNs_Human <- data.frame(All_Peaks_GRNs_Human)
write.table(All_Peaks_GRNs_Human[,c(1:3)],file="Human_PeakMat_Peaks.bed",sep="\t",quote=F,col.names=F,row.names=F)



####### globus_env #######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R



conda create -n globus_env python=3.9 -y
conda activate globus_env
mamba install -c conda-forge globus-cli

globus gcp create mapped "omb2_plyu3"

globus endpoint create --personal "omb2_plyu3"

./globusconnectpersonal -setup 8c716978-2e06-46a0-93b1-6eb3cf9f2b90

conda activate globus_env
/home/plyu3/globusconnectpersonal-3.2.6/globusconnectpersonal -start

/zp1/data/ShareData

globus ls fb5c6d1f-2c32-11f0-989b-02718ce60f75:/zp1/data/ShareData

globus endpoint search --filter-display-name

#######



#######
####### Next run bedtools #######
#######


####### for Fish !!! ########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR

cd /zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out

for cell in Cone Rod AC HC BC RGC MG; do
  bedtools bedpetobam \
    -i Zebrafish_${cell}_fragments_cl_bamGR_pe.bed \
    -g danRer11.sizes \
    > Zebrafish_${cell}_fragments_cl_bamGR_pe.bam \
  && samtools sort \
    -o Zebrafish_${cell}_fragments_cl_bamGR_pe_s.bam \
    Zebrafish_${cell}_fragments_cl_bamGR_pe.bam \
  && samtools index \
    Zebrafish_${cell}_fragments_cl_bamGR_pe_s.bam
done

#########------- Next for Mouse ######

cd /zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow

bedtools bedpetobam -i Mouse_Cone_fragments_cl_bamGR_pe.bed -g mm10.sizes \
  > Mouse_Cone_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Mouse_Cone_fragments_cl_bamGR_pe_s.bam Mouse_Cone_fragments_cl_bamGR_pe.bam \
&& samtools index Mouse_Cone_fragments_cl_bamGR_pe_s.bam

bedtools bedpetobam -i Mouse_Rod_fragments_cl_bamGR_pe.bed -g mm10.sizes \
  > Mouse_Rod_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Mouse_Rod_fragments_cl_bamGR_pe_s.bam Mouse_Rod_fragments_cl_bamGR_pe.bam \
&& samtools index Mouse_Rod_fragments_cl_bamGR_pe_s.bam

bedtools bedpetobam -i Mouse_AC_fragments_cl_bamGR_pe.bed -g mm10.sizes \
  > Mouse_AC_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Mouse_AC_fragments_cl_bamGR_pe_s.bam Mouse_AC_fragments_cl_bamGR_pe.bam \
&& samtools index Mouse_AC_fragments_cl_bamGR_pe_s.bam

bedtools bedpetobam -i Mouse_HC_fragments_cl_bamGR_pe.bed -g mm10.sizes \
  > Mouse_HC_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Mouse_HC_fragments_cl_bamGR_pe_s.bam Mouse_HC_fragments_cl_bamGR_pe.bam \
&& samtools index Mouse_HC_fragments_cl_bamGR_pe_s.bam

bedtools bedpetobam -i Mouse_BC_fragments_cl_bamGR_pe.bed -g mm10.sizes \
  > Mouse_BC_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Mouse_BC_fragments_cl_bamGR_pe_s.bam Mouse_BC_fragments_cl_bamGR_pe.bam \
&& samtools index Mouse_BC_fragments_cl_bamGR_pe_s.bam

bedtools bedpetobam -i Mouse_RGC_fragments_cl_bamGR_pe.bed -g mm10.sizes \
  > Mouse_RGC_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Mouse_RGC_fragments_cl_bamGR_pe_s.bam Mouse_RGC_fragments_cl_bamGR_pe.bam \
&& samtools index Mouse_RGC_fragments_cl_bamGR_pe_s.bam

bedtools bedpetobam -i Mouse_MG_fragments_cl_bamGR_pe.bed -g mm10.sizes \
  > Mouse_MG_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Mouse_MG_fragments_cl_bamGR_pe_s.bam Mouse_MG_fragments_cl_bamGR_pe.bam \
&& samtools index Mouse_MG_fragments_cl_bamGR_pe_s.bam


#########------- Next for Human ######

cd /zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq

bedtools bedpetobam -i Human_Cone_fragments_cl_bamGR_pe.bed -g hg38.sizes \
  > Human_Cone_fragments_cl_bamGR_pe.bam \
&& samtools sort -o Human_Cone_fragments_cl_bamGR_pe_s.bam Human_Cone_fragments_cl_bamGR_pe.bam \
&& samtools index Human_Cone_fragments_cl_bamGR_pe_s.bam


for cell in Rod AC HC BC RGC MG; do
  bedtools bedpetobam \
    -i Human_${cell}_fragments_cl_bamGR_pe.bed \
    -g hg38.sizes \
    > Human_${cell}_fragments_cl_bamGR_pe.bam \
  && samtools sort \
    -o Human_${cell}_fragments_cl_bamGR_pe_s.bam \
    Human_${cell}_fragments_cl_bamGR_pe.bam \
  && samtools index \
    Human_${cell}_fragments_cl_bamGR_pe_s.bam
done

########
######## Next run the TOBIAS #######-----------
########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

########
conda activate ArchR
R

TOBIAS
mamba install bioconda::tobias

cd /zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq

####wget ftp://ftp.ensembl.org/pub/release-*/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
####gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
seqs <- getSeq(genome, names(genome)[1:25])

####
####
####


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")

writeXStringSet(
  seqs,
  filepath = "Homo_sapiens.hg38.primary_assembly.fa",
  format   = "fasta"
)

cd /zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq

# 循环执行 TOBIAS ATACorrect
for cell in Rod Cone AC BC HC RGC MG; do
  TOBIAS ATACorrect \
    --read_shift -1 -1 \
    --bam "Human_${cell}_fragments_cl_bamGR_pe_s.bam" \
    --genome Homo_sapiens.hg38.primary_assembly.fa \
    --peaks Human_PeakMat_Peaks.bed \
    --outdir "Human_${cell}_res" \
    --cores 35
done

for cell in Rod RGC MG; do
  TOBIAS ATACorrect \
    --read_shift -1 -1 \
    --bam "Human_${cell}_fragments_cl_bamGR_pe_s.bam" \
    --genome Homo_sapiens.hg38.primary_assembly.fa \
    --peaks Human_PeakMat_Peaks.bed \
    --outdir "Human_${cell}_res" \
    --cores 35
done


######
###### Next run Mouse ########
######

library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
seqs <- getSeq(genome, names(genome))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
writeXStringSet(
  seqs,
  filepath = "Mus_musculus.mm10.primary_assembly.fa",
  format   = "fasta"
)

####
#### GRanges to bam #######
####

cd /zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow


####
####

for cell in Rod Cone AC BC HC RGC MG; do
  TOBIAS ATACorrect --read_shift -1 -1 \
    --bam Mouse_${cell}_fragments_cl_bamGR_pe_s.bam \
    --genome Mus_musculus.mm10.primary_assembly.fa \
    --peaks Mouse_PeakMat_Peaks.bed \
    --outdir Mouse_${cell}_res \
    --cores 30
done


####
####
####




##########
##########
##########

BSgenome.Drerio.UCSC.danRer11
library(BSgenome.Drerio.UCSC.danRer11)
cd /zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out
genome <- BSgenome.Drerio.UCSC.danRer11
seqs <- getSeq(genome, names(genome)[1:26])

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
writeXStringSet(
  seqs,
  filepath = "Danio_rerio.danRer11.primary_assembly.fa",
  format   = "fasta"
)

###### run #######
cd /zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out

for cell in Rod Cone AC BC HC RGC MG; do
  TOBIAS ATACorrect --read_shift -1 -1 \
    --bam Zebrafish_${cell}_fragments_cl_bamGR_pe_s.bam \
    --genome Danio_rerio.danRer11.primary_assembly.fa \
    --peaks Fish_PeakMat_Peaks.bed \
    --outdir Zebrafish_${cell}_res \
    --cores 30
done




############### Test the raw signal：##################
############### And the TOBIAS output： ###############

###### load the Mouse HC fragments #######
######

conda activate ArchR2 
R
library(ArchR)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
Mouse_ATAC_project <- readRDS("Project_all_clcl_Jul2")

#####
head(Mouse_ATAC_project@cellColData)
table(Mouse_ATAC_project$celltype)

#####
Celltypes = c("HC")
cell_id = rownames(Mouse_ATAC_project@cellColData)[which(Mouse_ATAC_project$celltype == Celltypes)]
Fragments_Need = getFragmentsFromProject(
        ArchRProj = Mouse_ATAC_project,
        cellNames = cell_id,
        verbose = FALSE,
        logFile = createLogFile("getFragmentsFromProject")
)

all_fragments = GRanges()
for(i in 1:length(Fragments_Need)){
        print(i)
        all_fragments <- c(all_fragments,Fragments_Need[[i]])
}

###### convert to bw file #######
###### 没啥问题 ###################

library(ArchR)
chr = seqnames(all_fragments)
s = start(all_fragments)
e = end(all_fragments)
GR_plus = GRanges(seqnames=chr,IRanges(s,s),strand='*')
GR_minus = GRanges(seqnames=chr,IRanges(e,e),strand='*')
GR <- c(GR_plus,GR_minus)

library(BSgenome.Mmusculus.UCSC.mm10)
seqinfo(GR) <- seqinfo(Mmusculus)[seqlevels(GR)]

cov_mm10 <- coverage(GR)
seqinfo(cov_mm10) <- seqinfo(GR)

###### 
export(cov_mm10, "HC_mm10_test.bw", format="BigWig")

######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######----------------------------
######----------------------------


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

######----------######

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_HC_res/Mouse_HC_fragments_cl_bamGR_pe_s_uncorrected.bw"
library(rtracklayer)
temp_bw = import.bw(file)
summary(width(temp_bw))
temp_bw[which(temp_bw$score > 4),]
saveRDS(Mouse_Rod_footprint_res,file="Mouse_Rod_footprint_res_2024")

#####
##### Next calculate the footprint score !!!! #############
#####

##### for Mouse: #######---------------------------------------------------------------------------------------------------------
cp  Mouse_total_PMWs_Res_Sep23 /zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/")
Mouse_total_PMWs_Res = readRDS("Mouse_total_PMWs_Res_Sep23")
Mouse_total_PMWs_Res_unlist = unlist(Mouse_total_PMWs_Res)
saveRDS(Mouse_total_PMWs_Res_unlist,file="Mouse_total_PMWs_Res_unlist")




##### for Zebrafish: #######---------------------------------------------------------------------------------------------------------
cp Fish_total_PMWs_Res_Sep23 /zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
Fish_total_PMWs_Res = readRDS("Fish_total_PMWs_Res_Sep23")
Fish_total_PMWs_Res_unlist = unlist(Fish_total_PMWs_Res)
saveRDS(Fish_total_PMWs_Res_unlist,file="Fish_total_PMWs_Res_unlist")

####matchtopeaks_list = Fish_total_PMWs_Res
####matchtopeaks_list <- as.list(matchtopeaks_list)
#####combined_gr <- do.call(c, matchtopeaks_list)

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_Rod_res/Zebrafish_Rod_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_Rod_footprint_res = Add_score_to_the_motifs(file,matchtopeaks_list)
saveRDS(Zebrafish_Rod_footprint_res,file="Zebrafish_Rod_footprint_res_2025")


####
### Mouse_Motif_matchtopeaks <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Motif_matchtopeaks_Jan13.rds")
####

#####
#####
#####


Add_score_to_the_motifs <- function(file,matchtopeaks_list){
    ##########
    library(rtracklayer)
    temp_bw = import.bw(file)
    res_list <- Check_normalized_Signal(x=matchtopeaks_list,temp_bw)
    ########## length(res_list) #######
    return(res_list)
}

Check_normalized_Signal <- function(x,temp_bw){
    library(GenomicRanges)
	library(rtracklayer)
	temp_bw = temp_bw
	### print(summary(width(temp_bw)))
    ### for each Hints, see the footprint score ####
    ### 将 ATAC GRanges 转为覆盖度 (RleList) #######
    ### not work ####
    ### x = x[[1]]
    ### x = x[which(x$score > 2)]
    ###
    flank_size = width(x)*3
    c_start <- start(x)
    c_end   <- end(x)
    lf_start <- c_start - flank_size
    lf_end   <- c_start - 1
    rf_start <- c_end + 1
    rf_end   <- c_end + flank_size
    ####
    GR_mid = GRanges(seqnames=seqnames(x),IRanges(c_start,c_end))
    GR_left = GRanges(seqnames=seqnames(x),IRanges(lf_start,lf_end))
    GR_right = GRanges(seqnames=seqnames(x),IRanges(rf_start,rf_end))
    ####
    Get_score <- function(x,y){
        library(GenomicRanges)
        res = findOverlaps(x,y)
        res = data.frame(res)
        res$score = y$score[res$subjectHits]
        #### add sums ####
        res_sum = tapply(res$score,res$queryHits,sum)
        ####
        x$score = 0
        x$score[as.numeric(names(res_sum))] = as.numeric(res_sum)
        return(x)
    }
    #####
    GR_mid = Get_score(GR_mid,temp_bw)
    GR_left = Get_score(GR_left,temp_bw)
    GR_right = Get_score(GR_right,temp_bw)
    #####
    GR_mid$score = GR_mid$score / c(width(x))
    GR_left$score = GR_left$score / c(width(x)*3)
    GR_right$score = GR_right$score / c(width(x)*3)
    #####
    GR_mid$left = GR_left$score
    GR_mid$right = GR_right$score
    #####
    return(GR_mid)
}




##### for Human: #######---------------------------------------------------------------------------------------------------------
cp Human_total_PMWs_Res_Sep23  /zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
Human_total_PMWs_Res = readRDS("Human_total_PMWs_Res_Sep23")
Human_total_PMWs_Res_unlist = unlist(Human_total_PMWs_Res)
saveRDS(Human_total_PMWs_Res_unlist,file="Human_total_PMWs_Res_unlist")




########## 
########## load the Motif tables and other informations #########
########## 

#########-------- for Mouse ------####
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
Mouse_ATAC_project <- readRDS("Project_all_clcl_Jul2")

Mouse_Gene_Anno = Mouse_ATAC_project@geneAnnotation$genes
All_Peaks_GRNs_Mouse <- readRDS("All_Peaks_GRNs_Mouse_Sep10")
Mouse_Peak_to_Gene = Find_total_genes_peak_fun(Gtf=Mouse_Gene_Anno,allpeaks=All_Peaks_GRNs_Mouse)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
saveRDS(Mouse_Peak_to_Gene,file="Mouse_Peak_to_Gene")

#######------ for Zebrafish ########


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
##### save(Project_all_clclclcl,file="Zebrafish_ATAC_FinalFinal_2024_Jun26")
load("Zebrafish_ATAC_FinalFinal_2024_Jun26")
Zebrafish_ATAC_project = Project_all_clclclcl

Zebrafish_Gene_Anno = Zebrafish_ATAC_project@geneAnnotation$genes

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
Fish_PeakMat = readRDS("Fish_aging_PeakMat_Aug28")
All_Peaks_GRNs_Fish = GRanges(rownames(Fish_PeakMat))
Fish_Peak_to_Gene = Find_total_genes_peak_fun(Gtf=Zebrafish_Gene_Anno,allpeaks=All_Peaks_GRNs_Fish)

saveRDS(Fish_Peak_to_Gene,file="Fish_Peak_to_Gene")



#######------ for Human ########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
Human_ATAC_Project <- readRDS(file="Human_ATAC_Final_Jun2024")
Human_Gene_Anno = Human_ATAC_Project@geneAnnotation$genes
All_Peaks_GRNs_Human <- readRDS("All_Peaks_GRNs_Human_Sep10")


Human_Peak_to_Gene = Find_total_genes_peak_fun(Gtf=Human_Gene_Anno,allpeaks=All_Peaks_GRNs_Human)
saveRDS(Human_Peak_to_Gene,file="Human_Peak_to_Gene")


Find_total_genes_peak_fun <- function(Gtf,allpeaks){
    ############
    all_genes = names(table(Gtf$symbol))
    ############
    library(parallel)
    ############
    res = mclapply(all_genes,Find_one_genes_peak_fun,Gtf=Gtf,allpeaks=allpeaks,mc.cores= 40)
    ############
    ####
    return(res)
}


Find_one_genes_peak_fun = function(genename,Gtf,allpeaks){
    ############
    library(GenomicRanges)
    ############
    Target_G = Gtf[which(Gtf$symbol == genename)][1]
    Other_G = Gtf[which(Gtf$symbol != genename)]
    ############
    Target_G_tss_GR = promoters(Target_G, upstream = 1000, downstream = 1000)
    Target_G_ext_GR = promoters(Target_G, upstream = 500000, downstream = 500000)
    ############
    ############
    TSS_peaks = as.character(allpeaks)[which(countOverlaps(allpeaks,Target_G_tss_GR) > 0)]
    if(length(TSS_peaks) == 0){
        TSS_peaks_tab = data.frame(peaks="NA",class="NA")
    }else{
        TSS_peaks_tab = data.frame(peaks=TSS_peaks,class="TSS")
    }
    Body_peaks = as.character(allpeaks)[which(countOverlaps(allpeaks,Target_G) > 0)]
    if(length(Body_peaks) == 0){
        Body_peaks_tab = data.frame(peaks="NA",class="NA")
    }else{
        Body_peaks_tab = data.frame(peaks=Body_peaks,class="Body")
    }
    Inter_peaks = as.character(allpeaks)[which(countOverlaps(allpeaks,Target_G_ext_GR) > 0)]
    if(length(Inter_peaks) == 0){
        Inter_peaks_tab = data.frame(peaks="NA",class="NA")
    }else{
        Inter_peaks_tab = data.frame(peaks=Inter_peaks,class="Inter")
    }
    ############
    ############
    Gene_table = rbind(TSS_peaks_tab,Body_peaks_tab,Inter_peaks_tab)
    K = which(Gene_table$peaks != "NA")
    if(length(K) == 0){
        return(data.frame(peaks="NA",class="NA",gene=genename))
    }
    Gene_table = Gene_table[K,]
    ############
    Gene_table = Gene_table[!duplicated(Gene_table$peaks),]
    ############
    ############ Next we will filterout peaks with overlap with Other genes TSS region and TSS body #####
    Gene_table_inter_peaks = GRanges(Gene_table$peaks[which(Gene_table$class == "Inter")])
    ############ get other G tss #####
    Other_G_TSS <- promoters(Other_G, upstream = 1000, downstream = 1000)
    ############ filterout with other gene body ####
    ############
    body_index = which(countOverlaps(Gene_table_inter_peaks,Other_G) > 0)
    tss_index = which(countOverlaps(Gene_table_inter_peaks,Other_G_TSS) > 0)
    ############
    remove_peaks = Gene_table_inter_peaks[c(body_index,tss_index)]
    ############
    if(length(remove_peaks) > 0){
        k = which(Gene_table$class == "Inter" & Gene_table$peaks %in% as.character(remove_peaks) == T)
        Gene_table = Gene_table[-k,]
    }
    ############
    if(dim(Gene_table)[1] == 0){
        return(data.frame(peaks="NA",class="NA",gene=genename))
    }
    Gene_table$gene=genename
    return(Gene_table)
}


######
###### Next Get_TSSandEnhancer_From_PtoG #######
######

######-------- for Zebrafish #############

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
Fish_Peak_to_Gene = readRDS("Fish_Peak_to_Gene")
Fish_Peak_to_Gene = do.call("rbind",Fish_Peak_to_Gene)
Fish_PtoG = readRDS("/home/plyu3/data/Fish_PtoG_Ann_Aug29")
### Peak_to_Gene_merge = Fish_Peak_to_Gene ###
### PtoG = Fish_PtoG
#####
#####
#####
Fish_network = readRDS("/home/plyu3/data/Fish_Network_cl_Sep26")
Fish_motif = readRDS("/home/plyu3/data/Fish_short_motif_Sep25")

#####

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
Fish_total_PMWs_Res_unlist = readRDS("Fish_total_PMWs_Res_unlist")
### PMWs_Res_unlist = Fish_total_PMWs_Res_unlist
Fish_total_PMWs_Res_unlist_cl = Filter_PMWs_regions(Fish_total_PMWs_Res_unlist,Fish_network,Fish_motif,Fish_PtoG)
saveRDS(Fish_total_PMWs_Res_unlist_cl,file="Fish_total_PMWs_Res_unlist_cl")

Fish_total_PMWs_Res_unlist_cl = readRDS("Fish_total_PMWs_Res_unlist_cl")

Filter_PMWs_regions <- function(PMWs_Res_unlist,Fish_network,Fish_motif,Fish_PtoG){
    ############
    allGenes = unique(Fish_network$Target)
    allTFs = unique(Fish_network$TF)
    ############
    allMotifs = unique(Fish_motif$ID[which(Fish_motif$Gene %in% allTFs == T)])
    ############
    PMWs_Res_unlist_cl1 = PMWs_Res_unlist[which(names(PMWs_Res_unlist) %in% allMotifs == T)]
    ############
    allpeaks = Fish_PtoG$Peak[which(Fish_PtoG$Gene %in% allGenes == T)]
    allpeaks = GRanges(allpeaks)
    ############
    index = which(countOverlaps(PMWs_Res_unlist_cl1,allpeaks) > 0)
    PMWs_Res_unlist_cl2 = PMWs_Res_unlist_cl1[index]
    ############
    index2 = which(PMWs_Res_unlist_cl2$score > 4)
    PMWs_Res_unlist_cl3 = PMWs_Res_unlist_cl2[index2]
    ############
    return(PMWs_Res_unlist_cl3)
}

######
###### Let see the speed of footprint !!! #########
######

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_Rod_res/Zebrafish_Rod_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_Rod_footprint_res = Add_score_to_the_motifs(file,Fish_total_PMWs_Res_unlist_cl)
saveRDS(Zebrafish_Rod_footprint_res,file="Zebrafish_Rod_footprint_res_2025")

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_Cone_res/Zebrafish_Cone_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_Cone_footprint_res = Add_score_to_the_motifs(file,Fish_total_PMWs_Res_unlist_cl)
saveRDS(Zebrafish_Cone_footprint_res,file="Zebrafish_Cone_footprint_res_2025")

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_AC_res/Zebrafish_AC_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_AC_footprint_res = Add_score_to_the_motifs(file,Fish_total_PMWs_Res_unlist_cl)
saveRDS(Zebrafish_AC_footprint_res,file="Zebrafish_AC_footprint_res_2025")

# RGC
file_RGC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_RGC_res/Zebrafish_RGC_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_RGC_footprint_res <- Add_score_to_the_motifs(file_RGC, Fish_total_PMWs_Res_unlist_cl)
saveRDS(Zebrafish_RGC_footprint_res, file = "Zebrafish_RGC_footprint_res_2025.rds")

# HC
file_HC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_HC_res/Zebrafish_HC_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_HC_footprint_res <- Add_score_to_the_motifs(file_HC, Fish_total_PMWs_Res_unlist_cl)
saveRDS(Zebrafish_HC_footprint_res, file = "Zebrafish_HC_footprint_res_2025.rds")

# MG
file_MG <- "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_MG_res/Zebrafish_MG_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_MG_footprint_res <- Add_score_to_the_motifs(file_MG, Fish_total_PMWs_Res_unlist_cl)
saveRDS(Zebrafish_MG_footprint_res, file = "Zebrafish_MG_footprint_res_2025.rds")

# BC
file_BC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out/Zebrafish_BC_res/Zebrafish_BC_fragments_cl_bamGR_pe_s_corrected.bw"
Zebrafish_BC_footprint_res <- Add_score_to_the_motifs(file_BC, Fish_total_PMWs_Res_unlist_cl)
saveRDS(Zebrafish_BC_footprint_res, file = "Zebrafish_BC_footprint_res_2025.rds")


######
###### Let us the mouse footprint !!! #######
######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
Mouse_Peak_to_Gene = readRDS("Mouse_Peak_to_Gene")
Mouse_Peak_to_Gene = do.call("rbind",Mouse_Peak_to_Gene)
Mouse_PtoG = readRDS("/home/plyu3/data/Mouse_PtoG_Ann_Sep10")
Mouse_network = readRDS("/home/plyu3/data/Mouse_Network_cl_Sep26")
Mouse_motif = readRDS("/home/plyu3/data/Mouse_short_motif_Sep25")
Mouse_total_PMWs_Res_unlist = readRDS("Mouse_total_PMWs_Res_unlist")
Mouse_total_PMWs_Res_unlist = unlist(Mouse_total_PMWs_Res_unlist)

Mouse_total_PMWs_Res_unlist_cl = Filter_PMWs_regions(Mouse_total_PMWs_Res_unlist,Mouse_network,Mouse_motif,Mouse_PtoG)
saveRDS(Mouse_total_PMWs_Res_unlist_cl,file="Mouse_total_PMWs_Res_unlist_cl")

Mouse_total_PMWs_Res_unlist_cl <- readRDS("Mouse_total_PMWs_Res_unlist_cl")
#######

# RGC
file_RGC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_RGC_res/Mouse_RGC_fragments_cl_bamGR_pe_s_corrected.bw"
Mouse_RGC_footprint_res <- Add_score_to_the_motifs(file_RGC, Mouse_total_PMWs_Res_unlist_cl)
saveRDS(Mouse_RGC_footprint_res, file = "Mouse_RGC_footprint_res_2025.rds")

# AC
file_AC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_AC_res/Mouse_AC_fragments_cl_bamGR_pe_s_corrected.bw"
Mouse_AC_footprint_res <- Add_score_to_the_motifs(file_AC, Mouse_total_PMWs_Res_unlist_cl)
saveRDS(Mouse_AC_footprint_res, file = "Mouse_AC_footprint_res_2025.rds")

# HC
file_HC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_HC_res/Mouse_HC_fragments_cl_bamGR_pe_s_corrected.bw"
Mouse_HC_footprint_res <- Add_score_to_the_motifs(file_HC, Mouse_total_PMWs_Res_unlist_cl)
saveRDS(Mouse_HC_footprint_res, file = "Mouse_HC_footprint_res_2025.rds")

# Rod
file_Rod <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_Rod_res/Mouse_Rod_fragments_cl_bamGR_pe_s_corrected.bw"
Mouse_Rod_footprint_res <- Add_score_to_the_motifs(file_Rod, Mouse_total_PMWs_Res_unlist_cl)
saveRDS(Mouse_Rod_footprint_res, file = "Mouse_Rod_footprint_res_2025.rds")

# Cone
file_Cone <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_Cone_res/Mouse_Cone_fragments_cl_bamGR_pe_s_corrected.bw"
Mouse_Cone_footprint_res <- Add_score_to_the_motifs(file_Cone, Mouse_total_PMWs_Res_unlist_cl)
saveRDS(Mouse_Cone_footprint_res, file = "Mouse_Cone_footprint_res_2025.rds")

# MG
file_MG <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_MG_res/Mouse_MG_fragments_cl_bamGR_pe_s_corrected.bw"
Mouse_MG_footprint_res <- Add_score_to_the_motifs(file_MG, Mouse_total_PMWs_Res_unlist_cl)
saveRDS(Mouse_MG_footprint_res, file = "Mouse_MG_footprint_res_2025.rds")

# BC
file_BC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow/Mouse_BC_res/Mouse_BC_fragments_cl_bamGR_pe_s_corrected.bw"
Mouse_BC_footprint_res <- Add_score_to_the_motifs(file_BC, Mouse_total_PMWs_Res_unlist_cl)
saveRDS(Mouse_BC_footprint_res, file = "Mouse_BC_footprint_res_2025.rds")


######
######
###### Let us see the Human footprint !!! #######
######
######


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
Human_Peak_to_Gene = readRDS("Human_Peak_to_Gene")
Human_Peak_to_Gene = do.call("rbind",Human_Peak_to_Gene)
Human_PtoG = readRDS("/home/plyu3/data/Human_PtoG_Ann_Sep10")
Human_network = readRDS("/home/plyu3/data/Human_Network_cl_Sep26")
Human_motif = readRDS("/home/plyu3/data/Human_short_motif_Sep25")
Human_total_PMWs_Res_unlist = readRDS("Human_total_PMWs_Res_unlist")

Human_total_PMWs_Res_unlist = unlist(Human_total_PMWs_Res_unlist)

Human_total_PMWs_Res_unlist_cl = Filter_PMWs_regions(Human_total_PMWs_Res_unlist,Human_network,Human_motif,Human_PtoG)
saveRDS(Human_total_PMWs_Res_unlist_cl,file="Human_total_PMWs_Res_unlist_cl")

length(Human_total_PMWs_Res_unlist)
length(Human_total_PMWs_Res_unlist_cl)

########
########
########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
Human_total_PMWs_Res_unlist_cl <- readRDS("Human_total_PMWs_Res_unlist_cl")

#######
#######
#######

# Rod
file = "/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq/Human_Rod_res/Human_Rod_fragments_cl_bamGR_pe_s_corrected.bw"
Human_Rod_footprint_res = Add_score_to_the_motifs(file,Human_total_PMWs_Res_unlist_cl)
saveRDS(Human_Rod_footprint_res,file="Human_Rod_footprint_res_2025")

# Cone
file_Cone <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq/Human_Cone_res/Human_Cone_fragments_cl_bamGR_pe_s_corrected.bw"
Human_Cone_footprint_res <- Add_score_to_the_motifs(file_Cone, Human_total_PMWs_Res_unlist_cl)
saveRDS(Human_Cone_footprint_res, file = "Human_Cone_footprint_res_2025.rds")

# AC
file_AC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq/Human_AC_res/Human_AC_fragments_cl_bamGR_pe_s_corrected.bw"
Human_AC_footprint_res <- Add_score_to_the_motifs(file_AC, Human_total_PMWs_Res_unlist_cl)
saveRDS(Human_AC_footprint_res, file = "Human_AC_footprint_res_2025.rds")

# BC
file_BC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq/Human_BC_res/Human_BC_fragments_cl_bamGR_pe_s_corrected.bw"
Human_BC_footprint_res <- Add_score_to_the_motifs(file_BC, Human_total_PMWs_Res_unlist_cl)
saveRDS(Human_BC_footprint_res, file = "Human_BC_footprint_res_2025.rds")

# HC
file_HC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq/Human_HC_res/Human_HC_fragments_cl_bamGR_pe_s_corrected.bw"
Human_HC_footprint_res <- Add_score_to_the_motifs(file_HC, Human_total_PMWs_Res_unlist_cl)
saveRDS(Human_HC_footprint_res, file = "Human_HC_footprint_res_2025.rds")

# RGC
file_RGC <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq/Human_RGC_res/Human_RGC_fragments_cl_bamGR_pe_s_corrected.bw"
Human_RGC_footprint_res <- Add_score_to_the_motifs(file_RGC, Human_total_PMWs_Res_unlist_cl)
saveRDS(Human_RGC_footprint_res, file = "Human_RGC_footprint_res_2025.rds")

# MG
file_MG <- "/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq/Human_MG_res/Human_MG_fragments_cl_bamGR_pe_s_corrected.bw"
Human_MG_footprint_res <- Add_score_to_the_motifs(file_MG, Human_total_PMWs_Res_unlist_cl)
saveRDS(Human_MG_footprint_res, file = "Human_MG_footprint_res_2025.rds")



###########
























#######----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




Get_TSSandEnhancer_From_PtoG <- function(Peak_to_Gene_merge,PtoG){
    ###############
    PtoG$index = paste(PtoG$Gene,PtoG$Peak)
    Peak_to_Gene_merge$index = paste(Peak_to_Gene_merge$gene,Peak_to_Gene_merge$peaks)
    ############### 这一步其实不需要 ###########
    Peak_to_Gene_merge_cl = Peak_to_Gene_merge
    ############### 不完全match也没事 ####
    k = which(PtoG$Peak %in% Peak_to_Gene_merge_cl$peaks == T)
    print(length(k))
    print(dim(PtoG))
    ###############
    m = match(Peak_to_Gene_merge_cl$index,PtoG$index)
    Peak_to_Gene_merge_cl$Correlation = as.numeric(PtoG$PtoG_cor[m])
    Peak_to_Gene_merge_cl$FDR = as.numeric(PtoG$PtoG_FDR[m])
    ###############
    k1_keep = which(Peak_to_Gene_merge_cl$class == "TSS")
    k2_keep = which(Peak_to_Gene_merge_cl$class == "Body" & Peak_to_Gene_merge_cl$Correlation > 0.25 & Peak_to_Gene_merge_cl$FDR < 0.01)
    k3_keep = which(Peak_to_Gene_merge_cl$class == "Inter" & Peak_to_Gene_merge_cl$Correlation > 0.25 & Peak_to_Gene_merge_cl$FDR < 0.01)
    ###############
    k_all = c(k1_keep,k2_keep,k3_keep)
    ###############
    Peak_to_Gene_merge_clcl = Peak_to_Gene_merge_cl[k_all,]
    ##### Peak_to_Gene_merge_clcl[which(Peak_to_Gene_merge_clcl$gene == "Thrb"),] ###
    ##### Peak_to_Gene_merge_clcl[] ####
    return(Peak_to_Gene_merge_clcl)
}

########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###########
########### 组装一下 GRNs ########----------------------------------------------------------------------------------------------------------------------------------------------------------------
###########
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)

####
#### 
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")

Fish_total_PMWs_Res_unlist_cl = readRDS("Fish_total_PMWs_Res_unlist_cl")
footprint = readRDS("Zebrafish_MG_footprint_res_2025.rds")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
Peak_to_Gene_add = readRDS("/home/plyu3/data/Fish_PtoG_Ann_Aug29")
motif_tf_table = readRDS("/home/plyu3/data/Fish_short_motif_Sep25")
total_GRNs_add = readRDS("/home/plyu3/data/Fish_Network_cl_Sep26")
Avg_Mat = readRDS("/home/plyu3/data/Zebrafish_celltype_exp_G_Sep26")
TAG = "MG"

#####
#####
#####
footprint = readRDS("Zebrafish_MG_footprint_res_2025.rds")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "MG"
Zebrafish_MG_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)


footprint = readRDS("Zebrafish_Rod_footprint_res_2025")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "Rod"
Zebrafish_Rod_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

footprint = readRDS("Zebrafish_Cone_footprint_res_2025")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "Cone"
Zebrafish_Cone_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)


footprint = readRDS("Zebrafish_AC_footprint_res_2025")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "AC"
Zebrafish_AC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)


footprint = readRDS("Zebrafish_BC_footprint_res_2025.rds")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "BC"
Zebrafish_BC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)


footprint = readRDS("Zebrafish_HC_footprint_res_2025.rds")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "HC"
Zebrafish_HC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)


footprint = readRDS("Zebrafish_MG_footprint_res_2025")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "MG"
Zebrafish_MG_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)


footprint = readRDS("Zebrafish_RGC_footprint_res_2025")
names(footprint) = names(Fish_total_PMWs_Res_unlist_cl)
TAG = "RGC"
Zebrafish_RGC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

#######
###########
#######

Zebrafish_GRNs_list = list(RGC=Zebrafish_RGC_GRNs,MG=Zebrafish_MG_GRNs,HC=Zebrafish_HC_GRNs,BC=Zebrafish_BC_GRNs,AC=Zebrafish_AC_GRNs,Cone=Zebrafish_Cone_GRNs,Rod=Zebrafish_Rod_GRNs)

####### see the dims ####

for(i in 1:length(Zebrafish_GRNs_list)){
    ##########
    tmp = Zebrafish_GRNs_list[[i]]
    tmptmp = tmp$TPG
    print(names(Zebrafish_GRNs_list)[i])
    print(dim(tmptmp))
}

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
saveRDS(Zebrafish_GRNs_list,file="Zebrafish_GRNs_list")

Zebrafish_GRNs_list <- readRDS("Zebrafish_GRNs_list")


############
#### for Human ############--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
############

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")

Human_total_PMWs_Res_unlist_cl = readRDS("Human_total_PMWs_Res_unlist_cl")
footprint = readRDS("Human_MG_footprint_res_2025.rds")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
Peak_to_Gene_add = readRDS("/home/plyu3/data/Human_PtoG_Ann_Sep10")
motif_tf_table = readRDS("/home/plyu3/data/Human_short_motif_Sep25")
total_GRNs_add = readRDS("/home/plyu3/data/Human_Network_cl_Sep26")
Avg_Mat = readRDS("/home/plyu3/data/Human_celltype_exp_G_Sep26")
TAG = "MG"

######
footprint = readRDS("Human_MG_footprint_res_2025.rds")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
TAG = "MG"
Human_MG_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

#######
footprint = readRDS("Human_HC_footprint_res_2025.rds")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
TAG = "HC"
Human_HC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

########
footprint = readRDS("Human_AC_footprint_res_2025.rds")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
TAG = "AC"
Human_AC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

########
footprint = readRDS("Human_BC_footprint_res_2025.rds")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
TAG = "BC"
Human_BC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

########
footprint = readRDS("Human_Cone_footprint_res_2025.rds")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
TAG = "Cone"
Human_Cone_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)


########
footprint = readRDS("Human_Rod_footprint_res_2025")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
TAG = "Rod"
Human_Rod_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)



########
footprint = readRDS("Human_RGC_footprint_res_2025.rds")
names(footprint) = names(Human_total_PMWs_Res_unlist_cl)
TAG = "RGC"
Human_RGC_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

#########
#########
#########

Human_GRNs_list = list(RGC=Human_RGC_GRNs,MG=Human_MG_GRNs,HC=Human_HC_GRNs,BC=Human_BC_GRNs,AC=Human_AC_GRNs,Cone=Human_Cone_GRNs,Rod=Human_Rod_GRNs)

####### see the dims ####

for(i in 1:length(Human_GRNs_list)){
    ##########
    tmp = Human_GRNs_list[[i]]
    tmptmp = tmp$TPG
    print(names(Human_GRNs_list)[i])
    print(dim(tmptmp))
}

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
saveRDS(Human_GRNs_list,file="Human_GRNs_list")

Human_GRNs_list = readRDS("Human_GRNs_list")

############
#### for Mouse ############--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
############

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")

Mouse_total_PMWs_Res_unlist_cl = readRDS("Mouse_total_PMWs_Res_unlist_cl")
footprint = readRDS("Mouse_MG_footprint_res_2025.rds")
names(footprint) = names(Mouse_total_PMWs_Res_unlist_cl)
Peak_to_Gene_add = readRDS("/home/plyu3/data/Mouse_PtoG_Ann_Sep10")
motif_tf_table = readRDS("/home/plyu3/data/Mouse_short_motif_Sep25")
total_GRNs_add = readRDS("/home/plyu3/data/Mouse_Network_cl_Sep26")
Avg_Mat = readRDS("/home/plyu3/data/Mouse_celltype_exp_G_Sep26")
TAG = "MG"

#######
#######
#######


######
footprint = readRDS("Mouse_MG_footprint_res_2025.rds")
names(footprint) = names(Mouse_total_PMWs_Res_unlist_cl)
TAG = "MG"
Mouse_MG_GRNs = Add_foot_print_to_peak(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat)

######
Mouse_Rod_footprint <- readRDS("Mouse_Rod_footprint_res_2025.rds")
names(Mouse_Rod_footprint) <- names(Mouse_total_PMWs_Res_unlist_cl)
Mouse_Rod_GRNs <- Add_foot_print_to_peak(Mouse_Rod_footprint, "Rod", Peak_to_Gene_add, motif_tf_table, total_GRNs_add, Avg_Mat)



Mouse_Cone_footprint <- readRDS("Mouse_Cone_footprint_res_2025.rds"); names(Mouse_Cone_footprint) <- names(Mouse_total_PMWs_Res_unlist_cl); Mouse_Cone_GRNs <- Add_foot_print_to_peak(Mouse_Cone_footprint, "Cone", Peak_to_Gene_add, motif_tf_table, total_GRNs_add, Avg_Mat)
Mouse_HC_footprint <- readRDS("Mouse_HC_footprint_res_2025.rds"); names(Mouse_HC_footprint) <- names(Mouse_total_PMWs_Res_unlist_cl); Mouse_HC_GRNs <- Add_foot_print_to_peak(Mouse_HC_footprint, "HC", Peak_to_Gene_add, motif_tf_table, total_GRNs_add, Avg_Mat)
Mouse_AC_footprint <- readRDS("Mouse_AC_footprint_res_2025.rds"); names(Mouse_AC_footprint) <- names(Mouse_total_PMWs_Res_unlist_cl); Mouse_AC_GRNs <- Add_foot_print_to_peak(Mouse_AC_footprint, "AC", Peak_to_Gene_add, motif_tf_table, total_GRNs_add, Avg_Mat)
Mouse_BC_footprint <- readRDS("Mouse_BC_footprint_res_2025.rds"); names(Mouse_BC_footprint) <- names(Mouse_total_PMWs_Res_unlist_cl); Mouse_BC_GRNs <- Add_foot_print_to_peak(Mouse_BC_footprint, "BC", Peak_to_Gene_add, motif_tf_table, total_GRNs_add, Avg_Mat)
Mouse_RGC_footprint <- readRDS("Mouse_RGC_footprint_res_2025.rds"); names(Mouse_RGC_footprint) <- names(Mouse_total_PMWs_Res_unlist_cl); Mouse_RGC_GRNs <- Add_foot_print_to_peak(Mouse_RGC_footprint, "RGC", Peak_to_Gene_add, motif_tf_table, total_GRNs_add, Avg_Mat)

Mouse_GRNs_list = list(RGC=Mouse_RGC_GRNs,MG=Mouse_MG_GRNs,HC=Mouse_HC_GRNs,BC=Mouse_BC_GRNs,AC=Mouse_AC_GRNs,Cone=Mouse_Cone_GRNs,Rod=Mouse_Rod_GRNs)

for(i in 1:length(Mouse_GRNs_list)){
    ##########
    tmp = Mouse_GRNs_list[[i]]
    tmptmp = tmp$TPG
    print(names(Mouse_GRNs_list)[i])
    print(dim(tmptmp))
}

#######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
saveRDS(Mouse_GRNs_list,file="Mouse_GRNs_list")

Mouse_GRNs_list = readRDS("Mouse_GRNs_list")
############




############
############
Add_foot_print_to_peak <- function(footprint,TAG,Peak_to_Gene_add,motif_tf_table,total_GRNs_add,Avg_Mat){
    ########## filter the footprint #####
    left_delta = footprint$left - footprint$score
    right_delta = footprint$right - footprint$score
    filter_foot_index = which(left_delta > 0.05 & right_delta > 0.05)
    ##########
    footprint_cl = footprint[filter_foot_index]
    ########## first we will filter motif by corr #####
    footprint_region = as.character(footprint_cl)
    footprint_name = names(footprint_cl)
    footprint_table = data.frame(Motif=names(footprint_cl),footprint=footprint_region)
    footprint_table = data.table(footprint_table)
    ########## merge footprint with peak ####
    findOverlaps_res = findOverlaps(GRanges(footprint_table$footprint),GRanges(Peak_to_Gene_add$Peak))
    ##########
    left_table = footprint_table[queryHits(findOverlaps_res),]
    right_table = Peak_to_Gene_add[subjectHits(findOverlaps_res),]
    merge_table = cbind(left_table,right_table)
    ########## Next merge tf motifs with footprint ####
    motif_tf_table_cl = motif_tf_table
    colnames(motif_tf_table_cl) = c("Motif","TF")
    motif_tf_table_cl = data.table(motif_tf_table_cl)
    ##########
    merge_table2 <- merge(
        x = merge_table,
        y = motif_tf_table_cl,
        by = c("Motif"),
        allow.cartesian = TRUE
    )
    ########## Next we will merge the gene-gene correaltion and importance ####
    total_GRNs_add = data.table(total_GRNs_add)
    ##########
    colnames(total_GRNs_add) = c("TF","Gene","GG_importance","GG_Corr")
    ##########
    merge_table3 <- merge(
        x = merge_table2,
        y = total_GRNs_add,
        by = c("TF","Gene")
    )
    ########### table(merge_table3$Class)
    ########### Next we will add the average expression of the TFs ################
    ####
    ####
    merge_table3$celltype = TAG
    Exp_G = Avg_Mat[[which(names(Avg_Mat) == TAG)]]
    #####
    #####
    ##### first we will filter PtoG_Cor ########
    ##### only positive ####
    #####
    k = which(merge_table3$PtoG_cor == "ND")
    merge_table3$PtoG_cor[k] = 0
    #### k1_filter = which(merge_table3$Class %in% c("Inter","Body") == T & as.numeric(merge_table3$PtoG_cor) < 0)
    #### merge_table3_f1 = merge_table3[-k1_filter,]
    #####
    k2_filter = which(merge_table3$GG_importance > 0.1 & abs(as.numeric(merge_table3$GG_Corr)) > 0.05)
    merge_table3_f2 = merge_table3[k2_filter,]
    ###### filter TF expression ######
    k3_filter = which(merge_table3_f2$TF %in% Exp_G == T)
    merge_table3_f3 = merge_table3_f2[k3_filter,]
    ###### merge_table3_f3[which(merge_table3_f3$TF == "Zic3"),] ######
    ######
    ###### output the TF-peak-Gene network and TF-Gene network ####
    ###### add pos and neg ####
    merge_table3_f3$Class = "pos"
    k_neg = which(merge_table3_f3$GG_Corr < 0)
    merge_table3_f3$Class[k_neg] = 'neg'
    #######
    merge_table3_f3_3col = merge_table3_f3[,c("TF","Peak","Gene","Class","celltype")]
    merge_table3_f3_2col = merge_table3_f3[,c("TF","Gene","Class","celltype")]
    ###
    merge_table3_f3_3col_index = paste(merge_table3_f3_3col$TF,merge_table3_f3_3col$Peak,merge_table3_f3_3col$Gene)
    merge_table3_f3_3col_cl = merge_table3_f3_3col[!duplicated(merge_table3_f3_3col_index),]
    ###
    merge_table3_f3_2col_index = paste(merge_table3_f3_2col$TF,merge_table3_f3_2col$Gene)
    merge_table3_f3_2col_cl = merge_table3_f3_2col[!duplicated(merge_table3_f3_2col_index),]
    ####
    return(list(Ori=merge_table3_f3,TPG=merge_table3_f3_3col_cl,TG=merge_table3_f3_2col_cl))
}




########--------------filter GRNs #####------------------------------------------------------------------------------------------------
########--------------filter GRNs #####------------------------------------------------------------------------------------------------
########--------------filter GRNs #####------------------------------------------------------------------------------------------------
########--------------filter GRNs #####------------------------------------------------------------------------------------------------

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)

########-----for Human samples------###

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
Human_GRNs_list = readRDS("Human_GRNs_list")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
Mouse_GRNs_list = readRDS("Mouse_GRNs_list")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
Zebrafish_GRNs_list <- readRDS("Zebrafish_GRNs_list")


########---------

Identify_expressed_TF <- function(count_list){
    ########
    out_list <- list()
    for(i in 1:length(count_list)){
        tmp = count_list[[i]][[1]]
        #### 
        tmp_rowmax = apply(tmp,1,max)
        #### tmp_rowmax[which(rownames(tmp) == 'NFIX')]
        k = which(tmp_rowmax > 0.5)
        print(length(k))
        ####
        exp_G = rownames(tmp)[k]
        out_list <- c(out_list,list(exp_G))
    }
    names(out_list) = names(count_list)
    return(out_list)
}

Filter_GRNs <- function(GRNs_list,Gene_list){
    #######
    celltype = names(GRNs_list)
    #######
    for(i in 1:length(celltype)){
        print(celltype[i])
        ######
        tmp_GRNs = GRNs_list[[celltype[i]]][[3]]
        tmp_Gene = Gene_list[[celltype[i]]]
        ######
        tmp_GRNs_cl = tmp_GRNs[which(tmp_GRNs$TF %in% tmp_Gene == T),]
        #######
        GRNs_list[[i]] = tmp_GRNs_cl
        print(dim(tmp_GRNs))
        print(dim(tmp_GRNs_cl))
    }
    return(GRNs_list)
}

###### load the Human pseudo count #######

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

Human_expressed = Identify_expressed_TF(Human_pseudo_count)
cell_types <- unique(sub("_[FM]$", "", names(Human_expressed)))

# 对每个细胞类型，合并 F 和 M 两个元素
Human_expressed_combined <- setNames(
  lapply(cell_types, function(ct) {
    c(Human_expressed[[paste0(ct, "_F")]],
      Human_expressed[[paste0(ct, "_M")]])
  }),
  cell_types
)

########
########
Human_GRNs_list_cl = Filter_GRNs(Human_GRNs_list,Human_expressed_combined)


########
########

Identify_expressed_TF <- function(count_list){
    ########
    out_list <- list()
    for(i in 1:length(count_list)){
        tmp = count_list[[i]]
        #### 
        tmp_rowmax = apply(tmp,1,max)
        #### tmp_rowmax[which(rownames(tmp) == 'NFIX')]
        k = which(tmp_rowmax > 0.5)
        print(length(k))
        ####
        exp_G = rownames(tmp)[k]
        out_list <- c(out_list,list(exp_G))
    }
    names(out_list) = names(count_list)
    return(out_list)
}

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
Mouse_pseudo_count <- readRDS("Mouse_pseudo_count_2025")
Mouse_expressed = Identify_expressed_TF(Mouse_pseudo_count)
Mouse_GRNs_list_cl = Filter_GRNs(Mouse_GRNs_list,Mouse_expressed)

setwd("/zp1/data/plyu3/Aging_2025_Zebrafish_Final")
Zebrafish_pseudo_count <- readRDS("Zebrafish_pseudo_count_2025")
Zebrafish_expressed = Identify_expressed_TF(Zebrafish_pseudo_count)
Zebrafish_GRNs_list_cl = Filter_GRNs(Zebrafish_GRNs_list,Zebrafish_expressed)


######## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######## repress Old Young activate Old Young ####
########
######## 4 classes ######
########
######## for each cell type, load the DEGs !! #######
######## for each cell type, load the aging clocks !!!! ######
######## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
names(Zebrafish_GRNs_list)
"RGC" "MG" "HC" "BC" "AC" "Cone" "Rod"
########
######## load the Zebrafish aging DEGs ########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########

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
Zebrafish_DOWN = kc[which(kc$Class=="DOWN"),]

Zebrafish_Old_DEGs = split(Zebrafish_UP$Gene,Zebrafish_UP$CT)
Zebrafish_Young_DEGs = split(Zebrafish_DOWN$Gene,Zebrafish_DOWN$CT)

########
######## load the Mouse aging DEGs 
########

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
Mouse_DOWN = kc[which(kc$Class=="DOWN"),]

#####
Mouse_Old_DEGs = split(Mouse_UP$Gene,Mouse_UP$CT)
Mouse_Young_DEGs = split(Mouse_DOWN$Gene,Mouse_DOWN$CT)

#####----------------------------------------------------------------------------------------------------------------------------------------------------------------
##### for Human load the human DEGs #######
#####----------------------------------------------------------------------------------------------------------------------------------------------------------------
#####----------------------------------------------------------------------------------------------------------------------------------------------------------------


setwd("/zp1/data/share/Human_aging_new")
load("Human_DEGs_Plot_Kmeans_order")

#######-----for M ######----------------------------------------------------------------------------------------------------------------------------------------------------------------

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
Human_DOWN = kc[which(kc$Class=="DOWN"),]


######
Human_Old_DEGs = split(Human_UP$Gene,Human_UP$CT)
Human_Young_DEGs = split(Human_DOWN$Gene,Human_DOWN$CT)

setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
save(Human_Old_DEGs,file="Human_Old_DEGs")
save(Human_Young_DEGs,file="Human_Young_DEGs")
save(Mouse_Old_DEGs,file="Mouse_Old_DEGs")
save(Mouse_Young_DEGs,file="Mouse_Young_DEGs")
save(Zebrafish_Old_DEGs,file="Zebrafish_Old_DEGs")
save(Zebrafish_Young_DEGs,file="Zebrafish_Young_DEGs")

######
###### Next we will load the clock genes #######
######




###### for Mouse ####------------------------------------------------------------------------------------------------------------------------
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

MG = Get_old_G_from_model(Mouse_MG_model$model)
RGC = Get_old_G_from_model(Mouse_RGC_model$model)
AC = Get_old_G_from_model(Mouse_AC_model$model)
HC = Get_old_G_from_model(Mouse_HC_model$model)
BC = Get_old_G_from_model(Mouse_BC_model$model)
RPE = Get_old_G_from_model(Mouse_RPE_model$model)
Rod = Get_old_G_from_model(Mouse_Rod_model$model)
Cone = Get_old_G_from_model(Mouse_Cone_model$model)
Microglia = Get_old_G_from_model(Mouse_Microglia_model$model)

G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

Mouse_Old_Clock = split(G,CT)

######

MG = Get_young_G_from_model(Mouse_MG_model$model)
RGC = Get_young_G_from_model(Mouse_RGC_model$model)
AC = Get_young_G_from_model(Mouse_AC_model$model)
HC = Get_young_G_from_model(Mouse_HC_model$model)
BC = Get_young_G_from_model(Mouse_BC_model$model)
RPE = Get_young_G_from_model(Mouse_RPE_model$model)
Rod = Get_young_G_from_model(Mouse_Rod_model$model)
Cone = Get_young_G_from_model(Mouse_Cone_model$model)
Microglia = Get_young_G_from_model(Mouse_Microglia_model$model)

G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

Mouse_Young_Clock = split(G,CT)

###### for Zebrafish ####------------------------------------------------------------------------------------------------------------------------


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

MG = Get_old_G_from_model(Zebrafish_MG_model$model)
RGC = Get_old_G_from_model(Zebrafish_RGC_model$model)
AC = Get_old_G_from_model(Zebrafish_AC_model$model)
HC = Get_old_G_from_model(Zebrafish_HC_model$model)
BC = Get_old_G_from_model(Zebrafish_BC_model$model)
RPE = Get_old_G_from_model(Zebrafish_RPE_model$model)
Rod = Get_old_G_from_model(Zebrafish_Rod_model$model)
Cone = Get_old_G_from_model(Zebrafish_Cone_model$model)
Microglia = Get_old_G_from_model(Zebrafish_Microglia_model$model)


G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

Zebrafish_Old_Clock = split(G,CT)

MG = Get_young_G_from_model(Zebrafish_MG_model$model)
RGC = Get_young_G_from_model(Zebrafish_RGC_model$model)
AC = Get_young_G_from_model(Zebrafish_AC_model$model)
HC = Get_young_G_from_model(Zebrafish_HC_model$model)
BC = Get_young_G_from_model(Zebrafish_BC_model$model)
RPE = Get_young_G_from_model(Zebrafish_RPE_model$model)
Rod = Get_young_G_from_model(Zebrafish_Rod_model$model)
Cone = Get_young_G_from_model(Zebrafish_Cone_model$model)
Microglia = Get_young_G_from_model(Zebrafish_Microglia_model$model)


G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

Zebrafish_Young_Clock = split(G,CT)

########
########------Next for Human ###########
########

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

MG = c(Get_old_G_from_model(Human_MG_model_F$model),Get_old_G_from_model(Human_MG_model_M$model))
RGC = c(Get_old_G_from_model(Human_RGC_model_F$model),Get_old_G_from_model(Human_RGC_model_M$model))
AC = c(Get_old_G_from_model(Human_AC_model_F$model),Get_old_G_from_model(Human_AC_model_M$model))
HC = c(Get_old_G_from_model(Human_HC_model_F$model),Get_old_G_from_model(Human_HC_model_M$model))
BC = c(Get_old_G_from_model(Human_BC_model_F$model),Get_old_G_from_model(Human_BC_model_M$model))
RPE = c(Get_old_G_from_model(Human_RPE_model_F$model),Get_old_G_from_model(Human_RPE_model_M$model))
Rod = c(Get_old_G_from_model(Human_Rod_model_F$model),Get_old_G_from_model(Human_Rod_model_M$model))
Cone = c(Get_old_G_from_model(Human_Cone_model_F$model),Get_old_G_from_model(Human_Cone_model_M$model))
Microglia = c(Get_old_G_from_model(Human_Microglia_model_F$model),Get_old_G_from_model(Human_Microglia_model_M$model))



G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

Human_Old_Clock = split(G,CT)

######
######
######

MG = c(Get_young_G_from_model(Human_MG_model_F$model),Get_young_G_from_model(Human_MG_model_M$model))
RGC = c(Get_young_G_from_model(Human_RGC_model_F$model),Get_young_G_from_model(Human_RGC_model_M$model))
AC = c(Get_young_G_from_model(Human_AC_model_F$model),Get_young_G_from_model(Human_AC_model_M$model))
HC = c(Get_young_G_from_model(Human_HC_model_F$model),Get_young_G_from_model(Human_HC_model_M$model))
BC = c(Get_young_G_from_model(Human_BC_model_F$model),Get_young_G_from_model(Human_BC_model_M$model))
RPE = c(Get_young_G_from_model(Human_RPE_model_F$model),Get_young_G_from_model(Human_RPE_model_M$model))
Rod = c(Get_young_G_from_model(Human_Rod_model_F$model),Get_young_G_from_model(Human_Rod_model_M$model))
Cone = c(Get_young_G_from_model(Human_Cone_model_F$model),Get_young_G_from_model(Human_Cone_model_M$model))
Microglia = c(Get_young_G_from_model(Human_Microglia_model_F$model),Get_young_G_from_model(Human_Microglia_model_M$model))



G = c(MG,RGC,AC,HC,BC,RPE,Rod,Cone,Microglia)
CT = rep(c("MG","RGC","AC","HC","BC","RPE","Rod","Cone","Microglia"),c(length(MG),length(RGC),length(AC),length(HC),length(BC),length(RPE),length(Rod),length(Cone),length(Microglia)))

Human_Young_Clock = split(G,CT)

######
######
######


setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
save(Human_Old_Clock,file="Human_Old_Clock")
save(Human_Young_Clock,file="Human_Young_Clock")
save(Mouse_Old_Clock,file="Mouse_Old_Clock")
save(Mouse_Young_Clock,file="Mouse_Young_Clock")
save(Zebrafish_Old_Clock,file="Zebrafish_Old_Clock")
save(Zebrafish_Young_Clock,file="Zebrafish_Young_Clock")
















###### for Human ####------------------------------------------------------------------------------------------------------------------------






#######
#######
#######

Make_the_barplots_for_Figure2 <- function(GRNs_list,UP_list,DOWN_list){
    ##########
    Res_table = list()
    for(i in 1:length(GRNs_list)){
        ######
        TAG = names(GRNs_list)[i]
        ######
        k = which(names(GRNs_list) == TAG)
        GRNs_tmp = GRNs_list[[k]]
        GRNs_tmp_3 = GRNs_tmp
        ######
        k2 = which(names(UP_list) == TAG)
        k3 = which(names(DOWN_list) == TAG)
        UP_G = UP_list[[k2]]
        DOWN_G = DOWN_list[[k3]]
        ######
        GRNs_tmp_3$AGE = "NODIFF"
        k4 = which(GRNs_tmp_3$Gene %in% UP_G == T)
        k5 = which(GRNs_tmp_3$Gene %in% DOWN_G == T)
        GRNs_tmp_3$AGE[k4] = "Old"
        GRNs_tmp_3$AGE[k5] = "Young"
        ######
        GRNs_tmp_3$combined_Tag = paste(GRNs_tmp_3$Class,GRNs_tmp_3$AGE,sep="_")
        ######
        RES = data.frame(table(GRNs_tmp_3$combined_Tag))
        ######
        RES$celltype = TAG
        ######
        k6 = grep("_NODIFF",RES$Var1)
        RES = RES[-k6,]
        ###
        Res_table <- c(Res_table,list(RES))
    }
    Res_table <- do.call("rbind",Res_table)
    ###########
    return(Res_table)
    ##########
}


########
########
GRNs_list = Zebrafish_GRNs_list_cl
UP_list = Zebrafish_Old
DOWN_list = Zebrafish_Young

Zebrafish_GRNs_RES = Make_the_barplots_for_Figure2(GRNs_list,UP_list,DOWN_list)

########
library(ggplot2)
########

Zebrafish_GRNs_RES$celltype = factor(Zebrafish_GRNs_RES$celltype,levels=c("MG","RGC","AC","HC","Rod","Cone","BC"))
Zebrafish_GRNs_RES$Var1 = factor(Zebrafish_GRNs_RES$Var1,levels=c("pos_Old","pos_Young","neg_Old","neg_Young"))

########
library(ggplot2)
ggplot(Zebrafish_GRNs_RES,aes(x=celltype,y=Freq,fill=Var1)) + geom_bar(stat = "identity",color="black") + theme_classic() + xlab("") + ylab("") + scale_y_continuous(expand = c(0, 0), limits=c(0,30000)) + theme(
    panel.border = element_rect(color = "black", fill = NA),  # 黑色框线
    panel.background = element_blank()) + scale_fill_manual(values = c(pos_Old = "#E41A1C",pos_Young = "#377EB8",neg_Old = "lightblue",neg_Young = "pink"))
ggsave("Zebrafish_GRNs_bar.png",height=3,width=4) 



########
######## load the Mouse DEGs ###################################---------------------------------------------------------------------------------------------------------
########

# setwd("/zp1/data/plyu3/Aging_2025_Mouse_Final")
# load("Mouse_DEGs_Plot_Kmeans_order")

# kc = Mouse_DEGs_Plot_Kmeans_order
# Upclusters = c(5,6,7,8)
# Downclusters = c(1,2,3,4)

# kc$CT = sapply(strsplit(kc$genes,split="__"),function(x) x[[1]])
# kc$Gene = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])

# kc$Class = "Unknown"
# kc$Class[which(kc$cluster %in% Upclusters)] = "UP"
# kc$Class[which(kc$cluster %in% Downclusters)] = "DOWN"

# #####
# Mouse_UP = kc[which(kc$Class=="UP"),]
# Mouse_UP_list = split(Mouse_UP,Mouse_UP$CT)

# Mouse_DOWN = kc[which(kc$Class=="DOWN"),]
# Mouse_DOWN_list = split(Mouse_DOWN,Mouse_DOWN$CT)


Mouse_Old = split(kc_index$G,kc_index$CT)
Mouse_Young = split(kc_index$G,kc_index$CT)



GRNs_list = Mouse_GRNs_list_cl
UP_list = Mouse_Old
DOWN_list = Mouse_Young

Mouse_GRNs_RES = Make_the_barplots_for_Figure2(GRNs_list,UP_list,DOWN_list)


########

########
library(ggplot2)
########

Mouse_GRNs_RES$celltype = factor(Mouse_GRNs_RES$celltype,levels=c("MG","RGC","AC","HC","Rod","Cone","BC"))
Mouse_GRNs_RES$Var1 = factor(Mouse_GRNs_RES$Var1,levels=c("pos_Old","pos_Young","neg_Old","neg_Young"))


library(ggplot2)
ggplot(Mouse_GRNs_RES,aes(x=celltype,y=Freq,fill=Var1)) + geom_bar(stat = "identity",color="black") + theme_classic() + xlab("") + ylab("") + scale_y_continuous(expand = c(0, 0), limits=c(0,30000)) + theme(
    panel.border = element_rect(color = "black", fill = NA),  # 黑色框线
    panel.background = element_blank()) + scale_fill_manual(values = c(pos_Old = "#E41A1C",pos_Young = "#377EB8",neg_Old = "lightblue",neg_Young = "pink"))
ggsave("Mouse_GRNs_bar.png",height=3,width=4) 




########




########
######## Next for Human #######
########


# setwd("/zp1/data/share/Human_aging_new")
# load("Human_DEGs_Plot_Kmeans_order")

# #######-----for M ######

# kc = Human_DEGs_Plot_Kmeans_order
# Upclusters = c(5,6,7,8,9)
# Downclusters = c(1,2,3,4)

# #######

# kc$CT = sapply(strsplit(kc$genes,split="_"),function(x) x[[1]])
# kc$Gene = sapply(strsplit(kc$genes,split="__"),function(x) x[[2]])

# kc$Class = "Unknown"
# kc$Class[which(kc$cluster %in% Upclusters)] = "UP"
# kc$Class[which(kc$cluster %in% Downclusters)] = "DOWN"

# #####
# Human_UP = kc[which(kc$Class=="UP"),]
# Human_UP_list = split(Human_UP,Human_UP$CT)

# Human_DOWN = kc[which(kc$Class=="DOWN"),]
# Human_DOWN_list = split(Human_DOWN,Human_DOWN$CT)

#######

########
########

Human_Old = split(kc_index$G,kc_index$CT)
Human_Young = split(kc_index$G,kc_index$CT)


GRNs_list = Human_GRNs_list_cl
UP_list = Human_Old
DOWN_list = Human_Young

Human_GRNs_RES = Make_the_barplots_for_Figure2(GRNs_list,UP_list,DOWN_list)

########
library(ggplot2)
########

Human_GRNs_RES$celltype = factor(Human_GRNs_RES$celltype,levels=c("MG","RGC","AC","HC","Rod","Cone","BC"))
Human_GRNs_RES$Var1 = factor(Human_GRNs_RES$Var1,levels=c("pos_Old","pos_Young","neg_Old","neg_Young"))

########
library(ggplot2)
ggplot(Human_GRNs_RES,aes(x=celltype,y=Freq,fill=Var1)) + geom_bar(stat = "identity",color="black") + theme_classic() + xlab("") + ylab("") + scale_y_continuous(expand = c(0, 0), limits=c(0,40000)) + theme(
    panel.border = element_rect(color = "black", fill = NA),  # 黑色框线
    panel.background = element_blank()) + scale_fill_manual(values = c(pos_Old = "#E41A1C",pos_Young = "#377EB8",neg_Old = "lightblue",neg_Young = "pink"))
ggsave("Human_GRNs_bar.png",height=3,width=4) 



##########
##########---------------------------------------------------------------------------------------------------------------------------------------------------------
########## 只计算 MG cells ####
##########---------------------------------------------------------------------------------------------------------------------------------------------------------
##########


Human_List_Add_DEGs = Make_the_barplots_for_Figure3(Human_GRNs_list_cl,Human_Old_DEGs,Human_Young_DEGs)
Human_List_Add_Clock = Make_the_barplots_for_Figure3(Human_GRNs_list_cl,Human_Old_Clock,Human_Young_Clock)

Zebrafish_List_Add_DEGs = Make_the_barplots_for_Figure3(Zebrafish_GRNs_list_cl,Zebrafish_Old_DEGs,Zebrafish_Young_DEGs)
Zebrafish_List_Add_Clock = Make_the_barplots_for_Figure3(Zebrafish_GRNs_list_cl,Zebrafish_Old_Clock,Zebrafish_Young_Clock)

Mouse_List_Add_DEGs = Make_the_barplots_for_Figure3(Mouse_GRNs_list_cl,Mouse_Old_DEGs,Mouse_Young_DEGs)
Mouse_List_Add_Clock = Make_the_barplots_for_Figure3(Mouse_GRNs_list_cl,Mouse_Old_Clock,Mouse_Young_Clock)


Make_the_barplots_for_Figure3 <- function(GRNs_list,UP_list,DOWN_list){
    ##########
    Res_table = list()
    for(i in 1:length(GRNs_list)){
        ######
        TAG = names(GRNs_list)[i]
        ######
        k = which(names(GRNs_list) == TAG)
        GRNs_tmp = GRNs_list[[k]]
        GRNs_tmp_3 = GRNs_tmp
        ######
        k2 = which(names(UP_list) == TAG)
        k3 = which(names(DOWN_list) == TAG)
        UP_G = UP_list[[k2]]
        DOWN_G = DOWN_list[[k3]]
        ######
        GRNs_tmp_3$AGE = "NODIFF"
        k4 = which(GRNs_tmp_3$Gene %in% UP_G == T)
        k5 = which(GRNs_tmp_3$Gene %in% DOWN_G == T)
        GRNs_tmp_3$AGE[k4] = "Old"
        GRNs_tmp_3$AGE[k5] = "Young"
        ######
        GRNs_tmp_3$combined_Tag = paste(GRNs_tmp_3$Class,GRNs_tmp_3$AGE,sep="_")
        ######
        ######
        ###
        Res_table <- c(Res_table,list(GRNs_tmp_3))
    }
    ###########
    names(Res_table) = names(GRNs_list)
    return(Res_table)
    ##########
}



#########
######### perform enrichment analysis ##########-------------------------------------------
#########


GRNs_list = Human_List_Add

Human_TF_enrich_DEGs = Perform_TF_enrich_FINAL(Human_List_Add_DEGs)
Human_TF_enrich_Clock = Perform_TF_enrich_FINAL(Human_List_Add_Clock)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
saveRDS(Human_TF_enrich,file="Human_TF_enrich")

Zebrafish_TF_enrich_DEGs = Perform_TF_enrich_FINAL(Zebrafish_List_Add_DEGs)
Zebrafish_TF_enrich_Clock = Perform_TF_enrich_FINAL(Zebrafish_List_Add_Clock)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
saveRDS(Zebrafish_TF_enrich,file="Zebrafish_TF_enrich")

Mouse_TF_enrich_DEGs = Perform_TF_enrich_FINAL(Mouse_List_Add_DEGs)
Mouse_TF_enrich_Clock = Perform_TF_enrich_FINAL(Mouse_List_Add_Clock)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
saveRDS(Mouse_TF_enrich,file="Mouse_TF_enrich")


#######
#######
library(writexl)
write_xlsx(Human_TF_enrich_DEGs, path = "Human_TF_enrich_DEGs.xlsx")
write_xlsx(Human_TF_enrich_Clock, path = "Human_TF_enrich_Clock.xlsx")
write_xlsx(Mouse_TF_enrich_DEGs, path = "Mouse_TF_enrich_DEGs.xlsx")
write_xlsx(Mouse_TF_enrich_Clock, path = "Mouse_TF_enrich_Clock.xlsx")
write_xlsx(Zebrafish_TF_enrich_DEGs, path = "Zebrafish_TF_enrich_DEGs.xlsx")
write_xlsx(Zebrafish_TF_enrich_Clock, path = "Zebrafish_TF_enrich_Clock.xlsx")





Perform_TF_enrich_FINAL <- function(GRNs_list){
    ######
    Total_res_tab = list()
    ######
    for(i in 1:length(GRNs_list)){
        print(i)
        #######
        tmp_res = list()
        #######
        tmp_GRNs = GRNs_list[[i]]
        #######
        tmp_all_TFs = names(table(tmp_GRNs$TF))
        for(j in 1:length(tmp_all_TFs)){
            tmp_TFs = tmp_all_TFs[j]
            ######
            tmp_table = data.frame(TF = tmp_TFs)
            ###### Next we will test for Old ######
            ###### pos_Old ####
            successes_in_population = length(which(tmp_GRNs$combined_Tag == "pos_Old"))
            successes_in_sample = length(which(tmp_GRNs$TF == tmp_TFs & tmp_GRNs$combined_Tag == "pos_Old"))
            population_size = length(tmp_GRNs$TF)
            sample_size = length(which(tmp_GRNs$TF == tmp_TFs))
            ######
            p_value <- phyper(successes_in_sample - 1, successes_in_population, population_size - successes_in_population, sample_size, lower.tail = FALSE)
            coverage = successes_in_sample / sample_size
            ######
            tmp_table$pos_Old__p = p_value
            tmp_table$pos_Old__cov = coverage
            ###### Young Young Young #####
            successes_in_population = length(which(tmp_GRNs$combined_Tag == "pos_Young"))
            successes_in_sample = length(which(tmp_GRNs$TF == tmp_TFs & tmp_GRNs$combined_Tag == "pos_Young"))
            population_size = length(tmp_GRNs$TF)
            sample_size = length(which(tmp_GRNs$TF == tmp_TFs))
            ######
            p_value <- phyper(successes_in_sample - 1, successes_in_population, population_size - successes_in_population, sample_size, lower.tail = FALSE)
            coverage = successes_in_sample / sample_size
            ######
            tmp_table$pos_Young__p = p_value
            tmp_table$pos_Young__cov = coverage
            #######-------------------------------------------------------
            successes_in_population = length(which(tmp_GRNs$combined_Tag == "neg_Old"))
            successes_in_sample = length(which(tmp_GRNs$TF == tmp_TFs & tmp_GRNs$combined_Tag == "neg_Old"))
            population_size = length(tmp_GRNs$TF)
            sample_size = length(which(tmp_GRNs$TF == tmp_TFs))
            ######
            p_value <- phyper(successes_in_sample - 1, successes_in_population, population_size - successes_in_population, sample_size, lower.tail = FALSE)
            coverage = successes_in_sample / sample_size
            ######
            tmp_table$neg_Old__p = p_value
            tmp_table$neg_Old__cov = coverage
            #######-------------------------------------------------------
            successes_in_population = length(which(tmp_GRNs$combined_Tag == "neg_Young"))
            successes_in_sample = length(which(tmp_GRNs$TF == tmp_TFs & tmp_GRNs$combined_Tag == "neg_Young"))
            population_size = length(tmp_GRNs$TF)
            sample_size = length(which(tmp_GRNs$TF == tmp_TFs))
            ######
            p_value <- phyper(successes_in_sample - 1, successes_in_population, population_size - successes_in_population, sample_size, lower.tail = FALSE)
            coverage = successes_in_sample / sample_size
            ######
            tmp_table$neg_Young__p = p_value
            tmp_table$neg_Young__cov = coverage
            ########
            tmp_res <- c(tmp_res,list(tmp_table))
            ######
        }
        tmp_res_merge = do.call("rbind",tmp_res)
        Total_res_tab = c(Total_res_tab,list(tmp_res_merge))
    }
    ####
    #### head(tmp_res_merge[order(tmp_res_merge$Young_p,decreasing=F),])
    ####
    names(Total_res_tab) = names(GRNs_list)
    ####
    return(Total_res_tab)
}



#############


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R
library(ArchR)


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Human_aging_scATACseq/Human_aging_scATACseq")
Human_TF_enrich <- readRDS(file="Human_TF_enrich")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/fish_aging_ATAC_2024May/ATAC_out")
Zebrafish_TF_enrich <- readRDS(file="Zebrafish_TF_enrich")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Aging_Mouse/Aging_arrow")
Mouse_TF_enrich <- readRDS(file="Mouse_TF_enrich")


#######
TAG1 = "MG"
TAG2 = "pos_Young"

H = Human_TF_enrich_Clock
M = Mouse_TF_enrich_Clock
Z = Zebrafish_TF_enrich_Clock

HMZ_ortholog_combined <- readRDS("/zp1/data/share/Human_aging_new/HMZ_ortholog_combined_2025")


Prepare_the_dot_plot <- function(H,M,Z,TAG1,TAG2){
    ######
    H_cl = H[[which(names(H) == TAG1)]]
    M_cl = M[[which(names(M) == TAG1)]]
    Z_cl = Z[[which(names(Z) == TAG1)]]
    ######
    ######
    H_clcl = H_cl[,c(1,grep(TAG2,colnames(H_cl)))]
    M_clcl = M_cl[,c(1,grep(TAG2,colnames(M_cl)))]
    Z_clcl = Z_cl[,c(1,grep(TAG2,colnames(Z_cl)))]
    ######
    H_cl1 = H_clcl$TF[which(H_clcl[,2] < 0.05)]
    M_cl1 = M_clcl$TF[which(M_clcl[,2] < 0.05)]
    Z_cl1 = Z_clcl$TF[which(Z_clcl[,2] < 0.05)]
    ######
    ###### find overlaps between HMZ #######
    ######
    HM_overlap = HMZ_ortholog_combined$HM[which(HMZ_ortholog_combined$HM$human %in% H_cl1 == T & HMZ_ortholog_combined$HM$mouse %in% M_cl1 == T),]
    HZ_overlap = HMZ_ortholog_combined$HZ[which(HMZ_ortholog_combined$HZ$human %in% H_cl1 == T & HMZ_ortholog_combined$HZ$zebrafish %in% Z_cl1 == T),]
    MZ_overlap = HMZ_ortholog_combined$MZ[which(HMZ_ortholog_combined$MZ$mouse %in% M_cl1 == T & HMZ_ortholog_combined$MZ$zebrafish %in% Z_cl1 == T),]
    ######
    HMZ_overlap = HMZ_ortholog_combined$HMZ[which(HMZ_ortholog_combined$HMZ$human %in% H_cl1 == T & HMZ_ortholog_combined$HMZ$mouse %in% M_cl1 == T & HMZ_ortholog_combined$HMZ$zebrafish %in% Z_cl1 == T),]
    ######
    ###### remove duplicates #######
    HM_overlap$tag = "NO"
    for(i in 1:length(HM_overlap$tag)){
        k = which(HM_overlap$human[i] %in% HMZ_overlap$human== T & HM_overlap$mouse[i] %in% HMZ_overlap$mouse == T)
        if(length(k) > 0){
            HM_overlap$tag[i] = "YES"
        }
    }
    HZ_overlap$tag = "NO"
    for(i in 1:length(HZ_overlap$tag)){
        k = which(HZ_overlap$human[i] %in% HMZ_overlap$human== T & HZ_overlap$zebrafish[i] %in% HMZ_overlap$zebrafish == T)
        if(length(k) > 0){
            HZ_overlap$tag[i] = "YES"
        }
    }
    MZ_overlap$tag = "NO"
    for(i in 1:length(MZ_overlap$tag)){
        k = which(MZ_overlap$mouse[i] %in% HMZ_overlap$mouse== T & MZ_overlap$zebrafish[i] %in% HMZ_overlap$zebrafish == T)
        if(length(k) > 0){
            MZ_overlap$tag[i] = "YES"
        }
    }
    ####
    MZ_overlap = MZ_overlap[MZ_overlap$tag=="NO",]
    HZ_overlap = HZ_overlap[HZ_overlap$tag=="NO",]
    HM_overlap = HM_overlap[HM_overlap$tag=="NO",]
    ###### covert to ggplot2 !!!! #####
    ######
    ###### y TF, x = specis, pvalue, cov ######
    ######
    res2 = list()
    if(dim(HMZ_overlap)[1] > 0){
        ####
        index = HMZ_overlap$HMZ_index
        res = data.frame(index=rep(index,each=3),sp=rep(c("H","M","Z"),dim(HMZ_overlap)[1]))
        res$genes = res$index
        res$pvalue = 1
        res$cov = 1
        ####
        HMZ_res = res
        HMZ_res$class = "HMZ"
        res2 = c(res2,list(HMZ_res))
    }
    #######
    #######
    #######
    if(dim(HM_overlap)[1] > 0){
        ####
        index = HM_overlap$HM_index
        res = data.frame(index=rep(index,each=2),sp=rep(c("H","M"),dim(HM_overlap)[1]))
        res$genes = res$index
        res$pvalue = 1
        res$cov = 1
        ####
        HM_res = res
        HM_res$class = "HM"
        res2 = c(res2,list(HM_res))
    }
    if(dim(HZ_overlap)[1] > 0){
        ####
        index = HZ_overlap$HZ_index
        res = data.frame(index=rep(index,each=2),sp=rep(c("H","Z"),dim(HZ_overlap)[1]))
        res$genes = res$index
        res$pvalue = 1
        res$cov = 1
        ####
        HZ_res = res
        HZ_res$class="HZ"
        res2 = c(res2,list(HZ_res))
    }
    if(dim(MZ_overlap)[1] > 0){
        ####
        index = MZ_overlap$MZ_index
        res = data.frame(index=rep(index,each=2),sp=rep(c("M","Z"),dim(MZ_overlap)[1]))
        res$genes = res$index
        res$pvalue = 1
        res$cov = 1
        ####
        MZ_res = res
        MZ_res$class="MZ"
        res2 = c(res2,list(MZ_res))
    }
    ###########
    ###########
    res_total = do.call('rbind',res2)
    ######
    for(i in 1:length(res_total$index)){
        #########
        genes_tmp = res_total$genes[i]
        tag1 = res_total$sp[i]
        tag2 = res_total$class[i]
        ########
        chars <- strsplit(tag2, "")[[1]]
        pos <- match(tag1, chars)
        ########
        res_total$genes[i] = strsplit(genes_tmp,split="__")[[1]][pos]
    }
    #####
    ##### #### #####
    #####
    for(i in 1:length(res_total$index)){
        #########
        genes_tmp = res_total$genes[i]
        tag1 = res_total$sp[i]
        ########
        if(tag1 == "H"){
            m1 = which(H_clcl$TF == genes_tmp)
            res_total$pvalue[i] = H_clcl[,2][m1]
            res_total$cov[i] = H_clcl[,3][m1]
        }
        if(tag1 == "M"){
            m1 = which(M_clcl$TF == genes_tmp)
            res_total$pvalue[i] = M_clcl[,2][m1]
            res_total$cov[i] = M_clcl[,3][m1]
            
        }
        if(tag1 == "Z"){
            m1 = which(Z_clcl$TF == genes_tmp)
            res_total$pvalue[i] = Z_clcl[,2][m1]
            res_total$cov[i] = Z_clcl[,3][m1] 
        }
    }
    ########
    res_total$sp = factor(res_total$sp,levels=c("Z","M","H"))
    res_total$logP = -log(res_total$pvalue)
    k = which(res_total$logP > 5)
    res_total$logP[k] = 5
    ####
    res_total$order = match(res_total$class,c("HMZ","HM","MZ","HZ"))
    return(res_total)
}


TAG1 = "Rod"
TAG2 = "pos_Old"

res_total = Prepare_the_dot_plot(H,M,Z,"MG","pos_Old")
res_total = Process_res2(res_total)

res_total = Prepare_the_dot_plot(H,M,Z,"Rod","pos_Old")
res_total = Process_res2(res_total)

res_total = Prepare_the_dot_plot(H,M,Z,"Rod","pos_Young")
res_total = Process_res2(res_total)

res_total = Prepare_the_dot_plot(H,M,Z,"Rod","neg_Young")
res_total = Process_res2(res_total)

res_total = Prepare_the_dot_plot(H,M,Z,"Rod","neg_Old")
res_total = Process_res2(res_total)
####
Process_res2 <- function(res_total){
    #######
    k1 = which(res_total$sp == "H")
    k2 = which(res_total$sp == "M")
    k3 = which(res_total$sp == "Z")
    #######
    res_total$sp = as.character(res_total$sp)
    res_total$sp[k1] = "Human"
    res_total$sp[k2] = "Mouse"
    res_total$sp[k3] = "Zebrafish"
    ########
    for(i in 1:length(res_total$index)){
        #########
        tmp = res_total$index[i]
        tmp = unlist(strsplit(tmp,split="__"))
        tmp = rev(tmp)
        tmp = paste(tmp,collapse="__")
        res_total$index[i] = tmp
    }
    #######
    res_total$sp = factor(res_total$sp,levels=c("Zebrafish","Mouse","Human"))
    return(res_total)
}


####
library(ggplot2)
ggplot(res_total,aes(x=sp,y=reorder(index,order,decreasing=T))) + geom_point(aes(size= cov,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='red', guide = guide_colorbar(reverse = FALSE, order = 1)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("HMZ_MG_overlap.png",height=4.5,width=3.5)


#####
library(ggplot2)
ggplot(res_total,aes(x=sp,y=reorder(index,order,decreasing=T))) + geom_point(aes(size= cov,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='blue', guide = guide_colorbar(reverse = FALSE, order = 1)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("HMZ_MG_overlap.png",height=4,width=3.5)
    

#####
library(ggplot2)
ggplot(res_total,aes(x=sp,y=reorder(index,order,decreasing=T))) + geom_point(aes(size= cov,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='pink', guide = guide_colorbar(reverse = FALSE, order = 1)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("HMZ_MG_overlap.png",height=4,width=3.5)


#####
library(ggplot2)
ggplot(res_total,aes(x=sp,y=reorder(index,order,decreasing=T))) + geom_point(aes(size= cov,color=logP)) + theme_classic() + scale_color_continuous(low='grey',high='lightblue', guide = guide_colorbar(reverse = FALSE, order = 1)) + scale_size_continuous(guide = guide_legend(reverse = FALSE, order = 2)) + xlab("") + ylab("") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),axis.line = element_line(color = "black")) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("HMZ_MG_overlap.png",height=4,width=3.5)




######
######
######

Output_GRNs:


for Zebrafish:

Zebrafish_List_Add_Rod = Zebrafish_List_Add$Rod
Zebrafish_List_Add_Rod_Plot = Zebrafish_List_Add_Rod[which(Zebrafish_List_Add_Rod$TF=="foxo1a"),]








######




    ######


}
























########
######## identify enriched TFs #############
########
######## How many TFs enriched in the MG cells ? ###############
########


























