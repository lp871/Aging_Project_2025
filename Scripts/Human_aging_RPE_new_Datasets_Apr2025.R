#############
#############
#############
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate /home/plyu3/anaconda3_new/envs/scvi-env
python

cd /zp1/data/plyu3/Human_New_Datasets

cd /zp1/data/share/Human_aging_new

cp /zp1/data/plyu3/Human_New_Datasets/Human_Retina_atlas.h5ad /zp1/data/share/Human_aging_new/
cp /zp1/data/plyu3/Human_New_Datasets/Human_Retina_atlas.h5ad /zp1/data/share/Human_aging_new/
cp /zp1/data/plyu3/Human_New_Datasets/Human_RPE_atlas.h5ad /zp1/data/share/Human_aging_new/

#############
#######
#############

import os 

os.chdir("/zp1/data/plyu3/Human_New_Datasets")

####### RPE altas ######
#######

import scanpy as sc

# 读取 h5ad 文件到 AnnData 对象中
adata = sc.read_h5ad("Human_RPE_atlas.h5ad")

#
os.chdir("/zp1/data/share/Human_aging_new")
adata2 = sc.read_h5ad("adata_majorclass_Astrocyte.h5ad")

#
adata2.var['index2'] = adata2.var.index

genes = adata2.var
File = "AllHuman_gene.tsv"
genes.to_csv(
        File,
        sep='\t', index=False, header=False
)

#


# 打印 AnnData 对象的简要信息
print(adata)

# 也可以输出更详细的信息，如细胞数和基因数：
print(f"细胞数：{adata.n_obs}")
print(f"基因数：{adata.n_vars}")

# 查看 obs 和 var 数据框的列信息
print("Observation (obs) 的列：", list(adata.obs.columns))
print("Variable (var) 的列：", list(adata.var.columns))

#
sc.pl.umap(adata, color='subclass')
sc.pl.umap(adata, color='subclass', save='Human_RPE_umap.png', show=False)

sc.pl.umap(adata, color='Age', save='AgeHuman_RPE_umap.png', show=False)

##
## 看一下每个年龄都有多少细胞，有多少样本 ######
##
import matplotlib.pyplot as plt

if 'Age' not in adata.obs.columns:
    print("数据中没有 'Age' 列，请检查 h5ad 文件。")
else:
    # 统计每个 Age 对应的细胞数，确保按照年龄大小排序（假设年龄为数值或可以排序）
    age_counts = adata.obs['Age'].value_counts().sort_index()
    
    # 创建条形图，横坐标为年龄，纵坐标为对应细胞数
    plt.figure(figsize=(10, 6))
    plt.bar(age_counts.index.astype(str), age_counts.values)
    plt.xlabel('Age')
    plt.ylabel('Single Cell Count')
    plt.title('Number of Single Cells by Age')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # 将图像保存为 PNG 文件，分辨率 dpi 可根据需要调整
    plt.savefig('age_barplot.png', dpi=300)
    plt.show()

######
######
######
if 'Age' not in adata.obs.columns or 'Donor' not in adata.obs.columns:
    print("数据中缺少 'Age' 或 'donor' 列，请检查 h5ad 文件。")
else:
    # 按 Age 分组，统计每个年龄下不重复 donor 的数量
    age_donor_counts = adata.obs.groupby('Age')['Donor'].nunique().sort_index()
    
    # 创建条形图
    plt.figure(figsize=(10, 6))
    plt.bar(age_donor_counts.index.astype(str), age_donor_counts.values)
    plt.xlabel('Age')
    plt.ylabel('Number of Donors')
    plt.title('Unique Donors per Age')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # 保存图像为 PNG 文件，dpi 可根据需要调整
    plt.savefig('age_donor_barplot.png', dpi=300)
    plt.show()


#######
#######
#######

#######
#######
import pandas as pd
# 替换为你的 Excel 文件路径
# 读取 Excel 文件，默认读取第一个 sheet
Sm_df = pd.read_excel("Human_aging_samples_list_2023.xlsx")
Big_df = pd.read_excel("Retina_snRNA_snATAC_donors.xlsx")

# 显示 DataFrame 的前几行，确认数据读取正确
print(df.head())

###### 比较这两组数据的 sampleid #########
###### 并且用 韦恩图 表示，图上表明 overlap 的个数 ########
######
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# 假设 Sm_df 和 Big_df 已经存在，如需从文件读取数据请自行补充相应代码
# 例如：
# Sm_df = pd.read_excel('Sm_df.xlsx')
# Big_df = pd.read_excel('Big_df.xlsx')

# 提取两个 DataFrame 中的 sampleid 转换为 set
set_sm = set(Sm_df['sampleid'])
set_big = set(Big_df['sampleid'])

# 计算各个区域的数量：
#   仅在 Sm_df 中的数量：set_sm - set_big
#   仅在 Big_df 中的数量：set_big - set_sm
#   两者交集中的数量：set_sm & set_big
only_sm = len(set_sm - set_big)
only_big = len(set_big - set_sm)
overlap = len(set_sm & set_big)

# 使用 venn2 绘制温恩图，并自动将数量标在图中
venn = venn2(subsets=(only_sm, only_big, overlap))

# 设置标题
plt.title('SampleID Overlap between Sm_df and Big_df')

# 保存图片为 PNG 文件，dpi 参数可根据需要调整
plt.savefig('sampleid_overlap_venn.png', dpi=100)

# 显示图形
plt.show()

########
######## 下一步 比较一下 donor id ##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########


plt.clf()

# 提取 donor 列转换成集合（确保每个 donor 只计数一次）
donor_set_sm = set(Sm_df['donor'])
donor_set_big = set(Big_df['donor'])

# 分别计算：
# 仅在 Sm_df 中的 donor 数量
only_sm = len(donor_set_sm - donor_set_big)
# 仅在 Big_df 中的 donor 数量
only_big = len(donor_set_big - donor_set_sm)
# 两个 DataFrame 都存在的 donor 数量
overlap = len(donor_set_sm & donor_set_big)

# 使用 venn2 绘制温恩图，subsets 参数的顺序为 (只在第一个集合中, 只在第二个集合中, 两个集合均有)
venn_diagram = venn2(subsets=(only_sm, only_big, overlap), set_labels=('Sm_df', 'Big_df'))

# 添加标题
plt.title("Donor Overlap between Sm_df and Big_df")

# 保存图形为 PNG 文件，dpi 参数可根据需要调整
plt.savefig("donor_overlap_venn.png", dpi=300)

#####
##### 打印一下每个年龄的 Donor 个数 ######
#####

#####
#####

df = Sm_df

donor_unique_counts = df.groupby("age")["donor"].nunique().sort_index()

plt.figure(figsize=(10, 6))
plt.bar(donor_unique_counts.index.astype(str), donor_unique_counts.values)
plt.xlabel("Age")
plt.ylabel("Unique Donor Count")
plt.title("Distribution of Unique Donors by Age")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("unique_donor_distribution.png", dpi=300)
plt.show()


########
########
########
df = Big_df

print(df.columns)

donor_unique_counts = df.groupby("age")["donor"].nunique().sort_index()

plt.figure(figsize=(10, 6))
plt.bar(donor_unique_counts.index.astype(str), donor_unique_counts.values)
plt.xlabel("Age")
plt.ylabel("Unique Donor Count")
plt.title("Distribution of Unique Donors by Age")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("unique_donor_distribution.png", dpi=300)
plt.show()







##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##########
########## 重新 DEGs and aging clocks ！！！！ #########
########## New human #######
##########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda env list
conda activate seurat4


cd /zp1/data/share/Human_aging_new

python
import anndata

# 1. 读入 .h5ad 文件
adata = anndata.read_h5ad('Human_Retina_atlas.h5ad')

adata.raw
# 2. 列出 obs 中所有的列名
print(adata.obs.columns.tolist())

# 3. （可选）查看 obs 的前几行，确认各列的含义
print(adata.obs.head())

ages = adata.obs['donor_age'].unique()
print("Unique donor ages:", ages)

for i in ages:
    print(i)


######## 下一步是 把 带有 pooled 去掉 #####
######## 下一步是 把 带有 pooled 去掉 #####

adata.obs['donor_age'] = adata.obs['donor_age'].astype(str)
####mask1 = adata.obs['donor_age'].str.contains('pooled [43 years,16 years]', case=False, na=False)
####mask2 = adata.obs['donor_age'].str.contains('pooled [43 years,16 years]', case=False, na=False)

adata_filtered = adata
print(adata_filtered.obs['donor_age'].value_counts())

#########
#########
print(adata_filtered.obs['majorclass'].value_counts())

import re

adata_filtered.raw

for mc in adata_filtered.obs['majorclass'].unique():
    # 构造掩码
    mask = adata_filtered.obs['majorclass'] == mc
    # 子集并拷贝
    adata_sub = adata_filtered[mask].copy()
    # 对类别名做简单清洗，去掉空格和特殊字符，确保文件名合法
    safe_mc = re.sub(r'[^0-9A-Za-z_-]', '_', str(mc))
    # 定义输出文件名
    out_fname = f"adata_majorclass_{safe_mc}.h5ad"
    # 写出文件
    adata_sub.write(out_fname)
    print(f"Saved subset for majorclass = '{mc}' → {out_fname}")

    ########
    ########

######### 输出 mtx cell gene ####
python

import os
import numpy as np
import pandas as pd
import anndata
from scipy import sparse
from scipy.io import mmwrite

##########

import os 
os.chdir("/zp1/data/share/Human_aging_new")

######
########## No RPE ##############
######

Rod_adata = anndata.read_h5ad("adata_majorclass_Rod.h5ad")

###### 查看一下 Rod_adata 的 dims ######
######
cell_counts = Rod_adata.obs['sample_id'].value_counts()
# 2. 打印结果
print(cell_counts)

###### we set as 100 #######
valid_samples = cell_counts[cell_counts >= 100].index

Rod_adata_filtered = Rod_adata[Rod_adata.obs['sample_id'].isin(valid_samples)].copy()

raw_mat = Rod_adata_filtered.raw.X
if issparse(raw_mat):
    total_nz = raw_mat.count_nonzero()   # 等价于 len(raw_mat.data)
else:
    total_nz = np.count_nonzero(raw_mat)
print(f"Total non-zero entries in raw.X: {total_nz}")

int32_max = np.iinfo(np.int32).max
print("Exceeds 2^31-1? ", nnz > int32_max)

### 查找 Rod_adata_filtered.raw.X 里面非零值的数目 
###

cell_counts = Rod_adata_filtered.obs['sample_id'].value_counts()
####
import numpy as np
import scanpy as sc


Rod_adata_filtered_sample = subsample_cells(Rod_adata_filtered,sample_key='sample_id', max_cells=30000, random_state=123)
raw_mat = Rod_adata_filtered_sample.raw.X
total_nz = raw_mat.count_nonzero()   # 等价于 len(raw_mat.data)

####
####

Rod_adata_filtered_sample.write("adata_majorclass_Rod_sample.h5ad")




######
######
######

def subsample_cells(adata, sample_key='sample_id', max_cells=30000, random_state=123):
    """
    对每个 sample_id：
      - 如果细胞数 > max_cells，则无放回地随机抽取 max_cells 个细胞；
      - 否则保留该 sample_id 下所有细胞。
    随机种子默认为 123，返回一个新的 AnnData 子集。
    """
    obs = adata.obs
    groups = obs.groupby(sample_key).indices  # dict: sample_id -> list of obs 索引
    ###
    rng = np.random.default_rng(random_state)
    selected_idx = []
    ####
    for sample, idxs in groups.items():
        n = len(idxs)
        if n > max_cells:
            # 无放回地抽样 max_cells 个
            pick = rng.choice(idxs, size=max_cells, replace=False)
        else:
            pick = np.array(idxs, dtype=int)
        selected_idx.append(pick)
    ####
    all_idx = np.concatenate(selected_idx)
    all_idx.sort()  # 保持原始顺序
    ####
    return adata[all_idx].copy()



###

export_adata_to_mtx("adata_majorclass_Rod_sample.h5ad","Human_Rod")
export_adata_to_mtx("adata_majorclass_AC.h5ad","Human_AC")
export_adata_to_mtx("adata_majorclass_BC.h5ad","Human_BC")
export_adata_to_mtx("adata_majorclass_RGC.h5ad","Human_RGC")
export_adata_to_mtx("adata_majorclass_MG.h5ad","Human_MG")
export_adata_to_mtx("adata_majorclass_Cone.h5ad","Human_Cone")
export_adata_to_mtx("adata_majorclass_HC.h5ad","Human_HC")
export_adata_to_mtx("adata_majorclass_Microglia.h5ad","Human_Microglia")

##
## export_adata_to_mtx("Human_RPE_atlas.h5ad","Human_RPE")
## h5ad_path = "adata_majorclass_AC.h5ad"

def export_adata_to_mtx(h5ad_path: str, tag: str):
    """
    从 h5ad 文件中导出：
      - matrix.mtx：稀疏表达矩阵
      - genes.tsv：基因名列表
      - metadata.tsv：细胞元信息

    参数：
      h5ad_path: 输入的 .h5ad 文件全路径
      out_dir: 输出目录（不存在则创建）
    """
    # 1. 读入 AnnData
    adata = anndata.read_h5ad(h5ad_path)
    # 3. 导出矩阵（MTX 格式）
    #### raw_mat = adata.raw.X
    #### nz_values = raw_mat.data
    #### print("前10个非零表达值：", nz_values[:10])
    ####
    X = adata.raw.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    File = tag + ".mtx"
    mmwrite(File, X)
    # 4. 导出基因名（单列，无表头）
    adata.var['index2'] = adata.var.index
    genes = adata.var
    File = tag + "_gene.tsv"
    genes.to_csv(
        File,
        sep='\t', index=False, header=False
    )
    # 5. 导出细胞元信息（TSV，第一列为细胞条码/索引）
    File = tag + "_meta.tsv"
    #### 去掉 assay ###
    adata.obs.drop('assay', axis=1, inplace=True)
    adata.obs.to_csv(
        File,
        sep='\t', index=True, header=True
    )



########
######## --------------------------------------------------------------------------------------------------------------
########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4

########
########
R
# file: create_seurat.R
library(Matrix)
library(Seurat)

#' 从 MTX、基因列表和元数据创建 Seurat 对象
#'
#' @param matrix_file  MTX 文件路径（e.g. "matrix.mtx"）
#' @param genes_file   基因名文件路径（TSV，无表头，单列）
#' @param metadata_file 细胞元信息文件路径（TSV，带表头，首列为细胞条码）
#' @return Seurat 对象

matrix_file="Human_Rod.mtx"
genes_file="Human_HC_gene.tsv"
metadata_file="Human_HC_meta.tsv"


readMM_big <- function(path_mtx) {
  # 1. 先把 header（前三个数字所在行）单独读出来
  con      <- gzfile(path_mtx, open = "r")
  on.exit(close(con), add = TRUE)
  # 跳过所有以 % 开头的注释行
  header   <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line)==0) stop("Unexpected EOF before header")
    if (!grepl("^\\s*%", line)) {
      header <- line
      break
    }
  }
  parts    <- strsplit(header, "\\s+")[[1]]
  nrow     <- as.integer(parts[1])
  ncol     <- as.integer(parts[2])
  # 注意：parts[3] 我们不再用，因为它太大会溢出
  
  # 2. 用 data.table::fread 直接读剩下的三列（i, j, x）
  #    fread 会把数字读成 double，没问题
  #    skip = number of lines we already consumed (including header)
  #    这里我们先把 con 关了，再让 fread 从头 skip
  close(con)
  skip_n <- sum(grepl("^\\s*%", readLines(path_mtx, n = 100))) + 1L
  dt <- fread(
    path_mtx,
    skip        = skip_n,
    header      = FALSE,
    col.names   = c("i","j","x"),
    colClasses  = c("integer","integer","numeric")
  )
  
  # 3. 用 sparseMatrix 重建
  #    Matrix::sparseMatrix 会自己检查 i,j 的范围是否在 [1,nrow]、[1,ncol]
  mat <- sparseMatrix(
    i    = dt$i,
    j    = dt$j,
    x    = dt$x,
    dims = c(nrow, ncol)
  )
  return(mat)
}


create_seurat_from_files <- function(matrix_file, genes_file, metadata_file) {
  # 1. 读取稀疏矩阵
  counts <- readMM(matrix_file)
  ### lines <- readLines(matrix_file, n = 10)
  # 2. 读取基因名
  genes <- read.table(
    genes_file,
    sep = "\t", header = FALSE,
    stringsAsFactors = FALSE
  )
  # 3. 读取细胞元信息
  metadata <- read.table(
    metadata_file,
    sep = "\t", header = TRUE
  )
  #####------------#####
  ##### see HC ####
  # 4. 赋予行名和列名
  counts = t(counts)
  rownames(counts) <- genes$V7
  colnames(counts) <- metadata$X
  # 5. 构建 Seurat 对象
  seu_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata
  )
  #####
  return(seu_obj)
}


create_seurat_from_filesbig <- function(matrix_file, genes_file, metadata_file) {
  # 1. 读取稀疏矩阵
  counts <- readMM_big(matrix_file)
  ### lines <- readLines(matrix_file, n = 10)
  # 2. 读取基因名
  genes <- read.table(
    genes_file,
    sep = "\t", header = FALSE,
    stringsAsFactors = FALSE
  )
  # 3. 读取细胞元信息
  metadata <- read.table(
    metadata_file,
    sep = "\t", header = TRUE
  )
  #####------------#####
  ##### see HC ####
  # 4. 赋予行名和列名
  counts = t(counts)
  rownames(counts) <- genes$V7
  colnames(counts) <- metadata$X
  # 5. 构建 Seurat 对象
  seu_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata
  )
  #####
  return(seu_obj)
}

#######
setwd("/zp1/data/share/Human_aging_new")

#######
Human_Rod_seurat = create_seurat_from_files("Human_Rod.mtx","Human_Rod_gene.tsv","Human_Rod_meta.tsv")
saveRDS(Human_Rod_seurat,file="Human_Rod_seurat.rds")

Human_AC_seurat = create_seurat_from_files("Human_AC.mtx","Human_AC_gene.tsv","Human_AC_meta.tsv")
saveRDS(Human_AC_seurat,file="Human_AC_seurat.rds")

Human_BC_seurat = create_seurat_from_files("Human_BC.mtx","Human_BC_gene.tsv","Human_BC_meta.tsv")
saveRDS(Human_BC_seurat,file="Human_BC_seurat.rds")

#### RPE 的 code 不同 #####
#Human_RPE_seurat = create_seurat_from_files("Human_RPE.mtx","Human_RPE_gene.tsv","Human_RPE_meta.tsv")
#saveRDS(Human_RPE_seurat,file="Human_RPE_seurat.rds")

Human_Cone_seurat = create_seurat_from_files("Human_Cone.mtx","Human_Cone_gene.tsv","Human_Cone_meta.tsv")
saveRDS(Human_Cone_seurat,file="Human_Cone_seurat.rds")

Human_RGC_seurat = create_seurat_from_files("Human_RGC.mtx","Human_RGC_gene.tsv","Human_RGC_meta.tsv")
saveRDS(Human_RGC_seurat,file="Human_RGC_seurat.rds")

Human_HC_seurat = create_seurat_from_files("Human_HC.mtx","Human_HC_gene.tsv","Human_HC_meta.tsv")
saveRDS(Human_HC_seurat,file="Human_HC_seurat.rds")

Human_MG_seurat = create_seurat_from_files("Human_MG.mtx","Human_MG_gene.tsv","Human_MG_meta.tsv")
saveRDS(Human_MG_seurat,file="Human_MG_seurat.rds")

Human_Microglia_seurat = create_seurat_from_files("Human_Microglia.mtx","Human_Microglia_gene.tsv","Human_Microglia_meta.tsv")
saveRDS(Human_Microglia_seurat,file="Human_Microglia_seurat.rds")


########--------------------------------------------------------------------------------------------------------------------------------------------------------
######## 接下来我们找一下之前做 差异基因 Human 的 code ！！！！ #############
########
########
######## harmony the gene names !!! #############
########----------------------------------------------------------------------------------------------------------------------------------------------------------------

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
# file: create_seurat.R

R
library(Matrix)
library(Seurat)

setwd("/zp1/data/share/Human_aging_new")
Human_RPE = readRDS("Human_RPE_seurat.rds")

HC_meta = read.table(
    "Human_HC_gene.tsv",
    sep = "\t", header = FALSE,
    stringsAsFactors = FALSE
)
Rod_meta = read.table(
    "Human_Rod_gene.tsv",
    sep = "\t", header = FALSE,
    stringsAsFactors = FALSE
)

########

all.equal(HC_meta$V2,Rod_meta$V2)
unique(HC_meta$V2)

########
########

Human_HC = readRDS("Human_HC_seurat.rds")

Human_HC = Human_HC_seurat
########------------------------------------------------------------------------------------------------------
library(rtracklayer)
# 1. 导入 GTF
setwd("/zp1/data/share/Human_aging_new")
gtf_gr <- import("Homo_sapiens.GRCh38.98.gtf.gz", format = "gtf")
mcols_df <- as.data.frame(gtf_gr)
mcols_df = mcols_df[which(mcols_df$type == 'gene'),]
mcols_df = mcols_df[,c("gene_name","gene_id")]
dim(genes)

save(mcols_df,file="mcols_df_2025")

#########
#########
#########


########-----------



########
########------------------------------------------------------------------------------------------------------
######## next for cl #####
########
########


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate seurat4
R
setwd("/zp1/data/share/Human_aging_new")
library(Seurat)

#### Human_RPE = readRDS("Human_RPE_seurat.rds")


Human_AC = readRDS("Human_AC_seurat.rds")
Human_AC_cl = convert_to_genes_NORPE(Human_AC,mcols_df)
saveRDS(Human_AC_cl,file="Human_AC_cl")

Human_HC = readRDS("Human_HC_seurat.rds")
Human_HC_cl = convert_to_genes_NORPE(Human_HC,mcols_df)
saveRDS(Human_HC_cl,file="Human_HC_cl")

Human_Rod = readRDS("Human_Rod_seurat.rds")
Human_Rod_cl = convert_to_genes_NORPE(Human_Rod,mcols_df)
saveRDS(Human_Rod_cl,file="Human_Rod_cl")

Human_Cone = readRDS("Human_Cone_seurat.rds")
Human_Cone_cl = convert_to_genes_NORPE(Human_Cone,mcols_df)
saveRDS(Human_Cone_cl,file="Human_Cone_cl")

Human_BC = readRDS("Human_BC_seurat.rds")
Human_BC_cl = convert_to_genes_NORPE(Human_BC,mcols_df)
saveRDS(Human_BC_cl,file="Human_BC_cl")

Human_Microglia = readRDS("Human_Microglia_seurat.rds")
Human_Microglia_cl = convert_to_genes_NORPE(Human_Microglia,mcols_df)
saveRDS(Human_Microglia_cl,file="Human_Microglia_cl")

Human_MG = readRDS("Human_MG_seurat.rds")
Human_MG_cl = convert_to_genes_NORPE(Human_MG,mcols_df)
saveRDS(Human_MG_cl,file="Human_MG_cl")

Human_RGC = readRDS("Human_RGC_seurat.rds")
Human_RGC_cl = convert_to_genes_NORPE(Human_RGC,mcols_df)
saveRDS(Human_RGC_cl,file="Human_RGC_cl")

#####
##### covert the genes ########
#####

setwd("/zp1/data/share/Human_aging_new")
load(file="mcols_df_2025")


#####
##### seurat_obj = Human_Microglia ###
#####

convert_to_genes_NORPE <- function(seurat_obj,mcols_df){
    ########
    rowN = rownames(seurat_obj)
    k = which(rowN %in% mcols_df$gene_id == T)
    print(length(k) / length(rowN))
    ########
    seurat_obj_cl = seurat_obj[k,]
    ########
    m = match(rownames(seurat_obj_cl),mcols_df$gene_id)
    ########
    rownames(seurat_obj_cl) = mcols_df$gene_name[m]
    ########
    return(seurat_obj_cl)
}


##### Get the overlapped Genes with RPE #######
setwd("/zp1/data/share/Human_aging_new")
Human_AC_cl = readRDS("Human_AC_cl")
Human_RPE_cl = readRDS("Human_RPE_seurat.rds")
total_genes = rownames(Human_AC_cl)
rpe_genes = rownames(Human_RPE_cl)
overlap_genes = total_genes[which(total_genes %in% rpe_genes == T)]

Human_RPE_clcl = Human_RPE_cl[which(rownames(Human_RPE_cl) %in% overlap_genes == T),]
dim(Human_RPE_clcl) #### 36381 ######
saveRDS(Human_RPE_clcl,file="Human_RPE_clcl")
#####
Human_AC_clcl = RemoveDuplicated_seurat(Human_AC_cl)
saveRDS(Human_AC_clcl,file="Human_AC_clcl")

seurat_obj = Human_AC_cl

RemoveDuplicated_seurat <- function(seurat_obj){
    ########
    mat = seurat_obj[["RNA"]]$counts
    ########
    print(dim(mat))
    duplicated_index = which(duplicated(rownames(mat)) == T)
    mat_cl = mat[-duplicated_index,]
    ########
    seurat_new = CreateSeuratObject(mat_cl)
    seurat_new@meta.data = seurat_obj@meta.data
    print(dim(seurat_new))
    ########
    return(seurat_new)
}

#####
#####





Human_BC_cl = readRDS("Human_BC_cl")
Human_BC_clcl = RemoveDuplicated_seurat(Human_BC_cl)
saveRDS(Human_BC_clcl,file="Human_BC_clcl")


Human_HC_cl = readRDS("Human_HC_cl")
Human_HC_clcl = RemoveDuplicated_seurat(Human_HC_cl)
saveRDS(Human_HC_clcl,file="Human_HC_clcl")


Human_Rod_cl = readRDS("Human_Rod_cl")
Human_Rod_clcl = RemoveDuplicated_seurat(Human_Rod_cl)
saveRDS(Human_Rod_clcl,file="Human_Rod_clcl")

Human_MG_cl = readRDS("Human_MG_cl")
Human_MG_clcl = RemoveDuplicated_seurat(Human_MG_cl)
saveRDS(Human_MG_clcl,file="Human_MG_clcl")

Human_RGC_cl = readRDS("Human_RGC_cl")
Human_RGC_clcl = RemoveDuplicated_seurat(Human_RGC_cl)
saveRDS(Human_RGC_clcl,file="Human_RGC_clcl")


Human_Cone_cl = readRDS("Human_Cone_cl")
Human_Cone_clcl = RemoveDuplicated_seurat(Human_Cone_cl)
saveRDS(Human_Cone_clcl,file="Human_Cone_clcl")

Human_Microglia_cl = readRDS("Human_Microglia_cl")
Human_Microglia_clcl = RemoveDuplicated_seurat(Human_Microglia_cl)
saveRDS(Human_Microglia_clcl,file="Human_Microglia_clcl")

library(Seurat)
Human_RPE_cl = readRDS("Human_RPE_seurat.rds")
Human_RPE_clcl = readRDS("Human_RPE_clcl")
saveRDS(Human_RPE_clcl,file="Human_RPE_clcl")

#####
##### 36381 ####
#####

#####
#####

#########-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######### Next we will see their average expression !!!! ######
#########-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#########
######### Seurat_Obj = Human_RPE_process1_M
#########

######### 这个要检查下 #####
#########

Prepare_Avg_Meta_ByDonor_MemOpt <- function(Seurat_Obj){
  library(Matrix)  # 确保 dgCMatrix 等函数可用

  # 1. 原始稀疏 counts 和元数据
  raw_counts    <- Seurat_Obj[["RNA"]]$counts
  meta          <- Seurat_Obj@meta.data

  # 2. 计算每个 cell 的 size factor（CPM 标准化因子）
  lib_sizes     <- colSums(raw_counts)
  size_factors  <- lib_sizes / 1e5

  # 3. age2 分组
  age2_fac      <- factor(meta$age2)
  ages          <- levels(age2_fac)
  n_gene        <- nrow(raw_counts)
  n_age         <- length(ages)

  # 4. 预分配最终矩阵：基因 × age2
  Mat_Avg       <- matrix(
    0,
    nrow = n_gene,
    ncol = n_age,
    dimnames = list(rownames(raw_counts), ages)
  )

  # 5. 两层循环：外层 age2，内层 donor
  for(i in seq_len(n_age)){
    age_i      <- ages[i]
    cells_i    <- which(age2_fac == age_i)
    donors_i   <- unique(meta$donor[cells_i])
    n_donors   <- length(donors_i)

    # 用于累加各 donor 贡献的向量
    col_accum  <- numeric(n_gene)

    for(d in donors_i){
      cells_j   <- cells_i[meta$donor[cells_i] == d]
      sf_j      <- size_factors[cells_j]

      # 只在这个 donor 的列上做归一化 + 求和
      submat    <- raw_counts[, cells_j, drop = FALSE]
      # CPM：每列除以各自的 sf_j
      submat[]  <- submat[] / rep(sf_j, each = nrow(submat))
      # donor 内平均 & 按 donor 数量归一
      col_accum <- col_accum + rowSums(submat) / (length(cells_j) * n_donors)
    }

    # 填入第 i 列
    Mat_Avg[, i] <- col_accum
  }

  # 6. 最后一次 log2(x+1)
  Mat_Avg_Log2 <- log2(Mat_Avg + 1)

  # 7. 构造 metadata
  bulk_meta <- data.frame(
    sample = ages,
    age    = ages,
    stringsAsFactors = FALSE
  )

  list(
    bulk_mat  = Mat_Avg_Log2,
    bulk_meta = bulk_meta
  )
}

Prepare_Avg_Meta_ByDonor_Fast2 <- function(Seurat_Obj){
  library(Matrix)

  # 1. 原始稀疏 counts 和 meta
  raw_counts <- Seurat_Obj[["RNA"]]$counts
  meta       <- Seurat_Obj@meta.data

  # 2. size factors
  sf         <- colSums(raw_counts) / 1e5

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

  list(
    bulk_mat  = mat_log2,   # 现在有 colnames 了
    bulk_meta = bulk_meta
  )
}


########## for RPE cells ###########

Human_RPE_process2_M = Prepare_Avg_Meta_ByDonor_MemOpt(Human_RPE_process1_M)
Human_RPE_process2_F = Prepare_Avg_Meta_ByDonor_MemOpt(Human_RPE_process1_F)
save(Human_RPE_process2_M,file="Human_RPE_process2_M")
save(Human_RPE_process2_F,file="Human_RPE_process2_F")

Human_RPE_DEGs_M = Identify_DEGs(Human_RPE_process2_M)
Human_RPE_DEGs_F = Identify_DEGs(Human_RPE_process2_F)
save(Human_RPE_DEGs_M,file="Human_RPE_DEGs_M")
save(Human_RPE_DEGs_F,file="Human_RPE_DEGs_F")

########### Tested genes: 18018; Significant (p<0.05): 4848
########### Tested genes: 18337; Significant (p<0.05): 3447


########## for AC cells ###########

Human_AC_process2_M = Prepare_Avg_Meta_ByDonor_Fast2(Human_AC_process1_M)
Human_AC_process2_F = Prepare_Avg_Meta_ByDonor_Fast2(Human_AC_process1_F)
save(Human_AC_process2_M,file="Human_AC_process2_M")
save(Human_AC_process2_F,file="Human_AC_process2_F")

Human_AC_DEGs_M = Identify_DEGs(Human_AC_process2_M)
Human_AC_DEGs_F = Identify_DEGs(Human_AC_process2_F)
save(Human_AC_DEGs_M,file="Human_AC_DEGs_M")
save(Human_AC_DEGs_F,file="Human_AC_DEGs_F")


########## for BC cells ###########

Human_BC_process2_M = Prepare_Avg_Meta_ByDonor_Fast2(Human_BC_process1_M)
Human_BC_process2_F = Prepare_Avg_Meta_ByDonor_Fast2(Human_BC_process1_F)
save(Human_BC_process2_M,file="Human_BC_process2_M")
save(Human_BC_process2_F,file="Human_BC_process2_F")

Human_BC_DEGs_M = Identify_DEGs(Human_BC_process2_M)
Human_BC_DEGs_F = Identify_DEGs(Human_BC_process2_F)
save(Human_BC_DEGs_M,file="Human_BC_DEGs_M")
save(Human_BC_DEGs_F,file="Human_BC_DEGs_F")

########## for Microglia cells ###########

Human_Microglia_process2_M = Prepare_Avg_Meta_ByDonor_MemOpt(Human_Microglia_process1_M)
Human_Microglia_process2_F = Prepare_Avg_Meta_ByDonor_MemOpt(Human_Microglia_process1_F)
save(Human_Microglia_process2_M,file="Human_Microglia_process2_M")
save(Human_Microglia_process2_F,file="Human_Microglia_process2_F")

Human_Microglia_DEGs_M = Identify_DEGs(Human_Microglia_process2_M)
Human_Microglia_DEGs_F = Identify_DEGs(Human_Microglia_process2_F)
save(Human_Microglia_DEGs_M,file="Human_Microglia_DEGs_M")
save(Human_Microglia_DEGs_F,file="Human_Microglia_DEGs_F")


########## for HC cells ###########

Human_HC_process2_M = Prepare_Avg_Meta_ByDonor_MemOpt(Human_HC_process1_M)
Human_HC_process2_F = Prepare_Avg_Meta_ByDonor_MemOpt(Human_HC_process1_F)
save(Human_HC_process2_M,file="Human_HC_process2_M")
save(Human_HC_process2_F,file="Human_HC_process2_F")

Human_HC_DEGs_M = Identify_DEGs(Human_HC_process2_M)
Human_HC_DEGs_F = Identify_DEGs(Human_HC_process2_F)
save(Human_HC_DEGs_M,file="Human_HC_DEGs_M")
save(Human_HC_DEGs_F,file="Human_HC_DEGs_F")


########## for MG cells ###########

Human_MG_process2_M = Prepare_Avg_Meta_ByDonor_MemOpt(Human_MG_process1_M)
Human_MG_process2_F = Prepare_Avg_Meta_ByDonor_MemOpt(Human_MG_process1_F)
save(Human_MG_process2_M,file="Human_MG_process2_M")
save(Human_MG_process2_F,file="Human_MG_process2_F")

Human_MG_DEGs_M = Identify_DEGs(Human_MG_process2_M)
Human_MG_DEGs_F = Identify_DEGs(Human_MG_process2_F)
save(Human_MG_DEGs_M,file="Human_MG_DEGs_M")
save(Human_MG_DEGs_F,file="Human_MG_DEGs_F")


########## for Cone cells ###########

Human_Cone_process2_M = Prepare_Avg_Meta_ByDonor_Fast2(Human_Cone_process1_M)
Human_Cone_process2_F = Prepare_Avg_Meta_ByDonor_Fast2(Human_Cone_process1_F)
save(Human_Cone_process2_M,file="Human_Cone_process2_M")
save(Human_Cone_process2_F,file="Human_Cone_process2_F")

Human_Cone_DEGs_M = Identify_DEGs(Human_Cone_process2_M)
Human_Cone_DEGs_F = Identify_DEGs(Human_Cone_process2_F)
save(Human_Cone_DEGs_M,file="Human_Cone_DEGs_M")
save(Human_Cone_DEGs_F,file="Human_Cone_DEGs_F")


########## for Rod cells ###########

Human_Rod_process2_M = Prepare_Avg_Meta_ByDonor_Fast2(Human_Rod_process1_M)
Human_Rod_process2_F = Prepare_Avg_Meta_ByDonor_Fast2(Human_Rod_process1_F)
save(Human_Rod_process2_M,file="Human_Rod_process2_M")
save(Human_Rod_process2_F,file="Human_Rod_process2_F")

Human_Rod_DEGs_M = Identify_DEGs(Human_Rod_process2_M)
Human_Rod_DEGs_F = Identify_DEGs(Human_Rod_process2_F)
save(Human_Rod_DEGs_M,file="Human_Rod_DEGs_M")
save(Human_Rod_DEGs_F,file="Human_Rod_DEGs_F")


########## for RGC cells ###########


Human_RGC_process2_M = Prepare_Avg_Meta_ByDonor_Fast2(Human_RGC_process1_M)
Human_RGC_process2_F = Prepare_Avg_Meta_ByDonor_Fast2(Human_RGC_process1_F)
save(Human_RGC_process2_M,file="Human_RGC_process2_M")
save(Human_RGC_process2_F,file="Human_RGC_process2_F")

Human_RGC_DEGs_M = Identify_DEGs(Human_RGC_process2_M)
Human_RGC_DEGs_F = Identify_DEGs(Human_RGC_process2_F)
save(Human_RGC_DEGs_M,file="Human_RGC_DEGs_M")
save(Human_RGC_DEGs_F,file="Human_RGC_DEGs_F")

#######  
####### Next we will plot the Heatmaps !!! ########
####### 




library(mgcv)


# 1. 相关计算函数（同你原来的一样，只改了命名）
Cor_FUN <- function(x, age) {
  vec <- as.numeric(x)
  cor(vec, age)
}

# 2. 平滑函数：加上对 s.table 空矩阵的保护，并用 tryCatch 捕获可能的 fit 错误
Smooth_FUN <- function(x, age) {
  vec <- as.numeric(x)
  df  <- data.frame(age = age, value = vec)
  
  if (length(unique(vec)) == 1) {
    # 常数序列，直接返回 NA p_value
    return(data.frame(age = age,
                      value = vec,
                      predict = vec,
                      p_value = NA,
                      row.names = NULL))
  }
  
  df <- df[order(df$age), ]
  rownames(df) <- NULL
  
  # 拦截模型拟合可能出错的情况
  fit <- tryCatch(
    gam(value ~ s(age, k = 5), data = df, select = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    df$predict <- NA
    df$p_value <- NA
    return(df)
  }
  
  s_tab <- summary(fit)$s.table
  # 如果没有平滑项，就给 NA
  if (is.null(s_tab) || nrow(s_tab) < 1) {
    pval <- NA
  } else {
    # 第 4 列通常是 p-value
    pval <- s_tab[1, 4]
  }
  
  df$predict <- predict(fit, newdata = df, type = "response")
  df$p_value <- pval
  df
}

# 3. 主函数：用 lapply 保证 Smooth_FUN 返回列表，不用动态扩容

## list_process2 = Human_RPE_process2_M

Identify_DEGs <- function(list_process2) {
  bulk_mat  <- list_process2$bulk_mat
  bulk_meta <- list_process2$bulk_meta
  bulk_meta$sample = as.character(as.numeric(bulk_meta$sample))
  bulk_meta$age = as.numeric(bulk_meta$age)
  # 强校验
  stopifnot(identical(colnames(bulk_mat), bulk_meta$sample))
  
  # 基因过滤
  keep <- rowMeans(bulk_mat) > 0.1
  bulk_mat <- bulk_mat[keep, , drop = FALSE]
  
  # 相关性
  corr_vals <- apply(bulk_mat, 1, Cor_FUN, age = bulk_meta$age)
  
  # 平滑拟合
  gene_list   <- rownames(bulk_mat)
  smooth_list <- lapply(seq_along(gene_list), function(i) {
    Smooth_FUN(bulk_mat[i, ], bulk_meta$age)
  })
  names(smooth_list) <- gene_list
  
  # 提取 p-values
  pvals <- sapply(smooth_list, function(df) df$p_value[1])
  
  # 汇总
  res <- data.frame(
    Gene  = gene_list,
    Cor   = corr_vals,
    pval  = pvals,
    stringsAsFactors = FALSE
  )
  degs <- subset(res, pval < 0.05)
  
  message("Tested genes: ", nrow(res),
          "; Significant (p<0.05): ", nrow(degs))
  degs
}



##########
########## Plot the Heatmaps !!!!!!!! ###############
##########

setwd("/zp1/data/share/Human_aging_new")

##########

load("Human_MG_process2_M")
load("Human_MG_process2_F")
load("Human_Rod_process2_M")
load("Human_Rod_process2_F")
load("Human_Cone_process2_M")
load("Human_Cone_process2_F")
load("Human_BC_process2_M")
load("Human_BC_process2_F")
load("Human_RGC_process2_M")
load("Human_RGC_process2_F")
load("Human_AC_process2_M")
load("Human_AC_process2_F")
load("Human_HC_process2_M")
load("Human_HC_process2_F")
load("Human_Microglia_process2_M")
load("Human_Microglia_process2_F")
load("Human_RPE_process2_M")
load("Human_RPE_process2_F")
load("Human_MG_DEGs_M")
load("Human_MG_DEGs_F")
load("Human_Rod_DEGs_M")
load("Human_Rod_DEGs_F")
load("Human_Cone_DEGs_M")
load("Human_Cone_DEGs_F")
load("Human_BC_DEGs_M")
load("Human_BC_DEGs_F")
load("Human_RGC_DEGs_M")
load("Human_RGC_DEGs_F")
load("Human_AC_DEGs_M")
load("Human_AC_DEGs_F")
load("Human_HC_DEGs_M")
load("Human_HC_DEGs_F")
load("Human_Microglia_DEGs_M")
load("Human_Microglia_DEGs_F")
load("Human_RPE_DEGs_M")
load("Human_RPE_DEGs_F")

#####
#####
#####


DEGs_List1 = list(Human_MG_process2_M[[1]],Human_MG_process2_F[[1]],Human_Rod_process2_M[[1]],Human_Rod_process2_F[[1]])
DEGs_List2 = list(Human_RGC_process2_M[[1]],Human_RGC_process2_F[[1]],Human_AC_process2_M[[1]],Human_AC_process2_F[[1]])
DEGs_List3 = list(Human_HC_process2_M[[1]],Human_HC_process2_F[[1]],Human_BC_process2_M[[1]],Human_BC_process2_F[[1]])
DEGs_List4 = list(Human_Cone_process2_M[[1]],Human_Cone_process2_F[[1]],Human_Microglia_process2_M[[1]],Human_Microglia_process2_F[[1]])
DEGs_List5 = list(Human_RPE_process2_M[[1]],Human_RPE_process2_F[[1]])

DEGs_List = c(DEGs_List1,DEGs_List2,DEGs_List3,DEGs_List4,DEGs_List5)
names(DEGs_List) = c("MG:M","MG:F","Rod:M","Rod:F","RGC:M","RGC:F","AC:M","AC:F","HC:M","HC:F","BC:M","BC:F","Cone:M","Cone:F","Microglia:M","Microglia:F","RPE:M","RPE:F")

#####
#####
#####
length(DEGs_List)
length(Avgs_List)


Avgs_List1 = list(Human_MG_DEGs_M,Human_MG_DEGs_F,Human_Rod_DEGs_M,Human_Rod_DEGs_F)
Avgs_List2 = list(Human_RGC_DEGs_M,Human_RGC_DEGs_F,Human_AC_DEGs_M,Human_AC_DEGs_F)
Avgs_List3 = list(Human_HC_DEGs_M,Human_HC_DEGs_F,Human_BC_DEGs_M,Human_BC_DEGs_F)
Avgs_List4 = list(Human_Cone_DEGs_M,Human_Cone_DEGs_F,Human_Microglia_DEGs_M,Human_Microglia_DEGs_F)
Avgs_List5 = list(Human_RPE_DEGs_M,Human_RPE_DEGs_F)

Avgs_List = c(Avgs_List1,Avgs_List2,Avgs_List3,Avgs_List4,Avgs_List5)
names(Avgs_List) = c("MG:M","MG:F","Rod:M","Rod:F","RGC:M","RGC:F","AC:M","AC:F","HC:M","HC:F","BC:M","BC:F","Cone:M","Cone:F","Microglia:M","Microglia:F","RPE:M","RPE:F")
col_index = c("15","25","37.5","50","60","70","80","90")

Merge_DEGs_genes_and_Avg_expression_Human <- function(DEGs_List,Avgs_List,col_index){
    #########
    #########
    for(i in 1:length(DEGs_List)){
        ####
        tmp = DEGs_List[[i]]
        rownames(tmp) = paste(names(DEGs_List)[i],rownames(tmp),sep=":")
        DEGs_List[[i]] = tmp
    }
    #########
    all_Matrix = do.call("rbind",DEGs_List)
    #########
    for(i in 1:length(Avgs_List)){
        ####
        tmp = Avgs_List[[i]]
        tmp$Gene = paste(names(Avgs_List)[i],rownames(tmp),sep=":")
        Avgs_List[[i]] = tmp
    }
    #########
    all_Avg = do.call("rbind",Avgs_List)
    all_Avg = all_Avg[which(all_Avg$pval < 0.05),]
    ######### all_Avg[grep("CLU",all_Avg$Gene),] #########
    ######### all_Avg[grep("APOE",all_Avg$Gene),]
    ######### all_Avg[grep("STAT",all_Avg$Gene),]
    all_Matrix_cl = all_Matrix[which(rownames(all_Matrix) %in% all_Avg$Gene == T),]
    #########
    all_Matrix_cl = all_Matrix_cl[,col_index]
    #########
    all_Matrix_cl_s = t(apply(all_Matrix_cl,1,scale))
    rownames(all_Matrix_cl_s) = rownames(all_Matrix_cl)
    colnames(all_Matrix_cl_s) = colnames(all_Matrix_cl)
    #########
    return(all_Matrix_cl_s)
}

Human_Res_Plot = Merge_DEGs_genes_and_Avg_expression_Human(DEGs_List,Avgs_List,col_index)


#######
k_up = which(Human_Res_Plot > 5)
k_down = which(Human_Res_Plot < -5)

Human_Res_Plot[k_up] = 5
Human_Res_Plot[k_down] = -5

kc = kmeans(Human_Res_Plot, 18, iter.max = 100)
kc_dat = data.frame(genes = names(kc$cluster),cluster=as.numeric(kc$cluster))
#### kc_dat[which(kc_dat$genes == 'MG:M:FRZB'),]

setwd("/zp1/data/share/Human_aging_new")
save(kc_dat,file="Human_Aging_kc_dat_2025Mar")

#####

######----------
library('ComplexHeatmap')
library('circlize')

celltypes <- c("MG","AC","Cone","RGC","HC","Rod","BC","RPE","Microglia","Astrocyte")
cols = c("#D11536","#F6BA00","#9EA220","#AAA9A9","#EF9000","#026AB1","#804537","#936DAD","#61BFB9","#EC6F64")

celltype = sapply(strsplit(rownames(Human_Res_Plot),split=":"),function(x) x[[1]])
sex = sapply(strsplit(rownames(Human_Res_Plot),split=":"),function(x) x[[2]])

anno_df = data.frame(
  celltype = celltype,
  sex = sex
)

col_list <- list(
  celltype = c("MG" = "#D11536", "RGC" = "#AAA9A9", "AC" = "#F6BA00","HC"="#EF9000","BC"="#804537","Rod"="#026AB1","Cone"="#9EA220","RPE"="#936DAD","Microglia"="#EC6F64"),
  sex = c("M"="#026AB1","F"="#D11536")
)

row_anno <- rowAnnotation(df = anno_df, col = col_list)

######----------
kc_dat_cl = kc_dat

row_sp = as.factor(kc_dat_cl$cluster)
colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")
col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))

### rownames(Human_Res_Plot)[grep("FRZB",rownames(Human_Res_Plot))]

labels = c("AC:F:STAT1","HC:F:STAT1","BC:F:STAT1","Cone:M:STAT1","Cone:F:STAT1","MG:M:FRZB")
at = match(labels,rownames(Human_Res_Plot))

setwd("/zp1/data/share/Human_aging_new")
png('Human_DEGs.png',height=7000,width=4500,res=72*12)
Heatmap(Human_Res_Plot, name = "XX",right_annotation = row_anno,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

#######
#######
#######

#### keep the cluster ######
rm = which(row_sp %in% c(1,5,7,10,11,16,18,8) == F)
row_sp_cl = as.character(row_sp[rm])
Human_Res_Plot_cl = Human_Res_Plot[rm,]

row_anno <- rowAnnotation(df = anno_df[rm,], col = col_list)

rownames(Human_Res_Plot_cl)[grep("STAT",rownames(Human_Res_Plot_cl))]
rownames(Human_Res_Plot_cl)[grep("APOE",rownames(Human_Res_Plot_cl))]

labels = c("RGC:F:STAT3","RPE:F:STAT3")
at = match(labels,rownames(Human_Res_Plot_cl))

setwd("/zp1/data/share/Human_aging_new")
png('Human_DEGs.png',height=7000,width=4500,res=72*12)
Heatmap(Human_Res_Plot_cl, name = "XX",right_annotation = row_anno,border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp_cl) + 
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

######
######


library(reshape2)

all_Matrix_cl_s_Plot = melt(Human_Res_Plot_cl)
all_Matrix_cl_s_Plot$class = kc_dat$cluster[match(all_Matrix_cl_s_Plot$Var1,kc_dat$genes)]

levels_tag = levels(as.factor(all_Matrix_cl_s_Plot$Var2))
all_Matrix_cl_s_Plot$Var2_1 = as.numeric(match(all_Matrix_cl_s_Plot$Var2,levels_tag))

summary(all_Matrix_cl_s_Plot$value)

all_Matrix_cl_s_Plot$class = factor(all_Matrix_cl_s_Plot$class,levels=c(12,13,14,15,17,2,3,4,6,9))

library(ggplot2)

ggplot(all_Matrix_cl_s_Plot, aes(x = Var2_1, y = value)) + 
  geom_point(size=0, position = position_jitter(width = 0.2, height = 0.2),alpha=0.01) + facet_wrap(~ class, ncol = 1) +
  stat_summary(fun = mean, geom = "line", color = "red", size = 1) +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

ggsave("Human_trend.png",height=12,width=1.5)


##########
########## 下面我们每个细胞系做一下GO #######
##########


Upclusters = c(12,3,6,9)
Downclusters = c(13,14,15,17)

kc_dat

Human_YO_Genes = Get_GOandKEGG_Human_UPDOWN(kc_dat,Upclusters,Downclusters)

save(Human_YO_Genes,file="Human_YO_Genes")

#######

Human_YO_GO = list()
for(i in 1:length(Human_YO_Genes)){
    print(i)
    tmp_res = Get_GOandKEGG_Human(Human_YO_Genes[[i]])
    ######
    Human_YO_GO <- c(Human_YO_GO,list(tmp_res))
}
names(Human_YO_GO) = names(Human_YO_Genes)

library(openxlsx)
write.xlsx(Human_YO_GO, file = "Human_GO_KEGG_by_YoungOld_Apr2025.xlsx")

########
########
########










