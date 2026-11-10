"""
SATURN: Mouse-Human Cross-Species Integration
完整流程脚本，分 4 个阶段：
  Step 1: 准备 Mouse h5ad
  Step 2: 准备 Human h5ad
  Step 3: 生成 protein embeddings (ESM2)
  Step 4: 生成 SATURN 输入文件 (CSV)
训练和下游分析在 run_saturn_train.sh 和 saturn_downstream.py 中
"""

import os
import scanpy as sc
import anndata as ad
import numpy as np
import scipy.sparse as sp
import pandas as pd

WORK_DIR = "/projects/hmz-aging/Zebrafish_Mouse_Human_inte/"
SATURN_DIR = "/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN"  # clone SATURN 到这里
SATURN_DATA_DIR = os.path.join(WORK_DIR, "saturn_data")
os.makedirs(SATURN_DATA_DIR, exist_ok=True)
os.chdir(WORK_DIR)

# ============================================================
# Step 1: 准备 Mouse h5ad
# SATURN 需要 raw counts 在 .X 里面，obs 里有 celltype 列
# ============================================================
print("=" * 60)
print("Step 1: 准备 Mouse h5ad")
print("=" * 60)

adata_ms = sc.read_h5ad("mouse_adata_decontx.h5ad")
print(f"Mouse 原始: {adata_ms.shape}")
print(f"obs columns: {list(adata_ms.obs.columns)}")

if "decontX_counts" in adata_ms.layers:
    adata_ms.X = adata_ms.layers["decontX_counts"].copy()
    print("使用 decontX_counts 作为 .X")
elif "counts" in adata_ms.layers:
    adata_ms.X = adata_ms.layers["counts"].copy()
    print("使用 counts 作为 .X")
else:
    print("使用默认 .X")

if sp.issparse(adata_ms.X):
    adata_ms.X = sp.csr_matrix(adata_ms.X, dtype=np.float32)
else:
    adata_ms.X = sp.csr_matrix(np.asarray(adata_ms.X, dtype=np.float32))

# SATURN 要求 celltype 标签列 -- 检查列名
print(f"\nMouse obs 列: {list(adata_ms.obs.columns)}")
if "celltype" in adata_ms.obs.columns:
    print(f"Mouse celltype 分布:\n{adata_ms.obs['celltype'].value_counts()}")
else:
    print("WARNING: 没有 'celltype' 列，请手动指定标签列名")

# 每个 celltype 抽样 2000 cells
sampled_idx = []
for ct, group in adata_ms.obs.groupby("celltype"):
    n = min(len(group), 2000)
    sampled_idx.append(group.sample(n=n, random_state=42).index)
sampled_idx = np.concatenate(sampled_idx)
adata_ms_sub = adata_ms[sampled_idx].copy()
print(f"\nMouse 抽样后: {adata_ms_sub.shape}")
print(adata_ms_sub.obs["celltype"].value_counts())

ms_out = os.path.join(SATURN_DATA_DIR, "mouse_saturn.h5ad")
adata_ms_sub.write(ms_out)
print(f"保存: {ms_out}")

# ============================================================
# Step 2: 准备 Human h5ad
# ============================================================
print("\n" + "=" * 60)
print("Step 2: 准备 Human h5ad")
print("=" * 60)

h5ad_dir = "/projects/hmz-aging/Human_New_Datasets/celltype_h5ad"
adata_list = []
for f in sorted(os.listdir(h5ad_dir)):
    if f.endswith("_clean4.h5ad"):
        ct = f.replace("_clean4.h5ad", "")
        path = os.path.join(h5ad_dir, f)
        a = sc.read_h5ad(path)
        a.obs["celltype"] = ct
        print(f"  {ct}: {a.shape}")
        adata_list.append(a)

adata_human = ad.concat(adata_list, merge="same")
adata_human.obs_names_make_unique()
print(f"\nHuman 合并后: {adata_human.shape}")
print(adata_human.obs["celltype"].value_counts())

if sp.issparse(adata_human.X):
    adata_human.X = sp.csr_matrix(adata_human.X, dtype=np.float32)
else:
    adata_human.X = sp.csr_matrix(np.asarray(adata_human.X, dtype=np.float32))

# 每个 celltype 抽样 2000 cells
sampled_idx_hu = []
for ct, group in adata_human.obs.groupby("celltype"):
    n = min(len(group), 2000)
    sampled_idx_hu.append(group.sample(n=n, random_state=42).index)
sampled_idx_hu = np.concatenate(sampled_idx_hu)
adata_hu_sub = adata_human[sampled_idx_hu].copy()
print(f"\nHuman 抽样后: {adata_hu_sub.shape}")
print(adata_hu_sub.obs["celltype"].value_counts())

hu_out = os.path.join(SATURN_DATA_DIR, "human_saturn.h5ad")
adata_hu_sub.write(hu_out)
print(f"保存: {hu_out}")

# ============================================================
# Step 3: Protein Embeddings
# SATURN 需要每个物种的 protein embedding (.pt 文件)
# 可以下载预计算的 或者自己用 ESM2 生成
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Protein Embeddings")
print("=" * 60)

EMBEDDING_DIR = "/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN/protein_embeddings/protein_embeddings_export/ESM2"

mouse_emb_path = os.path.join(EMBEDDING_DIR, "mouse_embedding.torch")
human_emb_path = os.path.join(EMBEDDING_DIR, "human_embedding.torch")

if os.path.exists(mouse_emb_path) and os.path.exists(human_emb_path):
    print("Protein embeddings 已存在")
else:
    print(f"WARNING: embedding 文件不存在，请先下载并解压到 {EMBEDDING_DIR}")
    print(f"  需要: {mouse_emb_path}")
    print(f"  需要: {human_emb_path}")

# 检查 embedding 文件中的基因覆盖率
import torch
for tag, emb_path, adata_check in [("Mouse", mouse_emb_path, adata_ms_sub),
                                     ("Human", human_emb_path, adata_hu_sub)]:
    if os.path.exists(emb_path):
        emb = torch.load(emb_path, map_location="cpu")
        emb_genes = set(k.lower() for k in emb.keys())
        adata_genes = set(g.lower() for g in adata_check.var_names)
        overlap = emb_genes & adata_genes
        print(f"{tag}: embedding 基因数={len(emb_genes)}, "
              f"adata 基因数={len(adata_genes)}, "
              f"重叠={len(overlap)} ({len(overlap)/len(adata_genes)*100:.1f}%)")
    else:
        print(f"{tag}: embedding 文件不存在 ({emb_path})")

# ============================================================
# Step 4: 生成 SATURN 输入 CSV 和 cell type mapping CSV
# ============================================================
print("\n" + "=" * 60)
print("Step 4: 生成输入文件")
print("=" * 60)

# 4a. Input data CSV (path, species, embedding_path)
run_csv = pd.DataFrame({
    "path": [ms_out, hu_out],
    "species": ["mouse", "human"],
    "embedding_path": [mouse_emb_path, human_emb_path]
})
run_csv_path = os.path.join(SATURN_DATA_DIR, "mouse_human_run.csv")
run_csv.to_csv(run_csv_path, index=False)
print(f"Input CSV 保存: {run_csv_path}")
print(run_csv)

# 4b. Cell type mapping CSV (用于训练时评分)
ms_cts = sorted(adata_ms_sub.obs["celltype"].unique())
hu_cts = sorted(adata_hu_sub.obs["celltype"].unique())
print(f"\nMouse celltypes: {ms_cts}")
print(f"Human celltypes: {hu_cts}")

# 自动匹配同名的 celltype
all_cts = sorted(set(ms_cts) | set(hu_cts))
ct_map_rows = []
for ct in all_cts:
    row = {
        "mouse_cell_type": ct if ct in ms_cts else np.nan,
        "human_cell_type": ct if ct in hu_cts else np.nan
    }
    ct_map_rows.append(row)
ct_map = pd.DataFrame(ct_map_rows)
ct_map_path = os.path.join(SATURN_DATA_DIR, "mouse_human_ct_map.csv")
ct_map.to_csv(ct_map_path, index=False)
print(f"\nCell type mapping CSV 保存: {ct_map_path}")
print(ct_map)

print("\n" + "=" * 60)
print("数据准备完成！")
print(f"所有文件在: {SATURN_DATA_DIR}")
print("=" * 60)
print(f"""
下一步:
  1. 确保 protein embeddings 存在:
     {mouse_emb_path}
     {human_emb_path}
  
  2. 运行 SATURN 训练:
     conda activate saturn_env
     bash run_saturn_train.sh
""")
