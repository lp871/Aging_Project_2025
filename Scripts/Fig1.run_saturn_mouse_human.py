"""
SATURN: Mouse-Human Cross-Species Integration
Full pipeline script in 4 stages:
  Step 1: Prepare mouse h5ad
  Step 2: Prepare human h5ad
  Step 3: Generate protein embeddings (ESM2)
  Step 4: Generate SATURN input files (CSV)
Training and downstream analysis are in run_saturn_train.sh and saturn_downstream.py
"""

import os
import scanpy as sc
import anndata as ad
import numpy as np
import scipy.sparse as sp
import pandas as pd

WORK_DIR = "/projects/hmz-aging/Zebrafish_Mouse_Human_inte/"
SATURN_DIR = "/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN"  # clone SATURN here
SATURN_DATA_DIR = os.path.join(WORK_DIR, "saturn_data")
os.makedirs(SATURN_DATA_DIR, exist_ok=True)
os.chdir(WORK_DIR)

# ============================================================
# Step 1: Prepare mouse h5ad
# SATURN requires raw counts in .X and a celltype column in obs
# ============================================================
print("=" * 60)
print("Step 1: Prepare mouse h5ad")
print("=" * 60)

adata_ms = sc.read_h5ad("mouse_adata_decontx.h5ad")
print(f"Mouse raw: {adata_ms.shape}")
print(f"obs columns: {list(adata_ms.obs.columns)}")

if "decontX_counts" in adata_ms.layers:
    adata_ms.X = adata_ms.layers["decontX_counts"].copy()
    print("Using decontX_counts as .X")
elif "counts" in adata_ms.layers:
    adata_ms.X = adata_ms.layers["counts"].copy()
    print("Using counts as .X")
else:
    print("Using default .X")

if sp.issparse(adata_ms.X):
    adata_ms.X = sp.csr_matrix(adata_ms.X, dtype=np.float32)
else:
    adata_ms.X = sp.csr_matrix(np.asarray(adata_ms.X, dtype=np.float32))

# SATURN requires a celltype label column -- check column names
print(f"\nMouse obs columns: {list(adata_ms.obs.columns)}")
if "celltype" in adata_ms.obs.columns:
    print(f"Mouse celltype distribution:\n{adata_ms.obs['celltype'].value_counts()}")
else:
    print("WARNING: no 'celltype' column; please specify the label column manually")

# Sample 2000 cells per celltype
sampled_idx = []
for ct, group in adata_ms.obs.groupby("celltype"):
    n = min(len(group), 2000)
    sampled_idx.append(group.sample(n=n, random_state=42).index)
sampled_idx = np.concatenate(sampled_idx)
adata_ms_sub = adata_ms[sampled_idx].copy()
print(f"\nMouse after subsampling: {adata_ms_sub.shape}")
print(adata_ms_sub.obs["celltype"].value_counts())

ms_out = os.path.join(SATURN_DATA_DIR, "mouse_saturn.h5ad")
adata_ms_sub.write(ms_out)
print(f"Saved: {ms_out}")

# ============================================================
# Step 2: Prepare human h5ad
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Prepare human h5ad")
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
print(f"\nHuman after merging: {adata_human.shape}")
print(adata_human.obs["celltype"].value_counts())

if sp.issparse(adata_human.X):
    adata_human.X = sp.csr_matrix(adata_human.X, dtype=np.float32)
else:
    adata_human.X = sp.csr_matrix(np.asarray(adata_human.X, dtype=np.float32))

# Sample 2000 cells per celltype
sampled_idx_hu = []
for ct, group in adata_human.obs.groupby("celltype"):
    n = min(len(group), 2000)
    sampled_idx_hu.append(group.sample(n=n, random_state=42).index)
sampled_idx_hu = np.concatenate(sampled_idx_hu)
adata_hu_sub = adata_human[sampled_idx_hu].copy()
print(f"\nHuman after subsampling: {adata_hu_sub.shape}")
print(adata_hu_sub.obs["celltype"].value_counts())

hu_out = os.path.join(SATURN_DATA_DIR, "human_saturn.h5ad")
adata_hu_sub.write(hu_out)
print(f"Saved: {hu_out}")

# ============================================================
# Step 3: Protein Embeddings
# SATURN requires a protein embedding (.pt file) for each species
# Download precomputed embeddings or generate them with ESM2
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Protein Embeddings")
print("=" * 60)

EMBEDDING_DIR = "/projects/hmz-aging/Zebrafish_Mouse_Human_inte/SATURN/protein_embeddings/protein_embeddings_export/ESM2"

mouse_emb_path = os.path.join(EMBEDDING_DIR, "mouse_embedding.torch")
human_emb_path = os.path.join(EMBEDDING_DIR, "human_embedding.torch")

if os.path.exists(mouse_emb_path) and os.path.exists(human_emb_path):
    print("Protein embeddings already exist")
else:
    print(f"WARNING: embedding files not found; please download and extract them to {EMBEDDING_DIR}")
    print(f"  Required: {mouse_emb_path}")
    print(f"  Required: {human_emb_path}")

# Check gene coverage of the embedding files
import torch
for tag, emb_path, adata_check in [("Mouse", mouse_emb_path, adata_ms_sub),
                                     ("Human", human_emb_path, adata_hu_sub)]:
    if os.path.exists(emb_path):
        emb = torch.load(emb_path, map_location="cpu")
        emb_genes = set(k.lower() for k in emb.keys())
        adata_genes = set(g.lower() for g in adata_check.var_names)
        overlap = emb_genes & adata_genes
        print(f"{tag}: embedding genes={len(emb_genes)}, "
              f"adata genes={len(adata_genes)}, "
              f"overlap={len(overlap)} ({len(overlap)/len(adata_genes)*100:.1f}%)")
    else:
        print(f"{tag}: embedding file not found ({emb_path})")

# ============================================================
# Step 4: Generate SATURN input CSV and cell type mapping CSV
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Generate input files")
print("=" * 60)

# 4a. Input data CSV (path, species, embedding_path)
run_csv = pd.DataFrame({
    "path": [ms_out, hu_out],
    "species": ["mouse", "human"],
    "embedding_path": [mouse_emb_path, human_emb_path]
})
run_csv_path = os.path.join(SATURN_DATA_DIR, "mouse_human_run.csv")
run_csv.to_csv(run_csv_path, index=False)
print(f"Input CSV saved: {run_csv_path}")
print(run_csv)

# 4b. Cell type mapping CSV (used for scoring during training)
ms_cts = sorted(adata_ms_sub.obs["celltype"].unique())
hu_cts = sorted(adata_hu_sub.obs["celltype"].unique())
print(f"\nMouse celltypes: {ms_cts}")
print(f"Human celltypes: {hu_cts}")

# Automatically match celltypes with identical names
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
print(f"\nCell type mapping CSV saved: {ct_map_path}")
print(ct_map)

print("\n" + "=" * 60)
print("Data preparation complete!")
print(f"All files are in: {SATURN_DATA_DIR}")
print("=" * 60)
print(f"""
Next steps:
  1. Make sure protein embeddings exist:
     {mouse_emb_path}
     {human_emb_path}
  
  2. Run SATURN training:
     conda activate saturn_env
     bash run_saturn_train.sh
""")
