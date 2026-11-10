"""
SATURN 下游分析 - Mouse x Human
训练完成后运行，生成 UMAP 和评分
"""

import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

WORK_DIR = "/projects/hmz-aging/Zebrafish_Mouse_Human_inte/"
SATURN_DATA = os.path.join(WORK_DIR, "saturn_data")
os.chdir(WORK_DIR)

# ============================================================
# Step 1: 找到 SATURN 输出的 h5ad 文件
# ============================================================
result_dir = os.path.join(SATURN_DATA, "saturn_results")
if not os.path.isdir(result_dir):
    result_dir = SATURN_DATA

h5ad_files = [f for f in os.listdir(result_dir) if f.endswith(".h5ad") and "pretrain" not in f and "ep_" not in f]
print("SATURN 输出文件:")
for f in sorted(h5ad_files):
    print(f"  {f}")

if not h5ad_files:
    print("没找到 SATURN 输出 h5ad，请检查训练是否完成")
    print(f"查找路径: {result_dir}")
    import sys; sys.exit(1)

saturn_h5ad = os.path.join(result_dir, sorted(h5ad_files)[-1])
print(f"\n使用: {saturn_h5ad}")

# ============================================================
# Step 2: 加载 SATURN 整合后的 adata
# ============================================================
adata = sc.read_h5ad(saturn_h5ad)
print(f"\nSATURN output: {adata.shape}")
print(f"obs columns: {list(adata.obs.columns)}")

# SATURN 输出的 .X 是 SATURN embeddings (latent space)
# .obs 中有:
#   - 'labels2': 原始 celltype 标签
#   - 'labels': species + celltype
#   - 'species': 物种名
#   - 'ref_labels': 参考标签
# .obsm 中有:
#   - 'macrogenes': macrogene values

print(f"\nSpecies 分布: {adata.obs['species'].value_counts().to_dict()}")
if "labels2" in adata.obs.columns:
    print(f"Cell types: {sorted(adata.obs['labels2'].unique())}")
elif "labels" in adata.obs.columns:
    print(f"Labels: {sorted(adata.obs['labels'].unique())}")

# ============================================================
# Step 3: UMAP 可视化
# ============================================================
print("\n计算 neighbors 和 UMAP...")
sc.pp.neighbors(adata, use_rep="X", n_neighbors=15)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)
sc.tl.leiden(adata, resolution=0.5)

color_map = {
    "MG": "#C44E52", "AC": "#E5C953", "Cone": "#BCBD22",
    "RGC": "#B0B0B0", "HC": "#E5933A", "Rod": "#4878A8",
    "BC": "#996633", "RPE": "#B89BC7", "Microglia": "#8ECFC9",
}

label_col = "labels2" if "labels2" in adata.obs.columns else "celltype"

# 3a. UMAP colored by celltype
fig, ax = plt.subplots(figsize=(7, 6))
celltypes = sorted(adata.obs[label_col].unique())
for ct in celltypes:
    mask = adata.obs[label_col] == ct
    c = color_map.get(ct, None)
    ax.scatter(
        adata.obsm["X_umap"][mask, 0],
        adata.obsm["X_umap"][mask, 1],
        c=c, s=8, alpha=0.6, label=ct, rasterized=True
    )
ax.legend(fontsize=9, frameon=False, loc="upper left", markerscale=2)
ax.set_title("SATURN: Mouse + Human (colored by celltype)", fontsize=12)
ax.axis("off")
plt.tight_layout()
plt.savefig("saturn_umap_celltype.png", dpi=300, bbox_inches="tight")
plt.close()
print("saved: saturn_umap_celltype.png")

# 3b. UMAP colored by species
fig, ax = plt.subplots(figsize=(7, 6))
species_colors = {"mouse": "#E74C3C", "human": "#2ECC71"}
for sp_name in ["mouse", "human"]:
    mask = adata.obs["species"] == sp_name
    ax.scatter(
        adata.obsm["X_umap"][mask, 0],
        adata.obsm["X_umap"][mask, 1],
        c=species_colors[sp_name], s=8, alpha=0.5,
        label=sp_name.capitalize(), rasterized=True
    )
ax.legend(fontsize=11, frameon=False, loc="upper left", markerscale=3)
ax.set_title("SATURN: Mouse + Human (colored by species)", fontsize=12)
ax.axis("off")
plt.tight_layout()
plt.savefig("saturn_umap_species.png", dpi=300, bbox_inches="tight")
plt.close()
print("saved: saturn_umap_species.png")

# 3c. UMAP split by species side-by-side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for i, sp_name in enumerate(["mouse", "human"]):
    ax = axes[i]
    mask_sp = adata.obs["species"] == sp_name
    adata_sp = adata[mask_sp]
    for ct in celltypes:
        mask_ct = adata_sp.obs[label_col] == ct
        c = color_map.get(ct, None)
        ax.scatter(
            adata_sp.obsm["X_umap"][mask_ct, 0],
            adata_sp.obsm["X_umap"][mask_ct, 1],
            c=c, s=8, alpha=0.6, label=ct, rasterized=True
        )
    ax.set_title(sp_name.capitalize(), fontsize=13)
    ax.axis("off")
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map.get(ct, "#888888"),
                       markersize=8, label=ct) for ct in celltypes]
fig.legend(handles=handles, loc="center right", fontsize=9, frameon=False)
plt.suptitle("SATURN: Mouse + Human (split by species)", fontsize=13)
plt.tight_layout(rect=[0, 0, 0.88, 0.95])
plt.savefig("saturn_umap_split.png", dpi=300, bbox_inches="tight")
plt.close()
print("saved: saturn_umap_split.png")

# ============================================================
# Step 4: 保存完整结果 (含 UMAP)
# ============================================================
adata.write(os.path.join(SATURN_DATA, "saturn_mouse_human_result.h5ad"))
print(f"\n完整结果保存: {os.path.join(SATURN_DATA, 'saturn_mouse_human_result.h5ad')}")

# ============================================================
# Step 5: Macrogene 分析 (可选)
# ============================================================
if "macrogenes" in adata.obsm:
    print("\n分析 Macrogenes...")
    macrogenes = adata.obsm["macrogenes"]
    print(f"Macrogene 矩阵: {macrogenes.shape}")

    mg_adata = sc.AnnData(X=macrogenes, obs=adata.obs.copy())
    sc.pp.neighbors(mg_adata, n_neighbors=15)
    sc.tl.umap(mg_adata, min_dist=0.3)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ct in celltypes:
        mask = mg_adata.obs[label_col] == ct
        c = color_map.get(ct, None)
        axes[0].scatter(mg_adata.obsm["X_umap"][mask, 0], mg_adata.obsm["X_umap"][mask, 1],
                        c=c, s=8, alpha=0.6, label=ct, rasterized=True)
    axes[0].set_title("Macrogene space (by celltype)", fontsize=12)
    axes[0].axis("off")
    for sp_name in ["mouse", "human"]:
        mask = mg_adata.obs["species"] == sp_name
        axes[1].scatter(mg_adata.obsm["X_umap"][mask, 0], mg_adata.obsm["X_umap"][mask, 1],
                        c=species_colors[sp_name], s=8, alpha=0.5, label=sp_name.capitalize(), rasterized=True)
    axes[1].set_title("Macrogene space (by species)", fontsize=12)
    axes[1].axis("off")
    axes[1].legend(fontsize=10, frameon=False)
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map.get(ct, "#888888"),
                           markersize=8, label=ct) for ct in celltypes]
    axes[0].legend(handles=handles, fontsize=8, frameon=False, loc="upper left")
    plt.tight_layout()
    plt.savefig("saturn_macrogene_umap.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("saved: saturn_macrogene_umap.png")

print("\n" + "=" * 60)
print("SATURN 下游分析完成！")
print("=" * 60)


##### ===== 186行: 展示整合质量 ===== #####
import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
os.chdir("/projects/hmz-aging/Zebrafish_Mouse_Human_inte/saturn_data/saturn_results/")
adata = sc.read_h5ad("test256_data_mouse_saturn_zebrafish_saturn_org_saturn_seed_0.h5ad")
adata.obs["labels2"] = adata.obs["labels2"].astype(str)
adata.obs["species"] = adata.obs["species"].astype(str)
X_emb = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
labels = adata.obs["labels2"].values
species = adata.obs["species"].values
color_map = {"MG": "#C44E52", "AC": "#E5C953", "Cone": "#BCBD22", "RGC": "#B0B0B0", "HC": "#E5933A", "Rod": "#4878A8", "BC": "#996633", "RPE": "#B89BC7", "Microglia": "#8ECFC9"}
ct_order = sorted(set(labels))
sc.settings.figdir = "./"

##### 图1: 混淆矩阵 #####
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score
ms_mask = species == "mouse"
zf_mask = species == "zebrafish"
clf = LogisticRegression(max_iter=5000, random_state=42)
clf.fit(X_emb[ms_mask], labels[ms_mask])
zf_pred = clf.predict(X_emb[zf_mask])
zf_true = labels[zf_mask]
acc = accuracy_score(zf_true, zf_pred)
print(f"Cross-species transfer accuracy: {acc*100:.1f}%")
cm = confusion_matrix(zf_true, zf_pred, labels=ct_order)
cm_pct = cm / cm.sum(axis=1, keepdims=True) * 100
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(cm_pct, cmap="Blues", vmin=0, vmax=100)
ax.set_xticks(range(len(ct_order)))
ax.set_yticks(range(len(ct_order)))
ax.set_xticklabels(ct_order, rotation=45, ha="right", fontsize=10)
ax.set_yticklabels(ct_order, fontsize=10)

for i in range(len(ct_order)):
    for j in range(len(ct_order)):
        val = cm_pct[i, j]
        ax.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=9, color="white" if val > 50 else "black")

ax.set_xlabel("Predicted (from Mouse classifier)", fontsize=11)
ax.set_ylabel("True Zebrafish celltype", fontsize=11)
ax.set_title(f"Mouse\u2192Zebrafish label transfer (accuracy={acc*100:.1f}%)", fontsize=12)
plt.colorbar(im, ax=ax, label="%")
plt.tight_layout()
plt.savefig("saturn_confusion_matrix.png", dpi=300, bbox_inches="tight")
plt.close()
print("saved: saturn_confusion_matrix.png")

##### 图2: Harmony UMAP #####
import harmonypy as hm
X_emb_np = np.array(X_emb, dtype=np.float64)
ho = hm.run_harmony(X_emb_np, adata.obs, "species", max_iter_harmony=20)
Z = ho.Z_corr; Z = Z.T if Z.shape[0] != adata.n_obs else Z; adata.obsm["X_harmony"] = Z; print(f"X_harmony shape: {Z.shape}")
sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=30, metric="cosine")
sc.tl.umap(adata, min_dist=0.3, spread=1.0)
sc.settings.figdir = "./"
sc.set_figure_params(figsize=(5, 5))
adata.uns["labels2_colors"] = [color_map.get(ct, "#888888") for ct in ct_order]
sc.pl.umap(adata, color="labels2", title="SATURN + Harmony (celltype)", frameon=True, save="_harmony_celltype.png")
sc.pl.umap(adata, color="species", palette={"mouse": "#C44E52", "zebrafish": "#A0A226"}, title="SATURN + Harmony (species)", frameon=True, save="_harmony_species.png")
print("saved: umap_harmony_celltype.png, umap_harmony_species.png")





######
###### 
###### 