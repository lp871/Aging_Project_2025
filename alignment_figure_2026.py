import os, scanpy as sc, numpy as np, pandas as pd, harmonypy as hm
from sklearn.neighbors import NearestNeighbors
os.chdir("/projects/hmz-aging/Zebrafish_Mouse_Human_inte/")
color_map = {"MG": "#C44E52", "AC": "#E5C953", "Cone": "#BCBD22", "RGC": "#B0B0B0", "HC": "#E5933A", "Rod": "#4878A8", "BC": "#996633", "RPE": "#B89BC7", "Microglia": "#8ECFC9"}

print("=== Loading Mouse-Zebrafish SATURN result ===")
adata_mz = sc.read_h5ad("saturn_data/saturn_results/test256_data_mouse_saturn_zebrafish_saturn_org_saturn_seed_0.h5ad")
adata_mz.obs["labels2"] = adata_mz.obs["labels2"].astype(str)
adata_mz.obs["species"] = adata_mz.obs["species"].astype(str)
X_mz = adata_mz.X if not hasattr(adata_mz.X, 'toarray') else adata_mz.X.toarray()
X_mz_np = np.array(X_mz, dtype=np.float64)
ho_mz = hm.run_harmony(X_mz_np, adata_mz.obs, "species", max_iter_harmony=20)
Z_mz = ho_mz.Z_corr
Z_mz = Z_mz.T if Z_mz.shape[0] != adata_mz.n_obs else Z_mz
print(f"Mouse-Zebrafish Harmony shape: {Z_mz.shape}")

ms_mask_mz = adata_mz.obs["species"].values == "mouse"
zf_mask_mz = adata_mz.obs["species"].values == "zebrafish"
H_ms_mz = Z_mz[ms_mask_mz]
H_zf = Z_mz[zf_mask_mz]
ms_labels_mz = adata_mz.obs["labels2"].values[ms_mask_mz]
zf_labels = adata_mz.obs["labels2"].values[zf_mask_mz]
nn_zf = NearestNeighbors(n_neighbors=20, metric="cosine")
nn_zf.fit(H_zf)
_, indices_ms2zf = nn_zf.kneighbors(H_ms_mz)
ms_to_zf = [{"ms_ct": ms_labels_mz[i], "zf_ct": pd.Series(zf_labels[indices_ms2zf[i]]).value_counts().index[0]} for i in range(len(ms_labels_mz))]
ms_zf_counts = pd.DataFrame(ms_to_zf).groupby(["ms_ct", "zf_ct"]).size().reset_index(name="count")
print("\nMouse -> Zebrafish alignment:")
print(ms_zf_counts.to_string(index=False))

print("\n=== Loading Mouse-Human SATURN result ===")
mh_dir = "saturn_data_mh/saturn_results/"
mh_files = sorted([f for f in os.listdir(mh_dir) if f.endswith(".h5ad") and "pretrain" not in f and "ep_" not in f])
print(f"Using: {mh_files[-1]}")
adata_mh = sc.read_h5ad(os.path.join(mh_dir, mh_files[-1]))
adata_mh.obs["labels2"] = adata_mh.obs["labels2"].astype(str)
adata_mh.obs["species"] = adata_mh.obs["species"].astype(str)
X_mh = adata_mh.X if not hasattr(adata_mh.X, 'toarray') else adata_mh.X.toarray()
X_mh_np = np.array(X_mh, dtype=np.float64)
ho_mh = hm.run_harmony(X_mh_np, adata_mh.obs, "species", max_iter_harmony=20)
Z_mh = ho_mh.Z_corr
Z_mh = Z_mh.T if Z_mh.shape[0] != adata_mh.n_obs else Z_mh
print(f"Mouse-Human Harmony shape: {Z_mh.shape}")

ms_mask_mh = adata_mh.obs["species"].values == "mouse"
hu_mask_mh = adata_mh.obs["species"].values == "human"
H_ms_mh = Z_mh[ms_mask_mh]
H_hu = Z_mh[hu_mask_mh]
ms_labels_mh = adata_mh.obs["labels2"].values[ms_mask_mh]
hu_labels = adata_mh.obs["labels2"].values[hu_mask_mh]
nn_hu = NearestNeighbors(n_neighbors=20, metric="cosine")
nn_hu.fit(H_hu)
_, indices_ms2hu = nn_hu.kneighbors(H_ms_mh)
ms_to_hu = [{"ms_ct": ms_labels_mh[i], "hu_ct": pd.Series(hu_labels[indices_ms2hu[i]]).value_counts().index[0]} for i in range(len(ms_labels_mh))]
ms_hu_counts = pd.DataFrame(ms_to_hu).groupby(["ms_ct", "hu_ct"]).size().reset_index(name="count")
print("\nMouse -> Human alignment:")
print(ms_hu_counts.to_string(index=False))

##### 保存为 CSV，供 R ggplot2 画冲击图 #####
ms_zf_out = ms_zf_counts.copy()
ms_zf_out.columns = ["Mouse", "Zebrafish", "count"]
ms_zf_out["pair"] = "Mouse_Zebrafish"
ms_hu_out = ms_hu_counts.copy()
ms_hu_out.columns = ["Mouse", "Human", "count"]
ms_hu_out["pair"] = "Mouse_Human"

ms_zf_out.to_csv("alignment_mouse_zebrafish.csv", index=False)
ms_hu_out.to_csv("alignment_mouse_human.csv", index=False)
print("\nsaved: alignment_mouse_zebrafish.csv")
print("saved: alignment_mouse_human.csv")

# 也保存一个合并的长格式表 (用于 ggalluvial)
long_df = []
for _, row in ms_zf_counts.iterrows():
    for _ in range(row["count"]): long_df.append({"Zebrafish": row["zf_ct"], "Mouse": row["ms_ct"]})

long_mz = pd.DataFrame(long_df)

long_df2 = []
for _, row in ms_hu_counts.iterrows():
    for _ in range(row["count"]): long_df2.append({"Mouse": row["ms_ct"], "Human": row["hu_ct"]})

long_mh = pd.DataFrame(long_df2)

long_mz.to_csv("alignment_long_mouse_zebrafish.csv", index=False)
long_mh.to_csv("alignment_long_mouse_human.csv", index=False)
print("saved: alignment_long_mouse_zebrafish.csv")
print("saved: alignment_long_mouse_human.csv")


##### R code (ggalluvial) #####
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
setwd("/projects/hmz-aging/Zebrafish_Mouse_Human_inte/")
mz <- read.csv("alignment_mouse_zebrafish.csv")
mh <- read.csv("alignment_mouse_human.csv")
alluvial_data <- expand.grid(Zebrafish = unique(mz$Zebrafish), Mouse = unique(c(mz$Mouse, mh$Mouse)), Human = unique(mh$Human), stringsAsFactors = FALSE)
alluvial_data <- alluvial_data %>% left_join(mz %>% select(Mouse, Zebrafish, count) %>% rename(count_zf = count), by = c("Mouse", "Zebrafish")) %>% left_join(mh %>% select(Mouse, Human, count) %>% rename(count_hu = count), by = c("Mouse", "Human")) %>% mutate(count_zf = replace_na(count_zf, 0), count_hu = replace_na(count_hu, 0)) %>% mutate(Freq = sqrt(count_zf * count_hu)) %>% filter(Freq > 0)
color_map <- c("MG" = "#C44E52", "AC" = "#E5C953", "Cone" = "#BCBD22", "RGC" = "#B0B0B0", "HC" = "#E5933A", "Rod" = "#4878A8", "BC" = "#996633", "RPE" = "#B89BC7", "Microglia" = "#8ECFC9")
real_cts <- sort(unique(c(as.character(mz$Zebrafish), as.character(mz$Mouse), as.character(mh$Mouse), as.character(mh$Human))))
spacer_size <- sum(alluvial_data$Freq) * 0.03
gap_names <- paste0("_gap", 1:(length(real_cts)-1))
spacer_rows <- data.frame(Zebrafish = gap_names, Mouse = gap_names, Human = gap_names, count_zf = 0, count_hu = 0, Freq = spacer_size, stringsAsFactors = FALSE)
alluvial_data <- rbind(alluvial_data, spacer_rows)
new_levels <- as.vector(rbind(real_cts, c(gap_names, NA)))[1:(2*length(real_cts)-1)]
alluvial_data$Zebrafish <- factor(alluvial_data$Zebrafish, levels = new_levels)
alluvial_data$Mouse <- factor(alluvial_data$Mouse, levels = new_levels)
alluvial_data$Human <- factor(alluvial_data$Human, levels = new_levels)
gap_colors <- setNames(rep("white", length(gap_names)), gap_names)
color_map_full <- c(color_map, gap_colors)
p <- ggplot(alluvial_data, aes(axis1 = Zebrafish, axis2 = Mouse, axis3 = Human, y = Freq)) + geom_alluvium(aes(fill = Mouse), width = 1/5, alpha = 0.3, knot.pos = 0.4) + geom_stratum(aes(fill = after_stat(stratum)), width = 1/5, color = NA) + geom_text(stat = "stratum", aes(label = ifelse(grepl("^_gap", after_stat(stratum)), "", as.character(after_stat(stratum)))), size = 3.5, fontface = "bold") + scale_x_discrete(limits = c("Zebrafish", "Mouse", "Human"), expand = c(0.1, 0.05)) + scale_fill_manual(values = color_map_full, na.value = "white") + labs(title = "Cross-species Cell Type Alignment", subtitle = "SATURN + Harmony kNN (k=20), Mouse as reference") + theme_minimal(base_size = 12) + theme(legend.position = "none", panel.grid = element_blank(), axis.text.y = element_blank(), axis.title = element_blank(), axis.text.x = element_text(face = "bold", size = 13), plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))
ggsave("alignment_sankey.png", p, width = 8, height = 7, dpi = 300)


########