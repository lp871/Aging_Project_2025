#####
#####  Table 1 for cell numbers ######
#####  按 时间点 × celltype 统计细胞数（宽表 + 长表）
#####

suppressMessages(library(Seurat))

## 工具：长表 -> 宽表（行=时间点，列=celltype，含 Total）
make_age_celltype_tables <- function(meta, age_col, species_label) {
    tab_long <- as.data.frame(table(
        age_tp   = meta[[age_col]],
        celltype = meta$celltype
    ), stringsAsFactors = FALSE)
    colnames(tab_long) <- c(age_col, "celltype", "n_cells")
    tab_long[[age_col]] <- as.numeric(as.character(tab_long[[age_col]]))
    tab_long$n_cells <- as.integer(tab_long$n_cells)
    tab_long$species <- species_label
    tab_long <- tab_long[order(tab_long[[age_col]], tab_long$celltype),
                         c("species", age_col, "celltype", "n_cells")]
    rownames(tab_long) <- NULL

    fml <- as.formula(paste("n_cells ~", age_col, "+ celltype"))
    tab_wide <- as.data.frame.matrix(xtabs(fml, data = tab_long))
    tab_wide <- cbind(setNames(data.frame(as.numeric(rownames(tab_wide))), age_col), tab_wide)
    rownames(tab_wide) <- NULL
    ct_cols <- setdiff(colnames(tab_wide), age_col)
    tab_wide$Total <- rowSums(tab_wide[, ct_cols, drop = FALSE])

    col_tot <- colSums(tab_wide[, setdiff(colnames(tab_wide), age_col), drop = FALSE])
    tab_wide <- rbind(tab_wide, c(setNames(NA, age_col), col_tot))
    tab_print <- tab_wide
    tab_print[[age_col]] <- as.character(tab_print[[age_col]])
    tab_print[[age_col]][nrow(tab_print)] <- "Total"

    list(long = tab_long, wide = tab_print)
}


########################################################################
##### Zebrafish
##### sample: <age>_Rep<n> → age_months
##### 数据: Zebrafish_seurat_clcl_2026.rds（已去 6_Rep3, 3_Rep2, 24_Rep1）
########################################################################
setwd("/projects/hmz-aging/Zebrafish_New_Datasets")
Zebrafish_seurat_clcl <- readRDS("Zebrafish_seurat_clcl_2026.rds")

zf_meta <- Zebrafish_seurat_clcl@meta.data
zf_meta$age_months <- as.numeric(sub("_Rep.*$", "", as.character(zf_meta$sample)))
if (any(is.na(zf_meta$age_months))) {
    stop("Zebrafish: 有 sample 无法解析年龄: ",
         paste(unique(zf_meta$sample[is.na(zf_meta$age_months)]), collapse = ", "))
}

zf_tabs <- make_age_celltype_tables(zf_meta, "age_months", "Zebrafish")
table1_zebrafish_long <- zf_tabs$long
table1_zebrafish <- zf_tabs$wide

message("\n=== Zebrafish Table 1 (wide: age_months × celltype) ===")
print(table1_zebrafish, row.names = FALSE)

write.csv(table1_zebrafish, "Table1_Zebrafish_cellnumbers_age_by_celltype.csv", row.names = FALSE)
write.csv(table1_zebrafish_long, "Table1_Zebrafish_cellnumbers_age_by_celltype_long.csv", row.names = FALSE)
message("已保存 Zebrafish Table1 (宽/长表)")


########################################################################
##### Mouse
##### sample: Wk<age>_Rep<n> → age_weeks；106 → 108（与 FIgure1_2026 一致）
##### 数据: Mouse_seurat_clcl_2026.rds
########################################################################
setwd("/projects/hmz-aging/Mouse_New_Datasets")
Mouse_seurat_clcl <- readRDS("Mouse_seurat_clcl_2026.rds")

ms_meta <- Mouse_seurat_clcl@meta.data
ms_age <- sub("^Wk", "", as.character(ms_meta$sample))
ms_age <- sub("_Rep.*$", "", ms_age)
ms_meta$age_weeks <- as.numeric(ms_age)
ms_meta$age_weeks[ms_meta$age_weeks == 106] <- 108
if (any(is.na(ms_meta$age_weeks))) {
    stop("Mouse: 有 sample 无法解析年龄: ",
         paste(unique(ms_meta$sample[is.na(ms_meta$age_weeks)]), collapse = ", "))
}

ms_tabs <- make_age_celltype_tables(ms_meta, "age_weeks", "Mouse")
table1_mouse_long <- ms_tabs$long
table1_mouse <- ms_tabs$wide

message("\n=== Mouse Table 1 (wide: age_weeks × celltype) ===")
print(table1_mouse, row.names = FALSE)

write.csv(table1_mouse, "Table1_Mouse_cellnumbers_age_by_celltype.csv", row.names = FALSE)
write.csv(table1_mouse_long, "Table1_Mouse_cellnumbers_age_by_celltype_long.csv", row.names = FALSE)
message("已保存 Mouse Table1 (宽/长表)")

message("\n对象: table1_zebrafish / table1_zebrafish_long / table1_mouse / table1_mouse_long")

########################################################################
##### Mouse RPE datasets (/projects/hmz-aging/RPE_smoking)
##### sample: <age>W_RPE → age_weeks（如 5W_RPE, 49W_RPE, 91W_RPE）
##### celltype 重命名后按 时间点 × celltype 统计
########################################################################
suppressMessages(library(dplyr))

setwd("/projects/hmz-aging/RPE_smoking")
Mouse_RPE <- readRDS("Mouse_RPE_RNA_Final_Jul16")

message("\n=== Mouse RPE: sample / celltype (raw) ===")
print(table(Mouse_RPE$sample))
print(table(Mouse_RPE$celltype))

Mouse_RPE$celltype <- dplyr::recode(
    as.character(Mouse_RPE$celltype),
    "RPE" = "Retinal pigment epithelial cells",
    "Stromal" = "Stromal cells",
    "VE" = "Vascular endothelial cells",
    "Microglia" = "Macrophage",
    "ciliary_epithelial" = "Ciliary epithelial cells"
)
message("\n=== Mouse RPE: celltype (recoded) ===")
print(table(Mouse_RPE$celltype))

rpe_meta <- Mouse_RPE@meta.data
rpe_meta$age_weeks <- as.numeric(sub("^([0-9]+)W_.*$", "\\1", as.character(rpe_meta$sample)))
if (any(is.na(rpe_meta$age_weeks))) {
    stop("Mouse RPE: 有 sample 无法解析年龄: ",
         paste(unique(rpe_meta$sample[is.na(rpe_meta$age_weeks)]), collapse = ", "))
}

rpe_tabs <- make_age_celltype_tables(rpe_meta, "age_weeks", "Mouse_RPE")
table1_mouse_rpe_long <- rpe_tabs$long
table1_mouse_rpe <- rpe_tabs$wide

message("\n=== Mouse RPE Table 1 (wide: age_weeks × celltype) ===")
print(table1_mouse_rpe, row.names = FALSE)

write.csv(table1_mouse_rpe,
          "Table1_Mouse_RPE_cellnumbers_age_by_celltype.csv",
          row.names = FALSE)
write.csv(table1_mouse_rpe_long,
          "Table1_Mouse_RPE_cellnumbers_age_by_celltype_long.csv",
          row.names = FALSE)
message("已保存 Mouse RPE Table1 (宽/长表) -> /projects/hmz-aging/RPE_smoking/")

message("\n对象追加: table1_mouse_rpe / table1_mouse_rpe_long")


########################################################################
##### Human clean5 Table1（下面整段在 python 中运行）
##### 主表: sample meta + AC/BC/RGC/Rod/Cone/HC/MG/Microglia + Total
##### RPE 单独表
##### meta: sampleid, donor_id, sex, age, age_numeric, age_interval_10yr, study_name
##### 规则: clean5 有 → 计数；原始有、clean5 去掉 → NA；从未有过 → 0
##### 输出:
#####   /projects/hmz-aging/Human_New_Datasets/celltype_h5ad/
#####     Table1_Human_clean5_sample_by_celltype.csv
#####     Table1_Human_RPE_clean5_sample.csv
########################################################################
# conda activate scprint
# python

conda activate 


import os
import numpy as np
import pandas as pd
import scanpy as sc

h5ad_dir = "/projects/hmz-aging/Human_New_Datasets/celltype_h5ad"
out_dir = h5ad_dir
major_cts = ["AC", "BC", "RGC", "Rod", "Cone", "HC", "MG", "Microglia"]

META_COLS = ["sampleid", "donor_id", "sex", "age", "age_numeric", "age_interval_10yr", "study_name"]


def parse_age_numeric(age_series):
    s = age_series
    if isinstance(s, pd.DataFrame):
        s = s.iloc[:, 0]
    s = pd.Series(s).astype(str)
    return pd.to_numeric(s.str.replace(r"[^0-9.]", "", regex=True), errors="coerce")


def add_age_interval_10yr(age_numeric):
    #### 与 Aging_Clocks_Functions.R / seurat_add_age_numeric_and_interval 一致 ####
    return pd.cut(
        age_numeric,
        bins=list(range(0, 100, 10)),
        labels=["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90"],
        right=False,
        include_lowest=True,
    ).astype(str).replace("nan", np.nan)


def normalize_obs_meta(obs):
    """统一原始 h5ad / clean5 的 sample 级 meta 列名。"""
    obs = obs.copy()
    obs.columns = [str(c) for c in obs.columns]
    #### 列名对齐：clean5=age；原始 retina=donor_age；RPE 原始可能 Gender/Age/Donor ####
    #### 注意 RPE.h5ad 常同时有 Age + donor_age，只能保留一列 age ####
    rename = {}
    if "donor_id" not in obs.columns and "Donor" in obs.columns:
        rename["Donor"] = "donor_id"
    if "sex" not in obs.columns and "Gender" in obs.columns:
        rename["Gender"] = "sex"
    if rename:
        obs = obs.rename(columns=rename)
    if "age" not in obs.columns:
        if "donor_age" in obs.columns:
            obs["age"] = obs["donor_age"]
        elif "Age" in obs.columns:
            obs["age"] = obs["Age"]
    #### 去掉重复列名（rename 撞名时可能出现）####
    obs = obs.loc[:, ~obs.columns.duplicated()].copy()
    if "sex" in obs.columns:
        obs["sex"] = obs["sex"].astype(str).str.lower()
    for col in ["sampleid", "donor_id", "sex", "age", "study_name"]:
        if col not in obs.columns:
            obs[col] = np.nan
        else:
            obs[col] = pd.Series(obs[col]).astype(str)
            obs[col] = obs[col].replace({"nan": np.nan, "None": np.nan, "": np.nan})
    #### age 统一成 "XX years"；纯数字则补后缀；>89 -> 90（与上游处理一致）####
    age = pd.Series(obs["age"]).astype(str)
    age = age.replace({">89 years": "90 years", ">89": "90 years"})
    num = parse_age_numeric(age)
    obs["age"] = [
        (f"{int(x) if float(x).is_integer() else x} years" if pd.notna(x) else a)
        for x, a in zip(num, age)
    ]
    obs["age_numeric"] = parse_age_numeric(obs["age"])
    obs["age_interval_10yr"] = add_age_interval_10yr(obs["age_numeric"])
    return obs




def load_sample_level(h5ad_path):
    """读 obs -> sample 级 counts + meta（每个 sampleid 一行 meta）。"""
    a = sc.read_h5ad(h5ad_path, backed="r")
    want = [c for c in ["sampleid", "donor_id", "sex", "age", "donor_age",
                        "study_name", "Donor", "Gender", "Age"] if c in a.obs.columns]
    if "sampleid" not in want:
        a.file.close()
        raise ValueError(f"{h5ad_path}: 缺少 sampleid")
    obs = normalize_obs_meta(a.obs[want].copy())
    a.file.close()
    counts = obs.groupby("sampleid", observed=True).size().astype(float)
    meta = (obs.groupby("sampleid", observed=True)
            .agg({
                "donor_id": "first",
                "sex": "first",
                "age": "first",
                "age_numeric": "first",
                "age_interval_10yr": "first",
                "study_name": "first",
            }))
    return counts, meta


def coalesce_meta(meta_prefer, meta_fallback, sample_ids):
    """优先用 clean5 meta，缺失再用原始 meta 补。"""
    rows = []
    for sid in sample_ids:
        row = {"sampleid": sid}
        src = None
        if sid in meta_prefer.index:
            src = meta_prefer.loc[sid]
        elif sid in meta_fallback.index:
            src = meta_fallback.loc[sid]
        for col in ["donor_id", "sex", "age", "age_numeric", "age_interval_10yr", "study_name"]:
            val = np.nan
            if src is not None and col in src.index and pd.notna(src[col]) and str(src[col]) not in ("nan", "None", ""):
                val = src[col]
            elif sid in meta_fallback.index and col in meta_fallback.columns:
                fb = meta_fallback.loc[sid, col]
                if pd.notna(fb) and str(fb) not in ("nan", "None", ""):
                    val = fb
            row[col] = val
        rows.append(row)
    return pd.DataFrame(rows).set_index("sampleid")


def build_table(celltypes, label):
    count_map, orig_set = {}, {}
    meta_c5_all, meta_orig_all = [], []
    for ct in celltypes:
        orig_f = os.path.join(h5ad_dir, f"{ct}.h5ad")
        c5_f = os.path.join(h5ad_dir, f"{ct}_clean5.h5ad")
        if not os.path.exists(c5_f):
            raise FileNotFoundError(f"缺少 clean5: {c5_f}")
        if not os.path.exists(orig_f):
            raise FileNotFoundError(f"缺少原始 h5ad: {orig_f}")
        print(f"\n[{label}] {ct}")
        c5_counts, c5_meta = load_sample_level(c5_f)
        orig_counts, orig_meta = load_sample_level(orig_f)
        count_map[ct] = c5_counts
        orig_set[ct] = set(orig_counts.index.astype(str))
        removed = sorted(orig_set[ct] - set(c5_counts.index.astype(str)))
        print(f"  clean5 samples={len(c5_counts)}, orig={len(orig_set[ct])}, removed={len(removed)}")
        if removed:
            print(f"  removed: {removed}")
        meta_c5_all.append(c5_meta)
        meta_orig_all.append(orig_meta)
    meta_c5 = pd.concat(meta_c5_all).groupby(level=0).first()
    meta_orig = pd.concat(meta_orig_all).groupby(level=0).first()
    all_samples = sorted(set().union(*orig_set.values()))
    sample_meta = coalesce_meta(meta_c5, meta_orig, all_samples)
    rows = []
    for sid in all_samples:
        row = sample_meta.loc[sid].to_dict()
        row["sampleid"] = sid
        for ct in celltypes:
            if sid in count_map[ct].index:
                row[ct] = int(count_map[ct].loc[sid])
            elif sid in orig_set[ct]:
                row[ct] = np.nan
            else:
                row[ct] = 0
        vals = [row[ct] for ct in celltypes]
        row["Total"] = int(np.nansum(vals)) if not all(pd.isna(v) for v in vals) else np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    df = df.sort_values(["donor_id", "sampleid"], kind="mergesort").reset_index(drop=True)
    df = df[META_COLS + celltypes + ["Total"]]
    return df


print("=" * 70)
print("Human Table1 from clean5 (+ sample metadata)")
print("=" * 70)

table_major = build_table(major_cts, "major")
out_major = os.path.join(out_dir, "Table1_Human_clean5_sample_by_celltype.csv")
table_major.to_csv(out_major, index=False, na_rep="NA")
print(f"\n已保存主表: {out_major}")
print(f"  rows={table_major.shape[0]}, cols={list(table_major.columns)}")
print(table_major.head(8).to_string(index=False))

table_rpe = build_table(["RPE"], "RPE")
out_rpe = os.path.join(out_dir, "Table1_Human_RPE_clean5_sample.csv")
table_rpe.to_csv(out_rpe, index=False, na_rep="NA")
print(f"\n已保存 RPE 表: {out_rpe}")
print(f"  rows={table_rpe.shape[0]}, cols={list(table_rpe.columns)}")
print(table_rpe.head(8).to_string(index=False))

#### donor 数量汇总 ####
n_donor_retina = table_major["donor_id"].nunique(dropna=True)
n_donor_rpe = table_rpe["donor_id"].nunique(dropna=True)
n_donor_union = pd.Series(
    pd.concat([table_major["donor_id"], table_rpe["donor_id"]], ignore_index=True)
).nunique(dropna=True)
n_donor_overlap = len(
    set(table_major["donor_id"].dropna().astype(str))
    & set(table_rpe["donor_id"].dropna().astype(str))
)
n_cells_retina = int(table_major["Total"].fillna(0).sum())
n_cells_rpe = int(table_rpe["Total"].fillna(0).sum())
n_cells_human = n_cells_retina + n_cells_rpe

age_retina = pd.to_numeric(table_major["age_numeric"], errors="coerce").dropna()
age_rpe = pd.to_numeric(table_rpe["age_numeric"], errors="coerce").dropna()
age_all = pd.concat([age_retina, age_rpe], ignore_index=True)

print("\n" + "=" * 70)
print("Donor / Cell / Age 汇总")
print("=" * 70)
print(f"  Retina (major 8 CT): {n_donor_retina} donors  ({table_major['sampleid'].nunique()} samples)  {n_cells_retina:,} cells")
print(f"  RPE:                 {n_donor_rpe} donors  ({table_rpe['sampleid'].nunique()} samples)  {n_cells_rpe:,} cells")
print(f"  两者并集 donors:     {n_donor_union}")
print(f"  两者交集 donors:     {n_donor_overlap}")
print(f"  Human 总计 cells:    {n_cells_human:,}  (= Retina {n_cells_retina:,} + RPE {n_cells_rpe:,})")
print(f"  Retina age range:    {age_retina.min():.0f} – {age_retina.max():.0f} years  (n_samples with age={len(age_retina)})")
print(f"  RPE age range:       {age_rpe.min():.0f} – {age_rpe.max():.0f} years  (n_samples with age={len(age_rpe)})")
print(f"  Human age range:     {age_all.min():.0f} – {age_all.max():.0f} years")
print(f"  Retina age_interval_10yr: {sorted(table_major['age_interval_10yr'].dropna().unique().tolist())}")
print(f"  RPE age_interval_10yr:    {sorted(table_rpe['age_interval_10yr'].dropna().unique().tolist())}")

print("\n=== Human Table1 完成 ===")


########################################################################
##### TableS1: combine all Table1 sheets into one Excel (first tab = Contents)
##### Requires: Table1 CSVs already generated; run this block in python
##### Output: /projects/hmz-aging/Human_New_Datasets/TableS1_cellnumbers.xlsx
########################################################################
# conda activate scprint
# python

import os
import pandas as pd

out_xlsx = "/projects/hmz-aging/Human_New_Datasets/TableS1_cellnumbers.xlsx"

#### sheet paths (wide tables; Human is sample-level) ####
sheet_files = [
    {
        "sheet": "Zebrafish",
        "file": "/projects/hmz-aging/Zebrafish_New_Datasets/Table1_Zebrafish_cellnumbers_age_by_celltype.csv",
        "description": "Zebrafish retina: cell numbers by age (months) x cell type.",
    },
    {
        "sheet": "Mouse_Retina",
        "file": "/projects/hmz-aging/Mouse_New_Datasets/Table1_Mouse_cellnumbers_age_by_celltype.csv",
        "description": "Mouse retina: cell numbers by age (weeks) x cell type.",
    },
    {
        "sheet": "Mouse_RPE",
        "file": "/projects/hmz-aging/RPE_smoking/Table1_Mouse_RPE_cellnumbers_age_by_celltype.csv",
        "description": "Mouse RPE/choroid: cell numbers by age (weeks) x cell type.",
    },
    {
        "sheet": "Human_Retina",
        "file": "/projects/hmz-aging/Human_New_Datasets/celltype_h5ad/Table1_Human_clean5_sample_by_celltype.csv",
        "description": "Human retina (clean5): per-sample metadata and cell numbers for AC, BC, RGC, Rod, Cone, HC, MG, and Microglia. NA indicates the sample was removed for that cell type during QC.",
    },
    {
        "sheet": "Human_RPE",
        "file": "/projects/hmz-aging/Human_New_Datasets/celltype_h5ad/Table1_Human_RPE_clean5_sample.csv",
        "description": "Human RPE (clean5): per-sample metadata and RPE cell numbers. NA indicates the sample was removed during QC.",
    },
]

#### load tables ####
tables = {}
for i, info in enumerate(sheet_files, start=1):
    if not os.path.exists(info["file"]):
        raise FileNotFoundError(
            f"Missing table file (run the corresponding Table1 section first): {info['file']}"
        )
    df = pd.read_csv(info["file"])
    tables[info["sheet"]] = df
    print(f"  [{i}] {info['sheet']}: {df.shape[0]} rows x {df.shape[1]} cols <- {info['file']}")

#### Contents tab ####
contents = pd.DataFrame({
    "Sheet": [x["sheet"] for x in sheet_files],
    "Description": [x["description"] for x in sheet_files],
    "N_rows": [tables[x["sheet"]].shape[0] for x in sheet_files],
    "N_cols": [tables[x["sheet"]].shape[1] for x in sheet_files],
    "Source_file": [x["file"] for x in sheet_files],
})

#### write Excel; Contents is the first sheet ####
with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    contents.to_excel(writer, sheet_name="Contents", index=False)
    for info in sheet_files:
        tables[info["sheet"]].to_excel(writer, sheet_name=info["sheet"], index=False)

print(f"\nSaved TableS1: {out_xlsx}")
print("Sheets:", ["Contents"] + [x["sheet"] for x in sheet_files])
print(contents.to_string(index=False))
print("\n=== TableS1 done ===")

#######
#######
#######


