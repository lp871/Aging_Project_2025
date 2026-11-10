########
######## 之前人的 aging clocks 是怎么处理的 ？？？ ####
########

########
######## Human 预处理（基于 clean5）########
######## clean5: X = raw counts；基因已跨 celltype 对齐；obs 列: sex, donor_id, sampleid, age, celltype ########
######## 流程：clean5 -> 每个 sampleid 内 KNN smoothing -> 加 subtype -> 导出 mtx -> R 建 Seurat ########
########

########################################################################
######## Python (scanpy): 每个 sampleid 内 KNN smoothing ########
########################################################################
conda activate scprint
python

import os
import gc
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


try:
    from tqdm.auto import tqdm
except ImportError:
    def tqdm(it, desc="", **kwargs):
        return it
try:
    from joblib import Parallel, delayed
except ImportError:
    Parallel = None
    def delayed(f):
        return f


def _smooth_one_sample(adata, sid, sample_col, max_k, n_pcs, min_cells, max_cells_per_sample, random_state, X_src, is_sparse):
    """处理单个 sample：KNN 用该 sample 内全部细胞，输出最多 max_cells_per_sample 个细胞（不放回抽样）。"""
    cell_mask = adata.obs[sample_col] == sid
    cell_indices = np.where(cell_mask)[0]
    n_cells = len(cell_indices)
    if is_sparse:
        sub_raw = X_src[cell_indices, :].toarray()
    else:
        sub_raw = np.asarray(X_src[cell_indices, :], dtype=np.float32)
    if n_cells < min_cells:
        return adata.obs.iloc[cell_indices].copy(), sub_raw.astype(np.float32), sub_raw.astype(np.float32), n_cells, 0
    adata_sub = adata[cell_mask, :].copy()
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    actual_n_pcs = min(n_pcs, n_cells - 1, adata_sub.n_vars - 1)
    if actual_n_pcs < 5:
        return adata.obs.iloc[cell_indices].copy(), sub_raw.astype(np.float32), sub_raw.astype(np.float32), n_cells, 0
    sc.tl.pca(adata_sub, n_comps=actual_n_pcs)
    pca_coords = adata_sub.obsm['X_pca'].copy()
    del adata_sub
    actual_k = max(1, min(max_k, n_cells - 1))
    knn = NearestNeighbors(n_neighbors=actual_k + 1, algorithm='auto', metric='euclidean')
    knn.fit(pca_coords)
    distances, indices = knn.kneighbors(pca_coords)
    del pca_coords
    row_i = np.repeat(np.arange(n_cells), indices.shape[1])
    col_j = indices.ravel()
    W = sparse.csr_matrix(
        (np.ones(len(row_i)), (row_i, col_j)),
        shape=(n_cells, n_cells), dtype=np.float32
    )
    sub_smoothed = (W @ sub_raw).astype(np.float32)
    del W
    if max_cells_per_sample is not None and n_cells > max_cells_per_sample:
        seed = (random_state + hash(str(sid))) % (2**32)
        rng = np.random.default_rng(seed)
        keep_idx = rng.choice(n_cells, max_cells_per_sample, replace=False)
    else:
        keep_idx = np.arange(n_cells)
    obs_kept = adata.obs.iloc[cell_indices[keep_idx]].copy()
    smoothed_kept = sub_smoothed[keep_idx, :]
    raw_kept = sub_raw[keep_idx, :]
    return obs_kept, smoothed_kept, raw_kept, n_cells, actual_k


def knn_smooth_by_sample(adata, sample_col='sampleid', n_neighbors=15, n_pcs=30, min_cells=50, progress_desc="", max_cells_per_sample=None, n_jobs=1, random_state=42):
    """
    对 adata 内的每个 sample 分别进行 KNN smoothing。
    KNN 始终在该 sample 内全部细胞上计算；若设 max_cells_per_sample，则每个 sample 最多保留该数量细胞（不放回抽样）。
    - adata.X: smoothed counts
    - adata.layers['raw']: 原始 raw counts
    - adata.layers['smoothed_counts']: smoothed counts (同 X)
    """
    print("="*70)
    print("KNN Smoothing by Sample")
    print("="*70)
    print(f"Input adata: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"Sample column: {sample_col}")
    max_k = int(n_neighbors)
    print(f"max_k: {max_k}, n_pcs: {n_pcs}, min_cells: {min_cells}")
    if max_cells_per_sample is not None:
        print(f"max_cells_per_sample: {max_cells_per_sample} (不放回抽样)")
        print(f"n_jobs: {n_jobs}")
    X_src = adata.X
    is_sparse = sparse.issparse(X_src)
    sample_ids = adata.obs[sample_col].unique()
    print(f"\nTotal samples: {len(sample_ids)}")
    if max_cells_per_sample is not None and n_jobs > 1 and Parallel is not None:
        n_jobs = min(n_jobs, len(sample_ids))
        def run_one_sample(sid):
            try:
                o, s, r, nc, k = _smooth_one_sample(adata, sid, sample_col, max_k, n_pcs, min_cells, max_cells_per_sample, random_state, X_src, is_sparse)
                return (o, s, r)
            except Exception:
                return (None, None, None)
        desc_par = f"{progress_desc} samples" if progress_desc else "samples"
        results = Parallel(n_jobs=n_jobs, backend="loky")(delayed(run_one_sample)(sid) for sid in tqdm(sample_ids, desc=desc_par))
        obs_list, smoothed_list, raw_list = [], [], []
        for (o, s, r) in results:
            if o is not None:
                obs_list.append(o)
                smoothed_list.append(s)
                raw_list.append(r)
        if not obs_list:
            raise RuntimeError("No samples processed successfully.")
        obs_out = pd.concat(obs_list, axis=0)
        X_smoothed = np.vstack(smoothed_list).astype(np.float32)
        X_raw = np.vstack(raw_list).astype(np.float32)
        adata_smoothed = sc.AnnData(X=sparse.csr_matrix(X_smoothed), obs=obs_out, var=adata.var.copy())
        adata_smoothed.layers['raw'] = X_raw
        adata_smoothed.layers['smoothed_counts'] = adata_smoothed.X
        print(f"\nOutput adata: {adata_smoothed.shape[0]} cells x {adata_smoothed.shape[1]} genes (after max_cells_per_sample)")
        return adata_smoothed
    # 串行分支（不抽样 / n_jobs=1）
    if is_sparse:
        smoothed_X = sparse.lil_matrix(adata.X.shape, dtype=np.float32)
    else:
        smoothed_X = np.zeros(adata.X.shape, dtype=np.float32)
    n_smoothed, n_skipped, k_used_list = 0, 0, []
    desc = f"{progress_desc} samples" if progress_desc else "samples"
    for sid in tqdm(sample_ids, desc=desc):
        cell_mask = adata.obs[sample_col] == sid
        cell_indices = np.where(cell_mask)[0]
        n_cells = len(cell_indices)
        if is_sparse:
            sub_raw = X_src[cell_indices, :].toarray()
        else:
            sub_raw = np.asarray(X_src[cell_indices, :], dtype=np.float32)
        if n_cells < min_cells:
            smoothed_X[cell_indices, :] = sub_raw
            n_skipped += 1
            del sub_raw
            gc.collect()
            continue
        adata_sub = adata[cell_mask, :].copy()
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)
        actual_n_pcs = min(n_pcs, n_cells - 1, adata_sub.n_vars - 1)
        if actual_n_pcs < 5:
            smoothed_X[cell_indices, :] = sub_raw
            n_skipped += 1
            del adata_sub, sub_raw
            gc.collect()
            continue
        sc.tl.pca(adata_sub, n_comps=actual_n_pcs)
        pca_coords = adata_sub.obsm['X_pca'].copy()
        del adata_sub
        gc.collect()
        actual_k = max(1, min(max_k, n_cells - 1))
        k_used_list.append((sid, n_cells, actual_k))
        knn = NearestNeighbors(n_neighbors=actual_k + 1, algorithm='auto', metric='euclidean')
        knn.fit(pca_coords)
        distances, indices = knn.kneighbors(pca_coords)
        del pca_coords
        row_i = np.repeat(np.arange(n_cells), indices.shape[1])
        col_j = indices.ravel()
        W = sparse.csr_matrix(
            (np.ones(len(row_i)), (row_i, col_j)),
            shape=(n_cells, n_cells), dtype=np.float32
        )
        sub_smoothed = (W @ sub_raw).astype(np.float32)
        smoothed_X[cell_indices, :] = sub_smoothed
        n_smoothed += 1
        del sub_raw, sub_smoothed, W
        gc.collect()
    print(f"\nSmoothed samples: {n_smoothed}")
    print(f"Skipped samples (< {min_cells} cells): {n_skipped}")
    if k_used_list:
        k_vals = [x[2] for x in k_used_list]
        print(f"Dynamic k: min={min(k_vals)}, max={max(k_vals)}, mean={np.mean(k_vals):.1f}")
    if is_sparse:
        smoothed_X = smoothed_X.tocsr()
    else:
        smoothed_X = sparse.csr_matrix(smoothed_X)
    adata_smoothed = sc.AnnData(X=smoothed_X, obs=adata.obs.copy(), var=adata.var.copy())
    adata_smoothed.layers['raw'] = adata.X.copy()
    adata_smoothed.layers['smoothed_counts'] = adata_smoothed.X
    print(f"\nOutput adata: {adata_smoothed.shape[0]} cells x {adata_smoothed.shape[1]} genes")
    return adata_smoothed


########################################################################
######## 批量处理所有 celltype 的 clean5 -> smoothed ########
########################################################################
def smooth_all_celltypes(h5ad_dir, celltypes, sample_col='sampleid', n_neighbors=15, n_pcs=30, min_cells=50, max_cells_per_sample=None, n_jobs=1, random_state=42):
    """批量对所有 celltype 的 clean5 进行 KNN smoothing，输出 *_clean5_smoothed.h5ad。"""
    for ct in celltypes:
        clean5_file = os.path.join(h5ad_dir, f"{ct}_clean5.h5ad")
        if not os.path.exists(clean5_file):
            print(f"\n{ct}: clean5 not found, skipping -> {clean5_file}")
            continue
        print(f"\n{'='*70}\nProcessing: {ct}\n{'='*70}")
        adata = sc.read_h5ad(clean5_file)
        adata.obs['celltype'] = ct  # 用文件名统一 celltype
        adata_smoothed = knn_smooth_by_sample(adata, sample_col=sample_col, n_neighbors=n_neighbors, n_pcs=n_pcs, min_cells=min_cells, progress_desc=ct, max_cells_per_sample=max_cells_per_sample, n_jobs=n_jobs, random_state=random_state)
        output_file = os.path.join(h5ad_dir, f"{ct}_clean5_smoothed.h5ad")
        adata_smoothed.write_h5ad(output_file)
        print(f"Saved: {output_file}")
        del adata, adata_smoothed
        gc.collect()
    print("\n=== All celltypes smoothed! ===")


h5ad_dir = "/projects/hmz-aging/Human_New_Datasets/celltype_h5ad"
celltypes = ['AC', 'BC', 'RGC', 'Rod', 'Cone', 'HC', 'MG', 'Microglia', 'RPE']

# 每 sample 内 KNN(k=15, n_pcs=30, min_cells=50) 平滑，每 sample 最多保留 100 细胞（不放回），5 核并行
smooth_all_celltypes(h5ad_dir, celltypes, sample_col='sampleid', n_neighbors=15, n_pcs=30, min_cells=50, max_cells_per_sample=100, n_jobs=5)


########################################################################
######## 验证：smooth 后每个细胞的 total counts >= smooth 前（抽 3 个 celltype）########
########################################################################
def check_smooth_total_counts(h5ad_dir, celltypes_to_check=None):
    """检查 smooth 后 adata.X 的 total counts（按细胞求和）是否 >= layers['raw'] 的 total counts。"""
    if celltypes_to_check is None:
        celltypes_to_check = ['AC', 'Microglia', 'RPE']
    for ct in celltypes_to_check:
        smooth_file = os.path.join(h5ad_dir, f"{ct}_clean5_smoothed.h5ad")
        if not os.path.exists(smooth_file):
            print(f"{ct}: {smooth_file} not found, skip")
            continue
        adata = sc.read_h5ad(smooth_file)
        X_smooth = adata.X
        X_raw = adata.layers['raw']
        total_smooth = np.array(X_smooth.sum(axis=1)).flatten() if sparse.issparse(X_smooth) else X_smooth.sum(axis=1)
        total_raw = np.array(X_raw.sum(axis=1)).flatten() if sparse.issparse(X_raw) else X_raw.sum(axis=1)
        diff = total_smooth - total_raw
        n_cells = adata.n_obs
        n_ge = (diff >= 0).sum()
        n_lt = (diff < 0).sum()
        print(f"\n=== {ct} (n_cells={n_cells}) ===")
        print(f"  total_smooth >= total_raw: {n_ge}/{n_cells} cells")
        if n_lt > 0:
            print(f"  total_smooth <  total_raw: {n_lt}/{n_cells} cells  <-- 异常")
        print(f"  total_raw   : min={total_raw.min():.0f}, max={total_raw.max():.0f}, mean={total_raw.mean():.0f}")
        print(f"  total_smooth: min={total_smooth.min():.0f}, max={total_smooth.max():.0f}, mean={total_smooth.mean():.0f}")
        print(f"  ratio (smooth/raw) mean: {(total_smooth / (total_raw + 1e-10)).mean():.2f}")
    print("\n=== check_smooth_total_counts done ===")

check_smooth_total_counts(h5ad_dir, celltypes_to_check=['AC', 'Microglia', 'RPE'])


########################################################################
######## Python: 将 smoothed adata 导出为 mtx + barcodes + genes + metadata ########
########################################################################
def export_smoothed_adata_to_mtx(h5ad_dir, output_base=None, celltypes=None):
    """将每个 celltype 的 *_clean5_smoothed.h5ad 导出为 mtx、barcodes.tsv、features.tsv、metadata.tsv，供 R 读入建 Seurat。"""
    from scipy.io import mmwrite
    if output_base is None:
        output_base = os.path.join(os.path.dirname(h5ad_dir), "smoothed_mtx_clean5")
    os.makedirs(output_base, exist_ok=True)
    if celltypes is None:
        celltypes = ['AC', 'BC', 'RGC', 'Rod', 'Cone', 'HC', 'MG', 'Microglia', 'RPE']
    for ct in celltypes:
        h5ad_file = os.path.join(h5ad_dir, f"{ct}_clean5_smoothed.h5ad")
        if not os.path.exists(h5ad_file):
            print(f"{ct}: {h5ad_file} not found, skip")
            continue
        out_dir = os.path.join(output_base, f"{ct}_clean5_smoothed")
        os.makedirs(out_dir, exist_ok=True)
        adata = sc.read_h5ad(h5ad_file)
        X = adata.X
        X = sparse.csr_matrix(X) if not sparse.issparse(X) else X.tocsr()
        mmwrite(os.path.join(out_dir, "matrix.mtx"), X.T)   # genes x cells
        pd.DataFrame(adata.obs_names, columns=["barcode"]).to_csv(os.path.join(out_dir, "barcodes.tsv"), sep="\t", index=False, header=False)
        pd.DataFrame(adata.var_names, columns=["gene"]).to_csv(os.path.join(out_dir, "features.tsv"), sep="\t", index=False, header=False)
        adata.obs.to_csv(os.path.join(out_dir, "metadata.tsv"), sep="\t", index=True)
        print(f"  {ct}: {out_dir}  (cells={adata.n_obs}, genes={adata.n_vars})")
    print(f"\nExport done. MTX 目录: {output_base}")

export_smoothed_adata_to_mtx(h5ad_dir, output_base="/projects/hmz-aging/Human_New_Datasets/smoothed_mtx_clean5")


########################################################################
######## R: 读 mtx 建 Seurat 并保存 Rds（供后续 elastic regression）########
########################################################################
conda activate ArchR
R

library(Matrix)
library(Seurat)
library(parallel)
mtx_base <- "/projects/hmz-aging/Human_New_Datasets/smoothed_mtx_clean5"
seurat_out_dir <- "/projects/hmz-aging/Human_New_Datasets/smoothed_seurat_clean5"
n_cores <- 5
dir.create(seurat_out_dir, showWarnings = FALSE, recursive = TRUE)

load_mtx_to_seurat <- function(mtx_dir) {
  mat <- readMM(file.path(mtx_dir, "matrix.mtx"))
  mat <- as(mat, "CsparseMatrix")
  genes <- read.table(file.path(mtx_dir, "features.tsv"), header = FALSE, stringsAsFactors = FALSE)$V1
  barcodes <- read.table(file.path(mtx_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)$V1
  rownames(mat) <- genes
  colnames(mat) <- barcodes
  meta <- read.table(file.path(mtx_dir, "metadata.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(meta) <- meta[[1]]
  CreateSeuratObject(counts = mat, meta.data = meta)
}

process_one_celltype <- function(ct) {
  mtx_dir <- file.path(mtx_base, paste0(ct, "_clean5_smoothed"))
  if (!dir.exists(mtx_dir)) { message(ct, ": mtx dir not found"); return(invisible(NULL)) }
  message("Reading: ", ct)
  seu <- load_mtx_to_seurat(mtx_dir)
  out_file <- file.path(seurat_out_dir, paste0(ct, "_clean5_smoothed.rds"))
  saveRDS(seu, out_file)
  message("  Saved: ", out_file, " (cells=", ncol(seu), ", genes=", nrow(seu), ")")
  invisible(out_file)
}

celltypes <- c("AC", "BC", "RGC", "Rod", "Cone", "HC", "MG", "Microglia", "RPE")
if (n_cores <= 1) {
  for (ct in celltypes) process_one_celltype(ct)
} else {
  mclapply(celltypes, process_one_celltype, mc.cores = n_cores)
}
message("=== All celltypes: mtx -> Seurat rds done ===")


####
####
####

/projects/hmz-aging/Human_New_Datasets/smoothed_seurat_clean5/<CT>_clean5_smoothed.rds
AC_clean5_smoothed.rds
BC_clean5_smoothed.rds
RGC_clean5_smoothed.rds
Rod_clean5_smoothed.rds
Cone_clean5_smoothed.rds
HC_clean5_smoothed.rds
MG_clean5_smoothed.rds
Microglia_clean5_smoothed.rds
RPE_clean5_smoothed.rds

########################################################################
######## Human Leave-One-Donor-Out (LODO)：折内计算 IPF 权重 ########
######## 平衡 donor + 10岁 age_interval + sex ########
######## interval 之间均衡优先：某 interval 只有单一 sex 时，该 sex 拿满整档权重 ########
########################################################################
conda activate ArchR
R

library(Matrix)
library(Seurat)
library(glmnet)
library(parallel)

setwd("/projects/hmz-aging/Human_New_Datasets")
source("Aging_Clocks_Functions.R")   # 复用 age interval / 老年 donor 下采样等辅助函数


######################################################################
######## 折内 IPF 权重：平衡 donor / age_interval / sex ########
######## 单性别 interval：该性别拿满 target_per_interval（interval 之间均衡优先）########
######################################################################
compute_balance_weights_ipf <- function(donor, sex, aint,
                                         target_per_donor = 1,
                                         target_per_interval = 1,
                                         max_iter = 100, tol = 1e-6) {
    donor <- as.character(donor); sex <- as.character(sex); aint <- as.character(aint)
    n <- length(donor)
    w <- rep(1.0, n)
    valid <- !is.na(donor) & donor != "" & !is.na(sex) & sex != "" & !is.na(aint) & aint != ""
    w[!valid] <- 0
    uniq_aint <- unique(aint[valid])
    for (iter in 1:max_iter) {
        # (a) 每个 donor 权重和 -> target_per_donor
        s_donor <- tapply(w, donor, sum, na.rm = TRUE)
        for (d in names(s_donor)) if (!is.na(s_donor[d]) && s_donor[d] > 0) w[donor == d] <- w[donor == d] * (target_per_donor / s_donor[d])
        # (b) 每个 interval 内：按“实际存在的 sex 数”平分 target_per_interval
        for (iv in uniq_aint) {
            idx <- aint == iv
            s_sex <- tapply(w[idx], sex[idx], sum, na.rm = TRUE)
            s_sex <- s_sex[!is.na(s_sex) & s_sex > 0]
            if (length(s_sex) == 0) next
            share <- target_per_interval / length(s_sex)   # 单性别 -> 全额 1.0；双性别 -> 各 0.5
            for (s in names(s_sex)) w[idx & sex == s] <- w[idx & sex == s] * (share / s_sex[s])
        }
        # 收敛检查：donor 边际 ~ target_per_donor 且 interval 边际 ~ target_per_interval
        sd_new <- tapply(w, donor, sum, na.rm = TRUE)
        si_new <- tapply(w, aint, sum, na.rm = TRUE)
        si_new <- si_new[names(si_new) %in% uniq_aint]
        if (length(si_new) > 0 &&
            max(abs(sd_new - target_per_donor), na.rm = TRUE) < tol &&
            max(abs(si_new - target_per_interval), na.rm = TRUE) < tol) break
    }
    w
}


######################################################################
######## 折内加权的 Elastic Net LODO ########
######## 每折只用训练集细胞现算权重，留出的 donor 完全不参与权重估计（无泄漏）########
######################################################################
elastic_net_LODO_perfold <- function(seu, donor_col = "donor_id", sex_col = "sex",
                                      age_col = "age", age_interval_col = "age_interval_10yr",
                                      alphas = c(0.1, 0.5, 1), ncores = 5, scale.factor = 1e4) {
    seu <- seurat_add_age_numeric_and_interval(seu, age_col = age_col)
    meta <- seu@meta.data
    if (any(is.na(meta$Age_numeric))) warning("Some ages could not be parsed to numeric")

    # 基因过滤 + depth normalization + log1p（与原 elastic_net_LODO 一致）
    counts <- GetAssayData(seu, layer = "counts")
    if (is.null(counts)) counts <- GetAssayData(seu, layer = "counts", assay = "RNA")
    mt_genes <- grep("^MT-|^mt-", rownames(counts), value = TRUE)
    counts <- counts[!rownames(counts) %in% mt_genes, , drop = FALSE]
    gene_pct <- Matrix::rowSums(counts > 0) / ncol(counts)
    counts <- counts[gene_pct >= 0.01, , drop = FALSE]
    counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
    col_sums <- Matrix::colSums(counts)
    counts_norm <- sweep(counts, 2, col_sums, FUN = "/") * scale.factor
    X <- Matrix::t(log1p(counts_norm))
    y <- meta$Age_numeric
    donors <- as.character(meta[[donor_col]])
    sexes  <- as.character(meta[[sex_col]])
    aints  <- as.character(meta[[age_interval_col]])
    cell_names <- rownames(meta)
    all_donors <- unique(donors)
    nfolds_inner <- 5

    process_one_donor <- function(test_donor) {
        idx_test  <- which(donors == test_donor)
        idx_train <- which(donors != test_donor)

        # —— 折内 IPF：只用训练集细胞现算权重 ——
        weights_train <- compute_balance_weights_ipf(donors[idx_train], sexes[idx_train], aints[idx_train])

        train_donors <- unique(donors[idx_train])
        n_folds <- min(nfolds_inner, length(train_donors))
        donor_to_fold <- setNames(sample(rep(seq_len(n_folds), length.out = length(train_donors))), train_donors)
        foldid <- donor_to_fold[donors[idx_train]]

        X_train <- X[idx_train, , drop = FALSE]; y_train <- y[idx_train]
        X_test  <- X[idx_test, , drop = FALSE]
        X_train_scaled <- scale(X_train)
        center <- attr(X_train_scaled, "scaled:center"); scale_val <- attr(X_train_scaled, "scaled:scale")
        scale_val[scale_val == 0] <- 1
        X_test_scaled <- sweep(sweep(X_test, 2, center, "-"), 2, scale_val, "/")

        best_cvm <- Inf; best_fit <- NULL; best_alpha <- NA_real_
        for (a in alphas) {
            fit <- tryCatch(cv.glmnet(X_train_scaled, y_train, alpha = a, foldid = foldid,
                                      weights = weights_train, standardize = FALSE), error = function(e) NULL)
            if (is.null(fit)) next
            cvm_min <- min(fit$cvm, na.rm = TRUE)
            if (cvm_min < best_cvm) { best_cvm <- cvm_min; best_fit <- fit; best_alpha <- a }
        }
        if (is.null(best_fit)) return(list(pred = NULL, alpha = NA_real_, lambda = NA_real_, ncoef = NA_integer_, coef = NULL))
        y_pred <- as.numeric(predict(best_fit, newx = X_test_scaled, s = "lambda.min"))
        coef_vec <- as.numeric(coef(best_fit, s = "lambda.min"))
        coef_df <- data.frame(gene = rownames(coef(best_fit, s = "lambda.min")), coef = coef_vec, stringsAsFactors = FALSE)
        coef_df <- coef_df[coef_df$coef != 0 & coef_df$gene != "(Intercept)", ]
        list(
            pred = data.frame(cell = cell_names[idx_test], actual = y[idx_test], predicted = y_pred,
                              donor = test_donor, alpha_chosen = best_alpha, stringsAsFactors = FALSE),
            alpha = best_alpha, lambda = best_fit$lambda.min, ncoef = sum(coef_vec != 0), coef = coef_df
        )
    }

    n_donors <- length(all_donors)
    if (ncores <= 1) {
        res_list <- vector("list", n_donors)
        pb <- utils::txtProgressBar(0, n_donors, style = 3)
        for (i in seq_len(n_donors)) { res_list[[i]] <- process_one_donor(all_donors[i]); utils::setTxtProgressBar(pb, i) }
        close(pb)
    } else if (requireNamespace("pbmcapply", quietly = TRUE)) {
        res_list <- pbmcapply::pbmclapply(all_donors, process_one_donor, mc.cores = ncores)
    } else {
        res_list <- mclapply(all_donors, process_one_donor, mc.cores = ncores)
    }

    pred_df <- do.call(rbind, lapply(res_list, function(r) r$pred))
    if (is.null(pred_df) || nrow(pred_df) == 0) return(list(predictions = NULL, metrics = NULL, model_params = NULL, coefs = NULL))
    alpha_vec  <- sapply(res_list, function(r) if (is.null(r$pred)) NA else r$alpha)
    lambda_vec <- sapply(res_list, function(r) if (is.null(r$pred)) NA else r$lambda)
    ncoef_vec  <- sapply(res_list, function(r) if (is.null(r$pred)) NA else r$ncoef)
    coef_list <- lapply(seq_along(res_list), function(i){
        if (!is.null(res_list[[i]]$coef) && nrow(res_list[[i]]$coef) > 0) { df <- res_list[[i]]$coef; df$donor <- all_donors[i]; df } else NULL
    })
    coef_df <- do.call(rbind, coef_list[!sapply(coef_list, is.null)])
    rmse <- sqrt(mean((pred_df$actual - pred_df$predicted)^2))
    mae  <- mean(abs(pred_df$actual - pred_df$predicted))
    pcc  <- cor(pred_df$actual, pred_df$predicted, use = "pairwise.complete.obs")
    r2   <- 1 - sum((pred_df$actual - pred_df$predicted)^2) / sum((pred_df$actual - mean(pred_df$actual))^2)
    list(predictions = pred_df,
         metrics = c(RMSE = rmse, MAE = mae, PCC = pcc, R2 = r2),
         model_params = data.frame(donor_left_out = all_donors, alpha_chosen = alpha_vec,
                                   lambda_min = lambda_vec, n_nonzero_coef = ncoef_vec, stringsAsFactors = FALSE),
         coefs = coef_df)
}


######################################################################
######## 主循环：对每个 celltype 跑折内加权 LODO 并保存 ########
######################################################################
seurat_dir <- "/projects/hmz-aging/Human_New_Datasets/smoothed_seurat_clean5"
out_dir    <- "/projects/hmz-aging/Human_New_Datasets/LODO_results_clean5_2026"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
celltypes <- c("AC", "BC", "RGC", "Rod", "Cone", "HC", "MG", "Microglia", "RPE")

for (ct in celltypes) {
    rds_file <- file.path(seurat_dir, paste0(ct, "_clean5_smoothed.rds"))
    if (!file.exists(rds_file)) { message(ct, ": rds not found, skip -> ", rds_file); next }
    message("\n======== ", ct, " ========")
    seu <- readRDS(rds_file)

    # 方案 A：不做老年段下采样，全量细胞交给折内 IPF 权重去偏（确定性、可复现、不丢数据）
    seu <- seurat_add_age_numeric_and_interval(seu, age_col = "age")
    # seu <- seurat_downsample_old_donors_10yr(seu, donor_col = "donor_id",
    #                                          age_interval_col = "age_interval_10yr", sex_col = "sex")

    res <- tryCatch(
        elastic_net_LODO_perfold(seu, donor_col = "donor_id", sex_col = "sex",
                                 age_col = "age", alphas = c(0.1, 0.5, 1), ncores = 5),
        error = function(e) { message(ct, " 出错: ", e$message); NULL })
    if (is.null(res)) next

    saveRDS(res, file.path(out_dir, paste0(ct, "_LODO_results_2026.rds")))
    write.csv(res$predictions, file.path(out_dir, paste0(ct, "_LODO_predictions_2026.csv")), row.names = FALSE)
    message(sprintf("  metrics: RMSE=%.3f MAE=%.3f PCC=%.3f R2=%.3f",
                    res$metrics["RMSE"], res$metrics["MAE"], res$metrics["PCC"], res$metrics["R2"]))
    message("  saved -> ", file.path(out_dir, paste0(ct, "_LODO_results_2026.rds")))
}
message("\n=== Human LODO（折内加权 / interval 之间均衡优先）完成 ===")


###
### ###
### ###
###


#####
##### 检查Human 的 LOOCV 和 Final 有没有跑完 ####
#####

# 输出目录
lodo_dir  <- "/projects/hmz-aging/Human_New_Datasets/LODO_results_clean5_2026"
final_dir <- "/projects/hmz-aging/Human_New_Datasets/FinalClock_clean5_2026"

lodo_celltypes  <- c("AC", "BC", "RGC", "Rod", "Cone", "HC", "MG")  # LODO 跑全部 9 类
final_celltypes <- c("AC", "BC", "RGC", "Rod", "Cone", "HC", "MG")                       # Final 只 major 7 类

# 小工具：返回文件是否存在 + 大小(MB) + 修改时间
file_stat <- function(fp) {
    if (!file.exists(fp)) return(data.frame(exists = FALSE, size_MB = NA_real_, mtime = as.POSIXct(NA)))
    info <- file.info(fp)
    data.frame(exists = TRUE, size_MB = round(info$size / 1024^2, 2), mtime = info$mtime)
}

# ---- 1) LODO 结果检查 ----
cat("\n========== Human LODO 结果检查 ==========\n")
lodo_tbl <- data.frame()
for (ct in lodo_celltypes) {
    rds <- file.path(lodo_dir, paste0(ct, "_LODO_results_2026.rds"))
    csv <- file.path(lodo_dir, paste0(ct, "_LODO_predictions_2026.csv"))
    s_rds <- file_stat(rds); s_csv <- file_stat(csv)
    cat(sprintf("%-10s  rds:%-5s (%6s MB)  csv:%-5s | %s\n",
                ct, s_rds$exists, ifelse(is.na(s_rds$size_MB), "-", s_rds$size_MB),
                s_csv$exists, ifelse(is.na(s_rds$mtime), "", format(s_rds$mtime, "%Y-%m-%d %H:%M"))))
    lodo_tbl <- rbind(lodo_tbl, data.frame(celltype = ct, rds = s_rds$exists,
                       csv = s_csv$exists, size_MB = s_rds$size_MB, mtime = s_rds$mtime))
}
cat(sprintf("LODO 完成: %d / %d\n", sum(lodo_tbl$rds), length(lodo_celltypes)))

# ---- 2) Final 模型检查 ----
cat("\n========== Human Final 模型检查（major 7 类）==========\n")
final_tbl <- data.frame()
for (ct in final_celltypes) {
    rds <- file.path(final_dir, paste0(ct, "_final_aging_clock_model_2026.rds"))
    csv <- file.path(final_dir, paste0(ct, "_final_model_coefs_2026.csv"))
    s_rds <- file_stat(rds); s_csv <- file_stat(csv)
    cat(sprintf("%-10s  rds:%-5s (%6s MB)  coefs:%-5s | %s\n",
                ct, s_rds$exists, ifelse(is.na(s_rds$size_MB), "-", s_rds$size_MB),
                s_csv$exists, ifelse(is.na(s_rds$mtime), "", format(s_rds$mtime, "%Y-%m-%d %H:%M"))))
    final_tbl <- rbind(final_tbl, data.frame(celltype = ct, rds = s_rds$exists,
                        coefs = s_csv$exists, size_MB = s_rds$size_MB, mtime = s_rds$mtime))
}
cat(sprintf("Final 完成: %d / %d\n", sum(final_tbl$rds), length(final_celltypes)))


############# Next ###############
#############




######################################################################
######## 读入 LOOCV 结果，画 violin plot（每个 celltype 一个分面）########
######## x = 真实年龄 interval, y = 预测年龄, 按年龄 viridis 配色, 内嵌白色箱线图 ####
######################################################################
library(ggplot2)
library(viridis)

base_dir <- "/projects/hmz-aging/Human_New_Datasets/LODO_results_clean5_2026"
violin_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")   # major 7 类，单行顺序

# 10 年年龄段（与训练时 age_interval_10yr 一致）
age_breaks <- seq(0, 90, by = 10)
age_labels <- c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90")

# 读入每个 celltype 的 LODO cell 级预测
pred_all <- list()
for (ct in violin_celltypes) {
    fp <- file.path(base_dir, ct, paste0(ct, "_LODO_results_2026.rds"))
    if (!file.exists(fp)) {
        # 结果直接存在 base_dir 下（非子文件夹）
        fp <- file.path(base_dir, paste0(ct, "_LODO_results_2026.rds"))
    }
    if (!file.exists(fp)) { message(ct, ": LODO 结果不存在，跳过 -> ", fp); next }
    res <- readRDS(fp)
    df <- res$predictions               # cell, actual(Age_numeric), predicted, donor
    df$celltype <- ct
    pred_all[[ct]] <- df[, c("celltype", "actual", "predicted")]
}
pred_df <- do.call(rbind, pred_all)

# celltype 顺序（只保留实际存在的）
present_ct <- intersect(violin_celltypes, unique(pred_df$celltype))
pred_df$celltype <- factor(pred_df$celltype, levels = present_ct)

# 连续真实年龄 -> 10 年 interval 作为 X 轴
pred_df$age_interval <- cut(pred_df$actual, breaks = age_breaks, labels = age_labels,
                            right = FALSE, include.lowest = TRUE)
pred_df <- pred_df[!is.na(pred_df$age_interval), ]
pred_df$age_interval <- factor(pred_df$age_interval,
                               levels = age_labels[age_labels %in% unique(as.character(pred_df$age_interval))])

p_violin <- ggplot(pred_df, aes(x = age_interval, y = predicted, fill = age_interval)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.9) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", lwd = 0.3) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "D") +
    facet_wrap(~ celltype, nrow = 1, scales = "free_x") +
    coord_cartesian(ylim = c(0, 95)) +
    labs(x = "Ages", y = "Predicted age") +
    theme_classic(base_size = 12) +
    theme(
        strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text       = element_text(face = "bold", size = 15),
        axis.text.x      = element_text(size = 16, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y      = element_text(size = 13, color = "black"),
        axis.title       = element_text(size = 18, face = "bold"),
        legend.position  = "none",
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

ncol_ct <- length(present_ct)
ggsave("Human_LODO_violin_by_celltype_2026.png", p_violin,
       width = max(4, 2.6 * ncol_ct), height = 3.6, dpi = 300)
message("Saved: Human_LODO_violin_by_celltype_2026.png")


######################################################################
######## 按画图顺序，逐个打印每个 celltype 的 RMSE / MAE / R2 / Pearson ########
######## （cell 级，基于连续真实年龄 vs 预测年龄现算）                     ########
######################################################################
metric_tbl <- data.frame()
for (ct in present_ct) {
    d <- pred_df[pred_df$celltype == ct & is.finite(pred_df$predicted), ]
    if (nrow(d) == 0) { message(ct, ": 无预测数据，跳过"); next }
    a <- d$actual; p <- d$predicted
    rmse    <- sqrt(mean((a - p)^2))
    mae     <- mean(abs(a - p))
    r2      <- 1 - sum((a - p)^2) / sum((a - mean(a))^2)
    pearson <- cor(a, p)
    cat(sprintf("%-5s  RMSE=%.3f  MAE=%.3f  R2=%.3f  Pearson=%.3f\n",
                ct, rmse, mae, r2, pearson))
    metric_tbl <- rbind(metric_tbl,
        data.frame(celltype = ct, RMSE = rmse, MAE = mae, R2 = r2, Pearson = pearson))
}
write.csv(metric_tbl, "Human_LODO_metrics_by_celltype_2026.csv", row.names = FALSE)
message("Saved: Human_LODO_metrics_by_celltype_2026.csv")


######################################################################
######## 诊断：核对 RMSE 是按连续年龄算的，且与存储 res$metrics 一致 ########
######## 同时给出 cell 级 vs donor(样本)级 RMSE 对比                    ########
######################################################################
cat("\n========== Human LODO RMSE 诊断 ==========\n")
for (ct in violin_celltypes) {
    fp <- file.path(base_dir, ct, paste0(ct, "_LODO_results_2026.rds"))
    if (!file.exists(fp)) fp <- file.path(base_dir, paste0(ct, "_LODO_results_2026.rds"))
    if (!file.exists(fp)) { message(ct, ": 结果不存在，跳过"); next }
    res <- readRDS(fp)
    d <- res$predictions
    d <- d[is.finite(d$actual) & is.finite(d$predicted), ]
    if (nrow(d) == 0) { message(ct, ": 无有效预测"); next }

    # 1) actual/predicted 的真实取值范围（确认 actual 是岁数而非 interval 索引）
    cat(sprintf("\n--- %s ---\n", ct))
    cat(sprintf("  actual   : n=%d, range=[%.1f, %.1f], unique=%d 个值\n",
                nrow(d), min(d$actual), max(d$actual), length(unique(d$actual))))
    cat(sprintf("  predicted: range=[%.2f, %.2f], mean=%.2f\n",
                min(d$predicted), max(d$predicted), mean(d$predicted)))

    # 2) cell 级 RMSE（连续年龄）——我重算 vs LODO 存储值
    rmse_cell <- sqrt(mean((d$actual - d$predicted)^2))
    stored <- if (!is.null(res$metrics)) as.numeric(res$metrics["RMSE"]) else NA_real_
    cat(sprintf("  cell 级 RMSE(重算)=%.3f | 存储 res$metrics RMSE=%.3f | 一致=%s\n",
                rmse_cell, stored,
                ifelse(is.na(stored), "无存储值", ifelse(abs(rmse_cell - stored) < 1e-3, "是", "否!!"))))

    # 3) donor(样本)级 RMSE：先对每个 donor 取预测均值，再和该 donor 真实年龄比
    if ("donor" %in% colnames(d)) {
        agg <- aggregate(cbind(actual, predicted) ~ donor, data = d, FUN = mean)
        rmse_donor <- sqrt(mean((agg$actual - agg$predicted)^2))
        cat(sprintf("  donor 级 RMSE=%.3f (n_donor=%d)\n", rmse_donor, nrow(agg)))
    }

    # 4) 基线参照：把所有细胞都预测成平均年龄时的 RMSE（= 年龄的 sd）
    cat(sprintf("  参照: 若全预测为平均年龄, RMSE=%.3f (即年龄 sd)\n",
                sqrt(mean((d$actual - mean(d$actual))^2))))
}






#######
####### 把三个物种的 final aging clocks 整理成1个 excel 表，包含多个 tabs 每个tab是一个 celltype 的 wetghts 哥 gene ####
#######
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# 三个物种 final clock 系数(csv: gene, coef)的路径构造函数
human_dir <- "/projects/hmz-aging/Human_New_Datasets/FinalClock_clean5_2026"
mouse_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
fish_dir  <- "/projects/hmz-aging/Zebrafish_New_Datasets"

coef_path <- list(
    Human     = function(ct) file.path(human_dir, paste0(ct, "_final_model_coefs_2026.csv")),
    Mouse     = function(ct) file.path(mouse_dir, ct, paste0(ct, "_final_model_coefs_2026.csv")),
    Zebrafish = function(ct) file.path(fish_dir,  ct, paste0(ct, "_final_model_coefs_2026.csv"))
)

excel_celltypes <- c("MG", "AC", "BC", "RGC", "Cone", "HC", "Rod")   # 三物种共有的 major 7 类
species_order   <- c("Human", "Mouse", "Zebrafish")

read_coef <- function(fp, species) {
    if (!file.exists(fp)) { message("  缺失: ", fp); return(NULL) }
    df <- read.csv(fp, stringsAsFactors = FALSE)
    if (!all(c("gene", "coef") %in% colnames(df))) { message("  列名异常: ", fp); return(NULL) }
    df <- df[df$coef != 0 & df$gene != "(Intercept)", c("gene", "coef")]
    if (nrow(df) == 0) return(NULL)
    df$species <- species
    df[order(-abs(df$coef)), c("species", "gene", "coef")]   # 按 |coef| 降序
}

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", halign = "center", fgFill = "#D9E1F2", border = "TopBottomLeftRight")

# 第一个 tab 预留为 Contents 目录（内容循环后再填）
addWorksheet(wb, "Contents")

# 每个物种 × 每个 celltype 单独一个 tab（共 species × celltype 个），tab 名如 Human_MG
contents <- data.frame()
for (sp in species_order) {
    for (ct in excel_celltypes) {
        df <- read_coef(coef_path[[sp]](ct), sp)
        sheet <- paste0(sp, "_", ct)
        if (is.null(df) || nrow(df) == 0) { message(sheet, ": 系数缺失，跳过该 tab"); next }
        df <- df[, c("gene", "coef")]   # 单物种 tab 不再需要 species 列
        addWorksheet(wb, sheet)
        writeData(wb, sheet, df, headerStyle = hs)
        freezePane(wb, sheet, firstRow = TRUE)
        setColWidths(wb, sheet, cols = 1:2, widths = c(18, 12))
        cat(sprintf("%-16s tab: %d genes\n", sheet, nrow(df)))
        contents <- rbind(contents, data.frame(
            Tab = sheet, Species = sp, CellType = ct, N_genes = nrow(df),
            Description = sprintf("%s %s final aging clock weights (gene, coef)", sp, ct),
            stringsAsFactors = FALSE))
    }
}

# 填 Contents：标题 + 目录表（Tab 名做成可点击的内部超链接）
writeData(wb, "Contents", "Final Aging Clock Weights — 3 Species (2026)",
          startRow = 1, startCol = 1)
addStyle(wb, "Contents", createStyle(textDecoration = "bold", fontSize = 14),
         rows = 1, cols = 1)
writeData(wb, "Contents", contents, startRow = 3, startCol = 1, headerStyle = hs)
# Tab 列写成内部超链接，点击跳到对应工作表
for (i in seq_len(nrow(contents))) {
    link <- sprintf("internal:'%s'!A1", contents$Tab[i])
    class(link) <- "hyperlink"
    writeData(wb, "Contents", x = link, startRow = 3 + i, startCol = 1)
}
freezePane(wb, "Contents", firstActiveRow = 4)
setColWidths(wb, "Contents", cols = 1:5, widths = c(16, 12, 10, 10, 55))

excel_out <- "Final_AgingClocks_weights_3species_2026.xlsx"
saveWorkbook(wb, excel_out, overwrite = TRUE)
message("Saved: ", excel_out, " (", length(names(wb)), " tabs, 含 Contents)")

