#!/usr/bin/env Rscript
# =============================================================================
# 2026final_1e5_sliceCleaned_OSKclean | SCT-kNN aging clocks → predict OSK
#
# Copied from train_xenium_SCT_knn_predict_OSK_2026final_1e5_sliceCleaned.R.
# The previous script is unchanged. This version uses the slice-cleaned OSK object for prediction.
#
# Aging Seurat (training, same as the previous script):
#   .../R1_4_seurat_Merged_dim_RNA_2024_sliceCleaned_colored.rds
# OSK Seurat (prediction, replaced in this version):
#   /projects/hmz-aging/Xenium5k_Mouse/OSK_New/OSK_Spatial_slice_clean/OSK_obj_cl_2026_sliceCleaned_OSKR3R4.rds
# For both Aging and OSK, SCTransform is run on RNA if SCT is absent.
#
# Clocks can reuse the previous sliceCleaned training results; OSK predictions go to a separate directory.
#
# Example:
#   # Prediction only: with existing Xmodel + OSK knn, predict and plot directly (no retraining or re-smoothing)
#   Rscript train_xenium_SCT_knn_predict_OSK_2026final_1e5_sliceCleaned_OSKclean.R --mode predict
#
#   # Full pipeline (training + prediction)
#   Rscript train_xenium_SCT_knn_predict_OSK_2026final_1e5_sliceCleaned_OSKclean.R --mode all --n_cores 20
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(FNN)
  library(glmnet)
  library(parallel)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}
has_flag <- function(flag) !is.na(match(flag, args))

# ---- paths ----
osk_root <- "/projects/hmz-aging/Xenium5k_Mouse/OSK_New"
proximity_root <- "/projects/hmz-aging/ProximityEffect"

xenium_seurat_cands <- c(
  "/path/to/Proximity_final/R1_4_seurat_Merged_dim_RNA_2024_sliceCleaned_colored.rds",
  file.path(
    proximity_root, "Spatial_slice_clean",
    "R1_4_seurat_Merged_dim_RNA_2024_sliceCleaned_colored.rds"
  )
)

mode <- tolower(get_arg("--mode", "all")) # train | predict | all
xenium_seurat <- get_arg("--xenium_seurat", NA_character_)
if (is.na(xenium_seurat) || !nzchar(xenium_seurat)) {
  hit <- xenium_seurat_cands[file.exists(xenium_seurat_cands)]
  xenium_seurat <- if (length(hit)) hit[[1]] else xenium_seurat_cands[[1]]
}
osk_seurat_cands <- c(
  file.path(
    osk_root, "OSK_Spatial_slice_clean",
    "OSK_obj_cl_2026_sliceCleaned_OSKR3R4.rds"
  ),
  file.path(
    osk_root, "OSK_Spatial_slice_clean",
    "OSK_obj_cl_2026_sliceCleaned_colored.rds"
  )
)
osk_seurat <- get_arg("--osk_seurat", NA_character_)
if (is.na(osk_seurat) || !nzchar(osk_seurat)) {
  hit_o <- osk_seurat_cands[file.exists(osk_seurat_cands)]
  osk_seurat <- if (length(hit_o)) hit_o[[1]] else osk_seurat_cands[[1]]
}
clock_dir <- get_arg(
  "--clock_dir",
  file.path(proximity_root, "Aging_Clocks_Xenium_SCT_knn_2026final_1e5_sliceCleaned")
)
out_dir <- get_arg(
  "--out_dir",
  file.path(osk_root, "OSK_predict_Xenium_SCT_knn_2026final_1e5_sliceCleaned_OSKclean")
)
xenium_knn_rds <- get_arg(
  "--xenium_knn_rds",
  file.path(clock_dir, "Xenium5k_Mouse_Obj_pseudobulk_list_2026final_1e5_sliceCleaned.rds")
)
osk_knn_rds <- get_arg(
  "--osk_knn_rds",
  file.path(out_dir, "OSK_obj_pseudobulk_list_SCT_knn_2026final_1e5_sliceCleaned_OSKclean.rds")
)

celltypes <- strsplit(
  get_arg("--celltypes", "MG,Rod,Cone,RGC,AC,HC,RPE,BC"),
  ",",
  fixed = TRUE
)[[1]]
# 2026final variant: depth-norm = 1e5 after kNN smoothing
scale_factor <- as.numeric(get_arg("--scale_factor", "1e5"))
n_cores <- as.integer(get_arg("--n_cores", "20"))
train_frac <- as.numeric(get_arg("--train_frac", "0.7"))
seed <- as.integer(get_arg("--seed", "42"))
ct_col_xenium <- get_arg("--ct_col_xenium", "celltype")
ct_col_osk <- get_arg("--ct_col_osk", "celltype")
sample_col_osk <- get_arg("--sample_col_osk", "sample")
condition_map_str <- get_arg(
  "--condition_map",
  "OSKR3=Ctrl,OSKR4=OSK,R3=Ctrl,R4=OSK"
)
ct_order <- strsplit(
  get_arg("--ct_order", "MG,RGC,AC,HC,Rod,Cone,BC,RPE"),
  ",",
  fixed = TRUE
)[[1]]

skip_xenium_knn <- has_flag("--skip_xenium_knn")
skip_osk_knn <- has_flag("--skip_osk_knn")
skip_train <- has_flag("--skip_train")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clock_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
  flush.console()
}

parse_condition_map <- function(s) {
  parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  kv <- strsplit(parts, "=", fixed = TRUE)
  setNames(
    vapply(kv, function(x) if (length(x) < 2) x[[1]] else paste(x[-1], collapse = "="), character(1)),
    vapply(kv, function(x) trimws(x[[1]]), character(1))
  )
}
cond_map <- parse_condition_map(condition_map_str)

ct_k_npcs <- function(ct) {
  # Consistent with OSK_clocks_New_2026 / Xenium_Aging_Clocks.R
  if (ct %in% c("Rod", "AC", "MG", "BC", "Cone", "RGC")) {
    list(k = 15L, npcs = 25L)
  } else if (ct %in% c("HC", "RPE")) {
    list(k = 2L, npcs = 5L)
  } else {
    list(k = 2L, npcs = 5L)
  }
}

ensure_rds <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) path else paste0(path, ".rds")
}

load_seurat_flexible <- function(path) {
  cands <- unique(c(path, ensure_rds(path), sub("\\.rds$", "", path)))
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("Seurat not found: ", path)
  log_msg("Loading: ", hit)
  obj <- readRDS(hit)
  if (!inherits(obj, "Seurat")) stop("Not Seurat: ", hit)
  obj
}

ensure_sct <- function(seu) {
  if (exists("JoinLayers", mode = "function") && "RNA" %in% Assays(seu)) {
    tryCatch({
      DefaultAssay(seu) <- "RNA"
      seu <- JoinLayers(seu, assay = "RNA")
      log_msg("JoinLayers(RNA) ok")
    }, error = function(e) {
      log_msg("JoinLayers skip: ", conditionMessage(e))
    })
  }
  if ("SCT" %in% Assays(seu)) {
    log_msg("SCT already present")
    return(seu)
  }
  if (!"RNA" %in% Assays(seu)) {
    stop("No SCT and no RNA. Assays: ", paste(Assays(seu), collapse = ", "))
  }
  log_msg("sliceCleaned has no SCT; running SCTransform on RNA ...")
  DefaultAssay(seu) <- "RNA"
  seu <- SCTransform(seu, verbose = FALSE, return.only.var.genes = FALSE)
  log_msg("SCTransform done | assays=", paste(Assays(seu), collapse = ","))
  seu
}

get_counts <- function(seu, assay) {
  out <- tryCatch(
    GetAssayData(seu, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(out)) {
    out <- tryCatch(
      GetAssayData(seu, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(out)) {
    a <- seu[[assay]]
    if (!is.null(a@layers) && length(a@layers) > 0) {
      out <- do.call(cbind, a@layers)
      rownames(out) <- rownames(seu)
      colnames(out) <- colnames(seu)
    }
  }
  if (is.null(out) || length(out) == 0) {
    stop("Cannot extract counts from assay=", assay)
  }
  if (!inherits(out, "dgCMatrix")) {
    out <- as(as(out, "CsparseMatrix"), "dgCMatrix")
  }
  out
}

# =============================================================================
# Core helpers (OSK_clocks_New_2026)
# =============================================================================

NormalizeAndPCASCT <- function(
    seu,
    scale.factor = 1e4,
    nfeatures = 2000,
    npcs = 30
) {
  # Original: DefaultAssay unchanged; after SCTransform the default is SCT → PCA/kNN computed on SCT
  # PseudoBulk then sums SCT@counts separately
  if ("SCT" %in% Assays(seu)) {
    DefaultAssay(seu) <- "SCT"
  }
  seu <- NormalizeData(
    seu,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor,
    verbose = FALSE
  )
  nfeatures_use <- min(nfeatures, nrow(seu))
  # Reuse SCT VariableFeatures if present; otherwise recompute
  if (!length(VariableFeatures(seu))) {
    seu <- FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = nfeatures_use,
      verbose = FALSE
    )
  }
  feats <- VariableFeatures(seu)
  if (!length(feats)) feats <- rownames(seu)[seq_len(min(nfeatures_use, nrow(seu)))]
  # Original: ScaleData(features=rownames(seu)); using VF is more robust
  seu <- ScaleData(seu, features = feats, verbose = FALSE)
  n_cells <- ncol(seu)
  max_npcs <- min(length(feats), n_cells) - 1L
  if (npcs > max_npcs) npcs <- max(1L, max_npcs)
  if (npcs < 1L) stop("Too few cells for PCA: ", n_cells)
  seu <- RunPCA(
    seu,
    features = feats,
    npcs = npcs,
    verbose = FALSE
  )
  seu
}

PseudoBulkByKNN_fast_SCT <- function(
    seu,
    k = 15,
    dims = 1:30,
    assay_out = "RNA",
    reduction = "pca"
) {
  emb <- Embeddings(seu, reduction)[, dims, drop = FALSE]
  knn <- get.knn(emb, k = k)
  nn_idx <- knn$nn.index
  ncells <- nrow(nn_idx)

  # self + k neighbors
  rows <- as.integer(c(t(nn_idx), seq_len(ncells)))
  cols <- rep(seq_len(ncells), each = k + 1L)
  A <- sparseMatrix(i = rows, j = cols, x = 1, dims = c(ncells, ncells))

  if (!"SCT" %in% Assays(seu)) {
    stop("PseudoBulkByKNN_fast_SCT needs SCT assay")
  }
  counts_mat <- get_counts(seu, "SCT")
  counts_mat <- counts_mat[, rownames(emb), drop = FALSE]
  pseudo_counts <- counts_mat %*% A
  colnames(pseudo_counts) <- colnames(counts_mat)

  meta <- seu@meta.data[colnames(pseudo_counts), , drop = FALSE]
  CreateSeuratObject(
    counts = pseudo_counts,
    meta.data = meta,
    assay = assay_out,
    project = "SCT_knn"
  )
}

smooth_by_ct_sample <- function(
    seu,
    celltypes,
    ct_col,
    sample_col,
    scale.factor = 1e4
) {
  if (!"SCT" %in% Assays(seu)) {
    stop("Need SCT. Assays: ", paste(Assays(seu), collapse = ", "))
  }
  if (!ct_col %in% colnames(seu@meta.data)) stop("Missing ", ct_col)
  if (!sample_col %in% colnames(seu@meta.data)) stop("Missing ", sample_col)

  seu$._ct_ <- as.character(seu[[ct_col]][, 1])
  seu$._sample_ <- as.character(seu[[sample_col]][, 1])

  present <- intersect(celltypes, unique(seu$._ct_))
  if (!length(present)) {
    stop("No celltypes. Available: ", paste(unique(seu$._ct_), collapse = ", "))
  }
  miss <- setdiff(celltypes, present)
  if (length(miss)) log_msg("WARNING skip CT: ", paste(miss, collapse = ", "))

  out <- list()
  for (ct in present) {
    kn <- ct_k_npcs(ct)
    log_msg("  SCT-kNN CT=", ct, " k=", kn$k, " npcs=", kn$npcs)
    cells_ct <- colnames(seu)[seu$._ct_ == ct]
    seu_ct <- seu[, cells_ct]
    by_s <- SplitObject(seu_ct, split.by = "._sample_")
    # ncol may return double in some Seurat/R versions; avoid vapply integer strictness
    by_s <- by_s[vapply(by_s, function(x) as.integer(ncol(x)), integer(1)) >= (kn$k + 2L)]
    if (!length(by_s)) {
      log_msg("    skip ", ct, ": no sample with enough cells")
      next
    }

    pb_list <- lapply(names(by_s), function(sn) {
      s <- by_s[[sn]]
      log_msg("    sample=", sn, " n=", ncol(s))
      s <- NormalizeAndPCASCT(s, scale.factor = scale.factor, npcs = kn$npcs)
      dims_use <- seq_len(min(kn$npcs, ncol(Embeddings(s, "pca"))))
      PseudoBulkByKNN_fast_SCT(s, k = kn$k, dims = dims_use, assay_out = "RNA")
    })
    names(pb_list) <- names(by_s)

    if (length(pb_list) == 1L) {
      merged <- pb_list[[1]]
    } else {
      merged <- merge(pb_list[[1]], pb_list[-1], add.cell.ids = names(pb_list))
    }
    merged$celltype <- ct
    if (!"sample" %in% colnames(merged@meta.data)) {
      if ("._sample_" %in% colnames(merged@meta.data)) {
        merged$sample <- merged$._sample_
      }
    }
    out[[ct]] <- merged
  }
  out
}

# Original: Elastic_net_one_round_Xenium
Elastic_net_one_round_Xenium <- function(
    seu,
    scale.factor = 1e4,
    train_frac = 0.7,
    alpha.vec = seq(0, 1, by = 0.1),
    nfolds = 5,
    ncores = 20,
    seed = 42
) {
  mat <- get_counts(seu, "RNA") # content = SCT-smoothed counts
  total_counts <- Matrix::colSums(mat)
  total_counts[total_counts == 0] <- 1
  mat_norm <- sweep(as.matrix(mat), 2, total_counts, "/") * scale.factor
  mat_norm <- log2(mat_norm + 1)

  if (!"age" %in% colnames(seu@meta.data)) {
    if ("sample" %in% colnames(seu@meta.data)) seu$age <- seu$sample
    else stop("Need meta$age or meta$sample")
  }
  ages <- as.numeric(as.character(seu$age[colnames(mat_norm)]))
  if (anyNA(ages)) {
    stop("Ages must be numeric 5/12/52/112. Got: ",
         paste(unique(as.character(seu$age)), collapse = ", "))
  }
  names(ages) <- colnames(mat_norm)

  set.seed(seed)
  train_idx <- integer(0)
  for (a in unique(ages)) {
    idx_age <- which(ages == a)
    n_train <- max(1L, floor(length(idx_age) * train_frac))
    train_idx <- c(train_idx, sample(idx_age, n_train))
  }
  train_idx <- unique(train_idx)
  test_idx <- setdiff(seq_along(ages), train_idx)
  if (!length(test_idx)) test_idx <- train_idx

  x_train <- t(mat_norm[, train_idx, drop = FALSE])
  y_train <- ages[train_idx]
  x_test <- t(mat_norm[, test_idx, drop = FALSE])
  y_test <- ages[test_idx]

  age_tab <- table(y_train)
  w_train <- as.numeric(1 / age_tab[as.character(y_train)])
  w_train <- w_train * length(w_train) / sum(w_train)

  cv_results <- mclapply(alpha.vec, function(a) {
    cv <- cv.glmnet(
      x = x_train, y = y_train,
      alpha = a, weights = w_train,
      type.measure = "mse",
      nfolds = min(nfolds, length(unique(y_train))),
      standardize = FALSE
    )
    list(alpha = a, cvm = min(cv$cvm), fit = cv)
  }, mc.cores = ncores)

  best <- Reduce(function(x, y) if (x$cvm < y$cvm) x else y, cv_results)
  fit <- best$fit
  preds_test <- as.numeric(predict(fit, newx = x_test, s = "lambda.min"))
  test_cor <- suppressWarnings(cor(y_test, preds_test))
  test_mse <- mean((y_test - preds_test)^2)
  log_msg("    alpha=", best$alpha,
          " test_cor=", signif(test_cor, 4),
          " test_mse=", signif(test_mse, 4))

  coefs <- as.matrix(coef(fit, s = "lambda.min"))
  coef_df <- data.frame(
    gene = rownames(coefs),
    coefficient = as.vector(coefs),
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[coef_df$gene != "(Intercept)", , drop = FALSE]

  list(
    model = fit,
    best_alpha = best$alpha,
    coefficients = coef_df,
    performance = list(test_mse = test_mse, test_cor = test_cor),
    predictions = data.frame(
      sample = colnames(mat_norm)[test_idx],
      actual = y_test,
      predicted = preds_test,
      set = "test",
      stringsAsFactors = FALSE
    )
  )
}

# Original: Prepare_Matrix_And_Meta_OSK: 1e4 + log2
Prepare_Matrix_And_Meta_OSK <- function(pb_list, scale.factor = 1e4) {
  results <- lapply(names(pb_list), function(ct) {
    tmp <- pb_list[[ct]]
    if ("RNA" %in% Assays(tmp)) {
      tryCatch(tmp[["RNA"]] <- JoinLayers(tmp[["RNA"]]), error = function(e) NULL)
    }
    mat <- get_counts(tmp, "RNA")
    meta <- tmp@meta.data
    keep <- intersect(c("sample", "celltype"), colnames(meta))
    meta <- meta[, keep, drop = FALSE]
    if (!"celltype" %in% colnames(meta)) meta$celltype <- ct
    if (!"sample" %in% colnames(meta)) stop("Missing sample for CT=", ct)
    meta$cellid <- rownames(meta)
    meta$celltype <- ct

    cs <- Matrix::colSums(mat)
    cs[cs == 0] <- 1
    # Original: factor = colSums/1e4; sweep(mat, 2, factor, "/")
    mat_log2 <- log2(sweep(as.matrix(mat), 2, cs / scale.factor, "/") + 1)
    list(Matrix = mat_log2, Meta = meta)
  })
  names(results) <- names(pb_list)
  results
}

# Original: PredictAgeFromModel
PredictAgeFromModel <- function(expr_mat, test_meta, fit_res) {
  coef_df <- coef(fit_res, s = fit_res$lambda.min)
  coef_df <- data.frame(
    gene = rownames(coef_df),
    coefficient = coef_df[, 1],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  mat_genes <- intersect(coef_df$gene, rownames(expr_mat))
  if (!length(mat_genes)) stop("No overlapping genes with model")
  coef_vec <- coef_df$coefficient[match(mat_genes, coef_df$gene)]
  intercept <- coef_df$coefficient[coef_df$gene == "(Intercept)"]
  if (!length(intercept)) intercept <- 0
  sub_mat <- expr_mat[mat_genes, , drop = FALSE]
  preds <- as.numeric(crossprod(as.matrix(sub_mat), coef_vec) + intercept)
  test_meta$preds <- preds
  test_meta$predict <- preds
  test_meta
}

map_xenium_age <- function(seu, sample1_col = "sample1") {
  if (!sample1_col %in% colnames(seu@meta.data)) {
    stop("Xenium needs ", sample1_col)
  }
  recode_map <- c(R1 = "5", R2 = "12", R3 = "52", R4 = "112")
  s1 <- as.character(seu[[sample1_col]][, 1])
  if (all(s1 %in% c("5", "12", "52", "112"))) {
    age <- s1
  } else if (all(s1 %in% names(recode_map))) {
    age <- unname(recode_map[s1])
  } else {
    stop("Cannot map sample1: ", paste(unique(s1), collapse = ", "))
  }
  seu$sample <- age
  seu$age <- as.numeric(age)
  seu
}

xmodel_path <- function(ct) {
  file.path(clock_dir, paste0("Mouse_", ct, "_Xmodel_1e5_sliceCleaned.rds"))
}

# =============================================================================
# TRAIN
# =============================================================================

run_train <- function() {
  log_msg("========== TRAIN (OSK_clocks_New_2026 style) ==========")
  if (skip_train) {
    log_msg("--skip_train: reuse Xmodels in ", clock_dir)
    return(invisible(NULL))
  }

  knn_path <- ensure_rds(xenium_knn_rds)
  if (skip_xenium_knn && (file.exists(knn_path) || file.exists(xenium_knn_rds))) {
    fp <- if (file.exists(knn_path)) knn_path else xenium_knn_rds
    log_msg("Load Xenium knn: ", fp)
    pb_list <- readRDS(fp)
  } else {
    seu_x <- load_seurat_flexible(xenium_seurat)
    if (!ct_col_xenium %in% colnames(seu_x@meta.data)) {
      stop("Xenium missing ", ct_col_xenium)
    }
    seu_x <- ensure_sct(seu_x)
    seu_x <- map_xenium_age(seu_x, "sample1")
    log_msg("Xenium n=", ncol(seu_x),
            " ages=", paste(sort(unique(seu_x$age)), collapse = ","))
    pb_list <- smooth_by_ct_sample(
      seu_x,
      celltypes = celltypes,
      ct_col = ct_col_xenium,
      sample_col = "sample",
      scale.factor = scale_factor
    )
    saveRDS(pb_list, knn_path)
    log_msg("Wrote ", knn_path)
  }

  for (ct in names(pb_list)) {
    s <- pb_list[[ct]]
    if (!"age" %in% colnames(s@meta.data)) {
      s$age <- as.numeric(as.character(s$sample))
    } else {
      s$age <- as.numeric(as.character(s$age))
    }
    pb_list[[ct]] <- s
  }

  for (ct in intersect(celltypes, names(pb_list))) {
    log_msg("Train Xmodel: ", ct, " n=", ncol(pb_list[[ct]]))
    model <- Elastic_net_one_round_Xenium(
      seu = pb_list[[ct]],
      scale.factor = scale_factor,
      train_frac = train_frac,
      ncores = n_cores,
      seed = seed
    )
    out_m <- xmodel_path(ct)
    saveRDS(model, out_m)
    log_msg("  saved ", out_m)
  }
  log_msg("========== TRAIN done ==========")
}

# =============================================================================
# PREDICT + boxplot only
# =============================================================================

run_predict <- function() {
  log_msg("========== PREDICT OSK (OSK_clocks_New_2026 style) ==========")

  models <- list()
  for (ct in celltypes) {
    fp <- xmodel_path(ct)
    if (!file.exists(fp)) {
      log_msg("WARNING missing: ", fp)
      next
    }
    models[[ct]] <- readRDS(fp)
    log_msg("  loaded ", fp)
  }
  if (!length(models)) {
    stop("No Mouse_*_Xmodel_1e5_sliceCleaned.rds under ", clock_dir,
         "\nRun --mode train first, or set --clock_dir")
  }

  knn_path <- ensure_rds(osk_knn_rds)
  # reuse existing OSK smooth if present (or --skip_osk_knn)
  if ((skip_osk_knn || file.exists(knn_path) || file.exists(osk_knn_rds)) &&
      (file.exists(knn_path) || file.exists(osk_knn_rds))) {
    fp <- if (file.exists(knn_path)) knn_path else osk_knn_rds
    log_msg("Reuse existing OSK knn/smooth: ", fp)
    pb_osk <- readRDS(fp)
  } else {
    log_msg("OSK knn not found; will smooth now -> ", knn_path)
    seu_o <- load_seurat_flexible(osk_seurat)
    if (!ct_col_osk %in% colnames(seu_o@meta.data)) {
      stop("OSK missing ", ct_col_osk)
    }
    seu_o <- ensure_sct(seu_o)
    sample_col_use <- sample_col_osk
    if (!sample_col_use %in% colnames(seu_o@meta.data)) {
      alt <- intersect(
        c("sample", "sample_id", "sample1"),
        colnames(seu_o@meta.data)
      )
      if (!length(alt)) stop("OSK missing sample column")
      sample_col_use <- alt[[1]]
      log_msg("Using sample_col=", sample_col_use)
    }
    smp <- as.character(seu_o[[sample_col_use]][, 1])
    seu_o$sample <- gsub("^R([1-4])$", "OSKR\\1", smp)
    log_msg("OSK n=", ncol(seu_o),
            " samples=", paste(sort(unique(seu_o$sample)), collapse = ","))

    pb_osk <- smooth_by_ct_sample(
      seu_o,
      celltypes = intersect(celltypes, names(models)),
      ct_col = ct_col_osk,
      sample_col = "sample",
      scale.factor = scale_factor
    )
    saveRDS(pb_osk, knn_path)
    log_msg("Wrote ", knn_path)
  }

  prepared <- Prepare_Matrix_And_Meta_OSK(pb_osk, scale.factor = scale_factor)

  pred_list <- lapply(intersect(names(prepared), names(models)), function(ct) {
    log_msg("Predict CT=", ct)
    mat <- prepared[[ct]]$Matrix
    meta <- prepared[[ct]]$Meta
    fit <- models[[ct]]$model
    res <- PredictAgeFromModel(mat, meta, fit)
    res$class <- ct
    res$celltype <- ct
    res
  })
  all_res <- do.call(rbind, pred_list)
  rownames(all_res) <- NULL

  # Original filtered OSKR3/OSKR4 first, then drew boxplots; here mapped to Ctrl/OSK
  all_res$condition <- unname(cond_map[as.character(all_res$sample)])
  d <- all_res[all_res$condition %in% c("Ctrl", "OSK"), , drop = FALSE]
  if (!nrow(d)) {
    stop("No Ctrl/OSK after map. samples=",
         paste(unique(all_res$sample), collapse = ","))
  }
  d$condition <- factor(d$condition, levels = c("Ctrl", "OSK"))
  d$age <- factor(as.character(d$sample), levels = c("OSKR3", "OSKR4"))
  # display labels to match prior figure style
  d$cond_lab <- factor(
    ifelse(d$condition == "Ctrl", "Control Diet 58wk", "OSK Diet 58wk"),
    levels = c("Control Diet 58wk", "OSK Diet 58wk")
  )

  ct_lev <- intersect(ct_order, unique(as.character(d$class)))
  ct_lev <- c(ct_lev, setdiff(unique(as.character(d$class)), ct_lev))
  d$class <- factor(as.character(d$class), levels = ct_lev)

  # one-sided Wilcoxon: H1 = OSK < Ctrl (rejuvenation)
  format_p_lab <- function(p) {
    if (is.na(p)) return("p = NA")
    if (p < 2.2e-16) return("p < 2.2e-16")
    if (p < 0.001) return(sprintf("p = %.2e", p))
    sprintf("p = %s", format(signif(p, 3), scientific = FALSE))
  }

  wil_rows <- list()
  log_msg("========== Wilcoxon one-sided: OSK < Ctrl ==========")
  for (ct in levels(d$class)) {
    dc <- d[d$class == ct, , drop = FALSE]
    x_ctrl <- dc$preds[dc$condition == "Ctrl"]
    x_osk <- dc$preds[dc$condition == "OSK"]
    if (length(x_ctrl) < 2L || length(x_osk) < 2L) {
      log_msg("  [", ct, "] skip (n too small)")
      next
    }
    wt <- wilcox.test(x_osk, x_ctrl, alternative = "less")
    pval <- wt$p.value
    log_msg(sprintf(
      "  [%s] Ctrl mean=%.3f (n=%d) | OSK mean=%.3f (n=%d) | Δ(OSK-Ctrl)=%+.3f | p(OSK<Ctrl)=%s",
      ct,
      mean(x_ctrl), length(x_ctrl),
      mean(x_osk), length(x_osk),
      mean(x_osk) - mean(x_ctrl),
      format.pval(pval, digits = 3)
    ))
    wil_rows[[ct]] <- data.frame(
      celltype = ct,
      n_Ctrl = length(x_ctrl),
      n_OSK = length(x_osk),
      mean_Ctrl = mean(x_ctrl),
      mean_OSK = mean(x_osk),
      delta_OSK_minus_Ctrl = mean(x_osk) - mean(x_ctrl),
      wilcox_p_OSK_less_Ctrl = pval,
      p_lab = format_p_lab(pval),
      stringsAsFactors = FALSE
    )
  }
  wil <- if (length(wil_rows)) do.call(rbind, wil_rows) else NULL
  if (!is.null(wil)) {
    wil_path <- file.path(
      out_dir, "OSK_wilcox_OSK_less_Ctrl_by_CT_1e5_sliceCleaned_OSKclean.tsv"
    )
    write.table(wil, wil_path, sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("Wrote wilcox TSV: ", wil_path)
  }

  # p-label position per facet (bottom-right-ish)
  p_anno <- NULL
  if (!is.null(wil)) {
    p_anno <- do.call(rbind, lapply(seq_len(nrow(wil)), function(i) {
      ct <- wil$celltype[i]
      dc <- d[d$class == ct, , drop = FALSE]
      yr <- range(dc$preds, na.rm = TRUE)
      data.frame(
        class = factor(ct, levels = levels(d$class)),
        x = 1.5,
        y = yr[1] + 0.08 * diff(yr),
        p_lab = wil$p_lab[i],
        stringsAsFactors = FALSE
      )
    }))
  }

  # Nature-style muted palette (print-friendly)
  fill_cols <- c(
    "Control Diet 58wk" = "#D9E2EC",
    "OSK Diet 58wk" = "#F1D6D0"
  )
  edge_cols <- c(
    "Control Diet 58wk" = "#3D5A80",
    "OSK Diet 58wk" = "#9B2226"
  )

  p <- ggplot(d, aes(x = cond_lab, y = preds, fill = cond_lab, color = cond_lab)) +
    facet_wrap(~ class, nrow = 2, scales = "free_y") +
    geom_boxplot(
      outlier.shape = NA, linewidth = 0.55, width = 0.48,
      fatten = 1.8
    ) +
    scale_fill_manual(values = fill_cols, guide = "none") +
    scale_color_manual(values = edge_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
    theme_classic(base_family = "Helvetica", base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size = 7, color = "#222222"),
      axis.text.y = element_text(size = 8, color = "#333333"),
      axis.title.y = element_text(size = 10, color = "#222222"),
      axis.line = element_line(color = "#222222", linewidth = 0.4),
      axis.ticks = element_line(color = "#222222", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 10, color = "#111111"),
      strip.background = element_blank(),
      panel.border = element_rect(color = "#222222", fill = NA, linewidth = 0.6),
      panel.spacing.x = unit(0.25, "cm"),
      panel.spacing.y = unit(0.35, "cm"),
      legend.position = "none",
      plot.margin = margin(6, 6, 4, 6)
    ) +
    xlab("") +
    ylab("Predicted")

  if (!is.null(p_anno)) {
    p <- p + geom_text(
      data = p_anno,
      aes(x = x, y = y, label = p_lab),
      inherit.aes = FALSE,
      size = 2.6,
      hjust = 0.5,
      color = "#111111",
      family = "Helvetica"
    )
  }

  n_ct <- length(levels(d$class))
  # narrower figure (leaner panels)
  plot_w <- max(5.8, 1.35 * ceiling(n_ct / 2))
  plot_h <- 4.4
  fp <- file.path(out_dir, "Mouse_Xenium_OSK_2026final_1e5_sliceCleaned_OSKclean.png")
  ggsave(fp, p, height = plot_h, width = plot_w, dpi = 600)
  log_msg("Boxplot: ", fp, " (", plot_w, " x ", plot_h, " in)")

  fp2 <- file.path(out_dir, "OSK_Ctrl_vs_OSK_boxplot_SCT_knn_2026final_1e5_sliceCleaned_OSKclean.png")
  ggsave(fp2, p, height = plot_h, width = plot_w, dpi = 600)
  log_msg("Boxplot: ", fp2)

  # also save predictions used for the plot
  pred_out <- file.path(
    out_dir, "OSK_predict_by_Xenium_SCT_knn_1e5_sliceCleaned_OSKclean.tsv"
  )
  write.table(d, pred_out, sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("Wrote predictions: ", pred_out)

  log_msg("========== PREDICT done ==========")
}

# =============================================================================
# Main
# =============================================================================

log_msg("========== 2026final SCT-kNN 1e5 | sliceCleaned aging + OSKclean start ==========")
log_msg("mode       : ", mode)
log_msg("xenium     : ", xenium_seurat)
log_msg("osk        : ", osk_seurat)
log_msg("clock_dir  : ", clock_dir)
log_msg("out_dir    : ", out_dir)
log_msg("preprocess : SCT-kNN sum -> depth ", scale_factor, " -> log2")
log_msg("NOTE       : aging=sliceCleaned; OSK=sliceCleaned_OSKR3R4; clocks reuse previous sliceCleaned dir")

if (!mode %in% c("train", "predict", "all")) {
  stop("--mode must be train | predict | all")
}

if (mode %in% c("train", "all")) run_train()
if (mode %in% c("predict", "all")) run_predict()

log_msg("========== all done ==========")
log_msg("Figure: ", file.path(out_dir, "Mouse_Xenium_OSK_2026final_1e5_sliceCleaned_OSKclean.png"))
