#!/usr/bin/env Rscript
# Per-slice, per-celltype spatial expression smooth → aging score → accel → viz.
#
# Pipeline:
#   1) For each panel × celltype: spatial SUM of ALL-gene raw RNA counts on (x,y)
#      among same-CT cells (default: small radius; optional: kNN); then depth-norm
#   2) LogNormalize on ALL genes (library size = full transcriptome)
#   3) AddModuleScore on full gene space (controls sampled from all genes);
#      aging = Up - Down
#   4) Within-slice relative acceleration (z): (score-mean)/sd | panel×CT×gene_set
#   5) Spatial map (white bg): INL=grey, RPE=dark green, MG=black,
#      Rod=aging color (Young blue → Old red); Cone not drawn
#
# Prerequisites:
#   summarize_4results_consistency.R  → consistency_genes_long_p0.05.tsv
#   cleaned Seurat with x/y/panel/celltype
#
# Example:
#   Rscript spatial_smooth_aging_score_accel.R
#   Rscript spatial_smooth_aging_score_accel.R --radius auto --radius_mult 3
#   Rscript spatial_smooth_aging_score_accel.R --radius 40
#   Rscript spatial_smooth_aging_score_accel.R --smooth_mode knn --k 50
#
# Plot only one slice + zoom to a region:
#   (auto-skips smooth if cell_spatialSmooth_aging_accel_long.tsv already exists)
#   Rscript spatial_smooth_aging_score_accel.R \
#     --panel R4_rep1 --xmin 2000 --xmax 4000 --ymin 1000 --ymax 3000
#   # force recompute scores:
#   Rscript spatial_smooth_aging_score_accel.R --force
#   # coords are native meta$x / meta$y (before --match_seurat swap)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(FNN)
  library(ggplot2)
  library(scales)
  library(ggnewscale) # dual color scales: MG + Rod on one map
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}
has_flag <- function(flag) !is.na(match(flag, args))

root <- "/projects/hmz-aging/ProximityEffect"
seurat_path <- get_arg(
  "--seurat",
  file.path(
    root, "Spatial_slice_clean",
    "R1_4_seurat_Merged_dim_RNA_2024_sliceCleaned_colored.rds"
  )
)
genes_long_path <- get_arg(
  "--genes_long",
  file.path(
    root, "Merged_4results_consistency_p0.05",
    "consistency_genes_long_p0.05.tsv"
  )
)
out_dir <- get_arg(
  "--out_dir",
  file.path(
    root, "Spatial_slice_clean",
    "spatial_smooth_aging_score_accel"
  )
)

# Only MG + Rod by default (faster; other CTs not scored/plotted)
celltypes <- strsplit(
  get_arg("--celltypes", "MG,Rod"),
  ",",
  fixed = TRUE
)[[1]]
plot_cts <- strsplit(get_arg("--plot_cts", "MG,Rod"), ",", fixed = TRUE)[[1]]
# only plot these (must be subset of celltypes)
plot_cts <- intersect(plot_cts, celltypes)
gene_sets <- strsplit(
  get_arg("--gene_set", "all4,exactly3"),
  ",",
  fixed = TRUE
)[[1]]

# smooth: radius (default) or knn — both SUM raw counts (not mean)
smooth_mode <- tolower(get_arg("--smooth_mode", "radius"))
if (!smooth_mode %in% c("radius", "knn")) {
  stop("--smooth_mode must be 'radius' or 'knn'")
}
# knn: SUM of self + k nearest same-CT cells in (x,y)
k_nn <- as.integer(get_arg("--k", "50"))
# radius: numeric (same units as meta$x/y) or "auto"
# auto = median(1st NN dist within slice×CT) * radius_mult  (smaller mult = tighter)
radius_arg <- get_arg("--radius", "auto")
radius_fixed <- if (identical(tolower(radius_arg), "auto")) {
  NA_real_
} else {
  as.numeric(radius_arg)
}
radius_mult <- as.numeric(get_arg("--radius_mult", "3"))
# safety cap for radius neighbor list (RANN max k)
max_neighbors <- as.integer(get_arg("--max_neighbors", "500"))
min_cells <- as.integer(get_arg("--min_cells", "5"))
min_genes <- as.integer(get_arg("--min_genes", "3"))
assay_use <- get_arg("--assay", "RNA")
n_cores <- as.integer(get_arg("--n_cores", "8"))
point_size <- as.numeric(get_arg("--point_size", "1.5"))
dpi <- as.integer(get_arg("--dpi", "300"))
slice_size_in <- as.numeric(get_arg("--slice_size_in", "8"))
match_seurat <- !has_flag("--no_match_seurat")
save_seu <- has_flag("--save_seurat")
color_lim <- get_arg("--color_lim", NULL)
if (!is.null(color_lim)) color_lim <- as.numeric(color_lim)

# ---- plot subset: one/several panels + optional xy crop ----
# --panel R4_rep1   or  --panels R4_rep1,R4_rep2
panels_arg <- get_arg("--panels", get_arg("--panel", NULL))
panels_plot <- if (is.null(panels_arg) || !nzchar(panels_arg)) {
  character(0)
} else {
  unique(strsplit(panels_arg, ",", fixed = TRUE)[[1]])
}
# crop box in native x/y (same units as Seurat meta); NA = no limit
parse_lim <- function(v) {
  if (is.null(v) || !nzchar(v) || identical(tolower(v), "na")) NA_real_ else as.numeric(v)
}
xmin <- parse_lim(get_arg("--xmin", NULL))
xmax <- parse_lim(get_arg("--xmax", NULL))
ymin <- parse_lim(get_arg("--ymin", NULL))
ymax <- parse_lim(get_arg("--ymax", NULL))
# skip smooth/score; only replot from existing long TSV
plot_only <- has_flag("--plot_only")
force_recompute <- has_flag("--force") || has_flag("--recompute")
score_tsv_arg <- get_arg("--score_tsv", NULL)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "figures_spatial")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "by_slice"), recursive = TRUE, showWarnings = FALSE)
gl_dir <- file.path(out_dir, "genelists_used")
dir.create(gl_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
  flush.console()
}

sample1_to_age_wk <- function(sample1) {
  map <- c(R1 = 5, R2 = 12, R3 = 52, R4 = 112)
  unname(map[as.character(sample1)])
}

get_counts <- function(seu, assay = "RNA") {
  DefaultAssay(seu) <- assay
  if (exists("JoinLayers", mode = "function")) {
    tryCatch(seu <- JoinLayers(seu, assay = assay), error = function(e) NULL)
  }
  out <- NULL
  if (exists("LayerData", mode = "function")) {
    out <- tryCatch(
      LayerData(seu, assay = assay, layer = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(out)) {
    out <- tryCatch(
      GetAssayData(seu, assay = assay, layer = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(out)) out <- tryCatch(seu[[assay]]$counts, error = function(e) NULL)
  if (is.null(out) || length(out) == 0) stop("Cannot extract counts")
  if (!inherits(out, "dgCMatrix")) out <- as(as(out, "CsparseMatrix"), "dgCMatrix")
  out
}

apply_smooth_A <- function(counts_sub, A) {
  out <- counts_sub %*% A
  colnames(out) <- colnames(counts_sub)
  rownames(out) <- rownames(counts_sub)
  out
}

# kNN smooth: SUM raw counts of self + k neighbors
smooth_counts_knn <- function(counts_sub, xy, k = 50L) {
  n <- ncol(counts_sub)
  if (n < 2) return(list(sm = counts_sub, radius = NA_real_, k = 0L, n_nbr_med = 0))
  k_use <- min(as.integer(k), n - 1L)
  if (k_use < 1L) {
    return(list(sm = counts_sub, radius = NA_real_, k = 0L, n_nbr_med = 0))
  }
  knn <- FNN::get.knn(xy, k = k_use)
  rows <- as.integer(c(t(knn$nn.index), seq_len(n)))
  cols <- rep(seq_len(n), each = k_use + 1L)
  A <- sparseMatrix(i = rows, j = cols, x = 1, dims = c(n, n))
  list(
    sm = apply_smooth_A(counts_sub, A),
    radius = NA_real_,
    k = k_use,
    n_nbr_med = k_use + 1L
  )
}

median_nn1_dist <- function(xy) {
  n <- nrow(xy)
  if (n < 2) return(NA_real_)
  d1 <- FNN::get.knn(xy, k = 1L)$nn.dist[, 1]
  as.numeric(stats::median(d1, na.rm = TRUE))
}

# Radius smooth: SUM raw counts of all same-CT cells within radius (includes self)
# Uses RANN radius search when available; else FNN all-k then filter.
smooth_counts_radius <- function(counts_sub, xy, radius, max_neighbors = 500L) {
  n <- ncol(counts_sub)
  if (n < 2 || !is.finite(radius) || radius <= 0) {
    return(list(sm = counts_sub, radius = radius, k = NA_integer_, n_nbr_med = 1))
  }
  max_k <- min(as.integer(max_neighbors), n)

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(
      data = xy, query = xy, k = max_k,
      searchtype = "radius", radius = radius
    )
    idx <- nn$nn.idx
    # RANN: unused slots are 0
    ok <- idx > 0L
  } else {
    # fallback: take up to max_k nearest, keep those with dist <= radius
    k_use <- min(max_k, n - 1L)
    if (k_use < 1L) {
      return(list(sm = counts_sub, radius = radius, k = NA_integer_, n_nbr_med = 1))
    }
    knn <- FNN::get.knn(xy, k = k_use)
    idx <- cbind(seq_len(n), knn$nn.index) # include self
    dists <- cbind(0, knn$nn.dist)
    ok <- dists <= radius
  }

  rows <- as.integer(idx[ok])
  cols <- as.integer(row(idx)[ok])
  # guarantee self (in case radius search missed due to numeric edge)
  rows <- c(rows, seq_len(n))
  cols <- c(cols, seq_len(n))
  A <- sparseMatrix(i = rows, j = cols, x = 1, dims = c(n, n))
  A@x <- rep(1, length(A@x)) # binarize (self may be duplicated)
  deg <- as.numeric(Matrix::colSums(A))
  # NO divide by degree: keep SUM of raw counts; depth-norm later
  list(
    sm = apply_smooth_A(counts_sub, A),
    radius = radius,
    k = NA_integer_,
    n_nbr_med = as.integer(stats::median(deg))
  )
}

sym_lim <- function(v, fixed = NULL) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(-1, 1))
  if (!is.null(fixed) && is.finite(fixed)) return(c(-abs(fixed), abs(fixed)))
  m <- as.numeric(stats::quantile(abs(v), 0.98, na.rm = TRUE))
  if (!is.finite(m) || m <= 0) m <- max(abs(v), na.rm = TRUE)
  if (!is.finite(m) || m <= 0) m <- 1
  c(-m, m)
}

# anatomical context colors (plot only; not scored)
col_inl <- "#9E9E9E"   # grey (INL = BC/AC/HC)
col_mg  <- "#000000"   # black
col_rpe <- "#145A32"   # dark green
pal_rod_age <- list(low = "#2166AC", mid = "#F2F2F2", high = "#B2182B") # Young→Old
inl_plot_cts <- c("BC", "AC", "HC")

log_msg("========== spatial_smooth_aging_score_accel start ==========")
log_msg("seurat : ", seurat_path)
log_msg("genes  : ", genes_long_path)
log_msg("out    : ", out_dir)
log_msg("CTs (process): ", paste(celltypes, collapse = ", "))
log_msg("plot CTs     : ", paste(plot_cts, collapse = ", "))
log_msg(
  "smooth_mode : ", smooth_mode,
  if (smooth_mode == "radius") {
    paste0(
      " | radius=", radius_arg,
      if (!is.finite(radius_fixed)) paste0(" (mult=", radius_mult, ")") else "",
      " | max_neighbors=", max_neighbors
    )
  } else {
    paste0(" | k=", k_nn)
  },
  " | gene_sets: ", paste(gene_sets, collapse = ", ")
)
log_msg("n_cores: ", n_cores, " | save_seurat: ", save_seu)
if (!length(plot_cts)) stop("No plot_cts left after intersect with celltypes")
if (length(panels_plot)) {
  log_msg("plot panels filter: ", paste(panels_plot, collapse = ", "))
}
if (any(is.finite(c(xmin, xmax, ymin, ymax)))) {
  log_msg(
    "plot crop (native x/y): xmin=", xmin, " xmax=", xmax,
    " ymin=", ymin, " ymax=", ymax
  )
}
log_msg("plot_only flag: ", plot_only, " | force_recompute: ", force_recompute)

crop_xy <- function(d) {
  if (is.null(d) || !nrow(d)) return(d)
  ok <- rep(TRUE, nrow(d))
  if (is.finite(xmin)) ok <- ok & is.finite(d$x) & d$x >= xmin
  if (is.finite(xmax)) ok <- ok & is.finite(d$x) & d$x <= xmax
  if (is.finite(ymin)) ok <- ok & is.finite(d$y) & d$y >= ymin
  if (is.finite(ymax)) ok <- ok & is.finite(d$y) & d$y <= ymax
  d[ok, , drop = FALSE]
}

# Auto: if score TSV already exists → reuse (plot only), unless --force/--recompute
score_tsv_default <- file.path(out_dir, "cell_spatialSmooth_aging_accel_long.tsv")
score_tsv_use <- if (!is.null(score_tsv_arg) && nzchar(score_tsv_arg)) {
  score_tsv_arg
} else {
  score_tsv_default
}
if (!isTRUE(force_recompute) && file.exists(score_tsv_use)) {
  if (!isTRUE(plot_only)) {
    log_msg(
      "Found existing scores → auto plot-only (skip smooth/score): ",
      score_tsv_use
    )
    log_msg("  To recompute, pass --force or --recompute")
  }
  plot_only <- TRUE
} else if (isTRUE(plot_only) && !file.exists(score_tsv_use)) {
  stop("plot_only requested but score TSV missing: ", score_tsv_use)
}

# ---------------------------------------------------------------------------
# plot_only: reload scores + context, skip smooth/score
# ---------------------------------------------------------------------------

if (isTRUE(plot_only)) {
  score_tsv <- score_tsv_use
  ctx_path <- file.path(out_dir, "plot_context_cells.tsv")
  if (!file.exists(score_tsv)) {
    stop("plot_only: score TSV not found: ", score_tsv)
  }
  if (!file.exists(ctx_path)) {
    log_msg("plot_only: no plot_context_cells.tsv — build from Seurat meta ...")
    if (!file.exists(seurat_path)) stop("Seurat not found for context: ", seurat_path)
    seu_tmp <- readRDS(seurat_path)
    meta_tmp <- seu_tmp@meta.data
    meta_tmp$cell_id <- rownames(meta_tmp)
    if (!"panel" %in% colnames(meta_tmp)) {
      if (all(c("sample1", "slice") %in% colnames(meta_tmp))) {
        meta_tmp$panel <- paste(meta_tmp$sample1, meta_tmp$slice, sep = "_")
      } else if (all(c("sample1", "rep") %in% colnames(meta_tmp))) {
        meta_tmp$panel <- paste(meta_tmp$sample1, meta_tmp$rep, sep = "_")
      } else {
        stop("Need panel or sample1+slice/rep for context")
      }
    }
    if (!"sample1" %in% colnames(meta_tmp)) {
      meta_tmp$sample1 <- sub("_.*$", "", meta_tmp$panel)
    }
    meta_tmp <- meta_tmp[
      is.finite(meta_tmp$x) & is.finite(meta_tmp$y) &
        meta_tmp$sample1 %in% c("R1", "R2", "R3", "R4"),
      ,
      drop = FALSE
    ]
    plot_bg <- meta_tmp[
      meta_tmp$celltype %in% c(inl_plot_cts, "RPE", "MG"),
      c("cell_id", "celltype", "panel", "x", "y"),
      drop = FALSE
    ]
    plot_bg$layer <- ifelse(
      plot_bg$celltype %in% inl_plot_cts, "INL",
      as.character(plot_bg$celltype)
    )
    plot_bg$panel <- as.character(plot_bg$panel)
    write.table(plot_bg, ctx_path, sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("Wrote ", ctx_path)
    rm(seu_tmp, meta_tmp)
    gc()
  }
  log_msg("plot_only load scores: ", score_tsv)
  score_long <- read.delim(score_tsv, stringsAsFactors = FALSE)
  plot_bg <- read.delim(ctx_path, stringsAsFactors = FALSE)
  plot_bg$panel <- as.character(plot_bg$panel)
  smooth_label <- get_arg("--smooth_label", "saved scores")
  log_msg(
    "  score rows=", nrow(score_long),
    " | context rows=", nrow(plot_bg),
    " | panels in scores=", length(unique(score_long$panel))
  )
} else {

if (!file.exists(seurat_path)) stop("Seurat not found: ", seurat_path)
if (!file.exists(genes_long_path)) {
  stop("Gene list not found: ", genes_long_path)
}

# ---------------------------------------------------------------------------
# Gene lists
# ---------------------------------------------------------------------------

glong <- read.delim(genes_long_path, stringsAsFactors = FALSE)
need_g <- c("celltype", "Gene", "agree_level", "direction")
miss_g <- setdiff(need_g, colnames(glong))
if (length(miss_g)) stop("genes_long missing: ", paste(miss_g, collapse = ", "))

genelists <- list()
genes_needed <- character(0)
for (ct in celltypes) {
  for (gs in gene_sets) {
    for (d in c("Up", "Down")) {
      key <- paste(ct, gs, d, sep = "__")
      genes <- unique(as.character(glong$Gene[
        glong$celltype == ct &
          glong$agree_level == gs &
          glong$direction == d
      ]))
      genes <- genes[!is.na(genes) & genes != ""]
      genelists[[key]] <- genes
      genes_needed <- union(genes_needed, genes)
      writeLines(genes, file.path(gl_dir, paste0(key, ".txt")))
      log_msg("genelist ", key, ": n=", length(genes))
    }
  }
}
log_msg("Union genes for scoring: ", length(genes_needed))
if (length(genes_needed) < 5) stop("Too few genes in union gene lists")

# ---------------------------------------------------------------------------
# Load Seurat meta + counts
# ---------------------------------------------------------------------------

log_msg("Loading Seurat ...")
t0 <- Sys.time()
seu <- readRDS(seurat_path)
if (exists("UpdateSeuratObject", mode = "function")) {
  tryCatch(seu <- UpdateSeuratObject(seu), error = function(e) NULL)
}
meta <- seu@meta.data
meta$cell_id <- rownames(meta)

need_m <- c("celltype", "x", "y")
miss_m <- setdiff(need_m, colnames(meta))
if (length(miss_m)) stop("meta missing: ", paste(miss_m, collapse = ", "))

if (!"panel" %in% colnames(meta)) {
  if (all(c("sample1", "slice") %in% colnames(meta))) {
    meta$panel <- paste(meta$sample1, meta$slice, sep = "_")
  } else if (all(c("sample1", "rep") %in% colnames(meta))) {
    meta$panel <- paste(meta$sample1, meta$rep, sep = "_")
  } else {
    stop("Need meta$panel or sample1+slice/rep")
  }
}
if (!"sample1" %in% colnames(meta)) {
  meta$sample1 <- sub("_.*$", "", meta$panel)
}
meta$age_wk <- sample1_to_age_wk(meta$sample1)

counts <- get_counts(seu, assay_use)
common <- intersect(meta$cell_id, colnames(counts))
meta <- meta[common, , drop = FALSE]
counts <- counts[, common, drop = FALSE]

# keep anatomical context cells for plotting (INL=BC/AC/HC, RPE, MG); no Cone
meta_all <- meta[
  is.finite(meta$x) & is.finite(meta$y) &
    meta$sample1 %in% c("R1", "R2", "R3", "R4"),
  ,
  drop = FALSE
]
plot_bg <- meta_all[
  meta_all$celltype %in% c(inl_plot_cts, "RPE", "MG"),
  c("cell_id", "celltype", "panel", "x", "y"),
  drop = FALSE
]
plot_bg$layer <- ifelse(
  plot_bg$celltype %in% inl_plot_cts, "INL",
  as.character(plot_bg$celltype)
)
plot_bg$panel <- as.character(plot_bg$panel)
ctx_path <- file.path(out_dir, "plot_context_cells.tsv")
write.table(plot_bg, ctx_path, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote ", ctx_path)
log_msg(
  "plot context cells: ",
  paste(sprintf("%s=%d", names(table(plot_bg$layer)), table(plot_bg$layer)), collapse = ", ")
)

meta <- meta[
  meta$celltype %in% celltypes &
    is.finite(meta$x) & is.finite(meta$y) &
    meta$sample1 %in% c("R1", "R2", "R3", "R4"),
  ,
  drop = FALSE
]
counts <- counts[, meta$cell_id, drop = FALSE]

# Smooth + depth-norm + AddModuleScore all need the FULL gene space:
#   - library size from all genes
#   - AddModuleScore samples control genes from all genes (expression bins)
genes_use <- intersect(genes_needed, rownames(counts))
log_msg(
  "all genes kept: ", nrow(counts),
  " | module genes available: ", length(genes_use)
)
if (length(genes_use) < 5) stop("Too few module genes in counts")

log_msg(
  "cells for smooth: ", nrow(meta),
  " | load min=",
  signif(as.numeric(difftime(Sys.time(), t0, units = "mins")), 3)
)

# ---------------------------------------------------------------------------
# Spatial smooth per panel × CT (parallel)
# ---------------------------------------------------------------------------

log_msg("Spatial smooth (per slice × CT, parallel) ...")
panels <- sort(unique(as.character(meta$panel)))
jobs <- list()
for (pn in panels) {
  for (ct in celltypes) {
    idx <- which(meta$panel == pn & meta$celltype == ct)
    if (length(idx) < min_cells) {
      log_msg("  SKIP ", pn, " ", ct, ": n=", length(idx))
      next
    }
    jobs[[paste(pn, ct, sep = "||")]] <- list(
      key = paste(pn, ct, sep = "||"),
      pn = pn,
      ct = ct,
      idx = idx
    )
  }
}
if (!length(jobs)) stop("No groups to smooth")

smooth_one_job <- function(job) {
  idx <- job$idx
  cells <- meta$cell_id[idx]
  xy <- as.matrix(meta[idx, c("x", "y"), drop = FALSE])
  storage.mode(xy) <- "double"
  rownames(xy) <- cells
  sub <- counts[, cells, drop = FALSE]

  if (identical(smooth_mode, "knn")) {
    res <- smooth_counts_knn(sub, xy, k = k_nn)
  } else {
    rad <- radius_fixed
    if (!is.finite(rad)) {
      nn1 <- median_nn1_dist(xy)
      rad <- nn1 * radius_mult
      if (!is.finite(rad) || rad <= 0) {
        # fallback: knn if auto radius fails
        res <- smooth_counts_knn(sub, xy, k = k_nn)
        return(list(
          key = job$key, sm = res$sm, meta = meta[idx, , drop = FALSE],
          radius = NA_real_, k = res$k, n_nbr_med = res$n_nbr_med,
          note = "auto-radius-failed→knn"
        ))
      }
    }
    res <- smooth_counts_radius(sub, xy, radius = rad, max_neighbors = max_neighbors)
  }
  list(
    key = job$key,
    sm = res$sm,
    meta = meta[idx, , drop = FALSE],
    radius = res$radius,
    k = res$k,
    n_nbr_med = res$n_nbr_med
  )
}

t_sm <- Sys.time()
if (.Platform$OS.type != "windows" && n_cores > 1L && length(jobs) > 1L) {
  res_list <- parallel::mclapply(
    jobs, smooth_one_job,
    mc.cores = min(n_cores, length(jobs))
  )
} else {
  res_list <- lapply(jobs, smooth_one_job)
}
# check errors
bad <- vapply(res_list, inherits, logical(1), what = "try-error")
if (any(bad)) stop("Smooth failed for some jobs")
rad_used <- vapply(res_list, function(x) {
  as.numeric(x$radius)[1]
}, numeric(1))
nbr_med <- vapply(res_list, function(x) {
  as.integer(x$n_nbr_med)[1]
}, integer(1))
log_msg(
  "  smooth done in ",
  signif(as.numeric(difftime(Sys.time(), t_sm, units = "secs")), 3),
  " sec | groups=", length(res_list)
)
if (identical(smooth_mode, "radius")) {
  log_msg(
    "  radius summary: median=",
    signif(stats::median(rad_used, na.rm = TRUE), 4),
    " | range=",
    signif(min(rad_used, na.rm = TRUE), 4),
    "-",
    signif(max(rad_used, na.rm = TRUE), 4),
    " | median_n_neighbors≈",
    stats::median(nbr_med, na.rm = TRUE)
  )
}
# for plot subtitle
smooth_label <- if (identical(smooth_mode, "knn")) {
  paste0("kNN k=", k_nn)
} else if (is.finite(radius_fixed)) {
  paste0("radius=", radius_fixed)
} else {
  paste0(
    "radius=auto(×", radius_mult, "), med≈",
    signif(stats::median(rad_used, na.rm = TRUE), 4)
  )
}

sm_mat <- do.call(cbind, lapply(res_list, function(x) x$sm))
meta_use <- do.call(rbind, lapply(res_list, `[[`, "meta"))
rownames(meta_use) <- meta_use$cell_id
meta_use <- meta_use[colnames(sm_mat), , drop = FALSE]
log_msg("Smoothed matrix (all genes): genes=", nrow(sm_mat), " cells=", ncol(sm_mat))

rm(seu, counts, res_list, jobs)
gc()

# ---------------------------------------------------------------------------
# Normalize + AddModuleScore on ALL genes
# ---------------------------------------------------------------------------

log_msg("CreateSeuratObject + NormalizeData (all genes) ...")
op <- options(Seurat.object.assay.version = "v3")
on.exit(options(op), add = TRUE)

seu_sm <- CreateSeuratObject(
  counts = sm_mat,
  meta.data = meta_use,
  assay = "RNA",
  project = "spatial_smooth"
)
rm(sm_mat)
gc()

seu_sm <- NormalizeData(
  seu_sm,
  normalization.method = "LogNormalize",
  verbose = FALSE
)

features_all <- rownames(seu_sm)
log_msg(
  "AddModuleScore on full gene space: genes=", length(features_all),
  " | module genes=", length(intersect(genes_use, features_all))
)
score_df <- seu_sm@meta.data
score_df$cell_id <- rownames(score_df)

for (ct in celltypes) {
  cells_ct <- colnames(seu_sm)[seu_sm$celltype == ct]
  if (!length(cells_ct)) next
  for (gs in gene_sets) {
    up_genes <- intersect(genelists[[paste(ct, gs, "Up", sep = "__")]], features_all)
    dn_genes <- intersect(genelists[[paste(ct, gs, "Down", sep = "__")]], features_all)
    log_msg(sprintf(
      "  score %s | %s | Up=%d Down=%d cells=%d",
      ct, gs, length(up_genes), length(dn_genes), length(cells_ct)
    ))

    up_col <- paste0("Up_", gs, "_", ct)
    dn_col <- paste0("Down_", gs, "_", ct)
    age_col <- paste0("aging_score_", gs)

    up_score <- rep(NA_real_, ncol(seu_sm))
    dn_score <- rep(NA_real_, ncol(seu_sm))
    names(up_score) <- colnames(seu_sm)
    names(dn_score) <- colnames(seu_sm)

    if (length(up_genes) >= min_genes) {
      nm <- gsub("[^A-Za-z0-9_]", "_", paste("tmpUp", gs, ct, sep = "_"))
      created <- paste0(nm, "1")
      if (created %in% colnames(seu_sm@meta.data)) seu_sm[[created]] <- NULL
      seu_sm <- AddModuleScore(
        seu_sm, features = list(up_genes), name = nm,
        assay = "RNA", search = FALSE
      )
      up_score[cells_ct] <- seu_sm@meta.data[cells_ct, created]
      seu_sm[[created]] <- NULL
    }
    if (length(dn_genes) >= min_genes) {
      nm <- gsub("[^A-Za-z0-9_]", "_", paste("tmpDn", gs, ct, sep = "_"))
      created <- paste0(nm, "1")
      if (created %in% colnames(seu_sm@meta.data)) seu_sm[[created]] <- NULL
      seu_sm <- AddModuleScore(
        seu_sm, features = list(dn_genes), name = nm,
        assay = "RNA", search = FALSE
      )
      dn_score[cells_ct] <- seu_sm@meta.data[cells_ct, created]
      seu_sm[[created]] <- NULL
    }

    aging <- up_score - dn_score
    aging[is.na(up_score) & !is.na(dn_score)] <- -dn_score[is.na(up_score) & !is.na(dn_score)]
    aging[!is.na(up_score) & is.na(dn_score)] <- up_score[!is.na(up_score) & is.na(dn_score)]

    # store only for this CT cells; others stay NA for this gene_set column
    # use long format instead to avoid CT overwrite
    score_df[[paste0("Up__", gs, "__", ct)]] <- up_score[score_df$cell_id]
    score_df[[paste0("Down__", gs, "__", ct)]] <- dn_score[score_df$cell_id]
    score_df[[paste0("aging__", gs, "__", ct)]] <- aging[score_df$cell_id]
  }
}

# long table
long_list <- list()
for (ct in celltypes) {
  for (gs in gene_sets) {
    uc <- paste0("Up__", gs, "__", ct)
    dc <- paste0("Down__", gs, "__", ct)
    ac <- paste0("aging__", gs, "__", ct)
    if (!ac %in% colnames(score_df)) next
    keep <- score_df$celltype == ct & is.finite(score_df[[ac]])
    if (!any(keep)) next
    long_list[[paste(ct, gs)]] <- data.frame(
      cell_id = score_df$cell_id[keep],
      celltype = ct,
      panel = score_df$panel[keep],
      sample1 = score_df$sample1[keep],
      age_wk = score_df$age_wk[keep],
      x = score_df$x[keep],
      y = score_df$y[keep],
      gene_set = gs,
      Up_score = score_df[[uc]][keep],
      Down_score = score_df[[dc]][keep],
      aging_score = score_df[[ac]][keep],
      stringsAsFactors = FALSE
    )
  }
}
score_long <- do.call(rbind, long_list)
rownames(score_long) <- NULL

# ---------------------------------------------------------------------------
# Within-slice acceleration z-score
# ---------------------------------------------------------------------------

log_msg("Within-slice accel z-score ...")
score_long$group_key <- paste(
  score_long$panel, score_long$celltype, score_long$gene_set, sep = "||"
)
g_mean <- tapply(score_long$aging_score, score_long$group_key, mean, na.rm = TRUE)
g_sd <- tapply(score_long$aging_score, score_long$group_key, sd, na.rm = TRUE)
g_n <- tapply(score_long$aging_score, score_long$group_key, function(x) sum(is.finite(x)))
score_long$slice_ct_mean <- as.numeric(g_mean[score_long$group_key])
score_long$slice_ct_sd <- as.numeric(g_sd[score_long$group_key])
score_long$n_in_group <- as.integer(g_n[score_long$group_key])
score_long$aging_accel <- (score_long$aging_score - score_long$slice_ct_mean) /
  score_long$slice_ct_sd
score_long$aging_accel[
  !is.finite(score_long$slice_ct_sd) |
    score_long$slice_ct_sd <= 0 |
    score_long$n_in_group < 2
] <- NA_real_

out_score <- file.path(out_dir, "cell_spatialSmooth_aging_accel_long.tsv")
write.table(score_long, out_score, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote ", out_score)

if (isTRUE(save_seu)) {
  saveRDS(seu_sm, file.path(out_dir, "seurat_spatialSmoothed_counts.rds"))
  log_msg("Wrote seurat_spatialSmoothed_counts.rds")
}
rm(seu_sm, sm_mat)
gc()

} # end !plot_only pipeline

# ---------------------------------------------------------------------------
# Spatial viz (white bg):
#   INL (BC/AC/HC) = grey, RPE = dark green, MG = black,
#   Rod = Young (blue) ↔ Old (red); Cone not drawn
# ---------------------------------------------------------------------------

log_msg("Plot INL/RPE/MG context + Rod aging color (white background; no Cone)")

layer_pal <- c(INL = col_inl, RPE = col_rpe, MG = col_mg)

panels_avail <- sort(unique(as.character(score_long$panel)))
panels_use <- if (length(panels_plot)) {
  miss_p <- setdiff(panels_plot, panels_avail)
  if (length(miss_p)) {
    log_msg("WARNING panels not in scores: ", paste(miss_p, collapse = ", "))
  }
  intersect(panels_plot, panels_avail)
} else {
  panels_avail
}
if (!length(panels_use)) stop("No panels left to plot")
log_msg("panels to plot: ", paste(panels_use, collapse = ", "))

for (gs in gene_sets) {
  for (pn in panels_use) {
    d_rod <- score_long[
      score_long$panel == pn &
        score_long$gene_set == gs &
        score_long$celltype == "Rod" &
        is.finite(score_long$aging_accel),
      ,
      drop = FALSE
    ]
    d_rod <- crop_xy(d_rod)
    if (!nrow(d_rod)) {
      log_msg("  skip ", pn, " | ", gs, " (no Rod in crop)")
      next
    }

    bg <- plot_bg[plot_bg$panel == pn, , drop = FALSE]
    bg <- crop_xy(bg)
    d_inl <- bg[bg$layer == "INL", , drop = FALSE]
    d_rpe <- bg[bg$layer == "RPE", , drop = FALSE]
    d_mg  <- bg[bg$layer == "MG", , drop = FALSE]

    # color scale from cropped Rod only (or fixed --color_lim)
    lim_rod <- sym_lim(d_rod$aging_accel, color_lim)
    age_lab <- d_rod$age_wk[1]

    # dummy rows for discrete legend (context layers)
    leg_ctx <- data.frame(
      layer = factor(c("INL", "RPE", "MG"), levels = c("INL", "RPE", "MG")),
      x = NA_real_, y = NA_real_,
      stringsAsFactors = FALSE
    )

    crop_tag <- ""
    if (any(is.finite(c(xmin, xmax, ymin, ymax)))) {
      crop_tag <- sprintf(
        "_x%.0f-%.0f_y%.0f-%.0f",
        ifelse(is.finite(xmin), xmin, NA),
        ifelse(is.finite(xmax), xmax, NA),
        ifelse(is.finite(ymin), ymin, NA),
        ifelse(is.finite(ymax), ymax, NA)
      )
      crop_tag <- gsub("NA", "all", crop_tag, fixed = TRUE)
    }

    # data extent in native x/y (printed on figure + log)
    all_xy <- rbind(
      if (nrow(d_rod)) d_rod[, c("x", "y"), drop = FALSE] else NULL,
      if (nrow(d_inl)) d_inl[, c("x", "y"), drop = FALSE] else NULL,
      if (nrow(d_rpe)) d_rpe[, c("x", "y"), drop = FALSE] else NULL,
      if (nrow(d_mg))  d_mg[,  c("x", "y"), drop = FALSE] else NULL
    )
    x_rng <- range(all_xy$x, na.rm = TRUE)
    y_rng <- range(all_xy$y, na.rm = TRUE)
    coord_lab <- sprintf(
      "native x=[%.1f, %.1f]  y=[%.1f, %.1f]",
      x_rng[1], x_rng[2], y_rng[1], y_rng[2]
    )
    log_msg("  ", pn, " | ", gs, " | ", coord_lab)

    make_plot <- function(swap_xy = FALSE) {
      p <- ggplot()
      add_ctx <- function(p, d, col, sz, al) {
        if (!nrow(d)) return(p)
        if (isTRUE(swap_xy)) {
          p + geom_point(
            data = d, aes(x = y, y = x),
            color = col, size = sz, stroke = 0, alpha = al
          )
        } else {
          p + geom_point(
            data = d, aes(x = x, y = y),
            color = col, size = sz, stroke = 0, alpha = al
          )
        }
      }
      p <- add_ctx(p, d_inl, col_inl, point_size * 0.7, 0.55)
      p <- add_ctx(p, d_rpe, col_rpe, point_size * 0.85, 0.85)
      p <- add_ctx(p, d_mg,  col_mg,  point_size * 0.9,  0.90)
      p <- p +
        geom_point(
          data = leg_ctx,
          aes(x = x, y = y, fill = layer),
          shape = 21, size = 3, stroke = 0, na.rm = TRUE
        ) +
        scale_fill_manual(
          values = layer_pal, name = "Context",
          breaks = c("INL", "RPE", "MG"), drop = FALSE
        )
      if (isTRUE(swap_xy)) {
        p <- p + geom_point(
          data = d_rod,
          aes(x = y, y = x, color = aging_accel),
          size = point_size * 0.95, stroke = 0, alpha = 0.9
        )
      } else {
        p <- p + geom_point(
          data = d_rod,
          aes(x = x, y = y, color = aging_accel),
          size = point_size * 0.95, stroke = 0, alpha = 0.9
        )
      }
      xlab <- if (isTRUE(swap_xy)) "y (native)" else "x"
      ylab <- if (isTRUE(swap_xy)) "x (native)" else "y"
      p <- p +
        scale_color_gradient2(
          low = pal_rod_age$low, mid = pal_rod_age$mid, high = pal_rod_age$high,
          midpoint = 0, limits = lim_rod, oob = squish,
          name = "Rod\nYoung → Old"
        ) +
        coord_fixed() +
        labs(
          title = paste0(pn, " | Rod aging accel | ", gs),
          subtitle = paste0(
            "age=", age_lab, " wk; smooth=", smooth_label, "; ",
            coord_lab,
            if (nzchar(crop_tag)) paste0("; crop", crop_tag) else ""
          ),
          x = xlab, y = ylab
        ) +
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey35"),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2),
          panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 11),
          legend.background = element_rect(fill = "white", color = NA)
        )
      if (isTRUE(swap_xy)) p <- p + scale_y_reverse()
      p
    }

    p2 <- make_plot(swap_xy = isTRUE(match_seurat))

    fp2 <- file.path(
      fig_dir, "by_slice",
      paste0(pn, "_INL_RPE_MG_Rod_accel_", gs, crop_tag, ".png")
    )
    ggsave(
      fp2, p2, width = slice_size_in, height = slice_size_in, dpi = dpi,
      bg = "white"
    )
    log_msg("  wrote ", basename(fp2), " | Rod n=", nrow(d_rod))
  }
}

log_msg("========== spatial_smooth_aging_score_accel done ==========")
