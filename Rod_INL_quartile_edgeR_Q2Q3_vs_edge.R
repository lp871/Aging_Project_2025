#!/usr/bin/env Rscript
# =============================================================================
# Rod → INL distance quartiles: edgeR DEG **per age (R1, R2, R3, R4 separately)**
#   Within each sample1:
#     Contrast 1: Q2Q3 vs Q1   (logFC>0 = higher in Q2+Q3)
#     Contrast 2: Q2Q3 vs Q4   (logFC>0 = higher in Q2+Q3)
#     Contrast 3: Q1 vs Q4     (logFC>0 = higher in Q1)
#
# Pseudobulk (per slice/panel):
#   Q1    = SUM Rod counts in Q1
#   Q2Q3  = SUM Rod counts in Q2 ∪ Q3
#   Q4    = SUM Rod counts in Q4
# Design (paired by panel within that age): ~ panel + group
#
# Rank keys: Q2Q3_vs_Q1__R1, Q2Q3_vs_Q4__R2, ...
# GSEA: Rod_INL_quartile_GSEA_Q2Q3_vs_edge.R
#
# Example:
#   Rscript Rod_INL_quartile_edgeR_Q2Q3_vs_edge.R
#   Rscript Rod_INL_quartile_edgeR_Q2Q3_vs_edge.R --sample1 R1,R4 --min_n 3
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(FNN)
  library(edgeR)
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
assign_path <- get_arg(
  "--assign",
  file.path(
    root, "Spatial_slice_clean",
    "Rod_INL_quartile_pseudobulk_aging_score",
    "Rod_INL_distance_quartile_assignment.tsv"
  )
)
out_dir <- get_arg(
  "--out_dir",
  file.path(root, "Spatial_slice_clean", "Rod_INL_quartile_edgeR_Q2Q3_vs_edge_byAge")
)
ortholog_path <- get_arg(
  "--ortholog",
  "/projects/hmz-aging/Zebrafish_New_Datasets/HMZ_ortholog_combined"
)

rod_ct <- get_arg("--rod_ct", "Rod")
inl_cts <- strsplit(get_arg("--inl_cts", "BC,AC,HC,MG"), ",", fixed = TRUE)[[1]]
k_inl <- as.integer(get_arg("--k_inl", "50"))
min_cells_q <- as.integer(get_arg("--min_cells_q", "10"))
min_n <- as.integer(get_arg("--min_n", "3")) # min paired panels per contrast
assay_use <- get_arg("--assay", "RNA")
force_reassign <- has_flag("--force_reassign")
# optional: restrict to age groups, e.g. R1,R4 or all
sample1_keep <- strsplit(get_arg("--sample1", "R1,R2,R3,R4"), ",", fixed = TRUE)[[1]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "DEG_edgeR"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "ranks_mouse"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "ranks_human"), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
  flush.console()
}

q_levels <- c(
  "Q1_0-25%(closest)",
  "Q2_25-50%",
  "Q3_50-75%",
  "Q4_75-100%(farthest)"
)

rod_distance_to_inl <- function(xy_rod, xy_inl, k = 50L) {
  xy_rod <- as.matrix(xy_rod); xy_inl <- as.matrix(xy_inl)
  storage.mode(xy_rod) <- "double"; storage.mode(xy_inl) <- "double"
  n_rod <- nrow(xy_rod); n_inl <- nrow(xy_inl)
  if (!n_rod) return(numeric(0))
  if (!n_inl) return(rep(Inf, n_rod))
  k_use <- min(as.integer(k), n_inl)
  nn <- FNN::get.knnx(data = xy_inl, query = xy_rod, k = k_use)
  if (k_use == 1L) as.numeric(nn$nn.dist[, 1]) else
    as.numeric(apply(nn$nn.dist, 1L, stats::median, na.rm = TRUE))
}

assign_distance_quartile <- function(dist) {
  d <- as.numeric(dist)
  out <- rep(NA_character_, length(d))
  ok <- is.finite(d)
  if (sum(ok) < 4L) {
    out[ok] <- q_levels[1]
    return(factor(out, levels = q_levels))
  }
  qs <- as.numeric(stats::quantile(d[ok], probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 7))
  eps <- max(.Machine$double.eps, diff(range(d[ok])) * 1e-9)
  if (qs[2] <= qs[1]) qs[2] <- qs[1] + eps
  if (qs[3] <= qs[2]) qs[3] <- qs[2] + eps
  br <- c(-Inf, qs[1], qs[2], qs[3], Inf)
  out[ok] <- as.character(cut(d[ok], breaks = br, labels = q_levels, include.lowest = TRUE, right = TRUE))
  factor(out, levels = q_levels)
}

get_counts <- function(seu, assay = "RNA") {
  DefaultAssay(seu) <- assay
  if (exists("JoinLayers", mode = "function")) {
    tryCatch(seu <- JoinLayers(seu, assay = assay), error = function(e) NULL)
  }
  out <- NULL
  if (exists("LayerData", mode = "function")) {
    out <- tryCatch(LayerData(seu, assay = assay, layer = "counts"), error = function(e) NULL)
  }
  if (is.null(out)) {
    out <- tryCatch(GetAssayData(seu, assay = assay, layer = "counts"), error = function(e) NULL)
  }
  if (is.null(out)) out <- tryCatch(seu[[assay]]$counts, error = function(e) NULL)
  if (is.null(out) || length(out) == 0) stop("Cannot extract counts")
  if (!inherits(out, "dgCMatrix")) out <- as(as(out, "CsparseMatrix"), "dgCMatrix")
  out
}

Build_Rank_from_DEG <- function(deg, eps = 1e-300) {
  pv <- deg$PValue
  pv[is.na(pv)] <- 1
  pv[pv < eps] <- eps
  stat <- sign(deg$logFC) * (-log10(pv))
  names(stat) <- deg$Gene
  stat <- stat[!is.na(stat)]
  stat <- tapply(stat, names(stat), function(x) x[which.max(abs(x))])
  sort(stat, decreasing = TRUE)
}

mouse_rank_to_human <- function(rank_vec, HM) {
  df <- data.frame(
    mouse = names(rank_vec),
    stat = as.numeric(rank_vec),
    stringsAsFactors = FALSE
  )
  mg <- merge(df, HM, by = "mouse")
  if (nrow(mg) == 0) return(numeric(0))
  stat <- tapply(mg$stat, mg$human, function(x) x[which.max(abs(x))])
  sort(stat, decreasing = TRUE)
}

write_rnk <- function(rank_vec, path) {
  df <- data.frame(
    gene = names(rank_vec),
    score = as.numeric(rank_vec),
    stringsAsFactors = FALSE
  )
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Paired edgeR: mid (Q2Q3) vs edge (Q1 or Q4)
# logFC > 0 = higher in Q2Q3
call_DEG_paired <- function(counts_mat, panels, group, mid_level = "Q2Q3",
                            edge_level, min_n = 3L) {
  keep <- group %in% c(mid_level, edge_level)
  counts <- round(as.matrix(counts_mat[, keep, drop = FALSE]))
  panels <- factor(panels[keep])
  group <- factor(group[keep], levels = c(edge_level, mid_level))
  # keep panels that have BOTH groups
  tab <- table(panels, group)
  ok_panels <- rownames(tab)[rowSums(tab > 0) == 2L]
  if (length(ok_panels) < min_n) {
    log_msg(
      "  insufficient paired panels for ", mid_level, " vs ", edge_level,
      ": n=", length(ok_panels), " (need >=", min_n, ")"
    )
    return(NULL)
  }
  keep2 <- as.character(panels) %in% ok_panels
  counts <- counts[, keep2, drop = FALSE]
  panels <- droplevels(factor(as.character(panels[keep2])))
  group <- droplevels(factor(as.character(group[keep2]), levels = c(edge_level, mid_level)))
  log_msg(
    "  edgeR ", mid_level, " vs ", edge_level,
    " | paired panels=", nlevels(panels),
    " | cols=", paste(sprintf("%s:%s", panels, group), collapse = ", ")
  )
  y <- DGEList(counts = counts)
  keep_g <- filterByExpr(y, group = group)
  y <- y[keep_g, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ panels + group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  # last coef = group mid vs edge
  qlf <- glmQLFTest(fit, coef = ncol(design))
  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  res <- data.frame(Gene = rownames(res), res, row.names = NULL)
  res$direction <- ifelse(res$logFC > 0, paste0("Up_", mid_level), paste0("Up_", edge_level))
  res$contrast <- paste0(mid_level, "_vs_", edge_level)
  res$n_panel <- nlevels(panels)
  res
}

# =============================================================================
log_msg("========== Rod INL quartile edgeR Q2Q3 vs edge BY AGE ==========")
log_msg("seurat : ", seurat_path)
log_msg("assign : ", assign_path)
log_msg("out    : ", out_dir)
log_msg("ages   : ", paste(sample1_keep, collapse = ", "), " (separate edgeR each)")

if (!file.exists(seurat_path)) stop("Seurat not found: ", seurat_path)
seu <- readRDS(seurat_path)
if (exists("UpdateSeuratObject", mode = "function")) {
  tryCatch(seu <- UpdateSeuratObject(seu), error = function(e) NULL)
}
meta <- seu@meta.data
meta$cell_id <- rownames(meta)
if (!"panel" %in% colnames(meta)) {
  if (all(c("sample1", "slice") %in% colnames(meta))) {
    meta$panel <- paste(meta$sample1, meta$slice, sep = "_")
  } else if (all(c("sample1", "rep") %in% colnames(meta))) {
    meta$panel <- paste(meta$sample1, meta$rep, sep = "_")
  } else stop("Need panel / sample1+slice/rep")
}
if (!"sample1" %in% colnames(meta)) meta$sample1 <- sub("_.*$", "", meta$panel)
meta$sample1 <- as.character(meta$sample1)
meta$panel <- as.character(meta$panel)

# ---- Rod quartile assignment ----
use_assign <- !isTRUE(force_reassign) && file.exists(assign_path)
if (use_assign) {
  log_msg("Reuse assignment: ", assign_path)
  rod_asg <- read.delim(assign_path, stringsAsFactors = FALSE)
} else {
  log_msg("Recompute Rod→INL quartiles (k=", k_inl, ") ...")
  meta_use <- meta[
    is.finite(meta$x) & is.finite(meta$y) &
      meta$celltype %in% c(rod_ct, inl_cts),
    ,
    drop = FALSE
  ]
  rod_asg_list <- list()
  for (pn in sort(unique(meta_use$panel))) {
    m <- meta_use[meta_use$panel == pn, , drop = FALSE]
    rod <- m[m$celltype == rod_ct, , drop = FALSE]
    inl <- m[m$celltype %in% inl_cts, , drop = FALSE]
    if (!nrow(rod) || !nrow(inl)) next
    dist <- rod_distance_to_inl(rod[, c("x", "y")], inl[, c("x", "y")], k = k_inl)
    qbin <- assign_distance_quartile(dist)
    rod_asg_list[[pn]] <- data.frame(
      cell_id = rod$cell_id, panel = pn,
      sample1 = as.character(rod$sample1),
      dist_to_inl = dist, dist_quartile = as.character(qbin),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rod_asg_list)) stop("No Rod quartile assignments")
  rod_asg <- do.call(rbind, rod_asg_list)
  rownames(rod_asg) <- NULL
  write.table(
    rod_asg, file.path(out_dir, "Rod_INL_distance_quartile_assignment.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
rod_asg$dist_quartile <- factor(as.character(rod_asg$dist_quartile), levels = q_levels)
rod_asg <- rod_asg[!is.na(rod_asg$dist_quartile), , drop = FALSE]
if (!"sample1" %in% colnames(rod_asg)) {
  rod_asg$sample1 <- sub("_.*$", "", rod_asg$panel)
}
rod_asg <- rod_asg[rod_asg$sample1 %in% sample1_keep, , drop = FALSE]
if (!nrow(rod_asg)) stop("No Rod cells after sample1 filter")

# group label: Q1 / Q2Q3 / Q4
rod_asg$group <- NA_character_
rod_asg$group[rod_asg$dist_quartile == q_levels[1]] <- "Q1"
rod_asg$group[rod_asg$dist_quartile %in% q_levels[2:3]] <- "Q2Q3"
rod_asg$group[rod_asg$dist_quartile == q_levels[4]] <- "Q4"
rod_asg <- rod_asg[!is.na(rod_asg$group), , drop = FALSE]
rod_asg$pb_id <- paste(rod_asg$panel, rod_asg$group, sep = "||")

n_tab <- as.data.frame(table(rod_asg$pb_id), stringsAsFactors = FALSE)
names(n_tab) <- c("pb_id", "n_cells")
keep_pb <- n_tab$pb_id[n_tab$n_cells >= min_cells_q]
rod_asg <- rod_asg[rod_asg$pb_id %in% keep_pb, , drop = FALSE]
log_msg(
  "pb units (panel×group) after min_cells_q=", min_cells_q, ": ",
  length(unique(rod_asg$pb_id))
)
write.table(
  n_tab, file.path(out_dir, "pseudobulk_n_cells_panel_x_group.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---- counts → pseudobulk ----
counts <- get_counts(seu, assay_use)
rm(seu); gc()
mt <- grep("^mt-", rownames(counts), ignore.case = TRUE)
if (length(mt)) {
  log_msg("Drop mt- genes: ", length(mt))
  counts <- counts[-mt, , drop = FALSE]
}
cells_use <- intersect(rod_asg$cell_id, colnames(counts))
rod_asg <- rod_asg[match(cells_use, rod_asg$cell_id), , drop = FALSE]
counts_use <- counts[, cells_use, drop = FALSE]

log_msg("Pseudobulk SUM (panel × group) ...")
pb_factor <- factor(rod_asg$pb_id, levels = sort(unique(rod_asg$pb_id)))
design_pb <- sparseMatrix(
  i = seq_len(length(pb_factor)),
  j = as.integer(pb_factor),
  x = 1,
  dims = c(length(pb_factor), nlevels(pb_factor))
)
colnames(design_pb) <- levels(pb_factor)
pb_mat <- as.matrix(counts_use %*% design_pb)
parts <- strsplit(colnames(pb_mat), "||", fixed = TRUE)
pb_meta <- data.frame(
  sample = colnames(pb_mat),
  panel = vapply(parts, `[`, character(1), 1L),
  group = vapply(parts, `[`, character(1), 2L),
  n_cells = as.integer(table(pb_factor))[colnames(pb_mat)],
  stringsAsFactors = FALSE
)
pb_meta$sample1 <- sub("_.*$", "", pb_meta$panel)
saveRDS(
  list(bulk_mat = pb_mat, bulk_meta = pb_meta),
  file.path(out_dir, "pseudobulk_Rod_panel_x_Qgroup.rds")
)
write.table(pb_meta, file.path(out_dir, "pseudobulk_meta.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Pseudobulk columns: ", ncol(pb_mat))

# ---- edgeR: each age (R1–R4) × three contrasts ----
# Naming A_vs_B: logFC>0 = higher in A (paired ~ panel + group)
contrasts <- list(
  Q2Q3_vs_Q1 = list(mid = "Q2Q3", edge = "Q1"),
  Q2Q3_vs_Q4 = list(mid = "Q2Q3", edge = "Q4"),
  Q1_vs_Q4 = list(mid = "Q1", edge = "Q4")
)

DEG <- list()
for (age in sample1_keep) {
  idx <- which(pb_meta$sample1 == age)
  if (!length(idx)) {
    log_msg("SKIP age ", age, ": no pseudobulk columns")
    next
  }
  mat_a <- pb_mat[, idx, drop = FALSE]
  meta_a <- pb_meta[idx, , drop = FALSE]
  log_msg("######## Age ", age, " | panels=", length(unique(meta_a$panel)),
          " | pb_cols=", ncol(mat_a), " ########")

  for (cn in names(contrasts)) {
    key <- paste(cn, age, sep = "__")
    log_msg("==== DEG: ", key, " ====")
    res <- call_DEG_paired(
      mat_a, meta_a$panel, meta_a$group,
      mid_level = contrasts[[cn]]$mid,
      edge_level = contrasts[[cn]]$edge,
      min_n = min_n
    )
    if (is.null(res)) next
    res$sample1 <- age
    res$contrast <- cn
    log_msg(sprintf(
      "  genes=%d | P<0.05 Up_Q2Q3=%d Up_%s=%d",
      nrow(res),
      sum(res$PValue < 0.05 & res$logFC > 0, na.rm = TRUE),
      contrasts[[cn]]$edge,
      sum(res$PValue < 0.05 & res$logFC < 0, na.rm = TRUE)
    ))
    DEG[[key]] <- res
    write.table(
      res,
      file.path(out_dir, "DEG_edgeR", paste0("DEG_edgeR_", key, "_Rod.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}
saveRDS(DEG, file.path(out_dir, "DEG_edgeR_Q2Q3_vs_edge_Rod_byAge.rds"))
if (!length(DEG)) stop("No DEG results (need >= min_n paired panels per age)")

# ---- ranks mouse → human ----
log_msg("Build ranks + map to human ...")
if (!file.exists(ortholog_path)) stop("Ortholog not found: ", ortholog_path)
load(ortholog_path)
if (!exists("HMZ_ortholog_combined")) stop("Expected HMZ_ortholog_combined")
HM <- unique(HMZ_ortholog_combined$HM[, c("human", "mouse")])
HM <- HM[HM$human != "" & HM$mouse != "" & !is.na(HM$human) & !is.na(HM$mouse), ]

rank_mouse <- list()
rank_human <- list()
for (key in names(DEG)) {
  rk_m <- Build_Rank_from_DEG(DEG[[key]])
  rk_h <- mouse_rank_to_human(rk_m, HM)
  rank_mouse[[key]] <- rk_m
  rank_human[[key]] <- rk_h
  write_rnk(rk_m, file.path(out_dir, "ranks_mouse", paste0("rank_mouse_", key, "_Rod.rnk")))
  write_rnk(rk_h, file.path(out_dir, "ranks_human", paste0("rank_human_", key, "_Rod.rnk")))
  log_msg(sprintf(
    "  %s: mouse=%d | human=%d",
    key, length(rk_m), length(rk_h)
  ))
}
saveRDS(rank_mouse, file.path(out_dir, "rank_mouse_Q2Q3_vs_edge_byAge.rds"))
saveRDS(rank_human, file.path(out_dir, "rank_human_Q2Q3_vs_edge_byAge.rds"))

# summary table
summ <- do.call(rbind, lapply(names(DEG), function(key) {
  r <- DEG[[key]]
  data.frame(
    key = key,
    sample1 = r$sample1[1],
    contrast = r$contrast[1],
    n_genes = nrow(r),
    n_P05_Up_Q2Q3 = sum(r$PValue < 0.05 & r$logFC > 0, na.rm = TRUE),
    n_P05_Up_edge = sum(r$PValue < 0.05 & r$logFC < 0, na.rm = TRUE),
    n_panel = r$n_panel[1],
    stringsAsFactors = FALSE
  )
}))
write.table(summ, file.path(out_dir, "DEG_summary_byAge.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

log_msg("========== DEG by age done ==========")
log_msg("Next: Rscript Rod_INL_quartile_GSEA_Q2Q3_vs_edge.R")
