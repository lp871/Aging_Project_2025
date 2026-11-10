#!/usr/bin/env Rscript
# Per-slice Rod → INL distance quartiles → pseudobulk → aging score.
#
# Pipeline:
#   1) Each panel: Rod→INL distance = median of k nearest INL (BC/AC/HC/MG)
#   2) Within panel: Q1–Q4 (0–25% closest … 75–100% farthest)
#   3) Pseudobulk = SUM raw RNA counts per panel × quartile
#   4) NormalizeData + AddModuleScore(Rod Up/Down gene lists)
#      aging_score = Up - Down
#   5) Print each slice's Q1–Q4 aging scores
#
# Example:
#   Rscript Rod_INL_quartile_pseudobulk_aging_score.R
#   Rscript Rod_INL_quartile_pseudobulk_aging_score.R --k_inl 50 --gene_set exactly3,all4

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(FNN)
  library(ggplot2)
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
    "Rod_INL_quartile_pseudobulk_aging_score"
  )
)

inl_cts <- strsplit(get_arg("--inl_cts", "BC,AC,HC,MG"), ",", fixed = TRUE)[[1]]
rod_ct <- get_arg("--rod_ct", "Rod")
k_inl <- as.integer(get_arg("--k_inl", "50"))
gene_sets <- strsplit(get_arg("--gene_set", "exactly3,all4"), ",", fixed = TRUE)[[1]]
min_cells_q <- as.integer(get_arg("--min_cells_q", "10"))
min_genes <- as.integer(get_arg("--min_genes", "3"))
assay_use <- get_arg("--assay", "RNA")
dpi <- as.integer(get_arg("--dpi", "300"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
gl_dir <- file.path(out_dir, "genelists_used")
dir.create(gl_dir, recursive = TRUE, showWarnings = FALSE)

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

# median distance to k nearest INL cells
rod_distance_to_inl <- function(xy_rod, xy_inl, k = 50L) {
  xy_rod <- as.matrix(xy_rod)
  xy_inl <- as.matrix(xy_inl)
  storage.mode(xy_rod) <- "double"
  storage.mode(xy_inl) <- "double"
  n_rod <- nrow(xy_rod)
  n_inl <- nrow(xy_inl)
  if (!n_rod) return(numeric(0))
  if (!n_inl) return(rep(Inf, n_rod))
  k_use <- min(as.integer(k), n_inl)
  if (k_use < 1L) return(rep(Inf, n_rod))
  nn <- FNN::get.knnx(data = xy_inl, query = xy_rod, k = k_use)
  if (k_use == 1L) {
    as.numeric(nn$nn.dist[, 1])
  } else {
    as.numeric(apply(nn$nn.dist, 1L, stats::median, na.rm = TRUE))
  }
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

log_msg("========== Rod INL-quartile pseudobulk aging score ==========")
log_msg("seurat : ", seurat_path)
log_msg("genes  : ", genes_long_path)
log_msg("out    : ", out_dir)
log_msg("INL    : ", paste(inl_cts, collapse = ", "))
log_msg("k_inl  : ", k_inl, " (median of k nearest INL distances)")
log_msg("sets   : ", paste(gene_sets, collapse = ", "))

if (!file.exists(seurat_path)) stop("Seurat not found: ", seurat_path)
if (!file.exists(genes_long_path)) stop("genes_long not found: ", genes_long_path)

# gene lists (Rod only)
glong <- read.delim(genes_long_path, stringsAsFactors = FALSE)
need_g <- c("celltype", "Gene", "agree_level", "direction")
miss_g <- setdiff(need_g, colnames(glong))
if (length(miss_g)) stop("genes_long missing: ", paste(miss_g, collapse = ", "))

genelists <- list()
for (gs in gene_sets) {
  for (d in c("Up", "Down")) {
    key <- paste(rod_ct, gs, d, sep = "__")
    genes <- unique(glong$Gene[
      glong$celltype == rod_ct &
        glong$agree_level == gs &
        glong$direction == d
    ])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    genelists[[key]] <- genes
    writeLines(genes, file.path(gl_dir, paste0(key, ".txt")))
    log_msg("genelist ", key, ": n=", length(genes))
  }
}

# load seurat
log_msg("Load Seurat ...")
seu <- readRDS(seurat_path)
if (exists("UpdateSeuratObject", mode = "function")) {
  tryCatch(seu <- UpdateSeuratObject(seu), error = function(e) NULL)
}
meta <- seu@meta.data
meta$cell_id <- rownames(meta)
need <- c("celltype", "x", "y")
miss <- setdiff(need, colnames(meta))
if (length(miss)) stop("meta missing: ", paste(miss, collapse = ", "))
if (!"panel" %in% colnames(meta)) {
  if (all(c("sample1", "slice") %in% colnames(meta))) {
    meta$panel <- paste(meta$sample1, meta$slice, sep = "_")
  } else if (all(c("sample1", "rep") %in% colnames(meta))) {
    meta$panel <- paste(meta$sample1, meta$rep, sep = "_")
  } else {
    stop("Need panel")
  }
}
if (!"sample1" %in% colnames(meta)) {
  meta$sample1 <- sub("_.*$", "", as.character(meta$panel))
}
meta$celltype <- as.character(meta$celltype)
meta$panel <- as.character(meta$panel)
meta <- meta[
  is.finite(meta$x) & is.finite(meta$y) &
    meta$celltype %in% c(rod_ct, inl_cts),
  ,
  drop = FALSE
]

counts <- get_counts(seu, assay_use)
common <- intersect(meta$cell_id, colnames(counts))
meta <- meta[common, , drop = FALSE]
counts <- counts[, common, drop = FALSE]
rm(seu)
gc()

# assign Rod quartiles per panel
panels <- sort(unique(meta$panel))
rod_asg_list <- list()
for (pn in panels) {
  m <- meta[meta$panel == pn, , drop = FALSE]
  rod <- m[m$celltype == rod_ct, , drop = FALSE]
  inl <- m[m$celltype %in% inl_cts, , drop = FALSE]
  if (!nrow(rod) || !nrow(inl)) {
    log_msg("SKIP assign ", pn, ": Rod=", nrow(rod), " INL=", nrow(inl))
    next
  }
  dist <- rod_distance_to_inl(rod[, c("x", "y")], inl[, c("x", "y")], k = k_inl)
  qbin <- assign_distance_quartile(dist)
  rod_asg_list[[pn]] <- data.frame(
    cell_id = rod$cell_id,
    panel = pn,
    sample1 = as.character(rod$sample1),
    dist_to_inl = dist,
    dist_quartile = qbin,
    stringsAsFactors = FALSE
  )
}
if (!length(rod_asg_list)) stop("No Rod quartile assignments")
rod_asg <- do.call(rbind, rod_asg_list)
rownames(rod_asg) <- NULL
rod_asg$dist_quartile <- factor(as.character(rod_asg$dist_quartile), levels = q_levels)
rod_asg$pb_id <- paste(rod_asg$panel, as.character(rod_asg$dist_quartile), sep = "||")

# filter small quartile groups
n_tab <- as.data.frame(table(rod_asg$pb_id), stringsAsFactors = FALSE)
names(n_tab) <- c("pb_id", "n_cells")
keep_pb <- n_tab$pb_id[n_tab$n_cells >= min_cells_q]
rod_asg <- rod_asg[rod_asg$pb_id %in% keep_pb, , drop = FALSE]
log_msg(
  "pseudobulk units (panel×Q): ", length(unique(rod_asg$pb_id)),
  " | min_cells_q=", min_cells_q
)
if (!nrow(rod_asg)) stop("No pb units after min_cells filter")

write.table(
  rod_asg,
  file.path(out_dir, "Rod_INL_distance_quartile_assignment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# pseudobulk SUM counts
log_msg("Pseudobulk SUM raw counts (panel × quartile) ...")
cells_use <- intersect(rod_asg$cell_id, colnames(counts))
rod_asg <- rod_asg[match(cells_use, rod_asg$cell_id), , drop = FALSE]
counts_use <- counts[, cells_use, drop = FALSE]
pb_factor <- factor(rod_asg$pb_id, levels = sort(unique(rod_asg$pb_id)))
design <- sparseMatrix(
  i = seq_len(length(pb_factor)),
  j = as.integer(pb_factor),
  x = 1,
  dims = c(length(pb_factor), nlevels(pb_factor))
)
colnames(design) <- levels(pb_factor)
pb_mat <- counts_use %*% design
colnames(pb_mat) <- levels(pb_factor)

parts <- strsplit(colnames(pb_mat), "||", fixed = TRUE)
pb_meta <- data.frame(
  pb_id = colnames(pb_mat),
  panel = vapply(parts, `[`, character(1), 1L),
  dist_quartile = vapply(parts, `[`, character(1), 2L),
  n_cells = as.integer(table(pb_factor))[colnames(pb_mat)],
  stringsAsFactors = FALSE
)
pb_meta$sample1 <- sub("_.*$", "", pb_meta$panel)
pb_meta$dist_quartile <- factor(pb_meta$dist_quartile, levels = q_levels)
rownames(pb_meta) <- pb_meta$pb_id
pb_meta <- pb_meta[colnames(pb_mat), , drop = FALSE]

saveRDS(
  list(counts = pb_mat, meta = pb_meta),
  file.path(out_dir, "Rod_quartile_pseudobulk_raw.rds")
)

# Normalize + AddModuleScore
log_msg("NormalizeData + AddModuleScore ...")
op <- options(Seurat.object.assay.version = "v3")
on.exit(options(op), add = TRUE)
seu_pb <- CreateSeuratObject(
  counts = pb_mat,
  meta.data = pb_meta,
  assay = "RNA",
  project = "RodQ"
)
seu_pb <- NormalizeData(
  seu_pb,
  normalization.method = "LogNormalize",
  verbose = FALSE
)
features_all <- rownames(seu_pb)

score_rows <- list()
for (gs in gene_sets) {
  up_genes <- intersect(genelists[[paste(rod_ct, gs, "Up", sep = "__")]], features_all)
  dn_genes <- intersect(genelists[[paste(rod_ct, gs, "Down", sep = "__")]], features_all)
  log_msg("  gene_set=", gs, " Up=", length(up_genes), " Down=", length(dn_genes))

  up_score <- rep(NA_real_, ncol(seu_pb))
  dn_score <- rep(NA_real_, ncol(seu_pb))
  names(up_score) <- colnames(seu_pb)
  names(dn_score) <- colnames(seu_pb)

  if (length(up_genes) >= min_genes) {
    nm <- gsub("[^A-Za-z0-9_]", "_", paste0("Up_", gs))
    created <- paste0(nm, "1")
    if (created %in% colnames(seu_pb@meta.data)) seu_pb[[created]] <- NULL
    seu_pb <- AddModuleScore(
      seu_pb, features = list(up_genes), name = nm,
      assay = "RNA", search = FALSE
    )
    up_score <- seu_pb@meta.data[[created]]
    names(up_score) <- colnames(seu_pb)
    seu_pb[[created]] <- NULL
  }
  if (length(dn_genes) >= min_genes) {
    nm <- gsub("[^A-Za-z0-9_]", "_", paste0("Down_", gs))
    created <- paste0(nm, "1")
    if (created %in% colnames(seu_pb@meta.data)) seu_pb[[created]] <- NULL
    seu_pb <- AddModuleScore(
      seu_pb, features = list(dn_genes), name = nm,
      assay = "RNA", search = FALSE
    )
    dn_score <- seu_pb@meta.data[[created]]
    names(dn_score) <- colnames(seu_pb)
    seu_pb[[created]] <- NULL
  }

  aging <- up_score - dn_score
  aging[is.na(up_score) & !is.na(dn_score)] <- -dn_score[is.na(up_score) & !is.na(dn_score)]
  aging[!is.na(up_score) & is.na(dn_score)] <- up_score[!is.na(up_score) & is.na(dn_score)]

  score_rows[[gs]] <- data.frame(
    pb_id = colnames(seu_pb),
    panel = seu_pb$panel,
    sample1 = seu_pb$sample1,
    dist_quartile = as.character(seu_pb$dist_quartile),
    n_cells = seu_pb$n_cells,
    gene_set = gs,
    Up_score = as.numeric(up_score),
    Down_score = as.numeric(dn_score),
    aging_score = as.numeric(aging),
    stringsAsFactors = FALSE
  )
}

score_df <- do.call(rbind, score_rows)
rownames(score_df) <- NULL
score_df$dist_quartile <- factor(score_df$dist_quartile, levels = q_levels)
score_df <- score_df[order(score_df$gene_set, score_df$panel, score_df$dist_quartile), ]

out_score <- file.path(out_dir, "Rod_quartile_pseudobulk_aging_score.tsv")
write.table(score_df, out_score, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote ", out_score)

# ---------------------------------------------------------------------------
# Print each slice's quartile aging scores
# ---------------------------------------------------------------------------

log_msg("")
log_msg("========== Aging score by slice × Rod→INL quartile ==========")
for (gs in gene_sets) {
  log_msg("---- gene_set: ", gs, " ----")
  d <- score_df[score_df$gene_set == gs, , drop = FALSE]
  for (pn in sort(unique(d$panel))) {
    dp <- d[d$panel == pn, , drop = FALSE]
    dp <- dp[order(dp$dist_quartile), , drop = FALSE]
    log_msg("  [", pn, "]")
    for (i in seq_len(nrow(dp))) {
      log_msg(
        sprintf(
          "    %-22s  n=%5d  aging_score=%+.5f  (Up=%+.4f Down=%+.4f)",
          as.character(dp$dist_quartile[i]),
          as.integer(dp$n_cells[i]),
          dp$aging_score[i],
          dp$Up_score[i],
          dp$Down_score[i]
        )
      )
    }
    # compact one-liner
    vals <- setNames(dp$aging_score, as.character(dp$dist_quartile))
    log_msg(
      "    summary: ",
      paste(sprintf("%s=%+.4f", names(vals), vals), collapse = " | ")
    )
  }
}
log_msg("==============================================================")

# wide table for easy reading
for (gs in gene_sets) {
  d <- score_df[score_df$gene_set == gs, , drop = FALSE]
  wide <- reshape(
    d[, c("panel", "dist_quartile", "aging_score", "n_cells")],
    idvar = "panel",
    timevar = "dist_quartile",
    direction = "wide"
  )
  write.table(
    wide,
    file.path(out_dir, paste0("Rod_quartile_aging_score_wide_", gs, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# plot: aging score ~ quartile, facet panel
q_pal <- c(
  "Q1_0-25%(closest)" = "#313695",
  "Q2_25-50%" = "#74ADD1",
  "Q3_50-75%" = "#FDAE61",
  "Q4_75-100%(farthest)" = "#A50026"
)
for (gs in gene_sets) {
  d <- score_df[score_df$gene_set == gs, , drop = FALSE]
  if (!nrow(d)) next
  p <- ggplot(d, aes(x = dist_quartile, y = aging_score, fill = dist_quartile)) +
    geom_col(width = 0.7, color = NA) +
    facet_wrap(~panel, scales = "free_y") +
    scale_fill_manual(values = q_pal, drop = FALSE) +
    labs(
      title = paste0("Rod pseudobulk aging score by INL-distance quartile | ", gs),
      subtitle = paste0(
        "dist=median of ", k_inl, " nearest INL; Q1=closest … Q4=farthest; ",
        "aging = Up-Down"
      ),
      x = NULL, y = "aging score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey35")
    )
  n_pn <- length(unique(d$panel))
  fp <- file.path(fig_dir, paste0("Rod_quartile_aging_score_bySlice_", gs, ".png"))
  ggsave(
    fp, p,
    width = max(10, 2.5 * ceiling(sqrt(n_pn))),
    height = max(7, 2.2 * ceiling(n_pn / ceiling(sqrt(n_pn)))),
    dpi = dpi
  )
  log_msg("Wrote ", fp)
}

log_msg("========== done ==========")
