#!/usr/bin/env Rscript
# edgeR DEG on cleaned spatial Seurat: R4 vs R1, per celltype.
# Pseudobulk unit = slice/panel (e.g. R1_rep1 … R4_rep4).
# Rank = sign(logFC) * (-log10(PValue)) → GSEA input (mouse + human).
# GSEA itself is a separate script: spatial_cleaned_GSEA_R4vsR1.R
#
# Example:
#   Rscript spatial_cleaned_DEG_R4vsR1.R
#
#   Rscript spatial_cleaned_DEG_R4vsR1.R \
#     --seurat /projects/hmz-aging/ProximityEffect/Spatial_slice_clean/R1_4_seurat_Merged_dim_RNA_2024_sliceCleaned_colored.rds \
#     --out_dir /projects/hmz-aging/ProximityEffect/Spatial_slice_clean/DEG_R4vsR1

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}

root <- "/projects/hmz-aging/ProximityEffect"
seurat_path <- get_arg(
  "--seurat",
  file.path(
    root, "Spatial_slice_clean",
    "R1_4_seurat_Merged_dim_RNA_2024_sliceCleaned_colored.rds"
  )
)
out_dir <- get_arg(
  "--out_dir",
  file.path(root, "Spatial_slice_clean", "DEG_R4vsR1")
)
ortholog_path <- get_arg(
  "--ortholog",
  "/projects/hmz-aging/Zebrafish_New_Datasets/HMZ_ortholog_combined"
)
min_n <- as.integer(get_arg("--min_n", "2"))

celltypes <- strsplit(
  get_arg("--celltypes", "AC,BC,Cone,HC,MG,RGC,RPE,Rod"),
  ",",
  fixed = TRUE
)[[1]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "ranks_mouse"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "ranks_human"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "DEG_edgeR"), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
  flush.console()
}

log_msg("========== spatial_cleaned_DEG_R4vsR1 start ==========")
log_msg("seurat : ", seurat_path)
log_msg("out_dir: ", out_dir)
log_msg("contrast: R4 (Old) vs R1 (Young); pseudobulk = slice/panel")
log_msg("celltypes: ", paste(celltypes, collapse = ", "))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

get_rna_counts <- function(seu) {
  if (!"RNA" %in% Assays(seu)) stop("Need RNA assay")
  DefaultAssay(seu) <- "RNA"
  if (exists("JoinLayers", mode = "function")) {
    tryCatch(seu <- JoinLayers(seu), error = function(e) NULL)
  }
  tryCatch(
    GetAssayData(seu, assay = "RNA", layer = "counts"),
    error = function(e) {
      tryCatch(
        GetAssayData(seu, assay = "RNA", slot = "counts"),
        error = function(e2) seu[["RNA"]]@counts
      )
    }
  )
}

resolve_panel <- function(meta) {
  if ("panel" %in% colnames(meta) && any(!is.na(meta$panel) & meta$panel != "")) {
    log_msg("Pseudobulk unit: panel")
    return(list(meta = meta, col = "panel"))
  }
  if ("slice" %in% colnames(meta) && "sample1" %in% colnames(meta)) {
    meta$panel <- paste(as.character(meta$sample1), as.character(meta$slice), sep = "_")
    log_msg("Pseudobulk unit: sample1_slice -> panel")
    return(list(meta = meta, col = "panel"))
  }
  if ("rep" %in% colnames(meta) && "sample1" %in% colnames(meta)) {
    meta$panel <- paste(as.character(meta$sample1), as.character(meta$rep), sep = "_")
    log_msg("Pseudobulk unit: sample1_rep -> panel")
    return(list(meta = meta, col = "panel"))
  }
  stop("Need meta$panel, or sample1+slice/rep for slice-level pseudobulk")
}

getPseudoBulkByCellPanel <- function(counts, meta, celltypes, panel_col = "panel") {
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  pb_list <- lapply(celltypes, function(ct) {
    cells_ct <- rownames(meta)[as.character(meta$celltype) == ct]
    cells_ct <- cells_ct[cells_ct %in% colnames(counts)]
    if (length(cells_ct) == 0) {
      log_msg("  skip empty: ", ct)
      return(NULL)
    }
    sub_counts <- counts[, cells_ct, drop = FALSE]
    sub_meta <- meta[cells_ct, , drop = FALSE]
    grp <- factor(sub_meta[[panel_col]], levels = sort(unique(sub_meta[[panel_col]])))
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- levels(grp)
    mat <- Matrix(sub_counts %*% Matrix(design, sparse = TRUE), sparse = TRUE)
    log_msg("  ", ct, ": cells=", length(cells_ct), " slices=", ncol(mat))
    mat
  })
  names(pb_list) <- celltypes
  pb_list[!vapply(pb_list, is.null, logical(1))]
}

# R4 vs R1; logFC > 0 = higher in R4
call_DEG_edgeR_R4vsR1 <- function(bulk_mat, bulk_meta, min_n = 2) {
  stopifnot(identical(colnames(bulk_mat), as.character(bulk_meta$sample)))
  keep <- bulk_meta$sample1 %in% c("R1", "R4")
  counts <- round(as.matrix(bulk_mat[, keep, drop = FALSE]))
  group <- factor(
    ifelse(bulk_meta$sample1[keep] == "R1", "R1", "R4"),
    levels = c("R1", "R4")
  )
  tab <- table(group)
  if (nlevels(droplevels(group)) < 2 || min(tab) < min_n) {
    log_msg(
      "  insufficient slices: ",
      paste(names(tab), as.integer(tab), collapse = ", ")
    )
    return(NULL)
  }
  log_msg(
    "  edgeR R4 vs R1 | slices: ",
    paste(names(tab), as.integer(tab), collapse = ", ")
  )
  y <- DGEList(counts = counts, group = group)
  keep_g <- filterByExpr(y, group = group)
  y <- y[keep_g, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~group) # intercept + R4
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  res <- data.frame(Gene = rownames(res), res, row.names = NULL)
  res$direction <- ifelse(res$logFC > 0, "Up_R4", "Up_R1")
  res
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
  # GSEA .rnk: gene \t score
  df <- data.frame(
    gene = names(rank_vec),
    score = as.numeric(rank_vec),
    stringsAsFactors = FALSE
  )
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ---------------------------------------------------------------------------
# 1. Load cleaned Seurat
# ---------------------------------------------------------------------------

if (!file.exists(seurat_path)) stop("Seurat not found: ", seurat_path)

log_msg("Loading cleaned Seurat ...")
t0 <- Sys.time()
seu <- readRDS(seurat_path)
log_msg(
  "  readRDS in ",
  signif(as.numeric(difftime(Sys.time(), t0, units = "mins")), 3),
  " min | cells=", ncol(seu)
)
if (!inherits(seu, "Seurat")) stop("Not a Seurat object")
if (exists("UpdateSeuratObject", mode = "function")) {
  tryCatch(seu <- UpdateSeuratObject(seu), error = function(e) NULL)
}
if (!"celltype" %in% colnames(seu@meta.data)) stop("Need celltype")
if (!"sample1" %in% colnames(seu@meta.data)) stop("Need sample1")

# keep R1/R4 only
seu <- subset(seu, subset = sample1 %in% c("R1", "R4"))
log_msg("After R1/R4 filter: cells=", ncol(seu))
print(table(seu$sample1, seu$celltype))

meta <- seu@meta.data
meta$cell_id <- rownames(meta)
pr <- resolve_panel(meta)
meta <- pr$meta
panel_col <- pr$col
seu$panel <- meta$panel

present_ct <- unique(as.character(seu$celltype))
celltypes <- intersect(celltypes, present_ct)
if (!length(celltypes)) stop("No requested celltypes present")
log_msg("Panels: ", paste(sort(unique(meta$panel)), collapse = ", "))

# ---------------------------------------------------------------------------
# 2. Pseudobulk by celltype × panel (RNA raw counts)
# ---------------------------------------------------------------------------

log_msg("Extract RNA raw counts ...")
counts <- get_rna_counts(seu)
mt <- grep("^mt-", rownames(counts), ignore.case = TRUE)
if (length(mt) > 0) {
  log_msg("Drop mt- genes: ", length(mt))
  counts <- counts[-mt, , drop = FALSE]
}

log_msg("Build slice-level pseudobulk ...")
pb_list <- getPseudoBulkByCellPanel(counts, meta, celltypes, panel_col = panel_col)

pb_input <- list()
for (ct in names(pb_list)) {
  mat <- as.matrix(pb_list[[ct]])
  units <- colnames(mat)
  unit_to_s1 <- tapply(as.character(meta$sample1), meta[[panel_col]], function(x) {
    ux <- unique(x)
    if (length(ux) != 1) stop("Mixed sample1 in panel")
    ux
  })
  s1 <- unname(unit_to_s1[units])
  if (anyNA(s1)) {
    s1 <- sub("_.*$", "", units)
  }
  bulk_meta <- data.frame(
    sample = units,
    sample1 = s1,
    stringsAsFactors = FALSE
  )
  # keep only R1/R4 columns
  keep_col <- bulk_meta$sample1 %in% c("R1", "R4")
  mat <- mat[, keep_col, drop = FALSE]
  bulk_meta <- bulk_meta[keep_col, , drop = FALSE]
  pb_input[[ct]] <- list(bulk_mat = mat, bulk_meta = bulk_meta)
  log_msg(
    "  ", ct, " panels: ",
    paste(bulk_meta$sample, collapse = ", ")
  )
}
saveRDS(pb_input, file.path(out_dir, "pseudobulk_by_slice_R1R4.rds"))

rm(seu, counts)
gc()

# ---------------------------------------------------------------------------
# 3. edgeR R4 vs R1 per celltype
# ---------------------------------------------------------------------------

log_msg("edgeR R4 vs R1 per celltype ...")
DEG <- list()
for (ct in names(pb_input)) {
  log_msg("==== DEG: ", ct, " ====")
  res <- call_DEG_edgeR_R4vsR1(
    pb_input[[ct]]$bulk_mat,
    pb_input[[ct]]$bulk_meta,
    min_n = min_n
  )
  if (is.null(res)) next
  log_msg(sprintf(
    "  %s: genes=%d | P<0.05 Up_R4=%d Up_R1=%d",
    ct, nrow(res),
    sum(res$PValue < 0.05 & res$logFC > 0, na.rm = TRUE),
    sum(res$PValue < 0.05 & res$logFC < 0, na.rm = TRUE)
  ))
  DEG[[ct]] <- res
  write.table(
    res,
    file.path(out_dir, "DEG_edgeR", paste0("DEG_edgeR_R4vsR1_", ct, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
saveRDS(DEG, file.path(out_dir, "DEG_edgeR_R4vsR1.rds"))
if (!length(DEG)) stop("No DEG results; check slice replicates")

# ---------------------------------------------------------------------------
# 4. Rank (mouse) → GSEA input
# ---------------------------------------------------------------------------

log_msg("Build mouse ranks: sign(logFC)*(-log10(PValue)) ...")
rank_mouse <- list()
for (ct in names(DEG)) {
  rk <- Build_Rank_from_DEG(DEG[[ct]])
  rank_mouse[[ct]] <- rk
  write_rnk(rk, file.path(out_dir, "ranks_mouse", paste0("rank_mouse_", ct, ".rnk")))
  # also tsv with header for convenience
  write.table(
    data.frame(gene = names(rk), score = as.numeric(rk)),
    file.path(out_dir, "ranks_mouse", paste0("rank_mouse_", ct, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  log_msg(sprintf("  %s: n_rank=%d", ct, length(rk)))
}
saveRDS(rank_mouse, file.path(out_dir, "rank_mouse_R4vsR1.rds"))

# ---------------------------------------------------------------------------
# 5. Mouse → human ortholog ranks (for human_aging_catalog GSEA)
# ---------------------------------------------------------------------------

log_msg("Map ranks mouse -> human ...")
if (!file.exists(ortholog_path)) stop("Ortholog not found: ", ortholog_path)
load(ortholog_path)
if (!exists("HMZ_ortholog_combined")) stop("Expected HMZ_ortholog_combined")
HM <- unique(HMZ_ortholog_combined$HM[, c("human", "mouse")])
HM <- HM[HM$human != "" & HM$mouse != "" & !is.na(HM$human) & !is.na(HM$mouse), ]

rank_human <- lapply(rank_mouse, mouse_rank_to_human, HM = HM)
names(rank_human) <- names(rank_mouse)

map_summary <- data.frame(
  celltype = names(rank_mouse),
  mouse = vapply(rank_mouse, length, integer(1)),
  human = vapply(rank_human, length, integer(1)),
  human_nonzero = vapply(rank_human, function(x) sum(x != 0), integer(1)),
  stringsAsFactors = FALSE
)
print(map_summary)
write.table(
  map_summary,
  file.path(out_dir, "ortholog_rank_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

for (ct in names(rank_human)) {
  rk <- rank_human[[ct]]
  write_rnk(rk, file.path(out_dir, "ranks_human", paste0("rank_human_", ct, ".rnk")))
  write.table(
    data.frame(gene = names(rk), score = as.numeric(rk)),
    file.path(out_dir, "ranks_human", paste0("rank_human_", ct, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
saveRDS(rank_human, file.path(out_dir, "rank_human_R4vsR1.rds"))
log_msg("GSEA rank inputs written under ranks_mouse/ and ranks_human/")
log_msg("Next: Rscript spatial_cleaned_GSEA_R4vsR1.R --rank_dir ", out_dir)
log_msg("========== spatial_cleaned_DEG_R4vsR1 done ==========")
