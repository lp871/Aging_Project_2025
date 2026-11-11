#!/usr/bin/env Rscript
# edgeR DEG on cleaned OSK Seurat: OSKR4 vs OSKR3, per celltype.
# Pseudobulk unit = slice/panel (e.g. OSKR3_rep1 … OSKR4_rep4).
# Rank = sign(logFC) * (-log10(PValue)) → GSEA input (mouse + human).
# GSEA: OSK_R4vsR3_fgsea_bubble_2026.R
#
# Ported from Proximity_final/spatial_cleaned_DEG_R4vsR1.R
#
# HPC:
#   conda activate seurat_env
#   cd /projects/hmz-aging/Xenium5k_Mouse/OSK_New
#   Rscript OSK_R4vsR3_DEG_rank_2026.R
#
#   Rscript OSK_R4vsR3_DEG_rank_2026.R \
#     --seurat OSK_Spatial_slice_clean/OSK_obj_cl_2026_sliceCleaned_OSKR3R4.rds \
#     --out_dir OSK_R4vsR3_Aging_GSEA_2026

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

root <- "/projects/hmz-aging/Xenium5k_Mouse/OSK_New"
seurat_path <- get_arg(
  "--seurat",
  file.path(
    root, "OSK_Spatial_slice_clean",
    "OSK_obj_cl_2026_sliceCleaned_OSKR3R4.rds"
  )
)
out_dir <- get_arg(
  "--out_dir",
  file.path(root, "OSK_R4vsR3_Aging_GSEA_2026")
)
ortholog_path <- get_arg(
  "--ortholog",
  "/projects/hmz-aging/Zebrafish_New_Datasets/HMZ_ortholog_combined"
)
min_n <- as.integer(get_arg("--min_n", "2"))

sample_r4 <- get_arg("--sample_r4", "OSKR4")
sample_r3 <- get_arg("--sample_r3", "OSKR3")

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

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- gsub("^R([1-4])$", "OSKR\\1", x)
  x
}

log_msg("========== OSK_R4vsR3_DEG_rank (edgeR) start ==========")
log_msg("seurat : ", seurat_path)
log_msg("out_dir: ", out_dir)
log_msg("contrast: ", sample_r4, " vs ", sample_r3, "; pseudobulk = slice/panel")
log_msg("celltypes: ", paste(celltypes, collapse = ", "))

# ---------------------------------------------------------------------------
# Helpers (from Proximity_final)
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
  s1 <- if ("sample1" %in% colnames(meta)) {
    as.character(meta$sample1)
  } else if ("sample" %in% colnames(meta)) {
    as.character(meta$sample)
  } else {
    stop("Need sample1 or sample")
  }
  if ("slice" %in% colnames(meta) && any(!is.na(meta$slice) & meta$slice != "")) {
    meta$panel <- paste(s1, as.character(meta$slice), sep = "_")
    log_msg("Pseudobulk unit: sample_slice -> panel")
    return(list(meta = meta, col = "panel"))
  }
  if ("rep" %in% colnames(meta) && any(!is.na(meta$rep) & meta$rep != "")) {
    meta$panel <- paste(s1, as.character(meta$rep), sep = "_")
    log_msg("Pseudobulk unit: sample_rep -> panel")
    return(list(meta = meta, col = "panel"))
  }
  stop("Need meta$panel, or sample(+1)+slice/rep for slice-level pseudobulk")
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

## OSKR4 vs OSKR3; logFC > 0 = higher in OSKR4
call_DEG_edgeR_R4vsR3 <- function(bulk_mat, bulk_meta, min_n = 2,
                                  s_r4 = "OSKR4", s_r3 = "OSKR3") {
  stopifnot(identical(colnames(bulk_mat), as.character(bulk_meta$sample)))
  keep <- bulk_meta$sample1 %in% c(s_r3, s_r4)
  counts <- round(as.matrix(bulk_mat[, keep, drop = FALSE]))
  group <- factor(
    ifelse(bulk_meta$sample1[keep] == s_r3, s_r3, s_r4),
    levels = c(s_r3, s_r4)
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
    "  edgeR ", s_r4, " vs ", s_r3, " | slices: ",
    paste(names(tab), as.integer(tab), collapse = ", ")
  )
  y <- DGEList(counts = counts, group = group)
  keep_g <- filterByExpr(y, group = group)
  y <- y[keep_g, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~group) # intercept + OSKR4
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  res <- data.frame(Gene = rownames(res), res, row.names = NULL)
  res$direction <- ifelse(res$logFC > 0, paste0("Up_", s_r4), paste0("Up_", s_r3))
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
  df <- data.frame(
    gene = names(rank_vec),
    score = as.numeric(rank_vec),
    stringsAsFactors = FALSE
  )
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

resolve_seurat <- function(path) {
  if (file.exists(path)) return(path)
  cands <- c(
    path,
    file.path(root, "OSK_Spatial_slice_clean",
              "OSK_obj_cl_2026_sliceCleaned_OSKR3R4.rds"),
    file.path(root, "OSK_Spatial_slice_clean",
              "OSK_obj_cl_2026_sliceCleaned_colored.rds"),
    file.path(root, "OSK_obj_cl_2026_withSlice.rds"),
    file.path(root, "OSK_obj_cl_2026.rds")
  )
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("Seurat not found: ", path)
  hit
}

# ---------------------------------------------------------------------------
# 1. Load cleaned Seurat
# ---------------------------------------------------------------------------

seurat_path <- resolve_seurat(seurat_path)
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

## unify sample1
if (!"sample1" %in% colnames(seu@meta.data)) {
  if ("sample" %in% colnames(seu@meta.data)) {
    seu$sample1 <- normalize_sample_id(seu$sample)
  } else {
    stop("Need sample or sample1")
  }
} else {
  seu$sample1 <- normalize_sample_id(seu$sample1)
}

seu <- subset(seu, subset = sample1 %in% c(sample_r3, sample_r4))
log_msg("After ", sample_r3, "/", sample_r4, " filter: cells=", ncol(seu))
print(table(seu$sample1, seu$celltype))

meta <- seu@meta.data
meta$cell_id <- rownames(meta)
meta$sample1 <- as.character(meta$sample1)
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
if (exists("JoinLayers", mode = "function")) {
  tryCatch(seu <- JoinLayers(seu, assay = "RNA"), error = function(e) NULL)
}
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
  if (anyNA(s1)) s1 <- sub("_.*$", "", units)
  s1 <- normalize_sample_id(s1)
  bulk_meta <- data.frame(
    sample = units,
    sample1 = s1,
    stringsAsFactors = FALSE
  )
  keep_col <- bulk_meta$sample1 %in% c(sample_r3, sample_r4)
  mat <- mat[, keep_col, drop = FALSE]
  bulk_meta <- bulk_meta[keep_col, , drop = FALSE]
  pb_input[[ct]] <- list(bulk_mat = mat, bulk_meta = bulk_meta)
  log_msg("  ", ct, " panels: ", paste(bulk_meta$sample, collapse = ", "))
}
saveRDS(pb_input, file.path(out_dir, "pseudobulk_by_slice_OSKR3R4.rds"))

rm(seu, counts)
gc()

# ---------------------------------------------------------------------------
# 3. edgeR OSKR4 vs OSKR3 per celltype
# ---------------------------------------------------------------------------

log_msg("edgeR ", sample_r4, " vs ", sample_r3, " per celltype ...")
DEG <- list()
for (ct in names(pb_input)) {
  log_msg("==== DEG: ", ct, " ====")
  res <- call_DEG_edgeR_R4vsR3(
    pb_input[[ct]]$bulk_mat,
    pb_input[[ct]]$bulk_meta,
    min_n = min_n,
    s_r4 = sample_r4,
    s_r3 = sample_r3
  )
  if (is.null(res)) next
  log_msg(sprintf(
    "  %s: genes=%d | P<0.05 Up_%s=%d Up_%s=%d",
    ct, nrow(res),
    sample_r4, sum(res$PValue < 0.05 & res$logFC > 0, na.rm = TRUE),
    sample_r3, sum(res$PValue < 0.05 & res$logFC < 0, na.rm = TRUE)
  ))
  DEG[[ct]] <- res
  write.table(
    res,
    file.path(out_dir, "DEG_edgeR", paste0("DEG_edgeR_OSKR4vsOSKR3_", ct, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
saveRDS(DEG, file.path(out_dir, "DEG_edgeR_OSKR4vsOSKR3.rds"))
if (!length(DEG)) stop("No DEG results; check slice replicates (need >= ", min_n, " per group)")

# ---------------------------------------------------------------------------
# 4. Rank (mouse)
# ---------------------------------------------------------------------------

log_msg("Build mouse ranks: sign(logFC)*(-log10(PValue)) ...")
rank_mouse <- list()
for (ct in names(DEG)) {
  rk <- Build_Rank_from_DEG(DEG[[ct]])
  rank_mouse[[ct]] <- rk
  write_rnk(rk, file.path(out_dir, "ranks_mouse", paste0("rank_mouse_", ct, ".rnk")))
  write.table(
    data.frame(gene = names(rk), score = as.numeric(rk)),
    file.path(out_dir, "ranks_mouse", paste0("rank_mouse_", ct, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  log_msg(sprintf("  %s: n_rank=%d", ct, length(rk)))
}
saveRDS(rank_mouse, file.path(out_dir, "rank_mouse_OSKR4vsOSKR3.rds"))

# ---------------------------------------------------------------------------
# 5. Mouse → human
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
saveRDS(rank_human, file.path(out_dir, "rank_human_OSKR4vsOSKR3.rds"))
## alias for GSEA script / older name
saveRDS(rank_human, file.path(out_dir, "OSK_R4vsR3_DEG_rank_human.rds"))

log_msg("GSEA rank inputs under ranks_mouse/ and ranks_human/")
log_msg("Next: conda activate clusterProfiler")
log_msg("  Rscript OSK_R4vsR3_fgsea_bubble_2026.R")
log_msg("========== OSK_R4vsR3_DEG_rank (edgeR) done ==========")
