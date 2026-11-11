#!/usr/bin/env Rscript
# fgsea on OSK OSKR4vsOSKR3 human ranks (from OSK_R4vsR3_DEG_rank_2026.R).
# Default: C1–C10 category_unions from human_aging_catalog.
# Optional: --run_terms for individual pathway terms.
#
# Ported from Proximity_final/spatial_cleaned_GSEA_R4vsR1.R
#
# Prerequisites:
#   conda activate seurat_env
#   Rscript OSK_R4vsR3_DEG_rank_2026.R
#   → writes rank_human_OSKR4vsOSKR3.rds
#
# HPC:
#   conda activate clusterProfiler
#   Rscript OSK_R4vsR3_fgsea_bubble_2026.R
#
#   # only redraw bubble (C6 → Homeostatic_signaling):
#   Rscript OSK_R4vsR3_fgsea_bubble_2026.R --replot_only
#
#   Rscript OSK_R4vsR3_fgsea_bubble_2026.R --run_terms

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}
has_flag <- function(flag) !is.na(match(flag, args))

root <- "/projects/hmz-aging/Xenium5k_Mouse/OSK_New"
rank_dir <- get_arg(
  "--rank_dir",
  file.path(root, "OSK_R4vsR3_Aging_GSEA_2026")
)
out_dir <- get_arg("--out_dir", rank_dir)
catalog_path <- get_arg(
  "--catalog",
  "/projects/hmz-aging/HMZ_GO_terms_Final/human_aging_catalog"
)
rank_rds <- get_arg(
  "--rank_rds",
  file.path(rank_dir, "rank_human_OSKR4vsOSKR3.rds")
)
run_terms <- has_flag("--run_terms")
replot_only <- has_flag("--replot_only")
min_size <- as.integer(get_arg("--minSize", "5"))
max_size_cat <- as.integer(get_arg("--maxSize_cat", "5000"))
max_size_term <- as.integer(get_arg("--maxSize_term", "2000"))
pval_cut <- as.numeric(get_arg("--pval_cut", "0.05"))
# Plot / run order (x-axis): match retina / Proximity figure
celltypes <- strsplit(
  get_arg("--celltypes", "MG,RGC,AC,HC,BC,Cone,Rod,RPE"),
  ",",
  fixed = TRUE
)[[1]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
  flush.console()
}

collapse_leading_edge <- function(le, sep = ";") {
  vapply(le, function(x) {
    if (is.null(x) || length(x) == 0) "" else paste(as.character(x), collapse = sep)
  }, character(1))
}

read_rnk <- function(path) {
  df <- read.table(
    path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
    col.names = c("gene", "score"), quote = "", comment.char = ""
  )
  rk <- as.numeric(df$score)
  names(rk) <- as.character(df$gene)
  rk[is.finite(rk)]
}

run_fgsea_one <- function(ranks, pathways, minSize = 5, maxSize = 5000) {
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  set.seed(123)
  res <- fgsea(
    pathways = pathways, stats = ranks,
    minSize = minSize, maxSize = maxSize, eps = 0
  )
  res[order(res$NES, decreasing = TRUE), ]
}

log_msg("========== OSK_R4vsR3_fgsea_bubble start ==========")
log_msg("rank_dir: ", rank_dir)
log_msg("out_dir : ", out_dir)
log_msg("catalog : ", catalog_path)
log_msg("run_terms: ", run_terms)
log_msg("replot_only: ", replot_only)

# C6 display label (same as aging Xenium GSEA)
cat_relabel <- c(
  C6_Intercellular_communication = "C6_Homeostatic_signaling"
)
relabel_categories <- function(df) {
  hit <- df$category %in% names(cat_relabel)
  if (any(hit)) {
    df$category[hit] <- unname(cat_relabel[df$category[hit]])
    log_msg(
      "Relabeled: ",
      paste(sprintf("%s→%s", names(cat_relabel), cat_relabel), collapse = "; ")
    )
  }
  df
}

make_bubble <- function(GSEA_df) {
  GSEA_df <- relabel_categories(GSEA_df)
  write.table(
    GSEA_df,
    file.path(out_dir, "GSEA_category_OSKR4vsOSKR3.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    GSEA_df[!is.na(GSEA_df$pval) & GSEA_df$pval < pval_cut, , drop = FALSE],
    file.path(out_dir, "GSEA_category_OSKR4vsOSKR3_enriched.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    GSEA_df,
    file.path(out_dir, "OSK_R4vsR3_GSEA_category.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  plot_df <- GSEA_df
  plot_df$neglogP <- -log10(pmax(plot_df$pval, 1e-300))
  cat_all <- unique(as.character(plot_df$category))
  cat_num <- as.numeric(sub("^C([0-9]+)_.*", "\\1", cat_all))
  cat_levels <- cat_all[order(cat_num)]
  ct_levels <- celltypes[celltypes %in% unique(plot_df$celltype)]
  if (length(ct_levels) == 0) ct_levels <- unique(as.character(plot_df$celltype))
  plot_df$category <- factor(plot_df$category, levels = rev(cat_levels))
  plot_df$celltype <- factor(plot_df$celltype, levels = ct_levels)
  plot_sig <- plot_df[plot_df$pval < pval_cut & !is.na(plot_df$pval), ]

  p <- ggplot(plot_sig, aes(x = celltype, y = category)) +
    geom_point(
      aes(size = neglogP, fill = NES),
      shape = 21, color = "grey30", stroke = 0.3
    ) +
    scale_fill_gradient2(
      low = "#3C5488", mid = "white", high = "#E64B35",
      midpoint = 0, name = "NES", limits = c(-3, 3), oob = squish
    ) +
    scale_size_continuous(range = c(4, 14), name = expression(-log[10] ~ italic(P))) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black")
    ) +
    labs(
      x = "", y = "",
      title = sprintf("OSK GSEA C1-C10 | OSKR4 vs OSKR3 (p < %s)", pval_cut)
    )
  ggsave(
    file.path(out_dir, "GSEA_category_OSKR4vsOSKR3_bubble.png"),
    p, width = 7, height = 7, dpi = 300, bg = "white"
  )
  ggsave(
    file.path(out_dir, "OSK_R4vsR3_GSEA_category_bubble.png"),
    p, width = 7, height = 7, dpi = 300, bg = "white"
  )
  ggsave(
    file.path(out_dir, "GSEA_category_OSKR4vsOSKR3_bubble.pdf"),
    p, width = 7, height = 7, bg = "white"
  )
  log_msg("Wrote GSEA_category_OSKR4vsOSKR3_bubble.png / .pdf")
  invisible(GSEA_df)
}

if (isTRUE(replot_only)) {
  tsv_cands <- c(
    file.path(out_dir, "GSEA_category_OSKR4vsOSKR3.tsv"),
    file.path(out_dir, "OSK_R4vsR3_GSEA_category.tsv"),
    file.path(rank_dir, "GSEA_category_OSKR4vsOSKR3.tsv")
  )
  tsv_hit <- tsv_cands[file.exists(tsv_cands)][1]
  if (is.na(tsv_hit)) stop("No GSEA TSV for --replot_only under ", out_dir)
  log_msg("Replot from: ", tsv_hit)
  GSEA_df <- read.delim(tsv_hit, stringsAsFactors = FALSE)
  make_bubble(GSEA_df)
  log_msg("========== replot done ==========")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Load ranks
# ---------------------------------------------------------------------------

rank_cands <- c(
  rank_rds,
  file.path(rank_dir, "rank_human_OSKR4vsOSKR3.rds"),
  file.path(rank_dir, "OSK_R4vsR3_DEG_rank_human.rds"),
  file.path(rank_dir, "rank_human_R4vsR1.rds")
)
rank_hit <- rank_cands[file.exists(rank_cands)][1]

if (!is.na(rank_hit)) {
  log_msg("Load ranks from: ", rank_hit)
  rank_human <- readRDS(rank_hit)
} else {
  rnk_dir <- file.path(rank_dir, "ranks_human")
  rnk_files <- list.files(rnk_dir, pattern = "^rank_human_.*\\.rnk$", full.names = TRUE)
  if (length(rnk_files) == 0) {
    stop(
      "No ranks found. Run DEG first:\n",
      "  conda activate seurat_env\n",
      "  Rscript OSK_R4vsR3_DEG_rank_2026.R\n",
      "Expected: ", rank_rds
    )
  }
  log_msg("Load ranks from ", length(rnk_files), " .rnk files in ", rnk_dir)
  rank_human <- lapply(rnk_files, read_rnk)
  names(rank_human) <- sub("^rank_human_|\\.rnk$", "", basename(rnk_files))
}

if (!is.null(celltypes) && length(celltypes) > 0) {
  keep <- intersect(celltypes, names(rank_human))
  if (length(keep) == 0) stop("No requested celltypes found in ranks")
  rank_human <- rank_human[keep]
}
log_msg("celltypes: ", paste(names(rank_human), collapse = ", "))

# ---------------------------------------------------------------------------
# Catalog
# ---------------------------------------------------------------------------

if (!file.exists(catalog_path)) stop("Catalog not found: ", catalog_path)
load(catalog_path)
if (!exists("human_aging_catalog")) stop("Expected human_aging_catalog")
category_sets <- human_aging_catalog$category_unions
log_msg("category_unions: ", length(category_sets), " sets")

# ---------------------------------------------------------------------------
# C1–C10 GSEA
# ---------------------------------------------------------------------------

GSEA_cat <- list()
for (ct in names(rank_human)) {
  ranks <- rank_human[[ct]]
  if (length(ranks) < 10) {
    log_msg("Skip ", ct, ": too few ranked genes (", length(ranks), ")")
    next
  }
  log_msg("==== fgsea category: ", ct, " (n=", length(ranks), ") ====")
  res <- run_fgsea_one(ranks, category_sets, minSize = min_size, maxSize = max_size_cat)
  res$celltype <- ct
  GSEA_cat[[ct]] <- res
}

if (length(GSEA_cat) == 0) stop("No GSEA category results")

GSEA_df <- do.call(rbind, lapply(GSEA_cat, function(x) {
  data.frame(
    celltype = x$celltype,
    category = x$pathway,
    NES = x$NES,
    ES = x$ES,
    pval = x$pval,
    padj = x$padj,
    size = x$size,
    nLeadingEdge = vapply(x$leadingEdge, length, integer(1)),
    leadingEdge = collapse_leading_edge(x$leadingEdge),
    stringsAsFactors = FALSE
  )
}))
rownames(GSEA_df) <- NULL
GSEA_df <- GSEA_df[order(GSEA_df$celltype, -GSEA_df$NES), ]
saveRDS(GSEA_cat, file.path(out_dir, "GSEA_category_OSKR4vsOSKR3.rds"))
GSEA_df <- make_bubble(GSEA_df)
log_msg("Wrote GSEA_category_OSKR4vsOSKR3.tsv / _enriched.tsv / .rds")

# ---------------------------------------------------------------------------
# Optional term-level GSEA
# ---------------------------------------------------------------------------

if (isTRUE(run_terms)) {
  term_sets <- human_aging_catalog$pathways
  term_anno <- human_aging_catalog$annotation
  log_msg("pathways (terms): ", length(term_sets), " sets")

  GSEA_term <- list()
  for (ct in names(rank_human)) {
    ranks <- rank_human[[ct]]
    if (length(ranks) < 10) next
    log_msg("==== fgsea term: ", ct, " ====")
    res <- run_fgsea_one(ranks, term_sets, minSize = min_size, maxSize = max_size_term)
    res$celltype <- ct
    GSEA_term[[ct]] <- res
  }

  term_df <- do.call(rbind, lapply(GSEA_term, function(x) {
    data.frame(
      celltype = x$celltype,
      term_id = x$pathway,
      NES = x$NES, ES = x$ES, pval = x$pval, padj = x$padj, size = x$size,
      nLeadingEdge = vapply(x$leadingEdge, length, integer(1)),
      leadingEdge = collapse_leading_edge(x$leadingEdge),
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(term_df) && nrow(term_df) > 0) {
    keep_anno <- intersect(
      c("term_id", "source", "gene_set", "category", "subcategory"),
      colnames(term_anno)
    )
    term_df <- merge(
      term_df, term_anno[, keep_anno, drop = FALSE],
      by = "term_id", all.x = TRUE
    )
    term_df <- term_df[order(term_df$celltype, term_df$category, -term_df$NES), ]
  }
  write.table(
    term_df,
    file.path(out_dir, "GSEA_term_OSKR4vsOSKR3.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    term_df[!is.na(term_df$pval) & term_df$pval < pval_cut, , drop = FALSE],
    file.path(out_dir, "GSEA_term_OSKR4vsOSKR3_enriched.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  saveRDS(GSEA_term, file.path(out_dir, "GSEA_term_OSKR4vsOSKR3.rds"))
  log_msg("Wrote GSEA_term_OSKR4vsOSKR3.tsv / _enriched.tsv / .rds")
}

log_msg("========== OSK_R4vsR3_fgsea_bubble done ==========")
