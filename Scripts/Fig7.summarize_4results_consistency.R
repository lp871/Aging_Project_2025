#!/usr/bin/env Rscript
# Count directional consistency across 4 results, per celltype, Up/Down separate.
#
# Sources / calls (pvalue < 0.05 for edgeR):
#   Xenium_clock : weight > 0 => Up; weight < 0 => Down
#   snRNA_clock  : weight > 0 => Up; weight < 0 => Down
#   Xenium_edgeR : PValue<0.05 & logFC>0 => Up; PValue<0.05 & logFC<0 => Down
#   snRNA_edgeR  : PValue<0.05 & logFC>0 => Up; PValue<0.05 & logFC<0 => Down
#
# Consistency:
#   all4 : all 4 sources call the same direction
#   exactly3 : exactly 3 call the same direction; the 4th MUST be missing (NA),
#              not the opposite direction
#
# Example:
#   Rscript summarize_4results_consistency.R
#   Rscript summarize_4results_consistency.R --merged /path/Merged_4results_clocks_edgeR.xlsx

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    install.packages("openxlsx", repos = "https://cloud.r-project.org")
  }
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}

root <- "/projects/hmz-aging/ProximityEffect"
merged_xlsx <- get_arg(
  "--merged",
  file.path(root, "Merged_4results_clocks_edgeR.xlsx")
)
out_dir <- get_arg(
  "--out_dir",
  file.path(root, "Merged_4results_consistency_p0.05")
)
pval_cut <- as.numeric(get_arg("--pval_cut", "0.05"))

# default: 5 major CTs requested for saving exactly3 lists
celltypes <- strsplit(
  get_arg("--celltypes", "MG,AC,BC,Cone,Rod"),
  ",",
  fixed = TRUE
)[[1]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
  flush.console()
}

log_msg("========== summarize_4results_consistency start ==========")
log_msg("merged   : ", merged_xlsx)
log_msg("out_dir  : ", out_dir)
log_msg("pval_cut : ", pval_cut)

if (!file.exists(merged_xlsx)) {
  stop(
    "Merged Excel not found: ", merged_xlsx,
    "\nRun merge_4results_clocks_edgeR_excel.R first."
  )
}

# +1 Up, -1 Down, NA no-call
call_from_weight <- function(w) {
  ifelse(is.na(w) | w == 0, NA_integer_, ifelse(w > 0, 1L, -1L))
}
call_from_edger <- function(logFC, pval, pcut) {
  ok <- !is.na(logFC) & !is.na(pval) & pval < pcut
  ifelse(!ok, NA_integer_, ifelse(logFC > 0, 1L, -1L))
}

# n_agree=4: all 4 == target
# n_agree=3: exactly 3 == target, and the remaining 1 is NA (not opposite)
genes_agree <- function(df, n_agree, direction = c("up", "down")) {
  direction <- match.arg(direction)
  target <- if (direction == "up") 1L else -1L
  opposite <- -target
  mat <- cbind(
    df$call_xe_w,
    df$call_sn_w,
    df$call_xe_e,
    df$call_sn_e
  )
  n_match <- rowSums(mat == target, na.rm = TRUE)
  n_opp <- rowSums(mat == opposite, na.rm = TRUE)
  n_na <- rowSums(is.na(mat))
  if (n_agree == 4L) {
    keep <- n_match == 4L & n_opp == 0L & n_na == 0L
  } else if (n_agree == 3L) {
    # 3 agree, 1 missing only (no opposite allowed)
    keep <- n_match == 3L & n_opp == 0L & n_na == 1L
  } else {
    stop("n_agree must be 3 or 4")
  }
  df$Gene[keep]
}

sheets <- getSheetNames(merged_xlsx)
sheets <- setdiff(sheets, "Contents")
cts_use <- intersect(celltypes, sheets)
if (length(cts_use) == 0) {
  # try all non-Contents sheets
  cts_use <- sheets
}
log_msg("celltypes in Excel: ", paste(cts_use, collapse = ", "))

summary_rows <- list()
gene_lists <- list()

cat("\n")
cat("============================================================\n")
cat(" Directional consistency (edgeR P < ", pval_cut, ")\n", sep = "")
cat(" Up  = clock weight>0 & edgeR logFC>0\n")
cat(" Down= clock weight<0 & edgeR logFC<0\n")
cat(" exactly3: 3 agree + 1 missing (NA); opposite not allowed\n")
cat("============================================================\n\n")

for (ct in cts_use) {
  df <- read.xlsx(merged_xlsx, sheet = ct)
  need <- c(
    "Gene",
    "Xenium_clock_weight", "snRNA_clock_weight",
    "Xenium_logFC", "Xenium_PValue",
    "snRNA_logFC", "snRNA_PValue"
  )
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) {
    log_msg("SKIP ", ct, ": missing cols ", paste(miss, collapse = ", "))
    next
  }

  df$call_xe_w <- call_from_weight(df$Xenium_clock_weight)
  df$call_sn_w <- call_from_weight(df$snRNA_clock_weight)
  df$call_xe_e <- call_from_edger(df$Xenium_logFC, df$Xenium_PValue, pval_cut)
  df$call_sn_e <- call_from_edger(df$snRNA_logFC, df$snRNA_PValue, pval_cut)

  n_call <- rowSums(!is.na(cbind(
    df$call_xe_w, df$call_sn_w, df$call_xe_e, df$call_sn_e
  )))

  up4 <- genes_agree(df, 4L, "up")
  dn4 <- genes_agree(df, 4L, "down")
  up3 <- genes_agree(df, 3L, "up")
  dn3 <- genes_agree(df, 3L, "down")

  cat("------------------------------------------------------------\n")
  cat(sprintf("[%s]\n", ct))
  cat(sprintf("  all 4 consistent  | Up: %d | Down: %d | total: %d\n",
              length(up4), length(dn4), length(up4) + length(dn4)))
  cat(sprintf("  exactly 3 agree   | Up: %d | Down: %d | total: %d\n",
              length(up3), length(dn3), length(up3) + length(dn3)))
  if (length(up4) > 0) {
    cat("  all4 Up genes: ", paste(up4, collapse = ", "), "\n", sep = "")
  }
  if (length(dn4) > 0) {
    cat("  all4 Down genes: ", paste(dn4, collapse = ", "), "\n", sep = "")
  }
  if (length(up3) > 0) {
    cat("  exactly3 Up genes: ", paste(up3, collapse = ", "), "\n", sep = "")
  }
  if (length(dn3) > 0) {
    cat("  exactly3 Down genes: ", paste(dn3, collapse = ", "), "\n", sep = "")
  }
  cat("\n")

  summary_rows[[ct]] <- data.frame(
    celltype = ct,
    n_all4_Up = length(up4),
    n_all4_Down = length(dn4),
    n_all4_total = length(up4) + length(dn4),
    n_exactly3_Up = length(up3),
    n_exactly3_Down = length(dn3),
    n_exactly3_total = length(up3) + length(dn3),
    n_genes_with_any_call = sum(n_call > 0),
    n_genes_with_4_calls = sum(n_call == 4),
    pval_cut = pval_cut,
    stringsAsFactors = FALSE
  )

  # long gene table for this CT
  add_block <- function(genes, level, direction) {
    if (length(genes) == 0) {
      return(data.frame(
        celltype = character(0), Gene = character(0),
        agree_level = character(0), direction = character(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      celltype = ct,
      Gene = genes,
      agree_level = level,
      direction = direction,
      stringsAsFactors = FALSE
    )
  }
  gene_lists[[ct]] <- rbind(
    add_block(up4, "all4", "Up"),
    add_block(dn4, "all4", "Down"),
    add_block(up3, "exactly3", "Up"),
    add_block(dn3, "exactly3", "Down")
  )
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

# overall totals
tot <- data.frame(
  celltype = "ALL_CT_SUM",
  n_all4_Up = sum(summary_df$n_all4_Up),
  n_all4_Down = sum(summary_df$n_all4_Down),
  n_all4_total = sum(summary_df$n_all4_total),
  n_exactly3_Up = sum(summary_df$n_exactly3_Up),
  n_exactly3_Down = sum(summary_df$n_exactly3_Down),
  n_exactly3_total = sum(summary_df$n_exactly3_total),
  n_genes_with_any_call = NA_integer_,
  n_genes_with_4_calls = NA_integer_,
  pval_cut = pval_cut,
  stringsAsFactors = FALSE
)

cat("============================================================\n")
cat(" SUMMARY TABLE\n")
cat("============================================================\n")
print(summary_df, row.names = FALSE)
cat("\n")
print(tot, row.names = FALSE)
cat("\n")

# save
sum_path <- file.path(out_dir, "consistency_summary_p0.05.csv")
write.csv(rbind(summary_df, tot), sum_path, row.names = FALSE)
log_msg("Wrote ", sum_path)

genes_long <- do.call(rbind, gene_lists)
genes_path <- file.path(out_dir, "consistency_genes_long_p0.05.tsv")
write.table(genes_long, genes_path, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote ", genes_path)

# ---- save exactly3 Up / Down only (requested CTs) ----
ex3 <- genes_long[genes_long$agree_level == "exactly3", , drop = FALSE]
ex3_up <- ex3[ex3$direction == "Up", , drop = FALSE]
ex3_dn <- ex3[ex3$direction == "Down", , drop = FALSE]

up_path <- file.path(out_dir, "exactly3_Up_genes_MG_AC_BC_Cone_Rod.tsv")
dn_path <- file.path(out_dir, "exactly3_Down_genes_MG_AC_BC_Cone_Rod.tsv")
write.table(ex3_up, up_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ex3_dn, dn_path, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote ", up_path, " n=", nrow(ex3_up))
log_msg("Wrote ", dn_path, " n=", nrow(ex3_dn))

# per-CT gene lists (one gene per line)
list_dir <- file.path(out_dir, "exactly3_gene_lists")
dir.create(list_dir, recursive = TRUE, showWarnings = FALSE)
for (ct in cts_use) {
  gu <- ex3_up$Gene[ex3_up$celltype == ct]
  gd <- ex3_dn$Gene[ex3_dn$celltype == ct]
  writeLines(gu, file.path(list_dir, paste0(ct, "_exactly3_Up.txt")))
  writeLines(gd, file.path(list_dir, paste0(ct, "_exactly3_Down.txt")))
  log_msg(sprintf("  %s exactly3 | Up=%d Down=%d", ct, length(gu), length(gd)))
}

# also Excel: summary + exactly3 Up/Down + per CT
wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb, "Summary", rbind(summary_df, tot))
addWorksheet(wb, "exactly3_Up")
writeData(wb, "exactly3_Up", ex3_up)
addWorksheet(wb, "exactly3_Down")
writeData(wb, "exactly3_Down", ex3_dn)

for (ct in names(gene_lists)) {
  addWorksheet(wb, substr(ct, 1, 31))
  writeData(wb, substr(ct, 1, 31), gene_lists[[ct]])
}
xlsx_path <- file.path(out_dir, "consistency_p0.05.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
log_msg("Wrote ", xlsx_path)

log_msg("========== summarize_4results_consistency done ==========")
