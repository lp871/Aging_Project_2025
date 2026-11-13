#!/usr/bin/env Rscript
########################################################################
# Re-run aging GSEA (fgsea) and SAVE leadingEdge genes (2026)
#
# Parameters identical to the earlier DEGs_revised_and_Rank*:
#   category: minSize=5, maxSize=5000, set.seed(123), eps=0
#   term:     minSize=5, maxSize=2000
#   gene sets: human_aging_catalog$category_unions / $pathways
#
# Species:
#   Mouse / Zebrafish: use ranks already mapped to human symbols
#   Human:             Human_DEG_rank.rds (by celltype_sex; ENSG->symbol)
#
# Outputs (in each species directory):
#   *_GSEA_withLE_2026.xlsx                  # multi-tab Excel (main output)
#   *_GSEA_category_withLE_2026.tsv          # statistics + leadingEdge (semicolon-joined)
#   *_GSEA_category_leadingEdge_long_2026.tsv# long table: one LE gene per row
#   *_GSEA_category_fgsea_2026.rds           # full fgsea list (with LE)
#   *_GSEA_term_withLE_2026.tsv / _leadingEdge_long / _fgsea.rds
#
# Excel tabs:
#   Contents | category_summary | category_LE_long
#            | term_summary | term_LE_long
#            | category_sig_p05 | term_sig_p05
#            | LE_<celltype> (LE per category, one tab per celltype)
#
# HPC:
#   conda activate clusterProfiler
#   Rscript Run_Aging_GSEA_withLeadingEdge_2026.R
#   Rscript Run_Aging_GSEA_withLeadingEdge_2026.R --species Mouse
#   Rscript Run_Aging_GSEA_withLeadingEdge_2026.R --species Mouse,Zebrafish
#   Rscript Run_Aging_GSEA_withLeadingEdge_2026.R --level category
########################################################################

suppressMessages({
    library(fgsea)
    library(dplyr)
})

if (!requireNamespace("openxlsx", quietly = TRUE)) {
    install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)

## args ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) == 0 || i[1] >= length(args)) return(default)
    args[i[1] + 1]
}
species_arg <- get_arg("--species", "Mouse,Zebrafish,Human")
level_arg   <- get_arg("--level", "both")   # category / term / both
species_run <- trimws(unlist(strsplit(species_arg, ",")))
level_run   <- match.arg(level_arg, c("category", "term", "both"))

catalog_path <- "/projects/hmz-aging/HMZ_GO_terms_Final/human_aging_catalog"
human_dir <- "/projects/hmz-aging/Human_New_Datasets"
mouse_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
fish_dir  <- "/projects/hmz-aging/Zebrafish_New_Datasets"
pb_root   <- file.path(human_dir, "pseudobulk_sample")


########################################################################
# helpers
########################################################################
run_fgsea_one <- function(ranks, pathways, minSize = 5, maxSize = 5000) {
    ranks <- ranks[is.finite(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    set.seed(123)
    res <- fgsea(pathways = pathways, stats = ranks,
                 minSize = minSize, maxSize = maxSize, eps = 0)
    res[order(res$NES, decreasing = TRUE), ]
}

## leadingEdge list -> "GENE1;GENE2;..."
le_to_str <- function(le) {
    if (is.null(le) || length(le) == 0) return(NA_character_)
    # fgsea may return list-column; each element is character vector
    if (is.list(le)) {
        return(vapply(le, function(x) {
            if (is.null(x) || length(x) == 0) return(NA_character_)
            paste(as.character(x), collapse = ";")
        }, character(1)))
    }
    paste(as.character(le), collapse = ";")
}

flatten_fgsea_with_le <- function(gsea_list, pathway_col = "pathway") {
    rows <- lapply(names(gsea_list), function(ct) {
        x <- gsea_list[[ct]]
        if (is.null(x) || nrow(x) == 0) return(NULL)
        data.frame(
            celltype = ct,
            pathway = as.character(x$pathway),
            NES = x$NES,
            ES = x$ES,
            pval = x$pval,
            padj = x$padj,
            size = x$size,
            leadingEdge = le_to_str(x$leadingEdge),
            n_leadingEdge = vapply(x$leadingEdge, length, integer(1)),
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)
    if (is.null(df) || nrow(df) == 0) return(df)
    rownames(df) <- NULL
    df
}

## long table: one gene per row
fgsea_list_to_le_long <- function(gsea_list) {
    rows <- list()
    for (ct in names(gsea_list)) {
        x <- gsea_list[[ct]]
        if (is.null(x) || nrow(x) == 0) next
        for (i in seq_len(nrow(x))) {
            genes <- x$leadingEdge[[i]]
            if (is.null(genes) || length(genes) == 0) next
            rows[[length(rows) + 1]] <- data.frame(
                celltype = ct,
                pathway = as.character(x$pathway[i]),
                NES = x$NES[i],
                ES = x$ES[i],
                pval = x$pval[i],
                padj = x$padj[i],
                size = x$size[i],
                gene = as.character(genes),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(rows) == 0) {
        return(data.frame(celltype = character(0), pathway = character(0),
                          NES = numeric(0), ES = numeric(0), pval = numeric(0),
                          padj = numeric(0), size = integer(0), gene = character(0),
                          stringsAsFactors = FALSE))
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

run_gsea_on_ranks <- function(rank_list, pathways, minSize, maxSize, label) {
    gsea_list <- list()
    for (ct in names(rank_list)) {
        ranks <- rank_list[[ct]]
        if (is.null(ranks) || length(ranks) < 10) {
            message("  skip (too few genes): ", ct)
            next
        }
        message(sprintf("  ==== %s fgsea: %s (n_rank=%d) ====", label, ct, length(ranks)))
        res <- run_fgsea_one(ranks, pathways, minSize = minSize, maxSize = maxSize)
        res$celltype <- ct
        gsea_list[[ct]] <- res
    }
    gsea_list
}

write_gsea_outputs <- function(gsea_list, out_prefix, pathway_name = "category",
                               anno = NULL) {
    ## full rds
    rds_fp <- paste0(out_prefix, "_fgsea_2026.rds")
    saveRDS(gsea_list, rds_fp)
    message("Saved: ", rds_fp)

    ## wide/summary with LE string
    df <- flatten_fgsea_with_le(gsea_list)
    if (is.null(df) || nrow(df) == 0) {
        message("  empty GSEA result, skip tsv")
        return(invisible(NULL))
    }
    names(df)[names(df) == "pathway"] <- pathway_name

    if (!is.null(anno) && pathway_name == "term_id") {
        keep <- intersect(c("term_id", "source", "gene_set", "category", "subcategory"),
                          colnames(anno))
        df <- merge(df, anno[, keep, drop = FALSE], by = "term_id", all.x = TRUE)
    }

    ## sort
    if ("celltype" %in% names(df) && "NES" %in% names(df)) {
        if ("category" %in% names(df) && pathway_name == "term_id") {
            df <- df[order(df$celltype, df$category, -df$NES), ]
        } else {
            df <- df[order(df$celltype, -df$NES), ]
        }
    }
    rownames(df) <- NULL

    tsv_fp <- paste0(out_prefix, "_withLE_2026.tsv")
    write.table(df, tsv_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved: ", tsv_fp, "  rows=", nrow(df))

    ## long LE
    le_long <- fgsea_list_to_le_long(gsea_list)
    names(le_long)[names(le_long) == "pathway"] <- pathway_name
    if (!is.null(anno) && pathway_name == "term_id" && nrow(le_long) > 0) {
        keep <- intersect(c("term_id", "source", "gene_set", "category", "subcategory"),
                          colnames(anno))
        le_long <- merge(le_long, anno[, keep, drop = FALSE], by = "term_id", all.x = TRUE)
    }
    le_fp <- paste0(out_prefix, "_leadingEdge_long_2026.tsv")
    write.table(le_long, le_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Saved: ", le_fp, "  rows=", nrow(le_long),
            "  unique_genes=", length(unique(le_long$gene)))
    invisible(list(df = df, le_long = le_long, gsea_list = gsea_list))
}

## ---- Excel: multi-tab workbook with Contents ----
safe_sheet_name <- function(x) {
    x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
    substr(x, 1, 31)
}

write_species_gsea_excel <- function(species, out_dir, category_out = NULL, term_out = NULL) {
    wb <- createWorkbook()
    hs <- createStyle(textDecoration = "bold", halign = "center",
                      fgFill = "#D9E1F2", border = "TopBottomLeftRight")
    title_style <- createStyle(textDecoration = "bold", fontSize = 14)

    sheet_specs <- list()
    push_sheet <- function(tab, description, df) {
        if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
            message("  skip empty sheet: ", tab)
            return(invisible(NULL))
        }
        tab31 <- safe_sheet_name(tab)
        ## avoid duplicate sheet names
        existing <- vapply(sheet_specs, `[[`, "", "Tab")
        if (tab31 %in% existing) {
            k <- 1
            while (paste0(substr(tab31, 1, 28), "_", k) %in% existing) k <- k + 1
            tab31 <- paste0(substr(tab31, 1, 28), "_", k)
        }
        sheet_specs[[length(sheet_specs) + 1]] <<- list(
            Tab = tab31,
            Description = description,
            N_rows = nrow(df),
            N_cols = ncol(df),
            df = df
        )
    }

    if (!is.null(category_out)) {
        push_sheet("category_summary",
                   "Category-level GSEA (C1-C10); leadingEdge as semicolon-separated genes",
                   category_out$df)
        push_sheet("category_LE_long",
                   "Category leadingEdge long table (one gene per row)",
                   category_out$le_long)
        if (!is.null(category_out$df) && nrow(category_out$df) > 0) {
            sig <- category_out$df[!is.na(category_out$df$pval) & category_out$df$pval < 0.05, ,
                                   drop = FALSE]
            push_sheet("category_sig_p05",
                       "Category GSEA with pval < 0.05",
                       sig)
        }
        ## per-celltype LE tabs (category)
        if (!is.null(category_out$le_long) && nrow(category_out$le_long) > 0 &&
            "celltype" %in% names(category_out$le_long)) {
            for (ct in unique(as.character(category_out$le_long$celltype))) {
                sub <- category_out$le_long[category_out$le_long$celltype == ct, , drop = FALSE]
                push_sheet(paste0("LE_", ct),
                           sprintf("Category leadingEdge genes for %s", ct),
                           sub)
            }
        }
    }

    if (!is.null(term_out)) {
        push_sheet("term_summary",
                   "Term-level GSEA; leadingEdge as semicolon-separated genes",
                   term_out$df)
        push_sheet("term_LE_long",
                   "Term leadingEdge long table (one gene per row)",
                   term_out$le_long)
        if (!is.null(term_out$df) && nrow(term_out$df) > 0) {
            sig <- term_out$df[!is.na(term_out$df$pval) & term_out$df$pval < 0.05, ,
                               drop = FALSE]
            push_sheet("term_sig_p05",
                       "Term GSEA with pval < 0.05",
                       sig)
        }
    }

    if (length(sheet_specs) == 0) {
        message("[", species, "] no sheets to write, skip Excel")
        return(invisible(NULL))
    }

    addWorksheet(wb, "Contents")
    contents <- data.frame(
        Tab = vapply(sheet_specs, `[[`, "", "Tab"),
        N_rows = vapply(sheet_specs, `[[`, 0L, "N_rows"),
        N_cols = vapply(sheet_specs, `[[`, 0L, "N_cols"),
        Description = vapply(sheet_specs, `[[`, "", "Description"),
        stringsAsFactors = FALSE
    )

    writeData(wb, "Contents",
              sprintf("%s Aging GSEA with leadingEdge (2026)", species),
              startRow = 1, startCol = 1)
    addStyle(wb, "Contents", title_style, rows = 1, cols = 1)
    writeData(wb, "Contents",
              sprintf("Catalog: human_aging_catalog | Generated: %s",
                      format(Sys.time(), "%Y-%m-%d %H:%M")),
              startRow = 2, startCol = 1)
    writeData(wb, "Contents", contents, startRow = 4, startCol = 1, headerStyle = hs)

    for (i in seq_along(sheet_specs)) {
        spec <- sheet_specs[[i]]
        addWorksheet(wb, spec$Tab)
        writeData(wb, spec$Tab, spec$df, headerStyle = hs)
        freezePane(wb, spec$Tab, firstRow = TRUE)
        setColWidths(wb, spec$Tab, cols = seq_len(min(ncol(spec$df), 14)), widths = "auto")

        link <- sprintf("internal:'%s'!A1", spec$Tab)
        class(link) <- "hyperlink"
        writeData(wb, "Contents", x = link, startRow = 4 + i, startCol = 1)
    }
    freezePane(wb, "Contents", firstActiveRow = 5)
    setColWidths(wb, "Contents", cols = 1:4, widths = c(28, 10, 10, 70))

    excel_out <- file.path(out_dir, sprintf("%s_GSEA_withLE_2026.xlsx", species))
    saveWorkbook(wb, excel_out, overwrite = TRUE)
    message(sprintf("Saved Excel: %s  (%d data tabs + Contents)",
                    excel_out, length(sheet_specs)))
    invisible(excel_out)
}


########################################################################
# load catalog
########################################################################
if (!file.exists(catalog_path)) stop("Cannot find catalog: ", catalog_path)
load(catalog_path)   # -> human_aging_catalog
stopifnot(exists("human_aging_catalog"))
category_sets <- human_aging_catalog$category_unions
term_sets     <- human_aging_catalog$pathways
term_anno     <- human_aging_catalog$annotation
message("[catalog] categories=", length(category_sets),
        " | terms=", length(term_sets))


########################################################################
# load ranks per species
########################################################################
load_mouse_ranks <- function() {
    cand <- c(
        file.path(mouse_dir, "Mouse_DEG_rank_human.rds"),
        file.path(mouse_dir, "Mouse_DEG_rank_human")
    )
    hit <- cand[file.exists(cand)]
    if (length(hit) == 0) stop("Cannot find Mouse_DEG_rank_human(.rds)")
    fp <- hit[[1]]
    message("[Mouse] load ranks: ", fp)
    if (grepl("\\.rds$", fp, ignore.case = TRUE)) {
        obj <- readRDS(fp)
    } else {
        env <- new.env(parent = emptyenv())
        load(fp, envir = env)
        if (exists("Mouse_DEG_rank_human", envir = env)) {
            obj <- get("Mouse_DEG_rank_human", envir = env)
        } else {
            objs <- ls(envir = env)
            if (length(objs) != 1) stop("RData contains more than one object: ", paste(objs, collapse = ","))
            obj <- get(objs[[1]], envir = env)
        }
    }
    obj
}

load_fish_ranks <- function() {
    fp <- file.path(fish_dir, "Zebrafish_DEG_rank_human.rds")
    if (!file.exists(fp)) {
        alt <- file.path(fish_dir, "Zebrafish_DEG_rank_human")
        if (file.exists(alt)) fp <- alt else stop("Cannot find Zebrafish_DEG_rank_human.rds")
    }
    message("[Zebrafish] load ranks: ", fp)
    if (grepl("\\.rds$", fp, ignore.case = TRUE)) {
        readRDS(fp)
    } else {
        env <- new.env(parent = emptyenv())
        load(fp, envir = env)
        if (exists("Zebrafish_DEG_rank_human", envir = env)) {
            get("Zebrafish_DEG_rank_human", envir = env)
        } else {
            get(ls(envir = env)[[1]], envir = env)
        }
    }
}

load_human_ranks <- function() {
    fp <- file.path(pb_root, "Human_DEG_rank.rds")
    if (!file.exists(fp)) stop("Cannot find Human_DEG_rank.rds: ", fp)
    message("[Human] load ranks: ", fp)
    x <- readRDS(fp)

    ## flatten nested list -> AC_male
    flatten_rank_list <- function(x) {
        if (length(x) == 0) return(x)
        if (is.numeric(x[[1]])) return(x)
        flat <- list()
        for (ct in names(x))
            for (sx in names(x[[ct]]))
                flat[[paste(ct, sx, sep = "_")]] <- x[[ct]][[sx]]
        flat
    }
    x <- flatten_rank_list(x)

    ## ENSG -> SYMBOL (same as the old script)
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        suppressMessages(library(org.Hs.eg.db))
        convert_rank_to_symbol <- function(rank_vec) {
            if (length(rank_vec) == 0) return(rank_vec)
            nm <- names(rank_vec)
            is_ens <- grepl("^ENSG", nm)
            if (!any(is_ens)) return(rank_vec)
            sym <- nm
            ens <- sub("\\..*$", "", nm[is_ens])
            sym[is_ens] <- AnnotationDbi::mapIds(
                org.Hs.eg.db, keys = ens, keytype = "ENSEMBL",
                column = "SYMBOL", multiVals = "first"
            )
            df <- data.frame(sym = sym, val = as.numeric(rank_vec),
                             stringsAsFactors = FALSE)
            df <- df[!is.na(df$sym), ]
            out <- tapply(df$val, df$sym, function(v) v[which.max(abs(v))])
            sort(out, decreasing = TRUE)
        }
        x <- lapply(x, convert_rank_to_symbol)
    } else {
        message("[Human] org.Hs.eg.db unavailable; can continue if ranks are already symbols")
    }
    x
}


########################################################################
# run
########################################################################
message("species = ", paste(species_run, collapse = ", "),
        " | level = ", level_run)

do_category <- level_run %in% c("category", "both")
do_term     <- level_run %in% c("term", "both")

## ---- Mouse ----
if ("Mouse" %in% species_run) {
    message("\n========== Mouse ==========")
    setwd(mouse_dir)
    ranks <- load_mouse_ranks()
    message("  celltypes: ", paste(names(ranks), collapse = ", "))

    mouse_cat <- NULL
    mouse_term <- NULL
    if (do_category) {
        gcat <- run_gsea_on_ranks(ranks, category_sets, 5, 5000, "Mouse-category")
        mouse_cat <- write_gsea_outputs(gcat, file.path(mouse_dir, "Mouse_GSEA_category"),
                                        pathway_name = "category")
    }
    if (do_term) {
        gterm <- run_gsea_on_ranks(ranks, term_sets, 5, 2000, "Mouse-term")
        mouse_term <- write_gsea_outputs(gterm, file.path(mouse_dir, "Mouse_GSEA_term"),
                                         pathway_name = "term_id", anno = term_anno)
    }
    write_species_gsea_excel("Mouse", mouse_dir, mouse_cat, mouse_term)
}

## ---- Zebrafish ----
if ("Zebrafish" %in% species_run) {
    message("\n========== Zebrafish ==========")
    setwd(fish_dir)
    ranks <- load_fish_ranks()
    message("  celltypes: ", paste(names(ranks), collapse = ", "))

    fish_cat <- NULL
    fish_term <- NULL
    if (do_category) {
        gcat <- run_gsea_on_ranks(ranks, category_sets, 5, 5000, "Zebrafish-category")
        fish_cat <- write_gsea_outputs(gcat, file.path(fish_dir, "Zebrafish_GSEA_category"),
                                       pathway_name = "category")
    }
    if (do_term) {
        gterm <- run_gsea_on_ranks(ranks, term_sets, 5, 2000, "Zebrafish-term")
        fish_term <- write_gsea_outputs(gterm, file.path(fish_dir, "Zebrafish_GSEA_term"),
                                        pathway_name = "term_id", anno = term_anno)
    }
    write_species_gsea_excel("Zebrafish", fish_dir, fish_cat, fish_term)
}

## ---- Human ----
if ("Human" %in% species_run) {
    message("\n========== Human ==========")
    setwd(human_dir)
    ranks <- load_human_ranks()
    message("  groups: ", paste(names(ranks), collapse = ", "))

    human_cat <- NULL
    human_term <- NULL
    if (do_category) {
        gcat <- run_gsea_on_ranks(ranks, category_sets, 5, 5000, "Human-category")
        human_cat <- write_gsea_outputs(gcat, file.path(human_dir, "Human_GSEA_category"),
                                        pathway_name = "category")
        ## Add sex / ct / mode (compatible with old tables)
        if (!is.null(human_cat$df)) {
            df <- human_cat$df
            df$sex  <- sub(".*_(male|female)$", "\\1", df$celltype)
            df$ct   <- sub("_(male|female)$", "", df$celltype)
            df$mode <- "noKmeans"
            human_cat$df <- df
            write.table(df, file.path(human_dir, "Human_GSEA_category_withLE_2026.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            message("Updated Human category tsv with sex/ct/mode")
        }
        if (!is.null(human_cat$le_long) && nrow(human_cat$le_long) > 0) {
            le <- human_cat$le_long
            le$sex  <- sub(".*_(male|female)$", "\\1", le$celltype)
            le$ct   <- sub("_(male|female)$", "", le$celltype)
            le$mode <- "noKmeans"
            human_cat$le_long <- le
            write.table(le, file.path(human_dir, "Human_GSEA_category_leadingEdge_long_2026.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
        }
    }
    if (do_term) {
        gterm <- run_gsea_on_ranks(ranks, term_sets, 5, 2000, "Human-term")
        human_term <- write_gsea_outputs(gterm, file.path(human_dir, "Human_GSEA_term"),
                                         pathway_name = "term_id", anno = term_anno)
        if (!is.null(human_term$df)) {
            df <- human_term$df
            df$sex  <- sub(".*_(male|female)$", "\\1", df$celltype)
            df$ct   <- sub("_(male|female)$", "", df$celltype)
            df$mode <- "noKmeans"
            human_term$df <- df
            write.table(df, file.path(human_dir, "Human_GSEA_term_withLE_2026.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
        }
        if (!is.null(human_term$le_long) && nrow(human_term$le_long) > 0) {
            le <- human_term$le_long
            le$sex  <- sub(".*_(male|female)$", "\\1", le$celltype)
            le$ct   <- sub("_(male|female)$", "", le$celltype)
            le$mode <- "noKmeans"
            human_term$le_long <- le
            write.table(le, file.path(human_dir, "Human_GSEA_term_leadingEdge_long_2026.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
        }
    }
    write_species_gsea_excel("Human", human_dir, human_cat, human_term)
}

message("\n=== Done: GSEA with leadingEdge ===")
message("Main output Excel:")
message("  /projects/hmz-aging/Mouse_New_Datasets/Mouse_GSEA_withLE_2026.xlsx")
message("  /projects/hmz-aging/Zebrafish_New_Datasets/Zebrafish_GSEA_withLE_2026.xlsx")
message("  /projects/hmz-aging/Human_New_Datasets/Human_GSEA_withLE_2026.xlsx")
message("tsv/rds are still written; Excel contains Contents + multiple tabs")
