#!/usr/bin/env Rscript
########################################################################
# 2026 Final clocks → TF enrichment → 种内 / 跨物种 overlap → 图 + Excel
#
# Old / Young 来自最新 Final aging clocks（不是 edgeR / kmeans）:
#   coef > 0  → Old
#   coef < 0  → Young
# Human 用合并版 FinalClock_clean5_2026，不用 bySex。
#
# GRN: 旧版 2025 lists（*_GRNs_list_clcl_2025），unique TF–Gene pairs
#   调控方向 Class = pos/neg，由 TF–Gene Pearson (GG_Cor / GG_Corr) 决定
#   标签 combined_Tag = pos_Old / pos_Young / neg_Old / neg_Young
#
# TF 富集: TF–Gene pair 水平超几何
# 选定 TF 点图暂不画。
#
# 另加 C1–C10 种内 overlap TF bubble（3 张）:
#   Human C2_down：各 CT 上 (LE∩edgeR)_male ∪ (LE∩edgeR)_female
#   Zebrafish C3 down
#   Mouse C5 up
#   三张都按全部 CT 平均 -log10(P) 取前 25
#
# 输出目录:
#   /projects/hmz-aging/Clock_TF_enrich_2026Aug/
#
# HPC:
#   conda activate clusterProfiler
#   Rscript Run_Clock_TF_enrich_overlap_2026Aug.R
########################################################################

suppressMessages({
    library(ggplot2)
})
if (!requireNamespace("openxlsx", quietly = TRUE)) {
    install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)

## ------------------------------------------------------------------
## 路径 / 参数
## ------------------------------------------------------------------
grn_dir         <- "/projects/hmz-aging/GRNs"
mouse_dir       <- "/projects/hmz-aging/Mouse_New_Datasets"
fish_dir        <- "/projects/hmz-aging/Zebrafish_New_Datasets"
human_final_dir <- "/projects/hmz-aging/Human_New_Datasets/FinalClock_clean5_2026"
out_dir         <- "/projects/hmz-aging/Clock_TF_enrich_2026Aug"
if (dir.exists(out_dir)) {
    old <- list.files(out_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    if (length(old)) {
        unlink(old, recursive = TRUE, force = TRUE)
        message("[clean] removed ", length(old), " items from ", out_dir)
    }
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

major_cts     <- c("MG", "RGC", "AC", "HC", "Rod", "Cone", "BC")
species_order <- c("Zebrafish", "Mouse", "Human")
p_cut         <- 0.05          # 显著 TF / overlap
min_shared_ct <- 3L            # 种内直方图：至少共享几个 CT
tags          <- c("pos_Old", "pos_Young", "neg_Old", "neg_Young")
bar_tags      <- c("neg_Young", "neg_Old", "pos_Young", "pos_Old",
                   "neg_NODIFF", "pos_NODIFF")

ortho_candidates <- c(
    file.path(fish_dir, "HMZ_ortholog_combined"),
    "/projects/hmz-aging/Zebrafish_New_Datasets/HMZ_ortholog_combined",
    file.path(mouse_dir, "HMZ_ortholog_combined"),
    "/projects/hmz-aging/HMZ_ortholog_combined",
    "/zp1/data/share/Human_aging_new/HMZ_ortholog_combined_2025"
)

## ------------------------------------------------------------------
## 工具
## ------------------------------------------------------------------
find_existing <- function(cands) {
    hit <- cands[file.exists(cands)]
    if (length(hit) == 0) NA_character_ else hit[[1]]
}

load_ortholog <- function(fp) {
    env <- new.env(parent = emptyenv())
    ok <- tryCatch({ load(fp, envir = env); TRUE }, error = function(e) FALSE)
    if (ok) {
        if (exists("HMZ_ortholog_combined", envir = env, inherits = FALSE))
            return(get("HMZ_ortholog_combined", envir = env))
        objs <- ls(envir = env)
        if (length(objs) == 1L) return(get(objs[[1]], envir = env))
    }
    obj <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (!is.null(obj)) return(obj)
    stop("无法读入 ortholog: ", fp)
}

load_r_object <- function(fp, prefer = NULL) {
    ## 无后缀的 2025 GRN 实际是 saveRDS，不是 load() 的 RData
    obj <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (!is.null(obj)) return(obj)
    env <- new.env(parent = emptyenv())
    load(fp, envir = env)
    if (!is.null(prefer) && exists(prefer, envir = env, inherits = FALSE))
        return(get(prefer, envir = env))
    objs <- ls(envir = env)
    if (length(objs) == 0) stop("空对象: ", fp)
    get(objs[[1]], envir = env)
}

## 最新 Final clocks：Human 合并版；Mouse/Fish 各 CT 目录
clock_rds_path <- function(species, ct) {
    if (species == "Human") {
        file.path(human_final_dir, paste0(ct, "_final_aging_clock_model_2026.rds"))
    } else {
        root <- if (species == "Mouse") mouse_dir else fish_dir
        file.path(root, ct, sprintf("%s_final_aging_clock_model_2026.rds", ct))
    }
}

clock_csv_path <- function(species, ct) {
    if (species == "Human") {
        file.path(human_final_dir, paste0(ct, "_final_model_coefs_2026.csv"))
    } else {
        root <- if (species == "Mouse") mouse_dir else fish_dir
        file.path(root, ct, sprintf("%s_final_model_coefs_2026.csv", ct))
    }
}

## coef>0 = Old, coef<0 = Young；优先 rds$coefs，缺则 csv
read_clock_old_young <- function(species, ct) {
    empty <- list(Old = character(0), Young = character(0), table = NULL, src = NA_character_)
    rds <- clock_rds_path(species, ct)
    csv <- clock_csv_path(species, ct)
    df <- NULL
    src <- NA_character_
    if (file.exists(rds)) {
        obj <- readRDS(rds)
        if (!is.null(obj$coefs)) {
            df <- as.data.frame(obj$coefs)
            src <- rds
        }
    }
    if (is.null(df) && file.exists(csv)) {
        df <- read.csv(csv, stringsAsFactors = FALSE)
        src <- csv
    }
    if (is.null(df) || !all(c("gene", "coef") %in% names(df))) return(empty)
    df$gene <- as.character(df$gene)
    df$coef <- as.numeric(df$coef)
    df <- df[df$gene != "(Intercept)" & !is.na(df$coef) & df$coef != 0, , drop = FALSE]
    if (nrow(df) == 0) return(empty)
    df$species <- species
    df$celltype <- ct
    df$direction <- ifelse(df$coef > 0, "Old", "Young")
    tab <- df[, c("species", "celltype", "gene", "coef", "direction"), drop = FALSE]
    list(Old = tab$gene[tab$direction == "Old"],
         Young = tab$gene[tab$direction == "Young"],
         table = tab, src = src)
}

## 读旧版 2025 GRN list（*_GRNs_list_clcl_2025）；收成 unique TF–Gene pairs
load_grn_list <- function(species) {
    cands <- c(
        file.path(grn_dir, sprintf("%s_GRNs_list_clcl_2025", species)),
        file.path(grn_dir, sprintf("%s_GRNs_list_clcl_2025.rds", species))
    )
    if (species == "Zebrafish") {
        cands <- c(
            cands,
            file.path(grn_dir, "Fish_GRNs_list_clcl_2025"),
            file.path(grn_dir, "Fish_GRNs_list_clcl_2025.rds")
        )
    }
    fp <- find_existing(cands)
    if (is.na(fp)) stop("找不到 ", species, " 旧 GRN: ", paste(cands, collapse = " | "))
    message("[GRN 2025] ", fp)
    normalize_grn_list(load_r_object(fp), species)
}

normalize_grn_list <- function(grn_list, species) {
    out <- list()
    for (ct in names(grn_list)) {
        x <- as.data.frame(grn_list[[ct]])
        if (!"Gene" %in% names(x) && "Target" %in% names(x)) x$Gene <- x$Target
        if (!all(c("TF", "Gene") %in% names(x))) {
            stop(species, "/", ct, " GRN 缺 TF/Gene: ",
                 paste(names(x), collapse = ","))
        }
        x$TF <- as.character(x$TF)
        x$Gene <- as.character(x$Gene)
        x <- x[!is.na(x$TF) & nzchar(x$TF) & !is.na(x$Gene) & nzchar(x$Gene), ]
        n0 <- nrow(x)
        ## 2025 clcl list 已有 Class=pos/neg；没有相关列也没关系
        if (!"Class" %in% names(x) || !all(x$Class %in% c("pos", "neg"))) {
            if (!"GG_Cor" %in% names(x) && "GG_Corr" %in% names(x)) x$GG_Cor <- x$GG_Corr
            if (!"GG_Cor" %in% names(x) && "Pearson" %in% names(x)) x$GG_Cor <- x$Pearson
            if (!"GG_Cor" %in% names(x)) {
                stop(species, "/", ct, " GRN 缺 Class/GG_Cor: ",
                     paste(names(x), collapse = ","))
            }
            x$GG_Cor <- as.numeric(as.character(x$GG_Cor))
            x$Class <- ifelse(!is.na(x$GG_Cor) & x$GG_Cor < 0, "neg", "pos")
        }
        x$Class <- as.character(x$Class)
        x <- x[!duplicated(paste(x$TF, x$Gene, sep = "\t")),
               c("TF", "Gene", "Class"), drop = FALSE]
        out[[ct]] <- x
        message(sprintf("  %s/%s: TF-Gene pairs=%d (from %d rows) | TFs=%d | Genes=%d | pos=%d neg=%d",
                        species, ct, nrow(x), n0,
                        length(unique(x$TF)), length(unique(x$Gene)),
                        sum(x$Class == "pos"), sum(x$Class == "neg")))
    }
    out
}

## 用 clock Old/Young 给每个 TF–Gene pair 打 AGE，再拼 combined_Tag
annotate_grn_with_clock <- function(grn_list, old_list, young_list) {
    out <- list()
    for (ct in names(grn_list)) {
        x <- grn_list[[ct]]
        old_g <- if (ct %in% names(old_list)) old_list[[ct]] else character(0)
        young_g <- if (ct %in% names(young_list)) young_list[[ct]] else character(0)
        x$AGE <- "NODIFF"
        x$AGE[x$Gene %in% old_g] <- "Old"
        x$AGE[x$Gene %in% young_g] <- "Young"
        x$combined_Tag <- paste(x$Class, x$AGE, sep = "_")
        out[[ct]] <- x
        tab <- table(x$combined_Tag)
        message(sprintf("  annotate %s: %s", ct,
                        paste(paste0(names(tab), "=", as.integer(tab)), collapse = " | ")))
    }
    out
}

## TF–Gene pair 水平超几何：每个 TF × pos_Old / pos_Young / neg_Old / neg_Young
enrich_tfs_pair <- function(grn_list) {
    res <- list()
    for (ct in names(grn_list)) {
        g <- grn_list[[ct]]
        tfs <- unique(g$TF)
        pop <- nrow(g)
        rows <- lapply(tfs, function(tf) {
            n <- sum(g$TF == tf)
            rec <- data.frame(TF = tf, stringsAsFactors = FALSE)
            for (tg in tags) {
                K <- sum(g$combined_Tag == tg)
                k <- sum(g$TF == tf & g$combined_Tag == tg)
                p <- phyper(k - 1, K, pop - K, n, lower.tail = FALSE)
                rec[[paste0(tg, "__p")]] <- p
                rec[[paste0(tg, "__cov")]] <- if (n > 0) k / n else NA_real_
                rec[[paste0(tg, "__k")]] <- k
                rec[[paste0(tg, "__K")]] <- K
            }
            rec$n_TF_pairs <- n
            rec$N_pairs <- pop
            rec
        })
        tab <- do.call(rbind, rows)
        tab$celltype <- ct
        res[[ct]] <- tab
        message(sprintf("  enrich %s: TFs=%d", ct, nrow(tab)))
    }
    res
}

enrich_list_to_df <- function(enrich_list, species) {
    df <- do.call(rbind, lapply(names(enrich_list), function(ct) {
        x <- enrich_list[[ct]]
        x$species <- species
        x
    }))
    rownames(df) <- NULL
    df
}

## 种内：每个 TF 在几个 CT 对 pos_Old / pos_Young 显著
within_overlap_mats <- function(enrich_list, pcut = 0.05) {
    cts <- names(enrich_list)
    tfs <- unique(unlist(lapply(enrich_list, function(x) x$TF)))
    mat_o <- matrix(0L, nrow = length(tfs), ncol = length(cts),
                    dimnames = list(tfs, cts))
    mat_y <- mat_o
    for (ct in cts) {
        tab <- enrich_list[[ct]]
        hit_o <- tab$TF[!is.na(tab$pos_Old__p) & tab$pos_Old__p < pcut]
        hit_y <- tab$TF[!is.na(tab$pos_Young__p) & tab$pos_Young__p < pcut]
        mat_o[hit_o, ct] <- 1L
        mat_y[hit_y, ct] <- 1L
    }
    to_df <- function(m) {
        d <- as.data.frame(m)
        d$TF <- rownames(d)
        d$sum <- rowSums(m)
        d[order(-d$sum, d$TF), ]
    }
    list(pos_Old = to_df(mat_o), pos_Young = to_df(mat_y))
}

## 跨物种：同一 CT、同一方向（pos_Old 或 pos_Young）的显著 TF，用 ortholog 对齐
prepare_cross_dot <- function(H, M, Z, ct, tag, ortho) {
    get_mat <- function(x) {
        if (is.null(x) || !ct %in% names(x)) return(NULL)
        mat <- x[[ct]]
        pcol <- paste0(tag, "__p")
        ccol <- paste0(tag, "__cov")
        if (!all(c("TF", pcol, ccol) %in% names(mat))) return(NULL)
        mat[, c("TF", pcol, ccol), drop = FALSE]
    }
    H_dat <- get_mat(H); M_dat <- get_mat(M); Z_dat <- get_mat(Z)
    if (is.null(H_dat) || is.null(M_dat) || is.null(Z_dat)) return(data.frame())
    H_sel <- H_dat$TF[H_dat[[2]] < p_cut]
    M_sel <- M_dat$TF[M_dat[[2]] < p_cut]
    Z_sel <- Z_dat$TF[Z_dat[[2]] < p_cut]
    selections <- list(H = H_sel, M = M_sel, Z = Z_sel)
    specs <- list(
        HMZ = list(df = ortho$HMZ, cols = c("human", "mouse", "zebrafish"),
                   idx = "HMZ_index", sp = c("H", "M", "Z")),
        HM  = list(df = ortho$HM,  cols = c("human", "mouse"),
                   idx = "HM_index",  sp = c("H", "M")),
        HZ  = list(df = ortho$HZ,  cols = c("human", "zebrafish"),
                   idx = "HZ_index",  sp = c("H", "Z")),
        MZ  = list(df = ortho$MZ,  cols = c("mouse", "zebrafish"),
                   idx = "MZ_index",  sp = c("M", "Z"))
    )
    lst <- lapply(names(specs), function(cl) {
        s <- specs[[cl]]
        df <- s$df
        if (is.null(df) || nrow(df) == 0) return(NULL)
        if (!s$idx %in% names(df)) {
            df[[s$idx]] <- apply(df[, s$cols, drop = FALSE], 1, paste, collapse = "__")
        }
        keep <- Reduce(`&`, lapply(seq_along(s$cols), function(i) {
            df[[s$cols[i]]] %in% selections[[s$sp[i]]]
        }))
        if (!any(keep)) return(NULL)
        ids <- df[[s$idx]][keep]
        do.call(rbind, lapply(s$sp, function(sp) {
            data.frame(index = ids, sp = sp, genes = ids, class = cl,
                       pvalue = NA_real_, cov = NA_real_,
                       stringsAsFactors = FALSE)
        }))
    })
    lst <- Filter(Negate(is.null), lst)
    if (length(lst) == 0) return(data.frame())
    res <- do.call(rbind, lst)
    ## HM/HZ/MZ 里已属于 HMZ 的对去掉，避免重复
    if (any(res$class == "HMZ")) {
        hmz_ids <- unique(res$index[res$class == "HMZ"])
        drop <- res$class != "HMZ" & vapply(seq_len(nrow(res)), function(i) {
            bits <- unlist(strsplit(res$index[i], "__"))
            any(vapply(hmz_ids, function(id) {
                all(bits %in% unlist(strsplit(id, "__")))
            }, logical(1)))
        }, logical(1))
        res <- res[!drop, , drop = FALSE]
    }
    mats <- list(H = H_dat, M = M_dat, Z = Z_dat)
    for (i in seq_len(nrow(res))) {
        sp <- res$sp[i]
        pos <- regexpr(sp, res$class[i])
        tf <- strsplit(res$genes[i], "__")[[1]][as.integer(pos)]
        row <- mats[[sp]][mats[[sp]]$TF == tf, , drop = FALSE]
        if (nrow(row) > 0) {
            res$pvalue[i] <- as.numeric(row[1, 2])
            res$cov[i] <- as.numeric(row[1, 3])
            res$genes[i] <- tf
        }
    }
    res$logP <- -log(pmax(res$pvalue, .Machine$double.xmin))
    res$order <- match(res$class, c("HMZ", "HM", "HZ", "MZ"))
    res$sp <- factor(res$sp, levels = c("Z", "M", "H"),
                     labels = c("Zebrafish", "Mouse", "Human"))
    ## y 轴标签 / 排序：zebrafish → mouse → human（对齐原 Process_res2）
    z <- m <- h <- character(nrow(res))
    new_idx <- character(nrow(res))
    for (i in seq_len(nrow(res))) {
        bits <- unlist(strsplit(res$index[i], "__"))
        sps <- strsplit(as.character(res$class[i]), "")[[1]]
        if (length(bits) != length(sps)) next
        names(bits) <- sps
        z[i] <- if ("Z" %in% names(bits)) unname(bits[["Z"]]) else ""
        m[i] <- if ("M" %in% names(bits)) unname(bits[["M"]]) else ""
        h[i] <- if ("H" %in% names(bits)) unname(bits[["H"]]) else ""
        new_idx[i] <- paste(c(z[i], m[i], h[i])[c(z[i], m[i], h[i]) != ""],
                            collapse = "__")
    }
    res$index <- new_idx
    res$z_gene <- z
    res$m_gene <- m
    res$h_gene <- h
    res$celltype <- ct
    res$tag <- tag
    rownames(res) <- NULL
    res
}

save_png_pdf <- function(plot, stem, w, h) {
    png <- file.path(out_dir, paste0(stem, ".png"))
    pdf <- file.path(out_dir, paste0(stem, ".pdf"))
    ggsave(png, plot, width = w, height = h, dpi = 300, bg = "white")
    tryCatch(ggsave(pdf, plot, width = w, height = h, device = cairo_pdf, bg = "white"),
             error = function(e) ggsave(pdf, plot, width = w, height = h, bg = "white"))
    message("  saved: ", png)
}

theme_bar <- function() {
    theme_classic(base_size = 12) +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank(),
            legend.title = element_blank()
        )
}

safe_sheet <- function(x) substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", x), 1, 31)

## 不再写 tsv：收集成 Excel tab
tsv_book <- list()
write_table <- function(df, tab, desc = NULL) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        message("  skip empty table: ", tab)
        return(invisible(NULL))
    }
    tab31 <- safe_sheet(tab)
    existing <- if (length(tsv_book)) vapply(tsv_book, `[[`, "", "Tab") else character(0)
    if (tab31 %in% existing) {
        k <- 1
        while (safe_sheet(paste0(substr(tab31, 1, 28), "_", k)) %in% existing) k <- k + 1
        tab31 <- safe_sheet(paste0(substr(tab31, 1, 28), "_", k))
    }
    tsv_book[[length(tsv_book) + 1]] <<- list(
        Tab = tab31,
        Description = if (is.null(desc)) tab else desc,
        N_rows = nrow(df), N_cols = ncol(df), df = df
    )
    message("  excel-tab: ", tab31, "  rows=", nrow(df))
}


########################################################################
# 0) ortholog
########################################################################
ortho_fp <- find_existing(ortho_candidates)
if (is.na(ortho_fp)) stop("找不到 HMZ_ortholog_combined")
HMZ_ortholog_combined <- load_ortholog(ortho_fp)
message("[Ortholog] ", ortho_fp)


########################################################################
# 1) 2026 Final clocks → Old / Young gene sets
########################################################################
message("\n========== 1) 2026 Final clocks → Old/Young gene sets ==========")
clock_genes <- list()
clock_tables <- list()
for (sp in species_order) {
    clock_genes[[sp]] <- list(Old = list(), Young = list())
    for (ct in major_cts) {
        oy <- read_clock_old_young(sp, ct)
        clock_genes[[sp]]$Old[[ct]] <- oy$Old
        clock_genes[[sp]]$Young[[ct]] <- oy$Young
        if (!is.null(oy$table)) clock_tables[[length(clock_tables) + 1]] <- oy$table
        message(sprintf("  %-10s %-4s  Old=%d  Young=%d  src=%s",
                        sp, ct, length(oy$Old), length(oy$Young),
                        if (is.na(oy$src)) "MISSING" else basename(oy$src)))
        if (is.na(oy$src)) {
            warning(sp, "/", ct, " 缺少 Final clock: ", clock_rds_path(sp, ct))
        }
    }
}
clock_df <- do.call(rbind, clock_tables)
rownames(clock_df) <- NULL
write_table(clock_df, "Clock_OldYoung_genesets",
            "2026 Final clocks；Old=coef>0；Young=coef<0。Human=FinalClock_clean5_2026")


########################################################################
# 2) GRNs + 标注
########################################################################
message("\n========== 2) 读 2025 GRNs 并用 2026 clock 标注 ==========")
grn_raw <- list()
grn_ann <- list()
bar_rows <- list()
for (sp in species_order) {
    grn_raw[[sp]] <- load_grn_list(sp)
    message("[annotate] ", sp)
    grn_ann[[sp]] <- annotate_grn_with_clock(
        grn_raw[[sp]], clock_genes[[sp]]$Old, clock_genes[[sp]]$Young
    )
    for (ct in names(grn_ann[[sp]])) {
        tab <- as.data.frame(table(grn_ann[[sp]][[ct]]$combined_Tag),
                             stringsAsFactors = FALSE)
        names(tab) <- c("combined_Tag", "n_pairs")
        if (nrow(tab) == 0) next
        tab$species <- sp
        tab$celltype <- ct
        n_all <- nrow(grn_ann[[sp]][[ct]])
        n_bar <- sum(tab$n_pairs)
        message(sprintf("  bar %s/%s: total GRNs=%d  stacked=%d",
                        sp, ct, n_all, n_bar))
        bar_rows[[length(bar_rows) + 1]] <- tab
    }
}
bar_df <- do.call(rbind, bar_rows)
rownames(bar_df) <- NULL
bar_df$combined_Tag <- factor(as.character(bar_df$combined_Tag), levels = bar_tags)
bar_df$celltype <- factor(bar_df$celltype, levels = major_cts)
write_table(bar_df, "GRN_clock_pair_counts",
            "全部 TF–Gene pair 数：4 类 aging + pos/neg_NODIFF（柱高=该 CT 总 GRN 数）")


########################################################################
# 3) TF enrichment
########################################################################
message("\n========== 3) TF–Gene pair 水平 TF enrichment ==========")
enrich_list <- list()
enrich_rows <- list()
for (sp in species_order) {
    message("[enrich] ", sp)
    enrich_list[[sp]] <- enrich_tfs_pair(grn_ann[[sp]])
    enrich_rows[[sp]] <- enrich_list_to_df(enrich_list[[sp]], sp)
}
enrich_df <- do.call(rbind, enrich_rows)
rownames(enrich_df) <- NULL
write_table(enrich_df, "TF_enrich_clock_all",
            "TF–Gene pair 水平超几何：每个 TF × 四个 combined_Tag 的 p / coverage / k")


########################################################################
# 4) 种内 overlap
########################################################################
message("\n========== 4) 种内 TF overlap（共享 CT 数） ==========")
within_list <- list()
within_long <- list()
for (sp in species_order) {
    within_list[[sp]] <- within_overlap_mats(enrich_list[[sp]], pcut = p_cut)
    for (dirn in c("pos_Old", "pos_Young")) {
        d <- within_list[[sp]][[dirn]]
        d$species <- sp
        d$direction <- dirn
        within_long[[length(within_long) + 1]] <- d
        write_table(d, sprintf("Within_%s_%s", sp, dirn),
                    sprintf("种内 %s %s：每个 TF 在哪些 CT 显著；sum=显著 CT 数", sp, dirn))
    }
}
within_df <- do.call(rbind, within_long)
rownames(within_df) <- NULL


########################################################################
# 5) 跨物种 overlap
########################################################################
message("\n========== 5) 跨物种 TF overlap（每个 cell type） ==========")
cross_rows <- list()
for (ct in major_cts) {
    for (tag in c("pos_Old", "pos_Young")) {
        tmp <- prepare_cross_dot(
            enrich_list$Human, enrich_list$Mouse, enrich_list$Zebrafish,
            ct, tag, HMZ_ortholog_combined
        )
        if (nrow(tmp) == 0) {
            message(sprintf("  %s / %s: no overlap", ct, tag))
            next
        }
        message(sprintf("  %s / %s: rows=%d classes=%s",
                        ct, tag, nrow(tmp),
                        paste(unique(tmp$class), collapse = ",")))
        cross_rows[[length(cross_rows) + 1]] <- tmp
        write_table(tmp, sprintf("Cross_%s_%s", ct, tag),
                    sprintf("跨物种 %s %s：HMZ/HM/HZ/MZ ortholog 对齐", ct, tag))
    }
}
cross_df <- if (length(cross_rows)) do.call(rbind, cross_rows) else data.frame()
if (nrow(cross_df) > 0) rownames(cross_df) <- NULL


########################################################################
# 6) 图
########################################################################
message("\n========== 6) 出图 ==========")
fill_cols <- c(
    pos_NODIFF = "#D9D9D9",
    neg_NODIFF = "#BDBDBD",
    pos_Old    = "#E98B8B",
    pos_Young  = "#99BADD",
    neg_Old    = "#D3EEF4",
    neg_Young  = "#FBE1E8"
)
stack_bottom_to_top <- c("neg_Young", "neg_Old", "pos_Young", "pos_Old",
                         "neg_NODIFF", "pos_NODIFF")
legend_breaks <- c("pos_Old", "pos_Young", "neg_Old", "neg_Young",
                   "pos_NODIFF", "neg_NODIFF")
legend_labels <- c("Activate Old", "Activate Young", "Repress Old", "Repress Young",
                   "Other Activate", "Other Repress")

## A) 三物种合并 barplot（对齐 GRNs_Plots）
if (nrow(bar_df) > 0) {
    bar_plot <- bar_df
    bar_plot$combined_Tag <- factor(as.character(bar_plot$combined_Tag),
                                    levels = stack_bottom_to_top)
    bar_plot$species <- factor(bar_plot$species, levels = species_order)
    bar_plot$celltype <- factor(bar_plot$celltype, levels = major_cts)
    bar_plot <- bar_plot[order(bar_plot$species, bar_plot$celltype, bar_plot$combined_Tag), ]
    stack_tot <- aggregate(n_pairs ~ celltype + species, data = bar_plot, FUN = sum)
    ymax_use <- tapply(stack_tot$n_pairs, stack_tot$species, max, na.rm = TRUE) * 1.05
    dummy <- do.call(rbind, lapply(names(ymax_use), function(sp) {
        data.frame(
            species = factor(sp, levels = species_order),
            celltype = factor(major_cts[1], levels = major_cts),
            n_pairs = as.numeric(ymax_use[[sp]]),
            stringsAsFactors = FALSE
        )
    }))
    p <- ggplot(bar_plot, aes(x = celltype, y = n_pairs, fill = combined_Tag)) +
        geom_col(color = "black", width = 0.8, size = 0.25) +
        geom_blank(data = dummy, aes(x = celltype, y = n_pairs), inherit.aes = FALSE) +
        facet_wrap(~species, nrow = 1, scales = "free_y") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(
            values = fill_cols,
            limits = stack_bottom_to_top,
            breaks = legend_breaks,
            labels = legend_labels,
            drop = FALSE, name = NULL
        ) +
        labs(x = NULL, y = NULL) +
        theme_classic(base_size = 11) +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size = unit(0.35, "cm")
        ) +
        guides(fill = guide_legend(nrow = 2))
    save_png_pdf(p, "HMZ_GRNs_clock_bar_2026Aug", 8.5, 3.2)
}

## B) TF number overlap barplot（对齐 GRNs_Plots / GRNs_for_aging_HMZ.R）
##    x = 该 TF 在几个 CT 显著；y = 这样的 TF 个数
##    合并图：上 Activate Old / 下 Activate Young；左→右 Zebrafish Mouse Human
##    另存每个物种 Old / Young 单图
## 种内 overlap 最高的 TF；不足 min_n 个则加上第二高一档
top_within_overlap_tfs <- function(overlap_df, min_n = 5L, min_overlap = 2L) {
    meta <- c("TF", "sum", "species", "direction")
    ct_cols <- setdiff(names(overlap_df), meta)
    d2 <- overlap_df[!is.na(overlap_df$sum) & overlap_df$sum >= min_overlap, ,
                     drop = FALSE]
    if (nrow(d2) == 0) {
        return(list(df = NULL, top_n = NA_integer_, extra_n = NA_integer_,
                    ct_cols = ct_cols))
    }
    sums <- sort(unique(as.integer(d2$sum)), decreasing = TRUE)
    top_n <- sums[1]
    hit <- d2[d2$sum == top_n, , drop = FALSE]
    extra_n <- NA_integer_
    if (nrow(hit) < min_n && length(sums) >= 2L) {
        extra_n <- sums[2]
        hit <- rbind(hit, d2[d2$sum == extra_n, , drop = FALSE])
    }
    hit <- hit[order(-hit$sum, hit$TF), ]
    rownames(hit) <- NULL
    list(df = hit, top_n = top_n, extra_n = extra_n, ct_cols = ct_cols)
}

print_top_within_overlap_tfs <- function(overlap_df, species, dirn,
                                         min_n = 5L, min_overlap = 2L) {
    res <- top_within_overlap_tfs(overlap_df, min_n = min_n,
                                  min_overlap = min_overlap)
    if (is.null(res$df) || nrow(res$df) == 0) {
        message(sprintf("  [top TFs] %s %s: 无 overlap>=%d 的 TF",
                        species, dirn, min_overlap))
        return(invisible(NULL))
    }
    n_top <- sum(res$df$sum == res$top_n)
    extra_txt <- if (!is.na(res$extra_n)) {
        sprintf("；不足%d个，接着打印共享 %d CT (n=%d)",
                min_n, res$extra_n, sum(res$df$sum == res$extra_n))
    } else ""
    message(sprintf("  [top TFs] %s %s: 最高 overlap = %d CT (n=%d)%s",
                    species, dirn, res$top_n, n_top, extra_txt))
    out_rows <- list()
    for (i in seq_len(nrow(res$df))) {
        cts <- res$ct_cols[as.integer(res$df[i, res$ct_cols]) == 1]
        message(sprintf("    %-16s  n_CT=%d  [%s]",
                        res$df$TF[i], res$df$sum[i], paste(cts, collapse = ",")))
        out_rows[[i]] <- data.frame(
            species = species, direction = dirn,
            TF = res$df$TF[i], n_CT = res$df$sum[i],
            celltypes = paste(cts, collapse = ","),
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, out_rows)
}

overlap_count_df <- function(overlap_df, species, class_label, min_overlap = 2) {
    df <- as.data.frame(table(overlap_df$sum), stringsAsFactors = FALSE)
    names(df) <- c("n_ct", "Freq")
    df$n_ct <- as.character(df$n_ct)
    df <- df[as.numeric(df$n_ct) >= min_overlap, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    data.frame(
        species = species, class = class_label,
        n_ct = df$n_ct, Freq = df$Freq,
        stringsAsFactors = FALSE
    )
}

overlap_parts <- list()
top_tf_parts <- list()
message("\n----- 种内 overlap 最高 TF（画 bar 时打印；不足5个打第二档） -----")
for (sp in species_order) {
    for (dirn in c("pos_Old", "pos_Young")) {
        cls <- if (dirn == "pos_Old") "Activate Old" else "Activate Young"
        ov <- within_list[[sp]][[dirn]]
        top_tab <- print_top_within_overlap_tfs(ov, sp, dirn, min_n = 5L,
                                                min_overlap = 2L)
        if (!is.null(top_tab) && nrow(top_tab) > 0)
            top_tf_parts[[length(top_tf_parts) + 1]] <- top_tab
        tab <- overlap_count_df(ov, sp, cls, min_overlap = 2)
        if (is.null(tab)) {
            message("  skip TF-overlap counts: ", sp, " ", dirn)
            next
        }
        overlap_parts[[length(overlap_parts) + 1]] <- tab
        write_table(tab, sprintf("%s_%s_CT_counts", sp,
                                if (dirn == "pos_Old") "old" else "young"),
                    sprintf("%s %s：共享 CT 数 → TF 个数（barplot 数据）", sp, cls))
        fill <- if (dirn == "pos_Old") "pink" else "lightblue"
        ymax <- max(tab$Freq) * 1.25
        p1 <- ggplot(tab, aes(x = factor(as.integer(n_ct)), y = Freq)) +
            geom_bar(stat = "identity", fill = fill, color = "black", width = 0.7) +
            geom_text(aes(label = Freq), vjust = -0.5, size = 5) +
            scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
            theme_classic(base_size = 14) +
            theme(axis.title = element_blank(),
                  axis.text.x = element_text(size = 18, color = "black"),
                  axis.text.y = element_text(size = 14, color = "black"))
        save_png_pdf(p1, sprintf("%s_%s_overlap_2026Aug", sp,
                                 if (dirn == "pos_Old") "old" else "young"),
                     3, 2.5)
    }
}
if (length(overlap_parts)) {
    odf <- do.call(rbind, overlap_parts)
    panel_levels <- c(
        "Activate Old | Zebrafish", "Activate Old | Mouse", "Activate Old | Human",
        "Activate Young | Zebrafish", "Activate Young | Mouse", "Activate Young | Human"
    )
    odf$panel <- factor(paste(odf$class, odf$species, sep = " | "), levels = panel_levels)
    odf$n_ct <- factor(odf$n_ct, levels = as.character(sort(unique(as.numeric(as.character(odf$n_ct))))))
    ymax_df <- aggregate(Freq ~ panel, data = odf, FUN = max)
    ymax_df$Freq <- ymax_df$Freq * 1.25
    dummy <- ymax_df
    dummy$n_ct <- factor(levels(odf$n_ct)[1], levels = levels(odf$n_ct))
    dummy$class <- gsub(" \\| .*", "", as.character(dummy$panel))
    p <- ggplot(odf, aes(x = n_ct, y = Freq, fill = class)) +
        geom_col(color = "black", width = 0.7) +
        geom_text(aes(label = Freq), vjust = -0.35, size = 5) +
        geom_blank(data = dummy, aes(x = n_ct, y = Freq, fill = class)) +
        facet_wrap(~panel, nrow = 2, ncol = 3, scales = "free") +
        scale_fill_manual(values = c("Activate Old" = "pink", "Activate Young" = "lightblue"),
                          guide = "none") +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = NULL, y = NULL) +
        theme_classic(base_size = 14) +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = 8),
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 12)
        )
    save_png_pdf(p, "HMZ_TF_overlap_bar_2026Aug", 6.5, 3.8)
} else {
    message("  [WARN] 没有 TF overlap 计数，跳过 HMZ_TF_overlap_bar")
}
if (length(top_tf_parts)) {
    top_tf_df <- do.call(rbind, top_tf_parts)
    rownames(top_tf_df) <- NULL
    write_table(top_tf_df, "Within_top_overlap_TFs",
                "画 overlap bar 时打印的种内最高 overlap TF；不足5个含第二档")
}

## C) 跨物种点图：每个 cell type × pos_Old / pos_Young
cross_colors <- c(pos_Old = "red", pos_Young = "blue")
if (nrow(cross_df) > 0) {
    for (ct in major_cts) {
        for (tag in c("pos_Old", "pos_Young")) {
            sub <- cross_df[cross_df$celltype == ct & cross_df$tag == tag, ]
            if (nrow(sub) == 0) {
                message("  skip conserved (empty): ", ct, " ", tag)
                next
            }
            hi <- cross_colors[[tag]]
            sub$logP[sub$logP > 5] <- 5
            ## 从上到下：HMZ → HM → HZ → MZ；组内按 fish → mouse → human
            sub$class <- factor(as.character(sub$class),
                                levels = c("HMZ", "HM", "HZ", "MZ"))
            ord <- unique(sub$index[order(as.integer(sub$class),
                                          sub$z_gene, sub$m_gene, sub$h_gene,
                                          sub$index)])
            sub$index <- factor(sub$index, levels = rev(ord))
            p <- ggplot(sub, aes(x = sp, y = index)) +
                geom_point(aes(size = cov, color = logP)) +
                scale_color_continuous(low = "grey", high = hi,
                                       limits = c(0, 5),
                                       name = expression(-log(italic(P))),
                                       guide = guide_colorbar(order = 1)) +
                scale_size_continuous(name = "Coverage",
                                      guide = guide_legend(order = 2)) +
                theme_classic(base_size = 12) +
                theme(
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
                    axis.line = element_line(color = "black"),
                    axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
                    axis.text.y = element_text(color = "black", size = 8, face = "italic"),
                    axis.title = element_blank()
                )
            stem <- sprintf("HMZ_%s_%s_conserved_2026Aug", ct,
                            if (tag == "pos_Old") "posOld" else "posYoung")
            save_png_pdf(p, stem, 3.5, 4)
        }
    }
}


########################################################################
# 7) C1–C10 → TF enrich bubble（3 张）
#    1 Human: C2_down = (LE∩edgeR)_male ∪ (LE∩edgeR)_female（每 CT；不用 C1）
#    2 Zebrafish: C3 down
#    3 Mouse: C5 up
########################################################################
message("\n========== 7) C1-C10 gene-set TF overlap bubbles ==========")
gs_dir <- "/projects/hmz-aging/Aging_GeneSets_LEedgeR_2026"
gs_fp <- file.path(gs_dir, "ALL_GeneSets_summary_2026.tsv")
if (!file.exists(gs_fp)) {
    alts <- list.files(gs_dir, pattern = "_summary\\.tsv$", full.names = TRUE)
    alts <- alts[!grepl("^ALL_", basename(alts))]
    if (length(alts)) {
        message("[GeneSets] 无 ALL summary，拼接: ", paste(basename(alts), collapse = ", "))
        gs_sum <- do.call(rbind, lapply(alts, function(f)
            read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)))
    } else {
        gs_sum <- NULL
        warning("找不到 ", gs_fp, "；跳过 C1-C10 bubble。请先跑 Build_Aging_GeneSets_LEedgeR_2026.R")
    }
} else {
    message("[GeneSets] ", gs_fp)
    gs_sum <- read.delim(gs_fp, stringsAsFactors = FALSE, check.names = FALSE)
}

parse_genes <- function(x) {
    if (length(x) == 0 || all(is.na(x))) return(character(0))
    unique(unlist(strsplit(as.character(x[!is.na(x) & nzchar(x)]), ";")))
}

gs_native <- function(species, label, celltype) {
    if (is.null(gs_sum)) return(character(0))
    rows <- gs_sum[gs_sum$species == species & gs_sum$label == label &
                       gs_sum$celltype == celltype, , drop = FALSE]
    if (nrow(rows) == 0) {
        rows <- gs_sum[gs_sum$species == species & gs_sum$label == label &
                           tolower(gs_sum$celltype) == tolower(celltype), ,
                       drop = FALSE]
    }
    if (nrow(rows) == 0) return(character(0))
    col <- if ("genes_native" %in% names(rows)) "genes_native" else "genes_human"
    parse_genes(rows[[col]])
}

## Human C2 基因集：见下方 h_gs（male/female 各自 C2_down 再 union）

enrich_tf_geneset <- function(grn_df, geneset, min_targets = 5L) {
    empty <- data.frame(TF = character(0), pvalue = numeric(0),
                        coverage = numeric(0), k = integer(0),
                        K = integer(0), n_TF = integer(0), N = integer(0),
                        stringsAsFactors = FALSE)
    if (is.null(grn_df) || nrow(grn_df) == 0 || length(geneset) == 0) return(empty)
    bg <- unique(as.character(grn_df$Gene))
    hit_bg <- intersect(bg, unique(as.character(geneset)))
    N <- length(bg)
    K <- length(hit_bg)
    if (N < 2 || K < 1) return(empty)
    tfs <- unique(as.character(grn_df$TF))
    rows <- lapply(tfs, function(tf) {
        tg <- unique(as.character(grn_df$Gene[grn_df$TF == tf]))
        n <- length(tg)
        if (n < min_targets) return(NULL)
        k <- length(intersect(tg, hit_bg))
        p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
        data.frame(TF = tf, pvalue = p, coverage = if (n > 0) k / n else NA_real_,
                   k = k, K = K, n_TF = n, N = N, stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(empty)
    do.call(rbind, rows)
}

## 按所有 cell type 的平均 -log10(P) 排序（缺测 CT 记 0），取 top_n
select_tfs_mean_logp <- function(enrich_df, top_n = 20L) {
    empty <- data.frame(TF = character(0), n_sig_CT = integer(0),
                        mean_neglogP = numeric(0), stringsAsFactors = FALSE)
    if (is.null(enrich_df) || nrow(enrich_df) == 0) return(empty)
    tfs <- unique(as.character(enrich_df$TF))
    rank_df <- do.call(rbind, lapply(tfs, function(tf) {
        nlp <- vapply(major_cts, function(ct) {
            hit <- enrich_df[enrich_df$TF == tf & enrich_df$celltype == ct, ]
            if (nrow(hit) == 0 || is.na(hit$pvalue[1])) return(0)
            as.numeric(-log10(max(hit$pvalue[1], 1e-300)))
        }, numeric(1))
        n_sig <- if ("k" %in% names(enrich_df)) {
            sum(vapply(major_cts, function(ct) {
                hit <- enrich_df[enrich_df$TF == tf & enrich_df$celltype == ct, ]
                nrow(hit) > 0 && !is.na(hit$pvalue[1]) && hit$pvalue[1] < p_cut &&
                    !is.na(hit$k[1]) && hit$k[1] > 0
            }, logical(1)))
        } else {
            sum(nlp > -log10(p_cut))
        }
        data.frame(TF = tf, n_sig_CT = as.integer(n_sig),
                   mean_neglogP = mean(nlp), stringsAsFactors = FALSE)
    }))
    rank_df <- rank_df[order(-rank_df$mean_neglogP, rank_df$TF), ]
    if (nrow(rank_df) > top_n) rank_df <- rank_df[seq_len(top_n), , drop = FALSE]
    message("  selected TFs by mean -log10(P) across all CT (", nrow(rank_df), "): ",
            paste(rank_df$TF, collapse = ", "))
    rank_df
}

## 跨 CT conserved：先按显著 CT 数，再按平均 -log10(P)
select_top_overlap_tfs <- function(enrich_df, top_n = 18L, pcut = 0.05,
                                   min_shared = 2L) {
    empty <- data.frame(TF = character(0), n_sig_CT = integer(0),
                        mean_neglogP = numeric(0), stringsAsFactors = FALSE)
    if (is.null(enrich_df) || nrow(enrich_df) == 0) return(empty)
    sig <- enrich_df[!is.na(enrich_df$pvalue) & enrich_df$pvalue < pcut &
                         enrich_df$k > 0, , drop = FALSE]
    if (nrow(sig) == 0) return(empty)
    n_ct <- tapply(sig$celltype, sig$TF, function(x) length(unique(x)))
    mean_nlp <- tapply(-log10(pmax(sig$pvalue, 1e-300)), sig$TF, mean)
    rank_df <- data.frame(
        TF = names(n_ct),
        n_sig_CT = as.integer(n_ct),
        mean_neglogP = as.numeric(mean_nlp[names(n_ct)]),
        stringsAsFactors = FALSE
    )
    rank_df <- rank_df[order(-rank_df$n_sig_CT, -rank_df$mean_neglogP, rank_df$TF), ]
    shared <- rank_df[rank_df$n_sig_CT >= min_shared, , drop = FALSE]
    if (nrow(shared) >= 8L) {
        rank_df <- shared
        message(sprintf("  conserved TFs (>=%d CT)=%d", min_shared, nrow(rank_df)))
    } else {
        message(sprintf("  conserved (>=%d CT) 只有 %d 个，改用全部显著 TF 按 CT 数排",
                        min_shared, nrow(shared)))
    }
    if (nrow(rank_df) > top_n) rank_df <- rank_df[seq_len(top_n), , drop = FALSE]
    message("  selected TFs (", nrow(rank_df), "): ", paste(rank_df$TF, collapse = ", "))
    rank_df
}

plot_ct_tf_bubble <- function(enrich_df, tf_rank, title, stem = NULL,
                              pcut = 0.05, do_save = TRUE,
                              y_text_size = 15,
                              fill_cols = c("#E8DAEF", "#AF7AC5", "#8E44AD", "#5B2C6F"),
                              square = FALSE) {
    if (is.null(tf_rank) || (is.data.frame(tf_rank) && nrow(tf_rank) == 0)) {
        message("  skip bubble (no TFs): ", if (is.null(stem)) title else stem)
        return(invisible(NULL))
    }
    tfs <- if (is.data.frame(tf_rank)) as.character(tf_rank$TF) else as.character(tf_rank)
    plot_df <- enrich_df[enrich_df$TF %in% tfs, , drop = FALSE]
    if (nrow(plot_df) == 0) return(invisible(NULL))
    plot_df$sig <- !is.na(plot_df$pvalue) & plot_df$pvalue < pcut &
        !is.na(plot_df$k) & plot_df$k > 0
    plot_df$neglogP <- -log10(pmax(plot_df$pvalue, 1e-300))
    pts <- plot_df[!is.na(plot_df$sig) & plot_df$sig, , drop = FALSE]
    if (nrow(pts) == 0) {
        message("  skip bubble (no sig points): ", if (is.null(stem)) title else stem)
        return(invisible(NULL))
    }
    pts$TF <- factor(pts$TF, levels = rev(tfs))
    pts$celltype <- factor(pts$celltype, levels = major_cts)
    pts$neglogP <- pmin(pmax(pts$neglogP, 0), 5)
    cov_vals <- pts$coverage[is.finite(pts$coverage) & pts$coverage >= 0]
    cov_max <- if (length(cov_vals)) max(cov_vals) else 0.01
    if (!is.finite(cov_max) || cov_max <= 0) cov_max <- 0.01
    p <- ggplot(pts, aes(x = celltype, y = TF)) +
        geom_point(aes(size = coverage, fill = neglogP),
                   shape = 21, color = "black", stroke = 0.35) +
        scale_fill_gradientn(
            colours = fill_cols,
            limits = c(0, 5),
            oob = scales::squish,
            name = expression(-log[10]~italic(P))
        ) +
        scale_size_continuous(name = "Coverage", range = c(3.4, 8.0),
                              limits = c(0, cov_max * 1.05)) +
        scale_x_discrete(drop = FALSE, limits = major_cts) +
        scale_y_discrete(drop = FALSE, limits = rev(tfs)) +
        theme_classic(base_size = 16) +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
            axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 15),
            axis.text.y = element_text(face = "italic", color = "black", size = y_text_size),
            axis.title = element_blank(),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 11),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.margin = margin(4, 4, 4, 4),
            aspect.ratio = if (isTRUE(square)) 1 else NULL
        ) +
        labs(title = title)
    if (isTRUE(square)) {
        panel_in <- max(3.2, 0.13 * max(length(tfs), length(major_cts)) + 1.2)
        fig_w <- panel_in + 1.8
        fig_h <- panel_in + 1.2
    } else {
        fig_w <- 4.6
        fig_h <- max(3.4, 0.16 * length(tfs) + 1.6)
    }
    if (isTRUE(do_save) && !is.null(stem)) {
        save_png_pdf(p, stem, fig_w, fig_h)
        if (is.data.frame(tf_rank) && nrow(tf_rank) > 0) {
            write_table(tf_rank, paste0(stem, "_TFs"),
                        paste(title, "| 跨 CT conserved TFs"))
        }
    }
    invisible(p)
}

run_geneset_enrich <- function(species, geneset_by_ct, title, tab,
                               rank_mode = c("overlap", "mean_logp"),
                               top_n = NULL) {
    message("[C-enrich] ", title)
    rows <- list()
    gs_n <- integer(length(major_cts))
    names(gs_n) <- major_cts
    for (ct in major_cts) {
        gs <- if (ct %in% names(geneset_by_ct)) geneset_by_ct[[ct]] else character(0)
        gs_n[ct] <- length(gs)
        grn <- if (!is.null(grn_raw[[species]]) && ct %in% names(grn_raw[[species]]))
            grn_raw[[species]][[ct]] else NULL
        message(sprintf("  %s/%s geneset=%d", species, ct, length(gs)))
        en <- enrich_tf_geneset(grn, gs)
        if (is.null(en) || nrow(en) == 0) next
        en$species <- species
        en$celltype <- ct
        rows[[ct]] <- en
    }
    if (!length(rows)) {
        message("  no enrich rows: ", title)
        return(invisible(NULL))
    }
    en_df <- do.call(rbind, rows)
    rownames(en_df) <- NULL
    write_table(en_df, tab, title)
    rank_mode <- match.arg(rank_mode)
    if (is.null(top_n)) top_n <- if (rank_mode == "mean_logp") 25L else 18L
    tf_rank <- if (rank_mode == "mean_logp") {
        select_tfs_mean_logp(en_df, top_n = top_n)
    } else {
        select_top_overlap_tfs(en_df, top_n = top_n, pcut = p_cut)
    }
    invisible(list(enrich = en_df, tfs = tf_rank, n_genes = gs_n))
}

if (!is.null(gs_sum)) {
    gs_size_rows <- list()

    ## 1) Human C2_down：male / female 各自 LE∩edgeR，再取并集（不用 C1）
    h_size_detail <- list()
    h_gs <- setNames(lapply(major_cts, function(ct) {
        g_m <- gs_native("Human", "C2_down", paste0(ct, "_male"))
        g_f <- gs_native("Human", "C2_down", paste0(ct, "_female"))
        g_u <- unique(c(g_m, g_f))
        message(sprintf("  Human C2_down %s: male=%d  female=%d  union=%d  intersect=%d",
                        ct, length(g_m), length(g_f), length(g_u),
                        length(intersect(g_m, g_f))))
        h_size_detail[[length(h_size_detail) + 1]] <<- data.frame(
            species = "Human", celltype = ct, label = "C2_down",
            n_male = length(g_m), n_female = length(g_f),
            n_union = length(g_u), n_intersect = length(intersect(g_m, g_f)),
            stringsAsFactors = FALSE
        )
        g_u
    }), major_cts)
    for (ct in major_cts) {
        gs_size_rows[[length(gs_size_rows) + 1]] <- data.frame(
            species = "Human", celltype = ct, label = "C2_down_mf_union",
            n_genes = length(h_gs[[ct]]), stringsAsFactors = FALSE
        )
    }
    if (length(h_size_detail)) {
        write_table(do.call(rbind, h_size_detail), "Human_C2_mf_sizes",
                    "Human C2_down：各 CT 上 male / female LE∩edgeR 大小及并集/交集大小")
    }
    h_gene_rows <- lapply(major_cts, function(ct) {
        g <- h_gs[[ct]]
        if (!length(g)) return(NULL)
        data.frame(celltype = ct, gene = g, stringsAsFactors = FALSE)
    })
    h_gene_rows <- Filter(Negate(is.null), h_gene_rows)
    if (length(h_gene_rows)) {
        write_table(do.call(rbind, h_gene_rows), "Human_C2_mf_genes",
                    "Human C2_down genes：(LE∩edgeR)_male ∪ (LE∩edgeR)_female（每 CT）")
    }
    h_res <- run_geneset_enrich(
        "Human", h_gs,
        "Human | C2 down | male ∪ female (LE∩edgeR)",
        "Human_C2_mf_enrich",
        rank_mode = "mean_logp", top_n = 25L
    )
    if (!is.null(h_res)) {
        plot_ct_tf_bubble(h_res$enrich, h_res$tfs,
                          "Human | C2 down | male ∪ female (LE∩edgeR)",
                          "Human_C2_LE_mf_bubble_2026Aug",
                          y_text_size = 9,
                          fill_cols = c("#F5EEF8", "#D2B4DE", "#BB8FCE", "#A569BD"),
                          square = TRUE)
    }

    ## 2) Zebrafish C3 down
    z_gs <- setNames(lapply(major_cts, function(ct)
        gs_native("Zebrafish", "C3_down", ct)), major_cts)
    for (ct in major_cts) {
        gs_size_rows[[length(gs_size_rows) + 1]] <- data.frame(
            species = "Zebrafish", celltype = ct, label = "C3_down",
            n_genes = length(z_gs[[ct]]), stringsAsFactors = FALSE
        )
    }
    z_res <- run_geneset_enrich(
        "Zebrafish", z_gs, "Zebrafish | C3 down", "Zebrafish_C3_enrich",
        rank_mode = "mean_logp", top_n = 25L
    )
    if (!is.null(z_res)) {
        plot_ct_tf_bubble(z_res$enrich, z_res$tfs,
                          "Zebrafish | C3 down",
                          "Zebrafish_C3_down_bubble_2026Aug",
                          y_text_size = 9,
                          fill_cols = c("#F5EEF8", "#D2B4DE", "#BB8FCE", "#A569BD"),
                          square = TRUE)
    }

    ## 3) Mouse C5 up
    m5_gs <- setNames(lapply(major_cts, function(ct)
        gs_native("Mouse", "C5_up", ct)), major_cts)
    for (ct in major_cts) {
        gs_size_rows[[length(gs_size_rows) + 1]] <- data.frame(
            species = "Mouse", celltype = ct, label = "C5_up",
            n_genes = length(m5_gs[[ct]]), stringsAsFactors = FALSE
        )
    }
    m5_res <- run_geneset_enrich("Mouse", m5_gs, "Mouse | C5 up", "Mouse_C5_up_enrich",
                                rank_mode = "mean_logp", top_n = 25L)
    if (!is.null(m5_res)) {
        plot_ct_tf_bubble(m5_res$enrich, m5_res$tfs,
                          "Mouse | C5 up",
                          "Mouse_C5_up_bubble_2026Aug",
                          y_text_size = 9,
                          fill_cols = c("#F5EEF8", "#D2B4DE", "#BB8FCE", "#A569BD"),
                          square = TRUE)
    }

    if (length(gs_size_rows)) {
        gs_size <- do.call(rbind, gs_size_rows)
        rownames(gs_size) <- NULL
        write_table(gs_size, "C1C10_geneset_sizes",
                    "各分析用的基因集大小（Human C2 = male∪female of LE∩edgeR down）")
    }
}


########################################################################
# 8) Excel：原 tsv 各一 tab + Contents
########################################################################
message("\n========== 8) Excel ==========")
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  border = "TopBottomLeftRight", halign = "center")
title_style <- createStyle(textDecoration = "bold", fontSize = 14)
sheet_specs <- tsv_book
push_sheet <- function(tab, desc, df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(invisible(NULL))
    tab31 <- safe_sheet(tab)
    existing <- if (length(sheet_specs)) vapply(sheet_specs, `[[`, "", "Tab") else character(0)
    if (tab31 %in% existing) {
        k <- 1
        while (safe_sheet(paste0(substr(tab31, 1, 28), "_", k)) %in% existing) k <- k + 1
        tab31 <- safe_sheet(paste0(substr(tab31, 1, 28), "_", k))
    }
    sheet_specs[[length(sheet_specs) + 1]] <<- list(
        Tab = tab31, Description = desc,
        N_rows = nrow(df), N_cols = ncol(df), df = df
    )
}

clock_sum <- aggregate(gene ~ species + celltype + direction, data = clock_df, FUN = length)
names(clock_sum)[names(clock_sum) == "gene"] <- "n_genes"
push_sheet("Clock_geneset_counts",
           "各物种×celltype 的 Old/Young clock 基因数",
           clock_sum)

sig_long <- do.call(rbind, lapply(tags, function(tg) {
    pcol <- paste0(tg, "__p")
    ccol <- paste0(tg, "__cov")
    kcol <- paste0(tg, "__k")
    sub <- enrich_df[!is.na(enrich_df[[pcol]]) & enrich_df[[pcol]] < p_cut, ]
    if (nrow(sub) == 0) return(NULL)
    data.frame(
        species = sub$species, celltype = sub$celltype, TF = sub$TF,
        tag = tg, pvalue = sub[[pcol]], coverage = sub[[ccol]],
        k = if (kcol %in% names(sub)) sub[[kcol]] else NA_integer_,
        n_TF_pairs = sub$n_TF_pairs, stringsAsFactors = FALSE
    )
}))
if (!is.null(sig_long) && nrow(sig_long) > 0) {
    sig_long <- sig_long[order(sig_long$species, sig_long$celltype, sig_long$tag, sig_long$pvalue), ]
    rownames(sig_long) <- NULL
    push_sheet("TF_enrich_SIG",
               sprintf("pvalue < %s 的显著 TF（四个 tag）", p_cut),
               sig_long)
}

for (sp in species_order) {
    push_sheet(paste0(sp, "_enrich"),
               paste(sp, "完整 TF enrich 表"),
               enrich_df[enrich_df$species == sp, ])
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
          "2026 Final clocks Old/Young → TF enrichment → overlap (2026Aug)",
          startRow = 1, startCol = 1)
addStyle(wb, "Contents", title_style, rows = 1, cols = 1)
writeData(wb, "Contents",
          sprintf("Gene sets: 2026 Final clocks (Human=FinalClock_clean5_2026) | GRNs: 2025 *_GRNs_list_clcl_2025 unique TF-Gene pairs | p<%s | %s",
                  p_cut, format(Sys.time(), "%Y-%m-%d %H:%M")),
          startRow = 2, startCol = 1)
writeData(wb, "Contents",
          "Old = coef>0; Young = coef<0. Clock enrich = TF-Gene pair hypergeometric on pos/neg × Old/Young. Category bubbles (2025 GRN): Human C2_down = (LE∩edgeR)_male ∪ (LE∩edgeR)_female; Zebrafish C3_down; Mouse C5_up.",
          startRow = 3, startCol = 1)
writeData(wb, "Contents", contents, startRow = 5, startCol = 1, headerStyle = hs)

for (i in seq_along(sheet_specs)) {
    spec <- sheet_specs[[i]]
    addWorksheet(wb, spec$Tab)
    writeData(wb, spec$Tab, spec$df, headerStyle = hs)
    freezePane(wb, spec$Tab, firstRow = TRUE)
    setColWidths(wb, spec$Tab, cols = seq_len(min(ncol(spec$df), 14)), widths = "auto")
    link <- sprintf("internal:'%s'!A1", spec$Tab)
    class(link) <- "hyperlink"
    writeData(wb, "Contents", x = link, startRow = 5 + i, startCol = 1)
}
freezePane(wb, "Contents", firstActiveRow = 6)
setColWidths(wb, "Contents", cols = 1:4, widths = c(28, 10, 10, 80))

excel_out <- file.path(out_dir, "Clock_TF_enrich_overlap_2026Aug.xlsx")
saveWorkbook(wb, excel_out, overwrite = TRUE)

message("\n=== Done ===")
message("输出目录: ", out_dir)
message("Excel:     ", excel_out)
message("图:")
message("  HMZ_GRNs_clock_bar_2026Aug.png")
message("  HMZ_TF_overlap_bar_2026Aug.png")
message("  {Zebrafish,Mouse,Human}_{old,young}_overlap_2026Aug.png")
message("  HMZ_{MG,RGC,AC,HC,Rod,Cone,BC}_{posOld,posYoung}_conserved_2026Aug.png")
message("  Human_C2_LE_mf_bubble_2026Aug.png")
message("  Zebrafish_C3_down_bubble_2026Aug.png")
message("  Mouse_C5_up_bubble_2026Aug.png")
