#!/usr/bin/env Rscript
########################################################################
# HMZ DEG 出图（2026）：不重跑 kmeans，接着之前算好的聚类结果
#
# Zebrafish: 只用 2026 revise（DEGs_revised_and_Rank）
#   读 Zebrafish_DEGs_Plot_neworder + Kmeans_order_neworder
#   Young 1–3 / Middle 4–5 / Old 6–9
#   不是第一季 Aging_DEGs.R（k=10, omb2）
# Mouse: DEGs_revised_and_Rank_Mouse
#   Young 1–4 / Middle 5–6 / Old 7–10
# Human: DEGs_revised_and_Rand_Human 第四轮 Kmeans4
#   Young 1–4 / Middle 5–6 / Old 7–9；bar 分 male / female
#   male/female Young-Old overlap：Common + Divergent 热图
#   四组 GO BP bubble（不要 barplot）
#   跨物种 Young/Old overlap：全部 cell type（与 Gene_clusters 同源）；
#     UpSet 图仍只出 MG/Rod/RGC；表含 HM/HZ/MZ/triple
#   另附 edgeR Old vs Young（每物种 × cell type；Human 按 sex）
#   输出全部写入 /projects/hmz-aging/HMZ_DEGs_plots_2026Aug/
#   Excel: HMZ_DEGs_tables_2026Aug.xlsx（Contents in English）
#
# HPC:
#   conda activate seurat_env
#   NCORES=10 Rscript HMZ_DEGs_from_Kmeans_plot_2026.R
#   默认 10 核；三物种热图、F/G/H 各最多 3 路并行
#   GO bubble 需要 clusterProfiler + org.Hs.eg.db；没有则跳过 GO 其余照出
########################################################################

req_pkgs <- c("ComplexHeatmap", "circlize", "ggplot2", "reshape2",
              "grid", "parallel", "openxlsx")
opt_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
pkg_ok <- function(p) requireNamespace(p, quietly = TRUE)
missing_req <- req_pkgs[!vapply(req_pkgs, pkg_ok, logical(1))]
missing_opt <- opt_pkgs[!vapply(opt_pkgs, pkg_ok, logical(1))]
if (length(missing_req)) {
    message("缺少以下 R 包，请先安装后再运行：")
    message(paste0("  - ", missing_req, collapse = "\n"))
}
if (length(missing_opt)) {
    message("可选包未安装（将跳过 GO，其余照出）：")
    message(paste0("  - ", missing_opt, collapse = "\n"))
}
if (length(missing_req)) {
    quit(save = "no", status = 1)
}
if (!length(missing_opt)) message("library check OK")
suppressMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(reshape2)
    library(grid)
})

fish_dir  <- "/projects/hmz-aging/Zebrafish_New_Datasets"
mouse_dir <- "/projects/hmz-aging/Mouse_New_Datasets"
human_dir <- "/projects/hmz-aging/Human_New_Datasets"
human_pb  <- file.path(human_dir, "pseudobulk_sample")
plot_dir  <- "/projects/hmz-aging/HMZ_DEGs_plots_2026Aug"
if (dir.exists(plot_dir)) {
    old <- list.files(plot_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    if (length(old)) {
        message("清空输出文件夹: ", plot_dir, " (", length(old), " 项)")
        unlink(old, recursive = TRUE, force = TRUE)
    }
}
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

ncores <- {
    v <- suppressWarnings(as.integer(Sys.getenv("NCORES", unset = NA_character_)))
    if (is.na(v) || v < 1L) v <- 10L
    nc <- tryCatch(parallel::detectCores(), error = function(e) v)
    if (is.na(nc) || nc < 1L) nc <- v
    min(as.integer(v), as.integer(nc))
}

safe_mclapply <- function(X, FUN, mc.cores = 1L) {
    n <- length(X)
    mc.cores <- max(1L, min(as.integer(mc.cores), n))
    if (n == 0L) return(list())
    if (mc.cores <= 1L || .Platform$OS.type != "unix") {
        return(lapply(X, FUN))
    }
    parallel::mclapply(X, FUN, mc.cores = mc.cores, mc.preschedule = FALSE)
}

celltype_colors <- c(
    MG = "#D11536", RGC = "#AAA9A9", AC = "#F6BA00", HC = "#EF9000",
    BC = "#804537", Rod = "#026AB1", Cone = "#9EA220",
    RPE = "#936DAD", Microglia = "#EC6F64"
)
sex_colors <- c(Male = "#3C5488", Female = "#E64B35")
col_fun <- colorRamp2(c(-2, -1, 0, 1, 2),
                      c("#026AB1", "lightblue", "white", "#EF9000", "#D11536"))

jobs <- list(
    Zebrafish = list(
        out_dir = fish_dir,
        ## 只用 2026 revise（DEGs_revised_and_Rank, k=12）
        ## 不要读第一季 Aging_DEGs.R 的 Zebrafish_DEGs_Plot / Kmeans_order（k=10）
        plot_raw = character(0),
        km_raw = character(0),
        plot_order = c(
            file.path(fish_dir, "Zebrafish_DEGs_Plot_neworder"),
            file.path(fish_dir, "Zebrafish_DEGs_Plot_neworder.rds")
        ),
        km_order = c(
            file.path(fish_dir, "Zebrafish_DEGs_Plot_Kmeans_order_neworder"),
            file.path(fish_dir, "Zebrafish_DEGs_Plot_Kmeans_order_neworder.rds")
        ),
        new_order_down = c(10, 8, 11),
        new_order_middle = c(7, 9),
        new_order_up = c(3, 2, 12, 4),
        rm_clusters = c(1, 5, 6),
        Down = c(1, 2, 3),
        Middle = c(4, 5),
        Up = c(6, 7, 8, 9),
        marks = c("Rod__hsp90ab1", "Rod__nr2e3", "MG__sparc", "Cone__gnat2",
                  "BC__aldoaa", "BC__gpia", "MG__adgra2", "BC__clocka",
                  "MG__b2m", "BC__arntl2", "MG__sash1a"),
        split_sex = FALSE,
        bar_ymax = 5500,
        ht_png_width = 4800,
        stem = "Zebrafish"
    ),
    Mouse = list(
        out_dir = mouse_dir,
        plot_raw = c(
            file.path(mouse_dir, "Mouse_DEGs_Plot_2026"),
            file.path(mouse_dir, "Mouse_DEGs_Plot_2026.rds"),
            file.path(mouse_dir, "Mouse_DEGs_Plot")
        ),
        km_raw = c(
            file.path(mouse_dir, "Mouse_DEGs_Plot_Kmeans_2026"),
            file.path(mouse_dir, "Mouse_DEGs_Plot_Kmeans_2026.rds"),
            file.path(mouse_dir, "Mouse_DEGs_Plot_Kmeans")
        ),
        plot_order = c(
            file.path(mouse_dir, "Mouse_DEGs_Plot_neworder_2026"),
            file.path(mouse_dir, "Mouse_DEGs_Plot_neworder_2026.rds")
        ),
        km_order = c(
            file.path(mouse_dir, "Mouse_DEGs_Plot_Kmeans_order_2026"),
            file.path(mouse_dir, "Mouse_DEGs_Plot_Kmeans_order_2026.rds")
        ),
        new_order_down = c(10, 5, 9, 6),
        new_order_middle = c(3, 2),
        new_order_up = c(8, 7, 1, 4),
        rm_clusters = integer(0),
        Down = c(1, 2, 3, 4),
        Middle = c(5, 6),
        Up = c(7, 8, 9, 10),
        marks = c("MG__Smad9", "MG__Galnt13", "Rod__Prdm1", "MG__Notch1",
                  "MG__Adgra2", "MG__Sash1", "MG__Aqp4", "MG__Stat1",
                  "MG__B2m", "MG__Zfhx4", "MG__Grik3", "MG__Relb",
                  "MG__Apoe", "MG__Clu", "MG__Glul", "Cone__Bmpr1b", "Cone__Pou6f2",
                  "Rod__Col4a3", "Rod__Pou6f2", "MG__A2m", "MG__C4b", "RGC__Nfib"),
        split_sex = FALSE,
        bar_ymax = 5000,
        stem = "Mouse"
    ),
    Human = list(
        out_dir = human_dir,
        plot_raw = c(
            file.path(human_pb, "Human_DEGs_Plot4.rds"),
            file.path(human_dir, "Human_DEGs_Plot4.rds"),
            file.path(human_pb, "Human_DEGs_Plot4")
        ),
        km_raw = c(
            file.path(human_pb, "Human_DEGs_Plot_Kmeans4.rds"),
            file.path(human_dir, "Human_DEGs_Plot_Kmeans4.rds"),
            file.path(human_pb, "Human_DEGs_Plot_Kmeans4")
        ),
        plot_order = c(
            file.path(human_pb, "Human_DEGs_Plot4_order.rds"),
            file.path(human_dir, "Human_DEGs_Plot4_order.rds")
        ),
        km_order = c(
            file.path(human_pb, "Human_DEGs_Plot_Kmeans4_order.rds"),
            file.path(human_dir, "Human_DEGs_Plot_Kmeans4_order.rds")
        ),
        ## 原始 Kmeans4: 2,6,9,5 → Young(1–4)；7,4 → Middle(5–6)；8,3,1 → Old(7–9)
        new_order_down = c(2, 6, 9, 5),
        new_order_middle = c(7, 4),
        new_order_up = c(8, 3, 1),
        rm_clusters = integer(0),
        Down = c(1, 2, 3, 4),
        Middle = c(5, 6),
        Up = c(7, 8, 9),
        marks = c("F_Cone__STAT1", "M_MG__BTN3A3", "M_RPE__BTN3A3", "M_RPE__APOE",
                  "F_AC__EBF1", "F_AC__BTN3A3", "F_RGC__ELOVL6", "M_RGC__CLOCK",
                  "M_RPE__DKK3", "F_RGC__CLOCK"),
        split_sex = TRUE,
        bar_ymax = NA_real_,
        ht_png_width = 5600,
        heatmap_col_mm = 8,
        stem = "Human"
    )
)


########################################################################
# helpers
########################################################################
find_existing <- function(cands) {
    hit <- cands[file.exists(cands)]
    if (length(hit) == 0) NA_character_ else hit[[1]]
}

load_r_object <- function(fp) {
    obj <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (!is.null(obj)) return(obj)
    env <- new.env(parent = emptyenv())
    load(fp, envir = env)
    objs <- ls(envir = env)
    if (length(objs) == 0) stop("空对象: ", fp)
    if (length(objs) == 1L) return(get(objs[[1]], envir = env))
    ## 优先带 Plot / Kmeans 的对象
    prefer <- objs[grepl("Plot|Kmeans|order", objs, ignore.case = TRUE)]
    if (length(prefer) == 1L) return(get(prefer[[1]], envir = env))
    as.list(env)
}

parse_row_meta <- function(rn) {
    has_mf <- grepl("^[MF]_", rn)
    sex <- ifelse(has_mf & grepl("^M_", rn), "Male",
           ifelse(has_mf & grepl("^F_", rn), "Female", NA_character_))
    rest <- ifelse(has_mf, sub("^[MF]_", "", rn), rn)
    ct <- sapply(strsplit(rest, "__", fixed = TRUE), function(x) x[[1]])
    list(celltype = unname(ct), sex = unname(sex))
}

mark_label <- function(r) {
    if (grepl("^[MF]_", r)) {
        sx <- ifelse(grepl("^M_", r), "M", "F")
        rest <- sub("^[MF]_", "", r)
        parts <- strsplit(rest, "__", fixed = TRUE)[[1]]
        if (length(parts) < 2) return(gsub("__", "_", r, fixed = TRUE))
        return(paste0(parts[1], "_", sx, "_", parts[2]))
    }
    gsub("__", "_", r, fixed = TRUE)
}

align_mat_km <- function(mat, km) {
    if (!all(c("genes", "cluster") %in% names(km))) {
        stop("kmeans 表缺 genes/cluster: ", paste(names(km), collapse = ","))
    }
    if (!identical(rownames(mat), km$genes)) {
        m <- match(km$genes, rownames(mat))
        if (any(is.na(m))) {
            m2 <- match(rownames(mat), km$genes)
            if (!any(is.na(m2))) {
                km <- km[m2, , drop = FALSE]
            } else {
                stop("矩阵行名与 kmeans$genes 对不齐")
            }
        } else {
            mat <- mat[m, , drop = FALSE]
        }
    }
    list(mat = mat, km = km)
}

Plot_dot <- function(Plot, Kmeans_order) {
    d <- melt(Plot)
    d$class <- Kmeans_order$cluster[match(d$Var1, Kmeans_order$genes)]
    ages <- unique(as.character(d$Var2))
    if (all(grepl("^[0-9.]+$", ages))) {
        age_levels <- as.character(sort(as.numeric(ages)))
        d$Var2 <- factor(d$Var2, levels = age_levels)
        d$Var2_1 <- as.numeric(d$Var2)
    } else {
        d$Var2_1 <- as.numeric(factor(d$Var2, levels = colnames(Plot)))
    }
    d$class <- factor(d$class, levels = sort(unique(Kmeans_order$cluster)))
    d
}

make_bar_data <- function(kc, Down, Middle, Up) {
    meta <- parse_row_meta(kc$genes)
    kc$cell_type <- meta$celltype
    index <- names(table(kc$cell_type))
    cell_type <- character(0)
    num_DEGs <- integer(0)
    direction <- character(0)
    for (ct in index) {
        cell_type <- c(cell_type, rep(ct, 3))
        num_DEGs <- c(
            num_DEGs,
            sum(kc$cell_type == ct & kc$cluster %in% Down),
            sum(kc$cell_type == ct & kc$cluster %in% Up),
            sum(kc$cell_type == ct & kc$cluster %in% Middle)
        )
        direction <- c(direction, c("Young", "Old", "Middle"))
    }
    data <- data.frame(cell_type = cell_type, num_DEGs = num_DEGs,
                       direction = direction, stringsAsFactors = FALSE)
    data$direction <- factor(data$direction, levels = c("Young", "Middle", "Old"))
    tot <- tapply(data$num_DEGs, data$cell_type, sum)
    data$cell_type <- factor(data$cell_type, levels = rev(names(tot)[order(tot)]))
    data
}

plot_bar <- function(data, title, ymax) {
    ggplot(data, aes(x = cell_type, y = num_DEGs, fill = direction)) +
        geom_col(position = "stack", width = 0.5) +
        scale_fill_manual(values = c(Young = "lightblue", Old = "pink",
                                     Middle = "darkgreen")) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, ymax * 1.05)) +
        labs(x = "", y = "Number of DEGs", title = title) +
        theme_classic(base_size = 14) +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            legend.position = "right"
        )
}

run_one <- function(species, job) {
    message("\n========== ", species, " ==========")
    dir.create(job$out_dir, showWarnings = FALSE, recursive = TRUE)
    new_order <- c(job$new_order_down, job$new_order_middle, job$new_order_up)

    plot_fp <- find_existing(job$plot_raw)
    km_fp <- find_existing(job$km_raw)
    already_ordered <- FALSE
    if (is.na(plot_fp) || is.na(km_fp)) {
        message("  无原始 kmeans，改用已重排文件")
        plot_fp <- find_existing(job$plot_order)
        km_fp <- find_existing(job$km_order)
        already_ordered <- TRUE
    }
    if (is.na(plot_fp) || is.na(km_fp)) {
        warning(species, ": 找不到 kmeans / Plot 文件，跳过")
        return(invisible(NULL))
    }
    message("  [Plot] ", plot_fp)
    message("  [Kmeans] ", km_fp)
    mat0 <- load_r_object(plot_fp)
    km0 <- load_r_object(km_fp)
    if (is.list(mat0) && !is.matrix(mat0) && !is.data.frame(mat0)) {
        hit <- names(mat0)[grepl("Plot", names(mat0), ignore.case = TRUE)]
        if (length(hit)) mat0 <- mat0[[hit[[1]]]]
    }
    if (is.list(km0) && !is.data.frame(km0)) {
        hit <- names(km0)[grepl("Kmeans|order", names(km0), ignore.case = TRUE)]
        if (length(hit)) km0 <- km0[[hit[[1]]]]
    }
    aligned <- align_mat_km(as.matrix(mat0), as.data.frame(km0, stringsAsFactors = FALSE))
    mat0 <- aligned$mat
    km0 <- aligned$km
    message(sprintf("  矩阵 %d x %d | cluster=%s",
                    nrow(mat0), ncol(mat0),
                    paste(sort(unique(km0$cluster)), collapse = ",")))

    if (!already_ordered) {
        keep <- !(km0$cluster %in% job$rm_clusters) & km0$cluster %in% new_order
        km <- km0[keep, , drop = FALSE]
        mat <- mat0[km$genes, , drop = FALSE]
        km$cluster <- match(km$cluster, new_order)
        ord <- order(km$cluster, km$genes)
        km <- km[ord, , drop = FALSE]
        mat <- mat[km$genes, , drop = FALSE]
        rownames(km) <- NULL
        message("  已按 new_order 重编号: down=", paste(job$new_order_down, collapse = ","),
                " | middle=", paste(job$new_order_middle, collapse = ","),
                " | up=", paste(job$new_order_up, collapse = ","))
    } else {
        km <- km0
        mat <- mat0
        message("  使用已重排文件，不再 match(new_order)")
    }
    message("  cluster 计数:")
    print(table(km$cluster))

    saveRDS(mat, file.path(plot_dir, paste0(job$stem, "_DEGs_Plot_neworder_2026Aug.rds")))
    saveRDS(km,  file.path(plot_dir, paste0(job$stem, "_DEGs_Plot_Kmeans_order_2026Aug.rds")))

    ## heatmap
    meta <- parse_row_meta(rownames(mat))
    if (isTRUE(job$split_sex) && any(!is.na(meta$sex))) {
        ra <- rowAnnotation(
            CellType = meta$celltype, Sex = meta$sex,
            col = list(CellType = celltype_colors, Sex = sex_colors),
            annotation_name_side = "top"
        )
    } else {
        ra <- rowAnnotation(
            CellType = meta$celltype,
            col = list(CellType = celltype_colors),
            annotation_name_side = "top"
        )
    }
    row_sp <- factor(km$cluster, levels = sort(unique(km$cluster)))
    mark_use <- job$marks[job$marks %in% rownames(mat)]
    at <- match(mark_use, rownames(mat))
    mark_labs <- vapply(mark_use, mark_label, character(1))
    ht <- Heatmap(
        mat, name = "z", right_annotation = ra, border = TRUE,
        use_raster = TRUE, raster_quality = 5,
        show_row_names = FALSE, show_column_names = TRUE,
        rect_gp = gpar(col = "white", lwd = 0),
        cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun,
        row_split = row_sp,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 12),
        width = unit(if (!is.null(job$heatmap_col_mm)) job$heatmap_col_mm else 5,
                     "mm") * ncol(mat),
        heatmap_legend_param = list(title = "z", border = "black",
                                    at = c(-2, -1, 0, 1, 2))
    )
    if (length(at)) {
        ht <- ht + rowAnnotation(
            link = anno_mark(at = at, labels = mark_labs,
                             labels_gp = gpar(fontsize = 12, fontface = "italic"))
        )
    }
    ht_file <- file.path(plot_dir, paste0(job$stem, "_DEGs_neworder_2026Aug.png"))
    png(ht_file, height = 7000,
        width = if (!is.null(job$ht_png_width)) job$ht_png_width else 4000,
        res = 72 * 12)
    ## 左边距加大，避免 row_split 编号（尤其是 7）被裁成一半
    draw(ht, padding = unit(c(2, 12, 2, 2), "mm"))
    dev.off()
    message("  saved: ", ht_file)

    ## trend
    dot <- Plot_dot(mat, km)
    p_trend <- ggplot(dot, aes(x = Var2_1, y = value)) +
        geom_point(size = 0, position = position_jitter(width = 0.2, height = 0.2),
                   alpha = 0.01) +
        stat_summary(fun = mean, geom = "line", color = "red", linewidth = 1) +
        facet_wrap(~ class, ncol = 1, strip.position = "right") +
        scale_x_discrete(breaks = NULL) +
        scale_y_continuous(breaks = c(-2, 0, 2)) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            strip.text = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank()
        ) +
        ylab("")
    tr_file <- file.path(plot_dir, paste0(job$stem, "_trend_2026Aug.png"))
    ggsave(tr_file, p_trend, height = 8, width = 1.2)
    message("  saved: ", tr_file)

    ## bar
    message(sprintf("  barplot 分组: Young=%s | Middle=%s | Old=%s",
                    paste(job$Down, collapse = ","),
                    paste(job$Middle, collapse = ","),
                    paste(job$Up, collapse = ",")))
    if (isTRUE(job$split_sex)) {
        sex_vec <- parse_row_meta(km$genes)$sex
        kc_M <- km[!is.na(sex_vec) & sex_vec == "Male", , drop = FALSE]
        kc_F <- km[!is.na(sex_vec) & sex_vec == "Female", , drop = FALSE]
        data_M <- make_bar_data(kc_M, job$Down, job$Middle, job$Up)
        data_F <- make_bar_data(kc_F, job$Down, job$Middle, job$Up)
        ymax <- if (is.finite(job$bar_ymax)) job$bar_ymax else
            max(tapply(data_M$num_DEGs, data_M$cell_type, sum),
                tapply(data_F$num_DEGs, data_F$cell_type, sum))
        ggsave(file.path(plot_dir, paste0(job$stem, "_data_barplot_male_2026Aug.png")),
               plot_bar(data_M, "Male", ymax), width = 7, height = 3)
        ggsave(file.path(plot_dir, paste0(job$stem, "_data_barplot_female_2026Aug.png")),
               plot_bar(data_F, "Female", ymax), width = 7, height = 3)
        write.table(data_M, file.path(plot_dir, paste0(job$stem, "_data_barplot_male_2026Aug.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        write.table(data_F, file.path(plot_dir, paste0(job$stem, "_data_barplot_female_2026Aug.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        message("  saved: male/female barplots")
    } else {
        data <- make_bar_data(km, job$Down, job$Middle, job$Up)
        ymax <- if (is.finite(job$bar_ymax)) job$bar_ymax else
            max(tapply(data$num_DEGs, data$cell_type, sum))
        ggsave(file.path(plot_dir, paste0(job$stem, "_data_barplot_2026Aug.png")),
               plot_bar(data, species, ymax), width = 7, height = 3)
        write.table(data, file.path(plot_dir, paste0(job$stem, "_data_barplot_2026Aug.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        message("  saved: ", paste0(job$stem, "_data_barplot_2026Aug.png"))
    }
    list(mat = mat, km = km)
}


########################################################################
# 跨物种 Young/Old overlap（基因全部转到人源）
#   表：所有 cell type（来自 Sheet A / Gene_clusters 的 km 行）
#   图：仍只画 Figure2 F/G/H = MG / Rod / RGC
########################################################################
ortho_candidates <- c(
    file.path(fish_dir, "HMZ_ortholog_combined"),
    "/projects/hmz-aging/Zebrafish_New_Datasets/HMZ_ortholog_combined",
    file.path(mouse_dir, "HMZ_ortholog_combined"),
    "/projects/hmz-aging/HMZ_ortholog_combined",
    "/zp1/data/share/Human_aging_new/HMZ_ortholog_combined_2025"
)
young_fill <- "lightblue"
old_fill   <- "pink"
point_col  <- "#026AB1"
fgh_cts <- c("MG", "Rod", "RGC")   # 只出 Figure2 panel
ct_order_hmz <- c("MG", "RGC", "AC", "HC", "Cone", "Rod", "BC",
                  "Microglia", "RPE")

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

norm_ortho <- function(ortho) {
    out <- ortho
    for (nm in names(out)) {
        if (is.data.frame(out[[nm]])) names(out[[nm]]) <- tolower(names(out[[nm]]))
    }
    out
}

native_symbol <- function(genes) sub(".*__", "", as.character(genes))

deg_symbols <- function(km, ct, clusters) {
    meta <- parse_row_meta(km$genes)
    keep <- !is.na(meta$celltype) & meta$celltype == ct & km$cluster %in% clusters
    unique(native_symbol(km$genes[keep]))
}

map_to_human <- function(genes, species, ortho) {
    ## between-species：全部转到人源 symbol，映不上的丢掉
    genes <- unique(as.character(genes))
    genes <- genes[nzchar(genes) & !is.na(genes)]
    if (!length(genes)) return(character(0))
    if (species == "Human") return(unique(toupper(genes)))
    if (species == "Mouse") {
        hm <- ortho$HM
        if (is.null(hm) || !all(c("human", "mouse") %in% names(hm))) return(character(0))
        mapped <- unique(as.character(hm$human[tolower(hm$mouse) %in% tolower(genes)]))
    } else {
        hz <- ortho$HZ
        if (is.null(hz) || !all(c("human", "zebrafish") %in% names(hz))) return(character(0))
        mapped <- unique(as.character(hz$human[tolower(hz$zebrafish) %in% tolower(genes)]))
    }
    mapped <- mapped[nzchar(mapped) & !is.na(mapped)]
    unique(toupper(mapped))
}

load_deg_pair <- function(species, deg_res) {
    if (is.list(deg_res) && is.list(deg_res[[species]]) &&
        all(c("mat", "km") %in% names(deg_res[[species]]))) {
        return(deg_res[[species]])
    }
    job <- jobs[[species]]
    plot_fp <- file.path(plot_dir, paste0(job$stem, "_DEGs_Plot_neworder_2026Aug.rds"))
    km_fp   <- file.path(plot_dir, paste0(job$stem, "_DEGs_Plot_Kmeans_order_2026Aug.rds"))
    if (!file.exists(plot_fp) || !file.exists(km_fp)) return(NULL)
    list(mat = load_r_object(plot_fp), km = load_r_object(km_fp))
}

## Sheet A 同源：从三物种 kmeans gene 表取 celltype 并集
cts_from_sheetA <- function(pairs) {
    cts <- unique(unlist(lapply(pairs, function(p) {
        if (is.null(p) || is.null(p$km)) return(character(0))
        meta <- parse_row_meta(p$km$genes)
        as.character(meta$celltype)
    })))
    cts <- unique(cts[!is.na(cts) & nzchar(cts)])
    c(intersect(ct_order_hmz, cts), setdiff(sort(cts), ct_order_hmz))
}

## 一对方向的 set → pairwise + triple gene 行 + 计数行
overlap_rows_one_dir <- function(ct, dirn, sets) {
    sets <- lapply(sets, function(x) unique(x[nzchar(x) & !is.na(x)]))
    Z <- sets$Z; M <- sets$M; H <- sets$H
    hm <- intersect(H, M)
    hz <- intersect(H, Z)
    mz <- intersect(M, Z)
    trip <- Reduce(intersect, list(H, M, Z))
    gene_parts <- list(
        if (length(hm)) data.frame(celltype = ct, direction = dirn, set = "HM",
                                   gene = sort(hm), stringsAsFactors = FALSE),
        if (length(hz)) data.frame(celltype = ct, direction = dirn, set = "HZ",
                                   gene = sort(hz), stringsAsFactors = FALSE),
        if (length(mz)) data.frame(celltype = ct, direction = dirn, set = "MZ",
                                   gene = sort(mz), stringsAsFactors = FALSE),
        if (length(trip)) data.frame(celltype = ct, direction = dirn, set = "triple",
                                     gene = sort(trip), stringsAsFactors = FALSE)
    )
    gene_parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, gene_parts)
    genes_df <- if (length(gene_parts)) do.call(rbind, gene_parts) else NULL
    stats_df <- data.frame(
        celltype = ct, direction = dirn,
        n_Z = length(Z), n_M = length(M), n_H = length(H),
        n_HM = length(hm), n_HZ = length(hz), n_MZ = length(mz),
        n_triple = length(trip),
        stringsAsFactors = FALSE
    )
    list(genes = genes_df, stats = stats_df, trip = trip,
         sets = list(Z = Z, M = M, H = H))
}

make_upset <- function(sets, fill) {
    sets <- lapply(sets, function(x) unique(x[nzchar(x) & !is.na(x)]))
    sets <- sets[lengths(sets) > 0]
    if (length(sets) < 2) return(NULL)
    m <- make_comb_mat(sets, mode = "distinct")
    set_ord <- intersect(c("Z", "M", "H"), set_name(m))
    cs <- comb_size(m)
    ss <- set_size(m)
    bar_gp <- gpar(fill = fill, col = "black", lwd = 0.7)
    name_gp <- gpar(fontsize = 10)
    num_gp  <- gpar(fontsize = 9)
    ax_gp   <- gpar(fontsize = 8)
    top_name <- paste0("Intersection", intToUtf8(160L), "size")
    ## 只加高顶部柱；点阵 width/height 不变
    top_h <- unit(6.6, "cm")
    right_w <- unit(2.6, "cm")
    top_bar <- anno_barplot(
        cs, gp = bar_gp, ylim = c(0, max(cs) * 1.22),
        add_numbers = TRUE, numbers_gp = num_gp,
        numbers_rot = 0, numbers_offset = unit(1.2, "mm"),
        axis_param = list(side = "left", gp = ax_gp, labels_rot = 0)
    )
    top_args <- list(
        top_bar,
        annotation_name_side = "left",
        annotation_name_rot = 90,
        annotation_name_gp = name_gp,
        annotation_name_offset = unit(10, "mm"),
        annotation_label = top_name,
        height = top_h
    )
    names(top_args)[1] <- top_name
    top_anno <- tryCatch(
        do.call(HeatmapAnnotation, top_args),
        error = function(e) {
            top_args$annotation_label <- NULL
            do.call(HeatmapAnnotation, top_args)
        }
    )
    right_anno <- tryCatch(
        upset_right_annotation(
            m, gp = bar_gp,
            add_numbers = TRUE, numbers_gp = num_gp,
            numbers_rot = 0, numbers_offset = unit(1.4, "mm"),
            ylim = c(0, max(ss) * 1.32),
            annotation_name_gp = name_gp,
            annotation_name_side = "bottom",
            annotation_name_rot = 0,
            annotation_name_offset = unit(6.5, "mm"),
            axis_param = list(side = "bottom", gp = ax_gp, labels_rot = 0),
            width = right_w
        ),
        error = function(e) {
            rowAnnotation(
                `Set size` = anno_barplot(
                    ss, gp = bar_gp, ylim = c(0, max(ss) * 1.32),
                    add_numbers = TRUE, numbers_gp = num_gp,
                    numbers_rot = 0, numbers_offset = unit(1.4, "mm"),
                    axis_param = list(side = "bottom", gp = ax_gp, labels_rot = 0)
                ),
                annotation_name_gp = name_gp,
                annotation_name_side = "bottom",
                annotation_name_rot = 0,
                annotation_name_offset = unit(6.5, "mm"),
                width = right_w
            )
        }
    )
    ht <- tryCatch(
        UpSet(
            m, set_order = set_ord, comb_order = order(-cs),
            pt_size = unit(1.8, "mm"), lwd = 1.1, comb_col = "black",
            bg_col = "#E6E6E6", bg_pt_col = "#CCCCCC",
            top_annotation = top_anno, right_annotation = right_anno,
            width = unit(5.0, "cm"),
            height = unit(13, "mm"),
            row_gap = unit(0.6, "mm"),
            show_heatmap_legend = FALSE,
            row_names_gp = gpar(fontsize = 13)
        ),
        error = function(e) {
            message("    UpSet 简化重试: ", conditionMessage(e))
            UpSet(
                m, set_order = set_ord, comb_order = order(-cs),
                pt_size = unit(2.4, "mm"), lwd = 1.2, comb_col = "black",
                top_annotation = top_anno, right_annotation = right_anno,
                row_names_gp = gpar(fontsize = 11)
            )
        }
    )
    ht
}

find_mat_rows <- function(mat, ct, gene, species) {
    rn <- rownames(mat)
    gene <- as.character(gene)[1]
    if (species == "Human") {
        rest <- sub("^[MF]_", "", rn)
        hit <- grepl(sprintf("^[MF]_%s__%s$", ct, gene), rn, ignore.case = TRUE) |
            (tolower(rest) == tolower(paste0(ct, "__", gene)))
    } else {
        hit <- tolower(rn) == tolower(paste0(ct, "__", gene))
    }
    if (!any(hit)) {
        meta <- parse_row_meta(rn)
        hit <- meta$celltype == ct & tolower(native_symbol(rn)) == tolower(gene)
    }
    rn[hit]
}

human_to_native <- function(hgene, species, ortho) {
    hgene <- toupper(as.character(hgene)[1])
    if (species == "Human") return(hgene)
    df <- if (species == "Mouse") ortho$HM else ortho$HZ
    sp_col <- if (species == "Mouse") "mouse" else "zebrafish"
    if (is.null(df) || !all(c("human", sp_col) %in% names(df))) return(character(0))
    hits <- as.character(df[[sp_col]][toupper(as.character(df$human)) == hgene])
    unique(hits[nzchar(hits) & !is.na(hits)])
}

first_in_mat <- function(mat, ct, genes, species) {
    for (g in genes) {
        if (length(find_mat_rows(mat, ct, g, species))) return(as.character(g)[1])
    }
    NA_character_
}

gene_age_cor <- function(mat, ct, gene, species) {
    rows <- find_mat_rows(mat, ct, gene, species)
    if (!length(rows)) return(NA_real_)
    ages <- suppressWarnings(as.numeric(colnames(mat)))
    if (all(is.na(ages))) ages <- seq_len(ncol(mat))
    age_vec <- rep(ages, each = length(rows))
    z_vec <- as.numeric(mat[rows, , drop = FALSE])
    ok <- is.finite(age_vec) & is.finite(z_vec)
    if (sum(ok) < 3L) return(NA_real_)
    suppressWarnings(cor(age_vec[ok], z_vec[ok], method = "pearson"))
}

pick_overlap_example <- function(trip_y, trip_o, ct, pairs, ortho) {
    cands <- unique(c(trip_o, trip_y))
    cands <- cands[nzchar(cands) & !is.na(cands)]
    if (!length(cands)) return(NULL)
    scored <- lapply(cands, function(hg) {
        nat_z <- first_in_mat(pairs$Zebrafish$mat, ct,
                              human_to_native(hg, "Zebrafish", ortho), "Zebrafish")
        nat_m <- first_in_mat(pairs$Mouse$mat, ct,
                              human_to_native(hg, "Mouse", ortho), "Mouse")
        nat_h <- first_in_mat(pairs$Human$mat, ct, hg, "Human")
        if (anyNA(c(nat_z, nat_m, nat_h))) return(NULL)
        cz <- gene_age_cor(pairs$Zebrafish$mat, ct, nat_z, "Zebrafish")
        cm <- gene_age_cor(pairs$Mouse$mat, ct, nat_m, "Mouse")
        ch <- gene_age_cor(pairs$Human$mat, ct, nat_h, "Human")
        if (anyNA(c(cz, cm, ch))) return(NULL)
        list(human = hg, Zebrafish = nat_z, Mouse = nat_m, Human = nat_h,
             cor_z = cz, cor_m = cm, cor_h = ch,
             cor_mean = mean(c(cz, cm, ch)))
    })
    scored <- Filter(Negate(is.null), scored)
    if (!length(scored)) return(NULL)
    means <- vapply(scored, `[[`, numeric(1), "cor_mean")
    scored[[which.max(means)]]
}

plot_gene_trend <- function(mat, ct, gene, species) {
    ## 与热图同一矩阵、同一数值：不重算 z；点色三物种相同
    empty <- ggplot() + theme_void()
    if (is.null(gene) || is.na(gene) || !nzchar(gene)) return(empty)
    rows <- find_mat_rows(mat, ct, gene, species)
    ttl <- paste(species, gene)
    if (!length(rows)) return(empty)
    ages <- colnames(mat)
    age_num <- suppressWarnings(as.numeric(ages))
    if (all(is.na(age_num))) age_num <- seq_along(ages)
    df <- do.call(rbind, lapply(seq_along(rows), function(i) {
        data.frame(age = age_num, z = as.numeric(mat[rows[i], ]),
                   stringsAsFactors = FALSE)
    }))
    df <- df[is.finite(df$age) & is.finite(df$z), , drop = FALSE]
    if (!nrow(df)) return(empty)
    ggplot(df, aes(x = age, y = z)) +
        geom_point(color = point_col, size = 1.7) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red",
                    fill = "grey75", linewidth = 0.7, alpha = 0.35) +
        scale_x_continuous(breaks = function(x) pretty(x, n = 4)) +
        scale_y_continuous(breaks = c(-2, 0, 2)) +
        labs(title = ttl,
             x = switch(species,
                        Zebrafish = "Age (Months)",
                        Mouse = "Age (Weeks)",
                        Human = "Age (Years)",
                        "Age"),
             y = "RNA z-score") +
        theme_classic(base_size = 12) +
        theme(
            plot.title = element_text(size = 12, hjust = 0.5, face = "italic",
                                      margin = margin(b = 2)),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.title.x = element_text(size = 11, margin = margin(t = 3)),
            axis.title.y = element_text(size = 11, margin = margin(r = 2)),
            axis.text = element_text(size = 10, color = "black"),
            legend.position = "none",
            aspect.ratio = 1,
            plot.margin = margin(2, 4, 2, 4)
        )
}

grab_on_page <- function(FUN, width, height) {
    fn <- tempfile(fileext = ".pdf")
    grDevices::pdf(fn, width = width, height = height)
    on.exit({
        if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
        unlink(fn)
    })
    grid.newpage()
    FUN()
    grid.grab(wrap = TRUE, wrap.grobs = TRUE)
}

draw_ct_panel <- function(g_y, g_o, g_z, g_m, g_h, newpage = TRUE) {
    if (newpage) grid.newpage()
    pushViewport(viewport(layout = grid.layout(
        2, 1, heights = unit(c(4.05, 2.45), "in")
    )))
    pushViewport(viewport(
        layout.pos.row = 1,
        layout = grid.layout(1, 2, widths = unit(c(0.5, 0.5), "npc"))
    ))
    pushViewport(viewport(layout.pos.col = 1))
    grid.draw(g_y)
    popViewport()
    pushViewport(viewport(layout.pos.col = 2))
    grid.draw(g_o)
    popViewport()
    popViewport()
    pushViewport(viewport(layout.pos.row = 2, layout = grid.layout(1, 3)))
    for (i in seq_along(list(g_z, g_m, g_h))) {
        g <- list(g_z, g_m, g_h)[[i]]
        pushViewport(viewport(layout.pos.col = i))
        pushViewport(viewport(width = unit(1, "snpc"), height = unit(1, "snpc")))
        grid.draw(g)
        popViewport()
        popViewport()
    }
    popViewport()
    popViewport()
}

plot_fgh_overlap <- function(deg_res = NULL) {
    message("\n========== 跨物种 Young/Old overlap（全部 cell type；图仅 F/G/H） ==========")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ortho_fp <- find_existing(ortho_candidates)
    if (is.na(ortho_fp)) stop("找不到 HMZ_ortholog_combined")
    message("  ortholog: ", ortho_fp)
    ortho <- norm_ortho(load_ortholog(ortho_fp))
    pairs <- lapply(c("Zebrafish", "Mouse", "Human"), function(sp) load_deg_pair(sp, deg_res))
    names(pairs) <- c("Zebrafish", "Mouse", "Human")
    miss <- names(pairs)[vapply(pairs, is.null, logical(1))]
    if (length(miss)) stop("缺 DEG 矩阵: ", paste(miss, collapse = ", "))

    ## cell type 并集 = Sheet A (Gene_clusters) 同源来源
    all_cts <- cts_from_sheetA(pairs)
    message("  SheetA celltypes (n=", length(all_cts), "): ",
            paste(all_cts, collapse = ", "))
    message("  UpSet 出图仅: ", paste(fgh_cts, collapse = ", "))

    one_ct <- function(ct, do_plot = FALSE) {
        message("  ---- ", ct, if (do_plot) " [plot]" else " [table]", " ----")
        mapped <- list()
        for (dirn in c("Young", "Old")) {
            mapped[[dirn]] <- list()
            for (sp in c("Zebrafish", "Mouse", "Human")) {
                cl <- if (dirn == "Young") jobs[[sp]]$Down else jobs[[sp]]$Up
                g <- deg_symbols(pairs[[sp]]$km, ct, cl)
                mapped[[dirn]][[sp]] <- map_to_human(g, sp, ortho)
                message(sprintf("    %s %s: native=%d  human=%d",
                                sp, dirn, length(g), length(mapped[[dirn]][[sp]])))
            }
        }
        sets_y <- list(Z = mapped$Young$Zebrafish, M = mapped$Young$Mouse,
                       H = mapped$Young$Human)
        sets_o <- list(Z = mapped$Old$Zebrafish,   M = mapped$Old$Mouse,
                       H = mapped$Old$Human)
        res_y <- overlap_rows_one_dir(ct, "Young", sets_y)
        res_o <- overlap_rows_one_dir(ct, "Old", sets_o)
        message(sprintf("    Young: n_H=%d n_M=%d n_Z=%d | HM=%d HZ=%d MZ=%d triple=%d",
                        res_y$stats$n_H, res_y$stats$n_M, res_y$stats$n_Z,
                        res_y$stats$n_HM, res_y$stats$n_HZ, res_y$stats$n_MZ,
                        res_y$stats$n_triple))
        message(sprintf("    Old:   n_H=%d n_M=%d n_Z=%d | HM=%d HZ=%d MZ=%d triple=%d",
                        res_o$stats$n_H, res_o$stats$n_M, res_o$stats$n_Z,
                        res_o$stats$n_HM, res_o$stats$n_HZ, res_o$stats$n_MZ,
                        res_o$stats$n_triple))
        if (length(res_y$trip))
            message("      Young triple: ", paste(sort(res_y$trip), collapse = ", "))
        if (length(res_o$trip))
            message("      Old triple: ", paste(sort(res_o$trip), collapse = ", "))

        gene_parts <- Filter(Negate(is.null), list(res_y$genes, res_o$genes))
        ov <- if (length(gene_parts)) do.call(rbind, gene_parts) else NULL
        st <- rbind(res_y$stats, res_o$stats)

        if (isTRUE(do_plot)) {
            ht_y <- make_upset(sets_y, young_fill)
            ht_o <- make_upset(sets_o, old_fill)
            eg <- pick_overlap_example(res_y$trip, res_o$trip, ct, pairs, ortho)
            if (is.null(eg)) {
                message("    example gene: no Z∩M∩H gene found in all three heatmaps")
                pZ <- plot_gene_trend(NULL, ct, NA_character_, "Zebrafish")
                pM <- plot_gene_trend(NULL, ct, NA_character_, "Mouse")
                pH <- plot_gene_trend(NULL, ct, NA_character_, "Human")
            } else {
                message(sprintf(
                    "    example gene: %s  (Z=%s  M=%s  H=%s)  mean cor=%.3f (Z=%.3f M=%.3f H=%.3f)",
                    eg$human, eg$Zebrafish, eg$Mouse, eg$Human,
                    eg$cor_mean, eg$cor_z, eg$cor_m, eg$cor_h))
                pZ <- plot_gene_trend(pairs$Zebrafish$mat, ct, eg$Zebrafish, "Zebrafish")
                pM <- plot_gene_trend(pairs$Mouse$mat, ct, eg$Mouse, "Mouse")
                pH <- plot_gene_trend(pairs$Human$mat, ct, eg$Human, "Human")
            }
            panel_tag <- switch(ct, MG = "F", Rod = "G", RGC = "H", "X")
            fp <- file.path(plot_dir,
                            sprintf("Figure2%s_%s_between_species_2026Aug.png",
                                    panel_tag, ct))
            g_y <- if (is.null(ht_y)) nullGrob() else grab_on_page(function() {
                draw(ht_y, newpage = TRUE, padding = unit(c(8, 8, 2, 2), "mm"))
            }, 4.35, 4.05)
            g_o <- if (is.null(ht_o)) nullGrob() else grab_on_page(function() {
                draw(ht_o, newpage = TRUE, padding = unit(c(8, 8, 2, 2), "mm"))
            }, 4.35, 4.05)
            g_z <- grab_on_page(function() print(pZ, newpage = TRUE), 2.05, 2.05)
            g_m <- grab_on_page(function() print(pM, newpage = TRUE), 2.05, 2.05)
            g_h <- grab_on_page(function() print(pH, newpage = TRUE), 2.05, 2.05)
            png(fp, width = 8.75, height = 6.65, units = "in", res = 300)
            tryCatch(draw_ct_panel(g_y, g_o, g_z, g_m, g_h, newpage = TRUE),
                     finally = {
                         if (dev.cur() != 1L) dev.off()
                     })
            message("  saved: ", fp)
        }
        list(overlap = ov, stats = st)
    }

    ct_out <- lapply(all_cts, function(ct) {
        tryCatch(one_ct(ct, do_plot = ct %in% fgh_cts), error = function(e) {
            message("  ERROR ", ct, ": ", conditionMessage(e))
            NULL
        })
    })
    names(ct_out) <- all_cts
    ok <- !vapply(ct_out, is.null, logical(1))
    if (any(!ok)) warning("跨物种 overlap 部分细胞类型失败: ",
                          paste(all_cts[!ok], collapse = ", "))
    overlap_rows <- lapply(ct_out[ok], `[[`, "overlap")
    stat_rows <- lapply(ct_out[ok], `[[`, "stats")
    overlap_rows <- Filter(function(x) !is.null(x) && nrow(x) > 0, overlap_rows)
    if (length(overlap_rows)) {
        ov_all <- do.call(rbind, overlap_rows)
        rownames(ov_all) <- NULL
        write.table(ov_all,
                    file.path(plot_dir, "HMZ_crossSp_overlap_genes_2026Aug.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        ## 兼容旧文件名：仅 triple
        trip_only <- ov_all[ov_all$set == "triple", c("celltype", "direction", "gene"),
                            drop = FALSE]
        if (nrow(trip_only)) {
            write.table(trip_only,
                        file.path(plot_dir, "HMZ_triple_overlap_genes_2026Aug.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
        }
    } else {
        ov_all <- NULL
    }
    if (length(stat_rows)) {
        st_all <- do.call(rbind, stat_rows)
        rownames(st_all) <- NULL
        write.table(st_all,
                    file.path(plot_dir, "HMZ_overlap_set_sizes_2026Aug.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        st_all <- NULL
    }
    invisible(list(triple = ov_all, sizes = st_all))
}


########################################################################
# Human male vs female：Young/Old 一致 (Common) 与不一致 (Divergent) 热图
# 行 = celltype__gene；左 Female / 右 Male；值 = 热图同一套 z-score
########################################################################
plot_human_sex_overlap <- function(mat, km, job) {
    message("\n========== Human male/female Young-Old overlap ==========")
    Down <- job$Down
    Up <- job$Up
    meta <- parse_row_meta(km$genes)
    tab <- data.frame(
        row = km$genes,
        key = sub("^[MF]_", "", km$genes),
        sex = meta$sex,
        ct  = meta$celltype,
        dir = ifelse(km$cluster %in% Down, "Young",
              ifelse(km$cluster %in% Up, "Old", "Middle")),
        stringsAsFactors = FALSE
    )
    tab <- tab[!is.na(tab$sex) & tab$dir %in% c("Young", "Old"), , drop = FALSE]
    f_y <- unique(tab$key[tab$sex == "Female" & tab$dir == "Young"])
    f_o <- unique(tab$key[tab$sex == "Female" & tab$dir == "Old"])
    m_y <- unique(tab$key[tab$sex == "Male"   & tab$dir == "Young"])
    m_o <- unique(tab$key[tab$sex == "Male"   & tab$dir == "Old"])
    groups <- list(
        `Young`            = intersect(f_y, m_y),
        `Old`              = intersect(f_o, m_o),
        `F:Old M:Young`    = intersect(f_o, m_y),
        `F:Young M:Old`    = intersect(f_y, m_o)
    )
    for (nm in names(groups)) {
        message(sprintf("  %s: %d", nm, length(groups[[nm]])))
    }
    gene_tab <- do.call(rbind, lapply(names(groups), function(nm) {
        k <- groups[[nm]]
        if (!length(k)) return(NULL)
        data.frame(class = nm, key = k,
                   celltype = sub("__.*$", "", k),
                   gene = sub("^[^_]+__", "", k),
                   stringsAsFactors = FALSE)
    }))
    if (!is.null(gene_tab) && nrow(gene_tab)) {
        fp <- file.path(plot_dir, "Human_sex_overlap_genes_2026Aug.tsv")
        write.table(gene_tab, fp, sep = "\t", quote = FALSE, row.names = FALSE)
        message("  tsv: ", fp)
    }

    pick_row <- function(keys, sex) {
        hit <- match(paste(sex, keys), paste(tab$sex, tab$key))
        tab$row[hit]
    }
    draw_pair_ht <- function(keys_a, lab_a, keys_b, lab_b, outfile, panel) {
        keys <- c(keys_a, keys_b)
        if (!length(keys)) {
            warning(panel, ": 无基因，跳过")
            return(invisible(NULL))
        }
        f_rn <- pick_row(keys, "Female")
        m_rn <- pick_row(keys, "Male")
        ok <- !is.na(f_rn) & !is.na(m_rn) & f_rn %in% rownames(mat) & m_rn %in% rownames(mat)
        keys <- keys[ok]; f_rn <- f_rn[ok]; m_rn <- m_rn[ok]
        keys_a <- keys_a[keys_a %in% keys]
        keys_b <- keys_b[keys_b %in% keys]
        keys <- c(keys_a, keys_b)
        f_rn <- pick_row(keys, "Female")
        m_rn <- pick_row(keys, "Male")
        ages <- colnames(mat)
        mat_F <- mat[f_rn, , drop = FALSE]
        mat_M <- mat[m_rn, , drop = FALSE]
        rownames(mat_F) <- keys
        rownames(mat_M) <- keys
        colnames(mat_F) <- paste0("F__", ages)
        colnames(mat_M) <- paste0("M__", ages)
        mat_both <- cbind(mat_F, mat_M)
        col_sp <- factor(c(rep("Female", length(ages)), rep("Male", length(ages))),
                         levels = c("Female", "Male"))
        row_sp <- factor(
            c(rep(sprintf("%s (%d genes)", lab_a, length(keys_a)), length(keys_a)),
              rep(sprintf("%s (%d genes)", lab_b, length(keys_b)), length(keys_b))),
            levels = c(sprintf("%s (%d genes)", lab_a, length(keys_a)),
                       sprintf("%s (%d genes)", lab_b, length(keys_b)))
        )
        ct <- sub("__.*$", "", keys)
        ra <- rowAnnotation(
            celltype = ct,
            col = list(celltype = celltype_colors),
            annotation_name_side = "top",
            show_annotation_name = TRUE,
            annotation_name_gp = gpar(fontsize = 13),
            annotation_legend_param = list(
                title_gp = gpar(fontsize = 13),
                labels_gp = gpar(fontsize = 12)
            )
        )
        ht <- Heatmap(
            mat_both, name = "z",
            col = col_fun,
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE,
            show_column_names = TRUE,
            column_labels = c(ages, ages),
            column_split = col_sp,
            row_split = row_sp,
            cluster_row_slices = FALSE,
            cluster_column_slices = FALSE,
            row_gap = unit(6, "mm"),
            column_gap = unit(6, "mm"),
            column_title_gp = gpar(fontsize = 18, fontface = "bold"),
            row_title_gp = gpar(fontsize = 14),
            row_title_rot = 90,
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 12),
            right_annotation = ra,
            border = TRUE,
            use_raster = TRUE, raster_quality = 5,
            rect_gp = gpar(col = NA),
            width = unit(96, "mm"),
            height = unit(118, "mm"),
            heatmap_legend_param = list(
                title = "z", border = "black",
                at = c(-2, -1, 0, 1, 2),
                title_gp = gpar(fontsize = 14),
                labels_gp = gpar(fontsize = 12),
                legend_height = unit(3.2, "cm")
            )
        )
        png(outfile, width = 7.2, height = 8.2, units = "in", res = 300)
        draw(ht, column_title = panel,
             column_title_gp = gpar(fontsize = 22, fontface = "bold"),
             padding = unit(c(2, 2, 2, 2), "mm"))
        dev.off()
        message("  saved: ", outfile)
    }
    draw_pair_ht(groups$Young, "Young", groups$Old, "Old",
                 file.path(plot_dir, "Human_sex_overlap_Common_2026Aug.png"),
                 "Common")
    draw_pair_ht(groups$`F:Old M:Young`, "F:Old M:Young",
                 groups$`F:Young M:Old`, "F:Young M:Old",
                 file.path(plot_dir, "Human_sex_overlap_Divergent_2026Aug.png"),
                 "Divergent")
    invisible(list(groups = groups, genes = gene_tab))
}


########################################################################
# Human sex overlap → GO BP bubble（top 10；x = −log10(p)，size = Count）
########################################################################
plot_human_sex_go <- function(groups, mat, out_dir) {
    message("\n========== Human sex overlap GO BP bubble ==========")
    if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
        warning("缺 clusterProfiler / org.Hs.eg.db，跳过 GO。可 conda activate clusterProfiler 再跑")
        return(invisible(NULL))
    }
    suppressMessages({
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(AnnotationDbi)
    })
    key_to_genes <- function(keys) unique(sub(".*__", "", keys))
    bg_sym <- unique(sub(".*__", "", rownames(mat)))
    bg_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = bg_sym,
                                       column = "ENTREZID", keytype = "SYMBOL",
                                       multiVals = "first")
    bg_entrez <- unique(as.character(na.omit(bg_entrez)))
    panel_spec <- list(
        list(key = "Young",         title = "Common Young",            fill = "#5B9BD5"),
        list(key = "Old",           title = "Common Old",              fill = "#E64B35"),
        list(key = "F:Old M:Young", title = "Divergent F:Old M:Young", fill = "#5B9BD5"),
        list(key = "F:Young M:Old", title = "Divergent F:Young M:Old", fill = "#E64B35")
    )
    wrap_term <- function(x, width = 42) paste(strwrap(x, width = width), collapse = "\n")
    parse_ratio <- function(r) {
        parts <- as.numeric(strsplit(as.character(r), "/")[[1]])
        if (length(parts) < 2 || !is.finite(parts[2]) || parts[2] == 0) return(NA_real_)
        parts[1] / parts[2]
    }
    run_one_go <- function(keys, title) {
        genes <- key_to_genes(keys)
        message(sprintf("  %s: keys=%d  unique genes=%d", title, length(keys), length(genes)))
        if (length(genes) < 5) {
            message("    基因太少，跳过")
            return(NULL)
        }
        entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes,
                                        column = "ENTREZID", keytype = "SYMBOL",
                                        multiVals = "first")
        entrez <- unique(as.character(na.omit(entrez)))
        message(sprintf("    Entrez=%d", length(entrez)))
        if (length(entrez) < 5) return(NULL)
        ego <- tryCatch(enrichGO(
            gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
            keyType = "ENTREZID", universe = bg_entrez,
            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
            readable = TRUE
        ), error = function(e) NULL)
        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
            message("    无显著 GO")
            return(NULL)
        }
        d <- as.data.frame(ego)
        d <- d[order(d$pvalue, d$p.adjust), ]
        d$GeneRatio_val <- vapply(d$GeneRatio, parse_ratio, numeric(1))
        d$neglog10_p <- -log10(pmax(d$pvalue, 1e-300))
        d$group <- title
        d
    }
    make_bubble <- function(d, title, fill) {
        df <- head(d, 10)
        df$lab <- vapply(df$Description, wrap_term, character(1))
        df <- df[order(df$neglog10_p, decreasing = FALSE), ]
        df$lab <- factor(df$lab, levels = unique(df$lab))
        ggplot(df, aes(x = neglog10_p, y = lab, size = Count, color = neglog10_p)) +
            geom_point(alpha = 0.9) +
            scale_color_gradient(low = "grey80", high = fill, name = expression(-log[10](p))) +
            scale_size_continuous(range = c(3, 9), name = "Count") +
            scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
            labs(title = title, x = expression(-log[10](p-value)), y = NULL) +
            theme_classic(base_size = 12) +
            theme(
                plot.title = element_text(size = 13, face = "bold", hjust = 0),
                axis.text.y = element_text(size = 10, color = "black", lineheight = 0.95),
                axis.text.x = element_text(size = 10, color = "black"),
                axis.title.x = element_text(size = 11),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
                legend.position = "right",
                plot.margin = margin(6, 8, 6, 6)
            )
    }
    go_tabs <- list()
    for (sp in panel_spec) {
        d <- run_one_go(groups[[sp$key]], sp$title)
        if (is.null(d) || nrow(d) == 0) next
        go_tabs[[sp$title]] <- d
        fp <- file.path(out_dir, sprintf("Human_sex_GO_bubble_%s_2026Aug.png",
                                         gsub("[^A-Za-z0-9]+", "_", sp$title)))
        ggsave(fp, make_bubble(d, sp$title, sp$fill),
               width = 6.4, height = max(3.6, 0.38 * min(10, nrow(d)) + 1.4),
               dpi = 300, bg = "white")
        message("  saved: ", fp)
    }
    if (length(go_tabs)) {
        all_go <- do.call(rbind, go_tabs)
        write.table(all_go, file.path(out_dir, "Human_sex_GO_enrich_2026Aug.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        saveRDS(go_tabs, file.path(out_dir, "Human_sex_GO_enrich_2026Aug.rds"))
    }
    invisible(go_tabs)
}


########################################################################
# 核对用：打印 Sheet A / Sheet C gene 个数
########################################################################
print_table2_gene_counts <- function(deg_res, fgh = NULL, edger = NULL) {
    message("\n################################################################")
    message("### Table2 核对：Sheet A  gene 个数 (Young / Middle / Old) ###")
    message("################################################################")
    for (sp in c("Zebrafish", "Mouse", "Human")) {
        x <- deg_res[[sp]]
        if (is.null(x) || is.null(x$km)) {
            message("\n[", sp, "] 缺失，跳过")
            next
        }
        job <- jobs[[sp]]
        ga <- km_to_gene_table(sp, x$km, job)
        ga <- ga[!is.na(ga$celltype) & nzchar(ga$celltype) &
                     ga$direction %in% c("Young", "Middle", "Old"), , drop = FALSE]
        message(sprintf("\n========== %s  (total rows=%d) ==========", sp, nrow(ga)))
        if (isTRUE(job$split_sex) && any(!is.na(ga$sex) & nzchar(ga$sex))) {
            for (sx in c("Male", "Female")) {
                g2 <- ga[!is.na(ga$sex) & ga$sex == sx, , drop = FALSE]
                message(sprintf("--- %s / %s (n=%d) ---", sp, sx, nrow(g2)))
                tab <- as.data.frame.matrix(table(g2$celltype, g2$direction))
                for (d in c("Young", "Middle", "Old")) {
                    if (!d %in% names(tab)) tab[[d]] <- 0L
                }
                tab <- tab[, c("Young", "Middle", "Old"), drop = FALSE]
                tab$Total <- rowSums(tab)
                tab <- tab[order(match(rownames(tab), ct_order_hmz), rownames(tab)), ,
                           drop = FALSE]
                print(tab)
                message(sprintf("  %s %s sum: Young=%d Middle=%d Old=%d Total=%d",
                                sp, sx,
                                sum(tab$Young), sum(tab$Middle), sum(tab$Old),
                                sum(tab$Total)))
            }
        } else {
            tab <- as.data.frame.matrix(table(ga$celltype, ga$direction))
            for (d in c("Young", "Middle", "Old")) {
                if (!d %in% names(tab)) tab[[d]] <- 0L
            }
            tab <- tab[, c("Young", "Middle", "Old"), drop = FALSE]
            tab$Total <- rowSums(tab)
            tab <- tab[order(match(rownames(tab), ct_order_hmz), rownames(tab)), ,
                       drop = FALSE]
            print(tab)
            message(sprintf("  %s sum: Young=%d Middle=%d Old=%d Total=%d",
                            sp,
                            sum(tab$Young), sum(tab$Middle), sum(tab$Old),
                            sum(tab$Total)))
        }
    }

    message("\n################################################################")
    message("### Table2 核对：Sheet C  跨物种 overlap gene 个数          ###")
    message("################################################################")
    if (is.null(fgh) || is.null(fgh$sizes) || !nrow(fgh$sizes)) {
        message("  (无 CrossSp_set_sizes)")
    } else {
        st <- fgh$sizes
        message("\n--- per celltype × direction (n_H/M/Z + pairwise + triple) ---")
        print(st, row.names = FALSE)
        message("\n--- triple gene 个数汇总 (celltype × Young/Old) ---")
        trip_w <- reshape2::dcast(st, celltype ~ direction, value.var = "n_triple",
                                  fill = 0)
        if (!"Young" %in% names(trip_w)) trip_w$Young <- 0
        if (!"Old" %in% names(trip_w)) trip_w$Old <- 0
        trip_w <- trip_w[, c("celltype", "Young", "Old"), drop = FALSE]
        trip_w <- trip_w[order(match(trip_w$celltype, ct_order_hmz),
                               trip_w$celltype), , drop = FALSE]
        print(trip_w, row.names = FALSE)
        message(sprintf("  triple sum: Young=%d Old=%d Total=%d",
                        sum(st$n_triple[st$direction == "Young"], na.rm = TRUE),
                        sum(st$n_triple[st$direction == "Old"], na.rm = TRUE),
                        sum(st$n_triple, na.rm = TRUE)))
        if (!is.null(fgh$triple) && nrow(fgh$triple)) {
            message("\n--- overlap gene 表行数 by set ---")
            print(as.data.frame(table(set = fgh$triple$set,
                                      direction = fgh$triple$direction),
                                stringsAsFactors = FALSE))
        }
    }
    message("\n################################################################")
    message("### Table2 核对：Sheet D  edgeR Old vs Young gene 个数       ###")
    message("################################################################")
    if (is.null(edger) || is.null(edger$summary) || !nrow(edger$summary)) {
        message("  (无 edgeR 结果)")
    } else {
        sm <- edger$summary
        message("\n--- edgeR_summary (PValue<0.05 / FDR<0.05) ---")
        print(sm, row.names = FALSE)
        for (sp in unique(sm$species)) {
            s2 <- sm[sm$species == sp, , drop = FALSE]
            message(sprintf(
                "  %s total P05: Up_Old=%d Down_Old=%d | FDR05: Up=%d Down=%d | genes_tested=%d",
                sp,
                sum(s2$n_P05_up, na.rm = TRUE),
                sum(s2$n_P05_down, na.rm = TRUE),
                sum(s2$n_FDR05_up, na.rm = TRUE),
                sum(s2$n_FDR05_down, na.rm = TRUE),
                sum(s2$n_genes, na.rm = TRUE)
            ))
        }
    }
    message("\n################################################################")
    invisible(NULL)
}

########################################################################
# Sheet D: edgeR Old vs Young（读已有结果，不重跑）
########################################################################
load_edger_object <- function(cands, prefer_name = NULL) {
    fp <- find_existing(cands)
    if (is.na(fp)) return(NULL)
    message("[edgeR load] ", fp)
    if (grepl("\\.rds$", fp, ignore.case = TRUE)) return(readRDS(fp))
    env <- new.env(parent = emptyenv())
    load(fp, envir = env)
    if (!is.null(prefer_name) && exists(prefer_name, envir = env, inherits = FALSE))
        return(get(prefer_name, envir = env))
    objs <- ls(envir = env)
    if (!length(objs)) return(NULL)
    get(objs[[1]], envir = env)
}

normalize_edger_df <- function(df, species, celltype, sex = NA_character_) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    need <- c("Gene", "logFC", "PValue")
    if (!all(need %in% names(df))) {
        warning(sprintf("edgeR %s/%s 缺列: %s", species, celltype,
                        paste(setdiff(need, names(df)), collapse = ",")))
        return(NULL)
    }
    out <- df
    if (!"direction" %in% names(out)) {
        out$direction <- ifelse(out$logFC > 0, "Up_Old", "Down_Old")
    }
    if (!"FDR" %in% names(out) && "FDR" %in% names(df)) out$FDR <- df$FDR
    if (!"FDR" %in% names(out) && "adj.P.Val" %in% names(out)) out$FDR <- out$adj.P.Val
    out$species  <- species
    out$celltype <- celltype
    out$sex      <- sex
    ## 列顺序：meta + 常见 edgeR 列 + 其余
    front <- c("species", "celltype", "sex", "Gene", "logFC", "logCPM", "F",
               "PValue", "FDR", "direction")
    keep <- c(intersect(front, names(out)), setdiff(names(out), front))
    out[, keep, drop = FALSE]
}

edger_one_summary <- function(df) {
    data.frame(
        species = df$species[1],
        celltype = df$celltype[1],
        sex = if (all(is.na(df$sex))) NA_character_ else as.character(df$sex[1]),
        n_genes = nrow(df),
        n_P05_up = sum(df$PValue < 0.05 & df$logFC > 0, na.rm = TRUE),
        n_P05_down = sum(df$PValue < 0.05 & df$logFC < 0, na.rm = TRUE),
        n_FDR05 = if ("FDR" %in% names(df)) sum(df$FDR < 0.05, na.rm = TRUE) else NA_integer_,
        n_FDR05_up = if ("FDR" %in% names(df))
            sum(df$FDR < 0.05 & df$logFC > 0, na.rm = TRUE) else NA_integer_,
        n_FDR05_down = if ("FDR" %in% names(df))
            sum(df$FDR < 0.05 & df$logFC < 0, na.rm = TRUE) else NA_integer_,
        stringsAsFactors = FALSE
    )
}

collect_edger_tables <- function() {
    message("\n========== 读入 edgeR Old vs Young ==========")
    message("  Zebrafish: Young<=3 mo, Old>=36 mo")
    message("  Mouse:     Young<=12 wk, Old>=100 wk")
    message("  Human:     Young<=40 yr, Old>=80 yr (by sex)")

    long_list <- list()
    push_df <- function(df) {
        if (is.null(df) || !nrow(df)) return(invisible(NULL))
        long_list[[length(long_list) + 1L]] <<- df
        invisible(NULL)
    }

    ## Zebrafish
    z <- load_edger_object(
        c(file.path(fish_dir, "Zebrafish_DEG_edgeR"),
          file.path(fish_dir, "Zebrafish_DEG_edgeR.rds"),
          file.path(fish_dir, "Zebrafish_DEG_edgeR_2026.rds")),
        prefer_name = "Zebrafish_DEG_edgeR"
    )
    if (is.null(z)) {
        warning("找不到 Zebrafish_DEG_edgeR，跳过鱼 edgeR")
    } else {
        for (ct in names(z)) {
            push_df(normalize_edger_df(z[[ct]], "Zebrafish", ct))
        }
    }

    ## Mouse
    m <- load_edger_object(
        c(file.path(mouse_dir, "Mouse_DEG_edgeR_2026.rds"),
          file.path(mouse_dir, "Mouse_DEG_edgeR_2026.RData"),
          file.path(mouse_dir, "Mouse_DEG_edgeR")),
        prefer_name = "Mouse_DEG_edgeR"
    )
    if (is.null(m)) {
        warning("找不到 Mouse_DEG_edgeR*，跳过鼠 edgeR；可先跑 Run_Mouse_edgeR_DEGs_2026.R")
    } else {
        for (ct in names(m)) {
            push_df(normalize_edger_df(m[[ct]], "Mouse", ct))
        }
    }

    ## Human by sex
    h_fp <- file.path(human_pb, "Human_DEG_edgeR_bySex.rds")
    h <- if (file.exists(h_fp)) {
        message("[edgeR load] ", h_fp)
        readRDS(h_fp)
    } else {
        warning("找不到 Human_DEG_edgeR_bySex.rds，跳过人 edgeR")
        NULL
    }
    if (!is.null(h)) {
        for (ct in names(h)) {
            for (sx in names(h[[ct]])) {
                sx_lab <- if (tolower(sx) %in% c("m", "male")) "Male" else
                    if (tolower(sx) %in% c("f", "female")) "Female" else as.character(sx)
                push_df(normalize_edger_df(h[[ct]][[sx]], "Human", ct, sex = sx_lab))
            }
        }
    }

    if (!length(long_list)) return(NULL)

    ## 对齐列后 rbind
    all_cols <- unique(unlist(lapply(long_list, names)))
    long_list <- lapply(long_list, function(d) {
        miss <- setdiff(all_cols, names(d))
        for (nm in miss) d[[nm]] <- NA
        d[, all_cols, drop = FALSE]
    })
    all_df <- do.call(rbind, long_list)
    rownames(all_df) <- NULL

    ## CT 排序
    all_df$celltype <- factor(all_df$celltype, levels = ct_order_hmz)
    all_df <- all_df[order(all_df$species, all_df$celltype, all_df$sex, all_df$PValue), ]
    all_df$celltype <- as.character(all_df$celltype)

    ## 按 species|celltype|sex 汇总
    keys <- paste(all_df$species, all_df$celltype, all_df$sex, sep = "||")
    summary_df <- do.call(rbind, lapply(split(all_df, keys), edger_one_summary))
    rownames(summary_df) <- NULL
    summary_df$celltype <- factor(summary_df$celltype, levels = ct_order_hmz)
    summary_df <- summary_df[order(summary_df$species, summary_df$celltype, summary_df$sex), ]
    summary_df$celltype <- as.character(summary_df$celltype)

    by_sp <- split(all_df, all_df$species)
    message(sprintf("  edgeR rows: Zebrafish=%d Mouse=%d Human=%d | summary rows=%d",
                    if (!is.null(by_sp$Zebrafish)) nrow(by_sp$Zebrafish) else 0L,
                    if (!is.null(by_sp$Mouse)) nrow(by_sp$Mouse) else 0L,
                    if (!is.null(by_sp$Human)) nrow(by_sp$Human) else 0L,
                    nrow(summary_df)))
    list(all = all_df, by_species = by_sp, summary = summary_df)
}
safe_sheet <- function(x) substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", x), 1, 31)

km_to_gene_table <- function(species, km, job) {
    meta <- parse_row_meta(km$genes)
    dirn <- ifelse(km$cluster %in% job$Down, "Young",
            ifelse(km$cluster %in% job$Up, "Old",
            ifelse(km$cluster %in% job$Middle, "Middle", "Other")))
    data.frame(
        species   = species,
        gene_id   = km$genes,
        celltype  = meta$celltype,
        sex       = meta$sex,
        gene      = sub(".*__", "", km$genes),
        cluster   = km$cluster,
        direction = dirn,
        stringsAsFactors = FALSE
    )
}

write_deg_excel <- function(deg_res, fgh, sex_out, go_tabs, fp, edger = NULL) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
        warning("openxlsx not installed; skip Excel")
        return(invisible(NULL))
    }
    library(openxlsx)
    sheets <- list()
    push <- function(tab, desc, df) {
        if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(invisible(NULL))
        sheets[[length(sheets) + 1]] <<- list(
            Tab = safe_sheet(tab), Description = desc, df = df,
            N_rows = nrow(df), N_cols = ncol(df)
        )
        invisible(NULL)
    }

    grp <- do.call(rbind, lapply(names(jobs), function(sp) {
        job <- jobs[[sp]]
        data.frame(
            species = sp,
            Young_clusters  = paste(job$Down, collapse = ","),
            Middle_clusters = paste(job$Middle, collapse = ","),
            Old_clusters    = paste(job$Up, collapse = ","),
            dropped_raw_kmeans = paste(job$rm_clusters, collapse = ","),
            stringsAsFactors = FALSE
        )
    }))
    push("Kmeans_groups",
         "Remapped cluster IDs used as Young / Middle / Old in heatmaps and barplots.",
         grp)

    gene_all <- list()
    bar_all  <- list()
    for (sp in names(deg_res)) {
        x <- deg_res[[sp]]
        if (is.null(x) || !is.list(x) || is.null(x$km)) next
        job <- jobs[[sp]]
        gt <- km_to_gene_table(sp, x$km, job)
        gene_all[[sp]] <- gt
        if (isTRUE(job$split_sex)) {
            sex_vec <- parse_row_meta(x$km$genes)$sex
            kc_M <- x$km[!is.na(sex_vec) & sex_vec == "Male", , drop = FALSE]
            kc_F <- x$km[!is.na(sex_vec) & sex_vec == "Female", , drop = FALSE]
            dM <- make_bar_data(kc_M, job$Down, job$Middle, job$Up)
            dF <- make_bar_data(kc_F, job$Down, job$Middle, job$Up)
            dM$species <- sp; dM$sex <- "Male"
            dF$species <- sp; dF$sex <- "Female"
            bar_all[[paste0(sp, "_M")]] <- dM
            bar_all[[paste0(sp, "_F")]] <- dF
        } else {
            d <- make_bar_data(x$km, job$Down, job$Middle, job$Up)
            d$species <- sp; d$sex <- NA_character_
            bar_all[[sp]] <- d
        }
    }
    if (length(gene_all)) {
        ga <- do.call(rbind, gene_all); rownames(ga) <- NULL
        push("Gene_clusters",
             "Each DEG row after kmeans remapping: species, cell type, gene, cluster, Young/Middle/Old.",
             ga)
        cl_n <- as.data.frame(table(ga$species, ga$cluster, ga$direction),
                              stringsAsFactors = FALSE)
        names(cl_n) <- c("species", "cluster", "direction", "n_genes")
        cl_n <- cl_n[cl_n$n_genes > 0, ]
        push("Cluster_gene_counts",
             "Number of DEGs per remapped cluster and Young/Middle/Old direction.",
             cl_n)
    }
    if (length(bar_all)) {
        bd <- do.call(rbind, bar_all); rownames(bd) <- NULL
        bd <- bd[, c("species", "sex", "cell_type", "direction", "num_DEGs")]
        push("DEG_barplot_counts",
             "Stacked-bar counts: number of Young / Middle / Old DEGs per cell type (Human split by sex).",
             bd)
    }
    if (!is.null(fgh) && is.list(fgh)) {
        push("CrossSp_overlap_genes",
             "Between-species overlap genes (human symbols) for ALL cell types in Gene_clusters. set = HM / HZ / MZ / triple; direction = Young (Down) or Old (Up). Same DEG lists as Sheet Gene_clusters.",
             fgh$triple)
        push("CrossSp_set_sizes",
             "Per celltype × Young/Old: n_H/M/Z after ortholog mapping, plus pairwise n_HM/n_HZ/n_MZ and n_triple. All cell types from Gene_clusters.",
             fgh$sizes)
    }
    if (!is.null(edger) && is.list(edger)) {
        push("edgeR_summary",
             "edgeR Old vs Young counts per species × cell type (Human split by sex). n_P05_* = PValue<0.05; n_FDR05_* = FDR<0.05. Up = higher in Old. Cutoffs: Zebrafish Young<=3mo Old>=36mo; Mouse Young<=12wk Old>=100wk; Human Young<=40yr Old>=80yr.",
             edger$summary)
        if (!is.null(edger$by_species$Zebrafish))
            push("edgeR_Zebrafish",
                 "Full edgeR table: Zebrafish Old vs Young per cell type (logFC>0 = Up_Old). All tested genes.",
                 edger$by_species$Zebrafish)
        if (!is.null(edger$by_species$Mouse))
            push("edgeR_Mouse",
                 "Full edgeR table: Mouse Old vs Young per cell type (logFC>0 = Up_Old). All tested genes.",
                 edger$by_species$Mouse)
        if (!is.null(edger$by_species$Human))
            push("edgeR_Human",
                 "Full edgeR table: Human Old vs Young per cell type × sex (logFC>0 = Up_Old). All tested genes.",
                 edger$by_species$Human)
    }
    if (!is.null(sex_out) && is.list(sex_out)) {
        gtab <- sex_out$genes
        if (!is.null(gtab) && nrow(gtab)) {
            gtab$class[gtab$class == "Young"] <- "Common_Young"
            gtab$class[gtab$class == "Old"] <- "Common_Old"
            gtab$class[gtab$class == "F:Old M:Young"] <- "Divergent_FOld_MYoung"
            gtab$class[gtab$class == "F:Young M:Old"] <- "Divergent_FYoung_MOld"
            push("Human_sex_overlap",
                 "Human male vs female Young/Old overlap. Common = same direction in both sexes; Divergent = opposite. Key = celltype__gene.",
                 gtab)
            cnt <- as.data.frame(table(gtab$class), stringsAsFactors = FALSE)
            names(cnt) <- c("class", "n_genes")
            push("Human_sex_overlap_n",
                 "Gene counts for Common Young, Common Old, and the two Divergent classes.",
                 cnt)
        }
    }
    if (!is.null(go_tabs) && length(go_tabs)) {
        go <- do.call(rbind, lapply(names(go_tabs), function(nm) {
            d <- go_tabs[[nm]]
            if (is.null(d) || !nrow(d)) return(NULL)
            keep <- intersect(c("group", "ID", "Description", "GeneRatio", "BgRatio",
                                "pvalue", "p.adjust", "qvalue", "geneID", "Count",
                                "GeneRatio_val", "neglog10_p"), names(d))
            d[, keep, drop = FALSE]
        }))
        rownames(go) <- NULL
        push("Human_sex_GO",
             "GO BP enrichment (enrichGO) for the four Human sex-overlap gene sets. Bubble plots use top 10 by p-value.",
             go)
    }

    if (!length(sheets)) {
        warning("No tables to write")
        return(invisible(NULL))
    }
    wb <- createWorkbook()
    hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
    title_style <- createStyle(textDecoration = "bold", fontSize = 14)
    addWorksheet(wb, "Contents")
    contents <- data.frame(
        Sheet = vapply(sheets, `[[`, "", "Tab"),
        N_rows = vapply(sheets, `[[`, 0L, "N_rows"),
        N_cols = vapply(sheets, `[[`, 0L, "N_cols"),
        Description = vapply(sheets, `[[`, "", "Description"),
        stringsAsFactors = FALSE
    )
    writeData(wb, "Contents",
              "HMZ DEG figures 2026Aug — source tables",
              startRow = 1, startCol = 1)
    addStyle(wb, "Contents", title_style, rows = 1, cols = 1)
    writeData(wb, "Contents",
              sprintf("Output dir: %s | %s", plot_dir, format(Sys.time(), "%Y-%m-%d %H:%M")),
              startRow = 2, startCol = 1)
    writeData(wb, "Contents",
              "Young/Old/Middle = remapped kmeans. Cross-species overlap (all CT in Gene_clusters) uses human orthologs (HM/HZ/MZ/triple). edgeR sheets = Old vs Young per species × cell type (Human by sex). Human sex overlap matches celltype__gene across male and female.",
              startRow = 3, startCol = 1)
    writeData(wb, "Contents", contents, startRow = 5, startCol = 1, headerStyle = hs)
    for (i in seq_along(sheets)) {
        spec <- sheets[[i]]
        addWorksheet(wb, spec$Tab)
        writeData(wb, spec$Tab, spec$df, headerStyle = hs)
        freezePane(wb, spec$Tab, firstActiveRow = 2)
        setColWidths(wb, spec$Tab, cols = seq_len(min(ncol(spec$df), 20)), widths = "auto")
    }
    freezePane(wb, "Contents", firstActiveRow = 6)
    setColWidths(wb, "Contents", cols = 1:4, widths = c(28, 10, 10, 90))
    saveWorkbook(wb, fp, overwrite = TRUE)
    message("Excel: ", fp)
    invisible(fp)
}


########################################################################
# run
########################################################################
sp_names <- c("Zebrafish", "Mouse", "Human")
message("并行核数 NCORES=", ncores, " （三物种同时出图）")
deg_res <- safe_mclapply(sp_names, function(sp) {
    tryCatch(run_one(sp, jobs[[sp]]), error = function(e) {
        warning(sp, " 失败: ", conditionMessage(e))
        NULL
    })
}, mc.cores = min(3L, ncores))
names(deg_res) <- sp_names
for (sp in sp_names) {
    if (inherits(deg_res[[sp]], "try-error")) {
        warning(sp, " 子进程失败: ", as.character(deg_res[[sp]]))
        deg_res[[sp]] <- NULL
    }
}
fgh <- tryCatch(plot_fgh_overlap(deg_res), error = function(e) {
    warning("F/G/H overlap 失败: ", conditionMessage(e))
    NULL
})
edger <- tryCatch(collect_edger_tables(), error = function(e) {
    warning("edgeR 读入失败: ", conditionMessage(e))
    NULL
})
sex_out <- NULL
go_tabs <- NULL
if (!is.null(deg_res$Human) && is.list(deg_res$Human)) {
    sex_out <- tryCatch(
        plot_human_sex_overlap(deg_res$Human$mat, deg_res$Human$km, jobs$Human),
        error = function(e) {
            warning("Human sex overlap 失败: ", conditionMessage(e))
            NULL
        })
    sex_groups <- if (is.list(sex_out) && !is.null(sex_out$groups)) sex_out$groups else sex_out
    if (!is.null(sex_groups)) {
        go_tabs <- tryCatch(plot_human_sex_go(sex_groups, deg_res$Human$mat, plot_dir),
                            error = function(e) {
                                warning("Human sex GO 失败: ", conditionMessage(e))
                                NULL
                            })
    }
}
tryCatch(write_deg_excel(deg_res, fgh, sex_out, go_tabs,
                         file.path(plot_dir, "HMZ_DEGs_tables_2026Aug.xlsx"),
                         edger = edger),
         error = function(e) warning("Excel 失败: ", conditionMessage(e)))

## 末尾集中打印，方便核对 Table2 gene 个数
tryCatch(print_table2_gene_counts(deg_res, fgh, edger),
         error = function(e) warning("打印 gene 个数失败: ", conditionMessage(e)))

message("\n=== Done ===")
message("全部输出: ", plot_dir)
message("Excel:    ", file.path(plot_dir, "HMZ_DEGs_tables_2026Aug.xlsx"))
