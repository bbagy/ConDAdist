#' Generate Ranked Lollipop Plots for DA and ConDA Results
#'
#' This function follows the same directory-first contract as `Go_volcanoPlot()`.
#' It reads CSV result files or bridge CSV files, recognizes the originating
#' DA tool, and generates ranked lollipop plots into `pdf/DA_plot` or
#' `pdf/ConDa_plot`.
#'
#' @param project Project name or identifier.
#' @param file_path Optional directory path to scan for CSV files.
#' @param files Optional filename pattern used together with `file_path`.
#' @param result Optional result directory returned by a Go DA-family tool or a
#'   ConDA bridge directory.
#' @param p_cutoff Numeric cutoff used to keep features. Adjusted p-values
#'   (\code{qval}) are used when available; otherwise raw \code{pval} is used.
#' @param mycols Optional two-color vector for negative / positive direction.
#' @param name Optional plot label used in output filenames.
#' @param font Base font size.
#' @param height Deprecated. Height is now determined automatically from the
#'   number of displayed features.
#' @param width Fixed panel-oriented base PDF width. Extra width for long
#'   feature labels is added automatically.
#'
#' @return The function saves PDF files and returns `NULL` invisibly.
#'
#' @param patchwork Logical. If \code{TRUE}, skip saving and return the plot object(s) for use with \code{Gg_patchwork()} or the \pkg{patchwork} package. Default \code{FALSE}.
#' @export
Go_lollipopPlot <- function(project,
                            file_path = NULL,
                            files = NULL,
                            result = NULL,
                            p_cutoff = 0.05,
                            mycols = NULL,
                            name = NULL,
                            font = 10,
                            height = NULL,
                            width = 2.2,
                            patchwork = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  if (!is.null(height)) {
    message("[Go_lollipopPlot] 'height' is deprecated and ignored. Plot height is auto-calculated.")
  }

  first_or_na <- function(x) {
    if (length(x) == 0 || all(is.na(x))) {
      return(NA_character_)
    }
    x[[1]]
  }

  build_optional_metric_note <- function(df) {
    first_metric <- function(candidates) {
      hits <- candidates[candidates %in% colnames(df)]
      if (length(hits) == 0) {
        return(NA_real_)
      }
      for (col in hits) {
        vals <- suppressWarnings(as.numeric(df[[col]]))
        vals <- vals[is.finite(vals) & !is.na(vals)]
        if (length(vals) > 0) {
          return(vals[1])
        }
      }
      NA_real_
    }

    auroc <- first_metric(c("auroc", "AUROC", "mean_auroc", "roc", "ROC"))
    auprc <- first_metric(c("auprc", "AUPRC", "mean_auprc", "prc", "PRC"))

    parts <- c(
      if (is.finite(auroc)) sprintf("ROC=%.3f", auroc) else NULL,
      if (is.finite(auprc)) sprintf("PRC=%.3f", auprc) else NULL
    )
    if (length(parts) == 0) {
      return("")
    }
    paste(parts, collapse = " | ")
  }

  calc_label_width_in <- function(labels, font_size_pt) {
    labels <- as.character(labels)
    labels <- labels[!is.na(labels)]
    if (length(labels) == 0) {
      return(1)
    }
    longest <- labels[which.max(nchar(labels))]
    width_in <- tryCatch(
      grid::convertWidth(
        grid::grobWidth(grid::textGrob(longest, gp = grid::gpar(fontsize = font_size_pt, fontface = "italic"))),
        "in",
        valueOnly = TRUE
      ),
      error = function(e) NA_real_
    )
    if (!is.finite(width_in)) {
      width_in <- max(1, nchar(longest) * 0.08)
    }
    width_in
  }

  fix_panel_width_grob <- function(plot_obj, panel_width_in) {
    gt <- ggplot2::ggplotGrob(plot_obj)
    panel_rows <- gt$layout[gt$layout$name == "panel", , drop = FALSE]
    if (nrow(panel_rows) > 0) {
      panel_cols <- unique(unlist(Map(seq.int, panel_rows$l, panel_rows$r)))
      gt$widths[panel_cols] <- grid::unit(panel_width_in / length(panel_cols), "in")
    }
    gt
  }

  build_taxonomy_plot_label <- function(df) {
    if (nrow(df) == 0) {
      return(character(0))
    }

    pick_first_valid <- function(candidates) {
      vals <- as.character(candidates)
      vals <- vals[!is.na(vals) & nzchar(trimws(vals)) & !vals %in% c("NA", "NA NA", "__", "s__", "g__", "f__", "o__", "c__", "p__", "k__")]
      if (length(vals) == 0) {
        return(NA_character_)
      }
      vals[1]
    }

    taxonomy_priority <- intersect(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "TaxaName"), colnames(df))
    asv_col <- if ("ASV" %in% colnames(df)) {
      "ASV"
    } else if ("Row.names" %in% colnames(df)) {
      "Row.names"
    } else {
      NULL
    }

    labels <- vapply(seq_len(nrow(df)), function(i) {
      vals <- vapply(taxonomy_priority, function(col) df[i, col], character(1))
      lbl <- pick_first_valid(vals)
      if (is.na(lbl) && !is.null(asv_col)) {
        lbl <- as.character(df[i, asv_col])
      }
      if (is.na(lbl) || !nzchar(trimws(lbl))) {
        lbl <- paste0("Feature_", i)
      }
      lbl
    }, character(1))

    labels
  }

  build_feature_label <- function(df) {
    if (nrow(df) == 0) {
      return(character(0))
    }
    if (all(c("KOName") %in% colnames(df))) {
      vals <- as.character(df$KOName)
      vals[is.na(vals) | !nzchar(trimws(vals))] <- if ("KO" %in% colnames(df)) as.character(df$KO[is.na(vals) | !nzchar(trimws(vals))]) else NA_character_
      vals[is.na(vals) | !nzchar(trimws(vals))] <- if ("ASV" %in% colnames(df)) as.character(df$ASV[is.na(vals) | !nzchar(trimws(vals))]) else paste0("Feature_", seq_len(nrow(df)))[is.na(vals) | !nzchar(trimws(vals))]
      return(vals)
    }
    if (all(c("symbol") %in% colnames(df))) {
      vals <- as.character(df$symbol)
      vals[is.na(vals) | !nzchar(trimws(vals))] <- if ("ASV" %in% colnames(df)) as.character(df$ASV[is.na(vals) | !nzchar(trimws(vals))]) else paste0("Feature_", seq_len(nrow(df)))[is.na(vals) | !nzchar(trimws(vals))]
      return(vals)
    }
    build_taxonomy_plot_label(df)
  }

  detect_tool <- function(filename1) {
    if (grepl("condadist", filename1)) return("condadist")
    if (grepl("wilcoxon", filename1)) return("wilcoxon")
    if (grepl("deseq2", filename1)) return("deseq2")
    if (grepl("aldex2", filename1)) return("aldex2")
    if (grepl("maaslin2", filename1)) return("maaslin2")
    if (grepl("maaslin", filename1)) return("maaslin2")
    if (grepl("corncob", filename1)) return("corncob")
    if (grepl("ancombc2", filename1)) return("ancom2")
    if (grepl("ancom2", filename1)) return("ancom2")
    NA_character_
  }

  first_existing <- function(df, candidates) {
    hits <- candidates[candidates %in% colnames(df)]
    if (length(hits) == 0) {
      return(NA_character_)
    }
    hits[1]
  }

  resolve_aldex_columns <- function(df) {
    est_col <- if ("diff.btw" %in% colnames(df)) {
      "diff.btw"
    } else {
      est_hits <- grep("\\.Est$", colnames(df), value = TRUE)
      est_hits <- est_hits[!grepl("kendall|spearman", est_hits, ignore.case = TRUE)]
      first_or_na(est_hits)
    }

    p_col <- if ("wi.ep" %in% colnames(df)) {
      "wi.ep"
    } else {
      p_hits <- grep("\\.pval$", colnames(df), value = TRUE)
      p_hits <- p_hits[!grepl("kendall|spearman", p_hits, ignore.case = TRUE)]
      first_or_na(p_hits)
    }

    q_col <- first_existing(df, c("wi.eBH", grep("\\.pval\\.padj$", colnames(df), value = TRUE), grep("\\.pval\\.holm$", colnames(df), value = TRUE)))
    list(effect = est_col, p = p_col, q = q_col)
  }

  build_plot_df <- function(df, tool, source_name = NULL) {
    baseline_col <- if ("baseline" %in% colnames(df)) "baseline" else if ("basline" %in% colnames(df)) "basline" else NULL
    basline <- if (!is.null(baseline_col)) unique(as.character(df[[baseline_col]]))[1] else "group1"
    smvar <- if ("smvar" %in% colnames(df)) unique(as.character(df$smvar))[1] else "group2"
    mvar <- if ("mvar" %in% colnames(df)) unique(as.character(df$mvar))[1] else "Group"

    conda_subtitle <- ""
    if (tool == "deseq2") {
      effect_col <- "log2FoldChange"
      p_col <- "pvalue"
      q_col <- "padj"
      title_tool <- "deseq2"
      file_tool <- "deseq2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "aldex2") {
      cols <- resolve_aldex_columns(df)
      effect_col <- cols$effect
      p_col <- cols$p
      q_col <- cols$q
      title_tool <- "aldex2"
      file_tool <- "aldex2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "ancom2") {
      effect_col <- "lfc_ancombc"
      p_col <- "pvalue_ancombc"
      q_col <- "qvalue_ancombc"
      title_tool <- "ancom2"
      file_tool <- "ancom2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "corncob") {
      effect_col <- "corncob_coef"
      p_col <- "corncob_pvalue"
      q_col <- "corncob_qvalue"
      title_tool <- "corncob"
      file_tool <- "corncob"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "maaslin2") {
      effect_col <- "maaslin2_coef"
      p_col <- "maaslin2_pvalue"
      q_col <- "maaslin2_qvalue"
      title_tool <- "maaslin2"
      file_tool <- "maaslin2"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "wilcoxon") {
      effect_col <- "log2FoldChange"
      p_col <- "pvalue"
      q_col <- "padj"
      title_tool <- "wilcoxon"
      file_tool <- "wilcoxon"
      out_subdir <- "DA_plot"
      score_label <- "Signed significance score (-log10 adj. p)"
    } else if (tool == "condadist") {
      effect_col <- "median_effect_size"
      p_col <- "fisher_combined_p"
      q_col <- "fisher_combined_q"
      # CSV 파일명(condadist.{sig}.(...).csv)에서 method/dist 시그니처 추출
      conda_sig <- if (!is.null(source_name)) {
        m <- regmatches(source_name, regexpr("^condadist\\.(.+?)(?=\\.\\()", source_name, perl = TRUE))
        if (length(m) == 1 && nzchar(m)) sub("^condadist\\.", "", m) else ""
      } else ""
      # 시그니처를 method 부분(대문자)과 distance 부분(소문자)으로 분리
      conda_parts <- if (nzchar(conda_sig)) strsplit(conda_sig, "\\.")[[1]] else character(0)
      conda_method_part <- if (length(conda_parts) >= 1) conda_parts[1] else ""
      conda_dist_parts  <- if (length(conda_parts) >= 2) conda_parts[-1] else character(0)
      # 이미지 title: DANMC 등 method abbreviation만 사용
      title_tool <- if (nzchar(conda_method_part)) conda_method_part else "ConDA-dist"
      # subtitle: dist 정보 (있을 때만)
      conda_subtitle <- if (length(conda_dist_parts) > 0) paste("dist:", paste(conda_dist_parts, collapse = " \u00b7 ")) else ""
      file_tool <- if (nzchar(conda_sig)) paste0("ConDA-dist.", conda_sig) else "ConDA-dist"
      out_subdir <- "ConDa_plot"
      score_label <- "Signed significance  (sign \u00d7 \u2212log\u2081\u2080 fisher q)"
    } else {
      message(sprintf("[Go_lollipopPlot] %s: required effect column not found for tool '%s'.", if (is.null(source_name)) "<unknown>" else source_name, tool))
      return(NULL)
    }

    if (!effect_col %in% colnames(df)) {
      message(sprintf("[Go_lollipopPlot] %s: no finite signed scores available after parsing.", if (is.null(source_name)) "<unknown>" else source_name))
      return(NULL)
    }

    feature_id <- if ("ASV" %in% colnames(df)) {
      as.character(df$ASV)
    } else if ("Row.names" %in% colnames(df)) {
      as.character(df$Row.names)
    } else if ("feature_id" %in% colnames(df)) {
      as.character(df$feature_id)
    } else {
      paste0("Feature_", seq_len(nrow(df)))
    }
    feature_label <- build_feature_label(df)
    effect <- suppressWarnings(as.numeric(df[[effect_col]]))
    pval <- if (p_col %in% colnames(df)) suppressWarnings(as.numeric(df[[p_col]])) else rep(NA_real_, nrow(df))
    qval <- if (!is.na(q_col) && q_col %in% colnames(df)) suppressWarnings(as.numeric(df[[q_col]])) else rep(NA_real_, nrow(df))

    base_score <- {
      use_q <- is.finite(qval) & !is.na(qval)
      score_in <- ifelse(use_q, qval, pval)
      -log10(pmax(score_in, 1e-300))
    }

    # condadist: final_score [0,1] used for point size; other tools: NA
    final_score_vals <- if (tool == "condadist" && "final_score" %in% colnames(df)) {
      suppressWarnings(as.numeric(df[["final_score"]]))
    } else {
      rep(NA_real_, nrow(df))
    }

    signed_score <- sign(effect) * base_score
    dir_label <- ifelse(
      !is.finite(effect) | is.na(effect),
      "NS",
      ifelse(effect > 0, smvar, ifelse(effect < 0, basline, "NS"))
    )

    small_n_warn <- {
      bc <- if ("bas.count" %in% colnames(df)) suppressWarnings(as.integer(unique(df$bas.count)[1])) else NA_integer_
      sc <- if ("smvar.count" %in% colnames(df)) suppressWarnings(as.integer(unique(df$smvar.count)[1])) else NA_integer_
      if (!is.na(bc) && !is.na(sc) && (bc <= 5L || sc <= 5L)) {
        "Warning: small sample size (n\u22645)"
      } else ""
    }
    wilcoxon_fallback_note <- if (identical(tool, "wilcoxon") && "notes" %in% colnames(df)) {
      note_vals <- unique(df$notes[!is.na(df$notes) & nzchar(df$notes)])
      if (length(note_vals) > 0) {
        m <- regmatches(note_vals[1], regexpr("wilcoxon fallback \\(([^\\s]+) skipped", note_vals[1]))
        orig <- if (length(m) > 0) sub("wilcoxon fallback \\(", "", sub(" skipped.*", "", m)) else NULL
        if (!is.null(orig) && nzchar(orig)) paste0("Wilcoxon fallback (", orig, " failed)") else "Wilcoxon fallback"
      } else "Wilcoxon fallback"
    } else ""
    final_subtitle <- paste(
      c(wilcoxon_fallback_note, conda_subtitle, build_optional_metric_note(df), small_n_warn)[nzchar(c(wilcoxon_fallback_note, conda_subtitle, build_optional_metric_note(df), small_n_warn))],
      collapse = " | "
    )

    out <- data.frame(
      feature_id = feature_id,
      feature_label = feature_label,
      effect = effect,
      pval = pval,
      qval = qval,
      base_score = base_score,
      signed_score = signed_score,
      final_score = final_score_vals,
      direction = dir_label,
      mvar = mvar,
      basline = basline,
      smvar = smvar,
      out_subdir = out_subdir,
      title_tool = title_tool,
      file_tool = file_tool,
      plot_subtitle = final_subtitle,
      score_label = score_label,
      stringsAsFactors = FALSE
    )

    out <- out[is.finite(out$signed_score) & !is.na(out$signed_score), , drop = FALSE]
    if (nrow(out) == 0) {
      message(sprintf("[Go_lollipopPlot] %s: no features passed p_cutoff <= %.3g.", if (is.null(source_name)) "<unknown>" else source_name, p_cutoff))
      return(NULL)
    }
    out$sig_value <- ifelse(is.finite(out$qval) & !is.na(out$qval), out$qval, out$pval)
    out$sig_type <- ifelse(is.finite(out$qval) & !is.na(out$qval), "Adj. p", "p")
    out <- out[is.finite(out$sig_value) & !is.na(out$sig_value) & out$sig_value <= p_cutoff, , drop = FALSE]
    if (nrow(out) == 0) {
      return(NULL)
    }
    out <- out[order(out$signed_score, decreasing = TRUE), , drop = FALSE]
    out$dirPadj <- ifelse(is.finite(out$qval) & !is.na(out$qval) & out$qval <= p_cutoff, TRUE, FALSE)
    out$feature_id <- factor(out$feature_id, levels = rev(unique(out$feature_id)))
    out
  }

  if (!is.null(dev.list())) dev.off()

  out_root <- file.path(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")))
  if (!dir.exists(out_root)) dir.create(out_root)
  out_pdf <- file.path(sprintf("%s_%s/pdf", project, format(Sys.Date(), "%y%m%d")))
  if (!dir.exists(out_pdf)) dir.create(out_pdf)

  if (!is.null(result)) {
    if (!is.character(result) || !dir.exists(result)) {
      stop("result must be a directory path returned by a Go DA-family tool or ConDA bridge directory.")
    }
    tool_files <- list.files(result, pattern = "\\.csv$", full.names = FALSE)
    if (length(tool_files) == 0) stop(sprintf("No CSV files found in: %s", result))
    file_list <- data.frame(path = result, file = tool_files, stringsAsFactors = FALSE)
  } else if (!is.null(file_path) && !is.null(files)) {
    filenames <- list.files(file_path, pattern = files)
    file_list <- data.frame(path = file_path, file = filenames, stringsAsFactors = FALSE)
  } else {
    tool_dirs <- c("deseq2", "aldex2", "ancom2", "corncob", "maaslin2")
    file_list <- do.call(rbind, lapply(tool_dirs, function(tool) {
      tool_path <- sprintf("%s_%s/table/%s", project, format(Sys.Date(), "%y%m%d"), tool)
      if (!dir.exists(tool_path)) return(NULL)
      tool_files <- list.files(tool_path, pattern = "\\.csv$", full.names = FALSE)
      if (length(tool_files) == 0) return(NULL)
      data.frame(path = tool_path, file = tool_files, stringsAsFactors = FALSE)
    }))
  }

  if (is.null(file_list) || nrow(file_list) == 0) {
    message("No DA result files found.")
    return(invisible(NULL))
  }

  print(file_list)

  plotlist_pw <- list()
  for (fn in seq_len(nrow(file_list))) {
    filename1 <- file.path(file_list$path[fn], file_list$file[fn])
    source_stub <- tools::file_path_sans_ext(basename(filename1))
    df <- utils::read.csv(filename1, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- gsub('^"|"$', "", colnames(df))

    tool <- detect_tool(filename1)
    if (is.na(tool)) {
      next
    }
    print(tool)

    plot_df <- build_plot_df(df, tool, source_name = basename(filename1))
    if (is.null(plot_df) || nrow(plot_df) == 0) {
      next
    }

    if (!is.null(mycols) && length(mycols) >= 2) {
      dircolors <- c(mycols[1], "grey70", mycols[2])
    } else {
      dircolors <- c("#f8766d", "grey70", "#7cae00")
    }
    bas.count <- if ("bas.count" %in% colnames(df)) unique(df$bas.count)[1] else NA
    smvar.count <- if ("smvar.count" %in% colnames(df)) unique(df$smvar.count)[1] else NA
    names(dircolors) <- c(plot_df$basline[1], "NS", plot_df$smvar[1])
    legend.labs <- c(
      if (is.na(bas.count)) plot_df$basline[1] else paste0(plot_df$basline[1], " (n=", bas.count, ")"),
      "NS",
      if (is.na(smvar.count)) plot_df$smvar[1] else paste0(plot_df$smvar[1], " (n=", smvar.count, ")")
    )
    names(legend.labs) <- c(plot_df$basline[1], "NS", plot_df$smvar[1])
    padj_shape <- c(19, 1)
    names(padj_shape) <- c(TRUE, FALSE)
    sig_shape_label <- if (all(plot_df$sig_type == "Adj. p")) {
      sprintf("Adj. p < %.3g", p_cutoff)
    } else {
      sprintf("p < %.3g", p_cutoff)
    }

    out_dir <- file.path(sprintf("%s_%s/pdf/%s", project, format(Sys.Date(), "%y%m%d"), plot_df$out_subdir[1]))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    table_name_token <- NULL
    if ("name_token" %in% colnames(df)) {
      table_name_vals <- unique(as.character(df$name_token))
      table_name_vals <- table_name_vals[!is.na(table_name_vals) & table_name_vals != "" & table_name_vals != "NA"]
      if (length(table_name_vals) > 0) {
        table_name_token <- gsub("^\\(|\\)$", "", table_name_vals[1])
      }
    }
    call_name_token <- if (is.null(name)) NULL else gsub("^\\(|\\)$", "", as.character(name))

    comparison_token <- NULL
    if ("comparison_token" %in% colnames(df)) {
      vals <- unique(as.character(df$comparison_token))
      vals <- vals[!is.na(vals) & vals != "" & vals != "NA"]
      if (length(vals) > 0) comparison_token <- gsub("^\\(|\\)$", "", vals[1])
    }
    if (is.null(comparison_token)) {
      name_token <- if (!is.null(call_name_token)) call_name_token else table_name_token
      comparison_token <- sprintf("%s.vs.%s%s",
                                  plot_df$basline[1],
                                  plot_df$smvar[1],
                                  if (is.null(name_token)) "" else paste(".", name_token, sep = ""))
    } else if (!is.null(call_name_token) &&
               !(comparison_token == call_name_token || endsWith(comparison_token, paste0(".", call_name_token)))) {
      comparison_token <- paste(comparison_token, call_name_token, sep = ".")
    }

    use_fscore_size <- any(!is.na(plot_df$final_score) & is.finite(plot_df$final_score))

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = signed_score, y = feature_id, color = direction)
    ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = signed_score, y = feature_id, yend = feature_id),
        linewidth = 0.7,
        alpha = 0.8
      ) +
      if (use_fscore_size) {
        ggplot2::geom_point(ggplot2::aes(shape = dirPadj, size = final_score))
      } else {
        ggplot2::geom_point(ggplot2::aes(shape = dirPadj), size = 3)
      } +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7, color = "grey40") +
      ggplot2::scale_color_manual(values = dircolors, labels = legend.labs, drop = FALSE) +
      ggplot2::scale_shape_manual(values = padj_shape, drop = FALSE) +
      if (use_fscore_size) {
        ggplot2::scale_size_continuous(name = "final_score", range = c(1.5, 6), limits = c(0, 1))
      } else NULL +
      ggplot2::scale_y_discrete(labels = stats::setNames(as.character(plot_df$feature_label), as.character(plot_df$feature_id))) +
      ggplot2::labs(
        title = sprintf("%s, %s (%s)", plot_df$mvar[1], plot_df$title_tool[1], sig_shape_label),
        subtitle = if (nzchar(plot_df$plot_subtitle[1])) plot_df$plot_subtitle[1] else NULL,
        x = plot_df$score_label[1],
        y = NULL,
        color = NULL,
        shape = sig_shape_label
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = font + 6),
        plot.title = ggplot2::element_text(size = font + 8),
        plot.subtitle = ggplot2::element_text(size = font + 5, color = "grey40"),
        axis.text.y = ggplot2::element_text(face = "italic"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.justification = c(0, 0),
        legend.box.just = "left",
        legend.box = "vertical",
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, -45),
        legend.key = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 2)
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(nrow = 1, byrow = TRUE, title.position = "top"),
        shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE, title.position = "left"),
        size  = if (use_fscore_size) ggplot2::guide_legend(nrow = 1, byrow = TRUE, title.position = "left") else NULL
      )

    pdf_file <- sprintf(
      "%s/%s.lollipop.%s.(%s).%s.(p=%s).%s.pdf",
      out_dir,
      plot_df$file_tool[1],
      plot_df$mvar[1],
      comparison_token,
      project,
      format(p_cutoff, scientific = FALSE, trim = TRUE),
      format(Sys.Date(), "%y%m%d")
    )

    max_lbl <- max(nchar(as.character(plot_df$feature_label)), na.rm = TRUE)
    plot_height <- max(3.8, min(8.5, 1.3 + 0.12 * nrow(plot_df) + 0.005 * max_lbl))
    plot_grob <- fix_panel_width_grob(p, width)
    grob_width <- tryCatch(grid::convertWidth(sum(plot_grob$widths), "in", valueOnly = TRUE), error = function(e) NA_real_)
    if (!is.finite(grob_width)) {
      grob_width <- width + max(1.1, min(5.5, calc_label_width_in(plot_df$feature_label, font + 6) + 0.25))
    }
    axis_label_width <- max(0, grob_width - width)
    message(sprintf("[Go_lollipopPlot] %s: %d features retained, panel width = %.2f, axis width = %.2f, total width = %.2f, auto height = %.2f",
                    plot_df$file_tool[1], nrow(plot_df), width, axis_label_width, grob_width, plot_height))

    plot_key <- gsub("[^A-Za-z0-9._-]+", "_",
                     if (!is.null(comparison_token) && nzchar(comparison_token)) {
                       paste(plot_df$file_tool[1], comparison_token, sep = "__")
                     } else {
                       source_stub
                     })
    plotlist_pw[[plot_key]] <- p
    if (!isTRUE(patchwork)) {
      ggplot2::ggsave(filename = pdf_file, plot = plot_grob, width = grob_width, height = plot_height, device = "pdf")
      print(normalizePath(pdf_file, winslash = "/", mustWork = FALSE))
    }
  }

  if (isTRUE(patchwork)) return(invisible(plotlist_pw))
  invisible(NULL)
}
