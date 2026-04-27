#' Generate QC plots for ConDA-dist results
Go_ConDaQCplot <- function(result,
                           output_dir,
                           x_method = c("deseq2", "median"),
                           label_col = c("auto", "TaxaName", "ShortName", "feature_id"),
                           top_n = 20,
                           width = 6,
                           height = 4.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }

  x_method <- match.arg(x_method)
  label_col <- match.arg(label_col)
  if (is.null(output_dir) || !nzchar(output_dir)) {
    stop("`output_dir` must be provided.")
  }

  final_scores <- result$final_scores
  da_table <- result$da_standardized$all_methods_standardized
  plot_taxonomy <- result$input_bundle$taxonomy %||% NULL
  final_scores <- Go_JoinPlotTaxonomy(final_scores, da_table, taxonomy = plot_taxonomy)
  prefix <- result$file_prefix %||% "ConDA-dist"
  plot_dir <- output_dir
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  volcano_x_col <- if (x_method == "median" || !paste0(x_method, "_effect_size") %in% colnames(final_scores)) {
    "combined_effect_rank"
  } else {
    paste0(x_method, "_effect_size")
  }
  final_scores$plot_label <- Go_SelectPlotLabel(final_scores, label_col = label_col)
  final_scores$plot_label[is.na(final_scores$plot_label) | !nzchar(final_scores$plot_label)] <- final_scores$feature_id[is.na(final_scores$plot_label) | !nzchar(final_scores$plot_label)]
  # feature_id를 breaks, plot_label을 display용 named vector로 유지
  # plot_label_unique는 volcano/scatter용 텍스트 라벨에만 사용
  final_scores$plot_label_unique <- Go_MakeUniqueLabels(final_scores$plot_label, final_scores$feature_id)
  q_col <- if ("combined_q" %in% colnames(final_scores)) "combined_q" else "cauchy_combined_q"
  p_col <- if ("combined_p" %in% colnames(final_scores)) "combined_p" else "cauchy_combined_p"
  final_scores$consensus_y <- -log10(pmax(final_scores[[q_col]], 1e-300))
  consensus_subtitle <- Go_BuildConsensusSupportSubtitle(final_scores, q_col = q_col)
  # High-contrast, colorblind-friendly palette for QC readability.
  class_cols <- c(
    Core_consensus = "#0072B2",
    Structure_driver = "#D55E00",
    Local_DA = "#009E73",
    Weak_signal = "#999999"
  )
  class_labels <- c(
    Core_consensus = "Core consensus",
    Structure_driver = "Community structure-associated",
    Local_DA = "DA-supported local signal",
    Weak_signal = "Weak signal"
  )
  class_note <- paste(
    "Core consensus: strong DA evidence with beta-structure support.",
    "Community structure-associated: comparatively stronger contribution to between-group community separation.",
    "DA-supported local signal: DA signal with more limited structure contribution.",
    "Weak signal: lower support across evidence layers."
  )
  overlap_cols <- c(`TRUE` = "#CC79A7", `FALSE` = "#D9D9D9")
  has_plotly <- requireNamespace("plotly", quietly = TRUE) &&
    requireNamespace("htmlwidgets", quietly = TRUE) &&
    requireNamespace("htmltools", quietly = TRUE)
  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#ece7df", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#f7f4ef", color = "#d9d3ca"),
      legend.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  p_volcano <- ggplot2::ggplot(
    final_scores,
    ggplot2::aes(
      x = .data[[volcano_x_col]],
      y = consensus_y,
      color = classification,
      size = beta_contribution_score,
      text = paste0(
        "Taxon: ", plot_label,
        "<br>Class: ", classification,
        "<br>Combined q: ", signif(.data[[q_col]], 3),
        "<br>Effect: ", signif(.data[[volcano_x_col]], 3),
        "<br>Beta contribution: ", signif(beta_contribution_score, 3)
      )
    )
  ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_manual(values = class_cols, labels = class_labels, drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(1.5, 4)) +
    ggplot2::labs(
      title = "ConDAdist",
      subtitle = paste("Consensus Volcano", consensus_subtitle, sep = " | "),
      x = volcano_x_col,
      y = "-log10(Combined q)",
      color = "Consensus class",
      size = "Beta contribution",
      caption = paste(
        "Higher points show stronger combined DA evidence across methods.",
        class_note
      )
    ) +
    base_theme

  p_scatter <- ggplot2::ggplot(
    final_scores,
    ggplot2::aes(
      x = cauchy_da_score,
      y = beta_contribution_score,
      color = classification,
      text = paste0(
        "Taxon: ", plot_label,
        "<br>Class: ", classification,
        "<br>Combined DA evidence: ", signif(cauchy_da_score, 3),
        "<br>Community structure contribution: ", signif(beta_contribution_score, 3),
        "<br>Combined q: ", signif(.data[[q_col]], 3)
      )
    )
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::scale_color_manual(values = class_cols, labels = class_labels, drop = FALSE) +
    ggplot2::labs(
      title = "ConDAdist",
      subtitle = "DA Evidence vs Community Structure Contribution",
      x = "Combined DA evidence",
      y = "Community structure contribution",
      color = "Consensus class",
      caption = paste(
        "Upper-right taxa show both stronger DA evidence and larger contribution to between-group beta-diversity separation.",
        class_note
      )
    ) +
    base_theme

  # feature_id → display label 매핑 (breaks = feature_id, labels = taxonomy name)
  label_map <- stats::setNames(final_scores$plot_label, final_scores$feature_id)

  overlap_df <- Go_BuildMethodOverlapLong(final_scores, da_table)
  p_overlap <- ggplot2::ggplot(
    overlap_df,
    ggplot2::aes(
      x = feature_id,
      y = method,
      fill = is_significant,
      text = paste0(
        "Taxon: ", plot_label,
        "<br>Method: ", method,
        "<br>Detected: ", is_significant
      )
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = overlap_cols) +
    ggplot2::scale_x_discrete(labels = label_map) +
    ggplot2::labs(
      title = "ConDAdist",
      subtitle = "Method Overlap",
      x = NULL,
      y = "Method",
      fill = "Detected",
      caption = "Filled cells indicate taxa detected as significant by each DA method. Use the HTML version to inspect taxa by hover."
    ) +
    base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  top_df <- final_scores[order(-final_scores$priority_score), , drop = FALSE]
  top_df <- utils::head(top_df, top_n)
  top_label_map <- stats::setNames(top_df$plot_label, top_df$feature_id)
  # feature_id를 factor로 사용 → 중복 없음, 에러 없음
  top_df$feature_id_fct <- factor(top_df$feature_id, levels = rev(top_df$feature_id))
  p_bar <- ggplot2::ggplot(
    top_df,
    ggplot2::aes(
      x = feature_id_fct,
      y = priority_score,
      fill = classification,
      text = paste0(
        "Taxon: ", plot_label,
        "<br>Class: ", classification,
        "<br>Priority score: ", signif(priority_score, 3),
        "<br>Combined q: ", signif(.data[[q_col]], 3),
        "<br>Beta contribution: ", signif(beta_contribution_score, 3)
      )
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = class_cols, labels = class_labels, drop = FALSE) +
    ggplot2::scale_x_discrete(labels = top_label_map) +
    ggplot2::labs(
      title = "ConDAdist",
      subtitle = paste(paste("Top", top_n, "Consensus Taxa"), consensus_subtitle, sep = " | "),
      x = NULL,
      y = "Priority score",
      fill = "Consensus class",
      caption = paste(
        "Priority score is a ranking index that combines DA evidence and beta-structure contribution.",
        class_note
      )
    ) +
    base_theme

  plots <- list(
    consensus_volcano = p_volcano,
    da_vs_beta_scatter = p_scatter,
    method_overlap_bubble = p_overlap,
    top_taxa_barplot = p_bar
  )

  file_paths <- c(
    consensus_volcano = file.path(plot_dir, paste0(prefix, ".consensus_volcano.pdf")),
    da_vs_beta_scatter = file.path(plot_dir, paste0(prefix, ".da_vs_beta_scatter.pdf")),
    method_overlap_bubble = file.path(plot_dir, paste0(prefix, ".method_overlap_bubble.pdf")),
    top_taxa_barplot = file.path(plot_dir, paste0(prefix, ".top_taxa_barplot.pdf"))
  )
  summary_html <- file.path(plot_dir, paste0(prefix, ".qc_summary_4panel.html"))

  for (nm in names(plots)) {
    ggplot2::ggsave(filename = file_paths[[nm]], plot = plots[[nm]], width = width, height = height, device = "pdf")
  }
  message("[ConDA] QC PDFs saved to: ", plot_dir)

  html_files <- NULL
  if (has_plotly) {
    saved_html <- tryCatch(
      Go_SaveQCPanelHTML(
        plots = plots,
        file = summary_html,
        plot_width = width * 100,
        plot_height = height * 100
      ),
      error = function(e) {
        message("[ConDA] QC HTML generation failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(saved_html) && file.exists(saved_html)) {
      html_files <- c(qc_summary_4panel = saved_html)
      message("[ConDA] QC HTML saved to: ", saved_html)
    }
  }

  invisible(list(
    plots = plots,
    files = file_paths,
    html_files = html_files
  ))
}

Go_BuildConsensusSupportSubtitle <- function(final_scores, q_col) {
  if (is.null(final_scores) || nrow(final_scores) == 0) {
    return("no features")
  }

  q_vals <- suppressWarnings(as.numeric(final_scores[[q_col]]))
  sig <- is.finite(q_vals) & !is.na(q_vals) & q_vals < 0.05
  core_n <- if ("classification" %in% colnames(final_scores)) {
    sum(final_scores$classification %in% "Core_consensus", na.rm = TRUE)
  } else {
    NA_integer_
  }
  med_methods <- if ("n_methods_significant" %in% colnames(final_scores) && any(sig)) {
    stats::median(final_scores$n_methods_significant[sig], na.rm = TRUE)
  } else {
    NA_real_
  }
  med_informative <- if ("n_informative_p" %in% colnames(final_scores) && any(sig)) {
    stats::median(final_scores$n_informative_p[sig], na.rm = TRUE)
  } else {
    NA_real_
  }

  parts <- c(
    sprintf("sig=%d/%d", sum(sig, na.rm = TRUE), nrow(final_scores)),
    if (is.finite(core_n)) sprintf("core=%d", core_n) else NULL,
    if (is.finite(med_methods)) sprintf("median method support=%.1f", med_methods) else NULL,
    if (is.finite(med_informative)) sprintf("median informative p=%.1f", med_informative) else NULL
  )
  paste(parts, collapse = " | ")
}

Go_SaveQCPanelHTML <- function(plots, file, plot_width = 900, plot_height = 600) {
  if (!requireNamespace("plotly", quietly = TRUE) ||
      !requireNamespace("htmltools", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    return(invisible(NULL))
  }

  widgets <- lapply(plots, function(p) {
    plotly::ggplotly(p, tooltip = "text", width = plot_width, height = plot_height)
  })
  plot_titles <- c(
    consensus_volcano = "Consensus Volcano",
    da_vs_beta_scatter = "DA vs Beta",
    method_overlap_bubble = "Method Overlap",
    top_taxa_barplot = "Top Taxa"
  )
  plot_notes <- c(
    consensus_volcano = "Higher points show stronger combined DA evidence across methods.",
    da_vs_beta_scatter = "Upper-right taxa show both stronger DA evidence and larger contribution to between-group beta-diversity separation.",
    method_overlap_bubble = "Filled cells indicate taxa detected as significant by each DA method. Use hover to inspect taxa.",
    top_taxa_barplot = "Priority score is a ranking index that combines DA evidence and beta-structure contribution."
  )
  widget_names <- names(widgets)
  tab_ids <- paste0("conda_tab_", seq_along(widgets))
  doc <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title("ConDAdist QC Summary"),
      htmltools::tags$style(htmltools::HTML("
        body { font-family: Arial, sans-serif; margin: 10px; }
        .conda-header { margin-bottom: 14px; }
        .conda-title { font-size: 24px; font-weight: 700; color: #12344d; margin-bottom: 4px; }
        .conda-subtitle { font-size: 13px; color: #5b6870; }
        .conda-tabs { display: flex; gap: 8px; margin-bottom: 12px; flex-wrap: wrap; }
        .conda-tab-btn {
          border: 1px solid #cfcfcf; background: #f7f7f7; color: #222;
          padding: 8px 12px; cursor: pointer; border-radius: 6px; font-size: 13px;
        }
        .conda-tab-btn.active { background: #0072B2; color: white; border-color: #0072B2; }
        .conda-panel { display: none; border: 1px solid #ddd; padding: 8px; background: #fff; width: fit-content; }
        .conda-panel.active { display: block; }
        .conda-panel-title { font-size: 18px; font-weight: 600; color: #12344d; margin: 2px 0 6px 0; }
        .conda-panel-note { font-size: 12px; color: #58636b; margin-bottom: 8px; max-width: 900px; }
      "))
    ),
    htmltools::tags$body(
      htmltools::tags$script(htmltools::HTML("
        function condaOpenTab(id, btn) {
          var panels = document.getElementsByClassName('conda-panel');
          for (var i = 0; i < panels.length; i++) { panels[i].classList.remove('active'); }
          var buttons = document.getElementsByClassName('conda-tab-btn');
          for (var j = 0; j < buttons.length; j++) { buttons[j].classList.remove('active'); }
          document.getElementById(id).classList.add('active');
          btn.classList.add('active');
        }
      ")),
      htmltools::tags$div(
        class = "conda-header",
        htmltools::tags$div(class = "conda-title", "ConDAdist"),
        htmltools::tags$div(
          class = "conda-subtitle",
          "QC summary of consensus DA evidence, community structure contribution, method overlap, and top-ranked taxa."
        )
      ),
      htmltools::tags$div(
        class = "conda-tabs",
        lapply(seq_along(widgets), function(i) {
          htmltools::tags$button(
            class = if (i == 1) "conda-tab-btn active" else "conda-tab-btn",
            onclick = sprintf("condaOpenTab('%s', this)", tab_ids[[i]]),
            plot_titles[[widget_names[[i]]]] %||% widget_names[[i]]
          )
        })
      ),
      htmltools::tagList(
        lapply(seq_along(widgets), function(i) {
          htmltools::tags$div(
            id = tab_ids[[i]],
            class = if (i == 1) "conda-panel active" else "conda-panel",
            htmltools::tags$div(class = "conda-panel-title", plot_titles[[widget_names[[i]]]] %||% widget_names[[i]]),
            htmltools::tags$div(class = "conda-panel-note", plot_notes[[widget_names[[i]]]] %||% ""),
            widgets[[i]]
          )
        })
      )
    )
  )
  # Step 1: save HTML + lib/ directly to the output path
  # (htmltools::save_html writes plotly.js into libdir = "lib/" next to the file)
  htmltools::save_html(doc, file = file)

  if (!file.exists(file)) {
    message("[ConDA] HTML output could not be saved.")
    return(invisible(NULL))
  }

  # Step 2: try pandoc to embed lib/ resources inline (single self-contained file)
  # If it succeeds, the lib/ dir becomes unnecessary but is left in place.
  # If it fails, the HTML + lib/ pair already works when opened locally.
  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    pandoc_info <- tryCatch(rmarkdown::find_pandoc(), error = function(e) NULL)
    pandoc_ver  <- tryCatch(numeric_version(pandoc_info$version), error = function(e) NULL)
    pandoc_opts <- if (!is.null(pandoc_ver) && pandoc_ver >= numeric_version("3.0")) {
      c("--embed-resources", "--standalone")
    } else {
      "--self-contained"
    }
    tmp_embedded <- paste0(file, ".tmp_embed.html")
    embedded_ok <- tryCatch({
      rmarkdown::pandoc_convert(input = file, output = tmp_embedded, options = pandoc_opts)
      file.exists(tmp_embedded) && file.size(tmp_embedded) > 1000L
    }, error = function(e) {
      message("[ConDA] pandoc embedding failed (", conditionMessage(e),
              "); HTML saved with companion lib/ directory.")
      FALSE
    })
    if (embedded_ok) {
      file.rename(tmp_embedded, file)   # replace with fully embedded version
    } else {
      unlink(tmp_embedded)
    }
  }

  invisible(file)
}

Go_BuildMethodOverlapLong <- function(final_scores, da_table) {
  methods <- unique(da_table$method[!is.na(da_table$method)])
  n <- nrow(final_scores)
  out <- lapply(methods, function(method) {
    sig_col <- paste0(method, "_is_significant")
    is_sig  <- final_scores[[sig_col]] %||% rep(FALSE, n)
    if (length(is_sig) != n) is_sig <- rep(FALSE, n)
    data.frame(
      feature_id        = final_scores$feature_id,
      plot_label        = final_scores$plot_label,
      plot_label_unique = final_scores$plot_label_unique,
      method            = method,
      is_significant    = is_sig,
      stringsAsFactors  = FALSE
    )
  })
  do.call(rbind, out)
}

Go_JoinPlotTaxonomy <- function(final_scores, da_table, taxonomy = NULL) {
  if (is.null(final_scores) || nrow(final_scores) == 0) {
    return(final_scores)
  }

  taxonomy_cols <- intersect(
    c("feature_id", "TaxaName", "ShortName", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Rank1"),
    union(colnames(da_table), colnames(taxonomy %||% data.frame()))
  )
  if (length(taxonomy_cols) <= 1) {
    return(final_scores)
  }

  tax_sources <- list()
  if (!is.null(taxonomy) && "feature_id" %in% colnames(taxonomy)) {
    tax_sources[[length(tax_sources) + 1L]] <- taxonomy[, intersect(taxonomy_cols, colnames(taxonomy)), drop = FALSE]
  }
  if (!is.null(da_table) && "feature_id" %in% colnames(da_table)) {
    tax_sources[[length(tax_sources) + 1L]] <- unique(da_table[, intersect(taxonomy_cols, colnames(da_table)), drop = FALSE])
  }
  if (length(tax_sources) == 0) {
    return(final_scores)
  }

  tax_map <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE, suffixes = c("", ".src")), tax_sources)
  out <- merge(final_scores, tax_map, by = "feature_id", all.x = TRUE, suffixes = c("", ".tax"))

  for (col in setdiff(taxonomy_cols, "feature_id")) {
    source_cols <- c(paste0(col, ".tax"), paste0(col, ".src"))
    for (src_col in source_cols) {
      if (!src_col %in% colnames(out)) {
        next
      }
      if (!col %in% colnames(out)) {
        out[[col]] <- out[[src_col]]
      } else {
        missing_idx <- is.na(out[[col]]) | !nzchar(as.character(out[[col]]))
        out[[col]][missing_idx] <- out[[src_col]][missing_idx]
      }
      out[[src_col]] <- NULL
    }
  }
  out
}

Go_SelectPlotLabel <- function(final_scores, label_col = "auto") {
  if (!identical(label_col, "auto")) {
    if (label_col %in% colnames(final_scores)) {
      return(as.character(final_scores[[label_col]]))
    }
    return(as.character(final_scores$feature_id))
  }

  taxonomy_label <- Go_BuildTaxonomyDisplayLabel(final_scores)
  if (any(Go_IsInformativeTaxonomyLabel(taxonomy_label, final_scores$feature_id))) {
    return(taxonomy_label)
  }

  preferred_cols <- c("TaxaName", "ShortName")
  for (col in preferred_cols) {
    if (!col %in% colnames(final_scores)) {
      next
    }
    vals <- as.character(final_scores[[col]])
    valid <- Go_IsInformativeTaxonomyLabel(vals, final_scores$feature_id)
    if (any(valid)) {
      return(vals)
    }
  }
  as.character(final_scores$feature_id)
}

Go_BuildTaxonomyDisplayLabel <- function(final_scores) {
  taxonomy_cols <- intersect(
    c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Rank1"),
    colnames(final_scores)
  )
  if (length(taxonomy_cols) == 0) {
    return(as.character(final_scores$feature_id))
  }

  labels <- apply(final_scores[, taxonomy_cols, drop = FALSE], 1, function(x) {
    vals <- as.character(x)
    vals <- vals[!is.na(vals) & nzchar(trimws(vals)) & vals != "__"]
    vals <- vals[!grepl("^(ASV|OTU|Feature)[._-]?[0-9A-Za-z_-]*$", vals, ignore.case = TRUE)]
    if (length(vals) == 0) {
      return(NA_character_)
    }
    vals[1]
  })
  as.character(labels)
}

Go_IsInformativeTaxonomyLabel <- function(labels, feature_ids) {
  labels <- as.character(labels)
  feature_ids <- as.character(feature_ids)
  valid <- !is.na(labels) & nzchar(trimws(labels))
  valid <- valid & labels != feature_ids
  valid <- valid & !grepl("^(ASV|OTU|Feature)[._-]?[0-9A-Za-z_-]*$", labels, ignore.case = TRUE)
  valid
}

Go_MakeUniqueLabels <- function(labels, feature_ids) {
  out <- as.character(labels)
  short_ids <- ifelse(
    nchar(feature_ids) > 16,
    paste0(substr(feature_ids, 1, 12), ".."),
    feature_ids
  )

  missing_idx <- is.na(out) | !nzchar(trimws(out))
  out[missing_idx] <- short_ids[missing_idx]

  dup_idx <- duplicated(out) | duplicated(out, fromLast = TRUE)
  if (any(dup_idx)) {
    dup_groups <- split(which(dup_idx), out[dup_idx])
    for (idx in dup_groups) {
      out[idx] <- paste0(out[idx], " [", seq_along(idx), "]")
    }
  }

  out
}
