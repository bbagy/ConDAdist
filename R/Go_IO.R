Go_RenameGroupColumns <- function(df, group_1, group_2) {
  g1 <- gsub("[^A-Za-z0-9_]", "_", group_1)
  g2 <- gsub("[^A-Za-z0-9_]", "_", group_2)
  rename_map <- c(
    mean_group1        = paste0("mean_", g1),
    mean_group2        = paste0("mean_", g2),
    prevalence_group1  = paste0("prevalence_", g1),
    prevalence_group2  = paste0("prevalence_", g2)
  )
  for (old in names(rename_map)) {
    if (old %in% colnames(df)) colnames(df)[colnames(df) == old] <- rename_map[[old]]
  }
  if ("direction" %in% colnames(df)) {
    df$direction <- gsub("up_in_group1", paste0("up_in_", g1), df$direction, fixed = TRUE)
    df$direction <- gsub("up_in_group2", paste0("up_in_", g2), df$direction, fixed = TRUE)
  }
  df
}

#' Export standard output tables
#' @export
Go_ExportResults <- function(output_dir, filtered_feature_table, standardized_da,
                             da_consensus, beta_summary,
                             beta_feature_contribution,
                             final_scores, file_prefix = NULL) {
  prefix <- if (is.null(file_prefix) || !nzchar(file_prefix)) "" else paste0(file_prefix, ".")
  files <- list(
    filtered_feature_table = file.path(output_dir, paste0(prefix, "filtered_feature_table.csv")),
    all_methods_standardized = file.path(output_dir, paste0(prefix, "all_methods_standardized.csv")),
    feature_consensus_summary = file.path(output_dir, paste0(prefix, "feature_consensus_summary.csv")),
    beta_summary = file.path(output_dir, paste0(prefix, "beta_summary.csv")),
    beta_feature_contribution = file.path(output_dir, paste0(prefix, "beta_feature_contribution.csv")),
    final_consensus_scores = file.path(output_dir, paste0(prefix, "final_consensus_scores.csv"))
  )

  filtered_feature_table_export <- as.data.frame(filtered_feature_table, stringsAsFactors = FALSE)
  filtered_feature_table_export$feature_id <- rownames(filtered_feature_table_export)
  filtered_feature_table_export <- filtered_feature_table_export[, c("feature_id", setdiff(colnames(filtered_feature_table_export), "feature_id")), drop = FALSE]

  utils::write.csv(filtered_feature_table_export, files$filtered_feature_table, row.names = FALSE)
  utils::write.csv(standardized_da, files$all_methods_standardized, row.names = FALSE)
  utils::write.csv(da_consensus, files$feature_consensus_summary, row.names = FALSE)
  utils::write.csv(beta_summary, files$beta_summary, row.names = FALSE)
  utils::write.csv(beta_feature_contribution, files$beta_feature_contribution, row.names = FALSE)
  utils::write.csv(final_scores, files$final_consensus_scores, row.names = FALSE)

  files
}

Go_ExportVolcanoBridge <- function(output_dir, da_table, final_scores,
                                   filtered_metadata, group_var, group_1, group_2,
                                   analysis_mode = NA_character_,
                                   methods,
                                   distances = NULL,
                                   single_output_dir = NULL,
                                   write_consensus = TRUE,
                                   file_prefix = NULL,
                                   name = NULL) {
  bridge_dir <- output_dir
  single_bridge_dir <- if (!is.null(single_output_dir)) single_output_dir else output_dir
  if (isTRUE(write_consensus)) {
    dir.create(bridge_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    dir.create(single_bridge_dir, recursive = TRUE, showWarnings = FALSE)
  }

  method_signature <- Go_MethodSignature(methods)
  bas.count <- sum(filtered_metadata[[group_var]] == group_1, na.rm = TRUE)
  smvar.count <- sum(filtered_metadata[[group_var]] == group_2, na.rm = TRUE)
  comparison_token <- paste0(group_1, ".vs.", group_2, if (is.null(name) || !nzchar(name)) "" else paste0(".", name))
  comparison_stub <- paste0("(", comparison_token, ")")

  files <- list()

  unique_methods <- unique(da_table$method[!is.na(da_table$method)])
  if (!isTRUE(write_consensus) && length(unique_methods) == 1) {
    method <- unique_methods[1]
    x <- da_table[da_table$method == method, , drop = FALSE]
    if (all(is.na(x$p_value)) && all(is.na(x$effect_size))) {
      message("[ConDA] Volcano output skipped for ", method, ": no native result was available.")
      return(list(dir = single_bridge_dir, files = files))
    }
    if (!"ASV" %in% colnames(x) || all(is.na(x$ASV) | !nzchar(x$ASV))) {
      x$ASV <- x$feature_id
    }
    x$basline <- group_1
    x$smvar <- group_2
    x$bas.count <- bas.count
    x$smvar.count <- smvar.count
    x$mvar <- group_var
    x$name_token <- if (is.null(name)) NA_character_ else as.character(name)
    x$comparison_token <- comparison_token
    x <- Go_RenameGroupColumns(x, group_1, group_2)

    # Fill missing q-values for tested features so plots retain their rows.
    # DESeq2 independent filtering can leave q-values undefined.
    .fill_q <- function(pv, qv) {
      na_q <- is.na(qv) & !is.na(pv)
      if (any(na_q)) {
        qv[na_q] <- stats::p.adjust(pv[na_q], method = "BH")
      }
      qv
    }

    bridge <- switch(
      method,
      deseq2 = {
        padj_bridge <- .fill_q(x$p_value, x$q_value)
        transform(
          x,
          log2FoldChange = coef,
          pvalue = p_value,
          padj = padj_bridge,
          deseq2.P = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
          deseq2.FDR = ifelse(padj_bridge < 0.05, ifelse(coef >= 0, "up", "down"), "NS")
        )
      },
      aldex2 = {
        q_bridge <- .fill_q(x$p_value, x$q_value)
        transform(x, diff.btw = coef, wi.ep = p_value, wi.eBH = q_bridge,
          aldex2.P   = ifelse(p_value < 0.05, ifelse(effect_size >= 0, "up", "down"), "NS"),
          aldex2.FDR = ifelse(q_bridge < 0.05, ifelse(effect_size >= 0, "up", "down"), "NS"))
      },
      ancombc2 = {
        q_bridge <- .fill_q(x$p_value, x$q_value)
        transform(x, lfc_ancombc = coef, pvalue_ancombc = p_value, qvalue_ancombc = q_bridge,
          ancom2.P   = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
          ancom2.FDR = ifelse(q_bridge < 0.05, ifelse(coef >= 0, "up", "down"), "NS"))
      },
      corncob = ,
      corncob_wald = ,
      corncob_lrt = {
        q_bridge <- .fill_q(x$p_value, x$q_value)
        transform(x, corncob_coef = coef, corncob_pvalue = p_value, corncob_qvalue = q_bridge,
          corncob.P   = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
          corncob.FDR = ifelse(q_bridge < 0.05, ifelse(coef >= 0, "up", "down"), "NS"))
      },
      wilcoxon = {
        q_bridge <- .fill_q(x$p_value, x$q_value)
        transform(x, log2FoldChange = coef, pvalue = p_value, padj = q_bridge,
          wilcoxon.P   = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
          wilcoxon.FDR = ifelse(q_bridge < 0.05, ifelse(coef >= 0, "up", "down"), "NS"))
      },
      x
    )

    tool_stub <- if (method == "ancombc2") "ancom2" else method
    file <- file.path(single_bridge_dir, paste0(tool_stub, ".", comparison_stub, ".csv"))
    utils::write.csv(bridge, file, row.names = FALSE)
    files[[tool_stub]] <- file
  }

  if (isTRUE(write_consensus)) {
    consensus <- final_scores
    if (!"ASV" %in% colnames(consensus) || all(is.na(consensus$ASV) | !nzchar(consensus$ASV))) {
      consensus$ASV <- consensus$feature_id
    }
    consensus$baseline <- group_1
    consensus$smvar <- group_2
    consensus$bas.count <- bas.count
    consensus$smvar.count <- smvar.count
    consensus$mvar <- group_var
    consensus$name_token <- if (is.null(name)) NA_character_ else as.character(name)
    consensus$comparison_token <- comparison_token
    consensus <- Go_RenameGroupColumns(consensus, group_1, group_2)
    consensus_p <- if ("combined_p" %in% colnames(consensus)) consensus$combined_p else consensus$cauchy_combined_p
    consensus_q <- if ("combined_q" %in% colnames(consensus)) consensus$combined_q else consensus$cauchy_combined_q
    # V1_JSD skeleton uses median_effect_size for up/down; V2_JSD uses combined_effect_rank.
    skel_vec <- if ("consensus_skeleton" %in% colnames(consensus)) unique(consensus$consensus_skeleton) else "v2"
    skel_vec <- skel_vec[!is.na(skel_vec)]
    use_v1 <- length(skel_vec) == 1 && identical(skel_vec, "v1") &&
              "median_effect_size" %in% colnames(consensus)
    direction_up <- if (use_v1) {
      consensus$median_effect_size >= 0
    } else {
      consensus$combined_effect_rank >= 0.5
    }
    consensus$condadist.P <- ifelse(
      consensus_p < 0.05,
      ifelse(direction_up, "up", "down"),
      "NS"
    )
    consensus$condadist.FDR <- ifelse(
      consensus_q < 0.05,
      ifelse(direction_up, "up", "down"),
      "NS"
    )
    dist_signature <- if (!is.null(distances) && length(distances) > 0)
      paste0(".", paste(distances, collapse = ".")) else ""
    consensus_file <- file.path(bridge_dir, paste0("condadist.", method_signature, dist_signature, ".", comparison_stub, ".csv"))
    utils::write.csv(consensus, consensus_file, row.names = FALSE)
    files$condadist <- consensus_file
  }

  list(dir = if (isTRUE(write_consensus)) bridge_dir else single_bridge_dir,
       files = files, analysis_mode = analysis_mode)
}
