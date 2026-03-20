#' Export standard output tables
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
                                   file_prefix = NULL,
                                   name = NULL) {
  bridge_dir <- output_dir
  dir.create(bridge_dir, recursive = TRUE, showWarnings = FALSE)

  method_signature <- Go_MethodSignature(methods)
  bas.count <- sum(filtered_metadata[[group_var]] == group_1, na.rm = TRUE)
  smvar.count <- sum(filtered_metadata[[group_var]] == group_2, na.rm = TRUE)
  comparison_token <- paste0(group_1, ".vs.", group_2, if (is.null(name) || !nzchar(name)) "" else paste0(".", name))
  comparison_stub <- paste0("(", comparison_token, ")")

  files <- list()

  unique_methods <- unique(da_table$method[!is.na(da_table$method)])
  if (length(unique_methods) == 1) {
    method <- unique_methods[1]
    x <- da_table[da_table$method == method, , drop = FALSE]
    if (!"ASV" %in% colnames(x) || all(is.na(x$ASV) | !nzchar(x$ASV))) {
      x$ASV <- x$feature_id
    }
    x$baseline <- group_1
    x$smvar <- group_2
    x$bas.count <- bas.count
    x$smvar.count <- smvar.count
    x$mvar <- group_var
    x$name_token <- if (is.null(name)) NA_character_ else as.character(name)
    x$comparison_token <- comparison_token

    bridge <- switch(
      method,
      deseq2 = transform(
        x,
        log2FoldChange = coef,
        pvalue = p_value,
        padj = q_value,
        deseq2.P = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
        deseq2.FDR = ifelse(q_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS")
      ),
      aldex2 = transform(
        x,
        diff.btw = coef,
        wi.ep = p_value,
        wi.eBH = q_value,
        aldex2.P = ifelse(p_value < 0.05, ifelse(effect_size >= 0, "up", "down"), "NS"),
        aldex2.FDR = ifelse(q_value < 0.05, ifelse(effect_size >= 0, "up", "down"), "NS")
      ),
      maaslin2 = transform(
        x,
        maaslin2_coef = coef,
        maaslin2_pvalue = p_value,
        maaslin2_qvalue = q_value,
        maaslin2.P = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
        maaslin2.FDR = ifelse(q_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS")
      ),
      ancombc2 = transform(
        x,
        lfc_ancombc = coef,
        pvalue_ancombc = p_value,
        qvalue_ancombc = q_value,
        ancom2.P = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
        ancom2.FDR = ifelse(q_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS")
      ),
      corncob = transform(
        x,
        corncob_coef = coef,
        corncob_pvalue = p_value,
        corncob_qvalue = q_value,
        corncob.P = ifelse(p_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS"),
        corncob.FDR = ifelse(q_value < 0.05, ifelse(coef >= 0, "up", "down"), "NS")
      ),
      x
    )

    tool_stub <- if (method == "ancombc2") "ancom2" else method
    file <- file.path(bridge_dir, paste0(tool_stub, ".", comparison_stub, ".volcano_bridge.csv"))
    utils::write.csv(bridge, file, row.names = FALSE)
    files[[tool_stub]] <- file
  }

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
  consensus$condadist.P <- ifelse(
    consensus$fisher_combined_p < 0.05,
    ifelse(consensus$median_effect_size >= 0, "up", "down"),
    "NS"
  )
  consensus$condadist.FDR <- ifelse(
    consensus$fisher_combined_q < 0.05,
    ifelse(consensus$median_effect_size >= 0, "up", "down"),
    "NS"
  )
  consensus_file <- file.path(bridge_dir, paste0("condadist.", method_signature, ".", comparison_stub, ".volcano_bridge.csv"))
  utils::write.csv(consensus, consensus_file, row.names = FALSE)
  files$condadist <- consensus_file

  list(dir = bridge_dir, files = files, analysis_mode = analysis_mode)
}
