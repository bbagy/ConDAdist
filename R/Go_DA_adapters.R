#' Filter low-information features
#'
#' @param feature_table Numeric feature matrix with features in rows.
#' @param metadata data.frame of aligned metadata.
#' @param prevalence Minimum fraction of samples with non-zero abundance.
#' @param abundance Minimum mean relative abundance threshold.
#' @param output_dir Directory for exported filtered table.
#'
#' @return A list with filtered feature table and metadata.
Go_FilterFeatures <- function(feature_table, metadata, prevalence = 0.1,
                              abundance = 1e-4, output_dir = NULL) {
  feature_table <- Go_AsMatrix(feature_table)
  relative_abundance <- sweep(
    feature_table,
    2,
    colSums(feature_table, na.rm = TRUE),
    "/"
  )
  relative_abundance[!is.finite(relative_abundance)] <- 0

  prevalence_vec <- rowMeans(feature_table > 0, na.rm = TRUE)
  abundance_vec <- rowMeans(relative_abundance, na.rm = TRUE)
  keep <- prevalence_vec >= prevalence & abundance_vec >= abundance

  filtered_table <- feature_table[keep, , drop = FALSE]

  if (!is.null(output_dir)) {
    utils::write.csv(
      filtered_table,
      file = file.path(output_dir, "filtered_feature_table.csv"),
      quote = FALSE
    )
  }

  list(
    feature_table = filtered_table,
    metadata = metadata,
    filter_summary = data.frame(
      n_features_input = nrow(feature_table),
      n_features_retained = nrow(filtered_table),
      prevalence = prevalence,
      abundance = abundance,
      filter_mode = "prevalence_abundance",
      stringsAsFactors = FALSE
    )
  )
}

#' Run all requested DA method adapters
Go_RunDAmethods <- function(feature_table, metadata, group_var, group_1, group_2,
                            random_effects = NULL, covariates = NULL, methods, method_controls = NULL,
                            alpha = 0.05) {
  adapters <- list(
    ancombc2 = Go_DA_ancombc2,
    aldex2 = Go_DA_aldex2,
    maaslin2 = Go_DA_maaslin,
    corncob = Go_DA_corncob,
    deseq2 = Go_DA_deseq2
  )

  results <- vector("list", length(methods))
  names(results) <- methods

  for (method in methods) {
    adapter <- adapters[[method]]
    if (is.null(adapter)) {
      results[[method]] <- Go_CreateAdapterFallback(
        feature_ids = rownames(feature_table),
        method = method,
        note = "Unsupported method."
      )
      next
    }

    results[[method]] <- tryCatch(
      adapter(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        random_effects = random_effects,
        covariates = covariates,
        control = method_controls[[method]] %||% if (identical(method, "maaslin2")) method_controls[["maaslin"]] else NULL,
        alpha = alpha
      ),
      error = function(e) {
        Go_CreateAdapterFallback(
          feature_ids = rownames(feature_table),
          method = method,
          note = conditionMessage(e)
        )
      }
    )
  }

  results
}

#' Standardize method-specific outputs into a shared schema
Go_StandardizeDA <- function(da_results, comparison, alpha = 0.05) {
  standardized <- lapply(names(da_results), function(method) {
    result <- da_results[[method]]
    Go_EnsureStandardDA(
      x = result,
      method = method,
      comparison = comparison,
      alpha = alpha
    )
  })

  standardized <- do.call(rbind, standardized)
  rownames(standardized) <- NULL

  list(
    all_methods_standardized = standardized
  )
}

#' ANCOM-BC2 adapter
Go_DA_ancombc2 <- function(feature_table, metadata, group_var, group_1, group_2,
                           random_effects = NULL,
                           covariates = NULL, control = NULL, alpha = 0.05) {
  control <- Go_GetDAMethodControls("ancombc2", control)
  native_check <- Go_ShouldUseNativeAdapter("ancombc2", feature_table, metadata, group_var, group_1, group_2)
  if (!isTRUE(native_check$ok)) {
    return(Go_BasicEffectAdapter(
      feature_table = feature_table,
      metadata = metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      covariates = covariates,
      method = "ancombc2",
      effect_type = "coef",
      notes = native_check$note,
      control = control,
      alpha = alpha
    ))
  }
  Go_RunNativeAdapter(
    method_name = "ANCOMBC2",
    package_names = c("ANCOMBC", "phyloseq"),
    native_fun = function() {
      prepared_base <- Go_PrepareDAInputs(
        feature_table, metadata, group_var, group_1, group_2,
        covariates = covariates, random_effects = random_effects
      )

      prepare_attempt <- function(prepared, cutoff = NULL, drop_zero_sum = FALSE) {
        attempt <- prepared

        if (isTRUE(drop_zero_sum)) {
          keep_samples <- colSums(attempt$feature_table, na.rm = TRUE) > 0
          if (any(keep_samples)) {
            attempt$feature_table <- attempt$feature_table[, keep_samples, drop = FALSE]
            attempt$metadata <- attempt$metadata[colnames(attempt$feature_table), , drop = FALSE]
          }
        }

        if (!is.null(cutoff) && is.finite(cutoff) && cutoff > 0) {
          legacy_filtered <- Go_FilterLegacyRelativeMean(
            feature_table = attempt$feature_table,
            metadata = attempt$metadata,
            cutoff = cutoff
          )
          attempt$feature_table <- legacy_filtered$feature_table
          attempt$metadata <- legacy_filtered$metadata
        }

        attempt$feature_table <- Go_RemoveZeroVarianceFeatures(attempt$feature_table)
        attempt
      }

      run_attempt <- function(prepared) {
        method_input <- Go_PrepareMethodInput(prepared, "ancombc2", control)
        ps <- method_input$phyloseq
        ANCOMBC::ancombc2(
          data = ps,
          taxa_are_rows = TRUE,
          fix_formula = method_input$full_formula_terms,
          rand_formula = Go_BuildRandomFormula(random_effects),
          group = prepared$temp_group_var,
          prv_cut = control$prv_cut,
          lib_cut = control$lib_cut,
          p_adj_method = control$p_adj_method,
          pseudo = control$pseudo,
          pseudo_sens = control$pseudo_sens,
          struc_zero = control$struc_zero,
          neg_lb = control$neg_lb,
          alpha = alpha,
          verbose = FALSE,
          global = control$global,
          pairwise = control$pairwise,
          dunnet = control$dunnet,
          trend = control$trend
        )
      }

      prepared <- prepare_attempt(prepared_base)
      fit <- tryCatch(
        run_attempt(prepared),
        error = function(e) e
      )
      used_cutoff <- NULL

      if (inherits(fit, "error")) {
        cutoff <- 0.001
        increment <- 0.0005
        final_cutoff <- 0.01

        while (cutoff <= final_cutoff) {
          prepared_retry <- prepare_attempt(
            prepared = prepared_base,
            cutoff = cutoff,
            drop_zero_sum = TRUE
          )
          fit_retry <- tryCatch(
            run_attempt(prepared_retry),
            error = function(e) e
          )
          if (!inherits(fit_retry, "error")) {
            prepared <- prepared_retry
            fit <- fit_retry
            used_cutoff <- cutoff
            break
          }
          cutoff <- cutoff + increment
        }
      }

      if (inherits(fit, "error")) {
        stop(conditionMessage(fit))
      }

      res <- fit$res
      coef_col <- Go_GetANCOMBCCoefColumn(res, "lfc_")
      p_col <- Go_GetANCOMBCCoefColumn(res, "p_")
      q_col <- Go_GetANCOMBCCoefColumn(res, "q_")

      control_note <- Go_ControlNote(control)
      if (!is.null(used_cutoff)) {
        control_note <- paste0(
          control_note,
          ", legacy_filter_cutoff=",
          format(used_cutoff, trim = TRUE, scientific = FALSE)
        )
      }

      out <- Go_CreateBaseDAResult(
        prepared, "ancombc2", "coef",
        Go_CombineNotes("Native ANCOMBC2 adapter", control_note)
      )
      out$ASV <- Go_MakeLegacyASVKey(out$feature_id)
      out <- Go_FillDAResult(
        base_result = out,
        feature_ids = res$taxon,
        coef = res[[coef_col]],
        effect_size = res[[coef_col]],
        p_value = res[[p_col]],
        q_value = res[[q_col]]
      )
      out$is_significant <- out$q_value < alpha
      out
    },
    fallback_fun = function(note) {
      Go_BasicEffectAdapter(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        covariates = covariates,
        method = "ancombc2",
        effect_type = "coef",
        notes = note,
        control = control,
        alpha = alpha
      )
    }
  )
}

#' ALDEx2 adapter
Go_DA_aldex2 <- function(feature_table, metadata, group_var, group_1, group_2,
                         random_effects = NULL,
                         covariates = NULL, control = NULL, alpha = 0.05) {
  control <- Go_GetDAMethodControls("aldex2", control)
  native_check <- Go_ShouldUseNativeAdapter("aldex2", feature_table, metadata, group_var, group_1, group_2)
  if (!isTRUE(native_check$ok)) {
    return(Go_BasicEffectAdapter(
      feature_table = feature_table,
      metadata = metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      covariates = covariates,
      method = "aldex2",
      effect_type = "aldex_effect",
      notes = native_check$note,
      control = control,
      alpha = alpha
    ))
  }
  Go_RunNativeAdapter(
    method_name = "ALDEx2",
    package_names = "ALDEx2",
    native_fun = function() {
      prepared <- Go_PrepareDAInputs(
        feature_table, metadata, group_var, group_1, group_2,
        covariates = covariates, random_effects = random_effects
      )
      method_input <- Go_PrepareMethodInput(prepared, "aldex2", control)
      has_covariates <- length(prepared$covariates) > 0

      if (!has_covariates) {
        conds <- as.character(method_input$condition_vector)
        conds[conds == "ref"] <- "0"
        conds[conds == "cmp"] <- "1"

        native_counts <- Go_PrepareCountMatrix(
          prepared$feature_table,
          min_count = 0,
          zero_replace = FALSE,
          zero_replace_value = 1
        )
        asv_matrix <- t(native_counts)
        if (!is.matrix(asv_matrix)) {
          asv_matrix <- as.matrix(asv_matrix)
        }

        set.seed(control$seed %||% 1L)
        fit_try <- try(ALDEx2::aldex(asv_matrix, conds, test = "t",
                                     mc.samples = control$mc_samples,
                                     denom = control$denom), silent = TRUE)
        if (inherits(fit_try, "try-error")) {
          set.seed(control$seed %||% 1L)
          fit <- ALDEx2::aldex(t(asv_matrix), conds, test = "t",
                               mc.samples = control$mc_samples,
                               denom = control$denom)
        } else {
          fit <- fit_try
        }

        fit_df <- as.data.frame(fit)
        feature_ids <- rownames(fit_df)

        out <- Go_CreateBaseDAResult(
          prepared, "aldex2", "aldex_effect",
          Go_CombineNotes("Native ALDEx2 adapter (t-test)", Go_ControlNote(control))
        )
        out <- Go_FillDAResult(
          base_result = out,
          feature_ids = feature_ids,
          coef = fit_df$diff.btw,
          effect_size = fit_df$effect,
          p_value = fit_df$wi.ep,
          q_value = fit_df$wi.eBH
        )
      } else {
        # GLM mode with covariates: model matrix approach
        native_counts <- Go_PrepareCountMatrix(
          prepared$feature_table,
          min_count = 0,
          zero_replace = control$zero_replace %||% FALSE,
          zero_replace_value = control$zero_replace_value %||% 0.5
        )
        design_terms <- c(prepared$covariates, prepared$temp_group_var)
        design_formula <- stats::as.formula(
          paste("~", paste(vapply(design_terms, Go_BacktickName, character(1)), collapse = " + "))
        )
        mod_matrix <- stats::model.matrix(design_formula, data = prepared$metadata)

        set.seed(control$seed %||% 1L)
        glm_fit <- ALDEx2::aldex(native_counts, mod_matrix, test = "glm",
                                 mc.samples = control$mc_samples,
                                 denom = control$denom,
                                 verbose = FALSE,
                                 useMC = control$use_mc)
        glm_df <- as.data.frame(glm_fit)

        # Extract columns for the group effect (.conda_groupcmp)
        group_col <- paste0(prepared$temp_group_var, "cmp")
        est_col  <- paste0(group_col, "Est")
        pval_col <- paste0(group_col, "pval")
        padj_col <- paste0(group_col, "pval.padj")
        if (!padj_col %in% colnames(glm_df)) {
          padj_col <- paste0(group_col, "pval.holm")
        }

        out <- Go_CreateBaseDAResult(
          prepared, "aldex2", "aldex_effect",
          Go_CombineNotes("Native ALDEx2 adapter (GLM)", Go_ControlNote(control))
        )
        out <- Go_FillDAResult(
          base_result = out,
          feature_ids = rownames(glm_df),
          coef = glm_df[[est_col]],
          effect_size = glm_df[[est_col]],
          p_value = glm_df[[pval_col]],
          q_value = glm_df[[padj_col]]
        )
      }
      out$is_significant <- out$q_value < alpha
      out
    },
    fallback_fun = function(note) {
      Go_BasicEffectAdapter(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        covariates = covariates,
        method = "aldex2",
        effect_type = "aldex_effect",
        notes = note,
        control = control,
        alpha = alpha
      )
    }
  )
}

#' Maaslin2 adapter
Go_DA_maaslin <- function(feature_table, metadata, group_var, group_1, group_2,
                          random_effects = NULL,
                          covariates = NULL, control = NULL, alpha = 0.05) {
  control <- Go_GetDAMethodControls("maaslin2", control)
  native_check <- Go_ShouldUseNativeAdapter("maaslin2", feature_table, metadata, group_var, group_1, group_2)
  if (!isTRUE(native_check$ok)) {
    return(Go_BasicEffectAdapter(
      feature_table = feature_table,
      metadata = metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      covariates = covariates,
      method = "maaslin2",
      effect_type = "coef",
      notes = native_check$note,
      control = control,
      alpha = alpha
    ))
  }
  Go_RunNativeAdapter(
    method_name = "Maaslin2",
    package_names = "Maaslin2",
    native_fun = function() {
      prepared <- Go_PrepareDAInputs(
        feature_table, metadata, group_var, group_1, group_2,
        covariates = covariates, random_effects = random_effects
      )
      method_input <- Go_PrepareMethodInput(prepared, "maaslin2", control)
      output_dir <- file.path(
        tempdir(),
        paste0("condadist_maaslin_", as.integer(stats::runif(1, 1, 1e9)))
      )
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      fit <- Maaslin2::Maaslin2(
        input_data = method_input$sample_feature_data,
        input_metadata = method_input$metadata,
        output = output_dir,
        min_abundance = control$min_abundance,
        min_prevalence = control$min_prevalence,
        normalization = control$normalization,
        transform = control$transform,
        analysis_method = control$analysis_method,
        max_significance = control$max_significance,
        fixed_effects = method_input$fixed_effects,
        random_effects = random_effects,
        standardize = control$standardize,
        plot_heatmap = FALSE,
        plot_scatter = FALSE,
        save_scatter = FALSE,
        reference = paste0(prepared$temp_group_var, ",ref")
      )
      res <- fit$results
      res <- res[res$metadata == prepared$temp_group_var & res$value == "cmp", , drop = FALSE]

      out <- Go_CreateBaseDAResult(
        prepared, "maaslin2", "coef",
        Go_CombineNotes("Native Maaslin2 adapter", Go_ControlNote(control))
      )
      out <- Go_FillDAResult(
        base_result = out,
        feature_ids = res$feature,
        coef = res$coef,
        effect_size = -res$coef,
        p_value = res$pval,
        q_value = res$qval
      )
      out$is_significant <- out$q_value < alpha
      out
    },
    fallback_fun = function(note) {
      Go_BasicEffectAdapter(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        covariates = covariates,
        method = "maaslin2",
        effect_type = "coef",
        notes = note,
        control = control,
        alpha = alpha
      )
    }
  )
}

#' corncob adapter
Go_DA_corncob <- function(feature_table, metadata, group_var, group_1, group_2,
                          random_effects = NULL,
                          covariates = NULL, control = NULL, alpha = 0.05) {
  control <- Go_GetDAMethodControls("corncob", control)
  native_check <- Go_ShouldUseNativeAdapter("corncob", feature_table, metadata, group_var, group_1, group_2)
  if (!isTRUE(native_check$ok)) {
    return(Go_BasicEffectAdapter(
      feature_table = feature_table,
      metadata = metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      covariates = covariates,
      method = "corncob",
      effect_type = "coef",
      notes = native_check$note,
      control = control,
      alpha = alpha
    ))
  }
  Go_RunNativeAdapter(
    method_name = "corncob",
    package_names = c("corncob", "phyloseq"),
    native_fun = function() {
      prepared0 <- Go_PrepareDAInputs(
        feature_table, metadata, group_var, group_1, group_2,
        covariates = covariates, random_effects = random_effects
      )

      run_attempt <- function(prepared, this_control, note_label) {
        if (!is.null(this_control$legacy_filter_cutoff) &&
            is.finite(this_control$legacy_filter_cutoff) &&
            this_control$legacy_filter_cutoff > 0) {
          legacy_filtered <- Go_FilterLegacyRelativeMean(
            feature_table = prepared$feature_table,
            metadata = prepared$metadata,
            cutoff = this_control$legacy_filter_cutoff
          )
          prepared$feature_table <- legacy_filtered$feature_table
          prepared$metadata <- legacy_filtered$metadata
        }

        filtered <- Go_FilterByFeatureThresholds(
          feature_table = prepared$feature_table,
          metadata = prepared$metadata,
          min_prevalence = this_control$min_prevalence %||% 0,
          min_total_count = this_control$min_total_count %||% 0
        )
        prepared$feature_table <- Go_RemoveZeroVarianceFeatures(filtered$feature_table)
        prepared$metadata <- filtered$metadata

        method_input <- Go_PrepareMethodInput(prepared, "corncob", this_control)
        ps <- method_input$phyloseq
        fit <- corncob::differentialTest(
          formula = method_input$full_formula,
          phi.formula = method_input$phi_formula,
          formula_null = method_input$null_formula,
          phi.formula_null = method_input$phi_null_formula,
          data = ps,
          test = "Wald",
          boot = this_control$boot,
          fdr_cutoff = this_control$fdr_cutoff,
          filter_discriminant = this_control$filter_discriminant,
          verbose = FALSE
        )

        coef_tab <- t(vapply(
          fit$all_models,
          Go_GetCorncobCoef,
          numeric(2),
          temp_group_var = prepared$temp_group_var
        ))
        out <- Go_CreateBaseDAResult(
          prepared, "corncob", "coef",
          Go_CombineNotes(note_label, Go_ControlNote(this_control))
        )
        out <- Go_FillDAResult(
          base_result = out,
          feature_ids = rownames(prepared$feature_table),
          coef = coef_tab[, "estimate"],
          effect_size = -coef_tab[, "estimate"],
          p_value = unname(fit$p[rownames(prepared$feature_table)]),
          q_value = unname(fit$p_fdr[rownames(prepared$feature_table)])
        )
        out$is_significant <- out$q_value < alpha
        out
      }

      first_try <- try(
        run_attempt(prepared0, control, "Native corncob adapter"),
        silent = TRUE
      )
      if (!inherits(first_try, "try-error")) {
        return(first_try)
      }

      retry_control <- utils::modifyList(
        control,
        list(
          min_prevalence = control$retry_min_prevalence %||% max(control$min_prevalence %||% 0, 0.1),
          min_total_count = control$retry_min_total_count %||% max(control$min_total_count %||% 0, 20),
          legacy_filter_cutoff = control$retry_legacy_filter_cutoff %||% max(control$legacy_filter_cutoff %||% 0, 0.001),
          phi_formula = control$retry_phi_formula %||% "~ 1",
          phi_null_formula = control$retry_phi_null_formula %||% "~ 1",
          filter_discriminant = TRUE
        )
      )
      second_try <- try(
        run_attempt(prepared0, retry_control, "Native corncob adapter (optimized retry)"),
        silent = TRUE
      )
      if (inherits(second_try, "try-error")) {
        stop(attr(second_try, "condition"))
      }
      second_try
    },
    fallback_fun = function(note) {
      Go_BasicEffectAdapter(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        covariates = covariates,
        method = "corncob",
        effect_type = "coef",
        notes = note,
        control = control,
        alpha = alpha
      )
    }
  )
}

#' DESeq2 adapter
Go_DA_deseq2 <- function(feature_table, metadata, group_var, group_1, group_2,
                         random_effects = NULL,
                         covariates = NULL, control = NULL, alpha = 0.05) {
  control <- Go_GetDAMethodControls("deseq2", control)
  native_check <- Go_ShouldUseNativeAdapter("deseq2", feature_table, metadata, group_var, group_1, group_2)
  if (!isTRUE(native_check$ok)) {
    return(Go_BasicEffectAdapter(
      feature_table = feature_table,
      metadata = metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      covariates = covariates,
      method = "deseq2",
      effect_type = "log2FC",
      notes = native_check$note,
      control = control,
      alpha = alpha
    ))
  }
  Go_RunNativeAdapter(
    method_name = "DESeq2",
    package_names = c("DESeq2", "S4Vectors"),
    native_fun = function() {
      prepared <- Go_PrepareDAInputs(
        feature_table, metadata, group_var, group_1, group_2,
        covariates = covariates, random_effects = random_effects
      )
      method_input <- Go_PrepareMethodInput(prepared, "deseq2", control)
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = method_input$counts,
        colData = method_input$metadata,
        design = method_input$design_formula
      )
      if (control$size_factors_type == "poscounts") {
        dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
      }
      dds <- Go_RunDESeqRobust(dds)
      res <- DESeq2::results(
        dds,
        contrast = c(prepared$temp_group_var, "cmp", "ref"),
        alpha = alpha
      )
      res_df <- as.data.frame(res)
      res_df$feature_id <- rownames(res_df)

      out <- Go_CreateBaseDAResult(
        prepared, "deseq2", "log2FC",
        Go_CombineNotes("Native DESeq2 adapter", Go_ControlNote(control))
      )
      out <- Go_FillDAResult(
        base_result = out,
        feature_ids = res_df$feature_id,
        coef = res_df$log2FoldChange,
        effect_size = res_df$log2FoldChange,
        p_value = res_df$pvalue,
        q_value = res_df$padj
      )
      out$is_significant <- out$q_value < alpha
      out
    },
    fallback_fun = function(note) {
      Go_BasicEffectAdapter(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        covariates = covariates,
        method = "deseq2",
        effect_type = "log2FC",
        notes = note,
        control = control,
        alpha = alpha
      )
    }
  )
}
