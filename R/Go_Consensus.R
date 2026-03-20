#' Summarize DA agreement across methods
Go_DAConsensus <- function(da_table, alpha = 0.05) {
  feature_ids <- unique(da_table$feature_id)

  out <- lapply(feature_ids, function(fid) {
    x <- da_table[da_table$feature_id == fid, , drop = FALSE]
    methods_run <- length(unique(x$method[!is.na(x$method)]))
    sig <- x$is_significant %in% TRUE
    support_score <- if (methods_run == 0) NA_real_ else sum(sig, na.rm = TRUE) / methods_run
    fisher_p <- Go_CombinePValuesFisher(x$p_value)

    direction_consistency <- Go_DirectionConsistency(x$direction[sig])
    effect_consistency <- Go_EffectConsistency(x$effect_size[sig])

    data.frame(
      feature_id            = fid,
      n_methods_run         = methods_run,
      n_methods_significant = sum(sig, na.rm = TRUE),
      DA_support_score      = support_score,
      fisher_combined_p     = fisher_p,
      direction_consistency = direction_consistency,
      effect_consistency    = effect_consistency,
      mean_effect_size      = mean(x$effect_size, na.rm = TRUE),
      median_effect_size    = stats::median(x$effect_size, na.rm = TRUE),
      min_method_q_value    = suppressWarnings(min(x$q_value, na.rm = TRUE)),
      stringsAsFactors      = FALSE
    )
  })

  out <- do.call(rbind, out)
  out$mean_effect_size[!is.finite(out$mean_effect_size)] <- NA_real_
  out$median_effect_size[!is.finite(out$median_effect_size)] <- NA_real_
  out$min_method_q_value[!is.finite(out$min_method_q_value)] <- NA_real_
  out$fisher_combined_q <- stats::p.adjust(out$fisher_combined_p, method = "BH")
  out$is_combined_significant <- out$fisher_combined_q < alpha
  out
}

#' Combine DA consensus and beta-diversity evidence
Go_FinalScore <- function(da_consensus, beta_contribution,
                          method_annotation = NULL,
                          beta_enabled = TRUE,
                          weights = c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15)) {
  required_names <- c("da", "beta", "direction", "effect")
  if (!all(required_names %in% names(weights))) {
    stop(
      "weights must be a named numeric vector with elements: ",
      paste(required_names, collapse = ", ")
    )
  }
  weights <- weights[required_names]
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("All weights must be finite non-negative numbers.")
  }
  total <- sum(weights)
  if (!isTRUE(all.equal(total, 1, tolerance = 1e-6))) {
    warning(
      "weights do not sum to 1 (sum = ", round(total, 4), "). ",
      "Scores will not be on a [0, 1] scale."
    )
  }

  merged <- merge(da_consensus, beta_contribution, by = "feature_id", all = TRUE)
  if (!is.null(method_annotation)) {
    merged <- merge(merged, method_annotation, by = "feature_id", all = TRUE)
  }

  if (!beta_enabled) {
    weights[["beta"]] <- 0
    remaining <- sum(weights[c("da", "direction", "effect")])
    if (remaining > 0) {
      weights[c("da", "direction", "effect")] <- weights[c("da", "direction", "effect")] / remaining
    }
  }

  merged$DA_support_score[is.na(merged$DA_support_score)]           <- 0
  merged$direction_consistency[is.na(merged$direction_consistency)] <- 0
  merged$effect_consistency[is.na(merged$effect_consistency)]       <- 0
  merged$delta_R2_score[is.na(merged$delta_R2_score)]               <- 0
  merged$simper_score[is.na(merged$simper_score)]                   <- 0
  merged$beta_contribution_score[is.na(merged$beta_contribution_score)] <- 0

  # Fisher q-value를 [0,1] 점수로 변환: q=0 → 1, q=1 → 0
  # rank 정규화로 데이터셋 의존성 제거
  fisher_q <- ifelse(
    is.na(merged$fisher_combined_q), 1, merged$fisher_combined_q
  )
  neg_log_q <- -log10(pmax(fisher_q, 1e-300))
  merged$fisher_da_score <- Go_NormalizeVector(neg_log_q)

  merged$core_final_score <-
    weights[["da"]]        * merged$fisher_da_score +
    weights[["beta"]]      * merged$beta_contribution_score +
    weights[["direction"]] * merged$direction_consistency +
    weights[["effect"]]    * merged$effect_consistency

  merged$final_score <- merged$core_final_score
  # Keep a ranking alias for backward compatibility with older exports/plots.
  merged$priority_score <- merged$final_score

  da_sig <- (
    !is.na(merged$is_combined_significant) & merged$is_combined_significant
  )
  beta_hi <- merged$beta_contribution_score >= 0.5

  merged$classification <- ifelse(
    merged$final_score >= 0.75 & da_sig, "Core_consensus",
    ifelse(
      da_sig & !beta_hi, "Local_DA",
      ifelse(!da_sig & beta_hi, "Structure_driver", "Weak_signal")
    )
  )

  if (!beta_enabled) {
    merged$analysis_mode <- if (max(merged$n_methods_run, na.rm = TRUE) <= 1) "single_method" else "da_only"
  } else {
    merged$analysis_mode <- "full"
  }

  merged[order(-merged$final_score, merged$feature_id), , drop = FALSE]
}

Go_BuildMethodAnnotation <- function(da_table) {
  methods <- unique(da_table$method[!is.na(da_table$method)])
  feature_ids <- unique(da_table$feature_id)

  out <- data.frame(
    feature_id = feature_ids,
    detected_methods = NA_character_,
    stringsAsFactors = FALSE
  )

  detected_list <- vector("list", length(feature_ids))
  names(detected_list) <- feature_ids

  for (method in methods) {
    method_df <- da_table[da_table$method == method, , drop = FALSE]
    rownames(method_df) <- method_df$feature_id
    idx <- match(feature_ids, method_df$feature_id)
    matched <- method_df[feature_ids[!is.na(idx)], , drop = FALSE]

    out[[paste0(method, "_detected")]] <- FALSE
    out[[paste0(method, "_is_significant")]] <- FALSE
    out[[paste0(method, "_p_value")]] <- NA_real_
    out[[paste0(method, "_q_value")]] <- NA_real_
    out[[paste0(method, "_effect_size")]] <- NA_real_
    out[[paste0(method, "_direction")]] <- NA_character_

    if (nrow(matched) == 0) {
      next
    }

    feat <- matched$feature_id
    rows <- match(feat, out$feature_id)
    out[[paste0(method, "_detected")]][rows] <- TRUE
    out[[paste0(method, "_is_significant")]][rows] <- matched$is_significant %in% TRUE
    out[[paste0(method, "_p_value")]][rows] <- matched$p_value
    out[[paste0(method, "_q_value")]][rows] <- matched$q_value
    out[[paste0(method, "_effect_size")]][rows] <- matched$effect_size
    out[[paste0(method, "_direction")]][rows] <- matched$direction

    sig_feat <- feat[matched$is_significant %in% TRUE]
    for (f in sig_feat) {
      detected_list[[f]] <- c(detected_list[[f]], method)
    }
  }

  out$detected_methods <- vapply(
    detected_list,
    function(x) {
      x <- unique(x)
      if (length(x) == 0) {
        return("")
      }
      paste(sort(x), collapse = ";")
    },
    character(1)
  )

  out
}
