# Consensus calculations and compatibility columns.

# ------------------------------------------------------------------------------
# Go_CombinePValuesFisher
#   X^2 = -2 * sum(log(p_i))  ~ chisq(df = 2k) under H0
# ------------------------------------------------------------------------------
#' @export
Go_CombinePValuesFisher <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & !is.na(p_values)]
  p_values <- p_values[p_values >= 0 & p_values <= 1]
  if (length(p_values) == 0) return(NA_real_)
  p_values <- pmax(p_values, 1e-300)
  stat <- -2 * sum(log(p_values))
  stats::pchisq(stat, df = 2 * length(p_values), lower.tail = FALSE)
}

# ------------------------------------------------------------------------------
# Go_CombinePValuesCauchy
#   T = sum(w_i * tan((0.5 - p_i) * pi))  ~ Cauchy(0,1) under H0
#   No independence assumption required.
# ------------------------------------------------------------------------------
#' @export
Go_CombinePValuesCauchy <- function(p_values, weights = NULL) {
  if (!is.null(weights) && length(weights) != length(p_values)) {
    stop("`weights` must have the same length as `p_values`.")
  }
  valid <- is.finite(p_values) & !is.na(p_values) & p_values >= 0 & p_values <= 1
  p_values <- p_values[valid]
  if (length(p_values) == 0) return(NA_real_)
  p_values <- pmin(pmax(p_values, 1e-15), 1 - 1e-15)
  n <- length(p_values)
  if (is.null(weights)) {
    weights <- rep(1 / n, n)
  } else {
    weights <- weights[valid]
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) weights <- rep(1, n)
    weights <- weights / sum(weights)
  }
  T_stat <- sum(weights * tan((0.5 - p_values) * pi))
  stats::pcauchy(T_stat, location = 0, scale = 1, lower.tail = FALSE)
}

# ------------------------------------------------------------------------------
# Go_CombinePValuesAdaptiveCauchy
#   Truncated/adaptive Cauchy:
#   - only p-values below `info_threshold` are treated as informative
#   - p ≈ 1 methods abstain instead of actively harming the combination
#   - keeps Cauchy dependence robustness while avoiding conservative-method drag
# ------------------------------------------------------------------------------
#' @export
Go_CombinePValuesAdaptiveCauchy <- function(p_values, weights = NULL, info_threshold = 0.5) {
  valid <- is.finite(p_values) & !is.na(p_values) & p_values >= 0 & p_values <= 1
  p_values <- p_values[valid]
  if (length(p_values) == 0) return(NA_real_)
  p_values <- pmin(pmax(p_values, 1e-15), 1 - 1e-15)

  if (is.null(weights)) {
    weights <- rep(1, length(p_values))
  } else {
    weights <- weights[valid]
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) {
      weights <- rep(1, length(p_values))
    }
  }

  keep <- p_values < info_threshold
  if (!any(keep)) return(1)

  Go_CombinePValuesCauchy(
    p_values = p_values[keep],
    weights = weights[keep]
  )
}

# ------------------------------------------------------------------------------
# Go_CombinePValuesFamilyPartialConjunction
#
# V5 configurable-panel rule:
#   1. Collapse related corncob Wald/LRT tests with Bonferroni min-p.
#   2. Treat other planned methods as separate model families.
#   3. Require all-but-one family support via partial conjunction.
#
# Missing/invalid planned tests are assigned p = 1 so the planned denominator
# remains fixed across features within a run.
# ------------------------------------------------------------------------------
Go_NormalizeFamilyMethod <- function(methods) {
  key <- tolower(gsub("[^[:alnum:]]+", "", as.character(methods)))
  aliases <- c(
    deseq2 = "deseq2",
    aldex2 = "aldex2",
    ancombc2 = "ancombc2",
    corncobwald = "corncob_wald",
    corncoblrt = "corncob_lrt"
  )
  out <- unname(aliases[key])
  out[is.na(out)] <- key[is.na(out)]
  out
}

Go_FamilyPartialConjunctionDetails <- function(p_values, methods,
                                                planned_methods = methods) {
  if (length(p_values) != length(methods)) {
    stop("`p_values` and `methods` must have the same length.")
  }

  method_key <- Go_NormalizeFamilyMethod(methods)
  planned_key <- unique(Go_NormalizeFamilyMethod(planned_methods))
  supported <- c("deseq2", "aldex2", "ancombc2", "corncob_wald", "corncob_lrt")
  if (length(planned_key) == 0L || any(!planned_key %in% supported)) {
    stop("`planned_methods` must contain supported DA methods.")
  }
  p_values <- suppressWarnings(as.numeric(p_values))
  valid <- is.finite(p_values) & p_values >= 0 & p_values <= 1

  extract_slot <- function(method_name) {
    idx <- which(method_key == method_name)
    if (length(idx) > 1L) {
      stop("Duplicate method rows for family partial conjunction: ", method_name)
    }
    if (length(idx) == 0L || !valid[idx]) {
      return(list(p = 1, estimable = FALSE))
    }
    list(p = p_values[idx], estimable = TRUE)
  }

  planned_families <- unique(ifelse(grepl("^corncob_", planned_key), "corncob", planned_key))
  family_p <- setNames(rep(1, length(planned_families)), planned_families)
  family_estimable <- setNames(rep(FALSE, length(planned_families)), planned_families)
  n_corncob_tests_estimable <- 0L
  for (family in planned_families) {
    members <- if (identical(family, "corncob")) {
      intersect(c("corncob_wald", "corncob_lrt"), planned_key)
    } else {
      family
    }
    slots <- lapply(members, extract_slot)
    member_p <- vapply(slots, `[[`, numeric(1), "p")
    member_estimable <- vapply(slots, `[[`, logical(1), "estimable")
    family_p[[family]] <- min(1, length(members) * min(member_p))
    family_estimable[[family]] <- any(member_estimable)
    if (identical(family, "corncob")) {
      n_corncob_tests_estimable <- sum(member_estimable)
    }
  }
  n_families <- length(family_p)
  h <- max(1L, n_families - 1L)
  combined_p <- min(1, (n_families - h + 1L) * sort(family_p)[[h]])

  diagnostic_family_p <- c(deseq2 = 1, aldex2 = 1, ancombc2 = 1, corncob = 1)
  diagnostic_family_p[names(family_p)] <- family_p

  list(
    combined_p = combined_p,
    family_p = diagnostic_family_p,
    n_families_planned = n_families,
    n_families_estimable = sum(family_estimable),
    partial_conjunction_h = h,
    n_corncob_tests_estimable = n_corncob_tests_estimable
  )
}

#' Combine DA method families using the V5 partial-conjunction rule
#' @param p_values Numeric p-values.
#' @param methods Method identifiers paired with `p_values`.
#' @param planned_methods Full method panel planned for the run. Missing planned
#'   tests are retained conservatively with p = 1.
#' @return A single raw partial-conjunction p-value.
#' @export
Go_CombinePValuesFamilyPartialConjunction <- function(p_values, methods,
                                                       planned_methods = methods) {
  Go_FamilyPartialConjunctionDetails(p_values, methods, planned_methods)$combined_p
}

# ------------------------------------------------------------------------------
# Go_CombinedEffectRank
#   Rank-normalise each method's effect_size to [0,1] within the feature set,
#   then average across methods.
#   0 = most negative, 0.5 = neutral, 1 = most positive.
# ------------------------------------------------------------------------------
#' @export
Go_CombinedEffectRank <- function(da_table, feature_ids) {
  methods <- unique(da_table$method[!is.na(da_table$method)])
  n <- length(feature_ids)

  rank_mat <- matrix(
    NA_real_,
    nrow = n, ncol = length(methods),
    dimnames = list(feature_ids, methods)
  )

  for (method in methods) {
    x <- da_table[da_table$method == method, , drop = FALSE]
    # Guard: take the first row per feature_id so rownames are unique.
    # Duplicate (method, feature_id) rows would otherwise cause x[feature_ids, ]
    # to return multiple rows per name, breaking the length assumption.
    x <- x[!duplicated(x$feature_id), , drop = FALSE]
    rownames(x) <- x$feature_id
    eff   <- x[feature_ids, "effect_size"]
    valid <- is.finite(eff)
    if (sum(valid) < 2) next
    r         <- rep(NA_real_, n)
    r[valid]  <- rank(eff[valid], ties.method = "average") / sum(valid)
    rank_mat[, method] <- r
  }

  rowMeans(rank_mat, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# Go_DAConsensus  (V5)
#
# In V4, the consensus "skeleton" is derived from p_combine:
#   - fisher                     => V1_JSD-style skeleton
#   - adaptive_cauchy            => V2_JSD-style skeleton
#   - family_partial_conjunction => V2_JSD-style skeleton
# ------------------------------------------------------------------------------
#' Summarize DA agreement across methods (V5)
#' @param da_table Standardized long-format DA result table.
#' @param alpha BH-adjusted significance threshold.
#' @param p_combine P-value combination rule.
#' @param planned_methods Full method panel planned for the run.
Go_DAConsensus <- function(da_table,
                           alpha = 0.05,
                           p_combine = c("family_partial_conjunction", "adaptive_cauchy",
                                         "fisher", "cauchy"),
                           planned_methods = unique(da_table$method)) {
  p_combine <- Go_ResolvePCombine(p_combine)
  v4_mode <- Go_ResolveV4Mode(p_combine)
  consensus_skeleton <- v4_mode$consensus_skeleton
  feature_ids <- unique(da_table$feature_id)

  combined_effect_rank <- Go_CombinedEffectRank(da_table, feature_ids)

  out <- lapply(feature_ids, function(fid) {
    x           <- da_table[da_table$feature_id == fid, , drop = FALSE]
    methods_run <- length(unique(x$method[!is.na(x$method)]))
    sig         <- x$is_significant %in% TRUE
    support_score <- if (methods_run == 0) NA_real_ else sum(sig, na.rm = TRUE) / methods_run

    family_details <- if (identical(p_combine, "family_partial_conjunction")) {
      Go_FamilyPartialConjunctionDetails(x$p_value, x$method, planned_methods)
    } else {
      NULL
    }
    combined_p <- switch(
      p_combine,
      fisher = Go_CombinePValuesFisher(x$p_value),
      adaptive_cauchy = Go_CombinePValuesAdaptiveCauchy(x$p_value),
      family_partial_conjunction = family_details$combined_p
    )
    n_p_informative <- sum(is.finite(x$p_value) & !is.na(x$p_value) & x$p_value > 0 & x$p_value < 0.5, na.rm = TRUE)

    direction_consistency <- Go_DirectionConsistency(x$direction[sig])

    effect_consistency <- if (identical(consensus_skeleton, "v2")) {
      # V2_JSD: 0 when no significant methods
      if (sum(sig, na.rm = TRUE) == 0) 0 else Go_EffectConsistency(x$effect_size[sig])
    } else {
      # V1_JSD: delegate fully to Go_EffectConsistency (empty -> 1, single -> 1)
      Go_EffectConsistency(x$effect_size[sig])
    }

    row <- data.frame(
      feature_id            = fid,
      n_methods_run         = methods_run,
      n_methods_significant = sum(sig, na.rm = TRUE),
      DA_support_score      = support_score,
      combined_p            = combined_p,
      cauchy_combined_p     = combined_p,
      n_informative_p       = n_p_informative,
      direction_consistency = direction_consistency,
      effect_consistency    = effect_consistency,
      combined_effect_rank  = unname(combined_effect_rank[fid]),
      min_method_q_value    = suppressWarnings(min(x$q_value, na.rm = TRUE)),
      p_combine_method      = p_combine,
      consensus_skeleton    = consensus_skeleton,
      stringsAsFactors      = FALSE
    )
    if (!is.null(family_details)) {
      row$family_deseq2_p <- unname(family_details$family_p[["deseq2"]])
      row$family_aldex2_p <- unname(family_details$family_p[["aldex2"]])
      row$family_ancombc2_p <- unname(family_details$family_p[["ancombc2"]])
      row$family_corncob_p <- unname(family_details$family_p[["corncob"]])
      row$n_families_planned <- family_details$n_families_planned
      row$n_families_estimable <- family_details$n_families_estimable
      row$partial_conjunction_h <- family_details$partial_conjunction_h
      row$n_corncob_tests_estimable <- family_details$n_corncob_tests_estimable
    }
    if (identical(consensus_skeleton, "v1")) {
      row$mean_effect_size   <- mean(x$effect_size, na.rm = TRUE)
      row$median_effect_size <- stats::median(x$effect_size, na.rm = TRUE)
    }
    row
  })

  out <- do.call(rbind, out)
  out$combined_effect_rank[!is.finite(out$combined_effect_rank)] <- NA_real_
  out$min_method_q_value[!is.finite(out$min_method_q_value)]     <- NA_real_
  if (identical(consensus_skeleton, "v1")) {
    out$mean_effect_size[!is.finite(out$mean_effect_size)]     <- NA_real_
    out$median_effect_size[!is.finite(out$median_effect_size)] <- NA_real_
  }
  out$combined_q              <- stats::p.adjust(out$combined_p, method = "BH")
  out$cauchy_combined_q       <- out$combined_q
  out$is_combined_significant <- out$combined_q < alpha

  if (identical(p_combine, "family_partial_conjunction")) {
    out$family_partial_conjunction_p <- out$combined_p
    out$family_partial_conjunction_q <- out$combined_q
  }

  # V1_JSD compatibility aliases (mirror combined_* columns so downstream
  # consumers that expect fisher_combined_p/q keep working).
  if (identical(p_combine, "fisher")) {
    out$fisher_combined_p <- out$combined_p
    out$fisher_combined_q <- out$combined_q
  }
  out
}

# ------------------------------------------------------------------------------
# Go_FinalScore  (V4)
# ------------------------------------------------------------------------------
#' Combine DA consensus and beta-diversity evidence (V4)
Go_FinalScore <- function(da_consensus, beta_contribution,
                          method_annotation = NULL,
                          beta_enabled = TRUE,
                          weights = c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15)) {
  required_names <- c("da", "beta", "direction", "effect")
  if (!all(required_names %in% names(weights))) {
    stop("weights must be a named numeric vector with elements: ",
         paste(required_names, collapse = ", "))
  }
  weights <- weights[required_names]
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("All weights must be finite non-negative numbers.")
  }
  total <- sum(weights)
  if (!isTRUE(all.equal(total, 1, tolerance = 1e-6))) {
    warning("weights do not sum to 1 (sum = ", round(total, 4), "). ",
            "Scores will not be on a [0, 1] scale.")
  }

  merged <- merge(da_consensus, beta_contribution, by = "feature_id", all = TRUE)
  if (!is.null(method_annotation)) {
    merged <- merge(merged, method_annotation, by = "feature_id", all = TRUE)
  }

  if (!beta_enabled) {
    weights[["beta"]] <- 0
    remaining <- sum(weights[c("da", "direction", "effect")])
    if (remaining > 0) {
      weights[c("da", "direction", "effect")] <-
        weights[c("da", "direction", "effect")] / remaining
    }
  }

  merged$DA_support_score[is.na(merged$DA_support_score)]           <- 0
  merged$direction_consistency[is.na(merged$direction_consistency)] <- 0
  merged$effect_consistency[is.na(merged$effect_consistency)]       <- 0
  if ("loo_separation_score" %in% colnames(merged)) {
    merged$loo_separation_score[is.na(merged$loo_separation_score)] <- 0
  }
  # V1_JSD compatibility: when beta contribution emitted SIMPER+delta_R2,
  # keep the downstream NA-fill behaviour on those columns.
  if ("simper_score" %in% colnames(merged)) {
    merged$simper_score[is.na(merged$simper_score)] <- 0
  }
  if ("delta_R2_score" %in% colnames(merged)) {
    merged$delta_R2_score[is.na(merged$delta_R2_score)] <- 0
  }
  merged$beta_contribution_score[is.na(merged$beta_contribution_score)] <- 0

  combined_q    <- ifelse(is.na(merged$combined_q), 1, merged$combined_q)
  neg_log_q     <- -log10(pmax(combined_q, 1e-300))
  merged$combined_da_score <- Go_NormalizeVector(neg_log_q)
  merged$cauchy_da_score   <- merged$combined_da_score

  # V1_JSD compatibility: when Fisher was used for p-value combination,
  # expose fisher_da_score as an alias for the combined DA score.
  p_combine_used <- unique(merged$p_combine_method)
  p_combine_used <- p_combine_used[!is.na(p_combine_used)]
  if (length(p_combine_used) == 1 && identical(p_combine_used, "fisher")) {
    merged$fisher_da_score <- merged$combined_da_score
  }
  if (length(p_combine_used) == 1 && identical(
    p_combine_used,
    "family_partial_conjunction"
  )) {
    merged$family_partial_conjunction_da_score <- merged$combined_da_score
  }

  merged$core_final_score <-
    weights[["da"]]        * merged$combined_da_score +
    weights[["beta"]]      * merged$beta_contribution_score +
    weights[["direction"]] * merged$direction_consistency +
    weights[["effect"]]    * merged$effect_consistency

  merged$final_score    <- merged$core_final_score
  merged$priority_score <- merged$final_score

  da_sig  <- !is.na(merged$is_combined_significant) & merged$is_combined_significant
  beta_hi <- merged$beta_contribution_score >= 0.5

  merged$classification <- ifelse(
    merged$final_score >= 0.75 & da_sig, "Core_consensus",
    ifelse(
      da_sig & !beta_hi, "Local_DA",
      ifelse(!da_sig & beta_hi, "Structure_driver", "Weak_signal")
    )
  )

  merged$analysis_mode <- if (!beta_enabled) {
    if (max(merged$n_methods_run, na.rm = TRUE) <= 1) "single_method" else "da_only"
  } else {
    "full"
  }

  merged[order(-merged$final_score, merged$feature_id), , drop = FALSE]
}

# ------------------------------------------------------------------------------
# Go_BuildMethodAnnotation  — unchanged from V1
# ------------------------------------------------------------------------------
Go_BuildMethodAnnotation <- function(da_table) {
  methods     <- unique(da_table$method[!is.na(da_table$method)])
  feature_ids <- unique(da_table$feature_id)

  out <- data.frame(
    feature_id       = feature_ids,
    detected_methods = NA_character_,
    stringsAsFactors = FALSE
  )

  detected_list <- vector("list", length(feature_ids))
  names(detected_list) <- feature_ids

  for (method in methods) {
    method_df <- da_table[da_table$method == method, , drop = FALSE]
    rownames(method_df) <- method_df$feature_id
    idx     <- match(feature_ids, method_df$feature_id)
    matched <- method_df[feature_ids[!is.na(idx)], , drop = FALSE]

    out[[paste0(method, "_detected")]]       <- FALSE
    out[[paste0(method, "_is_significant")]] <- FALSE
    out[[paste0(method, "_p_value")]]        <- NA_real_
    out[[paste0(method, "_q_value")]]        <- NA_real_
    out[[paste0(method, "_effect_size")]]    <- NA_real_
    out[[paste0(method, "_direction")]]      <- NA_character_

    if (nrow(matched) == 0) next

    feat <- matched$feature_id
    rows <- match(feat, out$feature_id)
    has_native_result <- is.finite(matched$p_value) | is.finite(matched$effect_size)
    out[[paste0(method, "_detected")]][rows]       <- has_native_result
    out[[paste0(method, "_is_significant")]][rows]  <- matched$is_significant %in% TRUE
    out[[paste0(method, "_p_value")]][rows]         <- matched$p_value
    out[[paste0(method, "_q_value")]][rows]         <- matched$q_value
    out[[paste0(method, "_effect_size")]][rows]     <- matched$effect_size
    out[[paste0(method, "_direction")]][rows]       <- matched$direction

    sig_feat <- feat[matched$is_significant %in% TRUE]
    for (f in sig_feat) detected_list[[f]] <- c(detected_list[[f]], method)
  }

  out$detected_methods <- vapply(
    detected_list,
    function(x) {
      x <- unique(x)
      if (length(x) == 0) return("")
      paste(sort(x), collapse = ";")
    },
    character(1)
  )
  out
}
