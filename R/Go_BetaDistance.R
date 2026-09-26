#' Compute modular beta-diversity distances
#'
#' @param feature_table Feature-by-sample count matrix.
#' @param metadata Sample metadata aligned to the feature-table columns.
#' @param group_var Metadata column defining the comparison groups.
#' @param group_1 Baseline group.
#' @param group_2 Target group.
#' @param distances Distance metrics to compute.
#' @param phy_tree Optional phylogenetic tree for UniFrac distances.
#' @param n_permutations Number of PERMANOVA permutations.
#' @param covariates Optional metadata columns included before the group term;
#'   the reported statistic is the marginal group partial R-squared.
#' @param strata Optional single metadata column defining permutation blocks.
#' @return A list of distance matrices and PERMANOVA summaries.
Go_BetaDistance <- function(feature_table, metadata, group_var, group_1, group_2,
                            distances = c("bray", "jaccard", "jsd"),
                            phy_tree = NULL,
                            n_permutations = 999L,
                            covariates = NULL,
                            strata = NULL) {
  feature_table <- Go_AsMatrix(feature_table)
  beta_design <- Go_PrepareBetaDesign(
    metadata = metadata,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    covariates = covariates,
    strata = strata
  )
  distances <- vapply(distances, Go_NormalizeDistanceName, character(1))
  calculators <- list(
    bray = Go_Dist_bray,
    jaccard = Go_Dist_jaccard,
    jsd = Go_Dist_jsd,
    unweighted_unifrac = Go_Dist_unweighted_unifrac,
    weighted_unifrac = Go_Dist_weighted_unifrac
  )

  matrices <- vector("list", length(distances))
  summaries <- vector("list", length(distances))
  names(matrices) <- distances
  names(summaries) <- distances

  for (metric in distances) {
    calculator <- calculators[[metric]]
    if (is.null(calculator)) {
      matrices[[metric]] <- NULL
      summaries[[metric]] <- data.frame(
        distance = metric,
        statistic = NA_real_,
        p_value = NA_real_,
        notes = "Unsupported distance metric.",
        stringsAsFactors = FALSE
      )
      next
    }

    dist_obj <- calculator(feature_table = feature_table, phy_tree = phy_tree)
    keep <- beta_design$keep
    group_factor <- beta_design$group_factor
    dist_mat <- as.matrix(dist_obj)[keep, keep, drop = FALSE]
    perm_result <- Go_RunPERMANOVA(
      dist_mat       = dist_mat,
      group_factor   = group_factor,
      n_permutations = n_permutations,
      covariate_data = beta_design$covariate_data,
      strata_factor  = beta_design$strata_factor
    )

    matrices[[metric]] <- dist_obj
    summaries[[metric]] <- data.frame(
      distance  = metric,
      statistic = perm_result$statistic,
      p_value   = perm_result$p_value,
      notes     = perm_result$notes,
      stringsAsFactors = FALSE
    )
  }

  list(
    distance_matrices = matrices,
    beta_summary = do.call(rbind, summaries)
  )
}

#' Estimate feature-level beta contributions
#'
#' @param feature_table Feature-by-sample count matrix.
#' @param metadata Sample metadata aligned to the feature-table columns.
#' @param group_var Metadata column defining the comparison groups.
#' @param group_1 Baseline group.
#' @param group_2 Target group.
#' @param beta_distances Output from \code{Go_BetaDistance()}.
#' @param phy_tree Optional phylogenetic tree for UniFrac distances.
#' @param n_beta_permutations Number of feature-level beta permutations.
#' @param method Contribution formula. \code{"loo_only"} (default, V2_JSD/V4 behaviour)
#'   uses only the leave-one-taxon-out separation delta. \code{"simper_loo"} (V1_JSD
#'   behaviour) averages SIMPER contribution with the LOO delta.
#' @param covariates Optional metadata columns used to calculate adjusted
#'   leave-one-feature-out changes in group partial R-squared.
#' @param strata Optional single metadata column defining permutation blocks.
#' @return A feature-level data frame of distance contributions and optional
#'   permutation p-values.
Go_BetaContribution <- function(feature_table, metadata, group_var, group_1, group_2,
                                beta_distances, phy_tree = NULL,
                                n_beta_permutations = 99L,
                                method = c("loo_only", "simper_loo"),
                                covariates = NULL,
                                strata = NULL) {
  method <- match.arg(method)
  feature_table <- Go_AsMatrix(feature_table)
  feature_ids   <- rownames(feature_table)
  beta_design <- Go_PrepareBetaDesign(
    metadata = metadata,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    covariates = covariates,
    strata = strata
  )
  keep          <- beta_design$keep
  ft            <- feature_table[, keep, drop = FALSE]
  group_factor  <- beta_design$group_factor

  full_scores <- vapply(
    names(beta_distances$distance_matrices),
    function(metric) {
      dm <- beta_distances$distance_matrices[[metric]]
      if (is.null(dm)) return(NA_real_)
      mat <- as.matrix(dm)
      Go_BetaGroupScore(
        mat[keep, keep, drop = FALSE],
        group_factor,
        covariate_data = beta_design$covariate_data
      )
    },
    numeric(1)
  )

  loo_separation_delta <- numeric(length(feature_ids))
  for (idx in seq_along(feature_ids)) {
    loo_separation_delta[idx] <- Go_LeaveOneTaxonOutScore(
      feature_table  = ft,
      target_feature = feature_ids[idx],
      group_factor   = group_factor,
      distances      = beta_distances$distance_matrices,
      full_scores    = full_scores,
      phy_tree       = phy_tree,
      covariate_data = beta_design$covariate_data
    )
  }
  loo_separation_score    <- Go_NormalizeVector(loo_separation_delta)

  if (identical(method, "simper_loo")) {
    if (ncol(beta_design$covariate_data) > 0L) {
      warning("SIMPER is not covariate-adjusted; adjusted beta ranking uses the LOO component only.")
    }
    simper_raw   <- Go_ComputeSIMPERContribution(feature_table = ft, group_factor = group_factor)
    simper_raw   <- simper_raw[feature_ids]
    simper_score <- Go_NormalizeVector(simper_raw)
    beta_contribution_score <- if (ncol(beta_design$covariate_data) > 0L) {
      loo_separation_score
    } else {
      rowMeans(cbind(simper_score, loo_separation_score), na.rm = TRUE)
    }
  } else {
    simper_raw   <- rep(NA_real_, length(feature_ids))
    simper_score <- rep(NA_real_, length(feature_ids))
    beta_contribution_score <- loo_separation_score
  }

  perm_plan <- Go_AdjustBetaPermutations(
    n_beta_permutations = n_beta_permutations,
    n_features = length(feature_ids)
  )

  # Restrict distance matrices to the retained samples.
  dist_matrices_sub <- lapply(beta_distances$distance_matrices, function(dm) {
    if (is.null(dm)) return(NULL)
    stats::as.dist(as.matrix(dm)[keep, keep, drop = FALSE])
  })

  if (perm_plan$n_permutations > 0L) {
    beta_perm_p <- Go_BetaPermutationPvalue(
      ft = ft,
      group_factor = group_factor,
      dist_matrices = dist_matrices_sub,
      observed_beta_score = beta_contribution_score,
      phy_tree = phy_tree,
      n_permutations = perm_plan$n_permutations,
      method = method,
      covariate_data = beta_design$covariate_data,
      strata_factor = beta_design$strata_factor
    )
  } else {
    beta_perm_p <- rep(NA_real_, length(feature_ids))
  }

  out <- data.frame(
    feature_id              = feature_ids,
    loo_separation_delta    = unname(loo_separation_delta),
    loo_separation_score    = unname(loo_separation_score),
    beta_contribution_score = unname(beta_contribution_score),
    beta_perm_p             = unname(beta_perm_p),
    beta_perm_q             = stats::p.adjust(unname(beta_perm_p), method = "BH"),
    beta_perm_note          = perm_plan$note,
    stringsAsFactors        = FALSE
  )

  if (identical(method, "simper_loo")) {
    # V1_JSD compatibility columns
    out$simper_contribution <- unname(simper_raw)
    out$simper_score        <- unname(simper_score)
    out$delta_R2            <- unname(loo_separation_delta)
    out$delta_R2_score      <- unname(loo_separation_score)
  }
  out
}

Go_Dist_jaccard <- function(feature_table, phy_tree = NULL) {
  feature_table <- Go_AsMatrix(feature_table)
  binary_table <- t(feature_table > 0) * 1
  if (requireNamespace("vegan", quietly = TRUE)) {
    return(vegan::vegdist(binary_table, method = "jaccard"))
  }
  stats::dist(binary_table, method = "manhattan")
}

Go_Dist_jsd <- function(feature_table, phy_tree = NULL) {
  feature_table <- Go_AsMatrix(feature_table)
  sample_mat <- t(feature_table)
  sample_mat <- sample_mat + 0.5
  sample_mat <- sweep(sample_mat, 1, rowSums(sample_mat), "/")
  sample_mat[!is.finite(sample_mat)] <- 0

  n <- nrow(sample_mat)
  P <- sample_mat
  out <- matrix(0, nrow = n, ncol = n, dimnames = list(rownames(P), rownames(P)))
  if (n < 2) {
    return(stats::as.dist(out))
  }

  # JSD(p,q) = H((p+q)/2) - (H(p)+H(q))/2; H(x) = -sum(x*log(x)).
  # The upstream pseudocount keeps P positive, so no zero guard is needed.
  H <- function(mat) -rowSums(mat * log(mat))
  Hp <- H(P)

  for (i in seq_len(n - 1L)) {
    idx <- (i + 1L):n
    Pi <- matrix(P[i, ], nrow = length(idx), ncol = ncol(P), byrow = TRUE)
    M  <- 0.5 * (Pi + P[idx, , drop = FALSE])
    Hm <- H(M)
    jsd <- Hm - 0.5 * Hp[i] - 0.5 * Hp[idx]
    d <- sqrt(pmax(jsd, 0))
    out[i, idx] <- d
    out[idx, i] <- d
  }

  stats::as.dist(out)
}

Go_Dist_aitchison <- function(feature_table, phy_tree = NULL) {
  Go_CompositionalFallbackDist(feature_table)
}

Go_Dist_bray <- function(feature_table, phy_tree = NULL) {
  feature_table <- Go_AsMatrix(feature_table)
  if (requireNamespace("vegan", quietly = TRUE)) {
    return(vegan::vegdist(t(feature_table), method = "bray"))
  }
  Go_CompositionalFallbackDist(feature_table)
}

Go_Dist_unweighted_unifrac <- function(feature_table, phy_tree = NULL) {
  Go_Dist_unifrac_common(feature_table = feature_table, phy_tree = phy_tree, weighted = FALSE)
}

Go_Dist_weighted_unifrac <- function(feature_table, phy_tree = NULL) {
  Go_Dist_unifrac_common(feature_table = feature_table, phy_tree = phy_tree, weighted = TRUE)
}

Go_Dist_unifrac_common <- function(feature_table, phy_tree, weighted) {
  feature_table <- Go_AsMatrix(feature_table)
  if (is.null(phy_tree)) {
    stop("UniFrac distances require `phy_tree`.")
  }
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq is required for UniFrac distances.")
  }
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(feature_table, taxa_are_rows = TRUE),
    phyloseq::phy_tree(phy_tree)
  )
  phyloseq::UniFrac(ps, weighted = weighted, normalized = TRUE, parallel = FALSE, fast = TRUE)
}

Go_PrepareBetaDesign <- function(metadata, group_var, group_1, group_2,
                                 covariates = NULL, strata = NULL) {
  metadata_rownames <- rownames(metadata)
  metadata <- data.frame(unclass(metadata), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(metadata) <- metadata_rownames
  covariates <- unique(as.character(covariates))
  covariates <- covariates[!is.na(covariates) & nzchar(covariates)]
  strata <- as.character(strata)
  strata <- strata[!is.na(strata) & nzchar(strata)]
  requested <- unique(c(group_var, covariates, strata))
  missing_columns <- setdiff(requested, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop("Unknown beta-model metadata column(s): ", paste(missing_columns, collapse = ", "), ".")
  }
  covariates <- unique(setdiff(covariates, c(group_var, strata)))
  if (length(strata) > 1L) {
    stop("`strata` must be NULL or one metadata column name.")
  }
  if (length(covariates) > 0L && !requireNamespace("vegan", quietly = TRUE)) {
    stop("vegan is required for covariate-adjusted distance analysis.")
  }

  group_vec <- metadata[[group_var]]
  keep <- group_vec %in% c(group_1, group_2)
  complete_columns <- unique(c(group_var, covariates, strata))
  if (length(complete_columns) > 0L) {
    keep <- keep & stats::complete.cases(metadata[, complete_columns, drop = FALSE])
  }
  group_factor <- factor(group_vec[keep], levels = c(group_1, group_2))
  if (any(table(group_factor) == 0L)) {
    stop("Both comparison groups require complete observations for the beta model.")
  }

  covariate_data <- metadata[keep, covariates, drop = FALSE]
  if (ncol(covariate_data) > 0L) {
    colnames(covariate_data) <- paste0(".cdd_cov_", seq_len(ncol(covariate_data)))
  }
  strata_factor <- if (length(strata) == 1L && nzchar(strata)) {
    factor(metadata[[strata]][keep])
  } else {
    NULL
  }
  list(
    keep = keep,
    group_factor = group_factor,
    covariate_data = covariate_data,
    strata_factor = strata_factor,
    covariates = covariates,
    strata = strata
  )
}

Go_BetaGroupScore <- function(dist_mat, group_factor, covariate_data = NULL) {
  if (is.null(covariate_data) || ncol(covariate_data) == 0L) {
    return(Go_GroupSeparationScore(dist_mat, group_factor))
  }
  Go_AdjustedGroupR2(dist_mat, group_factor, covariate_data)
}

Go_AdjustedGroupR2 <- function(dist_mat, group_factor, covariate_data) {
  if (!requireNamespace("vegan", quietly = TRUE)) return(NA_real_)
  tryCatch({
    dist_obj <- stats::as.dist(dist_mat)
    df <- data.frame(covariate_data, group = group_factor, check.names = FALSE)
    covariate_terms <- colnames(covariate_data)
    full_formula <- stats::reformulate(c(covariate_terms, "group"), response = "dist_obj")
    reduced_formula <- stats::reformulate(covariate_terms, response = "dist_obj")
    full_fit <- vegan::dbrda(full_formula, data = df)
    reduced_fit <- vegan::dbrda(reduced_formula, data = df)
    full_constrained <- full_fit$CCA$tot.chi %||% 0
    reduced_constrained <- reduced_fit$CCA$tot.chi %||% 0
    total_inertia <- full_fit$tot.chi
    if (!is.finite(total_inertia) || total_inertia <= 0) return(NA_real_)
    as.numeric((full_constrained - reduced_constrained) / total_inertia)
  }, error = function(e) NA_real_)
}

Go_RunPERMANOVA <- function(dist_mat, group_factor, n_permutations = 999L,
                            covariate_data = NULL, strata_factor = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    score <- Go_GroupSeparationScore(dist_mat, group_factor)
    return(list(
      statistic = score,
      p_value   = NA_real_,
      notes     = "vegan not available; centroid separation proxy reported."
    ))
  }

  min_group_size <- min(table(group_factor))
  if (min_group_size < 3L) {
    adjusted <- !is.null(covariate_data) && ncol(covariate_data) > 0L
    score <- if (adjusted) NA_real_ else Go_GroupSeparationScore(dist_mat, group_factor)
    return(list(
      statistic = score,
      p_value   = NA_real_,
      notes     = paste0(
        "PERMANOVA skipped: smallest group has ", min_group_size,
        " sample(s); ",
        if (adjusted) "no adjusted proxy was substituted." else "centroid separation proxy reported."
      )
    ))
  }

  result <- tryCatch({
    if (!is.null(covariate_data) && ncol(covariate_data) > 0L && n_permutations < 2L) {
      return(list(
        statistic = Go_AdjustedGroupR2(dist_mat, group_factor, covariate_data),
        p_value = NA_real_,
        notes = paste0(
          "Marginal distance-based R2 without permutation; covariates = ",
          ncol(covariate_data), "."
        )
      ))
    }
    dist_obj_for_adonis <- stats::as.dist(dist_mat)
    df <- data.frame(covariate_data, group = group_factor, check.names = FALSE)
    model_terms <- c(colnames(covariate_data), "group")
    model_formula <- stats::reformulate(model_terms, response = "dist_obj_for_adonis")
    fit <- if (is.null(strata_factor)) {
      vegan::adonis2(
        model_formula,
        data = df,
        permutations = n_permutations,
        by = "margin"
      )
    } else {
      vegan::adonis2(
        model_formula,
        data = df,
        permutations = n_permutations,
        by = "margin",
        strata = strata_factor
      )
    }
    group_row <- match("group", rownames(fit))
    if (is.na(group_row)) stop("PERMANOVA did not return the group term.")
    list(
      statistic = fit[["R2"]][group_row],
      p_value   = fit[["Pr(>F)"]][group_row],
      notes     = paste0(
        "Marginal PERMANOVA (adonis2); covariates = ", ncol(covariate_data),
        "; restricted permutation = ", !is.null(strata_factor),
        "; permutations = ", n_permutations, "."
      )
    )
  }, error = function(e) {
    adjusted <- !is.null(covariate_data) && ncol(covariate_data) > 0L
    score <- if (adjusted) NA_real_ else Go_GroupSeparationScore(dist_mat, group_factor)
    list(
      statistic = score,
      p_value   = NA_real_,
      notes     = paste0(
        "PERMANOVA failed (", conditionMessage(e),
        "); ",
        if (adjusted) "no adjusted proxy was substituted." else "centroid separation proxy reported."
      )
    )
  })

  result
}

Go_BetaPermutationPvalue <- function(ft, group_factor, dist_matrices,
                                     observed_beta_score, phy_tree = NULL,
                                     n_permutations = 99L,
                                     method = c("loo_only", "simper_loo"),
                                     covariate_data = NULL,
                                     strata_factor = NULL) {
  method <- match.arg(method)
  n_features <- length(observed_beta_score)

  if (is.null(n_permutations) || n_permutations < 1L) {
    return(rep(NA_real_, n_features))
  }

  # null_mat: permutation × feature
  null_mat <- matrix(NA_real_, nrow = n_permutations, ncol = n_features)

  ## Scoped seeding: same inputs must give the same beta_perm_p/q on every
  ## run (reproducibility), without changing the caller's own RNG stream --
  ## save/restore .Random.seed around the fixed internal seed.
  has_seed <- exists(".Random.seed", envir = .GlobalEnv)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(1L)

  for (perm_i in seq_len(n_permutations)) {
    perm_group <- Go_PermuteGroup(group_factor, strata_factor)

    perm_full_scores <- vapply(
      names(dist_matrices),
      function(metric) {
        dm <- dist_matrices[[metric]]
        if (is.null(dm)) return(NA_real_)
        Go_BetaGroupScore(
          as.matrix(dm), perm_group,
          covariate_data = covariate_data
        )
      },
      numeric(1)
    )

    perm_loo_delta <- vapply(
      rownames(ft),
      function(fid) {
        tryCatch(
          Go_LeaveOneTaxonOutScore(
            feature_table  = ft,
            target_feature = fid,
            group_factor   = perm_group,
            distances      = dist_matrices,
            full_scores    = perm_full_scores,
            phy_tree       = phy_tree,
            covariate_data = covariate_data
          ),
          error = function(e) NA_real_
        )
      },
      numeric(1)
    )
    perm_loo_score <- Go_NormalizeVector(perm_loo_delta)

    if (identical(method, "simper_loo") &&
        (is.null(covariate_data) || ncol(covariate_data) == 0L)) {
      perm_simper <- tryCatch(
        Go_ComputeSIMPERContribution(feature_table = ft, group_factor = perm_group),
        error = function(e) rep(NA_real_, n_features)
      )
      perm_simper_score <- Go_NormalizeVector(unname(perm_simper))
      null_mat[perm_i, ] <- rowMeans(
        cbind(perm_simper_score, perm_loo_score), na.rm = TRUE
      )
    } else {
      null_mat[perm_i, ] <- perm_loo_score
    }
  }

  # empirical p-value: P(null >= observed)
  vapply(
    seq_len(n_features),
    function(j) {
      obs  <- observed_beta_score[j]
      null <- null_mat[, j]
      null <- null[is.finite(null)]
      if (length(null) == 0 || !is.finite(obs)) return(NA_real_)
      (sum(null >= obs) + 1) / (length(null) + 1)
    },
    numeric(1)
  )
}

Go_PermuteGroup <- function(group_factor, strata_factor = NULL) {
  if (is.null(strata_factor)) return(sample(group_factor))
  out <- group_factor
  for (idx in split(seq_along(group_factor), strata_factor, drop = TRUE)) {
    out[idx] <- sample(group_factor[idx])
  }
  factor(out, levels = levels(group_factor))
}
