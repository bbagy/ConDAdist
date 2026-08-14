#' Compute modular beta-diversity distances
Go_BetaDistance <- function(feature_table, metadata, group_var, group_1, group_2,
                            distances = c("bray", "jaccard", "jsd"),
                            phy_tree = NULL,
                            n_permutations = 999L) {
  feature_table <- Go_AsMatrix(feature_table)
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
    group_vec <- metadata[[group_var]]
    keep <- group_vec %in% c(group_1, group_2)
    group_factor <- factor(group_vec[keep], levels = c(group_1, group_2))
    dist_mat <- as.matrix(dist_obj)[keep, keep, drop = FALSE]
    perm_result <- Go_RunPERMANOVA(
      dist_mat       = dist_mat,
      group_factor   = group_factor,
      n_permutations = n_permutations
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
#' @param method Contribution formula. \code{"loo_only"} (default, V2_JSD/V4 behaviour)
#'   uses only the leave-one-taxon-out separation delta. \code{"simper_loo"} (V1_JSD
#'   behaviour) averages SIMPER contribution with the LOO delta.
Go_BetaContribution <- function(feature_table, metadata, group_var, group_1, group_2,
                                beta_distances, phy_tree = NULL,
                                n_beta_permutations = 99L,
                                method = c("loo_only", "simper_loo")) {
  method <- match.arg(method)
  feature_table <- Go_AsMatrix(feature_table)
  feature_ids   <- rownames(feature_table)
  group_vec     <- metadata[[group_var]]
  keep          <- group_vec %in% c(group_1, group_2)
  ft            <- feature_table[, keep, drop = FALSE]
  group_factor  <- factor(group_vec[keep], levels = c(group_1, group_2))

  full_scores <- vapply(
    names(beta_distances$distance_matrices),
    function(metric) {
      dm <- beta_distances$distance_matrices[[metric]]
      if (is.null(dm)) return(NA_real_)
      mat <- as.matrix(dm)
      Go_GroupSeparationScore(mat[keep, keep, drop = FALSE], group_factor)
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
      phy_tree       = phy_tree
    )
  }
  loo_separation_score    <- Go_NormalizeVector(loo_separation_delta)

  if (identical(method, "simper_loo")) {
    simper_raw   <- Go_ComputeSIMPERContribution(feature_table = ft, group_factor = group_factor)
    simper_raw   <- simper_raw[feature_ids]
    simper_score <- Go_NormalizeVector(simper_raw)
    beta_contribution_score <- rowMeans(
      cbind(simper_score, loo_separation_score), na.rm = TRUE
    )
  } else {
    simper_raw   <- rep(NA_real_, length(feature_ids))
    simper_score <- rep(NA_real_, length(feature_ids))
    beta_contribution_score <- loo_separation_score
  }

  perm_plan <- Go_AdjustBetaPermutations(
    n_beta_permutations = n_beta_permutations,
    n_features = length(feature_ids)
  )

  # dist_matrices를 keep 샘플로 서브셋 후 전달 (차원 일치 보장)
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
      method = method
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

  ## JSD(p,q) = H((p+q)/2) - 0.5*H(p) - 0.5*H(q), with Shannon entropy
  ## H(x) = -sum(x*log(x)). Mathematically equivalent to the KL-based
  ## definition but lets each row of the i-vs-rest comparison be computed
  ## in one vectorized pass instead of a per-sample-pair R loop -- this
  ## distance gets recomputed thousands of times by the LOO/permutation
  ## machinery in Go_BetaContribution(), so the per-call cost compounds
  ## heavily (unlike bray/jaccard, which reuse compiled vegan::vegdist()).
  ## All entries of P are guaranteed > 0 by the +0.5 pseudocount above, so
  ## no zero-guard is needed here (verified: identical to the old
  ## zero-guarded KL-divergence loop to floating-point precision).
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

Go_RunPERMANOVA <- function(dist_mat, group_factor, n_permutations = 999L) {
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
    score <- Go_GroupSeparationScore(dist_mat, group_factor)
    return(list(
      statistic = score,
      p_value   = NA_real_,
      notes     = paste0(
        "PERMANOVA skipped: smallest group has ", min_group_size,
        " sample(s); centroid separation proxy reported."
      )
    ))
  }

  result <- tryCatch({
    dist_obj_for_adonis <- stats::as.dist(dist_mat)
    df  <- data.frame(group = group_factor)
    fit <- vegan::adonis2(
      dist_obj_for_adonis ~ group,
      data         = df,
      permutations = n_permutations,
      by           = "margin"
    )
    list(
      statistic = fit[["R2"]][1],
      p_value   = fit[["Pr(>F)"]][1],
      notes     = paste0("PERMANOVA (adonis2); permutations = ", n_permutations, ".")
    )
  }, error = function(e) {
    score <- Go_GroupSeparationScore(dist_mat, group_factor)
    list(
      statistic = score,
      p_value   = NA_real_,
      notes     = paste0(
        "PERMANOVA failed (", conditionMessage(e),
        "); centroid separation proxy reported."
      )
    )
  })

  result
}

Go_BetaPermutationPvalue <- function(ft, group_factor, dist_matrices,
                                     observed_beta_score, phy_tree = NULL,
                                     n_permutations = 99L,
                                     method = c("loo_only", "simper_loo")) {
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
    perm_group <- sample(group_factor)

    perm_full_scores <- vapply(
      names(dist_matrices),
      function(metric) {
        dm <- dist_matrices[[metric]]
        if (is.null(dm)) return(NA_real_)
        Go_GroupSeparationScore(as.matrix(dm), perm_group)
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
            phy_tree       = phy_tree
          ),
          error = function(e) NA_real_
        )
      },
      numeric(1)
    )
    perm_loo_score <- Go_NormalizeVector(perm_loo_delta)

    if (identical(method, "simper_loo")) {
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
