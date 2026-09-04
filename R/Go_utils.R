Go_NormalizeInputBundle <- function(feature_table, metadata = NULL, phy_tree = NULL) {
  if (inherits(feature_table, "phyloseq")) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("phyloseq package is required to use phyloseq inputs.")
    }

    ps <- feature_table
    otu <- phyloseq::otu_table(ps)
    mat <- methods::as(otu, "matrix")
    if (!phyloseq::taxa_are_rows(otu)) {
      mat <- t(mat)
    }

    if (is.null(metadata)) {
      metadata <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
    }
    if (is.null(phy_tree) && !is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
      phy_tree <- phyloseq::phy_tree(ps)
    }

    taxonomy <- Go_ExtractTaxonomyTable(ps)

    return(list(
      feature_table = Go_AsMatrix(mat),
      metadata = metadata,
      phy_tree = phy_tree,
      taxonomy = taxonomy,
      input_type = "phyloseq"
    ))
  }

  if (is.null(metadata)) {
    stop("metadata must be provided when feature_table is not a phyloseq object.")
  }

  list(
    feature_table = Go_AsMatrix(feature_table),
    metadata = metadata,
    phy_tree = phy_tree,
    taxonomy = NULL,
    input_type = "matrix"
  )
}

Go_ExtractTaxonomyTable <- function(ps) {
  tax <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
  if (is.null(tax)) {
    return(NULL)
  }
  tax_df <- as.data.frame(methods::as(tax, "matrix"), stringsAsFactors = FALSE)
  tax_df$feature_id <- rownames(tax_df)

  rank_cols <- intersect(
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Rank1"),
    colnames(tax_df)
  )

  tax_df$TaxaName <- apply(tax_df[, rank_cols, drop = FALSE], 1, function(x) {
    vals <- x[!is.na(x) & nzchar(trimws(x)) & x != "__"]
    if (length(vals) == 0) {
      return(NA_character_)
    }
    paste(vals, collapse = "; ")
  })

  short_cols <- intersect(
    c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Rank1"),
    colnames(tax_df)
  )
  tax_df$ShortName <- apply(tax_df[, short_cols, drop = FALSE], 1, function(x) {
    vals <- as.character(x)
    vals <- vals[!is.na(vals) & nzchar(trimws(vals)) & vals != "__"]
    vals <- vals[!grepl("^NA(\\s+NA)*$", vals)]
    vals <- unique(vals)
    if (length(vals) == 0) {
      return(NA_character_)
    }
    vals[1L]
  })
  # Append ASV_ID suffix if present: "Species (ASV000001)"
  if ("ASV_ID" %in% colnames(tax_df)) {
    has_asv <- !is.na(tax_df$ASV_ID) & nzchar(tax_df$ASV_ID)
    has_name <- !is.na(tax_df$ShortName)
    tax_df$ShortName[has_asv & has_name]  <- paste0(tax_df$ShortName[has_asv & has_name],
                                                     " (", tax_df$ASV_ID[has_asv & has_name], ")")
    tax_df$ShortName[has_asv & !has_name] <- tax_df$ASV_ID[has_asv & !has_name]
  }
  tax_df
}

Go_AppendTaxonomy <- function(df, taxonomy) {
  if (is.null(taxonomy) || nrow(df) == 0) {
    return(df)
  }
  out <- merge(df, taxonomy, by = "feature_id", all.x = TRUE)
  if (!"ShortName" %in% colnames(out)) {
    out$ShortName <- out$feature_id
  }
  out$ShortName[is.na(out$ShortName) | !nzchar(out$ShortName)] <- out$feature_id[is.na(out$ShortName) | !nzchar(out$ShortName)]
  if ("TaxaName" %in% colnames(out)) {
    out$TaxaName[is.na(out$TaxaName) | !nzchar(out$TaxaName)] <- out$ShortName[is.na(out$TaxaName) | !nzchar(out$TaxaName)]
  }
  out
}

Go_AssertInputs <- function(feature_table, metadata, group_var, group_1, group_2) {
  if (is.null(rownames(metadata))) {
    stop("metadata row names must contain sample identifiers.")
  }
  if (is.null(colnames(feature_table))) {
    stop("feature_table columns must contain sample identifiers.")
  }
  if (!group_var %in% colnames(metadata)) {
    stop("group_var was not found in metadata.")
  }
  groups <- unique(as.character(metadata[[group_var]]))
  if (!group_1 %in% groups || !all(group_2 %in% groups)) {
    stop("group_1 and group_2 must be present in metadata[[group_var]].")
  }
  invisible(TRUE)
}

Go_ResolveOrderedLevels <- function(metadata, group_var, orders = NULL) {
  groups <- unique(as.character(metadata[[group_var]]))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (is.null(orders) || length(orders) == 0) {
    return(groups)
  }
  ordered <- intersect(as.character(orders), groups)
  extra <- setdiff(groups, ordered)
  c(ordered, extra)
}

Go_BuildComparisonPlan <- function(metadata, group_var, group_1, group_2,
                                   orders = NULL, pairwise_all = FALSE) {
  ordered_levels <- Go_ResolveOrderedLevels(metadata = metadata, group_var = group_var, orders = orders)

  if (isTRUE(pairwise_all)) {
    if (length(ordered_levels) < 2) {
      stop("At least two levels are required in `group_var` to run pairwise comparisons.")
    }
    pair_mat <- utils::combn(ordered_levels, 2)
    return(data.frame(
      group_1 = pair_mat[1, ],
      group_2 = pair_mat[2, ],
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(group_1) || !nzchar(group_1)) {
    stop("`group_1` must be provided unless `pairwise_all = TRUE`.")
  }

  targets <- unique(as.character(group_2))
  targets <- targets[!is.na(targets) & nzchar(targets)]

  if (!group_1 %in% ordered_levels || !all(targets %in% ordered_levels)) {
    stop("`group_1` and `group_2` must be present in metadata[[group_var]].")
  }

  targets <- ordered_levels[ordered_levels %in% targets & ordered_levels != group_1]
  if (length(targets) == 0) {
    stop("No valid target groups remained after applying `orders`.")
  }

  data.frame(
    group_1 = rep(group_1, length(targets)),
    group_2 = targets,
    stringsAsFactors = FALSE
  )
}

Go_AlignInputs <- function(feature_table, metadata) {
  feature_table <- Go_AsMatrix(feature_table)
  sample_ids <- intersect(colnames(feature_table), rownames(metadata))
  if (length(sample_ids) == 0) {
    stop("No overlapping sample IDs between feature_table columns and metadata row names.")
  }
  list(
    feature_table = feature_table[, sample_ids, drop = FALSE],
    metadata = metadata[sample_ids, , drop = FALSE]
  )
}

Go_PrepareDAInputs <- function(feature_table, metadata, group_var, group_1, group_2,
                               covariates = NULL, random_effects = NULL) {
  feature_table <- Go_AsMatrix(feature_table)
  keep_cols <- c(group_var, covariates, random_effects)
  keep_cols <- unique(keep_cols[!is.na(keep_cols) & nzchar(keep_cols)])
  keep_cols <- intersect(keep_cols, colnames(metadata))
  if (!group_var %in% keep_cols) {
    stop("group_var must be present in metadata.")
  }

  keep_samples <- metadata[[group_var]] %in% c(group_1, group_2)
  md <- metadata[keep_samples, keep_cols, drop = FALSE]
  ft <- feature_table[, rownames(md), drop = FALSE]

  complete_rows <- stats::complete.cases(md)
  md <- md[complete_rows, , drop = FALSE]
  ft <- ft[, rownames(md), drop = FALSE]

  md$.conda_group <- factor(
    ifelse(md[[group_var]] == group_1, "ref", "cmp"),
    levels = c("ref", "cmp")
  )

  list(
    feature_table = ft,
    metadata = md,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    covariates = setdiff(keep_cols, group_var),
    random_effects = intersect(random_effects, colnames(md)),
    temp_group_var = ".conda_group",
    comparison = paste0(group_1, "_vs_", group_2)
  )
}

Go_GetDAMethodControls <- function(method, control = NULL) {
  method <- if (identical(method, "maaslin")) "maaslin2" else method
  method <- if (identical(method, "corncob")) "corncob_wald" else method
  defaults <- switch(
    method,
    ancombc2 = list(
      prv_cut = 0.1,
      lib_cut = 1000,
      p_adj_method = "BH",
      pseudo = 0,
      pseudo_sens = FALSE,
      struc_zero = TRUE,
      neg_lb = TRUE,
      global = TRUE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE
    ),
    aldex2 = list(
      mc_samples = 128L,
      denom = "iqlr",
      use_mc = FALSE,
      paired_test = FALSE,
      zero_replace = FALSE,
      zero_replace_value = 0.5,
      seed = 1L
    ),
    maaslin2 = list(
      min_abundance = 0,
      min_prevalence = 0,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 1,
      standardize = FALSE
    ),
    corncob_wald = list(
      phi_formula = "auto",
      phi_null_formula = "auto",
      test_type = "Wald",
      boot = FALSE,
      fdr_cutoff = 1,
      filter_discriminant = TRUE,
      min_prevalence = 0.05,
      min_total_count = 10,
      legacy_filter_cutoff = 0,
      retry_min_prevalence = 0.1,
      retry_min_total_count = 20,
      retry_legacy_filter_cutoff = 0.001,
      retry_phi_formula = "~ 1",
      retry_phi_null_formula = "~ 1"
    ),
    corncob_lrt = list(
      phi_formula = "auto",
      phi_null_formula = "auto",
      test_type = "LRT",
      boot = FALSE,
      fdr_cutoff = 1,
      filter_discriminant = TRUE,
      min_prevalence = 0.05,
      min_total_count = 10,
      legacy_filter_cutoff = 0,
      retry_min_prevalence = 0.1,
      retry_min_total_count = 20,
      retry_legacy_filter_cutoff = 0.001,
      retry_phi_formula = "~ 1",
      retry_phi_null_formula = "~ 1"
    ),
    deseq2 = list(
      min_count = 0,
      zero_replace = FALSE,
      zero_replace_value = 1,
      size_factors_type = "poscounts",
      lfc_shrink = TRUE,
      lfc_shrink_type = "ashr"
    ),
    list()
  )
  utils::modifyList(defaults, control %||% list())
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

Go_ResolvePCombine <- function(p_combine) {
  # V4 is a single version, but we keep "cauchy" as a deprecated alias to avoid
  # breaking older scripts that still pass p_combine = "cauchy".
  p <- match.arg(
    arg = p_combine,
    choices = c("family_partial_conjunction", "adaptive_cauchy", "fisher", "cauchy")
  )
  if (identical(p, "cauchy")) {
    warning("p_combine = \"cauchy\" is deprecated; using \"adaptive_cauchy\" instead.")
    p <- "adaptive_cauchy"
  }
  p
}

Go_DefaultPCombineForMethods <- function(methods) {
  resolved <- Go_ResolveMethods(methods)
  if (length(resolved) == 1L) {
    # A one-method run has nothing to combine. Retain the V2 skeleton so
    # single-method and one-method-plus-beta behavior stays backward compatible.
    return("adaptive_cauchy")
  }
  "family_partial_conjunction"
}

Go_ValidateFamilyPartialConjunctionPanel <- function(methods, p_combine) {
  if (!identical(p_combine, "family_partial_conjunction")) {
    return(invisible(TRUE))
  }
  resolved <- Go_ResolveMethods(methods)
  allowed <- Go_AllDAMethods()
  unsupported <- setdiff(resolved, allowed)
  if (length(unsupported) > 0L) {
    stop("Family partial conjunction does not support: ",
         paste(unsupported, collapse = ", "), ".")
  }
  invisible(TRUE)
}

Go_CDDPresetDefinitions <- function() {
  standard_weights <- c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15)
  broad_panel_config <- list(
    methods = c("deseq2", "aldex2", "ancombc2", "corncob_wald", "corncob_lrt"),
    distances = NULL,
    weights = standard_weights,
    p_combine = "family_partial_conjunction",
    calibration = "conservative_general_default"
  )
  list(
    ## "full_cdd" is the reader-facing preset name used in the manuscript.
    ## "broad_panel" is kept as an exact-duplicate alias, not deprecated, so
    ## every existing script and stored result keeps resolving identically.
    full_cdd = broad_panel_config,
    broad_panel = broad_panel_config
  )
}

Go_ResolveCDDPreset <- function(preset, methods = NULL, distances = NULL,
                                weights = NULL, p_combine = NULL,
                                supplied = list()) {
  preset <- match.arg(preset, c("full_cdd", "broad_panel", "custom"))
  definitions <- Go_CDDPresetDefinitions()
  component_names <- c("methods", "distances", "weights", "p_combine")
  explicitly_supplied <- component_names[vapply(
    component_names, function(x) isTRUE(supplied[[x]]), logical(1)
  )]

  if (!identical(preset, "custom")) {
    if (length(explicitly_supplied) > 0L) {
      stop("Preset `", preset, "` owns its configuration; remove explicit: ",
           paste(explicitly_supplied, collapse = ", "),
           ", or use preset = \"custom\".")
    }
    config <- definitions[[preset]]
  } else {
    config <- list(
      methods = methods %||% Go_AllDAMethods(),
      distances = distances,
      weights = weights %||% c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15),
      p_combine = p_combine,
      calibration = "custom_not_independently_calibrated"
    )
    if (is.null(config$p_combine)) {
      config$p_combine <- Go_DefaultPCombineForMethods(config$methods)
    }
  }
  config$preset <- preset
  config$engine_version <- "V5"
  config
}

Go_ResolveV4Mode <- function(p_combine) {
  p <- Go_ResolvePCombine(p_combine)
  if (identical(p, "fisher")) {
    return(list(
      consensus_skeleton = "v1",
      beta_contribution = "simper_loo"
    ))
  }
  list(
    consensus_skeleton = "v2",
    beta_contribution = "loo_only"
  )
}

Go_PrepareMethodInput <- function(prepared, method, control = list()) {
  ft <- prepared$feature_table
  md <- prepared$metadata

  out <- list(
    feature_table = ft,
    metadata = md,
    phyloseq = Go_CreatePhyloseq(ft, md),
    counts = Go_PrepareCountMatrix(ft),
    condition_vector = as.character(md[[prepared$temp_group_var]]),
    full_formula_terms = Go_BuildFormulaTerms(prepared, include_group = TRUE),
    fixed_effects = c(prepared$temp_group_var, prepared$covariates)
  )

  if (method == "aldex2") {
    out$counts <- Go_PrepareCountMatrix(
      ft,
      zero_replace = control$zero_replace,
      zero_replace_value = control$zero_replace_value
    )
  }

  if (method == "deseq2") {
    out$counts <- Go_PrepareCountMatrix(
      ft,
      min_count = control$min_count,
      zero_replace = control$zero_replace,
      zero_replace_value = control$zero_replace_value
    )
    design_terms <- c(prepared$covariates, prepared$temp_group_var)
    out$design_formula <- stats::as.formula(
      paste("~", paste(vapply(design_terms, Go_BacktickName, character(1)), collapse = " + "))
    )
  }

  if (method %in% c("maaslin", "maaslin2")) {
    out$sample_feature_data <- as.data.frame(t(ft))
  }

  if (method %in% c("corncob", "corncob_wald", "corncob_lrt")) {
    out$full_formula <- stats::as.formula(paste("~", out$full_formula_terms))
    out$null_formula <- stats::as.formula(paste("~", Go_BuildFormulaTerms(prepared, include_group = FALSE)))
    out$phi_formula <- if (identical(control$phi_formula, "auto")) {
      out$full_formula
    } else {
      stats::as.formula(control$phi_formula)
    }
    out$phi_null_formula <- if (identical(control$phi_null_formula, "auto")) {
      out$phi_formula
    } else {
      stats::as.formula(control$phi_null_formula)
    }
  }

  out
}

Go_DetectAbundanceType <- function(feature_table) {
  lib_sizes <- colSums(feature_table, na.rm = TRUE)
  mean_lib <- mean(lib_sizes, na.rm = TRUE)
  is_integer_like <- mean(abs(feature_table - round(feature_table)) < 1e-8, na.rm = TRUE) > 0.99
  prop_like <- mean(abs(lib_sizes - 1) < 0.05, na.rm = TRUE) > 0.8
  percent_like <- mean(abs(lib_sizes - 100) < 5, na.rm = TRUE) > 0.8
  if (prop_like || percent_like) "relative"
  else if (is_integer_like && mean_lib > 1000) "absolute"
  else "unknown"
}

Go_PrepareCountMatrix <- function(feature_table, min_count = 0,
                                  zero_replace = FALSE, zero_replace_value = 0.5) {
  abundance_type <- Go_DetectAbundanceType(feature_table)
  if (abundance_type == "relative") {
    lib_sizes <- colSums(feature_table, na.rm = TRUE)
    mean_lib <- mean(lib_sizes, na.rm = TRUE)
    scale_factor <- if (mean_lib <= 1.5) 1e4L else round(median(lib_sizes))
    feature_table <- feature_table * scale_factor
    message(sprintf(
      "[ConDA] Relative abundance detected. Converted to pseudo-counts (scale = %s).",
      format(scale_factor, big.mark = ",")
    ))
  }
  counts <- round(pmax(feature_table, 0))
  counts[counts < min_count] <- 0
  if (zero_replace) {
    counts[counts == 0] <- zero_replace_value
  }
  counts
}

Go_FilterLegacyRelativeMean <- function(feature_table, metadata, cutoff) {
  feature_table <- Go_AsMatrix(feature_table)
  keep_samples <- colSums(feature_table, na.rm = TRUE) > 1
  if (any(keep_samples)) {
    feature_table <- feature_table[, keep_samples, drop = FALSE]
    metadata <- metadata[colnames(feature_table), , drop = FALSE]
  }
  if (ncol(feature_table) == 0 || nrow(feature_table) == 0) {
    return(list(feature_table = feature_table, metadata = metadata))
  }

  rel <- sweep(feature_table, 2, colSums(feature_table, na.rm = TRUE), "/")
  rel[!is.finite(rel)] <- 0
  keep_taxa <- rowMeans(rel, na.rm = TRUE) >= cutoff
  if (all(!keep_taxa)) {
    keep_taxa[] <- TRUE
  }

  list(
    feature_table = feature_table[keep_taxa, , drop = FALSE],
    metadata = metadata
  )
}

Go_FilterByFeatureThresholds <- function(feature_table, metadata,
                                         min_prevalence = 0,
                                         min_total_count = 0) {
  feature_table <- Go_AsMatrix(feature_table)
  if (nrow(feature_table) == 0 || ncol(feature_table) == 0) {
    return(list(feature_table = feature_table, metadata = metadata))
  }

  prevalence <- rowMeans(feature_table > 0, na.rm = TRUE)
  total_count <- rowSums(feature_table, na.rm = TRUE)
  keep <- prevalence >= min_prevalence & total_count >= min_total_count
  if (all(!keep)) {
    keep[] <- TRUE
  }

  list(
    feature_table = feature_table[keep, , drop = FALSE],
    metadata = metadata
  )
}

Go_AsMatrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("Feature_", seq_len(nrow(x)))
  }
  x
}

Go_RemoveZeroVarianceFeatures <- function(feature_table) {
  feature_table <- Go_AsMatrix(feature_table)
  if (nrow(feature_table) == 0) {
    return(feature_table)
  }
  zero_var <- apply(feature_table, 1, function(x) stats::var(x, na.rm = TRUE) == 0)
  zero_var[is.na(zero_var)] <- FALSE
  feature_table[!zero_var, , drop = FALSE]
}

Go_StandardSchema <- function(feature_ids) {
  data.frame(
    feature_id = feature_ids,
    ASV = feature_ids,
    method = NA_character_,
    comparison = NA_character_,
    coef = NA_real_,
    effect_size = NA_real_,
    effect_type = NA_character_,
    p_value = NA_real_,
    q_value = NA_real_,
    direction = NA_character_,
    is_significant = NA,
    mean_group1 = NA_real_,
    mean_group2 = NA_real_,
    prevalence_group1 = NA_real_,
    prevalence_group2 = NA_real_,
    n_samples = NA_integer_,
    bas.count = NA_integer_,
    smvar.count = NA_integer_,
    notes = NA_character_,
    stringsAsFactors = FALSE
  )
}

Go_MakeLegacyASVKey <- function(feature_ids, max_chars = 100L) {
  feature_ids <- as.character(feature_ids)
  short_ids <- substr(feature_ids, 1, max_chars)
  if (anyDuplicated(short_ids)) {
    short_ids <- make.unique(short_ids)
  }
  short_ids
}

Go_BuildFormulaTerms <- function(prepared, include_group = TRUE) {
  terms <- prepared$covariates
  if (include_group) {
    terms <- c(terms, prepared$temp_group_var)
  }
  if (length(terms) == 0) {
    return("1")
  }
  paste(vapply(terms, Go_BacktickName, character(1)), collapse = " + ")
}

Go_BacktickName <- function(x) {
  if (make.names(x) == x) {
    return(x)
  }
  paste0("`", x, "`")
}

Go_CreateBaseDAResult <- function(prepared, method, effect_type, notes = NA_character_) {
  ft <- prepared$feature_table
  md <- prepared$metadata
  g <- md[[prepared$temp_group_var]]
  g1 <- g == "ref"
  g2 <- g == "cmp"

  out <- Go_StandardSchema(rownames(ft))
  out$method <- method
  out$comparison <- prepared$comparison
  out$effect_type <- effect_type
  out$mean_group1 <- rowMeans(ft[, g1, drop = FALSE], na.rm = TRUE)
  out$mean_group2 <- rowMeans(ft[, g2, drop = FALSE], na.rm = TRUE)
  out$prevalence_group1 <- rowMeans(ft[, g1, drop = FALSE] > 0, na.rm = TRUE)
  out$prevalence_group2 <- rowMeans(ft[, g2, drop = FALSE] > 0, na.rm = TRUE)
  out$n_samples <- ncol(ft)
  out$bas.count <- as.integer(sum(g1))
  out$smvar.count <- as.integer(sum(g2))
  out$notes <- notes
  out
}

Go_GetMinGroupSize <- function(metadata, group_var, group_1, group_2) {
  grp <- as.character(metadata[[group_var]])
  n1 <- sum(grp == as.character(group_1), na.rm = TRUE)
  n2 <- sum(grp == as.character(group_2), na.rm = TRUE)
  as.integer(min(n1, n2))
}

Go_ShouldUseNativeAdapter <- function(method, feature_table, metadata, group_var, group_1, group_2) {
  n_features <- nrow(feature_table)
  min_group_n <- Go_GetMinGroupSize(metadata, group_var, group_1, group_2)

  thresholds <- switch(
    method,
    ancombc2 = list(min_group_n = 3L, min_features = 10L),
    corncob = list(min_group_n = 3L, min_features = 10L),
    corncob_wald = list(min_group_n = 3L, min_features = 10L),
    corncob_lrt = list(min_group_n = 3L, min_features = 10L),
    deseq2 = list(min_group_n = 3L, min_features = 5L),
    maaslin2 = list(min_group_n = 2L, min_features = 2L),
    aldex2 = list(min_group_n = 2L, min_features = 2L),
    list(min_group_n = 2L, min_features = 2L)
  )

  if (min_group_n < thresholds$min_group_n) {
    return(list(
      ok = FALSE,
      note = paste0(
        "Native ", method, " skipped: smallest group has ", min_group_n,
        " sample(s); using robust fallback."
      )
    ))
  }
  if (n_features < thresholds$min_features) {
    return(list(
      ok = FALSE,
      note = paste0(
        "Native ", method, " skipped: only ", n_features,
        " feature(s) retained after filtering; using robust fallback."
      )
    ))
  }

  list(ok = TRUE, note = NA_character_)
}

Go_FillDAResult <- function(base_result, feature_ids, coef, effect_size, p_value,
                            q_value = NULL, notes = NULL) {
  idx <- match(feature_ids, base_result$feature_id)
  idx <- idx[!is.na(idx)]
  feature_ids <- feature_ids[!is.na(match(feature_ids, base_result$feature_id))]

  if (length(idx) == 0) {
    return(base_result)
  }

  coef <- Go_RecycleToLength(coef, length(idx))
  effect_size <- Go_RecycleToLength(effect_size, length(idx))
  p_value <- Go_RecycleToLength(p_value, length(idx))
  if (is.null(q_value)) {
    q_value <- stats::p.adjust(p_value, method = "BH")
  }
  q_value <- Go_RecycleToLength(q_value, length(idx))

  base_result$coef[idx] <- coef
  base_result$effect_size[idx] <- effect_size
  base_result$p_value[idx] <- p_value
  base_result$q_value[idx] <- q_value
  base_result$direction[idx] <- Go_DirectionFromEffect(effect_size)
  if (!is.null(notes)) {
    base_result$notes[idx] <- notes
  }

  base_result
}

Go_RecycleToLength <- function(x, n) {
  if (length(x) == n) {
    return(x)
  }
  rep(x, length.out = n)
}

Go_DirectionFromEffect <- function(effect_size) {
  ifelse(
    is.na(effect_size),
    NA_character_,
    ifelse(effect_size > 0, "up_in_group2", ifelse(effect_size < 0, "up_in_group1", "neutral"))
  )
}

Go_CreateAdapterFallback <- function(feature_ids, method, note) {
  out <- Go_StandardSchema(feature_ids)
  out$method <- method
  out$notes <- note
  out
}

Go_EnsureStandardDA <- function(x, method, comparison, alpha) {
  required <- colnames(Go_StandardSchema("x"))
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      x[[col]] <- NA
    }
  }
  extra_cols <- setdiff(colnames(x), required)
  x <- x[, c(required, extra_cols), drop = FALSE]
  x$method[is.na(x$method)] <- method
  x$comparison[is.na(x$comparison)] <- comparison
  x$is_significant <- ifelse(!is.na(x$q_value), x$q_value < alpha, FALSE)
  x
}

Go_BasicEffectAdapter <- function(feature_table, metadata, group_var, group_1, group_2,
                                  method, effect_type, notes = NULL,
                                  covariates = NULL, control = NULL, alpha = 0.05) {
  prepared <- Go_PrepareDAInputs(
    feature_table = feature_table,
    metadata = metadata,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    covariates = covariates
  )
  ft <- prepared$feature_table
  g1 <- prepared$metadata[[prepared$temp_group_var]] == "ref"
  g2 <- prepared$metadata[[prepared$temp_group_var]] == "cmp"

  effect_size <- log2(
    (rowMeans(ft[, g2, drop = FALSE], na.rm = TRUE) + 1e-08) /
      (rowMeans(ft[, g1, drop = FALSE], na.rm = TRUE) + 1e-08)
  )
  p_value <- apply(ft, 1, function(v) {
    tryCatch(
      stats::wilcox.test(v[g1], v[g2], exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  })
  q_value <- stats::p.adjust(p_value, method = "BH")

  out <- Go_CreateBaseDAResult(
    prepared = prepared,
    method = method,
    effect_type = effect_type,
    notes = Go_CombineNotes(notes, Go_ControlNote(control))
  )
  out <- Go_FillDAResult(
    base_result = out,
    feature_ids = rownames(ft),
    coef = effect_size,
    effect_size = effect_size,
    p_value = p_value,
    q_value = q_value
  )
  out$is_significant <- out$q_value < alpha
  out
}

Go_ControlNote <- function(control) {
  if (is.null(control) || length(control) == 0) {
    return(NULL)
  }
  paste0(
    "Controls: ",
    paste(
      paste(names(control), vapply(control, Go_ControlValueString, character(1)), sep = "="),
      collapse = ", "
    )
  )
}

Go_ControlValueString <- function(x) {
  if (length(x) > 1) {
    return(paste(x, collapse = "|"))
  }
  as.character(x)
}

Go_CombineNotes <- function(...) {
  vals <- c(...)
  vals <- vals[!vapply(vals, is.null, logical(1)) & !is.na(vals) & nzchar(vals)]
  paste(vals, collapse = " ; ")
}

Go_MethodAvailabilityNote <- function(pkg_name) {
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    paste0(pkg_name, " available. Adapter currently uses scaffold effect estimator and can be upgraded to native calls.")
  } else {
    paste0(pkg_name, " not installed. Adapter returned scaffold estimates so orchestration remains testable.")
  }
}

Go_CreatePhyloseq <- function(feature_table, metadata) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq is required for this adapter.")
  }
  phyloseq::phyloseq(
    phyloseq::otu_table(feature_table, taxa_are_rows = TRUE),
    phyloseq::sample_data(metadata)
  )
}

Go_NormalizeDistanceName <- function(x) {
  key <- tolower(x)
  aliases <- c(
    wunifrac = "weighted_unifrac",
    weighted_unifrac = "weighted_unifrac",
    weightedunifrac = "weighted_unifrac",
    unifrac = "unweighted_unifrac",
    unweighted_unifrac = "unweighted_unifrac",
    unweightedunifrac = "unweighted_unifrac",
    bray = "bray",
    jaccard = "jaccard",
    jsd = "jsd",
    jensen_shannon = "jsd",
    jensenshannon = "jsd"
  )
  if (key %in% names(aliases)) {
    return(unname(aliases[[key]]))
  }
  key
}

Go_AllDAMethods <- function() {
  c("deseq2", "aldex2", "ancombc2", "corncob_wald", "corncob_lrt")
}

Go_AllDistanceMetrics <- function() {
  c("bray", "jaccard", "jsd")
}

Go_ResolveMethods <- function(methods) {
  if (is.null(methods) || length(methods) == 0) {
    return(character(0))
  }
  methods <- unique(tolower(as.character(methods)))
  methods[methods == "maaslin"] <- "maaslin2"
  methods[methods == "corncob"] <- "corncob_lrt"
  methods
}

Go_ResolveDistances <- function(distances, phy_tree = NULL) {
  if (is.null(distances) || length(distances) == 0) {
    return(NULL)
  }

  out <- unique(vapply(distances, Go_NormalizeDistanceName, character(1)))
  supported <- c(Go_AllDistanceMetrics(), "unweighted_unifrac", "weighted_unifrac")
  unsupported <- setdiff(out, supported)
  if (length(unsupported) > 0) {
    stop(
      "Unsupported distance(s) for *_JSD line: ",
      paste(unsupported, collapse = ", "),
      ". Allowed distances are: ",
      paste(supported, collapse = ", "),
      "."
    )
  }
  phylo_metrics <- c("unweighted_unifrac", "weighted_unifrac")
  requested_phylo <- intersect(out, phylo_metrics)

  if (is.null(phy_tree)) {
    out <- setdiff(out, phylo_metrics)
    if (length(requested_phylo) > 0) {
      message(
        "No phylogenetic tree was found in psIN, so UniFrac distances were skipped: ",
        paste(requested_phylo, collapse = ", "),
        ". Recommended distances without a tree: bray, jaccard, jsd."
      )
    }
  }

  if (length(out) == 0) {
    return(NULL)
  }
  if (length(out) > 3) {
    message(
      "ConDA-dist currently allows at most 3 distances per run. ",
      "Please reduce `distances` to 3 or fewer representative metrics. ",
      "Recommended sets: c('bray','jaccard','jsd') without a tree, ",
      "or c('bray','jsd','unifrac') with a tree."
    )
    stop("Too many distances were requested: ", paste(out, collapse = ", "))
  }
  out
}


Go_CombinePValuesFisher <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & !is.na(p_values)]
  p_values <- p_values[p_values > 0 & p_values <= 1]
  if (length(p_values) == 0) {
    return(NA_real_)
  }
  stat <- -2 * sum(log(p_values))
  stats::pchisq(stat, df = 2 * length(p_values), lower.tail = FALSE)
}

Go_CreateEmptyBetaDistance <- function() {
  list(
    distance_matrices = list(),
    beta_summary = data.frame(
      distance = character(0),
      statistic = numeric(0),
      p_value = numeric(0),
      notes = character(0),
      stringsAsFactors = FALSE
    )
  )
}

Go_CreateEmptyBetaContribution <- function(feature_ids) {
  data.frame(
    feature_id              = feature_ids,
    loo_separation_delta    = NA_real_,
    loo_separation_score    = NA_real_,
    beta_contribution_score = NA_real_,
    beta_perm_p             = NA_real_,
    beta_perm_q             = NA_real_,
    stringsAsFactors        = FALSE
  )
}

Go_SanitizeName <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) {
    x <- "analysis"
  }
  x
}

Go_ComparisonLabel <- function(group_1, group_2) {
  paste0(Go_SanitizeName(group_1), ".vs.", Go_SanitizeName(group_2))
}

Go_MethodSignature <- function(methods) {
  methods <- unique(tolower(as.character(methods)))
  methods[methods == "maaslin"] <- "maaslin2"
  methods[methods == "corncob"] <- "corncob_lrt"
  ordered_methods <- c("deseq2", "aldex2", "ancombc2", "maaslin2", "corncob_wald", "corncob_lrt")
  methods <- ordered_methods[ordered_methods %in% methods]
  if (length(methods) <= 1) {
    return(methods[1] %||% "condadist")
  }
  map <- c(deseq2 = "D", aldex2 = "A", ancombc2 = "N", maaslin2 = "M", corncob_wald = "W", corncob_lrt = "L")
  paste0(unname(map[methods]), collapse = "")
}

Go_FilePrefix <- function(project, group_1, group_2, name = NULL, methods = NULL) {
  parts <- c(
    if (!is.null(methods) && length(methods) > 0) Go_MethodSignature(methods) else NULL,
    Go_ComparisonLabel(group_1, group_2),
    Go_SanitizeName(project)
  )
  if (!is.null(name) && nzchar(name)) {
    parts <- c(parts, Go_SanitizeName(name))
  }
  paste(parts, collapse = ".")
}

Go_UniquePath <- function(path) {
  if (!dir.exists(path) && !file.exists(path)) {
    return(path)
  }
  idx <- 1L
  repeat {
    candidate <- sprintf("%s_%02d", path, idx)
    if (!dir.exists(candidate) && !file.exists(candidate)) {
      return(candidate)
    }
    idx <- idx + 1L
  }
}

Go_path <- function(project, pdf = "yes", table = "yes", path = NULL) {
  if (is.null(project) || !nzchar(project)) {
    stop("Invalid input: 'project' cannot be NULL or empty.")
  }

  createDir <- function(dirPath, dirType) {
    if (!dir.exists(dirPath)) {
      dir.create(dirPath, recursive = TRUE, showWarnings = FALSE)
      cat(sprintf("%s directory created at: %s\n", dirType, dirPath))
    }
  }

  root_dir <- file.path(sprintf("%s_%s", Go_SanitizeName(project), format(Sys.Date(), "%y%m%d")))
  createDir(root_dir, "Main")

  dirs <- list(main = root_dir)

  if (!is.null(table) && tolower(table) == "yes") {
    table_dir <- file.path(root_dir, "table")
    createDir(table_dir, "Table")
    dirs$table <- table_dir

    conda_dist_dir <- file.path(table_dir, "ConDaDist")
    createDir(conda_dist_dir, "ConDaDist")
    dirs$conda_dist <- conda_dist_dir

    # These directories are created on demand by Go_ExportVolcanoBridge
    dirs$conda_dist_volcano <- file.path(table_dir, "ConDaDist_plot_Tab")
    dirs$conda_dist_single_volcano <- file.path(table_dir, "ConDaDist_plot_single_Tab")
  }

  if (!is.null(pdf) && tolower(pdf) == "yes") {
    pdf_dir <- file.path(root_dir, "pdf")
    createDir(pdf_dir, "PDF")
    dirs$pdf <- pdf_dir

    conda_plot_dir <- file.path(pdf_dir, "ConDaPlot")
    createDir(conda_plot_dir, "ConDaPlot")
    dirs$conda_plot <- conda_plot_dir
  }

  dirs
}

Go_FindGoToolsVolcanoPath <- function() {
  candidates <- c(
    file.path(getwd(), "..", "Gotools", "R", "Go_volcanoPlot_v1.R"),
    file.path(getwd(), "..", "..", "Gotools", "R", "Go_volcanoPlot_v1.R"),
    "/Users/heekukpark/Documents/Myscripts/Gotools/R/Go_volcanoPlot_v1.R"
  )
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || !nzchar(hit)) {
    return(NULL)
  }
  hit
}

Go_FindGoToolsMaaslinPath <- function() {
  candidates <- c(
    file.path(getwd(), "..", "Gotools", "R", "Go_Maaslin2_V2.R"),
    file.path(getwd(), "..", "..", "Gotools", "R", "Go_Maaslin2_V2.R"),
    "/Users/heekukpark/Documents/Myscripts/Gotools/R/Go_Maaslin2_V2.R"
  )
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || !nzchar(hit)) {
    return(NULL)
  }
  hit
}

Go_IsNativeMaaslinSingleMode <- function(methods, distances) {
  identical(Go_ResolveMethods(methods), "maaslin2") && is.null(distances)
}

Go_BuildRandomFormula <- function(random_effects) {
  random_effects <- unique(as.character(random_effects))
  random_effects <- random_effects[!is.na(random_effects) & nzchar(random_effects)]
  if (length(random_effects) == 0) {
    return(NULL)
  }
  paste(vapply(random_effects, function(x) sprintf("(1|%s)", Go_BacktickName(x)), character(1)), collapse = " + ")
}

Go_RunNativeMaaslinSingle <- function(psIN, project, group_var, group_1, group_2,
                                      covariates = NULL, random_effects = NULL,
                                      name = NULL) {
  maaslin_path <- Go_FindGoToolsMaaslinPath()
  if (is.null(maaslin_path)) {
    stop("Go_Maaslin2_V2.R was not found.")
  }
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq is required for native MaAsLin2 single mode.")
  }

  env <- new.env(parent = globalenv())
  sys.source(maaslin_path, envir = env)
  if (!exists("Go_Maaslin2", envir = env, inherits = FALSE)) {
    stop("Go_Maaslin2 was not loaded from Gotools.")
  }

  target_levels <- unique(as.character(group_2))
  keep_levels <- unique(c(group_1, target_levels))
  metadata_df <- as.data.frame(phyloseq::sample_data(psIN))
  keep_cols <- unique(c(group_var, covariates, random_effects))
  keep_cols <- keep_cols[!is.na(keep_cols) & nzchar(keep_cols) & keep_cols %in% colnames(metadata_df)]
  keep_samples <- rownames(metadata_df)[metadata_df[[group_var]] %in% keep_levels]
  ps_sub <- phyloseq::prune_samples(keep_samples, psIN)

  if (length(keep_cols) > 0) {
    meta_sub <- as.data.frame(phyloseq::sample_data(ps_sub))
    keep_complete <- stats::complete.cases(meta_sub[, keep_cols, drop = FALSE])
    ps_sub <- phyloseq::prune_samples(rownames(meta_sub)[keep_complete], ps_sub)
  }

  safe_tag <- function(x) {
    x <- gsub("[^A-Za-z0-9+._-]", "_", x)
    gsub("__+", "_", x)
  }
  tag_FE <- function(fx) sprintf("(FE=%s)", safe_tag(paste(fx, collapse = "+")))
  tag_RE <- function(re) {
    if (is.null(re) || length(re) == 0) "(RE=None)" else sprintf("(RE=%s)", safe_tag(paste(re, collapse = "+")))
  }

  fixed_effects <- c(group_var, covariates)
  orders <- keep_levels
  out_root <- sprintf("%s_%s", Go_SanitizeName(project), format(Sys.Date(), "%y%m%d"))
  out_dir <- file.path(
    out_root, "table", "MaAsLin2",
    sprintf(
      "%s.%s.%s",
      if (is.null(name) || !nzchar(name)) "MaAsLin2.Base" else sprintf("MaAsLin2.%s", safe_tag(name)),
      tag_FE(fixed_effects),
      tag_RE(random_effects)
    )
  )

  env$Go_Maaslin2(
    psIN = ps_sub,
    project = project,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    orders = orders,
    out_dir = out_dir,
    name = name,
    combination = NULL,
    global = FALSE
  )

  out_dir
}

Go_ExportNativeMaaslin2VolcanoBridge <- function(native_dir, bridge_dir, psIN,
                                                 group_var, group_1, group_2,
                                                 name = NULL) {
  results_file <- file.path(native_dir, "all_results.csv")
  if (!file.exists(results_file)) {
    stop("Native MaAsLin2 results were not found: all_results.csv")
  }

  df <- utils::read.csv(results_file, check.names = FALSE, stringsAsFactors = FALSE)
  df <- df[df$metadata == group_var, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No rows for group_var were found in native MaAsLin2 results.")
  }

  target_levels <- unique(as.character(group_2))
  df <- df[df$value %in% target_levels, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No requested target levels were found in native MaAsLin2 results.")
  }

  taxonomy <- NULL
  if (inherits(psIN, "phyloseq")) {
    taxonomy <- Go_ExtractTaxonomyTable(psIN)
  }

  dir.create(bridge_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list()

  for (target in unique(df$value)) {
    x <- df[df$value == target, , drop = FALSE]
    x$feature_id <- x$feature
    x$ASV <- x$feature
    x$maaslin2_coef <- x$coef
    x$maaslin2_pvalue <- x$pval
    x$maaslin2_qvalue <- x$qval
    x$maaslin2.P <- ifelse(x$pval < 0.05, ifelse(x$coef >= 0, "up", "down"), "NS")
    x$maaslin2.FDR <- ifelse(x$qval < 0.05, ifelse(x$coef >= 0, "up", "down"), "NS")
    x$basline <- group_1
    x$smvar <- target
    x$mvar <- group_var
    x$name_token <- if (is.null(name)) NA_character_ else as.character(name)
    x$comparison_token <- paste0(group_1, ".vs.", target, if (is.null(name) || !nzchar(name)) "" else paste0(".", name))
    if (!is.null(taxonomy)) {
      x <- Go_AppendTaxonomy(x, taxonomy)
    }
    comparison_stub <- paste0("(", x$comparison_token[1], ")")
    file <- file.path(bridge_dir, paste0("maaslin2.", comparison_stub, ".volcano_bridge.csv"))
    utils::write.csv(x, file, row.names = FALSE)
    files[[target]] <- file
  }

  list(dir = bridge_dir, files = files)
}

Go_RunVolcanoPlotBridge <- function(project, bridge_dir, comparison_name = NULL, name = NULL) {
  volcano_path <- Go_FindGoToolsVolcanoPath()
  if (is.null(volcano_path)) {
    stop("Go_volcanoPlot_v1.R was not found.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggplot2 and ggrepel are required to run Go_volcanoPlot.")
  }
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrepel))
  env <- new.env(parent = globalenv())
  sys.source(volcano_path, envir = env)
  if (!exists("Go_volcanoPlot", envir = env, inherits = FALSE)) {
    stop("Go_volcanoPlot was not loaded from Gotools.")
  }
  env$Go_volcanoPlot(
    project = project,
    result = bridge_dir,
    fc = 1,
    name = name %||% comparison_name,
    overlaps = 20,
    font = 4,
    height = 7,
    width = 7
  )
  list(
    output_dir = file.path(sprintf("%s_%s", Go_SanitizeName(project), format(Sys.Date(), "%y%m%d")), "pdf"),
    bridge_dir = bridge_dir
  )
}

Go_CreateComparisonDir <- function(conda_dist_dir, methods, group_1, group_2) {
  subdir <- file.path(conda_dist_dir, Go_MethodSignature(methods), Go_ComparisonLabel(group_1, group_2))
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  subdir
}

Go_CreateQCPlotDir <- function(conda_dist_dir, group_1, group_2) {
  qc_root <- file.path(conda_dist_dir, "QC_plot")
  dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
  qc_dir <- file.path(qc_root, Go_ComparisonLabel(group_1, group_2))
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  qc_dir
}

Go_RunNativeAdapter <- function(method_name, package_names, native_fun, fallback_fun) {
  missing_pkgs <- package_names[!vapply(package_names, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "[ConDA] ", method_name, " requires package(s) that are not installed: ",
      paste(missing_pkgs, collapse = ", "), ".\n",
      "Install with: BiocManager::install(c(",
      paste0('"', missing_pkgs, '"', collapse = ", "), "))"
    )
  }

  tryCatch(
    suppressWarnings(suppressMessages(native_fun())),
    error = function(e) {
      message("[ConDA] WARNING: Native ", method_name, " failed — result excluded from consensus. Reason: ", conditionMessage(e))
      fallback_fun(
        note = paste0("Native ", method_name, " adapter failed: ", conditionMessage(e))
      )
    }
  )
}

#' List all ConDA-dist dependencies with their install source
#' @export
Go_DependencyList <- function() {
  list(
    # Bioconductor
    bioc = c(
      "phyloseq",    # core input/output
      "DESeq2",      # DA method
      "S4Vectors",   # required by DESeq2 adapter
      "ALDEx2",      # DA method
      "ANCOMBC",     # DA method (includes ancombc2 + pulls in CVXR)
      "Maaslin2",    # DA method
      "BiocParallel" # used by several Bioc methods
    ),
    # CRAN
    cran = c(
      "corncob",    # DA method
      "vegan",      # beta-diversity distances
      "ggplot2",    # plotting
      "ggrepel",    # label repulsion in plots
      "rmarkdown",  # optional QC report export
      "plotly",     # optional interactive plots
      "htmlwidgets",
      "htmltools"
    )
  )
}

#' Check which ConDA-dist dependencies are missing
#'
#' Returns a named list with $bioc and $cran vectors of missing package names.
#' @export
Go_CheckDependencies <- function() {
  deps <- Go_DependencyList()
  list(
    bioc = deps$bioc[!vapply(deps$bioc, requireNamespace, logical(1), quietly = TRUE)],
    cran = deps$cran[!vapply(deps$cran, requireNamespace, logical(1), quietly = TRUE)]
  )
}

#' Install all missing ConDA-dist dependencies
#'
#' Call this once after loading ConDA-dist to ensure all required packages are
#' available. Bioconductor packages are installed via \code{BiocManager};
#' CRAN packages via \code{install.packages}.
#'
#' @param ask If \code{TRUE} (default in interactive sessions), prompt before
#'   installing. Set \code{FALSE} to install without prompting.
#'
#' @examples
#' \dontrun{
#' condadist_dependency()
#' }
#' @export
condadist_dependency <- function(ask = interactive()) {
  # ANCOMBC depends on CVXR which depends on clarabel (a Rust package).
  # clarabel must be compiled from source and requires the Rust toolchain.
  # Check for cargo (Rust package manager) before attempting installation.
  rust_packages <- c("clarabel", "CVXR")  # packages that need Rust to compile
  needs_rust <- any(!vapply(rust_packages, requireNamespace, logical(1), quietly = TRUE))
  if (needs_rust) {
    cargo_found <- nzchar(Sys.which("cargo"))
    if (!cargo_found) {
      # Try common Rust install locations not always on PATH in R sessions
      cargo_paths <- c(
        path.expand("~/.cargo/bin/cargo"),
        "/usr/local/bin/cargo",
        "/opt/homebrew/bin/cargo"
      )
      cargo_found <- any(file.exists(cargo_paths))
    }
    if (!cargo_found) {
      stop(
        "[ConDA] ANCOMBC requires the 'clarabel' package which must be compiled from Rust source.\n",
        "Rust toolchain (cargo) was not found on this system.\n\n",
        "Install Rust by running this command in your Terminal:\n",
        "  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh\n\n",
        "After installation, restart R and run condadist_dependency() again.\n",
        "If Rust is already installed, make sure 'cargo' is on your PATH."
      )
    }
  }

  deps <- Go_DependencyList()
  all_bioc <- deps$bioc
  all_cran <- deps$cran

  # Separate into: (1) not loadable + not in installed.packages → fresh install
  #                (2) in installed.packages but not loadable → broken, force reinstall
  installed_names <- rownames(utils::installed.packages())

  not_loadable_bioc <- all_bioc[!vapply(all_bioc, requireNamespace, logical(1), quietly = TRUE)]
  not_loadable_cran <- all_cran[!vapply(all_cran, requireNamespace, logical(1), quietly = TRUE)]

  fresh_bioc  <- not_loadable_bioc[!not_loadable_bioc %in% installed_names]
  broken_bioc <- not_loadable_bioc[ not_loadable_bioc %in% installed_names]
  fresh_cran  <- not_loadable_cran[!not_loadable_cran %in% installed_names]
  broken_cran <- not_loadable_cran[ not_loadable_cran %in% installed_names]

  all_needing_action <- c(fresh_bioc, broken_bioc, fresh_cran, broken_cran)

  if (length(all_needing_action) == 0) {
    message("[ConDA] All dependencies are installed.")
    return(invisible(TRUE))
  }

  if (length(c(fresh_bioc, fresh_cran)) > 0)
    message("[ConDA] Missing packages: ", paste(c(fresh_bioc, fresh_cran), collapse = ", "))
  if (length(c(broken_bioc, broken_cran)) > 0)
    message("[ConDA] Installed but broken (will force reinstall): ",
            paste(c(broken_bioc, broken_cran), collapse = ", "))

  if (ask) {
    answer <- readline("[ConDA] Install/repair packages now? [y/N] ")
    if (!tolower(trimws(answer)) %in% c("y", "yes")) {
      message("[ConDA] Installation cancelled.")
      return(invisible(FALSE))
    }
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("[ConDA] Installing BiocManager first...")
    utils::install.packages("BiocManager", quiet = TRUE)
  }

  # terra is not a direct ConDA-dist dependency.
  # It may be pulled in as a transitive dependency by BiocManager.
  # If it is already loadable, skip it entirely to avoid triggering
  # a source recompile that requires GDAL/OpenMP system libraries.
  # If it is missing, warn the user rather than attempting to compile.
  if (!requireNamespace("terra", quietly = TRUE)) {
    message(
      "[ConDA] Note: 'terra' is missing but is not a direct ConDA-dist dependency.\n",
      "  If installation later fails due to terra, install system libraries first:\n",
      "    brew install libomp\n",
      "  Then set ~/.R/Makevars:\n",
      "    LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp\n",
      "    CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp\n",
      "  Then retry: install.packages(\"terra\", repos = \"https://rspatial.r-universe.dev\")"
    )
  }

  # CVXR 1.0-11 must be installed before ANCOMBC.
  # Newer CVXR versions do not export 'solve', which causes ANCOMBC lazy-load failure.
  needs_ancombc <- "ANCOMBC" %in% c(fresh_bioc, broken_bioc)
  if (needs_ancombc) {
    cvxr_ok <- tryCatch({
      if (requireNamespace("CVXR", quietly = TRUE)) {
        # check if solve is exported by the installed version
        "solve" %in% getNamespaceExports("CVXR")
      } else {
        FALSE
      }
    }, error = function(e) FALSE)

    if (!cvxr_ok) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        message("[ConDA] Installing remotes (needed for CVXR version pinning)...")
        utils::install.packages("remotes", quiet = TRUE)
      }
      message("[ConDA] Installing CVXR 1.0-11 (required for ANCOMBC compatibility)...")
      remotes::install_version("CVXR", version = "1.0-11", quiet = TRUE,
                               repos = "https://cloud.r-project.org")
    }
  }

  if (length(fresh_cran) > 0) {
    message("[ConDA] Installing CRAN packages: ", paste(fresh_cran, collapse = ", "))
    utils::install.packages(fresh_cran, quiet = TRUE)
  }
  if (length(broken_cran) > 0) {
    message("[ConDA] Force reinstalling broken CRAN packages: ", paste(broken_cran, collapse = ", "))
    utils::install.packages(broken_cran, quiet = TRUE)
  }

  bioc_to_install <- c(fresh_bioc, broken_bioc)
  if (length(bioc_to_install) > 0) {
    force_flag <- length(broken_bioc) > 0
    message("[ConDA] Installing Bioconductor packages",
            if (force_flag) " (force = TRUE for broken ones)" else "", ": ",
            paste(bioc_to_install, collapse = ", "))
    BiocManager::install(bioc_to_install, ask = FALSE, update = FALSE,
                         force = force_flag)
  }

  still_missing <- c(
    all_bioc[!vapply(all_bioc, requireNamespace, logical(1), quietly = TRUE)],
    all_cran[!vapply(all_cran, requireNamespace, logical(1), quietly = TRUE)]
  )
  if (length(still_missing) > 0) {
    warning("[ConDA] The following packages still could not be loaded after installation: ",
            paste(still_missing, collapse = ", "),
            "\nTry running library(\"", still_missing[1], "\") manually to see the error.")
    return(invisible(FALSE))
  }

  message("[ConDA] All dependencies installed successfully.")
  invisible(TRUE)
}

Go_AdjustBetaPermutations <- function(n_beta_permutations, n_features) {
  if (is.null(n_beta_permutations) || n_beta_permutations < 1L) {
    return(list(n_permutations = 0L, note = "Beta permutation disabled."))
  }
  if (n_features >= 3000L) {
    return(list(
      n_permutations = 0L,
      note = paste0("Beta permutation skipped for ", n_features, " features to avoid excessive runtime.")
    ))
  }
  if (n_features >= 1500L) {
    return(list(
      n_permutations = min(as.integer(n_beta_permutations), 10L),
      note = paste0("Beta permutation reduced for ", n_features, " features.")
    ))
  }
  if (n_features >= 750L) {
    return(list(
      n_permutations = min(as.integer(n_beta_permutations), 25L),
      note = paste0("Beta permutation reduced for ", n_features, " features.")
    ))
  }
  list(
    n_permutations = as.integer(n_beta_permutations),
    note = paste0("Beta permutation used ", as.integer(n_beta_permutations), " iterations.")
  )
}

Go_MergeMethodControls <- function(base_controls = NULL, override_controls = NULL) {
  base_controls <- base_controls %||% list()
  override_controls <- override_controls %||% list()
  methods <- union(names(base_controls), names(override_controls))
  out <- stats::setNames(vector("list", length(methods)), methods)
  for (method in methods) {
    out[[method]] <- utils::modifyList(
      base_controls[[method]] %||% list(),
      override_controls[[method]] %||% list()
    )
  }
  out
}

Go_GetRetryPlan <- function(prevalence, abundance, method_controls = NULL) {
  stable_controls <- list(
    ancombc2 = list(struc_zero = FALSE, pseudo = 1, pseudo_sens = FALSE),
    aldex2 = list(zero_replace = TRUE, zero_replace_value = 0.5, mc_samples = 128L),
    maaslin2 = list(normalization = "TSS", transform = "LOG", min_prevalence = 0.05),
    corncob_wald = list(phi_formula = "~ 1", phi_null_formula = "~ 1", test_type = "Wald", filter_discriminant = TRUE, boot = FALSE),
    corncob_lrt = list(phi_formula = "~ 1", phi_null_formula = "~ 1", test_type = "LRT", filter_discriminant = TRUE, boot = FALSE),
    deseq2 = list(size_factors_type = "poscounts", min_count = 1)
  )

  list(
    list(
      attempt_id = 1L,
      label = "initial",
      prevalence = prevalence,
      abundance = abundance,
      method_controls = method_controls %||% list(),
      notes = "Initial analysis settings."
    ),
    list(
      attempt_id = 2L,
      label = "optimized_retry",
      prevalence = max(prevalence, 0.15),
      abundance = max(abundance, 5e-04),
      method_controls = Go_MergeMethodControls(method_controls, stable_controls),
      notes = "Retry with stricter feature filtering and method-stabilizing controls."
    )
  )
}

Go_ExtractMethodLevelStatus <- function(da_table) {
  if (is.null(da_table) || nrow(da_table) == 0) {
    return(data.frame())
  }
  split_notes <- split(da_table$notes, da_table$method)
  data.frame(
    method = names(split_notes),
    note = vapply(split_notes, function(x) paste(unique(stats::na.omit(x)), collapse = " | "), character(1)),
    stringsAsFactors = FALSE
  )
}

Go_ShouldRetryAnalysis <- function(filtered, da_table, methods, attempt_id, max_attempts) {
  if (attempt_id >= max_attempts) {
    return(FALSE)
  }
  if (is.null(filtered$feature_table) || nrow(filtered$feature_table) < 2) {
    return(TRUE)
  }
  if (is.null(da_table) || nrow(da_table) == 0) {
    return(TRUE)
  }

  method_status <- Go_ExtractMethodLevelStatus(da_table)
  if (nrow(method_status) == 0) {
    return(TRUE)
  }

  failed_native <- grepl("Native .* adapter failed", method_status$note)
  if (any(failed_native)) {
    return(TRUE)
  }

  unresolved <- vapply(methods, function(method) {
    this <- da_table[da_table$method == method, , drop = FALSE]
    if (nrow(this) == 0) {
      return(TRUE)
    }
    all(is.na(this$q_value))
  }, logical(1))
  any(unresolved)
}

Go_BuildOptimizationSummary <- function(attempt_plan, final_attempt_id, final_status,
                                        error_message = NA_character_, da_table = NULL) {
  method_status <- Go_ExtractMethodLevelStatus(da_table)
  data.frame(
    attempt_id = vapply(attempt_plan, `[[`, integer(1), "attempt_id"),
    attempt_label = vapply(attempt_plan, `[[`, character(1), "label"),
    attempt_notes = vapply(attempt_plan, `[[`, character(1), "notes"),
    selected_attempt = vapply(attempt_plan, function(x) x$attempt_id == final_attempt_id, logical(1)),
    final_status = final_status,
    error_message = error_message,
    method_status = if (nrow(method_status) == 0) NA_character_ else paste(
      paste(method_status$method, method_status$note, sep = ": "),
      collapse = " || "
    ),
    stringsAsFactors = FALSE
  )
}

Go_GetANCOMBCCoefColumn <- function(result_df, prefix, temp_group_var = ".conda_group") {
  target <- paste0(prefix, temp_group_var, "cmp")
  if (target %in% colnames(result_df)) {
    return(target)
  }

  matches <- grep(paste0("^", prefix), colnames(result_df), value = TRUE)
  matches <- matches[grepl("cmp$", matches)]
  matches <- matches[grepl(temp_group_var, matches, fixed = TRUE)]
  if (length(matches) > 0) {
    return(matches[1])
  }
  NA_character_
}

Go_GetCorncobCoef <- function(model_summary, temp_group_var = ".conda_group") {
  coef_table <- model_summary$coefficients
  if (is.null(coef_table)) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }
  row_id <- rownames(coef_table)
  target <- paste0("mu.", temp_group_var, "cmp")
  if (!target %in% row_id) {
    hits <- row_id[grepl("^mu\\.", row_id) & grepl("cmp$", row_id) & grepl(temp_group_var, row_id, fixed = TRUE)]
    if (length(hits) == 0) {
      return(c(estimate = NA_real_, p_value = NA_real_))
    }
    target <- hits[1]
  }
  c(
    estimate = coef_table[target, "Estimate"],
    p_value = coef_table[target, "Pr(>|t|)"]
  )
}

Go_RunDESeqRobust <- function(dds, sf_type = "ratio") {
  sf_chain <- unique(c(sf_type, "poscounts", "iterate"))
  last_error <- NULL

  .try_deseq <- function(dds_sf) {
    tryCatch(
      DESeq2::DESeq(dds_sf, quiet = TRUE),
      error = function(e) {
        if (!grepl("dispersion estimates", conditionMessage(e), fixed = TRUE)) {
          return(e)
        }
        tryCatch({
          dds2 <- DESeq2::estimateDispersionsGeneEst(dds_sf, quiet = TRUE)
          DESeq2::dispersions(dds2) <- S4Vectors::mcols(dds2)$dispGeneEst
          DESeq2::nbinomWaldTest(dds2, quiet = TRUE)
        }, error = function(e2) e2)
      }
    )
  }

  for (sf in sf_chain) {
    dds_sf <- tryCatch(
      DESeq2::estimateSizeFactors(dds, type = sf),
      error = function(e) e
    )
    if (inherits(dds_sf, "error")) {
      last_error <- dds_sf
      next
    }
    result <- .try_deseq(dds_sf)
    if (!inherits(result, "error")) {
      return(list(dds = result, sf_type_used = sf))
    }
    last_error <- result
  }

  stop(if (!is.null(last_error)) conditionMessage(last_error)
       else "DESeq2 failed with all size factor methods.")
}

Go_GroupSeparationScore <- function(dist_matrix, group_factor) {
  idx1 <- which(group_factor == levels(group_factor)[1])
  idx2 <- which(group_factor == levels(group_factor)[2])

  within_1 <- if (length(idx1) > 1) mean(dist_matrix[idx1, idx1][upper.tri(dist_matrix[idx1, idx1])], na.rm = TRUE) else 0
  within_2 <- if (length(idx2) > 1) mean(dist_matrix[idx2, idx2][upper.tri(dist_matrix[idx2, idx2])], na.rm = TRUE) else 0
  between <- mean(dist_matrix[idx1, idx2, drop = FALSE], na.rm = TRUE)

  score <- between - mean(c(within_1, within_2), na.rm = TRUE)
  if (!is.finite(score)) {
    score <- 0
  }
  score
}

Go_CompositionalFallbackDist <- function(feature_table) {
  feature_table <- Go_AsMatrix(feature_table)
  sample_mat <- t(feature_table)
  sample_mat <- sample_mat + 0.5
  clr <- log(sample_mat) - rowMeans(log(sample_mat))
  stats::dist(clr, method = "euclidean")
}

Go_ComputeSIMPERContribution <- function(feature_table, group_factor) {
  feature_table <- Go_AsMatrix(feature_table)
  out <- stats::setNames(rep(NA_real_, nrow(feature_table)), rownames(feature_table))
  if (length(unique(group_factor)) < 2) {
    return(out)
  }
  if (!requireNamespace("vegan", quietly = TRUE)) {
    approx_score <- rowMeans(feature_table, na.rm = TRUE)
    names(approx_score) <- rownames(feature_table)
    return(approx_score)
  }

  simper_fit <- tryCatch(
    vegan::simper(t(feature_table), group_factor, permutations = 0),
    error = function(e) NULL
  )
  if (is.null(simper_fit) || length(simper_fit) == 0) {
    return(out)
  }

  pair_means <- lapply(simper_fit, function(x) {
    if (!is.null(x$average)) {
      return(x$average)
    }
    if (!is.null(x$overall)) {
      return(x$overall)
    }
    NULL
  })
  pair_means <- pair_means[!vapply(pair_means, is.null, logical(1))]
  if (length(pair_means) == 0) {
    return(out)
  }

  taxa_union <- unique(unlist(lapply(pair_means, names)))
  mat <- matrix(NA_real_, nrow = length(taxa_union), ncol = length(pair_means),
    dimnames = list(taxa_union, NULL)
  )
  for (i in seq_along(pair_means)) {
    mat[names(pair_means[[i]]), i] <- pair_means[[i]]
  }
  out[rownames(mat)] <- rowMeans(mat, na.rm = TRUE)
  out
}

Go_LeaveOneTaxonOutScore <- function(feature_table, target_feature, group_factor,
                                     distances, full_scores, phy_tree = NULL,
                                     covariate_data = NULL) {
  if (length(distances) == 0) {
    return(0)
  }

  loo_score <- 0
  valid_n   <- 0

  for (metric in names(distances)) {
    dist_obj <- distances[[metric]]
    if (is.null(dist_obj)) next

    full_metric_score <- full_scores[[metric]]
    if (is.null(full_metric_score) || !is.finite(full_metric_score)) next

    reduced_table <- feature_table[setdiff(rownames(feature_table), target_feature), , drop = FALSE]
    reduced_dist <- switch(
      metric,
      bray               = Go_Dist_bray(reduced_table, phy_tree = phy_tree),
      jaccard            = Go_Dist_jaccard(reduced_table, phy_tree = phy_tree),
      jsd                = Go_Dist_jsd(reduced_table, phy_tree = phy_tree),
      unweighted_unifrac = Go_Dist_unweighted_unifrac(reduced_table, phy_tree = phy_tree),
      weighted_unifrac   = Go_Dist_weighted_unifrac(reduced_table, phy_tree = phy_tree),
      NULL
    )
    if (is.null(reduced_dist)) next

    loo_score <- loo_score + Go_BetaGroupScore(
      as.matrix(reduced_dist),
      group_factor,
      covariate_data = covariate_data
    )
    valid_n   <- valid_n + 1
  }

  if (valid_n == 0) {
    return(0)
  }

  full_score_sum <- sum(
    unlist(full_scores[intersect(names(distances), names(full_scores))]),
    na.rm = TRUE
  )
  (full_score_sum - loo_score) / valid_n
}

Go_NormalizeVector <- function(x, clip_probs = c(0.05, 0.95)) {
  n <- length(x)
  if (n == 0) return(x)
  finite_idx <- is.finite(x)
  n_finite <- sum(finite_idx)
  if (n_finite == 0) return(rep(0, n))
  result <- rep(0, n)
  vals <- x[finite_idx]
  ## Robust min-max: scale relative to the 5th-95th percentile range, then
  ## clamp to [0,1] -- not plain rank(x)/n (percentile rank), and not
  ## plain min-max either. Two artifacts were found in sequence:
  ##  1. percentile rank forces a uniform 0-1 distribution regardless of
  ##     x's true spread, so any downstream >=0.5 threshold always split
  ##     the data ~50/50 by construction (Structure_driver counts scaled
  ##     with dataset size, not real signal -- docs/20260813_stage3_...,
  ##     section D/L).
  ##  2. plain min-max, tried as the fix, turned out to be wrecked by
  ##     heavy-tailed inputs (e.g. -log10(q): one feature at q=1e-108 next
  ##     to a median significant feature at q=1e-7 crushes that median
  ##     feature's score to ~0.07, not "clearly significant" -- verified
  ##     on cdi_schubert, where this collapsed Core_consensus from 1164 to
  ##     7 features across 18 datasets, section L continuation). Clipping
  ##     to the 5th-95th percentile before scaling keeps a handful of
  ##     extreme values from setting the whole scale, while still letting
  ##     genuine 0-vs-100%-flat inputs collapse toward a single value
  ##     rather than being forced apart.
  bounds <- stats::quantile(vals, probs = clip_probs, na.rm = TRUE, names = FALSE, type = 7)
  lo <- bounds[1]; hi <- bounds[2]
  result[finite_idx] <- if (hi > lo) {
    pmin(pmax((vals - lo) / (hi - lo), 0), 1)
  } else {
    rep(0, n_finite)  # no variation in the bulk -- nothing stands out
  }
  result
}

Go_DirectionConsistency <- function(direction_vec) {
  direction_vec <- direction_vec[!is.na(direction_vec) & direction_vec != "neutral"]
  if (length(direction_vec) == 0) {
    return(0)
  }
  tab <- table(direction_vec)
  max(tab) / sum(tab)
}

Go_EffectConsistency <- function(effect_vec) {
  effect_vec <- effect_vec[is.finite(effect_vec)]
  if (length(effect_vec) <= 1) {
    return(1)
  }
  sign_consistency <- max(mean(effect_vec >= 0), mean(effect_vec <= 0))
  variability <- stats::sd(abs(effect_vec), na.rm = TRUE)
  sign_consistency * (1 / (1 + variability))
}

# ------------------------------------------------------------------------------
# Go_WeightSensitivity  (V2)
#   Assess rank stability of final_score across a grid of weight perturbations.
#   Returns a data.frame with one row per feature:
#     feature_id, mean_rank, sd_rank, rank_stability (1 - sd/max_sd)
# ------------------------------------------------------------------------------
#' Weight sensitivity analysis for ranking stability
#'
#' Perturbs the da/beta/direction/effect weights on a grid and reports how
#' consistently each feature is ranked across all weight combinations.
#'
#' @param da_consensus  Output of \code{Go_DAConsensus()}.
#' @param beta_contribution Output of \code{Go_BetaContribution()}.
#' @param n_grid Number of grid steps per weight dimension (default 5).
#' @param beta_enabled Logical; whether beta component is active (default TRUE).
#' @return A data.frame ordered by mean_rank with columns:
#'   feature_id, mean_rank, sd_rank, rank_stability.
#' @export
Go_WeightSensitivity <- function(da_consensus, beta_contribution,
                                 n_grid = 5L,
                                 beta_enabled = TRUE) {
  stopifnot(is.data.frame(da_consensus), is.data.frame(beta_contribution))

  # Generate a simplex grid over the four weight dimensions
  grid_vals <- seq(0, 1, length.out = n_grid)
  weight_grid <- expand.grid(
    da        = grid_vals,
    beta      = grid_vals,
    direction = grid_vals,
    effect    = grid_vals
  )
  # Keep only rows that sum > 0
  row_sums <- rowSums(weight_grid)
  weight_grid <- weight_grid[row_sums > 0, , drop = FALSE]
  # Normalise each row to sum 1
  weight_grid <- weight_grid / row_sums[row_sums > 0]

  if (!beta_enabled) {
    weight_grid$beta <- 0
    rs <- rowSums(weight_grid[, c("da", "direction", "effect")])
    zero_rows <- rs == 0
    weight_grid <- weight_grid[!zero_rows, , drop = FALSE]
    rs <- rs[!zero_rows]
    weight_grid[, c("da", "direction", "effect")] <-
      weight_grid[, c("da", "direction", "effect")] / rs
  }

  n_combos   <- nrow(weight_grid)
  feature_ids <- da_consensus$feature_id

  rank_mat <- matrix(NA_real_, nrow = length(feature_ids), ncol = n_combos,
                     dimnames = list(feature_ids, NULL))

  for (i in seq_len(n_combos)) {
    w <- as.numeric(weight_grid[i, ])
    names(w) <- c("da", "beta", "direction", "effect")
    # Go_FinalScore() returns rows sorted by -final_score, so we must align
    # by matching feature_id from the OUTPUT, not from da_consensus input.
    fs_out <- tryCatch(
      Go_FinalScore(
        da_consensus      = da_consensus,
        beta_contribution = beta_contribution,
        beta_enabled      = beta_enabled,
        weights           = w
      ),
      error = function(e) NULL
    )
    if (is.null(fs_out)) {
      scores <- rep(NA_real_, length(feature_ids))
    } else {
      scores <- fs_out$final_score[match(feature_ids, fs_out$feature_id)]
    }
    rank_mat[, i] <- rank(-scores, ties.method = "average", na.last = "keep")
  }

  mean_rank <- rowMeans(rank_mat, na.rm = TRUE)
  sd_rank   <- apply(rank_mat, 1, stats::sd, na.rm = TRUE)
  max_sd    <- max(sd_rank, na.rm = TRUE)
  rank_stability <- if (max_sd > 0) 1 - sd_rank / max_sd else rep(1, length(feature_ids))

  out <- data.frame(
    feature_id      = feature_ids,
    mean_rank       = mean_rank,
    sd_rank         = sd_rank,
    rank_stability  = rank_stability,
    stringsAsFactors = FALSE
  )
  out[order(out$mean_rank), , drop = FALSE]
}
