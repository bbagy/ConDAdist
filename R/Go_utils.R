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
    # 여러 rank에 같은 값이 반복되는 경우(예: Order=Family=Genus 모두 "Clostridia UCG-014")
    # 중복 제거 후 unique 값의 첫 번째를 사용
    vals <- unique(vals)
    if (length(vals) == 0) {
      return(NA_character_)
    }
    vals[1L]
  })
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
  defaults <- switch(
    method,
    ancombc2 = list(
      prv_cut = 0,
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
      denom = "all",
      use_mc = FALSE,
      paired_test = FALSE,
      zero_replace = TRUE,
      zero_replace_value = 0.5
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
    corncob = list(
      phi_formula = "auto",
      phi_null_formula = "auto",
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
      size_factors_type = "ratio"
    ),
    list()
  )
  utils::modifyList(defaults, control %||% list())
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
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

  if (method == "corncob") {
    out$full_formula <- stats::as.formula(paste("~", out$full_formula_terms))
    out$null_formula <- stats::as.formula(paste("~", Go_BuildFormulaTerms(prepared, include_group = FALSE)))
    out$phi_formula <- if (identical(control$phi_formula, "auto")) {
      out$full_formula
    } else {
      stats::as.formula(control$phi_formula)
    }
    out$phi_null_formula <- if (identical(control$phi_null_formula, "auto")) {
      out$null_formula
    } else {
      stats::as.formula(control$phi_null_formula)
    }
  }

  out
}

Go_PrepareCountMatrix <- function(feature_table, min_count = 0,
                                  zero_replace = FALSE, zero_replace_value = 0.5) {
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
  out$notes <- notes
  out
}

Go_GetMinGroupSize <- function(metadata, group_var, group_1, group_2) {
  keep <- metadata[[group_var]] %in% c(group_1, group_2)
  grp <- metadata[[group_var]][keep]
  if (length(grp) == 0) {
    return(0L)
  }
  as.integer(min(table(grp)))
}

Go_ShouldUseNativeAdapter <- function(method, feature_table, metadata, group_var, group_1, group_2) {
  n_features <- nrow(feature_table)
  min_group_n <- Go_GetMinGroupSize(metadata, group_var, group_1, group_2)

  thresholds <- switch(
    method,
    ancombc2 = list(min_group_n = 3L, min_features = 10L),
    corncob = list(min_group_n = 3L, min_features = 10L),
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
    aitchison = "aitchison"
  )
  if (key %in% names(aliases)) {
    return(unname(aliases[[key]]))
  }
  key
}

Go_AllDAMethods <- function() {
  c("ancombc2", "aldex2", "maaslin2", "corncob", "deseq2")
}

Go_AllDistanceMetrics <- function() {
  c("bray", "jaccard", "aitchison")
}

Go_ResolveMethods <- function(methods) {
  if (is.null(methods) || length(methods) == 0) {
    return(character(0))
  }
  methods <- unique(tolower(as.character(methods)))
  methods[methods == "maaslin"] <- "maaslin2"
  methods
}

Go_ResolveDistances <- function(distances, phy_tree = NULL) {
  if (is.null(distances) || length(distances) == 0) {
    return(NULL)
  }

  out <- unique(vapply(distances, Go_NormalizeDistanceName, character(1)))
  phylo_metrics <- c("unweighted_unifrac", "weighted_unifrac")
  requested_phylo <- intersect(out, phylo_metrics)

  if (is.null(phy_tree)) {
    out <- setdiff(out, phylo_metrics)
    if (length(requested_phylo) > 0) {
      message(
        "No phylogenetic tree was found in psIN, so UniFrac distances were skipped: ",
        paste(requested_phylo, collapse = ", "),
        ". Recommended distances without a tree: bray, jaccard, aitchison."
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
      "Recommended sets: c('bray','jaccard','aitchison') without a tree, ",
      "or c('bray','aitchison','unifrac') with a tree."
    )
    stop("Too many distances were requested: ", paste(out, collapse = ", "))
  }
  out
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
    feature_id = feature_ids,
    simper_contribution = NA_real_,
    simper_score = NA_real_,
    delta_R2 = NA_real_,
    delta_R2_score = NA_real_,
    beta_contribution_score = NA_real_,
    beta_perm_p = NA_real_,
    beta_perm_q = NA_real_,
    stringsAsFactors = FALSE
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
  ordered_methods <- c("deseq2", "aldex2", "ancombc2", "maaslin2", "corncob")
  methods <- ordered_methods[ordered_methods %in% methods]
  if (length(methods) <= 1) {
    return(methods[1] %||% "condadist")
  }
  map <- c(deseq2 = "D", aldex2 = "A", ancombc2 = "N", maaslin2 = "M", corncob = "C")
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

    conda_dist_volcano_dir <- file.path(table_dir, "ConDaDist_volcano")
    createDir(conda_dist_volcano_dir, "ConDaDist_volcano")
    dirs$conda_dist_volcano <- conda_dist_volcano_dir
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
    file.path(getwd(), "R_source", "Gotools", "R", "Go_volcanoPlot_v1.R"),
    file.path(getwd(), "..", "R_source", "Gotools", "R", "Go_volcanoPlot_v1.R"),
    "/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_volcanoPlot_v1.R"
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
    file.path(getwd(), "R_source", "Gotools", "R", "Go_Maaslin2_V2.R"),
    file.path(getwd(), "..", "R_source", "Gotools", "R", "Go_Maaslin2_V2.R"),
    "/Users/heekukpark/Dropbox/04_scripts/R_source/Gotools/R/Go_Maaslin2_V2.R"
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
  available <- all(vapply(package_names, requireNamespace, logical(1), quietly = TRUE))
  if (!available) {
    return(fallback_fun(
      note = paste0("Missing package(s): ", paste(package_names[!vapply(package_names, requireNamespace, logical(1), quietly = TRUE)], collapse = ", "))
    ))
  }

  tryCatch(
    suppressWarnings(suppressMessages(native_fun())),
    error = function(e) {
      fallback_fun(
        note = paste0("Native ", method_name, " adapter failed: ", conditionMessage(e))
      )
    }
  )
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
    corncob = list(phi_formula = "~ 1", phi_null_formula = "~ 1", filter_discriminant = TRUE, boot = FALSE),
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

Go_RunDESeqRobust <- function(dds) {
  tryCatch(
    DESeq2::DESeq(dds, quiet = TRUE),
    error = function(e) {
      if (!grepl("dispersion estimates", conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }
      dds <- DESeq2::estimateSizeFactors(dds)
      dds <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
      DESeq2::dispersions(dds) <- S4Vectors::mcols(dds)$dispGeneEst
      DESeq2::nbinomWaldTest(dds, quiet = TRUE)
    }
  )
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

Go_LeaveOneTaxonOutScore <- function(feature_table, target_feature, group_factor,
                                     distances, full_scores, phy_tree = NULL) {
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
      aitchison          = Go_Dist_aitchison(reduced_table, phy_tree = phy_tree),
      unweighted_unifrac = Go_Dist_unweighted_unifrac(reduced_table, phy_tree = phy_tree),
      weighted_unifrac   = Go_Dist_weighted_unifrac(reduced_table, phy_tree = phy_tree),
      NULL
    )
    if (is.null(reduced_dist)) next

    loo_score <- loo_score + Go_GroupSeparationScore(as.matrix(reduced_dist), group_factor)
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

Go_NormalizeVector <- function(x) {
  n <- length(x)
  if (n == 0) return(x)
  finite_idx <- is.finite(x)
  n_finite <- sum(finite_idx)
  if (n_finite == 0) return(rep(0, n))
  result <- rep(0, n)
  result[finite_idx] <- rank(x[finite_idx], ties.method = "average") / n_finite
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
