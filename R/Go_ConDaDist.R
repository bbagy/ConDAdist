#' Run ConDA-dist on a `phyloseq` object
#'
#' `Go_ConDaDist()` is the main orchestration entrypoint for `ConDA-dist`.
#' It supports both:
#'
#' - single-method DA runs with standardized exports
#' - multi-method consensus runs with optional beta-diversity evidence
#'
#' In single mode, the function tries to stay close to the older Go DA-family
#' workflow. In consensus mode, it combines multiple DA methods and optional
#' distance-based feature contribution scores into a final ranked table.
#'
#' @param psIN Standard `phyloseq` object used in the Go tool workflow.
#' @param group_var Column in `sample_data(psIN)` that defines the comparison
#'   groups.
#' @param group_1 Baseline comparison group.
#' @param group_2 One or more target groups compared against `group_1`.
#' @param project Project name used to create the dated output directory.
#' @param covariates Optional character vector of metadata covariates.
#' @param name Optional label appended to output file names.
#' @param random_effects Optional random-effect metadata variables used for
#'   native `MaAsLin2` single-mode runs and `ANCOMBC2` mixed-effects models.
#' @param methods Differential abundance methods to run. Defaults to all
#'   supported methods: `ancombc2`, `aldex2`, `maaslin2`, `corncob`, and
#'   `deseq2`.
#' @param distances Beta-diversity distances to compute. Defaults to the
#'   representative set `bray`, `jaccard`, and `aitchison`. At most 3 distances
#'   are allowed per run. Set to `NULL` to disable beta-diversity and run
#'   DA-only mode.
#' @param method_controls Optional named list of per-method control lists.
#' @param prevalence Minimum fraction of samples with non-zero abundance used in
#'   feature filtering.
#' @param abundance Minimum mean relative abundance threshold used in feature
#'   filtering.
#' @param alpha Significance cutoff used in consensus summaries.
#' @param n_permutations Number of PERMANOVA permutations.
#' @param n_beta_permutations Number of feature-level beta permutations.
#' @param weights Named numeric vector controlling final score weights.
#' @param qc_plot Generate QC plots automatically after analysis.
#' @param volcano_plot Generate volcano plots through Gotools
#'   `Go_volcanoPlot()` using bridge CSV files exported from `ConDA-dist`.
#' @param continue_on_error Keep running remaining pairwise comparisons even if
#'   one comparison fails.
#'
#' @return Output directory path. For a single comparison this is the
#'   comparison directory. For multiple comparisons this is the method-signature
#'   directory containing all comparison subdirectories.
#'
#' @examples
#' \dontrun{
#' # single-method mode
#' res_dir <- Go_ConDaDist(
#'   psIN = ps,
#'   group_var = "TreatmentGroup",
#'   group_1 = "Control",
#'   group_2 = "GLP-2",
#'   project = "DemoProj",
#'   methods = "deseq2",
#'   distances = NULL
#' )
#'
#' # consensus mode
#' res_dir <- Go_ConDaDist(
#'   psIN = ps,
#'   group_var = "TreatmentGroup",
#'   group_1 = "Control",
#'   group_2 = "GLP-2",
#'   project = "DemoProj",
#'   methods = c("deseq2", "aldex2", "ancombc2"),
#'   distances = c("bray", "jaccard", "aitchison")
#' )
#'
#' # native MaAsLin2 single-mode run
#' res_dir <- Go_ConDaDist(
#'   psIN = ps,
#'   group_var = "TreatmentGroup",
#'   group_1 = "Control",
#'   group_2 = c("GLP-2", "D7", "D14"),
#'   project = "DemoProj",
#'   methods = "maaslin2",
#'   distances = NULL,
#'   random_effects = c("SubjectID")
#' )
#' }
Go_ConDaDist <- function(
  psIN,
  group_var,
  group_1,
  group_2,
  project,
  covariates = NULL,
  name = NULL,
  random_effects = NULL,
  methods = Go_AllDAMethods(),
  distances = Go_AllDistanceMetrics(),
  method_controls = NULL,
  prevalence = 0.1,
  abundance = 1e-4,
  alpha = 0.05,
  n_permutations = 999L,
  n_beta_permutations = 99L,
  weights = c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15),
  qc_plot = TRUE,
  volcano_plot = FALSE,
  continue_on_error = TRUE
) {
  result <- Go_RunConDaDistMain(
    feature_table = psIN,
    metadata = NULL,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    project = project,
    covariates = covariates,
    name = name,
    random_effects = random_effects,
    methods = methods,
    distances = distances,
    method_controls = method_controls,
    prevalence = prevalence,
    abundance = abundance,
    alpha = alpha,
    n_permutations = n_permutations,
    n_beta_permutations = n_beta_permutations,
    weights = weights,
    qc_plot = qc_plot,
    volcano_plot = volcano_plot,
    continue_on_error = continue_on_error
  )
  invisible(if (is.list(result) && !is.null(result$return_dir)) result$return_dir else result)
}

Go_RunConDaDistMain <- function(
  feature_table,
  metadata = NULL,
  group_var,
  group_1,
  group_2,
  project,
  covariates = NULL,
  name = NULL,
  random_effects = NULL,
  methods = Go_AllDAMethods(),
  distances = Go_AllDistanceMetrics(),
  method_controls = NULL,
  prevalence = 0.1,
  abundance = 1e-4,
  alpha = 0.05,
  n_permutations = 999L,
  n_beta_permutations = 99L,
  weights = c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15),
  qc_plot = TRUE,
  volcano_plot = FALSE,
  continue_on_error = TRUE
) {
  if (is.null(project) || !nzchar(project)) {
    stop("`project` must be provided.")
  }
  methods <- Go_ResolveMethods(methods)
  normalized_bundle <- Go_NormalizeInputBundle(
    feature_table = feature_table,
    metadata = metadata
  )
  distances <- Go_ResolveDistances(distances, phy_tree = normalized_bundle$phy_tree)

  if (Go_IsNativeMaaslinSingleMode(methods, distances)) {
    message("[ConDA] Single MaAsLin2 mode: using native Go_Maaslin2 workflow with combined levels.")
    native_dir <- Go_RunNativeMaaslinSingle(
      psIN = feature_table,
      project = project,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      covariates = covariates,
      random_effects = random_effects,
      name = name
    )
    if (isTRUE(volcano_plot)) {
      output_layout <- Go_path(project = project, pdf = "no", table = "yes")
      bridge_out <- Go_ExportNativeMaaslin2VolcanoBridge(
        native_dir = native_dir,
        bridge_dir = output_layout$conda_dist_volcano,
        psIN = feature_table,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        name = name
      )
      message("[ConDA] Rendering volcano plots from native MaAsLin2 bridge.")
      Go_RunVolcanoPlotBridge(
        project = project,
        bridge_dir = bridge_out$dir,
        comparison_name = if (length(group_2) == 1) paste0(group_1, ".vs.", group_2) else name,
        name = name
      )
    }
    return(invisible(native_dir))
  }

  output_layout <- Go_path(project = project, pdf = "no", table = "yes")
  root_output_dir <- output_layout$main

  group_targets <- unique(as.character(group_2))
  comparison_results <- lapply(group_targets, function(group_2_target) {
    comparison_dir <- Go_CreateComparisonDir(output_layout$conda_dist, methods, group_1, group_2_target)
    file_prefix <- Go_FilePrefix(project = project, group_1 = group_1, group_2 = group_2_target, name = name, methods = methods)
    message("[ConDA] Starting comparison: ", group_1, " vs ", group_2_target)

    run_one <- function() {
      Go_RunSingleDAensemble(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2_target,
        project = project,
        name = name,
        random_effects = random_effects,
        covariates = covariates,
        methods = methods,
        distances = distances,
        method_controls = method_controls,
        prevalence = prevalence,
        abundance = abundance,
        alpha = alpha,
        n_permutations = n_permutations,
        n_beta_permutations = n_beta_permutations,
        weights = weights,
        qc_plot = qc_plot,
        volcano_plot = volcano_plot,
        volcano_bridge_root_dir = output_layout$conda_dist_volcano,
        output_dir = comparison_dir,
        file_prefix = file_prefix
      )
    }

    if (!isTRUE(continue_on_error)) {
      return(run_one())
    }

    tryCatch(
      run_one(),
      error = function(e) {
        Go_BuildFailedComparison(
          group_1 = group_1,
          group_2 = group_2_target,
          comparison_dir = comparison_dir,
          file_prefix = file_prefix,
          error_message = conditionMessage(e)
        )
      }
    )
  })
  names(comparison_results) <- paste0(group_1, ".vs.", group_targets)

  method_dir <- dirname(comparison_results[[1]]$comparison_dir)

  if (length(comparison_results) == 1) {
    single <- comparison_results[[1]]
    single$project <- project
    single$output_root_dir <- root_output_dir
    single$output_table_dir <- output_layout$table
    single$return_dir <- single$comparison_dir
    return(invisible(single))
  }

  invisible(list(
    project = project,
    output_root_dir = root_output_dir,
    output_table_dir = output_layout$table,
    comparisons = comparison_results,
    return_dir = method_dir
  ))
}

Go_BuildFailedComparison <- function(group_1, group_2, comparison_dir, file_prefix, error_message) {
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

  error_df <- data.frame(
    comparison = paste0(group_1, ".vs.", group_2),
    status = "failed",
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    error_df,
    file.path(comparison_dir, paste0(file_prefix, "run_status.csv")),
    row.names = FALSE
  )

  list(
    comparison = paste0(group_1, ".vs.", group_2),
    comparison_dir = comparison_dir,
    file_prefix = file_prefix,
    status = "failed",
    error_message = error_message,
    filtered = NULL,
    input_bundle = NULL,
    da_raw = NULL,
    da_standardized = NULL,
    beta_distances = Go_CreateEmptyBetaDistance(),
    beta_contribution = data.frame(),
    da_consensus = data.frame(),
    method_annotation = data.frame(),
    final_scores = data.frame(),
    optimization = error_df,
    exported = list(run_status = file.path(comparison_dir, paste0(file_prefix, "run_status.csv")))
  )
}

Go_RunSingleDAensemble <- function(
  feature_table,
  metadata = NULL,
  group_var,
  group_1,
  group_2,
  project,
  name = NULL,
  random_effects = NULL,
  covariates = NULL,
  methods,
  distances,
  method_controls = NULL,
  prevalence = 0.1,
  abundance = 1e-4,
  alpha = 0.05,
  n_permutations = 999L,
  n_beta_permutations = 99L,
  weights = c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15),
  qc_plot = TRUE,
  volcano_plot = FALSE,
  volcano_bridge_root_dir = NULL,
  output_dir,
  file_prefix = NULL
) {
  input_bundle <- Go_NormalizeInputBundle(
    feature_table = feature_table,
    metadata = metadata
  )

  methods <- Go_ResolveMethods(methods)
  distances <- Go_ResolveDistances(distances, phy_tree = input_bundle$phy_tree)

  Go_AssertInputs(
    feature_table = input_bundle$feature_table,
    metadata = input_bundle$metadata,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  attempt_plan <- Go_GetRetryPlan(
    prevalence = prevalence,
    abundance = abundance,
    method_controls = method_controls
  )
  result <- NULL
  final_error <- NULL

  for (attempt in attempt_plan) {
    message(
      "[ConDA] ", group_1, " vs ", group_2,
      " | attempt ", attempt$attempt_id, "/", length(attempt_plan),
      " (", attempt$label, ")"
    )
    attempt_result <- tryCatch(
      Go_RunSingleDAAttempt(
        input_bundle = input_bundle,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        random_effects = random_effects,
        covariates = covariates,
        methods = methods,
        distances = distances,
        method_controls = attempt$method_controls,
        prevalence = attempt$prevalence,
        abundance = attempt$abundance,
        alpha = alpha,
        n_permutations = n_permutations,
        n_beta_permutations = n_beta_permutations,
        weights = weights,
        output_dir = output_dir,
        file_prefix = file_prefix
      ),
      error = function(e) {
        list(status = "failed", error_message = conditionMessage(e))
      }
    )

    attempt_result$attempt_id <- attempt$attempt_id
    attempt_result$attempt_label <- attempt$label
    attempt_result$attempt_notes <- attempt$notes

    if (!identical(attempt_result$status, "completed")) {
      final_error <- attempt_result$error_message
      message("[ConDA] Attempt failed: ", final_error)
      next
    }

    retry_needed <- Go_ShouldRetryAnalysis(
      filtered = attempt_result$filtered,
      da_table = attempt_result$da_standardized$all_methods_standardized,
      methods = methods,
      attempt_id = attempt$attempt_id,
      max_attempts = length(attempt_plan)
    )

    result <- attempt_result
    if (!retry_needed) {
      message("[ConDA] Comparison completed with attempt ", attempt$attempt_id, ".")
      break
    }
    message("[ConDA] Retrying with optimized settings.")
  }

  if (is.null(result)) {
    failed <- Go_BuildFailedComparison(
      group_1 = group_1,
      group_2 = group_2,
      comparison_dir = output_dir,
      file_prefix = file_prefix,
      error_message = final_error %||% "Analysis failed after optimization retries."
    )
    failed$optimization <- Go_BuildOptimizationSummary(
      attempt_plan = attempt_plan,
      final_attempt_id = length(attempt_plan),
      final_status = "failed",
      error_message = final_error %||% "Analysis failed after optimization retries."
    )
    utils::write.csv(
      failed$optimization,
      file.path(output_dir, paste0(file_prefix, ".optimization_status.csv")),
      row.names = FALSE
    )
    failed$exported$optimization <- file.path(output_dir, paste0(file_prefix, ".optimization_status.csv"))
    return(invisible(failed))
  }

  optimization <- Go_BuildOptimizationSummary(
    attempt_plan = attempt_plan,
    final_attempt_id = result$attempt_id,
    final_status = result$status,
    da_table = result$da_standardized$all_methods_standardized
  )
  utils::write.csv(
    optimization,
    file.path(output_dir, paste0(file_prefix, ".optimization_status.csv")),
    row.names = FALSE
  )
  result$optimization <- optimization
  result$exported$optimization <- file.path(output_dir, paste0(file_prefix, ".optimization_status.csv"))

  if (isTRUE(qc_plot)) {
    qc_dir <- Go_CreateQCPlotDir(dirname(output_dir), group_1, group_2)
    message("[ConDA] Generating QC plots.")
    qc_out <- tryCatch(
      Go_ConDaQCplot(
        result = result,
        output_dir = qc_dir,
        label_col = "ShortName"
      ),
      error = function(e) {
        message("[ConDA] QC plot generation failed: ", conditionMessage(e))
        NULL
      }
    )
    result$qc_plot <- qc_out
    if (!is.null(qc_out)) {
      result$exported$qc_plot <- qc_out$files
      result$exported$qc_plot_html <- qc_out$html_files
    }
  }

  if (isTRUE(volcano_plot)) {
    message("[ConDA] Exporting volcano bridge tables.")
    bridge_out <- tryCatch(
      Go_ExportVolcanoBridge(
        output_dir = volcano_bridge_root_dir %||% file.path(output_dir, "volcano_bridge"),
        da_table = result$da_standardized$all_methods_standardized,
        final_scores = result$final_scores,
        filtered_metadata = result$filtered$metadata,
        group_var = group_var,
        group_1 = group_1,
        group_2 = group_2,
        analysis_mode = unique(result$final_scores$analysis_mode %||% NA_character_),
        methods = methods,
        file_prefix = file_prefix,
        name = name
      ),
      error = function(e) {
        message("[ConDA] Volcano bridge export failed: ", conditionMessage(e))
        NULL
      }
    )
    result$volcano_bridge <- bridge_out
    if (!is.null(bridge_out)) {
      result$exported$volcano_bridge <- bridge_out$files
      result$exported$volcano_bridge_dir <- bridge_out$dir
      message("[ConDA] Generating volcano plots through Gotools.")
      volcano_out <- tryCatch(
        Go_RunVolcanoPlotBridge(
          project = project,
          bridge_dir = bridge_out$dir,
          comparison_name = paste0(group_1, ".vs.", group_2),
          name = name
        ),
        error = function(e) {
          message("[ConDA] Volcano plot generation failed: ", conditionMessage(e))
          NULL
        }
      )
      result$volcano_plot <- volcano_out
      if (!is.null(volcano_out)) {
        result$exported$volcano_plot_dir <- volcano_out$output_dir
      }
    }
  }

  invisible(list(
    comparison = paste0(group_1, ".vs.", group_2),
    comparison_dir = output_dir,
    file_prefix = file_prefix,
    status = result$status,
    error_message = result$error_message,
    filtered = result$filtered,
    input_bundle = input_bundle,
    da_raw = result$da_raw,
    da_standardized = result$da_standardized,
    beta_distances = result$beta_distances,
    beta_contribution = result$beta_contribution,
    da_consensus = result$da_consensus,
    method_annotation = result$method_annotation,
    final_scores = result$final_scores,
    qc_plot = result$qc_plot %||% NULL,
    volcano_bridge = result$volcano_bridge %||% NULL,
    volcano_plot = result$volcano_plot %||% NULL,
    optimization = optimization,
    exported = result$exported
  ))
}

Go_RunSingleDAAttempt <- function(
  input_bundle,
  group_var,
  group_1,
  group_2,
  random_effects = NULL,
  covariates = NULL,
  methods,
  distances,
  method_controls = NULL,
  prevalence = 0.1,
  abundance = 1e-4,
  alpha = 0.05,
  n_permutations = 999L,
  n_beta_permutations = 99L,
  weights = c(da = 0.4, beta = 0.3, direction = 0.15, effect = 0.15),
  output_dir,
  file_prefix = NULL
) {
  message("[ConDA] Aligning samples and filtering features.")
  aligned <- Go_AlignInputs(
    feature_table = input_bundle$feature_table,
    metadata = input_bundle$metadata
  )
  single_method_mode <- length(methods) == 1 && (is.null(distances) || length(distances) == 0)
  if (isTRUE(single_method_mode)) {
    message("[ConDA] Single-method mode: using aligned input without additional feature filtering.")
    filtered <- list(
      feature_table = aligned$feature_table,
      metadata = aligned$metadata,
      filter_summary = data.frame(
        n_features_input = nrow(aligned$feature_table),
        n_features_retained = nrow(aligned$feature_table),
        prevalence = NA_real_,
        abundance = NA_real_,
        filter_mode = "aligned_only",
        stringsAsFactors = FALSE
      )
    )
    if (length(methods) == 1 && identical(methods, "ancombc2")) {
      method_controls <- method_controls %||% list()
      method_controls$ancombc2 <- utils::modifyList(
        method_controls$ancombc2 %||% list(),
        list(legacy_filter_cutoff = 0.001)
      )
      message("[ConDA] Single ANCOM mode: applying Go_Ancom2-style legacy cutoff guidance (0.001).")
    }
    if (length(methods) == 1 && identical(methods, "corncob")) {
      method_controls <- method_controls %||% list()
      method_controls$corncob <- utils::modifyList(
        method_controls$corncob %||% list(),
        list(
          min_prevalence = 0.05,
          min_total_count = 10,
          phi_formula = "auto",
          phi_null_formula = "auto",
          retry_min_prevalence = 0.1,
          retry_min_total_count = 20,
          retry_legacy_filter_cutoff = 0.001,
          retry_phi_formula = "~ 1",
          retry_phi_null_formula = "~ 1"
        )
      )
      message("[ConDA] Single corncob mode: applying stability-focused filtering and retry settings.")
    }
  } else {
    filtered <- Go_FilterFeatures(
      feature_table = aligned$feature_table,
      metadata = aligned$metadata,
      prevalence = prevalence,
      abundance = abundance,
      output_dir = NULL
    )
  }

  if (nrow(filtered$feature_table) < 2) {
    stop("Too few features retained after filtering.")
  }

  message("[ConDA] Running DA methods: ", paste(methods, collapse = ", "))
  da_raw <- Go_RunDAmethods(
    feature_table = filtered$feature_table,
    metadata = filtered$metadata,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    random_effects = random_effects,
    covariates = covariates,
    methods = methods,
    method_controls = method_controls,
    alpha = alpha
  )

  da_standardized <- Go_StandardizeDA(
    da_results = da_raw,
    comparison = paste0(group_1, "_vs_", group_2),
    alpha = alpha
  )
  da_standardized$all_methods_standardized <- Go_AppendTaxonomy(
    da_standardized$all_methods_standardized,
    input_bundle$taxonomy
  )

  beta_enabled <- !is.null(distances) && length(distances) > 0
  if (beta_enabled) {
    message("[ConDA] Computing beta distances: ", paste(distances, collapse = ", "))
    beta_distances <- Go_BetaDistance(
      feature_table = filtered$feature_table,
      metadata = filtered$metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      distances = distances,
      phy_tree = input_bundle$phy_tree,
      n_permutations = n_permutations
    )

    message("[ConDA] Estimating beta contribution.")
    beta_contribution <- Go_BetaContribution(
      feature_table = filtered$feature_table,
      metadata = filtered$metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      beta_distances = beta_distances,
      phy_tree = input_bundle$phy_tree,
      n_beta_permutations = n_beta_permutations
    )
  } else {
    beta_distances <- Go_CreateEmptyBetaDistance()
    beta_contribution <- Go_CreateEmptyBetaContribution(rownames(filtered$feature_table))
  }

  message("[ConDA] Building consensus tables.")
  da_consensus <- Go_DAConsensus(
    da_table = da_standardized$all_methods_standardized,
    alpha = alpha
  )
  method_annotation <- Go_BuildMethodAnnotation(
    da_table = da_standardized$all_methods_standardized
  )

  final_scores <- Go_FinalScore(
    da_consensus = da_consensus,
    beta_contribution = beta_contribution,
    method_annotation = method_annotation,
    beta_enabled = beta_enabled,
    weights = weights
  )
  final_scores <- Go_AppendTaxonomy(final_scores, input_bundle$taxonomy)

  message("[ConDA] Exporting result tables.")
  exported <- Go_ExportResults(
    output_dir = output_dir,
    filtered_feature_table = filtered$feature_table,
    standardized_da = da_standardized$all_methods_standardized,
    da_consensus = da_consensus,
    beta_summary = beta_distances$beta_summary,
    beta_feature_contribution = beta_contribution,
    final_scores = final_scores,
    file_prefix = file_prefix
  )

  list(
    status = "completed",
    error_message = NA_character_,
    filtered = filtered,
    da_raw = da_raw,
    da_standardized = da_standardized,
    beta_distances = beta_distances,
    beta_contribution = beta_contribution,
    da_consensus = da_consensus,
    method_annotation = method_annotation,
    final_scores = final_scores,
    exported = exported
  )
}
