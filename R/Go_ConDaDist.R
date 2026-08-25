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
#' Comparison modes:
#'
#' - `pairwise_all = FALSE`: run `group_1` against each value in `group_2`
#' - `pairwise_all = TRUE`: ignore `group_1` and `group_2`, and run all
#'   pairwise contrasts using `orders`
#'
#' @param psIN Standard `phyloseq` object used in the Go tool workflow.
#' @param group_var Column in `sample_data(psIN)` that defines the comparison
#'   groups.
#' @param group_1 Baseline comparison group. Can be `NULL` when
#'   `pairwise_all = TRUE`.
#' @param group_2 One or more target groups compared against `group_1`. Can be
#'   `NULL` when `pairwise_all = TRUE`.
#' @param orders Optional ordered levels for `group_var`. When
#'   `pairwise_all = TRUE`, all pairwise comparisons follow this order.
#' @param project Project name used to create the dated output directory.
#' @param covariates Optional character vector of metadata covariates.
#' @param name Optional label appended to output file names.
#' @param random_effects Optional random-effect metadata variables used for
#'   native `MaAsLin2` single-mode runs and `ANCOMBC2` mixed-effects models.
#' @param methods Differential abundance methods to run. Defaults to all
#'   supported methods. Available options:
#'   \itemize{
#'     \item \code{"deseq2"}       — DESeq2 (negative binomial, poscounts normalization)
#'     \item \code{"aldex2"}       — ALDEx2
#'     \item \code{"ancombc2"}     — ANCOM-BC2 (log-linear, bias-corrected)
#'     \item \code{"corncob_lrt"}  — corncob likelihood-ratio test
#'     \item \code{"corncob_wald"} — corncob Wald test (paired with LRT as one family)
#'   }
#' @param distances Beta-diversity distances to compute. At most 3 distances
#'   are allowed per run. Set to \code{NULL} to disable beta-diversity and run
#'   DA-only mode. Available options:
#'   \itemize{
#'     \item \code{"bray"}             — Bray-Curtis dissimilarity (recommended default)
#'     \item \code{"jaccard"}          — Jaccard distance (presence/absence)
#'     \item \code{"jsd"}              — Jensen-Shannon distance
#'     \item \code{"unifrac"}          — Unweighted UniFrac (requires phylogenetic tree)
#'     \item \code{"weighted_unifrac"} — Weighted UniFrac (requires phylogenetic tree)
#'   }
#'   Recommended without tree: \code{c("bray", "jaccard", "jsd")}.
#'   Recommended with tree: \code{c("bray", "jsd", "unifrac")}.
#' @param method_controls Optional named list of per-method control lists.
#' @param prevalence Minimum fraction of samples with non-zero abundance used in
#'   feature filtering.
#' @param abundance Minimum mean relative abundance threshold used in feature
#'   filtering.
#' @param alpha Significance cutoff used in consensus summaries.
#' @param n_permutations Number of PERMANOVA permutations.
#' @param n_beta_permutations Number of feature-level beta permutations.
#' @param weights Named numeric vector controlling final score weights.
#' @param p_combine DA p-value combination rule used in the consensus layer.
#'   One of \code{"family_partial_conjunction"} (default; V5 all-but-one
#'   method-family agreement), \code{"adaptive_cauchy"} (legacy,
#'   exploratory only), or \code{"fisher"}.
#'   When exactly one method is requested and this argument is omitted, CDD
#'   retains the prior V2 single-method skeleton because no cross-method
#'   combination is performed.
#'   \code{"cauchy"} is accepted for backward compatibility but is deprecated
#'   and mapped to \code{"adaptive_cauchy"}.
#' @param preset Configuration contract. \code{"broad_panel"} is the conservative
#'   general default (five tests, no distance); \code{"custom"} accepts explicit
#'   component arguments and is not
#'   independently calibrated. Supplying a component without `preset`
#'   automatically selects `custom` for backward compatibility.
#' @param qc_plot Generate QC plots automatically after analysis.
#' Volcano bridge tables and Gotools volcano plots are generated automatically.
#' @param pairwise_all If `TRUE`, ignore the baseline-vs-target pattern and run
#'   all pairwise comparisons across ordered `group_var` levels.
#' @param continue_on_error Keep running remaining pairwise comparisons even if
#'   one comparison fails.
#' @param filter_scope Whether prevalence and abundance filters are computed on
#'   each pairwise subset (`"pairwise"`) or on the full aligned data (`"global"`).
#'
#' @return Volcano/lollipop bridge table directory path. Single-method mode
#'   returns `table/ConDaDist_plot_single_Tab`; consensus mode returns
#'   `table/ConDaDist_plot_Tab`.
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
#'   methods = c("deseq2", "aldex2", "ancombc2", "corncob_wald", "corncob_lrt"),
#'   distances = c("bray", "jaccard", "jsd")
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
  group_1 = NULL,
  group_2 = NULL,
  orders = NULL,
  project,
  covariates = NULL,
  name = NULL,
  random_effects = NULL,
  # Solution1 DA tests: "deseq2", "aldex2", "ancombc2",
  # "corncob_wald", "corncob_lrt"
  methods = NULL,
  # distances: "bray", "jaccard", "jsd", "unifrac", "weighted_unifrac"
  # unifrac / weighted_unifrac require a phylogenetic tree in psIN
  # set NULL to run DA-only (no beta-diversity)
  distances = NULL,
  method_controls = NULL,
  prevalence = 0.1,
  abundance = 1e-4,
  alpha = 0.05,
  n_permutations = 999L,
  n_beta_permutations = 99L,
  weights = NULL,
  p_combine = NULL,
  preset = c("broad_panel", "custom"),
  qc_plot = TRUE,
  pairwise_all = FALSE,
  continue_on_error = TRUE,
  filter_scope = "pairwise"
) {
  supplied <- list(
    methods = !missing(methods), distances = !missing(distances),
    weights = !missing(weights), p_combine = !missing(p_combine)
  )
  if (missing(preset) && any(unlist(supplied))) preset <- "custom"
  config <- Go_ResolveCDDPreset(
    preset = preset, methods = methods, distances = distances,
    weights = weights, p_combine = p_combine, supplied = supplied
  )
  methods <- config$methods
  distances <- config$distances
  weights <- config$weights
  p_combine <- config$p_combine
  p_combine <- Go_ResolvePCombine(p_combine)
  result <- Go_RunConDaDistMain(
    feature_table = psIN,
    metadata = NULL,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    orders = orders,
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
    p_combine = p_combine,
    qc_plot = qc_plot,
    pairwise_all = pairwise_all,
    continue_on_error = continue_on_error,
    filter_scope = filter_scope
  )
  root_dir <- if (is.list(result) && !is.null(result$return_dir)) result$return_dir else result
  invisible(root_dir)
}

Go_RunConDaDistMain <- function(
  feature_table,
  metadata = NULL,
  group_var,
  group_1 = NULL,
  group_2 = NULL,
  orders = NULL,
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
  p_combine = c("family_partial_conjunction", "adaptive_cauchy", "fisher", "cauchy"),
  qc_plot = TRUE,
  pairwise_all = FALSE,
  continue_on_error = TRUE,
  filter_scope = "pairwise"
) {
  p_combine <- Go_ResolvePCombine(p_combine)
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
    output_layout <- Go_path(project = project, pdf = "no", table = "yes")
    bridge_out <- Go_ExportNativeMaaslin2VolcanoBridge(
      native_dir = native_dir,
      bridge_dir = output_layout$conda_dist_single_volcano,
      psIN = feature_table,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      name = name
    )
    return(invisible(normalizePath(output_layout$conda_dist_single_volcano, winslash = "/", mustWork = FALSE)))
  }

  output_layout <- Go_path(project = project, pdf = "no", table = "yes")
  root_output_dir <- output_layout$main

  comparison_plan <- Go_BuildComparisonPlan(
    metadata = normalized_bundle$metadata,
    group_var = group_var,
    group_1 = group_1,
    group_2 = group_2,
    orders = orders,
    pairwise_all = pairwise_all
  )

  comparison_results <- lapply(seq_len(nrow(comparison_plan)), function(i) {
    group_1_target <- comparison_plan$group_1[i]
    group_2_target <- comparison_plan$group_2[i]
    comparison_dir <- Go_CreateComparisonDir(output_layout$conda_dist, methods, group_1_target, group_2_target)
    file_prefix <- Go_FilePrefix(project = project, group_1 = group_1_target, group_2 = group_2_target, name = name, methods = methods)
    message("[ConDA] Starting comparison: ", group_1_target, " vs ", group_2_target)

    run_one <- function() {
      Go_RunSingleDAensemble(
        feature_table = feature_table,
        metadata = metadata,
        group_var = group_var,
        group_1 = group_1_target,
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
        p_combine = p_combine,
        qc_plot = qc_plot,
        volcano_bridge_root_dir = output_layout$conda_dist_volcano,
        single_volcano_bridge_root_dir = output_layout$conda_dist_single_volcano,
        output_dir = comparison_dir,
        file_prefix = file_prefix,
        filter_scope = filter_scope
      )
    }

    if (!isTRUE(continue_on_error)) {
      return(run_one())
    }

    tryCatch(
      run_one(),
      error = function(e) {
        Go_BuildFailedComparison(
          group_1 = group_1_target,
          group_2 = group_2_target,
          comparison_dir = comparison_dir,
          file_prefix = file_prefix,
          error_message = conditionMessage(e)
        )
      }
    )
  })
  names(comparison_results) <- paste0(comparison_plan$group_1, ".vs.", comparison_plan$group_2)

  method_dir <- dirname(comparison_results[[1]]$comparison_dir)
  return_dir <- if (length(methods) == 1 && length(distances) == 0) {
    output_layout$conda_dist_single_volcano
  } else {
    output_layout$conda_dist_volcano
  }

  if (length(comparison_results) == 1) {
    single <- comparison_results[[1]]
    single$project <- project
    single$output_root_dir <- root_output_dir
    single$output_table_dir <- output_layout$table
    single$return_dir <- normalizePath(return_dir, winslash = "/", mustWork = FALSE)
    return(invisible(single))
  }

  invisible(list(
    project = project,
    output_root_dir = root_output_dir,
    output_table_dir = output_layout$table,
    comparisons = comparison_results,
    return_dir = normalizePath(return_dir, winslash = "/", mustWork = FALSE)
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
  p_combine = c("family_partial_conjunction", "adaptive_cauchy", "fisher", "cauchy"),
  qc_plot = TRUE,
  volcano_bridge_root_dir = NULL,
  single_volcano_bridge_root_dir = NULL,
  output_dir,
  file_prefix = NULL,
  filter_scope = "pairwise"
) {
  p_combine <- Go_ResolvePCombine(p_combine)
  input_bundle <- Go_NormalizeInputBundle(
    feature_table = feature_table,
    metadata = metadata
  )

  methods <- Go_ResolveMethods(methods)
  distances <- Go_ResolveDistances(distances, phy_tree = input_bundle$phy_tree)
  single_method_mode <- length(methods) == 1 && (is.null(distances) || length(distances) == 0)

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
        p_combine = p_combine,
        output_dir = output_dir,
        file_prefix = file_prefix,
        filter_scope = filter_scope
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
    result$file_prefix  <- file_prefix
    result$input_bundle <- input_bundle
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

  message("[ConDA] Exporting volcano bridge tables.")
  bridge_out <- tryCatch(
    Go_ExportVolcanoBridge(
      output_dir = volcano_bridge_root_dir %||% file.path(output_dir, "volcano_bridge"),
      single_output_dir = single_volcano_bridge_root_dir %||% file.path(output_dir, "volcano_bridge"),
      da_table = result$da_standardized$all_methods_standardized,
      final_scores = result$final_scores,
      filtered_metadata = result$filtered$metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      analysis_mode = unique(result$final_scores$analysis_mode %||% NA_character_),
      methods = methods,
      distances = distances,
      write_consensus = !single_method_mode,
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
  p_combine = c("family_partial_conjunction", "adaptive_cauchy", "fisher", "cauchy"),
  output_dir,
  file_prefix = NULL,
  filter_scope = "pairwise"
) {
  p_combine <- Go_ResolvePCombine(p_combine)
  Go_ValidateFamilyPartialConjunctionPanel(methods, p_combine)
  v4_mode <- Go_ResolveV4Mode(p_combine)
  consensus_skeleton <- v4_mode$consensus_skeleton
  beta_contribution <- v4_mode$beta_contribution
  message("[TIMING] filter_start ", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"))
  message("[ConDA] Aligning samples and filtering features.")
  aligned <- Go_AlignInputs(
    feature_table = input_bundle$feature_table,
    metadata = input_bundle$metadata
  )
  single_method_mode <- length(methods) == 1 && (is.null(distances) || length(distances) == 0)

  filter_input <- aligned
  if (identical(filter_scope, "pairwise")) {
    pairwise_mask <- as.character(filter_input$metadata[[group_var]]) %in% c(group_1, group_2)
    filter_input$feature_table <- filter_input$feature_table[, pairwise_mask, drop = FALSE]
    filter_input$metadata <- filter_input$metadata[pairwise_mask, , drop = FALSE]
  }
  filtered <- Go_FilterFeatures(
    feature_table = filter_input$feature_table,
    metadata = filter_input$metadata,
    prevalence = prevalence,
    abundance = abundance,
    output_dir = NULL
  )
  filtered$filter_summary$filter_scope <- filter_scope

  if (isTRUE(single_method_mode)) {
    message("[ConDA] Single-method mode: filter_scope = \"", filter_scope, "\" (",
            nrow(filtered$feature_table), " features retained).")
    if (length(methods) == 1 && methods[[1]] %in% c("corncob", "corncob_wald", "corncob_lrt")) {
      corncob_key <- if (identical(methods[[1]], "corncob")) "corncob_wald" else methods[[1]]
      method_controls <- method_controls %||% list()
      method_controls[[corncob_key]] <- utils::modifyList(
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
        ),
        method_controls[[corncob_key]] %||% list()
      )
      message("[ConDA] Single ", corncob_key, " mode: applying stability-focused filtering and retry settings.")
    }
  }

  if (nrow(filtered$feature_table) < 2) {
    stop("Too few features retained after filtering.")
  }

  message("[TIMING] da_methods_start ", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"))
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
    message("[TIMING] da_methods_end/beta_distance_start ", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"))
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

    message("[TIMING] beta_distance_end/beta_contribution_start ", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"))
    message("[ConDA] Estimating beta contribution (method = ", beta_contribution, ").")
    beta_contribution_tbl <- Go_BetaContribution(
      feature_table = filtered$feature_table,
      metadata = filtered$metadata,
      group_var = group_var,
      group_1 = group_1,
      group_2 = group_2,
      beta_distances = beta_distances,
      phy_tree = input_bundle$phy_tree,
      n_beta_permutations = n_beta_permutations,
      method = beta_contribution
    )
  } else {
    beta_distances <- Go_CreateEmptyBetaDistance()
    beta_contribution_tbl <- Go_CreateEmptyBetaContribution(rownames(filtered$feature_table))
  }

  message("[TIMING] beta_contribution_end/consensus_start ", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"))
  message("[ConDA] Building consensus tables (p_combine = ", p_combine, ").")
  da_consensus <- Go_DAConsensus(
    da_table = da_standardized$all_methods_standardized,
    alpha = alpha,
    p_combine = p_combine,
    planned_methods = methods
  )
  method_annotation <- Go_BuildMethodAnnotation(
    da_table = da_standardized$all_methods_standardized
  )

  final_scores <- Go_FinalScore(
    da_consensus = da_consensus,
    beta_contribution = beta_contribution_tbl,
    method_annotation = method_annotation,
    beta_enabled = beta_enabled,
    weights = weights
  )
  final_scores <- Go_AppendTaxonomy(final_scores, input_bundle$taxonomy)

  message("[TIMING] consensus_end/export_start ", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"))
  message("[ConDA] Exporting result tables.")
  exported <- Go_ExportResults(
    output_dir = output_dir,
    filtered_feature_table = filtered$feature_table,
    standardized_da = da_standardized$all_methods_standardized,
    da_consensus = da_consensus,
    beta_summary = beta_distances$beta_summary,
    beta_feature_contribution = beta_contribution_tbl,
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
    beta_contribution = beta_contribution_tbl,
    da_consensus = da_consensus,
    method_annotation = method_annotation,
    final_scores = final_scores,
    exported = exported
  )
}
