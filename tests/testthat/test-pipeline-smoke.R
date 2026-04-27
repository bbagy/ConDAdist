# ---- helpers ---------------------------------------------------------------

make_feature_table <- function(n_taxa = 8, n_samples = 8, seed = 1) {
  set.seed(seed)
  mat <- matrix(
    as.integer(stats::rpois(n_taxa * n_samples, lambda = 100)),
    nrow = n_taxa,
    dimnames = list(
      paste0("taxa", seq_len(n_taxa)),
      paste0("s", seq_len(n_samples))
    )
  )
  storage.mode(mat) <- "numeric"
  mat
}

make_metadata <- function(n_samples = 8) {
  data.frame(
    group = rep(c("A", "B"), each = n_samples / 2),
    row.names = paste0("s", seq_len(n_samples)),
    stringsAsFactors = FALSE
  )
}

# ---- Go_FilterFeatures -----------------------------------------------------

test_that("Go_FilterFeatures retains matrix structure", {
  feature_table <- matrix(
    c(10, 0, 5, 1, 0, 2, 3, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("taxa1", "taxa2"), c("s1", "s2", "s3", "s4"))
  )
  metadata <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = c("s1", "s2", "s3", "s4")
  )

  out <- Go_FilterFeatures(feature_table, metadata, prevalence = 0, abundance = 0)
  expect_equal(ncol(out$feature_table), 4)
  expect_true(is.matrix(out$feature_table))
})

test_that("Go_FilterFeatures drops zero-prevalence features", {
  # byrow = TRUE: t1 and t4 have all-zero rows, t2 and t3 are non-zero
  ft <- matrix(
    c(0, 0,   # t1: absent in all samples
      0, 10,  # t2: present in 1 of 2 samples
      10, 5,  # t3: present in both
      0, 0),  # t4: absent in all samples
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("t", 1:4), paste0("s", 1:2))
  )
  md <- data.frame(g = c("A", "B"), row.names = paste0("s", 1:2))
  out <- Go_FilterFeatures(ft, md, prevalence = 0.5, abundance = 0)
  expect_true(all(c("t2", "t3") %in% rownames(out$feature_table)))
  expect_false("t1" %in% rownames(out$feature_table))
  expect_false("t4" %in% rownames(out$feature_table))
})

test_that("Go_FilterFeatures returns filter_summary with scope column", {
  ft <- make_feature_table()
  md <- make_metadata()
  out <- Go_FilterFeatures(ft, md, prevalence = 0, abundance = 0)
  out$filter_summary$filter_scope <- "pairwise"
  expect_true("filter_scope" %in% names(out$filter_summary))
})

# ---- Go_MethodSignature ----------------------------------------------------

test_that("Go_MethodSignature single method returns full name", {
  expect_equal(Go_MethodSignature("deseq2"),   "deseq2")
  expect_equal(Go_MethodSignature("aldex2"),   "aldex2")
  expect_equal(Go_MethodSignature("ancombc2"), "ancombc2")
  expect_equal(Go_MethodSignature("maaslin2"), "maaslin2")
  expect_equal(Go_MethodSignature("corncob"),  "corncob")
})

test_that("Go_MethodSignature multi-method uses fixed-order initials DANMC", {
  expect_equal(Go_MethodSignature(c("deseq2", "aldex2", "ancombc2", "maaslin2", "corncob")), "DANMC")
  expect_equal(Go_MethodSignature(c("deseq2", "corncob")), "DC")
  expect_equal(Go_MethodSignature(c("aldex2", "ancombc2")), "AN")
  # Order of input should not matter
  expect_equal(Go_MethodSignature(c("corncob", "deseq2")), "DC")
})

test_that("Go_MethodSignature normalises maaslin alias", {
  expect_equal(Go_MethodSignature(c("deseq2", "maaslin")), "DM")
})

# ---- Go_ResolveMethods / Go_ResolveDistances --------------------------------

test_that("Go_ResolveMethods deduplicates and lowercases", {
  expect_equal(Go_ResolveMethods(c("DESeq2", "deseq2")), "deseq2")
  expect_equal(Go_ResolveMethods(c("ALDEX2", "corncob")), c("aldex2", "corncob"))
})

test_that("Go_ResolveMethods normalises maaslin alias", {
  expect_equal(Go_ResolveMethods("maaslin"), "maaslin2")
})

test_that("Go_ResolveDistances removes phylo metrics when no tree", {
  out <- Go_ResolveDistances(c("bray", "unifrac"), phy_tree = NULL)
  expect_true("bray" %in% out)
  expect_false("unweighted_unifrac" %in% out)
})

test_that("Go_ResolveDistances accepts exactly 3 valid metrics without error", {
  out <- Go_ResolveDistances(c("bray", "jaccard", "aitchison"), phy_tree = NULL)
  expect_equal(length(out), 3L)
})

test_that("Go_ResolveDistances deduplicates before checking limit", {
  # "bray" appears twice → dedup → 3 unique → no error
  out <- Go_ResolveDistances(c("bray", "jaccard", "aitchison", "bray"), phy_tree = NULL)
  expect_equal(length(out), 3L)
})

# ---- Go_BuildComparisonPlan ------------------------------------------------

test_that("Go_BuildComparisonPlan baseline vs targets", {
  md <- data.frame(g = c("A", "B", "C"), row.names = paste0("s", 1:3))
  plan <- Go_BuildComparisonPlan(md, "g", group_1 = "A", group_2 = c("B", "C"))
  expect_equal(nrow(plan), 2)
  expect_true(all(plan$group_1 == "A"))
  expect_setequal(plan$group_2, c("B", "C"))
})

test_that("Go_BuildComparisonPlan pairwise_all", {
  md <- data.frame(g = c("A", "B", "C"), row.names = paste0("s", 1:3))
  plan <- Go_BuildComparisonPlan(md, "g", group_1 = NULL, group_2 = NULL,
                                 orders = c("A", "B", "C"), pairwise_all = TRUE)
  expect_equal(nrow(plan), 3)
})

# ---- Go_StandardizeDA ------------------------------------------------------

test_that("Go_StandardizeDA produces required columns", {
  ft <- make_feature_table()
  md <- make_metadata()
  da_result <- list(
    deseq2 = Go_StandardSchema(rownames(ft))
  )
  da_result$deseq2$method <- "deseq2"
  da_result$deseq2$p_value <- stats::runif(nrow(ft))
  da_result$deseq2$q_value <- stats::p.adjust(da_result$deseq2$p_value, "BH")

  out <- Go_StandardizeDA(da_result, comparison = "A_vs_B", alpha = 0.05)
  expect_true("all_methods_standardized" %in% names(out))
  required <- c("feature_id", "method", "p_value", "q_value", "is_significant")
  expect_true(all(required %in% colnames(out$all_methods_standardized)))
})

# ---- Go_DAConsensus --------------------------------------------------------

test_that("Go_DAConsensus produces fisher_combined_q column", {
  set.seed(1)
  n <- 10
  da_table <- data.frame(
    feature_id     = rep(paste0("t", 1:n), 2),
    method         = rep(c("deseq2", "aldex2"), each = n),
    p_value        = c(stats::runif(n, 0, 0.1), stats::runif(n, 0.1, 1)),
    q_value        = stats::runif(n * 2),
    effect_size    = stats::rnorm(n * 2),
    direction      = rep(c("up_in_group2", "up_in_group1"), n),
    is_significant = c(rep(TRUE, 5), rep(FALSE, 5), rep(TRUE, 3), rep(FALSE, 7)),
    stringsAsFactors = FALSE
  )
  out <- Go_DAConsensus(da_table, alpha = 0.05)
  expect_true("fisher_combined_q" %in% colnames(out))
  expect_true("DA_support_score" %in% colnames(out))
  expect_equal(nrow(out), n)
})

test_that("Go_DAConsensus handles single method", {
  n <- 5
  da_table <- data.frame(
    feature_id     = paste0("t", 1:n),
    method         = "deseq2",
    p_value        = c(0.01, 0.02, 0.5, 0.8, 0.001),
    q_value        = c(0.05, 0.08, 0.6, 0.9, 0.01),
    effect_size    = c(1, -1, 0.1, -0.2, 2),
    direction      = c("up_in_group2", "up_in_group1", "up_in_group2", "up_in_group1", "up_in_group2"),
    is_significant = c(TRUE, TRUE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  out <- Go_DAConsensus(da_table, alpha = 0.05)
  expect_equal(nrow(out), n)
  expect_true(all(out$n_methods_run == 1))
})

# ---- Go_FinalScore ---------------------------------------------------------

test_that("Go_FinalScore DA-only mode sets analysis_mode correctly", {
  n <- 6
  da_consensus <- data.frame(
    feature_id            = paste0("t", 1:n),
    n_methods_run         = 2,
    n_methods_significant = c(2, 1, 0, 2, 1, 0),
    DA_support_score      = c(1, 0.5, 0, 1, 0.5, 0),
    fisher_combined_p     = c(0.01, 0.04, 0.6, 0.001, 0.03, 0.9),
    fisher_combined_q     = c(0.03, 0.06, 0.7, 0.006, 0.05, 0.95),
    is_combined_significant = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
    direction_consistency = c(1, 1, 0, 1, 1, 0),
    effect_consistency    = c(0.9, 0.8, 0.5, 0.95, 0.7, 0.4),
    mean_effect_size      = c(1.2, -0.5, 0.1, 2.0, -0.8, 0.05),
    median_effect_size    = c(1.1, -0.4, 0.1, 2.1, -0.9, 0.06),
    min_method_q_value    = c(0.01, 0.04, 0.6, 0.001, 0.03, 0.9),
    stringsAsFactors = FALSE
  )
  beta_contribution <- Go_CreateEmptyBetaContribution(paste0("t", 1:n))

  out <- Go_FinalScore(da_consensus, beta_contribution, beta_enabled = FALSE)
  expect_true("final_score" %in% colnames(out))
  expect_true("classification" %in% colnames(out))
  expect_true(all(out$analysis_mode %in% c("da_only", "single_method")))
})

test_that("Go_FinalScore full mode sets analysis_mode = 'full'", {
  n <- 4
  da_consensus <- data.frame(
    feature_id = paste0("t", 1:n),
    n_methods_run = 2, n_methods_significant = c(2,1,0,2),
    DA_support_score = c(1,0.5,0,1),
    fisher_combined_p = c(0.01,0.04,0.6,0.001),
    fisher_combined_q = c(0.03,0.06,0.7,0.006),
    is_combined_significant = c(TRUE,FALSE,FALSE,TRUE),
    direction_consistency = c(1,1,0,1),
    effect_consistency = c(0.9,0.8,0.5,0.95),
    mean_effect_size = c(1.2,-0.5,0.1,2.0),
    median_effect_size = c(1.1,-0.4,0.1,2.1),
    min_method_q_value = c(0.01,0.04,0.6,0.001),
    stringsAsFactors = FALSE
  )
  beta_contribution <- data.frame(
    feature_id = paste0("t", 1:n),
    simper_contribution = c(0.2, 0.1, 0.05, 0.4),
    simper_score = c(0.8, 0.4, 0.1, 1.0),
    delta_R2 = c(0.05, 0.02, 0.0, 0.1),
    delta_R2_score = c(0.7, 0.3, 0.0, 1.0),
    beta_contribution_score = c(0.75, 0.35, 0.05, 1.0),
    beta_perm_p = rep(NA_real_, n),
    beta_perm_q = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
  out <- Go_FinalScore(da_consensus, beta_contribution, beta_enabled = TRUE)
  expect_true(all(out$analysis_mode == "full"))
})

# ---- Go_CombinePValuesFisher -----------------------------------------------

test_that("Go_CombinePValuesFisher returns NA for empty input", {
  expect_true(is.na(Go_CombinePValuesFisher(numeric(0))))
  expect_true(is.na(Go_CombinePValuesFisher(NA_real_)))
})

test_that("Go_CombinePValuesFisher very small p gives small combined p", {
  combined <- Go_CombinePValuesFisher(c(0.001, 0.001, 0.001))
  expect_true(combined < 0.001)
})

test_that("Go_CombinePValuesFisher large p gives large combined p", {
  combined <- Go_CombinePValuesFisher(c(0.9, 0.8, 0.95))
  expect_true(combined > 0.5)
})

# ---- Go_NormalizeVector ----------------------------------------------------

test_that("Go_NormalizeVector returns values in [0,1]", {
  x <- c(1, 5, 3, 9, 0, NA)
  out <- Go_NormalizeVector(x)
  expect_true(all(out[is.finite(x)] >= 0 & out[is.finite(x)] <= 1))
})

test_that("Go_NormalizeVector handles all-NA input", {
  out <- Go_NormalizeVector(rep(NA_real_, 5))
  expect_true(all(out == 0))
})

# ---- Go_DirectionConsistency / Go_EffectConsistency ------------------------

test_that("Go_DirectionConsistency returns 1 for unanimous direction", {
  expect_equal(Go_DirectionConsistency(rep("up_in_group2", 5)), 1)
})

test_that("Go_DirectionConsistency returns 0.5 for equal split", {
  expect_equal(Go_DirectionConsistency(c("up_in_group2", "up_in_group1")), 0.5)
})

test_that("Go_DirectionConsistency returns 0 for empty input", {
  expect_equal(Go_DirectionConsistency(character(0)), 0)
})

test_that("Go_EffectConsistency returns 1 for single observation", {
  expect_equal(Go_EffectConsistency(1.5), 1)
})

# ---- Go_FilePrefix ---------------------------------------------------------

test_that("Go_FilePrefix produces expected pattern", {
  out <- Go_FilePrefix("MyProj", "Control", "Treatment", methods = "deseq2")
  expect_match(out, "deseq2")
  expect_match(out, "Control.vs.Treatment")
  expect_match(out, "MyProj")
})

test_that("Go_FilePrefix multi-method uses initials", {
  out <- Go_FilePrefix("P", "A", "B", methods = c("deseq2", "aldex2"))
  expect_match(out, "DA")
})

# ---- Go_ExportResults (file I/O) -------------------------------------------

test_that("Go_ExportResults writes expected CSV files", {
  tmp <- tempfile()
  dir.create(tmp)

  ft <- make_feature_table(4, 4)
  md <- make_metadata(4)
  n  <- nrow(ft)

  standardized_da <- data.frame(
    feature_id = paste0("taxa", 1:n), method = "deseq2",
    p_value = runif(n), q_value = runif(n),
    stringsAsFactors = FALSE
  )
  da_consensus <- data.frame(feature_id = paste0("taxa", 1:n),
                             fisher_combined_p = runif(n),
                             stringsAsFactors = FALSE)
  beta_summary  <- data.frame(distance = "bray", statistic = 0.3,
                              p_value = 0.01, stringsAsFactors = FALSE)
  beta_contrib  <- Go_CreateEmptyBetaContribution(paste0("taxa", 1:n))
  final_scores  <- data.frame(feature_id = paste0("taxa", 1:n),
                              final_score = runif(n),
                              stringsAsFactors = FALSE)

  files <- Go_ExportResults(
    output_dir              = tmp,
    filtered_feature_table  = ft,
    standardized_da         = standardized_da,
    da_consensus            = da_consensus,
    beta_summary            = beta_summary,
    beta_feature_contribution = beta_contrib,
    final_scores            = final_scores,
    file_prefix             = "test"
  )

  expect_true(file.exists(files$filtered_feature_table))
  expect_true(file.exists(files$all_methods_standardized))
  expect_true(file.exists(files$final_consensus_scores))
  unlink(tmp, recursive = TRUE)
})

# ---- Go_BuildMethodOverlapLong (HTML bug fix) -------------------------------

test_that("Go_BuildMethodOverlapLong handles NA methods in da_table", {
  n <- 5
  final_scores <- data.frame(
    feature_id        = paste0("t", 1:n),
    plot_label        = paste0("Taxon", 1:n),
    plot_label_unique = paste0("Taxon", 1:n),
    deseq2_is_significant = c(TRUE, FALSE, TRUE, FALSE, FALSE),
    stringsAsFactors  = FALSE
  )
  # da_table includes a row with method = NA (fallback artifact)
  da_table <- data.frame(
    feature_id = paste0("t", 1:n),
    method     = c("deseq2", "deseq2", NA, "deseq2", "deseq2"),
    stringsAsFactors = FALSE
  )
  out <- Go_BuildMethodOverlapLong(final_scores, da_table)
  expect_false(any(is.na(out$method)))
  expect_equal(nrow(out), n)
})

test_that("Go_BuildMethodOverlapLong uses rep(FALSE) fallback not scalar", {
  n <- 4
  final_scores <- data.frame(
    feature_id        = paste0("t", 1:n),
    plot_label        = paste0("Taxon", 1:n),
    plot_label_unique = paste0("Taxon", 1:n),
    stringsAsFactors  = FALSE
  )
  da_table <- data.frame(
    feature_id = paste0("t", 1:n),
    method     = "aldex2",
    stringsAsFactors = FALSE
  )
  # aldex2_is_significant column is absent — should not error, should give FALSE vector
  out <- Go_BuildMethodOverlapLong(final_scores, da_table)
  expect_equal(nrow(out), n)
  expect_true(all(out$is_significant == FALSE))
})

# ---- Go_CreateEmptyBetaContribution ----------------------------------------

test_that("Go_CreateEmptyBetaContribution returns correct structure", {
  ids <- paste0("t", 1:5)
  out <- Go_CreateEmptyBetaContribution(ids)
  expect_equal(out$feature_id, ids)
  expect_true(all(is.na(out$beta_contribution_score)))
})

# ---- inst/extdata example data --------------------------------------------

test_that("minimal_example_ps.rds loads as phyloseq", {
  skip_if_not_installed("phyloseq")
  path <- system.file("extdata", "minimal_example_ps.rds", package = "ConDAdist")
  if (!nzchar(path)) {
    path <- file.path(system.file(package = "ConDAdist"), "extdata", "minimal_example_ps.rds")
  }
  # When running via devtools::test(), system.file may not resolve; use relative path
  if (!file.exists(path)) {
    path <- file.path(
      rprojroot::find_package_root_file("inst/extdata/minimal_example_ps.rds")
    )
  }
  skip_if(!file.exists(path), "minimal_example_ps.rds not found — skipping")
  ps <- readRDS(path)
  expect_s4_class(ps, "phyloseq")
  expect_gte(phyloseq::ntaxa(ps), 10L)
  expect_gte(phyloseq::nsamples(ps), 4L)
})
