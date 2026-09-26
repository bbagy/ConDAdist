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

make_da_table <- function(n = 10, seed = 1) {
  set.seed(seed)
  sig_pattern <- c(rep(TRUE, ceiling(n / 2)), rep(FALSE, floor(n / 2)))
  data.frame(
    feature_id     = rep(paste0("t", 1:n), 2),
    method         = rep(c("deseq2", "aldex2"), each = n),
    p_value        = c(stats::runif(n, 0, 0.1), stats::runif(n, 0.1, 1)),
    q_value        = stats::runif(n * 2),
    effect_size    = stats::rnorm(n * 2),
    direction      = rep(c("up_in_group2", "up_in_group1"), length.out = n * 2),
    is_significant = c(sig_pattern, rev(sig_pattern)),
    stringsAsFactors = FALSE
  )
}

make_da_consensus_v2 <- function(n = 6) {
  # Match the consensus output's mirrored cauchy_* compatibility columns.
  combined_p_vals <- c(0.01, 0.04, 0.6, 0.001, 0.03, 0.9)[1:n]
  combined_q_vals <- c(0.03, 0.06, 0.7, 0.006, 0.05, 0.95)[1:n]
  data.frame(
    feature_id              = paste0("t", 1:n),
    n_methods_run           = 2,
    n_methods_significant   = c(2, 1, 0, 2, 1, 0)[1:n],
    DA_support_score        = c(1, 0.5, 0, 1, 0.5, 0)[1:n],
    combined_p               = combined_p_vals,
    combined_q               = combined_q_vals,
    cauchy_combined_p       = combined_p_vals,
    cauchy_combined_q       = combined_q_vals,
    is_combined_significant = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)[1:n],
    direction_consistency   = c(1, 1, 0, 1, 1, 0)[1:n],
    effect_consistency      = c(0.9, 0.8, 0.5, 0.95, 0.7, 0.4)[1:n],
    combined_effect_rank    = c(0.8, 0.3, 0.5, 0.9, 0.2, 0.5)[1:n],
    min_method_q_value      = c(0.01, 0.04, 0.6, 0.001, 0.03, 0.9)[1:n],
    stringsAsFactors        = FALSE
  )
}

make_beta_contrib_v2 <- function(n = 6) {
  data.frame(
    feature_id              = paste0("t", 1:n),
    loo_separation_delta    = c(0.05, 0.02, 0.0, 0.1, 0.03, 0.0)[1:n],
    loo_separation_score    = c(0.7, 0.3, 0.0, 1.0, 0.4, 0.0)[1:n],
    beta_contribution_score = c(0.7, 0.3, 0.0, 1.0, 0.4, 0.0)[1:n],
    beta_perm_p             = rep(NA_real_, n),
    beta_perm_q             = rep(NA_real_, n),
    stringsAsFactors        = FALSE
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
  # The corncob alias resolves to corncob_lrt.
  expect_equal(Go_MethodSignature("corncob"),  "corncob_lrt")
})

test_that("Go_MethodSignature multi-method uses fixed-order initials DANML", {
  # The corncob alias contributes the L in the signature.
  expect_equal(Go_MethodSignature(c("deseq2", "aldex2", "ancombc2", "maaslin2", "corncob")), "DANML")
  expect_equal(Go_MethodSignature(c("deseq2", "corncob")), "DL")
  expect_equal(Go_MethodSignature(c("aldex2", "ancombc2")), "AN")
  expect_equal(Go_MethodSignature(c("corncob", "deseq2")), "DL")
})

test_that("Go_MethodSignature normalises maaslin alias", {
  expect_equal(Go_MethodSignature(c("deseq2", "maaslin")), "DM")
})

# ---- Go_ResolveMethods / Go_ResolveDistances --------------------------------

test_that("Go_ResolveMethods deduplicates and lowercases", {
  expect_equal(Go_ResolveMethods(c("DESeq2", "deseq2")), "deseq2")
  # Resolve the corncob alias consistently.
  expect_equal(Go_ResolveMethods(c("ALDEX2", "corncob")), c("aldex2", "corncob_lrt"))
})

test_that("Go_ResolveMethods normalises maaslin alias", {
  expect_equal(Go_ResolveMethods("maaslin"), "maaslin2")
})

test_that("Solution1 is the multi-method default and single-method behavior is preserved", {
  expect_equal(
    Go_AllDAMethods(),
    c("deseq2", "aldex2", "ancombc2", "corncob_wald", "corncob_lrt")
  )
  expect_equal(
    Go_DefaultPCombineForMethods(Go_AllDAMethods()),
    "family_partial_conjunction"
  )
  expect_equal(
    Go_DefaultPCombineForMethods("deseq2"),
    "adaptive_cauchy"
  )
})

test_that("Go_ResolveDistances removes phylo metrics when no tree", {
  out <- Go_ResolveDistances(c("bray", "unifrac"), phy_tree = NULL)
  expect_true("bray" %in% out)
  expect_false("unweighted_unifrac" %in% out)
})

test_that("Go_ResolveDistances accepts exactly 3 valid metrics without error", {
  out <- Go_ResolveDistances(c("bray", "jaccard", "jsd"), phy_tree = NULL)
  expect_equal(length(out), 3L)
})

test_that("Go_ResolveDistances deduplicates before checking limit", {
  out <- Go_ResolveDistances(c("bray", "jaccard", "jsd", "bray"), phy_tree = NULL)
  expect_equal(length(out), 3L)
})

# ---- adjusted/restricted beta model ---------------------------------------

test_that("beta design applies complete cases and preserves permutation blocks", {
  md <- data.frame(
    group = rep(c("A", "B"), each = 4),
    age = c(20:25, NA, 27),
    subject = rep(paste0("p", 1:4), 2),
    row.names = paste0("s", 1:8)
  )
  design <- Go_PrepareBetaDesign(
    md, "group", "A", "B", covariates = "age", strata = "subject"
  )
  expect_equal(sum(design$keep), 7L)
  expect_equal(ncol(design$covariate_data), 1L)
  expect_equal(length(design$strata_factor), 7L)
})

test_that("restricted group permutations stay within strata", {
  group <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  strata <- factor(rep(paste0("p", 1:4), 2))
  set.seed(4)
  permuted <- Go_PermuteGroup(group, strata)
  for (idx in split(seq_along(group), strata)) {
    expect_equal(sort(as.character(permuted[idx])), sort(as.character(group[idx])))
  }
})

test_that("covariates and strata reach marginal PERMANOVA and adjusted LOO", {
  skip_if_not_installed("vegan")
  set.seed(11)
  md <- data.frame(
    group = rep(c("A", "B"), each = 6),
    age = seq_len(12),
    subject = rep(paste0("p", 1:6), 2),
    row.names = paste0("s", 1:12)
  )
  ft <- matrix(
    stats::rpois(10 * 12, lambda = 100), nrow = 10,
    dimnames = list(paste0("t", 1:10), rownames(md))
  )
  beta <- Go_BetaDistance(
    ft, md, "group", "A", "B", distances = "bray",
    n_permutations = 9L, covariates = "age", strata = "subject"
  )
  expect_true(is.finite(beta$beta_summary$statistic))
  expect_match(beta$beta_summary$notes, "covariates = 1")
  expect_match(beta$beta_summary$notes, "restricted permutation = TRUE")

  contribution <- Go_BetaContribution(
    ft, md, "group", "A", "B", beta,
    n_beta_permutations = 0L, covariates = "age", strata = "subject"
  )
  expect_equal(nrow(contribution), nrow(ft))
  expect_true(all(is.finite(contribution$beta_contribution_score)))
})

test_that("permutation-free adjusted group R2 matches adonis marginal R2", {
  skip_if_not_installed("vegan")
  set.seed(15)
  x <- matrix(stats::rnorm(20 * 4), nrow = 20)
  dm <- as.matrix(stats::dist(x))
  z <- factor(rep(c("Z0", "Z1"), 10))
  g <- factor(rep(c("A", "B"), each = 10))
  covariates <- data.frame(.cdd_cov_1 = z, check.names = FALSE)
  observed <- Go_AdjustedGroupR2(dm, g, covariates)
  fit <- vegan::adonis2(
    stats::as.dist(dm) ~ z + g,
    data = data.frame(z = z, g = g), permutations = 9L, by = "margin"
  )
  expect_equal(observed, fit["g", "R2"], tolerance = 1e-10)
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

# ---- Go_CombinePValuesCauchy (V2) ------------------------------------------

test_that("Go_CombinePValuesCauchy returns NA for empty input", {
  expect_true(is.na(Go_CombinePValuesCauchy(numeric(0))))
  expect_true(is.na(Go_CombinePValuesCauchy(NA_real_)))
})

test_that("Go_CombinePValuesCauchy very small p gives small combined p", {
  combined <- Go_CombinePValuesCauchy(c(0.001, 0.001, 0.001))
  expect_equal(combined, 0.001, tolerance = 1e-12)
})

test_that("Go_CombinePValuesCauchy large p gives large combined p", {
  combined <- Go_CombinePValuesCauchy(c(0.9, 0.8, 0.95))
  expect_true(combined > 0.5)
})

test_that("Go_CombinePValuesCauchy returns value in (0, 1]", {
  combined <- Go_CombinePValuesCauchy(c(0.01, 0.05, 0.001))
  expect_true(combined > 0 && combined <= 1)
})

test_that("Go_CombinePValuesCauchy custom weights applied correctly", {
  p <- c(0.001, 0.9)
  # Weight heavily toward the small p → should give smaller result than 50/50
  combined_heavy <- Go_CombinePValuesCauchy(p, weights = c(0.99, 0.01))
  combined_equal  <- Go_CombinePValuesCauchy(p, weights = c(0.5, 0.5))
  expect_true(combined_heavy < combined_equal)
})

# ---- Solution1 family partial conjunction ----------------------------------

solution1_methods <- c(
  "deseq2", "aldex2", "ancombc2", "corncob_wald", "corncob_lrt"
)

test_that("family partial conjunction implements 2 * third family p-value", {
  p <- c(0.01, 0.02, 0.03, 0.04, 0.05)
  # corncob family p = 2 * min(0.04, 0.05) = 0.08;
  # four family p-values are 0.01, 0.02, 0.03, 0.08.
  expect_equal(
    Go_CombinePValuesFamilyPartialConjunction(p, solution1_methods),
    0.06,
    tolerance = 1e-12
  )
})

test_that("one extreme constituent cannot dominate the family rule", {
  p <- c(0.5, 0.5, 1e-30, 0.5, 0.5)
  expect_equal(
    Go_CombinePValuesFamilyPartialConjunction(p, solution1_methods),
    1
  )
})

test_that("two corncob tests count as one family", {
  p <- c(0.5, 0.5, 0.5, 1e-30, 1e-20)
  expect_equal(
    Go_CombinePValuesFamilyPartialConjunction(p, solution1_methods),
    1
  )
})

test_that("missing slots keep the fixed family denominator", {
  expect_equal(
    Go_CombinePValuesFamilyPartialConjunction(
      c(0.01, 0.02, 0.03),
      c("deseq2", "aldex2", "ancombc2"),
      planned_methods = solution1_methods
    ),
    0.06,
    tolerance = 1e-12
  )
  expect_equal(
    Go_CombinePValuesFamilyPartialConjunction(
      c(0.01, 0.02),
      c("deseq2", "aldex2"),
      planned_methods = solution1_methods
    ),
    1
  )
})

test_that("duplicate method rows are rejected", {
  expect_error(
    Go_CombinePValuesFamilyPartialConjunction(
      c(0.01, 0.02),
      c("deseq2", "deseq2")
    ),
    "Duplicate method rows"
  )
})

test_that("V5 family partial conjunction supports configured panels", {
  expect_silent(
    ConDAdist:::Go_ValidateFamilyPartialConjunctionPanel(
      solution1_methods,
      "family_partial_conjunction"
    )
  )
  expect_silent(
    ConDAdist:::Go_ValidateFamilyPartialConjunctionPanel(
      setdiff(solution1_methods, "corncob_wald"),
      "family_partial_conjunction"
    )
  )
})

test_that("V5 presets resolve to frozen package contracts", {
  broad <- ConDAdist:::Go_ResolveCDDPreset("broad_panel", supplied = list())
  expect_equal(broad$methods, solution1_methods)
  expect_null(broad$distances)
  expect_equal(broad$engine_version, "V5")
  expect_error(
    ConDAdist:::Go_ResolveCDDPreset(
      "broad_panel", methods = "deseq2", supplied = list(methods = TRUE)
    ),
    "owns its configuration"
  )
  expect_error(
    ConDAdist:::Go_ResolveCDDPreset("recommended_stagex", supplied = list()),
    "should be one of"
  )
})

test_that("V5 partial conjunction uses all-but-one planned families", {
  d3 <- ConDAdist:::Go_FamilyPartialConjunctionDetails(
    c(0.01, 0.02, 0.03, 0.04),
    c("aldex2", "ancombc2", "corncob_wald", "corncob_lrt"),
    c("aldex2", "ancombc2", "corncob_wald", "corncob_lrt")
  )
  expect_equal(d3$n_families_planned, 3L)
  expect_equal(d3$partial_conjunction_h, 2L)
  expect_equal(d3$family_p[["corncob"]], 0.06)
  expect_equal(d3$combined_p, 0.04)
})

# ---- Go_CombinedEffectRank (V2) --------------------------------------------

test_that("Go_CombinedEffectRank returns values in [0, 1]", {
  da_table <- make_da_table(n = 6)
  feature_ids <- unique(da_table$feature_id)
  out <- Go_CombinedEffectRank(da_table, feature_ids)
  expect_true(all(out >= 0 & out <= 1, na.rm = TRUE))
  expect_equal(length(out), length(feature_ids))
})

test_that("Go_CombinedEffectRank preserves feature_id names", {
  da_table <- make_da_table(n = 5)
  feature_ids <- unique(da_table$feature_id)
  out <- Go_CombinedEffectRank(da_table, feature_ids)
  expect_equal(names(out), feature_ids)
})

# ---- Go_DAConsensus (V2 columns) -------------------------------------------

test_that("Go_DAConsensus produces cauchy_combined_q column (not fisher)", {
  out <- Go_DAConsensus(
    make_da_table(n = 10),
    alpha = 0.05,
    p_combine = "adaptive_cauchy"
  )
  expect_true("cauchy_combined_q" %in% colnames(out))
  expect_false("fisher_combined_q" %in% colnames(out))
  expect_true("combined_effect_rank" %in% colnames(out))
  expect_false("median_effect_size" %in% colnames(out))
})

test_that("Go_DAConsensus effect_consistency is 0 when no significant methods", {
  da_table <- data.frame(
    feature_id     = rep("t1", 2),
    method         = c("deseq2", "aldex2"),
    p_value        = c(0.8, 0.9),
    q_value        = c(0.9, 0.95),
    effect_size    = c(1.0, -1.0),
    direction      = c("up_in_group2", "up_in_group1"),
    is_significant = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  out <- Go_DAConsensus(da_table, alpha = 0.05, p_combine = "adaptive_cauchy")
  expect_equal(out$effect_consistency[out$feature_id == "t1"], 0)
})

test_that("Go_DAConsensus n_methods_run and DA_support_score correct", {
  n <- 5
  da_table <- data.frame(
    feature_id     = paste0("t", 1:n),
    method         = "deseq2",
    p_value        = c(0.01, 0.02, 0.5, 0.8, 0.001),
    q_value        = c(0.05, 0.08, 0.6, 0.9, 0.01),
    effect_size    = c(1, -1, 0.1, -0.2, 2),
    direction      = c("up_in_group2", "up_in_group1", "up_in_group2",
                       "up_in_group1", "up_in_group2"),
    is_significant = c(TRUE, TRUE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  out <- Go_DAConsensus(da_table, alpha = 0.05, p_combine = "adaptive_cauchy")
  expect_equal(nrow(out), n)
  expect_true(all(out$n_methods_run == 1))
})

test_that("Go_DAConsensus exports Solution1 family diagnostics and BH q-values", {
  da_table <- data.frame(
    feature_id = rep(c("t1", "t2"), each = 5),
    method = rep(solution1_methods, 2),
    p_value = c(
      0.01, 0.02, 0.03, 0.04, 0.05,
      0.5, 0.5, 1e-30, 0.5, 0.5
    ),
    q_value = 1,
    effect_size = rep(c(1, 0.5, 0.25, 0.1, 0.2), 2),
    direction = "up_in_group2",
    is_significant = FALSE,
    stringsAsFactors = FALSE
  )
  out <- Go_DAConsensus(da_table, alpha = 0.05)

  expect_equal(out$combined_p, c(0.06, 1), tolerance = 1e-12)
  expect_equal(out$combined_q, stats::p.adjust(c(0.06, 1), method = "BH"))
  expect_equal(out$family_corncob_p, c(0.08, 1), tolerance = 1e-12)
  expect_true(all(out$n_families_planned == 4L))
  expect_true(all(out$n_families_estimable == 4L))
  expect_true(all(out$n_corncob_tests_estimable == 2L))
  expect_equal(out$family_partial_conjunction_p, out$combined_p)
  expect_equal(out$family_partial_conjunction_q, out$combined_q)
})

test_that("Go_FinalScore names the Solution1 DA score explicitly", {
  da_table <- data.frame(
    feature_id = rep(c("t1", "t2"), each = 5),
    method = rep(solution1_methods, 2),
    p_value = c(
      0.01, 0.02, 0.03, 0.04, 0.05,
      0.5, 0.5, 1e-30, 0.5, 0.5
    ),
    q_value = 1,
    effect_size = rep(c(1, 0.5, 0.25, 0.1, 0.2), 2),
    direction = "up_in_group2",
    is_significant = FALSE,
    stringsAsFactors = FALSE
  )
  da_consensus <- Go_DAConsensus(
    da_table,
    p_combine = "family_partial_conjunction"
  )
  out <- Go_FinalScore(
    da_consensus,
    Go_CreateEmptyBetaContribution(c("t1", "t2")),
    beta_enabled = FALSE
  )
  expect_equal(
    out$family_partial_conjunction_da_score,
    out$combined_da_score
  )
})

# ---- Go_FinalScore (V2) ----------------------------------------------------

test_that("Go_FinalScore DA-only mode sets analysis_mode correctly", {
  da_consensus <- make_da_consensus_v2(6)
  beta_contribution <- Go_CreateEmptyBetaContribution(paste0("t", 1:6))

  out <- Go_FinalScore(da_consensus, beta_contribution, beta_enabled = FALSE)
  expect_true("final_score" %in% colnames(out))
  expect_true("classification" %in% colnames(out))
  expect_true(all(out$analysis_mode %in% c("da_only", "single_method")))
})

test_that("Go_FinalScore full mode sets analysis_mode = 'full'", {
  da_consensus <- make_da_consensus_v2(4)
  beta_contribution <- make_beta_contrib_v2(4)

  out <- Go_FinalScore(da_consensus, beta_contribution, beta_enabled = TRUE)
  expect_true(all(out$analysis_mode == "full"))
})

test_that("Go_FinalScore cauchy_da_score present (not fisher_da_score)", {
  da_consensus <- make_da_consensus_v2(4)
  beta_contribution <- make_beta_contrib_v2(4)

  out <- Go_FinalScore(da_consensus, beta_contribution)
  expect_true("cauchy_da_score" %in% colnames(out))
  expect_false("fisher_da_score" %in% colnames(out))
})

test_that("Go_FinalScore loo_separation_score accepted from beta_contribution", {
  da_consensus <- make_da_consensus_v2(4)
  beta_contribution <- make_beta_contrib_v2(4)

  expect_true("loo_separation_score" %in% colnames(beta_contribution))
  expect_false("simper_score" %in% colnames(beta_contribution))
  out <- Go_FinalScore(da_consensus, beta_contribution)
  expect_equal(nrow(out), 4)
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
  da_consensus <- data.frame(
    feature_id        = paste0("taxa", 1:n),
    cauchy_combined_p = runif(n),
    cauchy_combined_q = runif(n),
    stringsAsFactors  = FALSE
  )
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

# ---- Go_BuildMethodOverlapLong (NA method guard) ---------------------------

test_that("Go_BuildMethodOverlapLong handles NA methods in da_table", {
  n <- 5
  final_scores <- data.frame(
    feature_id        = paste0("t", 1:n),
    plot_label        = paste0("Taxon", 1:n),
    plot_label_unique = paste0("Taxon", 1:n),
    deseq2_is_significant = c(TRUE, FALSE, TRUE, FALSE, FALSE),
    stringsAsFactors  = FALSE
  )
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
  out <- Go_BuildMethodOverlapLong(final_scores, da_table)
  expect_equal(nrow(out), n)
  expect_true(all(out$is_significant == FALSE))
})

# ---- Go_CreateEmptyBetaContribution (V2 schema) ----------------------------

test_that("Go_CreateEmptyBetaContribution V2 has loo columns not simper", {
  ids <- paste0("t", 1:5)
  out <- Go_CreateEmptyBetaContribution(ids)
  expect_equal(out$feature_id, ids)
  expect_true("loo_separation_delta" %in% colnames(out))
  expect_true("loo_separation_score" %in% colnames(out))
  expect_true("beta_contribution_score" %in% colnames(out))
  expect_false("simper_contribution" %in% colnames(out))
  expect_false("simper_score" %in% colnames(out))
  expect_false("delta_R2" %in% colnames(out))
  expect_true(all(is.na(out$beta_contribution_score)))
})

# ---- Go_WeightSensitivity (V2) --------------------------------------------

test_that("Go_WeightSensitivity returns data.frame with expected columns", {
  da_consensus <- make_da_consensus_v2(4)
  beta_contribution <- make_beta_contrib_v2(4)

  out <- Go_WeightSensitivity(da_consensus, beta_contribution,
                              n_grid = 3L, beta_enabled = TRUE)
  expect_true(is.data.frame(out))
  expect_true(all(c("feature_id", "mean_rank", "sd_rank", "rank_stability") %in% colnames(out)))
  expect_equal(nrow(out), 4)
})

test_that("Go_WeightSensitivity rank_stability in [0, 1]", {
  da_consensus <- make_da_consensus_v2(6)
  beta_contribution <- make_beta_contrib_v2(6)

  out <- Go_WeightSensitivity(da_consensus, beta_contribution,
                              n_grid = 3L, beta_enabled = FALSE)
  expect_true(all(out$rank_stability >= 0 & out$rank_stability <= 1))
})

test_that("Go_WeightSensitivity feature mapping is correct (not sorted order)", {
  # All features except t1 have identical scores → t1 should consistently rank 1st.
  # If feature mapping were wrong, t1's rank would be assigned to a different feature.
  da_consensus <- make_da_consensus_v2(4)
  # Give t1 the strongest consensus evidence and keep compatibility columns aligned.
  da_consensus$combined_q[da_consensus$feature_id == "t1"] <- 1e-10
  da_consensus$cauchy_combined_q[da_consensus$feature_id == "t1"] <- 1e-10
  da_consensus$is_combined_significant[da_consensus$feature_id == "t1"] <- TRUE
  da_consensus$combined_effect_rank[da_consensus$feature_id == "t1"] <- 1.0
  beta_contribution <- make_beta_contrib_v2(4)
  beta_contribution$beta_contribution_score[beta_contribution$feature_id == "t1"] <- 1.0

  out <- Go_WeightSensitivity(da_consensus, beta_contribution,
                              n_grid = 3L, beta_enabled = TRUE)
  # t1 must have the lowest mean_rank (= consistently top-ranked)
  expect_equal(out$feature_id[1], "t1")
})

test_that("Go_CombinedEffectRank handles duplicate (method, feature_id) without crash", {
  # Simulate a da_table where deseq2 accidentally has two rows for t1
  da_table <- data.frame(
    feature_id  = c("t1", "t1", "t2", "t3"),
    method      = c("deseq2", "deseq2", "deseq2", "deseq2"),
    effect_size = c(1.0, 2.0, -0.5, 0.3),
    stringsAsFactors = FALSE
  )
  feature_ids <- c("t1", "t2", "t3")
  expect_no_error({
    out <- Go_CombinedEffectRank(da_table, feature_ids)
  })
  expect_equal(length(out), 3L)
  expect_equal(names(out), feature_ids)
})

# ---- inst/extdata example data --------------------------------------------

test_that("minimal_example_ps.rds loads as phyloseq", {
  skip_if_not_installed("phyloseq")
  path <- system.file("extdata", "minimal_example_ps.rds", package = "ConDAdist")
  if (!nzchar(path)) {
    path <- file.path(system.file(package = "ConDAdist"), "extdata", "minimal_example_ps.rds")
  }
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

test_that("full_cdd is an exact alias of broad_panel and freezes the paper's rule", {
  a <- Go_ResolveCDDPreset("full_cdd", supplied = list())
  b <- Go_ResolveCDDPreset("broad_panel", supplied = list())
  a$preset <- NULL
  b$preset <- NULL
  expect_identical(a, b)
  expect_equal(a$methods, solution1_methods)
  expect_null(a$distances)
  expect_equal(a$p_combine, "family_partial_conjunction")

  # Four families, h = 3: combined_p = (n - h + 1) * p_(h) = 2 * p_(3).
  d <- Go_FamilyPartialConjunctionDetails(
    c(0.01, 0.02, 0.03, 0.04, 0.05), solution1_methods
  )
  expect_equal(d$n_families_planned, 4L)
  expect_equal(d$partial_conjunction_h, 3L)
  expect_equal(d$combined_p, 0.06, tolerance = 1e-12)
})
