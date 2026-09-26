# Average raw leave-one-taxon-out deltas across distances before applying
# percentile-clip normalization, which is nonlinear.

make_ft <- function(n_taxa = 10, n_samples = 10, seed = 42) {
  set.seed(seed)
  mat <- matrix(
    as.integer(stats::rpois(n_taxa * n_samples, lambda = 80)),
    nrow = n_taxa,
    dimnames = list(paste0("taxa", seq_len(n_taxa)), paste0("s", seq_len(n_samples)))
  )
  storage.mode(mat) <- "numeric"
  mat
}

make_md <- function(n_samples = 10) {
  data.frame(
    Group = rep(c("Control", "Treatment"), each = n_samples / 2),
    row.names = paste0("s", seq_len(n_samples)),
    stringsAsFactors = FALSE
  )
}

test_that("joint multi-metric Go_BetaContribution matches combine-raw-then-normalize-once", {
  ft <- make_ft()
  md <- make_md()
  group_factor <- factor(md$Group, levels = c("Control", "Treatment"))
  feature_ids <- rownames(ft)

  full_bray <- Go_GroupSeparationScore(as.matrix(Go_Dist_bray(ft)), group_factor)
  full_jsd <- Go_GroupSeparationScore(as.matrix(Go_Dist_jsd(ft)), group_factor)
  full_scores <- list(bray = full_bray, jsd = full_jsd)
  dist_matrices <- list(bray = Go_Dist_bray(ft), jsd = Go_Dist_jsd(ft))

  joint_delta <- vapply(feature_ids, function(fid) {
    Go_LeaveOneTaxonOutScore(
      feature_table = ft, target_feature = fid, group_factor = group_factor,
      distances = dist_matrices, full_scores = full_scores, phy_tree = NULL
    )
  }, numeric(1))
  expected_joint_normalized <- Go_NormalizeVector(joint_delta)

  bd_joint <- Go_BetaDistance(
    feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", distances = c("bray", "jsd"),
    phy_tree = NULL, n_permutations = 1L
  )
  bc_joint <- Go_BetaContribution(
    feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", beta_distances = bd_joint,
    phy_tree = NULL, n_beta_permutations = 0L, method = "loo_only"
  )
  actual_joint <- bc_joint$beta_contribution_score[match(feature_ids, bc_joint$feature_id)]

  expect_equal(unname(actual_joint), unname(expected_joint_normalized), tolerance = 1e-10)
})

test_that("joint contract is NOT equivalent to per-metric-normalize-then-average (regression guard)", {
  ft <- make_ft()
  md <- make_md()

  bd_bray <- Go_BetaDistance(feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", distances = "bray", phy_tree = NULL, n_permutations = 1L)
  bc_bray <- Go_BetaContribution(feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", beta_distances = bd_bray, phy_tree = NULL,
    n_beta_permutations = 0L, method = "loo_only")

  bd_jsd <- Go_BetaDistance(feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", distances = "jsd", phy_tree = NULL, n_permutations = 1L)
  bc_jsd <- Go_BetaContribution(feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", beta_distances = bd_jsd, phy_tree = NULL,
    n_beta_permutations = 0L, method = "loo_only")

  feature_ids <- rownames(ft)
  wrong_average <- 0.5 * (
    bc_bray$beta_contribution_score[match(feature_ids, bc_bray$feature_id)] +
    bc_jsd$beta_contribution_score[match(feature_ids, bc_jsd$feature_id)]
  )

  bd_joint <- Go_BetaDistance(feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", distances = c("bray", "jsd"), phy_tree = NULL, n_permutations = 1L)
  bc_joint <- Go_BetaContribution(feature_table = ft, metadata = md, group_var = "Group",
    group_1 = "Control", group_2 = "Treatment", beta_distances = bd_joint, phy_tree = NULL,
    n_beta_permutations = 0L, method = "loo_only")
  actual_joint <- bc_joint$beta_contribution_score[match(feature_ids, bc_joint$feature_id)]

  max_diff <- max(abs(actual_joint - wrong_average))
  expect_gt(max_diff, 1e-6)
})
