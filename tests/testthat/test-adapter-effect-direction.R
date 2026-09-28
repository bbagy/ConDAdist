## Regression test (2026-09-28): corncob and MaAsLin2 adapters negated their
## coefficients, reporting taxa increased in group_2 as up_in_group1.
test_that("adapters report a taxon increased in group_2 as up_in_group2", {
  set.seed(1)
  n <- 20; p <- 30; grp <- rep(c("A", "B"), each = n / 2)
  counts <- matrix(stats::rnbinom(p * n, mu = 200, size = 5), p, n,
                   dimnames = list(paste0("t", seq_len(p)), paste0("s", seq_len(n))))
  counts["t1", grp == "B"] <- counts["t1", grp == "B"] * 5
  md <- data.frame(Group = grp, row.names = colnames(counts))
  check <- function(fun, pkg) {
    skip_if_not_installed(pkg)
    res <- suppressMessages(fun(counts, md, "Group", "A", "B"))
    expect_gt(res$effect_size[res$feature_id == "t1"], 0)
    expect_equal(res$direction[res$feature_id == "t1"], "up_in_group2")
  }
  check(Go_DA_deseq2, "DESeq2")
  check(Go_DA_corncob_wald, "corncob")
  check(Go_DA_corncob_lrt, "corncob")
  check(Go_DA_maaslin, "Maaslin2")
})
