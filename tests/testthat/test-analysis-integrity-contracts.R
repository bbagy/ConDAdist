test_that("DA covariates and random effects retain separate roles", {
  ft <- matrix(seq_len(24), nrow = 3,
               dimnames = list(paste0("t", 1:3), paste0("s", 1:8)))
  md <- data.frame(
    group = rep(c("A", "B"), each = 4),
    age = seq_len(8),
    subject = rep(1:4, 2),
    row.names = colnames(ft)
  )
  out <- Go_PrepareDAInputs(
    ft, md, "group", "A", "B",
    covariates = "age", random_effects = "subject"
  )
  expect_equal(out$covariates, "age")
  expect_equal(out$random_effects, "subject")
  expect_false("subject" %in% out$covariates)
})

test_that("unknown or duplicated DA model roles fail loudly", {
  ft <- matrix(seq_len(12), nrow = 2,
               dimnames = list(c("t1", "t2"), paste0("s", 1:6)))
  md <- data.frame(group = rep(c("A", "B"), each = 3), age = 1:6,
                   row.names = colnames(ft))
  expect_error(
    Go_PrepareDAInputs(ft, md, "group", "A", "B", covariates = "BMI"),
    "Unknown DA-model metadata column"
  )
  expect_error(
    Go_PrepareDAInputs(
      ft, md, "group", "A", "B",
      covariates = "age", random_effects = "age"
    ),
    "both fixed covariates and random effects"
  )
})

test_that("complete-case cohort is fixed before feature filtering", {
  ft <- matrix(seq_len(24), nrow = 3,
               dimnames = list(paste0("t", 1:3), paste0("s", 1:8)))
  md <- data.frame(
    group = rep(c("A", "B"), each = 4),
    age = c(1:7, NA),
    row.names = colnames(ft)
  )
  out <- Go_PrepareAnalysisCohort(
    ft, md, "group", "A", "B", covariates = "age"
  )
  expect_equal(ncol(out$feature_table), 7)
  expect_equal(rownames(out$metadata), colnames(out$feature_table))
  expect_false(anyNA(out$metadata$age))
})

test_that("methods without mixed-model support are explicitly skipped", {
  ft <- matrix(seq_len(24), nrow = 3,
               dimnames = list(paste0("t", 1:3), paste0("s", 1:8)))
  md <- data.frame(
    group = rep(c("A", "B"), each = 4),
    subject = rep(1:4, 2),
    row.names = colnames(ft)
  )
  out <- suppressMessages(Go_RunDAmethods(
    ft, md, "group", "A", "B",
    random_effects = "subject", methods = "deseq2"
  ))$deseq2
  expect_true(all(is.na(out$p_value)))
  expect_match(unique(out$notes), "does not support random effects")
})

test_that("skipped rows are not labeled as detected native results", {
  skipped <- Go_CreateSkippedAdapterResult("t1", "ancombc2", "Skipped: test")
  out <- Go_BuildMethodAnnotation(skipped)
  expect_false(out$ancombc2_detected)
})

test_that("deterministically skipped methods do not trigger a futile retry", {
  skipped <- Go_CreateSkippedAdapterResult(
    c("t1", "t2"), "deseq2",
    "Skipped: deseq2 does not support random effects."
  )
  filtered <- list(feature_table = matrix(1, nrow = 2, ncol = 2))
  expect_false(Go_ShouldRetryAnalysis(
    filtered, skipped, methods = "deseq2", attempt_id = 1L, max_attempts = 2L
  ))
})

test_that("zero p-values remain strong evidence in legacy combiners", {
  expect_equal(Go_CombinePValuesFisher(c(0.01, 0.2)),
               stats::pchisq(-2 * sum(log(c(0.01, 0.2))), 4, lower.tail = FALSE))
  expect_lte(Go_CombinePValuesFisher(c(0, 0.5)), 1e-250)
  expect_true(is.finite(Go_CombinePValuesCauchy(c(0, 0.5))))
  expect_lt(Go_CombinePValuesAdaptiveCauchy(c(0, 0.8)), 1e-10)
})
