test_that("failed native adapters are skipped without Wilcoxon substitution", {
  ft <- matrix(c(1:10, 101:110), nrow = 1,
               dimnames = list("t1", paste0("s", 1:20)))
  md <- data.frame(group = rep(c("A", "B"), each = 10),
                   row.names = colnames(ft))
  methods <- c("ancombc2", "corncob_wald", "corncob_lrt")
  adapters <- list(Go_DA_ancombc2, Go_DA_corncob_wald, Go_DA_corncob_lrt)
  raw <- setNames(lapply(adapters, function(f) {
    suppressMessages(f(ft, md, "group", "A", "B"))
  }), methods)
  standardized <- Go_StandardizeDA(raw, "A.vs.B")$all_methods_standardized
  expect_equal(standardized$method, methods)
  expect_false(any(standardized$method == "wilcoxon"))
  expect_true(all(is.na(standardized$p_value)))
  expect_true(all(is.na(standardized$effect_size)))
  expect_true(all(grepl("skipped", standardized$notes, ignore.case = TRUE)))
  expect_false(any(standardized$is_significant))
})

test_that("a skipped single method does not create a volcano bridge", {
  skipped <- Go_CreateSkippedAdapterResult(
    feature_ids = c("t1", "t2"), method = "deseq2",
    note = "Skipped: test failure."
  )
  out_dir <- tempfile("condadist-skipped-")
  out <- suppressMessages(Go_ExportVolcanoBridge(
    output_dir = out_dir, single_output_dir = out_dir,
    da_table = skipped, final_scores = data.frame(),
    filtered_metadata = data.frame(group = c("A", "A", "B", "B")),
    group_var = "group", group_1 = "A", group_2 = "B",
    methods = "deseq2", write_consensus = FALSE
  ))
  expect_length(out$files, 0)
  expect_length(list.files(out_dir, pattern = "[.]csv$"), 0)
})
