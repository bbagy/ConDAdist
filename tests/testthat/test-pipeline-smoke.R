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
