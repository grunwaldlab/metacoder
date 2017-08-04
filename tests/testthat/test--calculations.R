library(metacoder)
context("Calculations")

test_that("Observation proportions", {
  # Get test data
  raw_data <- read.table("example_data/small_abund_mat.tsv", sep = "\t")
  data <- taxa::parse_tax_data(raw_data, class_cols = "tax", class_sep = ";")

  # Default usage
  result <- calc_obs_props(data, dataset = "tax_data", cols = c("sam_1", "sam_2"))
  expect_equal(sum(result$sam_1), 1)
  expect_equal(sum(result$sam_2), 1)
  expect_true("otu_id" %in% colnames(result))
  
  # Discarding other columns
  result <- calc_obs_props(data, dataset = "tax_data", cols = c("sam_1", "sam_2"),
                           keep_other_cols = FALSE)
  expect_equal(sum(result$sam_1), 1)
  expect_equal(sum(result$sam_2), 1)
  expect_true(! "otu_id" %in% colnames(result))
  expect_true("taxon_id" %in% colnames(result))
})
