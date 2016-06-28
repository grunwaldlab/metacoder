#| ## Testing column subset for `taxmap` objects
#|
library(metacoder)
context("Subsetting columns of `taxmap` objects")
#|
#| ### Select `taxon_data` columns
#|
#| ####  Code shared by tests
obj <- taxmap(taxon_ids = LETTERS[1:5], parent_ids = c(NA, 1, 2, 2, 1), 
                  obs_taxon_ids = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],
                                          other_col = 1:5,
                                          yet_another = 6:10, 
                                          stringsAsFactors = FALSE),
                  obs_data = data.frame(obs_attr = LETTERS[1:10],
                                         other_obs_col = 1:10,
                                         another_obs_col = 11:20, 
                                         stringsAsFactors = FALSE))
#|
#| ####  Selecting taxon_data with unquoted column names
test_that("Selecting taxon_data with unquoted column names works", {
  result <- select_taxa(obj, name)
  expect_s3_class(result, "taxmap")
  expect_false("other_col" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Selecting taxon_data using dplyr functions
test_that("Selecting taxon_data with dplyr functions works", {
  result <- select_taxa(obj, matches("name"))
  expect_s3_class(result, "taxmap")
  expect_false("other_col" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Selecting taxon_data with index
test_that("Selecting taxon_data with index works", {
  result <- select_taxa(obj, which("name" == colnames(obj$taxon_data)))
  expect_s3_class(result, "taxmap")
  expect_false("other_col" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Selecting taxon_data with multiple types
test_that("Selecting taxon_data with multiple types works", {
  result <- select_taxa(obj, which("name" == colnames(obj$taxon_data)),
                        other_col, 2)
  expect_s3_class(result, "taxmap")
  expect_false("yet_another" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids", "other_col") %in% colnames(result$taxon_data)))
})
#|
#|
#| ### Select `obs_data` columns
#|
#|
#| ####  Selecting obs_data with unquoted column names
test_that("Selecting obs_data with unquoted column names works", {
  result <- select_obs(obj, obs_attr)
  expect_s3_class(result, "taxmap")
  expect_false("other_obs_col" %in% colnames(result$obs_data))
  expect_true(all(c("obs_attr", "obs_taxon_ids") %in% colnames(result$obs_data)))
})
#|
#| ####  Selecting obs_data with  dplyr functions
test_that("Selecting obs_data with  dplyr functionsworks", {
  result <- select_obs(obj, matches("obs_attr"))
  expect_s3_class(result, "taxmap")
  expect_false("other_obs_col" %in% colnames(result$obs_data))
  expect_true(all(c("obs_attr", "obs_taxon_ids") %in% colnames(result$obs_data)))
})
#|
#| ####  Selecting obs_data with index
test_that("Selecting obs_data with index works", {
  result <- select_obs(obj, which("obs_attr" == colnames(obj$obs_data)))
  expect_s3_class(result, "taxmap")
  expect_false("other_obs_col" %in% colnames(result$obs_data))
  expect_true(all(c("obs_attr", "obs_taxon_ids") %in% colnames(result$obs_data)))
})
#|
#| ####  Selecting obs_data with multiple types
test_that("Selecting obs_data with multiple types works", {
  result <- select_obs(obj, which("obs_attr" == colnames(obj$obs_data)),
                        other_obs_col)
  expect_s3_class(result, "taxmap")
  expect_false("another_obs_col" %in% colnames(result$obs_data))
  expect_true(all(c("obs_attr", "obs_taxon_ids", "other_obs_col") %in% colnames(result$obs_data)))
})
