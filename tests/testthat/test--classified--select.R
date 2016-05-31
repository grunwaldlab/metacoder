#| ## Testing column subset for `classified` objects
#|
library(metacoder)
context("Subsetting columns of `classified` objects")
#|
#| ### Select `taxon_data` columns
#|
#| ####  Code shared by tests
obj <- classified(taxa = c(1, 2, 3, 4, 5), parents = c(NA, 1, 2, 2, 1), 
                  item_taxa = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],
                                          other_col = 1:5,
                                          yet_another = 6:10, 
                                          stringsAsFactors = FALSE),
                  item_data = data.frame(item_attr = LETTERS[1:10], other_item_col = 1:5,  stringsAsFactors = FALSE))
#|
#| ####  Selecting with unquoted column names
test_that("Selecting with unquoted column names works", {
  result <- select_taxa(obj, name)
  expect_s3_class(result, "classified")
  expect_false("other_col" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Selecting with quoted column names
test_that("Selecting with quoted column names works", {
  result <- select_taxa(obj, "name")
  expect_s3_class(result, "classified")
  expect_false("other_col" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Selecting with logical vector
test_that("Selecting with logical vector works", {
  result <- select_taxa(obj, "name" == colnames(obj$taxon_data))
  expect_s3_class(result, "classified")
  expect_false("other_col" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Selecting with index
test_that("Selecting with index works", {
  result <- select_taxa(obj, which("name" == colnames(obj$taxon_data)))
  expect_s3_class(result, "classified")
  expect_false("other_col" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Selecting with multiple types
test_that("Selecting with multiple types works", {
  result <- select_taxa(obj, which("name" == colnames(obj$taxon_data)),
                        other_col, "taxon_ids", 2)
  expect_s3_class(result, "classified")
  expect_false("yet_another" %in% colnames(result$taxon_data))
  expect_true(all(c("name", "taxon_ids", "parent_ids", "other_col") %in% colnames(result$taxon_data)))
})