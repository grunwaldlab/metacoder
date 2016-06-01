#| ## Testing `classified` object
#|
library(metacoder)
context("Creating `classified` object")
#|
#| ### Creating `classified` objects
#|
#| #### Minimal usage
test_that("Minimal usage works", {
  result_from_num <- classified(taxon_ids = numeric(0), parent_ids =  numeric(0))
  expect_s3_class(result_from_num, "classified")
  expect_equal(dim(result_from_num$taxon_data), c(0, 2))
  expect_equal(dim(result_from_num$item_data), c(0, 1))
  result_from_char <- classified(taxon_ids = character(0), parent_ids =  character(0))
  expect_identical(result_from_num, result_from_char)
})
#|
#| #### Basic use
test_that("Basic use works", {
  result_from_num <- classified(taxa = c(1, 2, 3, 4),
                                parent_ids = c(NA, 1, 2, 3),
                                item_taxon_ids = c(2, 2, 1, 1, 3),
                                taxon_ids = LETTERS[1:4])
  expect_s3_class(result_from_num, "classified")
  expect_equal(dim(result_from_num$taxon_data), c(4, 2))
  expect_equal(dim(result_from_num$item_data), c(5, 1))

  result_from_char <- classified(taxon_ids = LETTERS[1:4],
                                 parent_ids = c(NA, LETTERS[1:3]),
                                 item_taxon_ids = LETTERS[c(2, 2, 1, 1, 3)])
  expect_identical(result_from_char$taxon_data, result_from_num$taxon_data)  
  expect_identical(result_from_char$item_data, result_from_num$item_data)  
  
  result_from_mixed <- classified(taxon_ids = LETTERS[1:4],
                                  parent_ids = c(NA, LETTERS[1:3]),
                                  item_taxon_ids = c(2, 2, 1, 1, 3))
  expect_identical(result_from_mixed$taxon_data, result_from_num$taxon_data)  
  expect_identical(result_from_mixed$item_data, result_from_num$item_data)  
})
#|
#| #### print method
test_that("Print method works", {
  result <- classified(taxon_ids = LETTERS[1:4],
                       parent_ids = c(NA, LETTERS[1:3]),
                       item_taxon_ids = c(2, 2, 1, 1, 3))
  expect_output(print(result), "`classified` object with")
  expect_output(print(result), "Source: local data frame")
  expect_output(print(result), "taxon_data")
  expect_output(print(result), "A, B, C, D")
})
