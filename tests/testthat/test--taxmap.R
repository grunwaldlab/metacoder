#| ## Testing `taxmap` object
#|
library(metacoder)
context("Creating `taxmap` object")
#|
#| ### Creating `taxmap` objects
#|
#| #### Minimal usage
test_that("Minimal usage works", {
  result_from_num <- taxmap(taxon_ids = numeric(0), supertaxon_ids =  numeric(0))
  expect_s3_class(result_from_num, "taxmap")
  expect_equal(dim(result_from_num$taxon_data), c(0, 2))
  expect_equal(dim(result_from_num$obs_data), c(0, 1))
  result_from_char <- taxmap(taxon_ids = character(0), supertaxon_ids =  character(0))
  expect_identical(result_from_num, result_from_char)
})
#|
#| #### Basic use
test_that("Basic use works", {
  result_from_num <- taxmap(taxa = c(1, 2, 3, 4),
                                supertaxon_ids = c(NA, 1, 2, 3),
                                obs_taxon_ids = c(2, 2, 1, 1, 3),
                                taxon_ids = LETTERS[1:4])
  expect_s3_class(result_from_num, "taxmap")
  expect_equal(dim(result_from_num$taxon_data), c(4, 2))
  expect_equal(dim(result_from_num$obs_data), c(5, 1))

  result_from_char <- taxmap(taxon_ids = LETTERS[1:4],
                                 supertaxon_ids = c(NA, LETTERS[1:3]),
                                 obs_taxon_ids = LETTERS[c(2, 2, 1, 1, 3)])
  expect_identical(result_from_char$taxon_data, result_from_num$taxon_data)  
  expect_identical(result_from_char$obs_data, result_from_num$obs_data)  
  
  result_from_mixed <- taxmap(taxon_ids = LETTERS[1:4],
                                  supertaxon_ids = c(NA, LETTERS[1:3]),
                                  obs_taxon_ids = c(2, 2, 1, 1, 3))
  expect_identical(result_from_mixed$taxon_data, result_from_num$taxon_data)  
  expect_identical(result_from_mixed$obs_data, result_from_num$obs_data)  
})
#|
#| #### print method
test_that("Print method works", {
  result <- taxmap(taxon_ids = LETTERS[1:4],
                   supertaxon_ids = c(NA, LETTERS[1:3]),
                   obs_taxon_ids = c(2, 2, 1, 1, 3))
  expect_output(print(result), "`taxmap` object with")
  expect_output(print(result), "A tibble")
  expect_output(print(result), "taxon_data")
  expect_output(print(result), "A, B, C, D")
})
