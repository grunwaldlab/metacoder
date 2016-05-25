#| ## Testing `classified` object
#|
library(metacoder)
context("Creating `classified` object")
#|
#| ### Creating `classified` objects
#|
#| #### Minimal usage
test_that("Minimal usage works", {
  result_from_num <- classified(taxa = numeric(0), parents =  numeric(0))
  expect_s3_class(result_from_num, "classified")
  expect_equal(result_from_num$taxa, character(0))
  expect_equal(dim(result_from_num$taxon_data), c(0, 2))
  expect_equal(dim(result_from_num$item_data), c(0, 1))
  result_from_char <- classified(taxa = character(0), parents =  character(0))
  expect_identical(result_from_num, result_from_char)
})
#|
#| #### Basic use
test_that("Basic use works", {
  result_from_num <- classified(taxa = c(1, 2, 3, 4),
                                parents = c(NA, 1, 2, 3),
                                item_taxa = c(2, 2, 1, 1, 3))
  expect_s3_class(result_from_num, "classified")
  expect_equal(result_from_num$taxa, as.character(c(1, 2, 3, 4)))
  expect_equal(dim(result_from_num$taxon_data), c(4, 2))
  expect_equal(dim(result_from_num$item_data), c(5, 1))
  expect_equal(class(result_from_num$taxa), "character")
  
  result_from_offset <- classified(taxa = c(1, 2, 3, 4) + 7,
                                   parents = c(NA, 1, 2, 3),
                                   item_taxa = c(2, 2, 1, 1, 3))
  expect_identical(result_from_offset$taxon_data, result_from_num$taxon_data)  
  expect_identical(result_from_offset$item_data, result_from_num$item_data)  
  
  result_from_char <- classified(taxa = letters[1:4],
                                 parents = c(NA, letters[1:3]),
                                 item_taxa = letters[c(2, 2, 1, 1, 3)])
  expect_identical(result_from_char$taxon_data, result_from_num$taxon_data)  
  expect_identical(result_from_char$item_data, result_from_num$item_data)  
  
  result_from_mixed <- classified(taxa = letters[1:4],
                                 parents = c(NA, letters[1:3]),
                                 item_taxa = c(2, 2, 1, 1, 3))
  expect_identical(result_from_mixed$taxon_data, result_from_num$taxon_data)  
  expect_identical(result_from_mixed$item_data, result_from_num$item_data)  
})


