#| ## Testing `extract_taxonomy` internal functions
#|
library(metacoder)
context("'extract_taxonomy' internal functions")
#|
#| ### Validating regex-key pairs
#|
#| #### Create test data
test_regex = "z fsasdfasfsfad f other(.*)random(.*)stuff sdf f \n sdff"
test_key = c("taxon_name", "taxon_id")
test_options <- c("taxon_id", "taxon_name", "taxon_info", "class", "item_id", "item_info")
#| 
#| #### Test capture group counting
test_that("capture group counting works", {
  expect_equal(count_capture_groups(test_regex), 2)
  expect_equal(count_capture_groups(""), 0)
})
#| 
#| #### Test validation
test_that("validating regex-key pairs works", {
  expect_equal(validate_regex_key_pair(regex = test_regex, key = test_key, key_options = test_options),
               test_key)
  expect_error(validate_regex_key_pair(regex = "", key = test_key, key_options = test_options))
})



