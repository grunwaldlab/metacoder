#| ## Testing `extract_taxonomy` internal functions
#|
library(metacoder)
context("'extract_taxonomy' internal functions")
#|
#| ### Validating that input is matched be regex
#|
#| #### Create test data
test_regex <- "match"
test_input <- c("matching", "not", "also a match")
#| 
#| #### Test
test_that("Input validation works", {
  expect_equal(validate_regex_match(test_input[1], test_regex),
               test_input[1])
  expect_error(validate_regex_match(test_input, test_regex),
               "could not be matched by the regex ")
  expect_warning(validate_regex_match(test_input, test_regex, mismatch_action = "warn"),
                 "They will be excluded.")
  expect_equal(validate_regex_match(test_input, test_regex, mismatch_action = "allow"),
               test_input[c(1, 3)])
})
#|
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
               setNames(test_key, test_key))
  expect_error(validate_regex_key_pair(regex = "", key = test_key, key_options = test_options))
})
test_that("default naming of keys works", {
  expect_true(all(names(validate_regex_key_pair(regex = test_regex, key = test_key, key_options = test_options)) == test_key))
  expect_true(all(names(validate_regex_key_pair(regex = test_regex, key = c(x = "taxon_name", "taxon_id"), key_options = test_options)) == c("x", "taxon_id")))
})

#|


