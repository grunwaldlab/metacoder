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
  vigilance = "error"
  expect_error(validate_regex_match(test_input, test_regex),
               "could not be matched by the regex ")
  vigilance = "warning"
  expect_warning(validate_regex_match(test_input, test_regex),
                 "could not be matched by the regex ")
  vigilance = "none"
  expect_equal(validate_regex_match(test_input, test_regex),
               test_input[c(1, 3)])
})
#|
#|
#| ### Validating regex-key pairs
#|
#| #### Create test data
test_regex = "z fsasdfasfsfad f other(.*)random(.*)stuff sdf f \n sdff"
test_key = c("name", "taxon_id")
#| 
#| #### Test capture group counting
test_that("capture group counting works", {
  expect_equal(count_capture_groups(test_regex), 2)
  expect_equal(count_capture_groups(""), 0)
})
