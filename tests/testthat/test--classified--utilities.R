#| ## Testing utility methods for `classified` objects
#|
library(metacoder)
context("Utility methods for `classified` objects")
#|
#| ### Getting supertaxa
#|
#| ####  Code shared by tests
obj <- classified(taxon_ids = LETTERS[1:5], parent_ids = c(NA, 1, 2, 2, 1), 
                  item_taxon_ids = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],  stringsAsFactors = FALSE),
                  item_data = data.frame(item_attr = letters[1:10],  stringsAsFactors = FALSE))
original_plot <- plot(obj, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
#|
#| ####  Immediate supertaxa
test_that("Getting immediate supertaxa works", {
  result <- supertaxa(obj, recursive = FALSE, simplify = FALSE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "list")
  expect_equal(names(result), obj$taxon_data$taxon_ids)
  expect_equal(result$A, character(0))
  expect_equal(result$B, "A")
  result <- supertaxa(obj, recursive = FALSE, simplify = FALSE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(sapply(result, function(x) x[1]) ==  names(result)))
  result <- supertaxa(obj, recursive = FALSE, simplify = TRUE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "character")
  result <- supertaxa(obj, recursive = FALSE, simplify = TRUE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(result %in% obj$taxon_data$taxon_ids))
})
#|
#| ####  All supertaxa
test_that("Getting all supertaxa works", {
  result <- supertaxa(obj, recursive = TRUE, simplify = FALSE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "list")
  expect_equal(names(result), obj$taxon_data$taxon_ids)
  expect_equal(result$A, character(0))
  expect_equal(result$B, "A")
  result <- supertaxa(obj, recursive = TRUE, simplify = FALSE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(sapply(result, function(x) x[1]) ==  names(result)))
  result <- supertaxa(obj, recursive = TRUE, simplify = TRUE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "character")
  result <- supertaxa(obj, recursive = TRUE, simplify = TRUE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(result %in% obj$taxon_data$taxon_ids))
})
#|
#| ### Getting subtaxa
#|
#| ####  Code shared by tests
test_that("Getting immediate subtaxa works", {
  result <- subtaxa(obj, recursive = FALSE, simplify = FALSE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "list")
  expect_equal(names(result), obj$taxon_data$taxon_ids)
  expect_equal(result$E, character(0))
  expect_equal(result$B, c("C", "D"))
  result <- subtaxa(obj, recursive = FALSE, simplify = FALSE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(sapply(result, function(x) x[1]) ==  names(result)))
  result <- subtaxa(obj, recursive = FALSE, simplify = TRUE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "character")
  result <- subtaxa(obj, recursive = FALSE, simplify = TRUE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(result %in% obj$taxon_data$taxon_ids))
})
test_that("Getting all subtaxa works", {
  result <- subtaxa(obj, recursive = TRUE, simplify = FALSE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "list")
  expect_equal(names(result), obj$taxon_data$taxon_ids)
  expect_equal(result$A, LETTERS[2:5])
  expect_equal(result$E, character(0))
  result <- subtaxa(obj, recursive = TRUE, simplify = FALSE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(sapply(result, function(x) x[1]) ==  names(result)))
  result <- subtaxa(obj, recursive = TRUE, simplify = TRUE, include_input = FALSE, index = FALSE, na = FALSE)
  expect_equal(class(result), "character")
  result <- subtaxa(obj, recursive = TRUE, simplify = TRUE, include_input = TRUE, index = FALSE, na = FALSE)
  expect_true(all(result %in% obj$taxon_data$taxon_ids))
})
