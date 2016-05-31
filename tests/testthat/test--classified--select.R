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
                  taxon_data = data.frame(name = letters[1:5], other_col = 1:5, stringsAsFactors = FALSE),
                  item_data = data.frame(item_attr = LETTERS[1:10], other_item_col = 1:5,  stringsAsFactors = FALSE))
#|
#| ####  Selecting with unquoted column names
test_that("Selecting with unquoted column names works", {
  result <- select_taxa(obj, name)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_s3_class(result, "classified")
  expect_equivalent(result$taxon_data$name, c("a", "b"))
  expect_equivalent(item_counts(result), c(10, 7))
})

