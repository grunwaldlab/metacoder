#| ## Testing sorting methods for `classified` objects
#|
library(metacoder)
context("Sorting `classified` objects")
#|
#| ### Sorting taxa
#|
#| ####  Code shared by tests
obj <- classified(taxa = c(1, 2, 3, 4, 5), parents = c(NA, 1, 2, 2, 1), 
                  item_taxa = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],  stringsAsFactors = FALSE),
                  item_data = data.frame(item_attr = LETTERS[1:10],  stringsAsFactors = FALSE))
#|
#| ####  Taxon sorting 
test_that("Taxon data sorting works", {
  result <- arrange_taxa(obj, desc(name))
  expect_s3_class(result, "classified")
  expect_equivalent(result$taxon_data$name, rev(letters[1:5]))
})

#|
#| ####  Item sorting 
test_that("Item data sorting with taxon_funcs works", {
  result <- arrange_items(obj, desc(item_attr))
  expect_s3_class(result, "classified")
  expect_equivalent(result$item_data$item_attr, rev(LETTERS[1:10]))
})
