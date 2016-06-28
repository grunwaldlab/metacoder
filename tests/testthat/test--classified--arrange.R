#| ## Testing sorting methods for `taxmap` objects
#|
library(metacoder)
context("Sorting `taxmap` objects")
#|
#| ### Sorting taxa
#|
#| ####  Code shared by tests
obj <- taxmap(taxon_ids = c(1, 2, 3, 4, 5), supertaxon_ids = c(NA, 1, 2, 2, 1), 
                  obs_taxon_ids = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],  stringsAsFactors = FALSE),
                  obs_data = data.frame(obs_attr = LETTERS[1:10],  stringsAsFactors = FALSE))
#|
#| ####  Taxon sorting 
test_that("Taxon data sorting works", {
  result <- arrange_taxa(obj, desc(name))
  expect_s3_class(result, "taxmap")
  expect_equivalent(result$taxon_data$name, rev(letters[1:5]))
})

#|
#| ####  observation sorting 
test_that("observation data sorting with taxon_funcs works", {
  result <- arrange_obs(obj, desc(obs_attr))
  expect_s3_class(result, "taxmap")
  expect_equivalent(result$obs_data$obs_attr, rev(LETTERS[1:10]))
})
