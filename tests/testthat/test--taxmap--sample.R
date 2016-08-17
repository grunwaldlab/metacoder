#| ## Testing column addition for `taxmap` objects
#|
library(metacoder)
context("Sampling`taxmap` objects")
#|
#| ### Add `taxon_data` columns
#|
#| ####  Code shared by tests
obj <- taxmap(taxon_ids = LETTERS[1:5], supertaxon_ids = c(NA, 1, 2, 2, 1), 
                  obs_taxon_ids = rep(c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4), 10),
                  taxon_data = data.frame(name = letters[1:5],
                                          stringsAsFactors = FALSE),
                  obs_data = data.frame(obs_attr = rep(LETTERS[1:10], 10),
                                         stringsAsFactors = FALSE))
original_plot <- heat_tree(obj, node_label = paste(taxon_ids, n_obs), node_color = n_obs, layout = "fr")
#|
#| ####  Sampling observations by taxon_weights
test_that("Sampling observations by taxon_weights works", {
  result <- sample_n_obs(obj, size = 10, taxon_weight = 1 / n_obs, replace = FALSE, use_supertaxa = TRUE)
  expect_s3_class(result, "taxmap")
  expect_equal(nrow(result$obs_data), 10)
})
