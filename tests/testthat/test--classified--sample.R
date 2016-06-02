#| ## Testing column addition for `classified` objects
#|
library(metacoder)
context("Sampling`classified` objects")
#|
#| ### Add `taxon_data` columns
#|
#| ####  Code shared by tests
obj <- classified(taxon_ids = LETTERS[1:5], parent_ids = c(NA, 1, 2, 2, 1), 
                  item_taxon_ids = rep(c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4), 10),
                  taxon_data = data.frame(name = letters[1:5],
                                          stringsAsFactors = FALSE),
                  item_data = data.frame(item_attr = rep(LETTERS[1:10], 10),
                                         stringsAsFactors = FALSE))
original_plot <- plot(obj, vertex_label = paste(taxon_ids, item_counts), vertex_color = item_counts, layout = "fr")
#|
#| ####  Sampling items by taxon_weights
test_that("Sampling items by taxon_weights works", {
  result <- sample_n_items(obj, size = 10, taxon_weight = 1 / item_counts, replace = FALSE, use_supertaxa = TRUE)
  expect_s3_class(result, "classified")
})
