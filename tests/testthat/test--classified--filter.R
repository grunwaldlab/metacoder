#| ## Testing filtering methods for `classified` objects
#|
library(metacoder)
context("filtering `classified` objects")
#|
#| ### Filtering taxa
#|
#| ####  Code shared by tests
obj <- classified(taxa = c(1, 2, 3, 4, 5), parents = c(NA, 1, 2, 2, 1), 
                  item_taxa = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],  stringsAsFactors = FALSE),
                  item_data = data.frame(item_attr = LETTERS[1:10],  stringsAsFactors = FALSE))
original_plot <- plot(obj, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
#|
#| ####  Taxon filtering with NSE
test_that("Taxon filtering with non-standard evaluation works", {
  result <- filter_taxa(obj, taxon_ranks < 3, item_counts > 1, 
                        subtaxa = FALSE, supertaxa = TRUE,
                        taxonless = FALSE, reassign = TRUE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_s3_class(result, "classified")
  expect_equivalent(result$taxon_data$name, c("a", "b"))
  expect_equivalent(unname(item_counts(result)), c(10, 7))
})
#|
#| ####  Taxon filtering with taxon_ids
test_that("Taxon filtering with taxon_ids works", {
  result <- filter_taxa(obj, c("1", "2"), 
                        subtaxa = FALSE, supertaxa = TRUE,
                        taxonless = FALSE, reassign = TRUE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
   expect_equivalent(result$taxon_data$name, c("a", "b"))
  expect_equivalent(unname(item_counts(result)), c(10, 7))
})
#|
#| ####  Taxon filtering with taxon_data indexes
test_that("Taxon filtering with taxon_data indexes works", {
  result <- filter_taxa(obj, c(1, 2), 
                        subtaxa = FALSE, supertaxa = TRUE,
                        taxonless = FALSE, reassign = TRUE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$taxon_data$name, c("a", "b"))
  expect_equivalent(unname(item_counts(result)), c(10, 7))
})
#|
#| ####  Taxon filtering with data stored in variables
test_that("Taxon filtering with data stored in variables", {
  input <- c(1, 2)
  result <- filter_taxa(obj, input, 
                        subtaxa = FALSE, supertaxa = TRUE,
                        taxonless = FALSE, reassign = TRUE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$taxon_data$name, c("a", "b"))
  expect_equivalent(unname(item_counts(result)), c(10, 7))
})
#|
#| ####  Removing items
test_that("Taxon filtering: removing items works", {
  result <- filter_taxa(obj, item_counts > 1, 
                        subtaxa = FALSE, supertaxa = TRUE,
                        taxonless = FALSE, reassign = FALSE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$taxon_data$name, c("a", "b", "c", "d"))
  expect_equivalent(unname(item_counts(result)), c(9, 7, 3, 2))
})
#|
#| ####  Adding NA to filtered items item_taxon_ids
test_that("Taxon filtering: NA items works", {
  result <- filter_taxa(obj, item_counts > 1, 
                        subtaxa = FALSE, supertaxa = TRUE,
                        taxonless = TRUE, reassign = FALSE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$taxon_data$name, c("a", "b", "c", "d"))
  expect_length(sum(is.na(item_data(result))), 1)
})
#|
#| ####  Removing supertaxa works
test_that("Taxon filtering: removing supertaxa works", {
  result <- filter_taxa(obj, taxon_ranks > 1, 
                        subtaxa = FALSE, supertaxa = FALSE,
                        taxonless = FALSE, reassign = TRUE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$taxon_data$name, c("b", "c", "d", "e"))
  expect_equal(sum(is.na(result$taxon_data$parent_ids)), 2)
})


#|
#| ### Filtering items
#|
#| ####  Removing items without removing taxa
test_that("Item filtering: filtering by logical vector", {
  original_plot
  result <- filter_items(obj, item_attr %in% LETTERS[1:5], itemless = TRUE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$item_data$item_attr, LETTERS[1:5])
  expect_equal(nrow(taxon_data(result)), nrow(taxon_data(obj)))
})
test_that("Item filtering: filtering by item data index", {
  original_plot
  result <- filter_items(obj, 1:5, itemless = TRUE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$item_data$item_attr, LETTERS[1:5])
  expect_equal(nrow(taxon_data(result)), nrow(taxon_data(obj)))
})
#|
#| ####  Removing items while removing taxa
test_that("Item filtering: filtering by logical vector", {
  original_plot
  result <- filter_items(obj, item_attr %in% LETTERS[1:5], itemless = FALSE)
  plot(result, vertex_label = paste(name, item_counts), vertex_color = item_counts, layout = "fr")
  expect_equivalent(result$item_data$item_attr, LETTERS[1:5])
  expect_equal(nrow(taxon_data(result)), 3)
})
