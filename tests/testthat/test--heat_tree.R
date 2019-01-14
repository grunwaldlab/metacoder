library(metacoder)
library(testthat)
context("Tree plotting")


x <- parse_tax_data(hmp_otus[1:5, ], class_cols = "lineage", class_sep = ";",
                       class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                       class_regex = "^(.+)__(.+)$")


test_that("basic tree plotting works", {
  expect_true(!is.null(heat_tree(x,
                                 node_size = n_obs,
                                 node_color = n_obs,
                                 node_label = taxon_names,
                                 tree_label = taxon_names,
                                 layout = "fruchterman-reingold")))
})

test_that("A single taxon can be plotted", {
  expect_true(!is.null(heat_tree(filter_taxa(x, 1),
                                 node_size = n_obs,
                                 node_color = n_obs,
                                 node_label = taxon_names,
                                 tree_label = taxon_names,
                                 layout = "fruchterman-reingold")))
})

test_that("Plotting with infinite inputs works ", {
  x$data$x <- n_obs(x)
  x$data$x[1] <- Inf
  x$data$x[2] <- -Inf

  expect_warning(heat_tree(x,
                           node_label = taxon_names,
                           node_size = x,
                           node_color = x))
})

x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                   class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                   class_regex = "^(.+)__(.+)$")


test_that("Argument validation works right ", {
  expect_message(heat_tree(x, node_size = unname(n_obs)), 
                 "Assuming values are in the same order as taxa.")
  expect_message(heat_tree(x, node_size = stats::setNames(unname(n_obs), letters[seq_along(x$taxon_ids())])), 
                 "is named, but not by taxon IDs")
  expect_error(heat_tree(x, node_size = unname(n_obs)[1:5]), 
               "has no taxon IDs and is a different length")
  expect_error(heat_tree(x, node_size = c(not_taxon = 1, 2)), 
               "is not named by taxon ids")
  expect_warning(heat_tree(x, node_size = stats::setNames(1:3, c("b", "c", NA))), 
                 "has NAs in its taxon IDs")
  expect_warning(heat_tree(x, node_size = c("b" = NA)), 
                 "has NAs. This might cause odd behavior")
  expect_message(heat_tree(x, node_label = stats::setNames(1:3, c("b", "c", NA))), 
                 "has NAs in its taxon IDs")
  expect_error(heat_tree(x, node_size = rep(Inf, 24)), 
               "has no finite, non-NA values")
  expect_error(heat_tree(x, node_size = letters[1:24]), 
               "is not numeric")
  expect_error(heat_tree(x, node_size_range = 1), 
               "must be of length 2")
  expect_error(heat_tree(x, node_size_range = c(2, 1)), 
               "is greater than its max")
  expect_error(heat_tree(x, node_size_trans = "invalid"), 
               "must be a function or the name of a built-in transformation function")
  expect_error(heat_tree(x, node_color_range = c("invalid")), 
               "must be hex color codes or a name returned by")
  expect_error(heat_tree(x, node_color_range = character(0)), 
               "has no values")
  expect_error(heat_tree(x, overlap_avoidance  = 1:3), 
               "has no values")
  expect_error(heat_tree(x, layout = "invalid"), 
               "argument must be one of the followin")
  expect_error(heat_tree(x, initial_layout = "invalid"), 
               "argument must be one of the followin")
  
})


if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}