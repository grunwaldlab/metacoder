library(metacoder)
library(testthat)
context("Tree plotting")


my_data <- parse_tax_data(hmp_otus[1:5, ], class_cols = "lineage", class_sep = ";",
                       class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                       class_regex = "^(.+)__(.+)$")


test_that("basic tree plotting works", {
  expect_true(!is.null(heat_tree(my_data,
                                 node_size = n_obs,
                                 node_color = n_obs,
                                 node_label = taxon_names,
                                 tree_label = taxon_names,
                                 layout = "fruchterman-reingold")))
})

test_that("A single taxon can be plotted", {
  expect_true(!is.null(heat_tree(filter_taxa(my_data, 1),
                                 node_size = n_obs,
                                 node_color = n_obs,
                                 node_label = taxon_names,
                                 tree_label = taxon_names,
                                 layout = "fruchterman-reingold")))
})

test_that("Plotting with infinite inputs works ", {
  my_data$data$x <- n_obs(my_data)
  my_data$data$x[1] <- Inf
  my_data$data$x[2] <- -Inf

  expect_warning(heat_tree(my_data,
                           node_label = taxon_names,
                           node_size = x,
                           node_color = x))
})


if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}