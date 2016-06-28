library(metacoder)
context("Tree plotting")

test_that("basic tree plotting works", {
  expect_true(!is.null(plot(contaminants,
                            node_size = n_obs,
                            node_color = n_obs,
                            node_label = name,
                            tree_label = name,
                            layout = "fruchterman-reingold")))
})
test_that("A single taxon can be plotted", {
  expect_true(!is.null(plot(filter_taxa(contaminants, 1, subtaxa = FALSE),
                            node_size = n_obs,
                            node_color = n_obs,
                            node_label = name,
                            tree_label = name,
                            layout = "fruchterman-reingold")))
})
