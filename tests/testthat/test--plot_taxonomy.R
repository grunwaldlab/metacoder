library(metacoder)
context("Tree plotting")

test_that("basic tree plotting works", {
  expect_true(!is.null(plot(contaminants,
                            vertex_size = item_counts,
                            vertex_color = item_counts,
                            vertex_label = name,
                            tree_label = name,
                            layout = "fruchterman-reingold")))
})
test_that("A single taxon can be plotted", {
  expect_true(!is.null(plot(filter_taxa(contaminants, 1, subtaxa = FALSE),
                            vertex_size = item_counts,
                            vertex_color = item_counts,
                            vertex_label = name,
                            tree_label = name,
                            layout = "fruchterman-reingold")))
})
