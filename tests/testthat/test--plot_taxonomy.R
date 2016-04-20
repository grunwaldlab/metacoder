library(metacoder)
context("Tree plotting")

test_that("basic tree plotting works", {
  expect_true(!is.null(plot(contaminants,
                            vertex_size = item_count,
                            vertex_color = item_count,
                            vertex_label = name,
                            tree_label = name,
                            layout = "fruchterman-reingold")))
})
