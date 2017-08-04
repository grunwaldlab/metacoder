# library(metacoder)
# context("Tree plotting")
# 
# test_that("basic tree plotting works", {
#   expect_true(!is.null(heat_tree(contaminants,
#                                  node_size = n_obs,
#                                  node_color = n_obs,
#                                  node_label = name,
#                                  tree_label = name,
#                                  layout = "fruchterman-reingold")))
# })
# test_that("A single taxon can be plotted", {
#   expect_true(!is.null(heat_tree(filter_taxa(contaminants, 1, subtaxa = FALSE),
#                                  node_size = n_obs,
#                                  node_color = n_obs,
#                                  node_label = name,
#                                  tree_label = name,
#                                  layout = "fruchterman-reingold")))
# })
# 
# test_that("Plotting with infinite inputs works ", {
#   data <- contaminants
#   data$taxon_data$x <- n_obs(data)
#   data$taxon_data$x[1] <- Inf
#   data$taxon_data$x[2] <- -Inf
#   
#   expect_warning(heat_tree(data, node_label = name, node_size = x, node_color = x))
# })
