#| ## Testing dplyr style methods for `classified` objects
#|
library(metacoder)
context("dplyr-style methods for `classified` objects")
#|
#| ### Filter
#|
#| #### Basic usage
test_that("Basic usage works", {
  obj <- classified(taxa = c(1, 2, 3), parents = c(NA, 1, 2), item_taxa = c(2, 2, 1, 1))
  result <- filter_taxa(obj, taxon_ranks < 2,
                        item_counts > 1,
                        c("1", "3"),
                        1:2)
})
#|
#| #### Basic use
