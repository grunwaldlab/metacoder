#| ## Testing column addition for `classified` objects
#|
library(metacoder)
context("Adding columns to `classified` objects")
#|
#| ### Add `taxon_data` columns
#|
#| ####  Code shared by tests
obj <- classified(taxa = c(1, 2, 3, 4, 5), parents = c(NA, 1, 2, 2, 1), 
                  item_taxa = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],
                                          stringsAsFactors = FALSE),
                  item_data = data.frame(item_attr = LETTERS[1:10],
                                         stringsAsFactors = FALSE))
#|
#| ####  Adding taxon_data columns
test_that("Adding taxon_data columns works", {
  result <- mutate_taxa(obj, new_name = toupper(name))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "name", "taxon_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Adding taxon_data columns without using column names
test_that("Adding taxon_data columns without using column names works", {
  stored_in_var = letters[1:5]
  result <- mutate_taxa(obj, new_name = 1:5, new_name2 = stored_in_var)
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "new_name2", "taxon_ids") %in% colnames(result$taxon_data)))
  expect_equal(result$taxon_data$new_name2, stored_in_var)
})
#|
#| ####  Adding taxon_data columns by referencing new columns
test_that("Adding taxon_data columns by referencing new columns works", {
  result <- mutate_taxa(obj, new_name = toupper(name),
                        newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "name", "taxon_ids") %in% colnames(result$taxon_data)))
})

#|
#| ### Add `item_data` columns
#|
#| ####  Adding item_data columns
test_that("Adding item_data columns works", {
  result <- mutate_items(obj, new_name = tolower(item_attr))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "item_attr", "item_taxon_ids") %in% colnames(result$item_data)))
})
#|
#| ####  Adding item_data columns without using column names
test_that("Adding item_data columns without using column names works", {
  stored_in_var = letters[1:10]
  result <- mutate_items(obj, new_name = 1:10, new_name2 = stored_in_var)
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "new_name2", "item_taxon_ids") %in% colnames(result$item_data)))
  expect_equal(result$item_data$new_name2, stored_in_var)
})
#|
#| ####  Adding item_data columns by referencing new columns
test_that("Adding item_data columns by referencing new columns works", {
  result <- mutate_items(obj, new_name = toupper(item_attr),
                         newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "item_attr", "item_taxon_ids") %in% colnames(result$item_data)))
})

#|
#| ### Replace `taxon_data` columns
#|
#| ####  Replacing taxon_data columns
test_that("Replacing taxon_data columns works", {
  result <- transmute_taxa(obj, new_name = toupper(name))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name","taxon_ids") %in% colnames(result$taxon_data)))
  expect_false(all(c("name") %in% colnames(result$taxon_data)))
})
#|
#| ####  Replacing taxon_data columns without using column names
test_that("Replacing taxon_data columns without using column names works", {
  stored_in_var = letters[1:5]
  result <- transmute_taxa(obj, new_name = 1:5, new_name2 = stored_in_var)
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "new_name2", "taxon_ids") %in% colnames(result$taxon_data)))
  expect_false(all(c("name") %in% colnames(result$taxon_data)))
  expect_equal(result$taxon_data$new_name2, stored_in_var)
})
#|
#| ####  Replacing taxon_data columns by referencing new columns
test_that("Replacing taxon_data columns by referencing new columns works", {
  result <- transmute_taxa(obj, new_name = toupper(name),
                           newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "newest_name", "taxon_ids") %in% colnames(result$taxon_data)))
  expect_false(all(c("name") %in% colnames(result$taxon_data)))
})
#|
#| ### Replace `item_data` columns
#|
#| ####  Replacing item_data columns
test_that("Replacing item_data columns works", {
  result <- transmute_items(obj, new_name = toupper(item_attr))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name","item_taxon_ids") %in% colnames(result$item_data)))
  expect_false(all(c("item_attr") %in% colnames(result$item_data)))
})
#|
#| ####  Replacing item_data columns without using column names
test_that("Replacing item_data columns without using column names works", {
  stored_in_var = letters[1:10]
  result <- transmute_items(obj, new_name = 1:10, new_name2 = stored_in_var)
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "new_name2", "item_taxon_ids") %in% colnames(result$item_data)))
  expect_false(all(c("item_attr") %in% colnames(result$item_data)))
  expect_equal(result$item_data$new_name2, stored_in_var)
})
#|
#| ####  Replacing item_data columns by referencing new columns
test_that("Replacing item_data columns by referencing new columns works", {
  result <- transmute_items(obj, new_name = toupper(item_attr),
                           newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "classified")
  expect_true(all(c("new_name", "newest_name", "item_taxon_ids") %in% colnames(result$item_data)))
  expect_false(all(c("item_attr") %in% colnames(result$item_data)))
})
