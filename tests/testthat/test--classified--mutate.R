#| ## Testing column addition for `taxmap` objects
#|
library(metacoder)
context("Adding columns to `taxmap` objects")
#|
#| ### Add `taxon_data` columns
#|
#| ####  Code shared by tests
obj <- taxmap(taxon_ids = LETTERS[1:5], parent_ids = c(NA, 1, 2, 2, 1), 
                  obs_taxon_ids = c(2, 2, 1, 1, 3, 4, 5, 3, 3, 4),
                  taxon_data = data.frame(name = letters[1:5],
                                          stringsAsFactors = FALSE),
                  obs_data = data.frame(obs_attr = LETTERS[1:10],
                                         stringsAsFactors = FALSE))
#|
#| ####  Adding taxon_data columns
test_that("Adding taxon_data columns works", {
  result <- mutate_taxa(obj, new_name = toupper(name))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "name", "taxon_ids") %in% colnames(result$taxon_data)))
})
#|
#| ####  Adding taxon_data columns without using column names
test_that("Adding taxon_data columns without using column names works", {
  stored_in_var = letters[1:5]
  result <- mutate_taxa(obj, new_name = 1:5, new_name2 = stored_in_var)
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "new_name2", "taxon_ids") %in% colnames(result$taxon_data)))
  expect_equal(result$taxon_data$new_name2, stored_in_var)
})
#|
#| ####  Adding taxon_data columns by referencing new columns
test_that("Adding taxon_data columns by referencing new columns works", {
  result <- mutate_taxa(obj, new_name = toupper(name),
                        newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "name", "taxon_ids") %in% colnames(result$taxon_data)))
})

#|
#| ### Add `obs_data` columns
#|
#| ####  Adding obs_data columns
test_that("Adding obs_data columns works", {
  result <- mutate_obs(obj, new_name = tolower(obs_attr))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "obs_attr", "obs_taxon_ids") %in% colnames(result$obs_data)))
})
#|
#| ####  Adding obs_data columns without using column names
test_that("Adding obs_data columns without using column names works", {
  stored_in_var = letters[1:10]
  result <- mutate_obs(obj, new_name = 1:10, new_name2 = stored_in_var)
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "new_name2", "obs_taxon_ids") %in% colnames(result$obs_data)))
  expect_equal(result$obs_data$new_name2, stored_in_var)
})
#|
#| ####  Adding obs_data columns by referencing new columns
test_that("Adding obs_data columns by referencing new columns works", {
  result <- mutate_obs(obj, new_name = toupper(obs_attr),
                         newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "obs_attr", "obs_taxon_ids") %in% colnames(result$obs_data)))
})

#|
#| ### Replace `taxon_data` columns
#|
#| ####  Replacing taxon_data columns
test_that("Replacing taxon_data columns works", {
  result <- transmute_taxa(obj, new_name = toupper(name))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name","taxon_ids") %in% colnames(result$taxon_data)))
  expect_false(all(c("name") %in% colnames(result$taxon_data)))
})
#|
#| ####  Replacing taxon_data columns without using column names
test_that("Replacing taxon_data columns without using column names works", {
  stored_in_var = letters[1:5]
  result <- transmute_taxa(obj, new_name = 1:5, new_name2 = stored_in_var)
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "new_name2", "taxon_ids") %in% colnames(result$taxon_data)))
  expect_false(all(c("name") %in% colnames(result$taxon_data)))
  expect_equal(result$taxon_data$new_name2, stored_in_var)
})
#|
#| ####  Replacing taxon_data columns by referencing new columns
test_that("Replacing taxon_data columns by referencing new columns works", {
  result <- transmute_taxa(obj, new_name = toupper(name),
                           newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "newest_name", "taxon_ids") %in% colnames(result$taxon_data)))
  expect_false(all(c("name") %in% colnames(result$taxon_data)))
})
#|
#| ### Replace `obs_data` columns
#|
#| ####  Replacing obs_data columns
test_that("Replacing obs_data columns works", {
  result <- transmute_obs(obj, new_name = toupper(obs_attr))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name","obs_taxon_ids") %in% colnames(result$obs_data)))
  expect_false(all(c("obs_attr") %in% colnames(result$obs_data)))
})
#|
#| ####  Replacing obs_data columns without using column names
test_that("Replacing obs_data columns without using column names works", {
  stored_in_var = letters[1:10]
  result <- transmute_obs(obj, new_name = 1:10, new_name2 = stored_in_var)
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "new_name2", "obs_taxon_ids") %in% colnames(result$obs_data)))
  expect_false(all(c("obs_attr") %in% colnames(result$obs_data)))
  expect_equal(result$obs_data$new_name2, stored_in_var)
})
#|
#| ####  Replacing obs_data columns by referencing new columns
test_that("Replacing obs_data columns by referencing new columns works", {
  result <- transmute_obs(obj, new_name = toupper(obs_attr),
                           newest_name = paste0(new_name, "!"))
  expect_s3_class(result, "taxmap")
  expect_true(all(c("new_name", "newest_name", "obs_taxon_ids") %in% colnames(result$obs_data)))
  expect_false(all(c("obs_attr") %in% colnames(result$obs_data)))
})
