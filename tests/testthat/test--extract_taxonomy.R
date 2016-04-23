#| ## Testing `extract_taxonomy`
#|
library(metacoder)
context("Extracting taxonomic information")
#|
#| ### Create dummy data
#|
#| This data has every kind of field that `extract_taxonomy` can interprete, but only one field will be used in each test.
#|
test_data = c("item_id: FJ712037.1 - taxon_name: Panthera leo - taxon_id: 9689 - class_name: Carnivora;Feliformia;Felidae;Pantherinae;Panthera - class_id: 33554;379583;9681;338153;9688",
              "item_id: KC879292.1 - taxon_name: Panthera tigris - taxon_id: 9694 - class_name: Carnivora;Feliformia;Felidae;Pantherinae;Panthera - class_id: 33554;379583;9681;338153;9688",
              "item_id: HW243304.1 - taxon_name: Ursus americanus - taxon_id: 9643 - class_name: Carnivora;Caniformia;Ursidae;Ursus - class_id: 33554;379584;9632;9639") # oh my!
test_regex <- "item_id: (.*) - taxon_name: (.*) - taxon_id: (.*) - class_name: (.*) - class_id: (.*)"
#|
#| ### 
test_that("Exracting by item_id works", {
  expect_true(TRUE)
  extract_taxonomy(test_data, key = "item_id", regex = "item_id: (.*?) -")
})

