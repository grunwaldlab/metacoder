#| ## Testing `extract_taxonomy`
#|
library(metacoder)
context("Extracting taxonomic information")
#|
#| ### Create dummy data
#|
#| This data has every kind of field that `extract_taxonomy` can interprete, but only one field will be used in each test.
#|
test_data = c("item_id: FJ712037.1 - name: Panthera leo - taxon_id: 9689 - class_name: Carnivora;Feliformia;Felidae;Pantherinae;Panthera - class_id: 33554;379583;9681;338153;9688",
              "item_id: KC879292.1 - name: Panthera tigris - taxon_id: 9694 - class_name: Carnivora;Feliformia;Felidae;Pantherinae;Panthera - class_id: 33554;379583;9681;338153;9688",
              "item_id: HW243304.1 - name: Ursus americanus - taxon_id: 9643 - class_name: Carnivora;Caniformia;Ursidae;Ursus - class_id: 33554;379584;9632;9639") # oh my!
test_regex <- "item_id: (.*) - name: (.*) - taxon_id: (.*) - class_name: (.*) - class_id: (.*)"

check_for_internet <- function() {
  if (! is.character(RCurl::getURL("www.google.com"))) {
    skip("Internet not available")
  }
}

#|
#| ### Exracting by item_id
test_that("Exracting by item_id works", {
  check_for_internet()
  result <- extract_taxonomy(test_data, key = "item_id", regex = "item_id: (.*?) -", database = "ncbi")
  expect_s3_class(result, "classified")
  expect_true("Eukaryota" %in% result$taxon_data$name)
})
test_that("Invalid IDs cause understandable errors", {
  check_for_internet()
  expect_error(extract_taxonomy(c("FJ712037.1", "notvalid", "HW243304.1"),
                                key = "item_id", regex = "(.*)", database = "ncbi",
                                vigilance = "error"),
               "3 item IDs failed to return classifications")
})

#|
#| ### Exracting by classification names
test_that("Exracting by classification names works", {
  result <- extract_taxonomy(test_data, key = "class", regex = "class_name: (.*?) -", 
                             class_key = "name", class_regex = "(.*)", class_sep = ";")
  expect_s3_class(result, "classified")
  expect_true("Caniformia" %in% result$taxon_data$name)
})
test_that("Looking up IDs for classification names works", {
  check_for_internet()
  result <- extract_taxonomy(c("Ascomycota_sp|AY773457|SH189849.06FU|reps|k__Fungi;p__Ascomycota;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__Ascomycota_sp", 
                               "Myrmecridium_schulzeri|EU041774|SH189850.06FU|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Incertae_sedis;f__Incertae_sedis;g__Myrmecridium;s__Myrmecridium_schulzeri", 
                               "Myrmecridium_sp|JX156014|SH189851.06FU|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Incertae_sedis;f__Incertae_sedis;g__Myrmecridium;s__Myrmecridium_sp"),
                             regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
                             key = c(seq_name = "item_info", sequence_id = "item_info",
                                     other_id = "item_info", "class"),
                             class_regex = "^(.*)__(.*)$",
                             class_key = c(unite_rank = "taxon_info", "name"),
                             class_sep = ";",
                             database = "ncbi")
  expect_s3_class(result, "classified")
  expect_equal(result$taxon_data$ncbi_id[1], "4751")
})




#|
#| ### Exracting by classification IDs
test_that("Exracting by classification IDs works", {
  result <- extract_taxonomy(test_data, key = "class", regex = "class_id: (.*?)$", 
                             class_key = "taxon_id", class_regex = "(.*)", class_sep = ";")
  expect_s3_class(result, "classified")
  expect_true("379583" %in% result$taxon_data$taxon_id)
})

#|
#| ### Exracting by taxon name
test_that("Exracting by taxon name works", {
  check_for_internet()
  result <- extract_taxonomy(test_data, key = "name", regex = "name: (.*?) - ", database = "ncbi")
  expect_s3_class(result, "classified")
  expect_true("379583" %in% result$taxon_data$ncbi_id)
})

#|
#| ### Exracting by taxon ID
test_that("Exracting by taxon ID works", {
  check_for_internet()
  result <- extract_taxonomy(test_data, key = "taxon_id", regex = "taxon_id: (.*?) - ", database = "ncbi")
  expect_s3_class(result, "classified")
  expect_true("Eukaryota" %in% result$taxon_data$name)
})

#|
#| ### Taxon info columns
test_that("Taxon info columns from key are added", {
  result <- extract_taxonomy(c("taxon_id: 9688 - class_name: Pantherinae;Panthera - class_id: 33554;379583;9681;338153;9688",
                               "taxon_id: 9688 - class_name: Pantherinae;Panthera - class_id: 33554;379583;9681;338153;9688",
                               "taxon_id: 9639 - class_name: Caniformia;Ursidae;Ursus - class_id: 33554;379584;9632;9639"), 
                             key = c("taxon_info", "class", my_custom_name = "taxon_info"), 
                             regex = "taxon_id: (.*?) - class_name: (.*) - class_id: (.*)", 
                             class_key = "name", class_regex = "(.*)", class_sep = ";")
  expect_s3_class(result, "classified")
  expect_true(all(c("taxon_info", "my_custom_name") %in% colnames(result$taxon_data)))
  expect_equal(result$taxon_data$taxon_info, c(NA, NA, NA, "9639", "9688"))
})
test_that("Taxon info columns from class are added", {
  result <- extract_taxonomy(c("taxon_id: 9688 - class_name: Pantherinae-a-1;Panthera-b-2 - class_id: 33554;379583;9681;338153;9688",
                               "taxon_id: 9688 - class_name: Pantherinae-a-1;Panthera-b-2 - class_id: 33554;379583;9681;338153;9688",
                               "taxon_id: 9639 - class_name: Caniformia-a-1;Ursidae-b-2;Ursus-c-3 - class_id: 33554;379584;9632;9639"), 
                             key = "class", 
                             regex = "- class_name: (.*) -", 
                             class_key = c("name", "taxon_info", my_custom_name = "taxon_info"),
                             class_regex = "^(.*)-(.*)-(.*)$", class_sep = ";")
  expect_s3_class(result, "classified")
  expect_true(all(c("taxon_info", "my_custom_name") %in% colnames(result$taxon_data)))
  expect_equal(result$taxon_data$taxon_info, c("a", "a", "b", "c", "b"))
})
test_that("Taxon info columns from both key and class are added", {
  result <- extract_taxonomy(c("taxon_id: 9688 - class_name: Pantherinae-a-1;Panthera-b-2 - class_id: 33554;379583;9681;338153;9688",
                               "taxon_id: 9688 - class_name: Pantherinae-a-1;Panthera-b-2 - class_id: 33554;379583;9681;338153;9688",
                               "taxon_id: 9639 - class_name: Caniformia-a-1;Ursidae-b-2;Ursus-c-3 - class_id: 33554;379584;9632;9639"), 
                             key = c("taxon_info", "class", my_custom_name = "taxon_info"), 
                             regex = "taxon_id: (.*?) - class_name: (.*) - class_id: (.*)", 
                             class_key = c("name", "taxon_info", my_custom_name = "taxon_info"),
                             class_regex = "^(.*)-(.*)-(.*)$", class_sep = ";")
  expect_s3_class(result, "classified")
  expect_true(all(c("taxon_info_1", "my_custom_name_1", "taxon_info_2", "my_custom_name_2") %in% colnames(result$taxon_data)))
  expect_equal(result$taxon_data$taxon_info_1, c("a", "a", "b", "c", "b"))
  expect_equal(result$taxon_data$taxon_info_2, c(NA, NA, NA, "9639", "9688"))
})

#|
#| ### Item info columns
test_that("Item info columns are added", {
  result <- extract_taxonomy(test_data,
                             key = c("item_info", my_custom_name = "item_info", "class", "item_info"),
                             regex = "item_id: (.*?) - name.* taxon_id: (.*?) - class_name: (.*) - class_id: (.*)", 
                             class_key = "name", class_regex = "(.*)", class_sep = ";")
  expect_s3_class(result, "classified")
  expect_equal(result$item_data$my_custom_name, c("9689", "9694", "9643"))
  expect_equal(result$item_data$item_info_1, c("FJ712037.1", "KC879292.1", "HW243304.1"))
})


#|
#| ### Invalid keys give warnings
test_that("Invalid keys give warnings", {
  expect_error(extract_taxonomy(test_data,
                                 key = c("item_info", "invalid", "class", "item_info"),
                                 regex = "item_id: (.*?) - name.* taxon_id: (.*?) - class_name: (.*) - class_id: (.*)", 
                                 class_key = "name", class_regex = "(.*)", class_sep = ";"),
               'Invalid key value "invalid" given.')
  
})
#|
#| ### Only specified keys can be duplicated
test_that("Only specified keys can be duplicated", {
  expect_error(extract_taxonomy(test_data,
                                key = c("item_info", "class", "class", "item_info"),
                                regex = "item_id: (.*?) - name.* taxon_id: (.*?) - class_name: (.*) - class_id: (.*)", 
                                class_key = "name", class_regex = "(.*)", class_sep = ";"),
               "have been used more than once:")
})

