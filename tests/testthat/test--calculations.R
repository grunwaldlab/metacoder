library(metacoder)
library(testthat)
context("Calculations")

# Make test data set
x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                   class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                   class_regex = "^(.+)__(.+)$")

test_that("Counting the number of samples with reads", {
  # Count samples with reads
  result <- calc_n_samples(x, data = "tax_data")
  expect_equal(colnames(result), c("taxon_id", "n_samples"))
  expect_equivalent(unlist(result[1, "n_samples"]), 17)
  
  # Return a vector instead of a table
  result <- calc_n_samples(x, data = "tax_data", drop = TRUE)
  expect_true(is.vector(result))
  
  # Only use some columns
  result <- calc_n_samples(x, data = "tax_data", cols = hmp_samples$sample_id[1:5])
  expect_equal(colnames(result), c("taxon_id", "n_samples"))
  
  # Return a count for each treatment
  result <- calc_n_samples(x, data = "tax_data", groups = hmp_samples$body_site)
  expect_equal(colnames(result), c("taxon_id", unique(hmp_samples$body_site)))
  
  # Rename output columns 
  result <- calc_n_samples(x, data = "tax_data", groups = hmp_samples$body_site,
                           out_names = c("A", "B", "C", "D", "E"))
  expect_equal(colnames(result), c("taxon_id", c("A", "B", "C", "D", "E")))
  
  # Cols can be factors
  result <- calc_n_samples(x, data = "tax_data", cols = hmp_samples$sample_id[1])
  expect_equal(result, 
               calc_n_samples(x, data = "tax_data", cols = as.factor(hmp_samples$sample_id[1])))
  
  # Cols can be logical
  expect_equal(result, 
               calc_n_samples(x, data = "tax_data",
                              cols = colnames(x$data$tax_data) == hmp_samples$sample_id[1]))
  
  # Cols can be numeric
  expect_equal(result, 
               calc_n_samples(x, data = "tax_data",
                              cols = which(colnames(x$data$tax_data) == hmp_samples$sample_id[1])))
})


test_that("Observation proportions", {
  # Calculate proportions for all numeric columns
  result <- calc_obs_props(x, "tax_data")
  expect_equal(colnames(x$data$tax_data)[-(1:3)], colnames(result)[-1])
  expect_true(all(result$`700016050` == x$data$tax_data$`700016050` / sum(x$data$tax_data$`700016050`)))
  
  # Calculate proportions for a subset of columns
  col_subset <- c("700035949", "700097855", "700100489")
  result <- calc_obs_props(x, "tax_data", cols = col_subset)
  expect_equal(col_subset, colnames(result)[-1])
  
  result <- calc_obs_props(x, "tax_data", cols = 4:6)
  expect_equal(col_subset, colnames(result)[-1])
  
  result <- calc_obs_props(x, "tax_data",
                           cols = startsWith(colnames(x$data$tax_data), "70001"))
  expect_equal(colnames(x$data$tax_data)[startsWith(colnames(x$data$tax_data), "70001")],
               colnames(result)[-1])
  
  # Including all other columns in ouput
  expect_warning(result <- calc_obs_props(x, "tax_data", other_cols = TRUE))
  expect_true(all(c("otu_id", "lineage") %in% colnames(result)))
  
  # Inlcuding specific columns in output
  result <- calc_obs_props(x, "tax_data", cols = col_subset,
                           other_cols = 2:3)
  expect_true(all(c("otu_id", "lineage") %in% colnames(result)))
  
  # Rename output columns
  result <- calc_obs_props(x, "tax_data", cols = col_subset,
                           out_names = c("a", "b", "c"))
  expect_equal(colnames(result), c("taxon_id", "a", "b", "c"))
  
  # Cols can be factors
  result <- calc_obs_props(x, data = "tax_data", cols = hmp_samples$sample_id[1])
  expect_equal(result, 
               calc_obs_props(x, data = "tax_data", cols = as.factor(hmp_samples$sample_id[1])))
  
  # Cols can be logical
  expect_equal(result, 
               calc_obs_props(x, data = "tax_data",
                              cols = colnames(x$data$tax_data) == hmp_samples$sample_id[1]))
  
  # Cols can be numeric
  expect_equal(result, 
               calc_obs_props(x, data = "tax_data",
                              cols = which(colnames(x$data$tax_data) == hmp_samples$sample_id[1])))
})


test_that("Summing counts per taxon", {
  # Calculate the taxon abundance for each numeric column (i.e. sample)
  result <- calc_taxon_abund(x, "tax_data")
  expect_equivalent(sum(x$data$tax_data$`700035949`), result$`700035949`[1])
  
  # Calculate the taxon abundance for a subset of columns
  expect_equal(calc_taxon_abund(x, "tax_data", cols = 4:5), 
               calc_taxon_abund(x, "tax_data", cols = c("700035949", "700097855")))
  
  # Calculate the taxon abundance for groups of columns (e.g. treatments)
  #  Note that we do not need to use the "cols" option for this since all
  #  numeric columns are samples in this dataset. If there were numeric columns
  #  that were not samples present in hmp_samples, the "cols" would be needed.
  result <- calc_taxon_abund(x, "tax_data", groups = hmp_samples$sex)
  expect_equal(colnames(result), c("taxon_id", "female", "male"))
  
  # Rename the output columns
  total_counts <- sum(x$data$tax_data[, hmp_samples$sample_id])
  result <- calc_taxon_abund(x, "tax_data", groups = hmp_samples$sex,
                             out_names = c("Women", "Men"))
  expect_equal(colnames(result), c("taxon_id", "Women", "Men"))
  expect_equal(total_counts, sum(result[1, c("Women", "Men")]))
  
  # Geting a total for all columns 
  result <- calc_taxon_abund(x, "tax_data", cols = hmp_samples$sample_id,
                             groups = rep("total", nrow(hmp_samples)))
  expect_equivalent(total_counts, result$total[1])
  
  # Cols can be factors
  result <- calc_taxon_abund(x, data = "tax_data", cols = hmp_samples$sample_id[1])
  expect_equal(result, 
               calc_taxon_abund(x, data = "tax_data", cols = as.factor(hmp_samples$sample_id[1])))
  
  # Cols can be logical
  expect_equal(result, 
               calc_taxon_abund(x, data = "tax_data",
                                cols = colnames(x$data$tax_data) == hmp_samples$sample_id[1]))
  
  # Cols can be numeric
  expect_equal(result, 
               calc_taxon_abund(x, data = "tax_data",
                                cols = which(colnames(x$data$tax_data) == hmp_samples$sample_id[1])))
})


test_that("Comparing groups of samples", {
  # Convert counts to proportions
  x$data$otu_table <- calc_obs_props(x, data = "tax_data", cols = hmp_samples$sample_id)
  
  # Get per-taxon counts
  x$data$tax_table <- calc_taxon_abund(x, data = "otu_table", cols = hmp_samples$sample_id)
  
  # Calculate difference between groups
  expect_warning(x$data$diff_table <- compare_groups(x, data = "tax_table",
                                                     cols = hmp_samples$sample_id,
                                                     groups = hmp_samples$body_site))
  expect_equal(nrow(x$data$diff_table),
               ncol(combn(length(unique(hmp_samples$body_site)), 2)) * nrow(x$data$tax_table))
  
})


test_that("Rarefying observation counts", {
  # Rarefy all numeric columns
  result <- rarefy_obs(x, "tax_data")
  expect_equal(length(unique(colSums(result[, hmp_samples$sample_id]))), 1)
  
  # Cols can be factors
  result <- rarefy_obs(x, data = "tax_data", cols = hmp_samples$sample_id[1])
  expect_equal(result, 
               rarefy_obs(x, data = "tax_data", cols = as.factor(hmp_samples$sample_id[1])))
  
  # Cols can be logical
  expect_equal(result, 
               rarefy_obs(x, data = "tax_data",
                          cols = colnames(x$data$tax_data) == hmp_samples$sample_id[1]))
  
  # Cols can be numeric
  expect_equal(result, 
               rarefy_obs(x, data = "tax_data",
                          cols = which(colnames(x$data$tax_data) == hmp_samples$sample_id[1])))
})


test_that("Converting low counts to zero", {
  # Calculate proportions for all numeric columns
  result <- zero_low_counts(x, "tax_data")
  expect_equal(sum(result[, hmp_samples$sample_id] == 1), 0)
  
  # Cols can be factors
  result <- zero_low_counts(x, data = "tax_data", cols = hmp_samples$sample_id[1])
  expect_equal(result, 
               zero_low_counts(x, data = "tax_data", cols = as.factor(hmp_samples$sample_id[1])))
  
  # Cols can be logical
  expect_equal(result, 
               zero_low_counts(x, data = "tax_data",
                               cols = colnames(x$data$tax_data) == hmp_samples$sample_id[1]))
  
  # Cols can be numeric
  expect_equal(result, 
               zero_low_counts(x, data = "tax_data",
                               cols = which(colnames(x$data$tax_data) == hmp_samples$sample_id[1])))
})
