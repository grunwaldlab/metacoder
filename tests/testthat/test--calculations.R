library(metacoder)
context("Calculations")

# Make test data set
x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                   class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                   class_regex = "^(.+)__(.+)$")

test_that("Counting the number of samples with reads", {
  # Count samples with reads
  calc_n_samples(x, dataset = "tax_data")
  
  # Return a vector instead of a table
  calc_n_samples(x, dataset = "tax_data", drop = TRUE)
  
  # Only use some columns
  calc_n_samples(x, dataset = "tax_data", cols = hmp_samples$sample_id[1:5])
  
  # Return a count for each treatment
  calc_n_samples(x, dataset = "tax_data", groups = hmp_samples$body_site)
  
  # Rename output columns 
  calc_n_samples(x, dataset = "tax_data", groups = hmp_samples$body_site,
                 out_names = c("A", "B", "C", "D", "E"))
  
  # Add results to input table
  calc_n_samples(x, dataset = "tax_data", append = TRUE)
})


test_that("Observation proportions", {
  # Calculate proportions for all numeric columns
  calc_obs_props(x, "tax_data")
  
  # Calculate proportions for a subset of columns
  calc_obs_props(x, "tax_data", cols = c("700035949", "700097855", "700100489"))
  calc_obs_props(x, "tax_data", cols = 4:6)
  calc_obs_props(x, "tax_data", cols = startsWith(colnames(x$data$tax_data), "70001"))
  
  # Including all other columns in ouput
  calc_obs_props(x, "tax_data", other_cols = TRUE)
  
  # Inlcuding specific columns in output
  calc_obs_props(x, "tax_data", cols = c("700035949", "700097855", "700100489"),
                 other_cols = 2:3)
  
  # Rename output columns
  calc_obs_props(x, "tax_data", cols = c("700035949", "700097855", "700100489"),
                 out_names = c("a", "b", "c"))
})


test_that("Summing counts per taxon", {
  # Calculate the taxon abundance for each numeric column (i.e. sample)
  calc_taxon_abund(x, "tax_data")
  
  # Calculate the taxon abundance for a subset of columns
  calc_taxon_abund(x, "tax_data", cols = 4:5)
  calc_taxon_abund(x, "tax_data", cols = c("700035949", "700097855"))
  calc_taxon_abund(x, "tax_data", cols = startsWith(colnames(x$data$tax_data), "70001"))
  
  # Calculate the taxon abundance for groups of columns (e.g. treatments)
  #  Note that we do not need to use the "cols" option for this since all
  #  numeric columns are samples in this dataset. If there were numeric columns
  #  that were not samples present in hmp_samples, the "cols" would be needed.
  calc_taxon_abund(x, "tax_data", groups = hmp_samples$sex)
  calc_taxon_abund(x, "tax_data", groups = hmp_samples$body_site)
  
  # The above example using the "cols" option, even though not needed in this case
  calc_taxon_abund(x, "tax_data", cols = hmp_samples$sample_id,
                   groups = hmp_samples$sex)
  
  # Rename the output columns
  calc_taxon_abund(x, "tax_data", cols = hmp_samples$sample_id[1:10],
                   out_names = letters[1:10])
  calc_taxon_abund(x, "tax_data", groups = hmp_samples$sex,
                   out_names = c("Women", "Men"))
  
  # Geting a total for all columns 
  calc_taxon_abund(x, "tax_data", cols = hmp_samples$sample_id,
                   groups = rep("total", nrow(hmp_samples)))
})


test_that("Comparing groups of samples", {
  # Convert counts to proportions
  x$data$otu_table <- calc_obs_props(x, dataset = "tax_data", cols = hmp_samples$sample_id)
  
  # Get per-taxon counts
  x$data$tax_table <- calc_taxon_abund(x, dataset = "otu_table", cols = hmp_samples$sample_id)
  
  # Calculate difference between groups
  x$data$diff_table <- compare_treatments(x, dataset = "tax_table",
                                          cols = hmp_samples$sample_id,
                                          groups = hmp_samples$body_site)
  
})


test_that("Rarefying observation counts", {
  # Rarefy all numeric columns
  rarefy_obs(x, "tax_data")
  
  # Rarefy a subset of columns
  rarefy_obs(x, "tax_data", cols = c("700035949", "700097855", "700100489"))
  rarefy_obs(x, "tax_data", cols = 4:6)
  rarefy_obs(x, "tax_data", cols = startsWith(colnames(x$data$tax_data), "70001"))
  
  # Including all other columns in ouput
  rarefy_obs(x, "tax_data", other_cols = TRUE)
  
  # Inlcuding specific columns in output
  rarefy_obs(x, "tax_data", cols = c("700035949", "700097855", "700100489"),
             other_cols = 2:3)
  
  # Rename output columns
  rarefy_obs(x, "tax_data", cols = c("700035949", "700097855", "700100489"),
             out_names = c("a", "b", "c"))
})


test_that("Converting low counts to zero", {
  # Calculate proportions for all numeric columns
  calc_obs_props(x, "tax_data")
  
  # Calculate proportions for a subset of columns
  calc_obs_props(x, "tax_data", cols = c("700035949", "700097855", "700100489"))
  calc_obs_props(x, "tax_data", cols = 4:6)
  calc_obs_props(x, "tax_data", cols = startsWith(colnames(x$data$tax_data), "70001"))
  
  # Including all other columns in ouput
  calc_obs_props(x, "tax_data", other_cols = TRUE)
  
  # Inlcuding specific columns in output
  calc_obs_props(x, "tax_data", cols = c("700035949", "700097855", "700100489"),
                 other_cols = 2:3)
  
  # Rename output columns
  calc_obs_props(x, "tax_data", cols = c("700035949", "700097855", "700100489"),
                 out_names = c("a", "b", "c"))
})
