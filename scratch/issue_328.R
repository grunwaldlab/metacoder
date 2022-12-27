library(metacoder)

x <- parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                    class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                    class_regex = "^(.+)__(.+)$")

meta <- hmp_samples[hmp_samples$body_site %in% c('Nose', 'Throat'), ]
# meta <- hmp_samples

# Convert counts to proportions
x$data$otu_table <- calc_obs_props(x, data = "tax_data", cols = meta$sample_id)

# Get per-taxon counts
x$data$tax_table <- calc_taxon_abund(x, data = "otu_table", cols = meta$sample_id)

# Calculate difference between treatments
x$data$diff_table <- compare_groups(x, data = "tax_table",
                                    cols = meta$sample_id,
                                    groups = meta$body_site)

# Plot results (might take a few minutes)
heat_tree_matrix(x,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions")
