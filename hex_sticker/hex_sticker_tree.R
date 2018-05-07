library(metacoder)
x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                   class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                   class_regex = "^(.+)__(.+)$")

x$taxa[[1]]$name$name <- "metacoder"

x %>%
  filter_taxa(n_supertaxa <= 3) %>%
  heat_tree(node_size = n_obs,
            node_color = n_obs,
            # node_color_range = c("#EEEEEE", "#999999"),
            # node_label = ifelse(is_root, taxon_names, NA),
            # node_label_size_range = c(.2, .2),
            make_node_legend = FALSE,
            output_file = "hex_sticker/hex_sticker_tree.svg")
