## ----, echo=FALSE, warning=FALSE, message=FALSE--------------------------
library(ape)
library(metacoder)
library(knitr)
opts_chunk$set(dev='svg', fig.width = 7.1, fig.height = 7.1, cache = TRUE, warning=FALSE, message=FALSE)

## ----parse_data----------------------------------------------------------
library(metacoder)
file_path <- system.file("extdata", "unite_general_release.fasta", package = "metacoder")
sequences <- ape::read.FASTA(file_path) # Reads an example FASTA file
names(sequences)[1] # Print an example of the sequence header format
data <- extract_taxonomy(names(sequences), # Parse the sequence headers
                         regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
                         key = c(name = "item_info", seq_id = "item_info",
                                 other_id = "item_info", "class_name"),
                         database = "none")
items <- data$items # Rename parts of the result for later convenience
taxa <- data$taxa # Rename parts of the result for later convenience

## ----all_data------------------------------------------------------------
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name)

## ------------------------------------------------------------------------
index <- taxonomic_sample(root_id = "1", item_ids = items$taxon_id, taxon_ids = taxa$taxon_id,
                          parent_ids = taxa$parent_id, ranks = taxa$rank, max_counts = c("o" = 20, "s" = 5),
                          min_counts = c("s" = 3))
sampled_items <- items[index, ] # Subsample sequence metadata
sampled_taxa_index <- restrict_taxonomy(taxa = taxa$taxon_id, # Find which taxa are still used in the subset
                                         parents = taxa$parent_id,
                                         subset = sampled_items$taxon_id)
sampled_taxa <- taxa[sampled_taxa_index, ] # Subsample the taxa to only those in the subset of items 
sampled_taxa$item_count <- get_taxon_count(sampled_taxa$taxon_id, # Count how many items are in each taxon
                                           sampled_taxa$parent_id,
                                           sampled_items$taxon_id)

## ----subset_1------------------------------------------------------------
plot_taxonomy(sampled_taxa$taxon_id, sampled_taxa$parent_id,
              size = sampled_taxa$item_count,
              vertex_color = sampled_taxa$item_count,
              vertex_label = ifelse(sampled_taxa$rank %in% c("s", "o"), sampled_taxa$item_count, NA),
              line_label = ifelse(sampled_taxa$rank == "o", sampled_taxa$name, NA))

## ----, results='hide', eval=FALSE----------------------------------------
#  sequence_data <- ncbi_taxon_sample(name = "fungi", target_level = "phylum",
#                             max_counts = c(phylum = 30),
#                             entrez_query = "18S[All Fields] AND 28S[All Fields]",
#                             min_length = 600, max_length = 10000)
#  ncbi_fungi_sample <- extract_taxonomy(sequence_data$gi_no, regex = "(.*)", key = "item_id")

## ------------------------------------------------------------------------
taxa <- ncbi_fungi_sample$taxa # Rename parts of the result for later convenience
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = ifelse(taxa$rank == "phylum", taxa$name, taxa$item_count),
              min_label_size = .01)

