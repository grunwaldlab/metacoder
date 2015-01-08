## ----, echo=FALSE, warning=FALSE, message=FALSE--------------------------
library(ape)
library(metacoder)
library(knitr)
opts_chunk$set(dev='svg', fig.width = 7.1, fig.height = 7.1, cache = TRUE, warning = FALSE, message = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  file_path <- system.file("extdata", "ncbi_basidiomycetes.fasta", package = "metacoder")
#  sequences <- ape::read.FASTA(file_path)
#  genbank_ex_data <- extract_taxonomy(names(sequences),
#                              regex = "^.*\\|(.*)\\|.*\\|(.*)\\|(.*)$",
#                              key = c("item_id", acc_no = "item_info", desc = "item_info"))

## ----, echo=FALSE, comment=NA--------------------------------------------
cat(genbank_ex_data$items$input[1])

## ------------------------------------------------------------------------
taxa <- genbank_ex_data$taxa
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name)

## ----eval=FALSE----------------------------------------------------------
#  file_path <- system.file("extdata", "unite_general_release.fasta", package = "metacoder")
#  sequences <- ape::read.FASTA(file_path)
#  unite_ex_data <- extract_taxonomy(names(sequences),
#                              regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#                              key = c(name = "taxon_name", seq_id = "item_id",
#                                      other_id = "item_info", tax_string = "item_info"))

## ----, echo=FALSE, comment=NA--------------------------------------------
cat(unite_ex_data$items$input[1])

## ------------------------------------------------------------------------
taxa <- unite_ex_data$taxa
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name)

## ----eval=FALSE----------------------------------------------------------
#  file_path <- system.file("extdata", "unite_general_release.fasta", package = "metacoder")
#  sequences <- ape::read.FASTA(file_path)
#  unite_ex_data_2 <- extract_taxonomy(names(sequences),
#                              regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#                              key = c(name = "item_info", seq_id = "item_info",
#                                      other_id = "item_info", "class_name"))

## ------------------------------------------------------------------------
taxa <- unite_ex_data_2$taxa
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name)

## ----eval=FALSE----------------------------------------------------------
#  file_path <- system.file("extdata", "unite_general_release.fasta", package = "metacoder")
#  sequences <- ape::read.FASTA(file_path)
#  unite_ex_data_3 <- extract_taxonomy(names(sequences),
#                              regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#                              key = c(name = "item_info", seq_id = "item_info",
#                                      other_id = "item_info", "class_name"),
#                              database = "none")

## ------------------------------------------------------------------------
taxa <- unite_ex_data_3$taxa
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name)

## ----eval=FALSE----------------------------------------------------------
#  file_path <- system.file("extdata", "ITSoneDB_ITS1_GBandHMM.fasta", package = "metacoder")
#  sequences <- ape::read.FASTA(file_path)
#  its1_ex_data <- extract_taxonomy(names(sequences),
#                           regex = "^(.*)\\|(.*)\\|(.*)\\|(.*)$",
#                           key = c("item_id", taxon_name = "taxon_info",
#                                   "taxon_id", description = "item_info"))

## ----, echo=FALSE, comment=NA--------------------------------------------
cat(its1_ex_data$items$input[1])

## ------------------------------------------------------------------------
taxa <- its1_ex_data$taxa
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name)

## ----eval=FALSE----------------------------------------------------------
#  file_path <- system.file("extdata", "pr2_stramenopiles_gb203.fasta", package = "metacoder")
#  sequences <- ape::read.FASTA(file_path)
#  pr2_ex_data <- extract_taxonomy(names(sequences),
#                          regex = "^(.*\\..*?)\\|(.*)$",
#                          key = c("item_id", "class_name"),
#                          class_tax_sep = "|",
#                          database = "none")

## ----, echo=FALSE, comment=NA--------------------------------------------
cat(pr2_ex_data$items$input[1])

## ------------------------------------------------------------------------
taxa <- pr2_ex_data$taxa
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name)

## ----eval=FALSE----------------------------------------------------------
#  file_path <- system.file("extdata", "rdp_current_Archaea_unaligned.fa", package = "metacoder")
#  sequences <- ape::read.FASTA(file_path)
#  rdp_ex_data <- extract_taxonomy(names(sequences),
#                          regex = "^(.*?) (.*)\\tLineage=(.*)",
#                          key = c(id = "item_info", description = "item_info", "class_name"),
#                          class_tax_sep = ";",
#                          class_rank_sep = ";",
#                          class_rank_rev = TRUE)

## ----, echo=FALSE, comment=NA--------------------------------------------
cat(rdp_ex_data$items$input[212])

## ------------------------------------------------------------------------
taxa <- rdp_ex_data$taxa
plot_taxonomy(taxa$taxon_id, taxa$parent_id,
              size = taxa$item_count,
              vertex_color = taxa$item_count,
              vertex_label = taxa$name,
              min_label_size = .0215)

