## ----, echo=FALSE, warning=FALSE, message=FALSE--------------------------
library(ape)
library(metacoder)
library(knitr)
opts_chunk$set(eval = FALSE)

## ------------------------------------------------------------------------
#  sequences <- ape::read.FASTA("inst/extdata/ncbi_basidiomycetes.fasta")
#  cat(names(sequences)[1])
#  genbank <- extract_taxonomy(names(sequences)[1:3],
#                              regex = "^.*\\|(.*)\\|.*\\|(.*)\\|(.*)$",
#                              key = c("item_id", acc_no = "item_info", desc = "item_info"))

## ------------------------------------------------------------------------
#  sequences <- ape::read.FASTA("inst/extdata/unite_general_release.fasta")
#  cat(names(sequences)[1])
#  unite_1 <- extract_taxonomy(names(sequences)[1:3],
#                              regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#                              key = c(name = "taxon_name", seq_id = "item_id",
#                                      other_id = "item_info", tax_string = "taxon_info"),
#                              database = "ncbi")

## ------------------------------------------------------------------------
#  sequences <- ape::read.FASTA("inst/extdata/unite_general_release.fasta")
#  cat(names(sequences)[1])
#  unite_2 <- extract_taxonomy(names(sequences)[1:3],
#                              regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#                              key = c(name = "item_info", seq_id = "item_info",
#                                      other_id = "item_info", "class_name"))

## ------------------------------------------------------------------------
#  sequences <- ape::read.FASTA("inst/extdata/unite_general_release.fasta")
#  cat(names(sequences)[1])
#  unite_3 <- extract_taxonomy(names(sequences)[1:3],
#                              regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#                              key = c(name = "item_info", seq_id = "item_info",
#                                      other_id = "item_info", "class_name"),
#                              database = "none")

## ------------------------------------------------------------------------
#  sequences <- ape::read.FASTA("inst/extdata/ITSoneDB_ITS1_GBandHMM.fasta")
#  cat(names(sequences)[1])
#  its1 <- extract_taxonomy(names(sequences)[1:3],
#                           regex = "^(.*)\\|(.*)\\|(.*)\\|(.*)$",
#                           key = c("item_id", taxon_name = "taxon_info",
#                                   "taxon_id", description = "item_info"))

## ------------------------------------------------------------------------
#  sequences <- ape::read.FASTA("inst/extdata/pr2_stramenopiles_gb203.fasta")
#  cat(names(sequences)[1])
#  pr2 <- extract_taxonomy(names(sequences)[3:6],
#                          regex = "^(.*\\..*?)\\|(.*)$",
#                          key = c("item_id", "class_name"),
#                          class_tax_sep = "|",
#                          database = "none")

## ------------------------------------------------------------------------
#  sequences <- ape::read.FASTA("inst/extdata/rdp_release11_3_Archaea_unaligned.fa")
#  cat(names(sequences)[212])
#  rdp <- extract_taxonomy(names(sequences)[212:214],
#                          regex = "^(.*?) (.*)\\tLineage=(.*)",
#                          key = c(id = "item_info", description = "item_info", "class_name"),
#                          class_tax_sep = ";",
#                          class_rank_sep = ";",
#                          class_rank_rev = TRUE)

