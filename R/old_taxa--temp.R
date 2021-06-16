
# file_path <- system.file("extdata", "unite_general_release.fasta", package = "metacoder")
# sequences <- ape::read.FASTA(file_path)
# library(taxa) # The parsers in taxa are used
# unite_ex_data_2 <- extract_tax_data(names(sequences)[1],
#                                     regex = "^(.*)\\|(.*)\\|(.*)\\|.*\\|(.*)$",
#                                     key = c(seq_name = "info", seq_id = "info",
#                                             other_id = "info", my_class = "class"),
#                                     class_regex = "^(.*)__(.*)$",
#                                     class_key = c(unite_rank = "info", my_name = "taxon_name"),
#                                     class_sep = ";",
#                                     database = "ncbi")
