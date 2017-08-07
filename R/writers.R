#' Write an immitation of the Greengenes databse
#' 
#' Attempts to save taxonomic and sequence information of a taxmap object in the
#' Greengenes output format. If the taxmap object was created using
#' [parse_greengenes()], then it should be able to replicate the format
#' exactly.
#' 
#' The taxonomy output file has a format like:
#' 
#' \preformatted{
#' 228054  k__Bacteria; p__Cyanobacteria; c__Synechococcophycideae; o__Synech...
#' 844608  k__Bacteria; p__Cyanobacteria; c__Synechococcophycideae; o__Synech...
#' ...
#' }
#' 
#' The optional sequence file has a format like:
#' 
#' \preformatted{
#' >1111886 
#' AACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTA...
#' >1111885 
#' AGAGTTTGATCCTGGCTCAGAATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGTACGAGAAATCCCGAGC...
#' ...
#' }
#' 
#' @param tax_file (\code{character} of length 1) The file path to save the
#'   taxonomy file.
#' @param seq_file (\code{character} of length 1) The file path to save the
#'   sequence fasta file. This is optional.
#' @param tax_names (\code{character} named by taxon ids) The names of taxa
#' @param ranks (\code{character} named by taxon ids) The ranks of taxa 
#' @param ids (\code{character} named by taxon ids) Sequence ids
#' @param sequences (\code{character} named by taxon ids) Sequences
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family writers
#'   
#' @export
write_greengenes <- function(obj, tax_file = NULL, seq_file = NULL,
                             tax_names = taxon_names, ranks = gg_rank, ids = gg_id,
                             sequences = gg_seq) {
  # non-standard argument evaluation
  my_data_used_func <- obj$data_used # needed to avoid errors when testing for some reason
  data_used <- eval(substitute(my_data_used_func(obj, tax_names, ranks, ids, sequences)))
  tax_names <- rlang::eval_tidy(rlang::enquo(tax_names), data = data_used)
  ranks <- rlang::eval_tidy(rlang::enquo(ranks), data = data_used)
  ids <- rlang::eval_tidy(rlang::enquo(ids), data = data_used)
  sequences <- rlang::eval_tidy(rlang::enquo(sequences), data = data_used)
  
  # Check that at least one output file is specified
  if (is.null(tax_file) && is.null(seq_file)) {
    stop('"tax_file" and/or "seq_file" must be specified for anything to be written.',
         call. = FALSE)
  }
  
  # Create taxonomy file
  if (!is.null(tax_file)) {
    tax_content <- vapply(seq_len(length(ids)),
                          FUN.VALUE = character(1),
                          function(i) {
                            class_ids <- rev(supertaxa(obj, names(ids[i]), value = "taxon_ids",
                                                       include_input = TRUE, simplify = TRUE))
                            my_names <- tax_names[class_ids]
                            my_ranks <- ranks[class_ids]
                            if (length(class_ids) < 7) {
                              my_names <- c(my_names, rep("", 7 - length(class_ids)))
                              default_ranks <- c("k", "p", "c", "o", "f", "g", "s")
                              my_ranks <- c(my_ranks, default_ranks[(length(class_ids) + 1):7])
                            }
                            paste0(ids[i], "\t", paste(my_ranks, my_names, sep = "__", collapse = "; "))
                          })
    writeLines(tax_content, tax_file)
  }
  
  # Create sequence file
  if (!is.null(seq_file)) {
    seq_content <- paste0(">", ids, "\n", sequences)
    seq_content <- seq_content[order(ids, decreasing = TRUE)]
    writeLines(seq_content, seq_file)
  }
  
}