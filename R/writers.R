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
#' @param obj A taxmap object
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
                            class_ids <- rev(obj$supertaxa(names(ids[i]), value = "taxon_ids",
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


#' Write an immitation of the RDP FASTA databse
#' 
#' Attempts to save taxonomic and sequence information of a taxmap object in the
#' RDP FASTA format. If the taxmap object was created using
#' [parse_rdp()], then it should be able to replicate the format
#' exactly.
#' 
#' The output file has a format like:
#' 
#' \preformatted{
#' >S000448483 Sparassis crispa; MBUH-PIRJO&ILKKA94-1587/ss5	Lineage=Root;rootrank;Fun...
#' ggattcccctagtaactgcgagtgaagcgggaagagctcaaatttaaaatctggcggcgtcctcgtcgtccgagttgtaa
#' tctggagaagcgacatccgcgctggaccgtgtacaagtctcttggaaaagagcgtcgtagagggtgacaatcccgtcttt
#' ...
#' }
#' 
#' @param obj A taxmap object
#' @param file (\code{character} of length 1) The file path to save the
#'   sequence fasta file. This is optional.
#' @param tax_names (\code{character} named by taxon ids) The names of taxa
#' @param ranks (\code{character} named by taxon ids) The ranks of taxa 
#' @param ids (\code{character} named by taxon ids) Sequence ids
#' @param info (\code{character} named by taxon ids) Info associated with
#'   sequences. In the example output shown here, this field corresponds to
#'   "Sparassis crispa; MBUH-PIRJO&ILKKA94-1587/ss5"
#' @param sequences (\code{character} named by taxon ids) Sequences
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family writers
#'   
#' @export
write_rdp <- function(obj, file, tax_names = taxon_names,
                      ranks = rdp_rank, ids = rdp_id, info = seq_name, 
                      sequences = rdp_seq) {
  # non-standard argument evaluation
  my_data_used_func <- obj$data_used # needed to avoid errors when testing for some reason
  data_used <- eval(substitute(my_data_used_func(obj, tax_names, ranks, ids,
                                                 sequences, info)))
  tax_names <- rlang::eval_tidy(rlang::enquo(tax_names), data = data_used)
  ranks <- rlang::eval_tidy(rlang::enquo(ranks), data = data_used)
  ids <- rlang::eval_tidy(rlang::enquo(ids), data = data_used)
  info <- rlang::eval_tidy(rlang::enquo(info), data = data_used)
  sequences <- rlang::eval_tidy(rlang::enquo(sequences), data = data_used)

  # Create sequence file
  headers <- vapply(seq_len(length(ids)),
                        FUN.VALUE = character(1),
                        function(i) {
                          class_ids <- rev(obj$supertaxa(names(ids[i]), value = "taxon_ids",
                                                     include_input = TRUE, simplify = TRUE))
                          my_names <- tax_names[class_ids]
                          my_ranks <- ranks[class_ids]
                          my_seq_name <- info[names(info) %in% class_ids]
                          paste0(ids[i], " ", my_seq_name, "\tLineage=", paste(my_names, my_ranks, sep = ";", collapse = ";"))
                        })
  seqinr::write.fasta(as.list(sequences), headers, file, as.string = TRUE, nbchar = 80)
}


#' Write an immitation of the Mothur taxonomy file
#' 
#' Attempts to save taxonomic information of a taxmap object in the
#' mothur `*.taxonomy` format. If the taxmap object was created using
#' [parse_mothur_taxonomy()], then it should be able to replicate the format
#' exactly.
#' 
#' The output file has a format like:
#' 
#' \preformatted{
#' AY457915	Bacteria(100);Firmicutes(99);Clostridiales(99);Johnsone...
#' AY457914	Bacteria(100);Firmicutes(100);Clostridiales(100);Johnso...
#' AY457913	Bacteria(100);Firmicutes(100);Clostridiales(100);Johnso...
#' AY457912	Bacteria(100);Firmicutes(99);Clostridiales(99);Johnsone...
#' AY457911	Bacteria(100);Firmicutes(99);Clostridiales(98);Ruminoco...
#' }
#' 
#' or...
#' 
#' \preformatted{
#' AY457915	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457914	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457913	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457912	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457911	Bacteria;Firmicutes;Clostridiales;Ruminococcus_et_rel.;...
#' }
#' 
#' @param obj A taxmap object
#' @param file (\code{character} of length 1) The file path to save the
#'   sequence fasta file. This is optional.
#' @param tax_names (\code{character} named by taxon ids) The names of taxa
#' @param ids (\code{character} named by taxon ids) Sequence ids
#' @param scores TBD
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family writers
#'   
#' @export
write_mothur_taxonomy <- function(obj, file, tax_names = taxon_names,
                                  ids = sequence_id, scores = score) {
  # non-standard argument evaluation
  my_data_used_func <- obj$data_used # needed to avoid errors when testing for some reason
  data_used <- eval(substitute(my_data_used_func(obj, tax_names, ids, scores)))
  tax_names <- rlang::eval_tidy(rlang::enquo(tax_names), data = data_used)
  ids <- rlang::eval_tidy(rlang::enquo(ids), data = data_used)
  if (length(data_used) > 2) {
    scores <- rlang::eval_tidy(rlang::enquo(scores), data = data_used)
  } else {
    scores <- NULL
  }

  # Create sequence file
  output <- vapply(seq_len(length(ids)),
                    FUN.VALUE = character(1),
                    function(i) {
                      class_ids <- rev(obj$supertaxa(names(ids[i]), value = "taxon_ids",
                                                 include_input = TRUE, simplify = TRUE))
                      my_names <- tax_names[class_ids]
                      my_scores <- scores[obj$data$class_data$input_index == i][class_ids]
                      if (is.null(scores)) {
                        line <- paste0(ids[i], "\t", paste(my_names, collapse = ";"))
                      } else {
                        line <- paste0(ids[i], "\t", paste(my_names, "(", my_scores, ")", collapse = ";", sep = ""))
                      }
                      paste0(line, ";")
                    })
  writeLines(output, file)
}


#' Write an immitation of the UNITE general FASTA databse
#' 
#' Attempts to save taxonomic and sequence information of a taxmap object in the
#' UNITE general FASTA format. If the taxmap object was created using
#' [parse_unite_general()], then it should be able to replicate the format
#' exactly.
#' 
#' The output file has a format like:
#' 
#' \preformatted{
#' >Glomeromycota_sp|KJ484724|SH523877.07FU|reps|k__Fungi;p__Glomeromycota;c__unid...
#' ATAATTTGCCGAACCTAGCGTTAGCGCGAGGTTCTGCGATCAACACTTATATTTAAAACCCAACTCTTAAATTTTGTAT...
#' ...
#' }
#' 
#' @param obj A taxmap object
#' @param file (\code{character} of length 1) The file path to save the
#'   sequence fasta file. This is optional.
#' @param tax_names (\code{character} named by taxon ids) The names of taxa
#' @param ranks (\code{character} named by taxon ids) The ranks of taxa 
#' @param sequences (\code{character} named by taxon ids) Sequences
#' @param seq_name (\code{character} named by taxon ids) Name of sequences.
#'   Usually a taxon name.
#' @param ids (\code{character} named by taxon ids) UNITE sequence ids
#' @param gb_acc (\code{character} named by taxon ids) Genbank accession
#'   numbers
#' @param type (\code{character} named by taxon ids) What type of sequence it
#'   is. Usually "rep" or "ref".
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family writers
#'   
#' @export
write_unite_general <- function(obj, file, tax_names = taxon_names,
                                ranks = unite_rank, sequences = unite_seq,
                                seq_name = organism, ids = unite_id,
                                gb_acc = acc_num, type = unite_type) {
  # non-standard argument evaluation
  my_data_used_func <- obj$data_used # needed to avoid errors when testing for some reason
  data_used <- eval(substitute(my_data_used_func(obj, tax_names, ranks, ids,
                                                 rdp_id, seq_name, gb_acc,
                                                 type, sequences)))
  tax_names <- rlang::eval_tidy(rlang::enquo(tax_names), data = data_used)
  ranks <- rlang::eval_tidy(rlang::enquo(ranks), data = data_used)
  ids <- rlang::eval_tidy(rlang::enquo(ids), data = data_used)
  acc_num <- rlang::eval_tidy(rlang::enquo(gb_acc), data = data_used)
  seq_name <- rlang::eval_tidy(rlang::enquo(seq_name), data = data_used)
  type <- rlang::eval_tidy(rlang::enquo(type), data = data_used)
  sequences <- rlang::eval_tidy(rlang::enquo(sequences), data = data_used)
  
  # Create sequence file
  headers <- vapply(seq_len(length(ids)),
                    FUN.VALUE = character(1),
                    function(i) {
                      class_ids <- rev(obj$supertaxa(names(ids[i]), value = "taxon_ids",
                                                 include_input = TRUE, simplify = TRUE))
                      my_names <- tax_names[class_ids]
                      my_ranks <- ranks[class_ids]
                      paste(sep = "|", seq_name[i], acc_num[i], ids[i], type[i], 
                            paste(my_ranks, my_names, sep = "__", collapse = ";"))
                    })
  
  seq_content <- paste0(">", headers, "\n",  toupper(sequences))
  writeLines(seq_content, file)
}

#' Write an immitation of the SILVA FASTA databse
#' 
#' Attempts to save taxonomic and sequence information of a taxmap object in the
#' SILVA FASTA format. If the taxmap object was created using
#' [parse_silva_fasta()], then it should be able to replicate the format
#' exactly.
#' 
#' The output file has a format like:
#' 
#' \preformatted{
#' >GCVF01000431.1.2369 Bacteria;Proteobacteria;Gammaproteobacteria;Oceanospiril...
#' CGUGCACGGUGGAUGCCUUGGCAGCCAGAGGCGAUGAAGGACGUUGUAGCCUGCGAUAAGCUCCGGUUAGGUGGCAAACA
#' ACCGUUUGACCCGGAGAUCUCCGAAUGGGGCAACCCACCCGUUGUAAGGCGGGUAUCACCGACUGAAUCCAUAGGUCGGU
#' ...
#' }
#' 
#' @param obj A taxmap object
#' @param file (\code{character} of length 1) The file path to save the
#'   sequence fasta file. This is optional.
#' @param tax_names (\code{character} named by taxon ids) The names of taxa
#' @param other_names (\code{character} named by taxon ids) Alternate names
#'   of taxa. Will be added after the primary name.
#' @param ids (\code{character} named by taxon ids) Sequence ids
#' @param start (\code{character}) The start position of the
#' sequence.
#' @param end (\code{character}) The end position of the
#' sequence.
#' @param sequences (\code{character} named by taxon ids) Sequences
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family writers
#'   
#' @export
write_silva_fasta <- function(obj, file, tax_names = taxon_names, 
                        other_names = other_name, ids = ncbi_id,
                        start = start_pos, end = end_pos,
                        sequences = silva_seq) {
  # non-standard argument evaluation
  my_data_used_func <- obj$data_used # needed to avoid errors when testing for some reason
  data_used <- eval(substitute(my_data_used_func(obj, tax_names, other_names, ids,
                                                 start, end, sequences, info)))
  tax_names <- rlang::eval_tidy(rlang::enquo(tax_names), data = data_used)
  other_names <- rlang::eval_tidy(rlang::enquo(other_names), data = data_used)
  ids <- rlang::eval_tidy(rlang::enquo(ids), data = data_used)
  start <- rlang::eval_tidy(rlang::enquo(start), data = data_used)
  end <- rlang::eval_tidy(rlang::enquo(end), data = data_used)
  sequences <- rlang::eval_tidy(rlang::enquo(sequences), data = data_used)
  
  # Create sequence file
  headers <- vapply(seq_len(length(ids)),
                    FUN.VALUE = character(1),
                    function(i) {
                      class_ids <- rev(obj$supertaxa(names(ids[i]), value = "taxon_ids",
                                                 include_input = TRUE, simplify = TRUE))
                      my_names <- tax_names[class_ids]
                      my_other <- other_names[class_ids]
                      my_names <- ifelse(my_other == "",
                                         my_names,
                                         paste0(my_names, " (", my_other, ")"))
                      if (length(my_names) >= 2) { # collapse species and genus
                        last_two <- c(length(my_names) - 1, length(my_names))
                        my_names <- c(my_names[-last_two], paste0(my_names[last_two], collapse = " "))
                      }
                      paste0(ids[i], ".", start[i], ".", end[i], " ",
                             paste(my_names, collapse = ";"))
                    })
  seqinr::write.fasta(as.list(sequences), headers, file, as.string = TRUE, nbchar = 80)
}