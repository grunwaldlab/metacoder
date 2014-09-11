#===================================================================================================
#' Execute EMBOSS Primerseach
#' 
#' @param seq_path A character vector of length 1. The path to the fasta file containing reference
#'   sequences to search for primer matches in.
#' @param primer_path A character vector of length 1. The path to the file containing primer pairs
#'   to match. The file should be whitespace-delimited with 3 columns: primer name, first primer
#'   sequence, and second primer sequence. 
#' @param mismatch An integer vector of length 1. The percentage of mismatches allowed.
#' @param output_path A character vector of length 1. Where the output of primersearch is saved.
#' @param program_path A character vector of length 1. The location of the primersearch binary.
#'   Ideally, it should be in your system's search path.
#' @param dont_run If TRUE, the command is generated, but not executed. This could be useful if you
#'   want to execute the command yourself.
#' @return The command generated as a character vector of length 1. 
#' @seealso \code{\link{parse_primersearch}}
#' @export
run_primersearch <- function(seq_path, primer_path, mismatch = 5, output_path = tempfile(),
                             program_path = 'primersearch', dont_run = FALSE, ...) {
  extra_args <- as.list(match.call(expand.dots=F))$...
  extra_args_string <- paste(names(extra_args), extra_args, collapse = " ", sep = " ")
  command <- gettextf('%s -seqall %s -infile %s -mismatchpercent %s -outfile %s',
                      program_path, seq_path, primer_path, mismatch, output_path)
  if (nchar(extra_args_string) > 0) {
    command <- paste(command, extra_args_string)
  }
  system(command)
  return(output_path)
}


#===================================================================================================
#' Parse EMBOSS primersearch output
#' 
#' Parses the output file from EMBOSS primersearch into a data.frame with rows corresponding to 
#' predicted amplicons and their associated information.
#' @param file_path The path to a primersearch output file.
#' @return A data frame with each row corresponding to amplicon data
#' @seealso \code{\link{run_primersearch}}
#' @importFrom stringr str_match_all
#' @importFrom plyr name_rows
#' @export
parse_primersearch <- function(file_path) {
  # Split output into chunks for each primer--------------------------------------------------------
  raw_output <- readLines(file_path)
  primer_indexes <- grep("Primer name ", raw_output, fixed=TRUE, value=FALSE)
  primer_chunk_id <- findInterval(seq_along(raw_output), primer_indexes)
  primer_chunks <- vapply(split(raw_output, primer_chunk_id)[-1],
                          paste, character(1), collapse = "\n")
  names(primer_chunks) <- str_match(primer_chunks, "Primer name ([^\n]*)")[,2]
  # Extract amplicon data from each chunk and combine ----------------------------------------------
  pattern <- paste("Amplimer ([0-9]+)",
                   "\tSequence: ([^\n]*)",
                   "\t([^\n]*)",
                   "\t([^\n]+) hits forward strand at ([0-9]+) with ([0-9]+) mismatches",
                   "\t([^\n]+) hits reverse strand at \\[([0-9]+)\\] with ([0-9]+) mismatches",
                   "\tAmplimer length: ([0-9]+) bp", sep = '\n')
  primer_data <- do.call(rbind, str_match_all(primer_chunks, pattern))[,-1]
  # Reformat amplicon data -------------------------------------------------------------------------
  primer_data <- name_rows(as.data.frame(primer_data))
  colnames(primer_data) <- c("amplimer", "sequence", "info", "forward_primer", "forward_index",
                             "forward_mismatch",  "reverse_primer", "reverse_index",
                             "reverse_mismatch", "length", "primer_pair")
  primer_data <- primer_data[c("primer_pair", "amplimer", "length", "sequence", "info",
                               "forward_primer", "forward_index", "forward_mismatch",
                               "reverse_primer", "reverse_index", "reverse_mismatch")]
  for (i in seq_along(primer_data)) primer_data[[i]] <- as.character(primer_data[[i]])
  numeric_cols <- c("amplimer", "length","forward_index", "forward_mismatch",
                    "reverse_index", "reverse_mismatch")
  for (col in numeric_cols) primer_data[[col]] <- as.numeric(primer_data[[col]])
  return(primer_data)
} 


#===================================================================================================
#' Use EMBOSS primersearch for in silico PCR 
#' 
#' @param sequence A list of character vectors of DNA sequence.
#' @param forward A character vector or list of primer sequences. If named and pair_name is not set, 
#'   the names are used to construct the primer pair names.
#' @param reverse A character vector or list of primer sequences. If named and pair_name is not set, 
#'   the names are used to construct the primer pair names.
#' @param pair_name A character vector of names for primer pairs.
#' @param mismatch An integer vector of length 1. The percentage of mismatches allowed.
#' @param ... Additional arguments are passed to \code{\link{run_primersearch}}.
#' @return Output from \code{\link{parse_primersearch}} (A dataframe)
#' @importFrom ape write.dna
#' @export
primersearch <- function(sequence, forward, reverse, pair_name = NULL, mismatch = 5, ...) {
  # Write fasta file for primersearch input---------------------------------------------------------
  sequence_path <- tempfile("primersearch_sequence_input_", fileext=".fasta")
  on.exit(file.remove(sequence_path))
  write.dna(sequence, sequence_path, format="fasta", nbcol=-1, colsep="")
  # Write primer file for primersearch input--------------------------------------------------------
  if (is.null(names(forward))) names(forward) <- seq_along(forward)
  if (is.null(names(reverse))) names(reverse) <- seq_along(reverse)
  if (is.null(pair_name)) pair_name <- paste(names(forward), names(reverse), sep="__and__")
  forward <- unlist(lapply(forward, paste, collapse=""))
  reverse <- unlist(lapply(reverse, paste, collapse=""))
  primer_path <- tempfile("primersearch_primer_input_", fileext=".txt")
  on.exit(file.remove(primer_path))
  write.table(cbind(pair_name, forward, reverse), primer_path,
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  # Run and parse primersearch ---------------------------------------------------------------------
  output_path <- run_primersearch(sequence_path, primer_path, mismatch = mismatch, ...)
  on.exit(file.remove(output_path))
  parse_primersearch(output_path)
}