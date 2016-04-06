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
#' @param ... Additional arguments are passed to \code{primersearch}.
#' @return The command generated as a character vector of length 1.
#' @seealso \code{\link{parse_primersearch}}
#' 
#' @keywords internal
run_primersearch <- function(seq_path, primer_path, mismatch = 5, output_path = tempfile(),
                             program_path = 'primersearch', dont_run = FALSE, ...) {
  extra_args <- as.list(match.call(expand.dots=F))$...
  if (Sys.info()['sysname'] == "Windows") {
    arguments <- c("-seqall", seq_path,
                 "-infile", primer_path,
                 "-mismatchpercent", mismatch,
                 "-outfile", output_path,
                 as.character(extra_args))
    system2(program_path, arguments, stdout = TRUE, stderr = TRUE)
  } else {
    extra_args_string <- paste(names(extra_args), extra_args, collapse = " ", sep = " ")
    command <- gettextf('%s -seqall %s -infile %s -mismatchpercent %s -outfile %s',
                        program_path, seq_path, primer_path, mismatch, output_path)
    if (nchar(extra_args_string) > 0) command <- paste(command, extra_args_string)
    system(command)
    
  }
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
#' 
#' @keywords internal
parse_primersearch <- function(file_path) {
  # Split output into chunks for each primer--------------------------------------------------------
  raw_output <- readLines(file_path)
  primer_indexes <- grep("Primer name ", raw_output, fixed=TRUE, value=FALSE)
  primer_chunk_id <- findInterval(seq_along(raw_output), primer_indexes)
  primer_chunks <- vapply(split(raw_output, primer_chunk_id)[-1],
                          paste, character(1), collapse = "\n")
  names(primer_chunks) <- stringr::str_match(primer_chunks, "Primer name ([^\n]*)")[,2]
  # Extract amplicon data from each chunk and combine ----------------------------------------------
  pattern <- paste("Amplimer ([0-9]+)",
                   "\tSequence: ([^\n]*)",
                   "\t([^\n]*)",
                   "\t([^\n]+) hits forward strand at ([0-9]+) with ([0-9]+) mismatches",
                   "\t([^\n]+) hits reverse strand at \\[([0-9]+)\\] with ([0-9]+) mismatches",
                   "\tAmplimer length: ([0-9]+) bp", sep = '\n')
  primer_data <- do.call(rbind, stringr::str_match_all(primer_chunks, pattern))[,-1]
  # Reformat amplicon data -------------------------------------------------------------------------
  primer_data <- plyr::name_rows(as.data.frame(primer_data))
  colnames(primer_data) <- c("amplimer", "seq_id", "name", "forward_primer", "forward_index",
                             "forward_mismatch",  "reverse_primer", "reverse_index",
                             "reverse_mismatch", "length", "primer_pair")
  primer_data <- primer_data[c("primer_pair", "amplimer", "length", "seq_id", "name",
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
#' @param sequence A list of character vectors of DNA sequence or a file path to FASTA file 
#'   containing sequences.
#' @param forward A character vector or list of primer sequences. If named and pair_name is not set, 
#'   the names are used to construct the primer pair names.
#' @param reverse A character vector or list of primer sequences. If named and pair_name is not set, 
#'   the names are used to construct the primer pair names.
#' @param seq_name Names of sequences
#' @param pair_name A character vector of names for primer pairs.
#' @param mismatch An integer vector of length 1. The percentage of mismatches allowed.
#' @param ... Additional arguments are passed to \code{\link{run_primersearch}}.
#' @return Output from \code{\link{parse_primersearch}} (A dataframe)
#' 
#' @export
#' @rdname primersearch
primersearch <- function(...) {
  UseMethod("primersearch")
}

#' @method primersearch default
#' @export
#' @rdname primersearch
primersearch.default <- function(sequence, forward, reverse,
                         seq_name = NULL, pair_name = NULL, mismatch = 5, ...) {
  # Read input file if supplied --------------------------------------------------------------------
  if (class(sequence) == "DNAbin") {
    if (is.null(seq_name)) seq_name <- names(sequence)
  } else if (is.atomic(sequence) && all(file.exists(sequence))) {
    if (length(sequence) == 1) {
      sequence <- as.character(ape::read.dna(sequence, format = "fasta"))
    } else {
      stop("Only one input file can be used currently.")
    }
  } else {
    if (is.atomic(sequence)) sequence <- lapply(as.character(sequence), seqinr::s2c)
    sequence <- ape::as.DNAbin(sequence)
  }
  # Write temporary fasta file for primersearch input ----------------------------------------------
  if (!is.null(seq_name)) names(sequence) <- seq_name
  names(sequence) <- paste(seq_along(sequence), names(sequence))
  sequence_path <- tempfile("primersearch_sequence_input_", fileext=".fasta")
  on.exit(file.remove(sequence_path))
  ape::write.dna(sequence, sequence_path, format="fasta", nbcol=-1, colsep="")    
  # Write primer file for primersearch input -------------------------------------------------------
  if (is.null(names(forward))) names(forward) <- seq_along(forward)
  if (is.null(names(reverse))) names(reverse) <- seq_along(reverse)
  if (is.null(pair_name)) pair_name <- paste(names(forward), names(reverse), sep="__and__")
  forward <- unlist(lapply(forward, paste, collapse = ""))
  reverse <- unlist(lapply(reverse, paste, collapse = ""))
  primer_path <- tempfile("primersearch_primer_input_", fileext=".txt")
  on.exit(file.remove(primer_path))
  write.table(cbind(pair_name, forward, reverse), primer_path,
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  # Run and parse primersearch ---------------------------------------------------------------------
  output_path <- run_primersearch(sequence_path, primer_path, mismatch = mismatch, ...)
  on.exit(file.remove(output_path))
  output <- parse_primersearch(output_path)
  # Extract amplicon sequence ---------------------------------------------------------------------
  output$seq_id <- as.numeric(output$seq_id)
  output$amp_start <- output$forward_index + nchar(output$forward_primer)
  output$amp_end <- vapply(sequence[output$seq_id], length, numeric(1)) -
    (output$reverse_index + nchar(output$reverse_primer)) + 1
  output$amplicon <- unlist(Map(function(seq, start, end) paste(seq[start:end], collapse = ""),
                         as.character(sequence[output$seq_id]), output$amp_start, output$amp_end))
  return(output)
}


#' @param embed If \code{TRUE}, embed output in input object and return
#' 
#' @method primersearch classified
#' @export
#' @rdname primersearch
primersearch.classified <- function(classified_data, embed = TRUE, ...) {
  if (is.null(classified_data$item_data$sequence)) {
    stop('"primersearch" requires a column in "item_data" called "sequence" when using an object of class "classified"')
  }
  result <- primersearch(sequence = classified_data$item_data$sequence,
                         seq_name = rownames(classified_data$item_data),
                         ...)
  
  if (embed) {
    classified_data$item_data[ , colnames(result)] <- NA
    classified_data$item_data[result$name, colnames(result)] <- result
    classified_data$item_data$amplified <- ! is.na(classified_data$item_data$seq_id)
    classified_data$taxon_funcs <- c(classified_data$taxon_funcs,
                                     count_amplified = function(obj, taxon) {
                                       sum(obj$item_data$amplified, na.rm = TRUE)
                                     },
                                     prop_amplified = function(obj, taxon) {
                                       sum(obj$item_data$amplified, na.rm = TRUE) / nrow(obj$item_data)
                                     })
    output <- classified_data
  } else {
    output <- result
  }
  return(output)
}