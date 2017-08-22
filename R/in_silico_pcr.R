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
#' 
#' 
#' @return The command generated as a character vector of length 1.
#' 
#' @seealso \code{\link{parse_primersearch}}
#' 
#' @keywords internal
run_primersearch <- function(seq_path, primer_path, mismatch = 5, output_path = tempfile(),
                             program_path = 'primersearch', dont_run = FALSE, ...) {
  # Check if primersearch is installed...
  primersearch_is_installed()
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
  primer_data <- stringr::str_match_all(primer_chunks, pattern)
  primer_data <- as.data.frame(cbind(rep(names(primer_chunks), vapply(primer_data, nrow, numeric(1))),
                       do.call(rbind, primer_data)[, -1]), stringsAsFactors = FALSE)
  # Reformat amplicon data -------------------------------------------------------------------------
  colnames(primer_data) <- c("pair_name", "amplimer", "seq_id", "name", "f_primer", "f_index",
                             "f_mismatch",  "r_primer", "r_index", "r_mismatch", "length")
  primer_data <- primer_data[, c("seq_id", "pair_name", "amplimer", "length", 
                               "f_primer", "f_index", "f_mismatch",
                               "r_primer", "r_index", "r_mismatch")]
  numeric_cols <- c("amplimer", "length","f_index", "f_mismatch",
                    "r_index", "r_mismatch", "seq_id")
  for (col in numeric_cols) primer_data[, col] <- as.numeric(primer_data[, col]) 
  return(primer_data)
} 


#' @rdname primersearch
#' @export
primersearch <- function(input, forward, reverse, mismatch = 5, ...) {
  UseMethod("primersearch")
}

#===================================================================================================
#' Use EMBOSS primersearch for in silico PCR
#' 
#' A pair of primers are aligned against a set of sequences.
#' The location of the best hits, quality of match, and predicted amplicons are returned.
#' Requires the EMBOSS tool kit (\url{http://emboss.sourceforge.net/}) to be installed.
#' 
#' @param input (\code{character})
#' @param forward (\code{character} of length 1) The forward primer sequence
#' @param reverse (\code{character} of length 1) The reverse primer sequence
#' @param mismatch An integer vector of length 1. The percentage of mismatches allowed.
#' @param ... Unused.
#' 
#' @return An object of type \code{\link{taxmap}}
#' 
#' @section Installing EMBOSS:
#' 
#' The command-line tool "primersearch" from the EMBOSS tool kit is needed to use this function.
#' How you install EMBOSS will depend on your operating system:
#' 
#' \strong{Linux:}
#' 
#' Open up a terminal and type:
#' 
#' \code{sudo apt-get install emboss}
#' 
#' \strong{Mac OSX:}
#' 
#' The easiest way to install EMBOSS on OSX is to use \href{http://brew.sh/}{homebrew}.
#' After installing homebrew, open up a terminal and type:
#' 
#' \code{brew install homebrew/science/emboss}
#' 
#' \strong{Windows:}
#' 
#' There is an installer for Windows here:
#' 
#' ftp://emboss.open-bio.org/pub/EMBOSS/windows/mEMBOSS-6.5.0.0-setup.exe
#' 
#' NOTE: This has not been tested by us yet.
#' 
#' @examples
#' \dontrun{
#' result <- primersearch(rdp_ex_data, 
#'                        forward = c("U519F" = "CAGYMGCCRCGGKAAHACC"),
#'                        reverse = c("Arch806R" = "GGACTACNSGGGTMTCTAAT"),
#'                        mismatch = 10)
#'                        
#' heat_tree(result, 
#'           node_size = n_obs,
#'           node_label = name,
#'           node_color = prop_amplified,
#'           node_color_range = c("red", "yellow", "green"),
#'           node_color_trans = "linear",
#'           node_color_interval = c(0, 1),
#'           layout = "fruchterman-reingold")
#' }
#' 
#' @method primersearch character
#' @rdname primersearch
#' @export
primersearch.character <- function(input, forward, reverse, mismatch = 5, ...) {
  
  # Write temporary fasta file for primersearch input ----------------------------------------------
  sequence_path <- tempfile("primersearch_sequence_input_", fileext = ".fasta")
  on.exit(file.remove(sequence_path))
  writeLines(text = paste0(">", seq_along(input), "\n", input),
             con = sequence_path)
    
  # Write primer file for primersearch input -------------------------------------------------------
  name_primer <- function(primer) {
    if (is.null(names(primer))) {
      to_be_named <- seq_along(primer)
    } else {
      to_be_named <- which(is.na(names(primer)) | names(primer) == "")
    }
    names(primer)[to_be_named] <- seq_along(primer)[to_be_named]
    return(primer)
  }
  forward <- name_primer(forward)
  reverse <- name_primer(reverse)
  pair_name <- paste(names(forward), names(reverse), sep = "_")
  primer_path <- tempfile("primersearch_primer_input_", fileext = ".txt")
  on.exit(file.remove(primer_path))
  primer_table <- as.data.frame(stringsAsFactors = FALSE,
                                cbind(pair_name, forward, reverse))
  utils::write.table(primer_table, primer_path,
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  
  # Run and parse primersearch ---------------------------------------------------------------------
  output_path <- run_primersearch(sequence_path, primer_path, mismatch = mismatch)
  on.exit(file.remove(output_path))
  output <- parse_primersearch(output_path)
  
  # Extract amplicon input ---------------------------------------------------------------------
  output$f_primer <- ifelse(vapply(output$f_primer, grepl, x = forward, FUN.VALUE = logical(1)), forward, reverse)
  output$r_primer <- ifelse(vapply(output$r_primer, grepl, x = reverse, FUN.VALUE = logical(1)), reverse, forward)
  output$r_index <- vapply(input[output$seq_id], nchar, numeric(1)) - output$r_index + 1
  output$amplicon <- unlist(Map(function(seq, start, end) substr(seq, start, end),
                                input[output$seq_id], output$f_index, output$r_index)) 
  return(output)
}


#' @method primersearch taxmap
#' 
#' @param sequence_col (\code{character} of length 1) The name of the column in \code{obs_data} that has the input sequences.
#' @param result_cols (\code{character}) The names of columns to include in the output.
#' By default, all output columns are included.
#' 
#' @rdname primersearch
#' @export
primersearch.taxmap <- function(input, forward, reverse, mismatch = 5,
                                    sequence_col = "sequence", result_cols = NULL, ...) {
  if (is.null(input$obs_data[[sequence_col]])) {
    stop(paste0('`sequence_col` "', sequence_col, '" does not exist. Check the input or change the value of the `sequence_col` option.'))
  }
  result <- primersearch(input = input$obs_data[[sequence_col]],
                         forward = forward, reverse = reverse, mismatch = mismatch)
  seq_id <- result$seq_id
  result <- result[, colnames(result) != "seq_id", drop = FALSE]
  pair_name <- result$pair_name
  if (!is.null(result_cols)) {
    result <- result[, result_cols, drop = FALSE]
  }
  
  overwritten_cols <-  colnames(input$obs_data)[colnames(input$obs_data) %in% colnames(result)]
  if (length(overwritten_cols) > 0) {
    warning(paste0('The following obs_data columns will be overwritten by primersearch:\n',
                   paste0(collapse = "\n", "    ", overwritten_cols)))
  }
  input$obs_data[ , colnames(result)] <- NA
  input$obs_data[seq_id, colnames(result)] <- result
  input$obs_data$amplified <- ! is.na(input$obs_data$length)
  input$taxon_funcs <- c(input$taxon_funcs,
                         count_amplified = function(obj, subset = obj$taxon_data$taxon_ids) {
                           vapply(obs(obj, subset), function(x) sum(obj$obs_data$amplified[x]), numeric(1))
                         },
                         prop_amplified = function(obj, subset = obj$taxon_data$taxon_ids) {
                           vapply(obs(obj, subset), function(x) sum(obj$obs_data$amplified[x]) / length(x), numeric(1))
                         })
  output <- input
  return(output)
}



#' Test if primersearch is installed
#' 
#' Test if primersearch is installed
#' 
#' @param must_be_installed (\code{logical} of length 1)
#' If \code{TRUE}, throw an error if primersearch is not installed.
#' 
#' @return \code{logical} of length 1
#' 
#' @keywords internal
primersearch_is_installed <- function(must_be_installed = TRUE) {
  test_result <- tryCatch(system2("primersearch", "--version", stdout = TRUE, stderr = TRUE),
                          error = function(e) e)
  is_installed <- grepl(pattern = "^EMBOSS", test_result)
  if (must_be_installed && ! is_installed) {
    stop("'primersearch' could not be found and is required for this function. Check that the EMBOSS tool kit is installed and is in the program search path. Type '?primersearch' for information on installing EMBOSS.")
  }
  return(invisible(is_installed))
}