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
#' A pair of primers are aligned against a set of sequences.
#' The location of the best hits, quality of match, and predicted amplicons are returned.
#' Requires the EMBOSS tool kit (\url{http://emboss.sourceforge.net/}) to be installed.
#' 
#' @param input A list of character vectors of DNA sequence, a \code{\link[ape]{DNAbin}} object,
#' or an object of type \code{\link{classified}}
#' @param forward (\code{character} of length 1) The forward primer sequence
#' @param reverse (\code{character} of length 1) The reverse primer sequence
#' @param seq_name (\code{character}) Names of sequences
#' @param pair_name (\code{character} of length 1) The name of the primer pair.
#' @param mismatch An integer vector of length 1. The percentage of mismatches allowed.
#' @param ... Additional arguments are passed to \code{\link{run_primersearch}}.
#' 
#' @return An object of type \code{\link{classified}}
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
#' NOTE: This has not been tested by us yet.
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
#'                        forward = "CAGYMGCCRCGGKAAHACC",
#'                        reverse = "GGACTACNSGGGTMTCTAAT",
#'                        pair_name = "U519F_Arch806R",
#'                        mismatch = 10)
#'                        
#' plot(result, 
#'      vertex_size = item_count,
#'      vertex_label = name,
#'      vertex_color = prop_amplified,
#'      vertex_color_range = c("red", "yellow", "green"),
#'      vertex_color_trans = "radius",
#'      vertex_color_interval = c(0, 1),
#'      layout = "fruchterman-reingold")
#' }
#' 
#' @export
#' @rdname primersearch
primersearch <- function(...) {
  UseMethod("primersearch")
}

#' @method primersearch default
#' @export
#' @rdname primersearch
primersearch.default <- function(input, forward, reverse,
                         seq_name = NULL, pair_name = NULL, mismatch = 5, ...) {
  # Read input file if supplied --------------------------------------------------------------------
  if (class(input) == "DNAbin") {
    if (is.null(seq_name)) seq_name <- names(input)
  } else if (is.atomic(input) && all(file.exists(input))) {
    if (length(input) == 1) {
      input <- as.character(ape::read.dna(input, format = "fasta"))
    } else {
      stop("Only one input file can be used currently.")
    }
  } else {
    if (is.atomic(input)) input <- lapply(as.character(input), seqinr::s2c)
    input <- ape::as.DNAbin(input)
  }
  # Write temporary fasta file for primersearch input ----------------------------------------------
  if (!is.null(seq_name)) names(input) <- seq_name
  names(input) <- paste(seq_along(input), names(input))
  sequence_path <- tempfile("primersearch_sequence_input_", fileext = ".fasta")
  on.exit(file.remove(sequence_path))
  ape::write.dna(input, sequence_path, format = "fasta", nbcol = -1, colsep = "")    
  # Write primer file for primersearch input -------------------------------------------------------
  if (is.null(names(forward))) names(forward) <- seq_along(forward)
  if (is.null(names(reverse))) names(reverse) <- seq_along(reverse)
  if (is.null(pair_name)) pair_name <- paste(names(forward), names(reverse), sep = "__and__")
  forward <- unlist(lapply(forward, paste, collapse = ""))
  reverse <- unlist(lapply(reverse, paste, collapse = ""))
  primer_path <- tempfile("primersearch_primer_input_", fileext = ".txt")
  on.exit(file.remove(primer_path))
  utils::write.table(cbind(pair_name, forward, reverse), primer_path,
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  # Run and parse primersearch ---------------------------------------------------------------------
  output_path <- run_primersearch(sequence_path, primer_path, mismatch = mismatch, ...)
  on.exit(file.remove(output_path))
  output <- parse_primersearch(output_path)
  # Extract amplicon input ---------------------------------------------------------------------
  output$seq_id <- as.numeric(output$seq_id)
  output$amp_start <- output$forward_index + nchar(output$forward_primer)
  output$amp_end <- vapply(input[output$seq_id], length, numeric(1)) -
    (output$reverse_index + nchar(output$reverse_primer)) + 1
  output$amplicon <- unlist(Map(function(seq, start, end) paste(seq[start:end], collapse = ""),
                         as.character(input[output$seq_id]), output$amp_start, output$amp_end))
  return(output)
}


#' @param embed If \code{TRUE}, embed output in input object and return
#' 
#' @method primersearch classified
#' @export
#' @rdname primersearch
primersearch.classified <- function(input, embed = TRUE, ...) {
  if (is.null(input$item_data$sequence)) {
    stop('"primersearch" requires a column in "item_data" called "sequence" when using an object of class "classified"')
  }
  result <- primersearch(input = input$item_data$sequence,
                         seq_name = rownames(input$item_data),
                         ...)
  
  if (embed) {
    overwritten_cols <-  colnames(input$item_data)[colnames(input$item_data) %in% colnames(result)]
    if (length(overwritten_cols) > 0) {
      warning(paste0('The following item_data columns will be overwritten by primersearch:\n',
                     paste0(collapse = "\n", "    ", overwritten_cols)))
    }
    input$item_data[ , colnames(result)] <- NA
    input$item_data[result$name, colnames(result)] <- result
    input$item_data$amplified <- ! is.na(input$item_data$length)
    input$taxon_funcs <- c(input$taxon_funcs,
                           count_amplified = function(obj, subset = obj$taxon_data$taxon_ids) {
                             vapply(items(obj, subset), function(x) sum(obj$item_data$amplified[x]), numeric(1))
                           },
                           prop_amplified = function(obj, subset = obj$taxon_data$taxon_ids) {
                             vapply(items(obj, subset), function(x) sum(obj$item_data$amplified[x]) / length(x), numeric(1))
                           })
    output <- input
  } else {
    output <- result
  }
  return(output)
}



#' Test if primersearch is installed
#' 
#' Test if primersearch is installed
#' 
#' @param must_be_installed (\code{logical} of length 1) If \code{TRUE}, throw an error if
#' primersearch is not installed.
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