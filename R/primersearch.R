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
#' @param ... Additional arguments are passed to \code{primersearch}.
#' 
#' 
#' @return The command generated as a character vector of length 1.
#' 
#' @seealso \code{\link{parse_primersearch}}
#' 
#' @keywords internal
run_primersearch <- function(seq_path, primer_path, mismatch = 5,
                             output_path = tempfile(),
                             program_path = 'primersearch', ...) {
  # Check if primersearch is installed...
  primersearch_is_installed()
  extra_args <- as.list(match.call(expand.dots = FALSE))$...
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
  # Split output into chunks for each primer
  raw_output <- readLines(file_path)
  primer_indexes <- grep("Primer name ", raw_output, fixed = TRUE, value = FALSE)
  primer_chunk_id <- findInterval(seq_along(raw_output), primer_indexes)
  primer_chunks <- vapply(split(raw_output, primer_chunk_id)[-1],
                          paste, character(1), collapse = "\n")
  names(primer_chunks) <- stringr::str_match(primer_chunks, "Primer name ([^\n]*)")[,2]
  # Extract amplicon data from each chunk and combine
  pattern <- paste("Amplimer ([0-9]+)",
                   "\tSequence: ([^\n]*)",
                   "\t([^\n]*)",
                   "\t([^\n]+) hits forward strand at ([0-9]+) with ([0-9]+) mismatches",
                   "\t([^\n]+) hits reverse strand at \\[([0-9]+)\\] with ([0-9]+) mismatches",
                   "\tAmplimer length: ([0-9]+) bp", sep = '\n')
  primer_data <- stringr::str_match_all(primer_chunks, pattern)
  primer_data <- as.data.frame(cbind(rep(names(primer_chunks), vapply(primer_data, nrow, numeric(1))),
                                     do.call(rbind, primer_data)[, -1]), stringsAsFactors = FALSE)
  # Reformat amplicon data
  colnames(primer_data) <- c("pair_name", "amplimer", "input", "name", "f_primer", "f_start",
                             "f_mismatch",  "r_primer", "r_start", "r_mismatch")
  primer_data <- primer_data[, c("input",  "f_primer", "f_start", "f_mismatch",
                                 "r_primer", "r_start", "r_mismatch")]
  numeric_cols <- c("f_start", "f_mismatch", "r_start", "r_mismatch", "input")
  for (col in numeric_cols) primer_data[, col] <- as.numeric(primer_data[, col]) 
  return(primer_data)
} 


#' Use EMBOSS primersearch for in silico PCR
#' 
#' A pair of primers are aligned against a set of sequences.
#' The location of the best hits, quality of match, and predicted amplicons are returned.
#' Requires the EMBOSS tool kit (\url{http://emboss.sourceforge.net/}) to be installed.
#' 
#' It can be confusing how the primer sequence relates to the binding sites on a
#' reference database sequence. A simplified diagram can help. For example, if
#' the top strand below (5' -> 3') is the database sequence, the forward primer
#' has the same sequence as the target region, since it will bind to the other
#' strand (3' -> 5') during PCR and extend on the 3' end. However, the reverse
#' primer must bind to the database strand, so it will have to be the complement
#' of the reference sequence. It also has to be reversed to make it in the
#' standard 5' -> 3' orientation. Therefore, the reverse primer must be the
#' reverse comlement of its binding site on the reference sequence.
#' \preformatted{
#' Primer 1: 5' AAGTACCTTAACGGAATTATAG 3'
#' Primer 2: 5' GCTCCACCTACGAAACGAAT   3'
#'  
#'                                <- TAAGCAAAGCATCCACCTCG 5'
#' 5' ...AAGTACCTTAACGGAATTATAG......ATTCGTTTCGTAGGTGGAGC... 3'
#' 
#' 3' ...TTCATGGAATTGCCTTAATATC......TAAGCAAAGCATCCACCTCG... 5'
#'    5' AAGTACCTTAACGGAATTATAG ->
#'}
#' However, a database might have either the top or the bottom strand as a
#' reference sequence. Since one implies the sequence of the other, either is
#' valid, but this is another source of confusion. If we take the diagram above
#' and rotate it 180 degrees, it would mean the same thing, but which primer we would
#' want to call "forward" and which we would want to call "reverse" would
#' change. Databases of a single locus (e.g. Greengenes) will likly have a
#' convention for which strand will be present, so relative to this convention,
#' there is a distinct "forward" and "reverse". However, computers dont know
#' about this convention, so the "forward" primer is whichever primer has the
#' same sequence as its binding region in the database (as opposed to the
#' reverse complement). For this reason, primersearch will redefine which primer
#' is "forward" and which is "reverse" based on how it binds the reference
#' sequence. See the example code for a demonstration of this.
#' 
#' @inheritParams parse_seq_input
#' @param forward (\code{character} of length 1) The forward primer sequence
#' @param reverse (\code{character} of length 1) The reverse primer sequence
#' @param mismatch An integer vector of length 1. The percentage of mismatches allowed.
#' 
#' @return A table with one row per predicted amplicon with the following info:
#' 
#' \preformatted{
#'            (f_primer)
#'    5' AAGTACCTTAACGGAATTATAG ->        (r_primer)
#'                                <- TAAGCAAAGCATCCACCTCG 5'
#' 5' ...AAGTACCTTAACGGAATTATAG......ATTCGTTTCGTAGGTGGAGC... 3'
#'       ^                    ^      ^                  ^
#'    f_start              f_end   r_rtart             r_end
#'      
#'       |--------------------||----||------------------|
#'              f_match       amplicon       r_match  
#'       |----------------------------------------------|
#'                            product
#'                            
#'   
#' f_mismatch: The number of mismatches on the forward primer
#' r_mismatch: The number of mismatches on the reverse primer
#' input: The index of the input sequence
#' }
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
#' @examples
#' \dontrun{
#' 
#' ### Dummy test data set ###
#' 
#' primer_1_site <- "AAGTACCTTAACGGAATTATAG"
#' primer_2_site <- "ATTCGTTTCGTAGGTGGAGC"
#' amplicon <- "NNNAGTGGATAGATAGGGGTTCTGTGGCGTTTGGGAATTAAAGATTAGAGANNN"
#' seq_1 <- paste0("AA", primer_1_site, amplicon, primer_2_site, "AAAA")
#' seq_2 <- rev_comp(seq_1)
#' f_primer <- "ACGTACCTTAACGGAATTATAG" # Note the "C" mismatch at position 2
#' r_primer <- rev_comp(primer_2_site)
#' seqs <- c(a = seq_1, b = seq_2)
#' 
#' result <- primersearch(seqs, forward = f_primer, reverse = r_primer)
#' 
#' 
#' ### Real data set ###
#' 
#' # Get example FASTA file
#' fasta_path <- system.file(file.path("extdata", "silva_subset.fa"),
#'                           package = "metacoder")
#' 
#' # Parse the FASTA file as a taxmap object
#' obj <- parse_silva_fasta(fasta_path)
#' 
#' # Simulate PCR with primersearch
#' pcr_result <- primersearch(obj$data$silva_seq, 
#'                            forward = c("U519F" = "CAGYMGCCRCGGKAAHACC"),
#'                            reverse = c("Arch806R" = "GGACTACNSGGGTMTCTAAT"),
#'                            mismatch = 10)
#' 
#' # Add result to input table 
#' #  NOTE: We want to add a function to handle running pcr on a
#' #        taxmap object directly, but we are still trying to figure out
#' #        the best way to implement it. For now, do the following:
#' obj$data$pcr <- pcr_result
#' obj$data$pcr$taxon_id <- obj$data$tax_data$taxon_id[pcr_result$input]
#' 
#' # Visualize which taxa were amplified
#' #  This work because only amplicons are returned by `primersearch`
#' n_amplified <- obj$obs_apply("pcr",
#'                              function(x) length(unique(x)),
#'                              value = "input",
#'                              simplify = TRUE)
#' prop_amped <- n_amplified / obj$n_obs()
#' heat_tree(obj,
#'           node_label = taxon_names, 
#'           node_color = prop_amped, 
#'           node_color_range = c("grey", "red", "purple", "green"),
#'           node_color_trans = "linear",
#'           node_color_axis_label = "Proportion amplified",
#'           node_size = n_obs,
#'           node_size_axis_label = "Number of sequences",
#'           layout = "da", 
#'           initial_layout = "re")
#' 
#' }
#' 
#' @export
primersearch <- function(input = NULL, file = NULL, forward, reverse, mismatch = 5) {
  
  # Read sequence info
  input <- parse_seq_input(input = input, file = file)
  
  # Write temporary fasta file for primersearch input
  sequence_path <- tempfile("primersearch_sequence_input_", fileext = ".fasta")
  on.exit(file.remove(sequence_path))
  writeLines(text = paste0(">", seq_along(input), "\n", input),
             con = sequence_path)
  
  # Write primer file for primersearch input
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
  
  # Run and parse primersearch
  output_path <- run_primersearch(sequence_path, primer_path, mismatch = mismatch)
  on.exit(file.remove(output_path))
  output <- dplyr::as_tibble(parse_primersearch(output_path))
  
  # Standardize primer regex
  #   primersearch seems to use some non-standard regex-like thing
  fix_regex <- function(x) {
    gsub(x, pattern = "?", replacement = ".", fixed = TRUE)
  }
  output$f_primer <- fix_regex(output$f_primer)
  output$r_primer <- fix_regex(output$r_primer)
  
  # Convert from regex to ambiguity codes
  #   I think most people will expect the ambiguity codes they supplied
  output$f_primer <- ifelse(vapply(output$f_primer, grepl, x = forward, FUN.VALUE = logical(1)),
                            forward, reverse)
  output$r_primer <- ifelse(vapply(output$r_primer, grepl, x = reverse, FUN.VALUE = logical(1)),
                            reverse, forward)
  
  # Make reverse primer position relative to start of the sequence
  #   primersearch returns the index of the start on the reverse complement
  output$r_start <- vapply(input[output$input], nchar, numeric(1)) - output$r_start - nchar(output$r_primer) + 2
  
  # Find the end index of the primer binding site
  output$f_end <- output$f_start + nchar(output$f_primer) - 1
  output$r_end <- output$r_start + nchar(output$r_primer) - 1
  
  # Extract primer matching region
  output$f_match <- unlist(Map(function(seq, start, end) substr(seq, start, end),
                               input[output$input], output$f_start, output$f_end)) 
  output$r_match <- unlist(Map(function(seq, start, end) substr(seq, start, end),
                               input[output$input], output$r_start, output$r_end))
  
  # Reverse complement matching region if the antisense strand is supplied
  output$f_match <- ifelse(output$f_primer == forward, output$f_match, rev_comp(output$f_match))
  output$r_match <- ifelse(output$f_primer == forward, output$r_match, rev_comp(output$r_match))
  
  # Extract amplicon input
  output$amplicon <- unlist(Map(function(seq, start, end) substr(seq, start, end),
                                input[output$input], output$f_end + 1, output$r_start - 1)) 
  if (! "amplicon" %in% colnames(output)) { # For empty tables
    output$amplicon <- character(0)
  }
  output$product <- paste0(output$f_primer, output$amplicon, rev_comp(output$r_primer))
  
  # Change order of output table
  output <- output[, c("input", "f_primer", "r_primer", "f_mismatch", "r_mismatch",
                       "f_start", "f_end", "r_start", "r_end",
                       "f_match", "r_match", "amplicon", "product")]
  
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


#' Read sequences in an unknown format
#'
#' Read sequences in an unknown format. This is meant to parse the sequence
#' input arguments of functions like \code{\link{primersearch}}.
#' 
#' @param input (\code{character}) One of the following: 
#' \describe{
#'   \item{A character vector of sequences}{See the example below for what this
#'   looks like. The parser \code{\link{read_fasta}} produces output like this.}
#'   \item{A list of character vectors}{Each vector should have one base per element.}
#'   \item{A "DNAbin" object}{This is the result of parsers like
#'   \code{\link[ape]{read.FASTA}}.}
#'   \item{A list of "SeqFastadna" objects}{This is the result of parsers like
#'   \code{\link[seqinr]{read.fasta}}.}
#'   Either "input" or "file" must be supplied but not both.
#' }
#' @param file The path to a FASTA file containing sequences to use. Either
#'   "input" or "file" must be supplied but not both.
#' 
#' @return A named character vector of sequences
#' 
#' @keywords internal
parse_seq_input <- function(input = NULL, file = NULL) {
  # Check parameters
  if (sum(! c(is.null(file), is.null(input))) != 1) {
    stop(call. = FALSE,
         "Either `file` or `input` must be supplied, but not both.")
  }
  
  if (! is.null(file) && (! is.character(file) || length(file) != 1)) {
    stop(call. = FALSE,
         "`file` must be a character vector of length 1 that is a valid path to a file.")
  }

  # Convert to common format
  if (! is.null(file)) {
    result <- read_fasta(file)
  } else if (length(input) == 0 || class(input) == "character") {
    result <- input
  } else if (class(input) == "DNAbin") {
    result <- toupper(vapply(as.character(input), paste, character(1), collapse = ""))
  } else if (class(input[[1]]) == "SeqFastadna" || class(input) == "list") {
    result <- vapply(input, paste, character(1), collapse = "")
  } else {
    stop(paste0('Could not parse sequence information of class "', class(input), '".'),
         call. = FALSE)
  }
  
  return(result)
}