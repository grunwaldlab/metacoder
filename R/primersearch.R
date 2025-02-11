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


#' Execute EMBOSS Primersearch
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
    system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
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
#' @export
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
                                     do.call(rbind, primer_data)[, -1, drop = FALSE]), stringsAsFactors = FALSE)
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
#' Requires the EMBOSS tool kit (\url{https://emboss.sourceforge.net/}) to be installed.
#' 
#' It can be confusing how the primer sequence relates to the binding sites on a
#' reference database sequence. A simplified diagram can help. For example, if
#' the top strand below (5' -> 3') is the database sequence, the forward primer
#' has the same sequence as the target region, since it will bind to the other
#' strand (3' -> 5') during PCR and extend on the 3' end. However, the reverse
#' primer must bind to the database strand, so it will have to be the complement
#' of the reference sequence. It also has to be reversed to make it in the
#' standard 5' -> 3' orientation. Therefore, the reverse primer must be the
#' reverse complement of its binding site on the reference sequence.
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
#' change. Databases of a single locus (e.g. Greengenes) will likely have a
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
#' The easiest way to install EMBOSS on OSX is to use \href{https://brew.sh/}{homebrew}.
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
#' \donttest{
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
#' result <- primersearch_raw(seqs, forward = f_primer, reverse = r_primer)
#' 
#' 
#' ### Real data set ###
#' 
#' # Get example FASTA file
#' fasta_path <- system.file(file.path("extdata", "silva_subset.fa"),
#'                           package = "metacoder")
#' 
#' # Parse the FASTA file as a taxmap object
#' obj <- parse_silva_fasta(file = fasta_path)
#' 
#' # Simulate PCR with primersearch
#' pcr_result <- primersearch_raw(obj$data$tax_data$silva_seq, 
#'                                forward = c("U519F" = "CAGYMGCCRCGGKAAHACC"),
#'                                reverse = c("Arch806R" = "GGACTACNSGGGTMTCTAAT"),
#'                                mismatch = 10)
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
#' n_amplified <- unlist(obj$obs_apply("pcr",
#'     function(x) length(unique(obj$data$tax_data$input[x]))))
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
#' }
#' 
#' @export
primersearch_raw <- function(input = NULL, file = NULL, forward, reverse, mismatch = 5) {
  
  # Read sequence info
  input <- parse_seq_input(input = input, file = file, output_format = "DNAbin", u_to_t = TRUE)
  
  # Write temporary fasta file for primersearch input
  sequence_path <- tempfile("primersearch_sequence_input_", fileext = ".fasta")
  on.exit(file.remove(sequence_path))
  # writeLines(text = paste0(">", seq_along(input), "\n", input),
  #            con = sequence_path)
  ape::write.FASTA(stats::setNames(input, seq_along(input)), file = sequence_path)
  
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
  output$r_start <- vapply(input[output$input], length, numeric(1)) - output$r_start - nchar(output$r_primer) + 2
  
  # Find the end index of the primer binding site
  output$f_end <- output$f_start + nchar(output$f_primer) - 1
  output$r_end <- output$r_start + nchar(output$r_primer) - 1
  
  # Extract primer matching region
  output$f_match <- unlist(Map(function(seq_index, start, end) paste0(as.character(input[seq_index])[[1]][start:end], collapse = ""),
                               output$input, output$f_start, output$f_end)) 
  output$r_match <- unlist(Map(function(seq_index, start, end) paste0(as.character(input[seq_index])[[1]][start:end], collapse = ""),
                               output$input, output$r_start, output$r_end))
  
  # Reverse complement matching region if the antisense strand is supplied
  output$f_match <- ifelse(output$f_primer == forward, output$f_match, rev_comp(output$f_match))
  output$r_match <- ifelse(output$f_primer == forward, output$r_match, rev_comp(output$r_match))
  
  # Extract amplicon input
  output$amplicon <- unlist(Map(function(seq_index, start, end) paste0(as.character(input[seq_index])[[1]][start:end], collapse = ""),
                                output$input, output$f_end + 1, output$r_start - 1)) 
  if (! "amplicon" %in% colnames(output)) { # For empty tables
    output$amplicon <- character(0)
  }
  output$product <- paste0(output$f_primer, output$amplicon, rev_comp(output$r_primer))
  
  # Change order of output table
  output <- output[, c("input", "f_primer", "r_primer", "f_mismatch", "r_mismatch",
                       "f_start", "f_end", "r_start", "r_end",
                       "f_match", "r_match", "amplicon", "product")]
  
  # Make all sequences upper case
  cols_with_seqs <- c("f_primer", "r_primer", "f_match", "r_match", "amplicon", "product")
  output[cols_with_seqs] <- lapply(output[cols_with_seqs], toupper)
  
  return(output)
}



#' Use EMBOSS primersearch for in silico PCR
#' 
#' A pair of primers are aligned against a set of sequences. A
#' \code{\link{taxmap}} object with two tables is returned: a table with
#' information for each predicted amplicon, quality of match, and predicted
#' amplicons, and a table with per-taxon amplification statistics. Requires the
#' EMBOSS tool kit (\url{https://emboss.sourceforge.net/}) to be installed.
#' 
#' It can be confusing how the primer sequence relates to the binding sites on a
#' reference database sequence. A simplified diagram can help. For example, if
#' the top strand below (5' -> 3') is the database sequence, the forward primer
#' has the same sequence as the target region, since it will bind to the other
#' strand (3' -> 5') during PCR and extend on the 3' end. However, the reverse
#' primer must bind to the database strand, so it will have to be the complement
#' of the reference sequence. It also has to be reversed to make it in the
#' standard 5' -> 3' orientation. Therefore, the reverse primer must be the
#' reverse complement of its binding site on the reference sequence.
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
#' change. Databases of a single locus (e.g. Greengenes) will likely have a
#' convention for which strand will be present, so relative to this convention,
#' there is a distinct "forward" and "reverse". However, computers dont know
#' about this convention, so the "forward" primer is whichever primer has the
#' same sequence as its binding region in the database (as opposed to the
#' reverse complement). For this reason, primersearch will redefine which primer
#' is "forward" and which is "reverse" based on how it binds the reference
#' sequence. See the example code in \code{\link{primersearch_raw}} for a
#' demonstration of this.
#' 
#' @param obj A \code{\link{taxmap}} object.
#' @param seqs The sequences to do in silico PCR on. This can be any variable in
#'   \code{obj$data} listed in \code{all_names(obj)} or an external variable. If
#'   an external variable (i.e. not in \code{obj$data}), it must be named by
#'   taxon IDs or have the same length as the number of taxa in \code{obj}.
#'   Currently, only character vectors are accepted.
#' @param clone If \code{TRUE}, make a copy of the input object and add on the results (like most R
#'   functions). If \code{FALSE}, the input will be changed without saving the result, which uses less RAM.
#' @inheritParams primersearch_raw
#' 
#' @return A copy of the input \code{\link{taxmap}} object with two tables added. One table contains amplicon information with one row per predicted amplicon with the following info:
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
#' }
#' 
#'  \describe{
#'     \item{taxon_id:}{The taxon IDs for the sequence.} 
#'     \item{seq_index:}{The index of the input sequence.} 
#'     \item{f_primer:}{The sequence of the forward primer.} 
#'     \item{r_primer:}{The sequence of the reverse primer.} 
#'     \item{f_mismatch:}{The number of mismatches on the forward primer.} 
#'     \item{r_mismatch:}{The number of mismatches on the reverse primer.} 
#'     \item{f_start:}{The start location of the forward primer.} 
#'     \item{f_end:}{The end location of the forward primer.} 
#'     \item{r_start:}{The start location of the reverse primer.}  
#'     \item{r_end:}{The end location of the reverse primer.} 
#'     \item{f_match:}{The sequence matched by the forward primer.} 
#'     \item{r_match:}{The sequence matched by the reverse primer.} 
#'     \item{amplicon:}{The sequence amplified by the primers, not including the primers.} 
#'     \item{product:}{The sequence amplified by the primers including the primers. This simulates a real PCR product.} 
#'   }
#'   
#' The other table contains per-taxon information about the PCR, with one row per taxon. It has the following columns:
#' 
#'  \describe{
#'     \item{taxon_ids:}{Taxon IDs.}
#'     \item{query_count:}{The number of sequences used as input.}
#'     \item{seq_count:}{The number of sequences that had at least one amplicon.}
#'     \item{amp_count:}{The number of amplicons. Might be more than one per sequence.}
#'     \item{amplified:}{If at least one sequence of that taxon had at least one amplicon.}
#'     \item{multiple:}{If at least one sequences had at least two amplicons.}
#'     \item{prop_amplified:}{The proportion of sequences with at least one amplicon.}
#'     \item{med_amp_len:}{The median amplicon length.}
#'     \item{min_amp_len:}{The minimum amplicon length.}
#'     \item{max_amp_len:}{The maximum amplicon length.}
#'     \item{med_prod_len:}{The median product length.}
#'     \item{min_prod_len:}{The minimum product length.}
#'     \item{max_prod_len:}{The maximum product length.}
#'   }
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
#' The easiest way to install EMBOSS on OSX is to use \href{https://brew.sh/}{homebrew}.
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
#' \donttest{
#' # Get example FASTA file
#' fasta_path <- system.file(file.path("extdata", "silva_subset.fa"),
#'                           package = "metacoder")
#' 
#' # Parse the FASTA file as a taxmap object
#' obj <- parse_silva_fasta(file = fasta_path)
#' 
#' # Simulate PCR with primersearch
#' # Have to replace Us with Ts in sequences since primersearch
#' #   does not understand Us.
#' obj <- primersearch(obj,
#'                     gsub(silva_seq, pattern = "U", replace = "T"), 
#'                     forward = c("U519F" = "CAGYMGCCRCGGKAAHACC"),
#'                     reverse = c("Arch806R" = "GGACTACNSGGGTMTCTAAT"),
#'                     mismatch = 10)
#'                            
#' # Plot what did not ampilify                          
#' obj %>%
#'   filter_taxa(prop_amplified < 1) %>%
#'   heat_tree(node_label = taxon_names, 
#'             node_color = prop_amplified, 
#'             node_color_range = c("grey", "red", "purple", "green"),
#'             node_color_trans = "linear",
#'             node_color_axis_label = "Proportion amplified",
#'             node_size = n_obs,
#'             node_size_axis_label = "Number of sequences",
#'             layout = "da", 
#'             initial_layout = "re")
#' } 
#' 
#' @importFrom rlang .data
#' @export
primersearch <- function(obj, seqs, forward, reverse, mismatch = 5, clone = TRUE) {
  # Non-standard argument evaluation
  data_used <- eval(substitute(obj$data_used(seqs)))
  sequences <- lazyeval::lazy_eval(lazyeval::lazy(seqs), data = data_used)
  
  # Make sure sequences are associated with taxon IDs
  if (is.null(names(sequences))) {  # seqs is not named
    if (length(sequences) == length(obj$taxa)) { # seqs is same length as number of taxa and not named
      message("`seq` is unnamed, so I will asssume it is in the same order as the taxa:\n  ",
              limited_print(type = "silent", prefix = "  ", obj$taxon_ids()), "\n",
              "If it is in a different order, name it by taxon IDs.")
      names(sequences) <- obj$taxon_ids()
    } else {
      stop(call. = FALSE,
           "`seqs`` is unnamed and of a different length (", length(sequences),
           ") than the number of taxa (", length(obj$taxa), "). ",
           "`seqs` must be named by taxon IDs or the same length as the number of taxa.")
    }
  } else { # seqs is named
    name_is_id <- names(sequences) %in% obj$taxon_ids()
    if (! all(name_is_id)) { # seq is named, but not by taxon ids
      stop(call. = FALSE,
           sum(! name_is_id), " of ", length(name_is_id), " taxon ids in `seqs` are invalid:\n",
           limited_print(type = "silent", prefix = "  ", names(sequences)[! name_is_id]), "\n",
           "Check that `seqs` is named by taxon IDs present in `taxon_ids(obj)`.")
      
    }
  }
  
  # Make copy of input object to construct output
  if (clone) {
    output <- obj$clone(deep = TRUE)
  } else {
    output <- obj
  }

  # Run primer search
  if ("amplicons" %in% names(output$data)) {
    warning(call. = FALSE,
            'The existing dataset "amplicons" will be overwritten.')
  }
  output$data$amplicons <- primersearch_raw(input = sequences, forward = forward,
                                            reverse = reverse, mismatch = mismatch) %>%
    dplyr::mutate(taxon_id = names(sequences)[.data$input]) %>%
    dplyr::rename(seq_index = .data$input) %>%
    dplyr::select(taxon_id , everything()) 
  
  # Make per-taxon table
  if ("tax_amp_stats" %in% names(output$data)) {
    warning(call. = FALSE,
            'The existing dataset "tax_amp_stats" will be overwritten.')
  }
  output$data$tax_amp_stats <- dplyr::tibble("taxon_id" = output$taxon_ids(),
                                             "query_count" = vapply(output$obs(sequences), length, numeric(1)),
                                             "seq_count" = vapply(output$obs("amplicons"),
                                                                    FUN.VALUE = numeric(1),
                                                                    FUN = function(i) length(unique(output$data$amplicons$seq_index[i]))),
                                             "amp_count" = vapply(output$obs("amplicons"), length, numeric(1)))
  output$data$tax_amp_stats$amplified <- output$data$tax_amp_stats$amp_count > 0
  
  # Check for multiple amplicons per sequence
  amp_per_seq_data <- output$data$amplicons %>%
    dplyr::group_by(.data$seq_index) %>%
    dplyr::count() %>%
    dplyr::mutate(taxon_id = names(sequences)[.data$seq_index],
                  multiple = .data$n > 1)
  output$data$tax_amp_stats$multiple <- unlist(output$obs_apply(amp_per_seq_data, function(i) any(amp_per_seq_data$multiple)))
  
  # Calculate proportion amplified
  output$data$tax_amp_stats$prop_amplified <- output$data$tax_amp_stats$seq_count / output$data$tax_amp_stats$query_count
  
  # Calculate amplicon length stats
  output$mutate_obs("tax_amp_stats",
                    med_amp_len = unlist(output$obs_apply("amplicons", value = "amplicon", func = function(s) {
                      if (length(s) == 0) {
                        return(NA_real_)
                      } else {
                        return(stats::median(nchar(s)))
                      }
                    })),
                    min_amp_len = unlist(output$obs_apply("amplicons", value = "amplicon", func = function(s) {
                      if (length(s) == 0) {
                        return(NA_real_)
                      } else {
                        return(min(nchar(s)))
                      }
                    })),
                    max_amp_len = unlist(output$obs_apply("amplicons", value = "amplicon", func = function(s) {
                      if (length(s) == 0) {
                        return(NA_real_)
                      } else {
                        return(max(nchar(s)))
                      }
                    })),
                    med_prod_len = unlist(output$obs_apply("amplicons", value = "product", func = function(s) {
                      if (length(s) == 0) {
                        return(NA_real_)
                      } else {
                        return(stats::median(nchar(s)))
                      }
                    })),
                    min_prod_len = unlist(output$obs_apply("amplicons", value = "product", func = function(s) {
                      if (length(s) == 0) {
                        return(NA_real_)
                      } else {
                        return(min(nchar(s)))
                      }
                    })),
                    max_prod_len = unlist(output$obs_apply("amplicons", value = "product", func = function(s) {
                      if (length(s) == 0) {
                        return(NA_real_)
                      } else {
                        return(max(nchar(s)))
                      }
                    }))
  )
  
  # Calculate consensus sequences
  # TODO
  
  return(output)
}
