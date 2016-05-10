#===================================================================================================
#' Parse summary.seqs output
#' 
#' Extract statistics from the command line output of mothur command \code{summary.seqs} and
#' return the resutls in a \code{data.frame}
#' 
#' @param text The text output of \code{summary.seqs}
#' @param file The path to saved output of \code{summary.seqs}
#' 
#' @return A \code{data.frame} of statistics
#' 
#' @keywords internal
parse_summary_seqs <- function(text = NULL, file = NULL) {
  # Parse arguments
  if (sum(missing(text), missing(file)) != 1)
    stop("Either 'text' or 'file' must be specified, but not both")
  if (!missing(file)) text = readLines(file)
  # Split into lines
  text <- unlist(strsplit(text, split = "\n", fixed = TRUE))
  # Subset the lines with statistics 
  text <- text[17:25]
  # Remove extra tab on header line
  text[1] <- gsub("\\t\\t", "\\\t", text[1])
  # Add missing tab onto 'mean' line
  text[length(text)] <- paste0(text[length(text)], "\t")
  # Parse into a data.frame
  utils::read.table(text = text, sep = "\t", header = TRUE, row.names = 1)
}




#' Parse mothur classification summary file
#' 
#' Parse mothur classification summary file
#' 
#' @param file_path (\code{character} of length 1)
#' The file path to the input file.
#' 
#' @return \code{\link{classified}}
#' 
#' @keywords internal
parse_mothur_summary <- function(file_path, unclassified = FALSE) {
  
  # Read file
  content <- readLines(file_path)
  
  # Parse header to make key
  header <- strsplit(content[[1]], split = "\t")[[1]]
  key <- c("taxon_info", "class", rep("taxon_info", length(header) - 2))
  key_names <- header
  key_names[2] <- ""
  names(key) <- key_names
  
  # Make regex
  regex <- paste0("^", paste0(collapse = "\t", rep("(.*?)", length(header))), "$")
  
  # Remove 'unclassified' rows
  if (! unclassified) {
    unclassified_rows <- grepl(content, pattern = "^(.*?)\\t(.*?)\\tunclassified\\t")
    content <- content[! unclassified_rows]
    message(paste0("Removed ", sum(unclassified_rows), " unclassified rows."))
  }
  
  # Extract taxonomic data
  extract_taxonomy(content[-1],
                   key = key,
                   regex = regex,
                   class_key = "name",
                   class_sep = "\\.",
                   return_input = FALSE,
                   return_match = FALSE)
}