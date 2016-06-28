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
#' @param untaxmap (\code{logical} of length 1)
#' If \code{FALSE}, remove any untaxmap rows.
#' 
#' @return \code{\link{taxmap}}
#' 
#' @export
parse_mothur_summary <- function(file_path, untaxmap = FALSE) {
  
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
  
  # Remove 'untaxmap' rows
  if (! untaxmap) {
    untaxmap_rows <- grepl(content, pattern = "^(.*?)\\t(.*?)\\tuntaxmap\\t")
    content <- content[! untaxmap_rows]
    message(paste0("Removed ", sum(untaxmap_rows), " untaxmap rows."))
  }
  
  # Extract taxonomic data
  result <- extract_taxonomy(content[-1],
                             key = key,
                             regex = regex,
                             class_key = "name",
                             class_sep = "\\.",
                             return_input = FALSE,
                             return_match = FALSE)
  
  # Add 'all' calculated column
  result$taxon_funcs <- c(result$taxon_funcs,
                          list(all = function(obj, subset = obj$taxon_data$taxon_ids) {
                            sample_cols <- header[6:length(header)]
                            sample_cols <- sample_cols[sample_cols %in% colnames(obj$taxon_data)]
                            apply(obj$taxon_data[subset, sample_cols], MARGIN = 1, sum)
                          }))
  
  return(result)
}