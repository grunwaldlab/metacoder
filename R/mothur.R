#===================================================================================================
#' Parse summary.seqs output
#' 
#' Extract statistics from the command line output of mothur command \code{summary.seqs} and
#' return the results in a \code{data.frame}
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
