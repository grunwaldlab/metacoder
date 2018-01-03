#' Apply a function to chunks of a file
#' 
#' Reads a file in chunks, applies a function to each of them, and returns to results of the function calls.
#' 
#' @param file_path (\code{character} of length 1) The path to a file to read. 
#' @param func (\code{function}) The function to run on each chunk of the file.
#' @param buffer_size (\code{numeric} of length 1) The number of lines in each chunk
#' @param simplify (\code{logical} of length 1) If \code{TRUE}, then the result is simplified to a vector.
#' @param skip (\code{numeric} of length 1) Where to start reading the file.
#' 
#' @return \code{list} of results of \code{func}
#' 
#' @keywords internal
read_lines_apply <- function(file_path, func, buffer_size = 1000, simplify = FALSE, skip = 0) {
  result <- list()
  while (length(chunk <- readr::read_lines(file_path, skip = skip, n_max = buffer_size)) > 0) {
    skip <- skip + buffer_size
    result <- c(result, list(func(chunk)))
  }
  if (simplify) {
    result <- unlist(result, recursive = FALSE)
  }
  return(result)
}


#' Get line numbers of FASTA headers
#' 
#' Get line numbers of FASTA headers without reading whole fasta file into RAM.
#' 
#' @param file_path (\code{character} of length 1) The path to a file to read.
#' @param buffer_size (\code{numeric} of length 1) The number of lines in each chunk.
#' @param return_headers (\code{logical} of length 1) If \code{TRUE}, name the result with the headers.
#' 
#' @return \code{numeric}
#' 
#' @keywords internal
fasta_headers <- function(file_path, buffer_size = 1000, return_headers = TRUE) {
  extract_headers <- function(chunk) {
    is_header <- grepl(pattern = "^>", chunk)
    output <- current_pos + which(is_header)
    if (return_headers) {
      headers <- chunk[is_header]
      headers <- gsub(pattern = "^>", replacement = "", headers)
      headers <- gsub(pattern = "\\n$", replacement = "", headers)
      names(output) <- headers
    }
    return(output)
  }
  current_pos <- 0
  read_lines_apply(file_path, extract_headers, buffer_size = buffer_size)
}



#' Read a FASTA file
#'
#' Reads a FASTA file. This is the FASTA parser for metacoder. It simply tries
#' to read a FASTA file into a named character vector with minimal fuss. It does
#' not do any checks for valid characters etc. Other FASTA parsers you might
#' want to consider include \code{\link[ape]{read.FASTA}} or
#' \code{\link[seqinr]{read.fasta}}.
#'
#' @param file_path (\code{character} of length 1) The path to a file to read.
#'
#' @return named \code{character} vector
#'
#' @examples
#'
#' # Get example FASTA file
#' fasta_path <- system.file(file.path("extdata", "silva_subset.fa"),
#'                           package = "metacoder")
#'
#' # Read fasta file
#' my_seqs <- read_fasta(fasta_path)
#'
#' @export
read_fasta <- function(file_path) {
  # Read raw string
  raw_data <- readr::read_file(file_path)
  
  # Return an empty vector an a warning if no sequences are found
  if (raw_data == "") {
    warning(paste0("No sequences found in the file: ", file_path))
    return(character(0))
  }
  
  # Find location of every header start 
  split_data <- stringr::str_split(raw_data, pattern = "\n>", simplify = TRUE)
  
  # Split the data for each sequence into lines
  split_data <- stringr::str_split(split_data, pattern = "\n")
  
  # The first lines are headers, so remvove those
  headers <- vapply(split_data, FUN = `[`, FUN.VALUE = character(1), 1)
  split_data <- lapply(split_data, FUN = `[`, -1)
  
  # Remove the > from the first sequence. The others were removed by the split
  headers[1] <- sub(headers[1], pattern = "^>", replacement = "")
  
  # Combine multiple lines into single sequences
  seqs <- vapply(split_data, FUN = paste0, FUN.VALUE = character(1), collapse = "")
  
  # Combine and return results 
  return(stats::setNames(seqs, headers))
}
