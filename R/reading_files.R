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
#' Read a FASTA file
#' 
#' @param file_path (\code{character} of length 1) The path to a file to read.
#' @param subset (\code{numeric}) Indexes of entries to return.
#' If not \code{NULL}, the file will first be indexed without loading the whole file into RAM.
#' 
#' @return names \code{character}
#' 
#' @export
read_fasta <- function(file_path, subset = NULL) {
  
}
