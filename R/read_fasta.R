#' Apply a function to chunks of a file
#' 
#' Reads a file in chunks, applies a function to each of them, and returns to results of the function calls.
#' 
#' @param file_path (\code{character} of length 1) The path to a file to read. 
#' @param func (\code{function}) The function to run on each chunk of the file.
#' @param buffer_size (\code{numeric} of length 1) The number of lines in each chunk
#' @param simplify (\code{logical} of length 1) If \code{TRUE}, then the result is simplified to a vector.
#' 
#' @return \code{list} of results of \code{func}
#' 
#' @keywords internal
read_lines_apply <- function(file_path, func, buffer_size = 1000, simplify = FALSE) {
  skip <- 0
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