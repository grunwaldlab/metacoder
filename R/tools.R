#' Revere complement sequences
#'
#' Make the reverse complement of one or more sequences stored as a character
#' vector. This is a wrapper for \code{\link[seqinr]{comp}} for character
#' vectors instead of lists of character vectors with one value per letter.
#' IUPAC ambiguity codes are handled and the upper/lower case is preserved.
#' 
#' @param seqs A character vector with one element per sequence.
#' 
#' @family sequence transformations
#' 
#' @examples 
#' 
#' rev_comp(c("aagtgGGTGaa", "AAGTGGT"))
#' 
#' @export
rev_comp <- function(seqs) {
  # Handle zero length inputs
  if (length(seqs) == 0) {
    return(seqs)
  }
  
  # Capture names to restore later
  seq_names <- names(seqs)
  
  # Split into a list of vectors 
  seq_list <- strsplit(seqs, "")
  
  # Reverse
  seq_list <- lapply(seq_list, rev)
  
  # Get case of input so it can be restored later
  seq_case <- lapply(seq_list, grepl, pattern = "^[[:upper:]]+$")
  
  # Find complement
  seq_list <- lapply(seq_list, seqinr::comp, ambiguous = TRUE)
  
  # Restore case
  seq_list <- lapply(seq_len(length(seq_list)), function(i) {
    ifelse(seq_case[[i]], toupper(seq_list[[i]]), tolower(seq_list[[i]]))
  })
  
  # Paste back together
  result <- vapply(seq_list, FUN = paste0, FUN.VALUE = character(1), collapse = "")
  
  # Restore names
  names(result) <- seq_names
  
  return(result)
}


#' Find complement of sequences
#'
#' Find the complement of one or more sequences stored as a character
#' vector. This is a wrapper for \code{\link[seqinr]{comp}} for character
#' vectors instead of lists of character vectors with one value per letter.
#' IUPAC ambiguity code are handled and the upper/lower case is preserved.
#' 
#' @param seqs A character vector with one element per sequence.
#' 
#' @family sequence transformations
#' 
#' @examples 
#' 
#' complement(c("aagtgGGTGaa", "AAGTGGT"))
#' 
#' @export
complement <- function(seqs) {
  # Handle zero length inputs
  if (length(seqs) == 0) {
    return(seqs)
  }
  
  # Capture names to restore later
  seq_names <- names(seqs)
  
  # Split into a list of vectors 
  seq_list <- strsplit(seqs, "")
  
  # Get case of input so it can be restored later
  seq_case <- lapply(seq_list, grepl, pattern = "^[[:upper:]]+$")
  
  # Find complement
  seq_list <- lapply(seq_list, seqinr::comp, ambiguous = TRUE)
  
  # Restore case
  seq_list <- lapply(seq_len(length(seq_list)), function(i) {
    ifelse(seq_case[[i]], toupper(seq_list[[i]]), tolower(seq_list[[i]]))
  })
  
  # Paste back together
  result <- vapply(seq_list, FUN = paste0, FUN.VALUE = character(1), collapse = "")
  
  # Restore names
  names(result) <- seq_names
  
  return(result)
}


#' Reverse sequences
#'
#' Find the reverse of one or more sequences stored as a character
#' vector. This is a wrapper for \code{\link{rev}} for character
#' vectors instead of lists of character vectors with one value per letter.
#' 
#' @param seqs A character vector with one element per sequence.
#' 
#' @family sequence transformations
#' 
#' @examples 
#' 
#' reverse(c("aagtgGGTGaa", "AAGTGGT"))
#' 
#' @export
reverse <- function(seqs) {
  # Handle zero length inputs
  if (length(seqs) == 0) {
    return(seqs)
  }
  
  # Capture names to restore later
  seq_names <- names(seqs)
  
  # Split into a list of vectors 
  seq_list <- strsplit(seqs, "")
  
  # Reverse
  seq_list <- lapply(seq_list, rev)
  
  # Paste back together
  result <- vapply(seq_list, FUN = paste0, FUN.VALUE = character(1), collapse = "")
  
  # Restore names
  names(result) <- seq_names
  
  return(result)
}


#' Capitalize
#' 
#' Make the first letter uppercase
#' 
#' @param text Some text
#' 
#' @keywords internal
capitalize <- function(text) {
  paste0(toupper(substr(text, 1, 1)), 
         substr(text, vapply(text, FUN.VALUE = numeric(1), function(x) min(c(2, nchar(x)))), nchar(text)))
}
