#' Format taxon subset value
#' 
#' Format an input to a \code{subset} option on functions like \code{subset}.
#' Converts logical and numeric vectors into taxon ids
#' 
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param index
#' 
#' @return \code{character}
#' 
#' @keywords internal
format_taxon_subset <- function(obj, index) {
  if (is.character(index)) {
    index <- vapply(index, FUN.VALUE = numeric(1),
                    function(x) {
                      result = which(x == obj$taxon_data$taxon_ids)
                      if (length(result) == 0) {
                        result = as.numeric(NA)
                      }
                      return(result)
                    })
  }
  if (is.logical(index)) {
    index <- which(index)
  }
  index <- unname(index)
  return(index)
}


#' Remove the redundant taxon names
#' 
#' Remove the names of parent taxa in the begining of their children's names in a \code{classified} object.
#' This is useful for remove genus names in species binomials.
#' 
#' @param obj a \code{classified} object
#' @param name_col (\code{character} of length 1)
#' The name of a column in \code{obj$taxon_data}
#' 
#' @return \code{character} 
#' 
#' @export
remove_redundant_names <- function(obj, name_col) {
  obj$taxon_data[, name_col] <- vapply(supertaxa(obj, recursive = FALSE, include_input = TRUE), 
                                       function(x) gsub(obj$taxon_data[x[1], name_col],
                                                        pattern = paste0("^", obj$taxon_data[x[2], name_col], "[_ ]+"),
                                                        replacement = ""),
                                       character(1))
  return(obj)
}
