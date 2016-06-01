#' Format taxon subset value
#' 
#' Format an input to a \code{subset} option on functions like \code{\link{supertaxa}}.
#' Converts logical and numeric vectors into taxon ids
#' 
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param index
#' If a \code{character}, then it should be values of \code{taxon_ids}.
#' If a \code{numeric}, then it should be row indexes of \code{taxon_data}.
#' If a \code{logical}, then it should correspond to rows of \code{taxon_data}.
#' 
#' @return \code{character}
#' 
#' @keywords internal
format_taxon_subset <- function(obj, index) {
  if (is.character(index)) {
    if (any(! index %in% obj$taxon_data$taxon_ids)) {
      stop("All `index` must be in `obj$taxon_data$taxon_ids`.")
    }
  }
  if (is.logical(index) || is.numeric(index)) {
    index <- obj$taxon_data$taxon_ids[index]
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
                                       function(x) gsub(obj$taxon_data[obj$taxon_data$taxon_ids == x[1], name_col],
                                                        pattern = paste0("^", obj$taxon_data[obj$taxon_data$taxon_ids == x[2], name_col], "[_ ]+"),
                                                        replacement = ""),
                                       character(1))
  return(obj)
}


#' Get names of taxon_data in an unevaluated expression
#' 
#' Get names of taxon_data in an unevaluated expression
#' 
#' @param obj a \code{classified} object
#' @param ... unevaluated expression
#' 
#' @return \code{character}
#' 
#' @keywords internal
taxon_data_cols_used <- function(obj, ...) {
  names_used <- unlist(lapply(lazyeval::lazy_dots(...), function(x) as.character(x$expr)))
  names_used[names_used %in% taxon_data_colnames(obj)]
}

#' Get names of item_data in an unevaluated expression
#' 
#' Get names of item_data in an unevaluated expression
#' 
#' @param obj a \code{classified} object
#' @param ... unevaluated expression
#' 
#' @return \code{character}
#' 
#' @keywords internal
item_data_cols_used <- function(obj, ...) {
  names_used <- unlist(lapply(lazyeval::lazy_dots(...), function(x) as.character(x$expr)))
  names_used[names_used %in% item_data_colnames(obj)]
}


#' Get column names of taxon_data
#' 
#' Get column names of taxon_data without calculating columns
#' 
#' @param obj a \code{classified} object
#' 
#' @return \code{character}
#' 
#' @export
taxon_data_colnames <- function(obj) {
  c(colnames(obj$taxon_data), names(obj$taxon_funcs))
}

#' Get column names of item_data
#' 
#' Get column names of item_data without calculating columns
#' 
#' @param obj a \code{classified} object
#' 
#' @return \code{character}
#' 
#' @export
item_data_colnames <- function(obj) {
  c(colnames(obj$item_data), names(obj$item_funcs))
}

