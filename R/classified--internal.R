#' Format taxon subset value
#' 
#' Format an input to a \code{subset} option on functions like \code{\link{supertaxa}}.
#' Converts logical and \code{taxon_ids} into indexes of \code{taxon_data}.
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
  if (is.null(index)) {
    output <- stats::setNames(1:nrow(obj$taxon_data), obj$taxon_data$taxon_ids)
  } else {
    if (is.null(names(index))) {
      names(index) <- index
    } 
    if (is.numeric(index)) {
      output <- index
      my_names <- names(index)
    } else if (is.character(index)) {
      output <- match(index, obj$taxon_data$taxon_ids)
      my_names <- names(index)
    } else if (is.logical(index)) {
      output <- which(index)
      my_names <- names(index)[output]
    } else {
      stop("Invalid subset value.")
    }
    names(output) <- my_names
  }
  return(output)
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
  my_supertaxa <- supertaxa(obj, recursive = FALSE, include_input = TRUE)
  has_parent <- vapply(my_supertaxa, length, numeric(1)) > 1
  obj$taxon_data[has_parent, name_col] <- vapply(my_supertaxa[has_parent], 
                                       function(x) gsub(obj$taxon_data[[name_col]][obj$taxon_data$taxon_ids == x[1]],
                                                        pattern = paste0("^", obj$taxon_data[[name_col]][obj$taxon_data$taxon_ids == x[2]], "[_ ]+"),
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

