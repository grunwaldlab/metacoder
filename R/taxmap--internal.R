#' Format taxon subset value
#' 
#' Format an input to a \code{subset} option on functions like \code{\link{supertaxa}}.
#' Converts logical and \code{taxon_ids} into indexes of \code{taxon_data}.
#' 
#' @param obj (\code{taxmap})
#' The \code{taxmap} object containing taxon information to be queried.
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
      output <- obj$taxon_data$taxon_ids[index]
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
#' Remove the names of parent taxa in the begining of their children's names in a \code{taxmap} object.
#' This is useful for remove genus names in species binomials.
#' 
#' @param obj a \code{taxmap} object
#' @param name_col (\code{character} of length 1)
#' The name of a column in \code{obj$taxon_data}
#' @param all_supertaxa (\code{logical} of length 1) If \code{TRUE}, check all supertaxa for redundant names instead of just the one immediate supertaxa.
#' 
#' @return \code{character} 
#' 
#' @export
remove_redundant_names <- function(obj, name_col, all_supertaxa = TRUE) {
  my_supertaxa <- supertaxa(obj, recursive = all_supertaxa, include_input = TRUE, simplify = FALSE, index = TRUE, na = FALSE)
  has_parent <- vapply(my_supertaxa, length, numeric(1)) > 1
  
  make_pattern <- function(x) {
    paste0("^", paste0(collapse = "|", paste0(obj$taxon_data[[name_col]][x[-1]], "[_ +-]+")))
  }
  
  obj$taxon_data[has_parent, name_col] <- vapply(my_supertaxa[has_parent], 
                                                 function(x) gsub(obj$taxon_data[[name_col]][x[1]],
                                                                  pattern = make_pattern(x),
                                                                  replacement = ""),
                                                 character(1))
  return(obj)
}


#' Get names of taxon_data in an unevaluated expression
#' 
#' Get names of taxon_data in an unevaluated expression
#' 
#' @param obj a \code{taxmap} object
#' @param ... unevaluated expression
#' 
#' @return \code{character}
#' 
#' @keywords internal
taxon_data_cols_used <- function(obj, ...) {
  decompose <- function(x) {
    if (class(x) %in% c("call", "(")) {
      return(lapply(1:length(x), function(i) decompose(x[[i]])))
    } else {
      return(as.character(x))
    }
  }
  
  expressions <- lapply(lazyeval::lazy_dots(...), function(x) x$expr)
  if (length(expressions) == 0) {
    return(character(0))
  } else {
    names_used <- unlist(lapply(1:length(expressions), function(i) decompose(expressions[[i]])))
    return(unique(names_used[names_used %in% taxon_data_colnames(obj)]))
  }
}

#' Get names of obs_data in an unevaluated expression
#' 
#' Get names of obs_data in an unevaluated expression
#' 
#' @param obj a \code{taxmap} object
#' @param ... unevaluated expression
#' 
#' @return \code{character}
#' 
#' @keywords internal
obs_data_cols_used <- function(obj, ...) {
  decompose <- function(x) {
    if (class(x) %in% c("call", "(")) {
      return(lapply(1:length(x), function(i) decompose(x[[i]])))
    } else {
      return(as.character(x))
    }
  }
  
  expressions <- lapply(lazyeval::lazy_dots(...), function(x) x$expr)
  if (length(expressions) == 0) {
    return(character(0))
  } else {
    names_used <- unlist(lapply(1:length(expressions), function(i) decompose(expressions[[i]])))
    return(unique(names_used[names_used %in% obs_data_colnames(obj)]))
  }
}


#' Get column names of taxon_data
#' 
#' Get column names of taxon_data without calculating columns
#' 
#' @param obj a \code{taxmap} object
#' 
#' @return \code{character}
#' 
#' @export
taxon_data_colnames <- function(obj) {
  c(colnames(obj$taxon_data), names(obj$taxon_funcs))
}

#' Get column names of obs_data
#' 
#' Get column names of obs_data without calculating columns
#' 
#' @param obj a \code{taxmap} object
#' 
#' @return \code{character}
#' 
#' @export
obs_data_colnames <- function(obj) {
  c(colnames(obj$obs_data), names(obj$obs_funcs))
}

