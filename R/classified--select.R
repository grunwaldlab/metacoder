#' Subset columns in a \code{\link{classified}} object
#' 
#' Subsets \code{taxon_data} columns in a \code{\link{classified}} object.
#' Takes and returns a \code{\link{classified}} object.
#' Any column name that appears in \code{taxon_data(.data)} can be used as if it was a vector on its own.
#' @param .data \code{\link{classified}}
#' @param ... One or more column names to return in the new object.
#' This can be one of three things:
#' \describe{
#'    \item{\code{character}}{The name of a column in \code{taxon_data}}
#'    \item{\code{unquoted column name}}{The name of a column in \code{taxon_data} typed as if it was a varaible on its own.}
#'    \item{\code{logical}}{A \code{TRUE}/\code{FALSE} vector of length equal to the number of columns in \code{taxon_data}}
#'    \item{\code{numeric}}{Indexes of columns in \code{taxon_data}}
#'  }
#'  
#' @return An object of type \code{\link{classified}}
#' 
#' @export
select_taxa <- function(.data, ...) {
  # non-standard argument evaluation ---------------------------------------------------------------
  dots <- lazyeval::lazy_dots(...)
  
  process_one <- function(unevaluated) {
    evaluated <- lazyeval::lazy_eval(unevaluated, data = .data$taxon_data)
    if (length(as.character(unevaluated$expr)) == 1 && as.character(unevaluated$expr) %in% colnames(.data$taxon_data)) { # column name
      return(colnames(.data$taxon_data) == as.character(unevaluated$expr))
    } else if (is.logical(evaluated)) { # logical vector
      if (length(evaluated) == ncol(.data$taxon_data)) {
        return(evaluated)
      } else {
        stop("Incorrect number of columns implied by logical vector.")
      }
    } else if (is.numeric(evaluated)) { # index
      if (all(evaluated %in% 1:ncol(.data$taxon_data))) {
        return(1:ncol(.data$taxon_data) == evaluated)
      } else {
        stop("Index supplied does not correspond to a column index.")
      }
    } else {
      stop("Could not parse column selector.")
    }
  }
  
  selected <- Reduce(`|`, c(lapply(dots, process_one),
                            list(c(TRUE, TRUE, rep(FALSE, ncol(.data$taxon_data) - 2)))))
  
  
  # Select columns ---------------------------------------------------------------------------------
  .data$taxon_data <- .data$taxon_data[ , selected, drop = FALSE]
  
  return(.data)
}




#' Subset columns in a \code{\link{classified}} object
#' 
#' Subsets \code{item_data} columns in a \code{\link{classified}} object.
#' Takes and returns a \code{\link{classified}} object.
#' Any column name that appears in \code{item_data(.data)} can be used as if it was a vector on its own.
#' @param .data \code{\link{classified}}
#' @param ... One or more column names to return in the new object.
#' This can be one of three things:
#' \describe{
#'    \item{\code{character}}{The name of a column in \code{item_data}}
#'    \item{\code{unquoted column name}}{The name of a column in \code{item_data} typed as if it was a varaible on its own.}
#'    \item{\code{logical}}{A \code{TRUE}/\code{FALSE} vector of length equal to the number of columns in \code{item_data}}
#'    \item{\code{numeric}}{Indexes of columns in \code{item_data}}
#'  }
#'  
#' @return An object of type \code{\link{classified}}
#' 
#' @export
select_items <- function(.data, ...) {
  # non-standard argument evaluation ---------------------------------------------------------------
  dots <- lazyeval::lazy_dots(...)
  
  process_one <- function(unevaluated) {
    evaluated <- lazyeval::lazy_eval(unevaluated, data = .data$item_data)
    if (length(as.character(unevaluated$expr)) == 1 && as.character(unevaluated$expr) %in% colnames(.data$item_data)) { # column name
      return(colnames(.data$item_data) == as.character(unevaluated$expr))
    } else if (is.logical(evaluated)) { # logical vector
      if (length(evaluated) == ncol(.data$item_data)) {
        return(evaluated)
      } else {
        stop("Incorrect number of columns implied by logical vector.")
      }
    } else if (is.numeric(evaluated)) { # index
      if (all(evaluated %in% 1:ncol(.data$item_data))) {
        return(1:ncol(.data$item_data) == evaluated)
      } else {
        stop("Index supplied does not correspond to a column index.")
      }
    } else {
      stop("Could not parse column selector.")
    }
  }
  
  selected <- Reduce(`|`, c(lapply(dots, process_one),
                            list(c(TRUE, rep(FALSE, ncol(.data$item_data) - 1)))))
  
  
  # Select columns ---------------------------------------------------------------------------------
  .data$item_data <- .data$item_data[ , selected, drop = FALSE]
  
  return(.data)
}