#' Sort columns of \code{\link{classified}} objects
#' 
#' Sort columns of \code{taxon_data} in \code{\link{classified}} objects.
#' Any column name that appears in \code{taxon_data(.data)} can be used as if it was a vector on its own.
#' See \link[dplyr]{arrange} for more details.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more column names to sort on.
#' Newly created columns can be referenced in the same function call.
#' 
#' @return An object of type \code{\link{classified}}
#' 
#' @export
arrange_taxa <- function(.data, ...) {
  my_taxon_data <- taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))
  dots <- lazyeval::lazy_eval(lazyeval::lazy_dots(...), data = my_taxon_data)
  .data$taxon_data <- dplyr::arrange_(.data$taxon_data,
                                      .dots = lapply(dots, function(x) ~ x))
  return(.data) 
} 


#' Sort columns of \code{\link{classified}} objects
#' 
#' Sort columns of \code{item_data} in \code{\link{classified}} objects.
#' Any column name that appears in \code{item_data(.data)} can be used as if it was a vector on its own.
#' See \link[dplyr]{arrange} for more details.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more column names to sort on.
#' Newly created columns can be referenced in the same function call.
#' 
#' @return An object of type \code{\link{classified}}
#' 
#' @export
arrange_items <- function(.data, ...) {
  my_item_data <- item_data(.data, col_subset = item_data_cols_used(.data, ...))
  dots <- lazyeval::lazy_eval(lazyeval::lazy_dots(...), data = my_item_data)
  .data$item_data <- dplyr::arrange_(.data$item_data,
                                     .dots = lapply(dots, function(x) ~ x))
  return(.data) 
} 
