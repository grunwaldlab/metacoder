#' Sort columns of \code{\link{taxmap}} objects
#' 
#' Sort columns of \code{taxon_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{taxon_data(.data)} can be used as if it was a vector on its own. See
#' \link[dplyr]{arrange} for more details.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to sort on. Newly created columns can be referenced in the
#'   same function call.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @export
arrange_taxa <- function(.data, ...) {
  my_taxon_data <- taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))
  unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
                   names(my_taxon_data), my_taxon_data)
  .data$taxon_data <- dplyr::arrange(.data$taxon_data, ...)
  return(.data) 
} 


#' Sort columns of \code{\link{taxmap}} objects
#' 
#' Sort columns of \code{item_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{item_data(.data)} can be used as if it was a vector on its own. See
#' \link[dplyr]{arrange} for more details.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to sort on. Newly created columns can be referenced in the
#'   same function call.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @export
arrange_items <- function(.data, ...) {
  my_item_data <- item_data(.data, col_subset = item_data_cols_used(.data, ...))
  unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
                   names(my_item_data), my_item_data)
  .data$item_data <- dplyr::arrange(.data$item_data, ...)
  return(.data) 
} 
