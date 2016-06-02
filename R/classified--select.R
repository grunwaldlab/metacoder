#' Subset columns in a \code{\link{classified}} object
#' 
#' Subsets \code{taxon_data} columns in a \code{\link{classified}} object. Takes and returns a
#' \code{\link{classified}} object. Any column name that appears in \code{taxon_data(.data)} can be
#' used as if it was a vector on its own. See \code{\link[dplyr]{select}} for more information.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more column names to return in the new object. This can be one of three things:
#'   \describe{ \item{\code{expression with unquoted column name}}{The name of a column in
#'   \code{taxon_data} typed as if it was a varaible on its own.} \item{\code{numeric}}{Indexes of
#'   columns in \code{taxon_data}} } To match column names with a character vector, use 
#'   \code{matches("my_col_name")}. To match a logical vector, convert it to a column index using
#'   \code{\link{which}}.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
select_taxa <- function(.data, ...) {
  .data$taxon_data <- dplyr::bind_cols(.data$taxon_data[ , c("taxon_ids", "parent_ids"), drop = FALSE],
                                       dplyr::select(.data$taxon_data, ...))
  return(.data)
}


#' Subset columns in a \code{\link{classified}} object
#' 
#' Subsets \code{item_data} columns in a \code{\link{classified}} object. Takes and returns a
#' \code{\link{classified}} object. Any column name that appears in \code{item_data(.data)} can be
#' used as if it was a vector on its own. See \code{\link[dplyr]{select}} for more information.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more column names to return in the new object. This can be one of three things:
#'   \describe{ \item{\code{expression with unquoted column name}}{The name of a column in
#'   \code{taxon_data} typed as if it was a varaible on its own.} \item{\code{numeric}}{Indexes of
#'   columns in \code{taxon_data}} } To match column names with a character vector, use 
#'   \code{matches("my_col_name")}. To match a logical vector, convert it to a column index using
#'   \code{\link{which}}.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
select_items <- function(.data, ...) {
  .data$item_data <- dplyr::bind_cols(.data$item_data[ , c("item_taxon_ids"), drop = FALSE],
                                      dplyr::select(.data$item_data, ...))
  return(.data)
}