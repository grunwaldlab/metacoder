#' Add columns to \code{\link{taxmap}} objects
#' 
#' Add columns to the \code{taxon_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{taxon_data(.data)} can be used as if it was a vector on its own. See
#' \code{\link[dplyr]{mutate}} for inspiration and more information.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to add to the new object. Newly created columns can be
#'   referenced in the same function call.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @export
mutate_taxa <- function(.data, ...) {
  my_taxon_data <- taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))
  unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
                   names(my_taxon_data), my_taxon_data)
  .data$taxon_data <- dplyr::mutate(.data$taxon_data, ...)
  return(.data) 
} 



#' Add columns to \code{\link{taxmap}} objects
#' 
#' Add columns to the \code{obs_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{obs_data(.data)} can be used as if it was a vector on its own. See
#' \code{\link[dplyr]{mutate}} for inspiration and more information.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to add to the new object. Newly created columns can be
#'   referenced in the same function call.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @export
mutate_obs <- function(.data, ...) {
  my_obs_data <- obs_data(.data, col_subset = obs_data_cols_used(.data, ...))
  unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
                   names(my_obs_data), my_obs_data)
  .data$obs_data <- dplyr::mutate(.data$obs_data, ...)
  return(.data) 
} 


#' Replace columns in \code{\link{taxmap}} objects
#' 
#' Replace columns of \code{taxon_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{taxon_data(.data)} can be used as if it was a vector on its own. See
#' \code{\link[dplyr]{transmute}} for inspiration and more information.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to add to the new object. Newly created columns can be
#'   referenced in the same function call.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @export
transmute_taxa <- function(.data, ...) {
  my_taxon_data <- taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))
  unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
                   names(my_taxon_data), my_taxon_data)
  .data$taxon_data <- dplyr::bind_cols(.data$taxon_data[ , c("taxon_ids", "parent_ids"), drop = FALSE],
                                       dplyr::transmute(.data$taxon_data, ...))
  return(.data) 
} 


#' Replace columns in \code{\link{taxmap}} objects
#' 
#' Replace columns of \code{obs_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{obs_data(.data)} can be used as if it was a vector on its own. See
#' \code{\link[dplyr]{transmute}} for inspiration and more information.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to add to the new object. Newly created columns can be
#'   referenced in the same function call.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @export
transmute_obs <- function(.data, ...) {
  my_obs_data <- obs_data(.data, col_subset = obs_data_cols_used(.data, ...))
  unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
                   names(my_obs_data), my_obs_data)
  .data$obs_data <- dplyr::bind_cols(.data$obs_data[ , c("obs_taxon_ids"), drop = FALSE],
                                      dplyr::transmute(.data$obs_data, ...))
  return(.data) 
} 