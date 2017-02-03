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
#' @examples 
#' # Add one or more taxon columns
#' mutate_taxa(unite_ex_data_3, x = 1, y = x+2)
#' 
#' @export
mutate_taxa <- function(.data, ...) {
  data_used <- taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))
  # unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
  #                  names(my_taxon_data), my_taxon_data)
  calculated_cols <- names(.data$taxon_funcs)
  modified_data <- dplyr::mutate(data_used, ...)
  new_cols <- colnames(modified_data)[! colnames(modified_data) %in% colnames(data_used)]
  .data$taxon_data <- dplyr::bind_cols(.data$taxon_data[ , ! colnames(.data$taxon_data) %in% new_cols], modified_data[, new_cols])
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
#' @examples 
#' # Add one or more observation columns
#' mutate_obs(unite_ex_data_3, x = 1, y = x+2)
#'   
#' @export
mutate_obs <- function(.data, ...) {
  data_used <- obs_data(.data, col_subset = obs_data_cols_used(.data, ...))
  # unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
  #                  names(my_obs_data), my_obs_data)
  calculated_cols <- names(.data$obs_data)
  modified_data <- dplyr::mutate(data_used, ...)
  new_cols <- colnames(modified_data)[! colnames(modified_data) %in% colnames(data_used)]
  .data$obs_data <- dplyr::bind_cols(.data$obs_data[ , ! colnames(.data$obs_data) %in% new_cols], modified_data[, new_cols])
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
#' @examples 
#' # Replace all taxon columns with new columns
#' transmute_taxa(unite_ex_data_3, x = 1, y = x+2)
#'   
#' @export
transmute_taxa <- function(.data, ...) {
  data_used <- taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))
  # unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
  #                  names(my_taxon_data), my_taxon_data)
  calculated_cols <- names(.data$taxon_funcs)
  modified_data <- dplyr::mutate(data_used, ...)
  new_cols <- colnames(modified_data)[! colnames(modified_data) %in% colnames(data_used)]
  .data$taxon_data <- dplyr::bind_cols(.data$taxon_data[ , c("taxon_ids", "supertaxon_ids"), drop = FALSE],
                                       modified_data[, new_cols])
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
#' @examples 
#' # Replace all observation columns with new columns
#' transmute_obs(unite_ex_data_3, x = 1, y = x+2)
#'   
#' @export
transmute_obs <- function(.data, ...) {
  data_used <- obs_data(.data, col_subset = obs_data_cols_used(.data, ...))
  # unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
  #                  names(my_obs_data), my_obs_data)
  calculated_cols <- names(.data$obs_data)
  modified_data <- dplyr::mutate(data_used, ...)
  new_cols <- colnames(modified_data)[! colnames(modified_data) %in% colnames(data_used)]
  .data$obs_data <- dplyr::bind_cols(.data$obs_data[ , c("obs_taxon_ids"), drop = FALSE],
                                     modified_data[, new_cols])
  return(.data)
} 