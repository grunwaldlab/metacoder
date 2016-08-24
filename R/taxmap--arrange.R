#' Sort columns of \code{\link{taxmap}} objects
#' 
#' Sort columns of \code{taxon_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{taxon_data(.data)} can be used as if it was a vector on its own. See
#' \link[dplyr]{arrange} for more details.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to sort on.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#' 
#' @examples
#' # Sort by taxon name alphabetically
#' arrange_taxa(unite_ex_data_3, name)
#' # Reverse order of sort
#' arrange_taxa(unite_ex_data_3, desc(name))
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
#' Sort columns of \code{obs_data} in \code{\link{taxmap}} objects. Any column name that
#' appears in \code{obs_data(.data)} can be used as if it was a vector on its own. See
#' \link[dplyr]{arrange} for more details.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more column names to sort on.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @examples
#' # Sort observations by sequence name alphabetically
#' arrange_obs(unite_ex_data_3, seq_name)
#' # Reverse order of sort
#' arrange_obs(unite_ex_data_3, desc(seq_name))
#' 
#' @export
arrange_obs <- function(.data, ...) {
  my_obs_data <- obs_data(.data, col_subset = obs_data_cols_used(.data, ...))
  unused <- mapply(function(name, value) assign(name, value, envir = parent.frame(2)),
                   names(my_obs_data), my_obs_data)
  .data$obs_data <- dplyr::arrange(.data$obs_data, ...)
  return(.data) 
} 
