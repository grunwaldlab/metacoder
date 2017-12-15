#' Get numeric columns from taxmap table
#' 
#' If columns are specified by the user, parse them and check that they are numeric.
#' If not, return all numeric columns.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj}.
#' @param cols The names/indexes of columns in \code{dataset} to use. By
#'   default, all numeric columns are used. Takes one of the following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All/No columns will used.}
#'     \item{Character vector:}{The names of columns to use}
#'     \item{Numeric vector:}{The indexes of columns to use}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Use the columns
#'   corresponding to \code{TRUE} values.}
#'   }
#'   
#' @keywords internal
get_numeric_cols <- function(obj, dataset, cols = NULL) {
  # Get input table
  input <- get_taxmap_table(obj, dataset)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(input, is.numeric, logical(1)))
    my_print("No `cols` specified, so using all numeric columns:\n  ", 
             limited_print(names(cols), type = "silent"))
  }
  
  # Parse user input for columns
  cols <- get_taxmap_cols(obj = obj, dataset = dataset, cols = cols)
  
  # Check that count columns are numeric
  col_is_num <- vapply(input[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:\n  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  return(cols)
}


#' Run some function to produce new columns.
#'
#' For a given table in a taxmap object, run some function to produce new columns.
#' This function handles all of the option parsing and formatting of the result.
#'
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj}.
#' @param func The function to apply. Should accept and return a table.
#' @param cols The names/indexes of columns in \code{dataset} to use. By
#'   default, all numeric columns are used. Takes one of the following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All/No columns will used.}
#'     \item{Character vector:}{The names of columns to use}
#'     \item{Numeric vector:}{The indexes of columns to use}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Use the columns
#'   corresponding to \code{TRUE} values.}
#'   }
#' @param other_cols Preserve in the output non-target columns present in the
#'   input data. New columns will always be on the end. The
#'   "taxon_id" column will always be preserved in the front. Takes one of the
#'   following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All non-target columns will be preserved or not.}
#'     \item{Character vector:}{The names of columns to preserve}
#'     \item{Numeric vector:}{The indexes of columns to preserve}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Preserve the columns
#'   corresponding to \code{TRUE} values.}
#'   }
#' @param out_names If supplied, rename the output proportion columns. Must be
#'   the same length as \code{cold}.
#'
#' @return A tibble
#' 
#' @keywords internal
do_calc_on_num_cols <- function(obj, dataset, func, cols = NULL,
                                other_cols = FALSE, out_names = NULL) {
  # Get input table
  input <- get_taxmap_table(obj, dataset)
  
  # Parse columns to use
  cols <- get_numeric_cols(obj, dataset, cols)
  
  # Check that new column names are the same length as calculation columns
  if (is.null(out_names)) {
    out_names <- cols
  } else if (length(out_names) != length(cols)) {
    stop(paste0('The `out_names` option (length = ', length(out_names), 
                ') must be the same length as the `cols` used (length = ',
                length(cols), ').'))
  }
  
  # Find other columns
  #   These might be added back to the output later
  cols_to_keep <- get_taxmap_other_cols(obj, dataset, cols, other_cols)
  
  # Calculate proportions
  result <- func(input[, cols])
  colnames(result) <- out_names
  
  # Add back other columns if specified
  result <- cbind(input[, cols_to_keep], result)
  
  # Convert to tibble
  result <- dplyr::as_tibble(result)
  
  return(result)
}
