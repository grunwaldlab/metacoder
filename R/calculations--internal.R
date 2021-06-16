#' Get numeric columns from taxmap table
#' 
#' If columns are specified by the user, parse them and check that they are numeric.
#' If not, return all numeric columns.
#' 
#' @param obj A taxmap object
#' @param data The name of a table in \code{obj}.
#' @param cols The names/indexes of columns in \code{data} to use. By
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
get_numeric_cols <- function(obj, data, cols = NULL) {
  # Get input table
  input <- get_taxmap_table(obj, data)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(input, is.numeric, logical(1)))
    if (length(cols) > 0) {
      my_print("No `cols` specified, so using all numeric columns:\n  ", 
               limited_print(names(cols), type = "silent"))
    } else {
      my_print("No `cols` specified and no numeric columns can be found.")
    }
  }
  
  # Parse user input for columns
  cols <- get_taxmap_cols(obj = obj, data = data, cols = cols)
  
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
#' For a given table in a taxmap object, run some function to produce new
#' columns. This function handles all of the option parsing and formatting of
#' the result.
#'
#' @param obj A \code{\link{taxmap}} object
#' @param data The name of a table in \code{obj$data}.
#' @param func The function to apply. Should have the following form:
#'   \code{function(count_table, cols = cols, groups = groups)} and return a table.
#' @param cols The columns in \code{data} to use. By
#'   default, all numeric columns are used. Takes one of the following inputs:
#'   \describe{
#'   \item{TRUE/FALSE:}{All/No columns will used.}
#'   \item{Character vector:}{The names of columns to use} \item{Numeric vector:}{The indexes of
#'   columns to use}
#'   \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Use the columns corresponding to \code{TRUE} values.} }
#' @param groups Group multiple columns per treatment/group. This should be a
#'   vector of group IDs (e.g. character, integer) the same length as
#'   \code{cols} that defines which samples go in which group. When used, there
#'   will be one column in the output for each unique value in \code{groups}.
#' @param other_cols Preserve in the output non-target columns present in the
#'   input data. New columns will always be on the end. The "taxon_id" column
#'   will be preserved in the front. Takes one of the following inputs:
#'   \describe{
#'   \item{NULL:}{No columns will be added back, not even the taxon id column.}
#'   \item{TRUE/FALSE:}{All/None of the non-target columns will be preserved.}
#'   \item{Character vector:}{The names of columns to preserve}
#'   \item{Numeric vector:}{The indexes of columns to preserve}
#'   \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Preserve the columns corresponding to \code{TRUE} values.}}
#' @param out_names The names of count columns in the output. Must be the same
#'   length and order as \code{cols} (or \code{unique(groups)}, if \code{groups} is used).
#'
#' @return A tibble
#'
#' @keywords internal
do_calc_on_num_cols <- function(obj, data, func, cols = NULL, groups = NULL,
                                other_cols = FALSE, out_names = NULL) {
  
  # Warn if groups is used with no cols specified
  if (is.null(cols) && !is.null(groups)) {
    message('NOTE: Using the "groups" option without the "cols" option can yeild incorrect results if the column order is different from the group order.\n')
  }
  
  # Get input table
  input <- get_taxmap_table(obj, data)
  
  # Parse columns to use
  cols <- get_numeric_cols(obj, data, cols)
  
  # Check that cols, groups and output names make sense
  groups <- check_option_groups(groups, cols)
  if (is.null(groups)) {
    if (is.null(out_names)) { # groups and out_names are NULL
      out_names <- colnames(input[cols])
    } else { # groups is NULL, but out_names set
      if (length(out_names) != length(cols)) {
        stop(call. = FALSE,
             "The length of `cols` (", length(cols),
             ") and `out_names` (", length(out_names),
             ") are not equal.")
      }
    }
    groups <- seq_along(cols)
  } else {
    if (length(groups) != length(cols)) {
      stop(call. = FALSE,
           "`groups` (", length(groups),
           ") must be the same length as `cols` (", length(cols), ").")
    }
    if (is.null(out_names)) { # groups is set, but out_names is NULL
      if (is.numeric(groups)) {
        warning(call. = FALSE,
                "Numeric groups used without supplying 'out_names'. This will result in numeric column names.")
      }
      out_names <- unique(groups)
    } else { # groups and out_names are both set
      if (length(out_names) != length(unique(groups))) {
        stop(call. = FALSE,
             "The length of 'unique(groups)' and 'out_names' are not equal")       
      }
    }
  }
  
  # Check that out_names is a character
  if (! is.null(out_names) && is.numeric(out_names)){
    warning(call. = FALSE,
            "  `out_names` is numeric. This will result in numeric column names.")
  }
  
  # Do calculaton
  if (length(cols) < 1) {
    warning(call. = FALSE,
            "  No cols specified. No calculation will be done.")
    result <- NULL
  } else {
    result <- func(input[, cols], cols = cols, groups = groups)
    result <- result[, as.character(unique(groups)), drop = FALSE] 
    colnames(result) <- out_names
  }
  
  # Add back other columns if specified
  if (! is.null(other_cols)) {
    cols_to_keep <- get_taxmap_other_cols(obj, data, cols, other_cols)
    if (is.null(result)) {
      result <- input[, cols_to_keep]
    } else {
      result <- cbind(input[, cols_to_keep], result)
    }
  }
  
  # Convert to tibble
  result <- dplyr::as_tibble(result)
  
  return(result)
}


#' Check option: groups
#' 
#' This option is used in a few of the calculation functions
#' 
#' @param groups The groups option to check
#' @param cols The cols option, if applicable
#'
#' @keywords internal
check_option_groups <- function(groups, cols = NULL) {
  # Do not do checks if NULL
  if (is.null(groups)) {
    return(groups)
  }
  
  # Check that groups and cols are the same length
  if (! is.null(cols)) {
    if (length(groups) == 1) {
      groups <- rep(groups, length(cols))
    } else if (length(groups) != length(cols)) {
      stop(call. = FALSE,
           "`groups` (", length(groups),
           ") must be length 1 or the same length as `cols` (", length(cols), ").")
    }
  }
  
  # Check for odd values
  if (any(is.na(groups))) {
    warning(call. = FALSE, 
            paste0("NA's detected in `groups` option. This might cause problems. Indexes with NA:\n  ",
                   limited_print(which(is.na(groups)), type = "silent")))
  }
  if (any(groups == "")) {
    warning(call. = FALSE, 
            paste0("Empty values ('') detected in `groups` option. This might cause problems. Indexes with empty values:\n  ",
                   limited_print(which(groups == ""), type = "silent")))
  }
  if (any(is.nan(groups))) {
    warning(call. = FALSE, 
            paste0("NaN's detected in `groups` option. This might cause problems. Indexes with NaN:\n  ",
                   limited_print(which(is.nan(groups)), type = "silent")))
  }
  if (any(is.infinite(groups))) {
    warning(call. = FALSE, 
            paste0("Infinite values detected in `groups` option. This might cause problems. Indexes with infinite values:\n  ",
                   limited_print(which(is.infinite(groups)), type = "silent")))
  }
  
  
  return(groups)
}
