#' Create an instance of \code{taxmap}
#'
#' Create an instance of \code{taxmap} containing observations taxmap by a taxonomy.
#'
#' @param taxon_ids (\code{character}) 
#' These are unique identifiers for taxa.
#' They will be coerced into characters.
#' @param supertaxon_ids (\code{character} OR (\code{numeric}))
#' Supertaxa of \code{taxa}.
#' If a \code{character} vector, then these should be in the same format as \code{taxon_ids}.
#' If a \code{numeric} vector, then it is interpreted as the indexes of \code{taxon_ids}.
#' Taxa without parents should be \code{NA}.
#' @param taxa (\code{character})
#' Objects representing taxa.
#' Currently, these can be anything, but this might change in the future.
#' @param obs_taxon_ids (\code{character} OR (\code{numeric}))
#' Taxon assignments of observations.
#' Supertaxa of \code{taxa}.
#' If a \code{character} vector, then these should be in the same format as \code{taxon_ids}.
#' If a \code{numeric} vector, then it is interpreted as the indexes of \code{taxon_ids}.
#' @param taxon_data (\code{data.frame})
#' A table with rows pretaining to \code{taxa}
#' @param obs_data A (\code{data.frame})
#' A table with rows pretaining to \code{obs_taxa}
#' @param taxon_funcs (\code{list} of named \code{function}s)
#' These the values produced by these functions will be accessible as a column in \code{taxon_data}.
#' The first parameter of each function should be a single \code{taxmap} object.
#' @param obs_funcs (\code{list} of named \code{function}s)
#' These the values produced by these functions will be accessible as a column in \code{obs_data}.
#' The first parameter of each function should be a single \code{taxmap} object.
#' 
#' @return An object of type \code{taxmap}
#'
#' @export
taxmap <- function(taxon_ids, supertaxon_ids,
                       taxa = taxon_ids,
                       obs_taxon_ids = numeric(0),
                       taxon_data = NULL, obs_data = NULL,
                       taxon_funcs = list(n_obs = n_obs,
                                          n_obs_1 = n_obs_1,
                                          n_supertaxa = n_supertaxa,
                                          n_subtaxa = n_subtaxa,
                                          n_subtaxa_1 = n_subtaxa_1,
                                          hierarchies = hierarchies),
                       obs_funcs = list()) {
  # Validate `taxon_ids` ---------------------------------------------------------------------------
  # Coerce into character vector
  taxon_ids <- as.character(taxon_ids)
  # Check that `taxon_ids` is the same length as `taxa`
  if (length(taxa) != length(taxon_ids)) {
    stop("'taxon_ids' must be the same length as 'taxa'")
  }
  # Check that taxa are unique
  if (length(unique(taxon_ids)) != length(taxon_ids)) { stop("'taxon_ids' must be unique") }
  # Name `taxa` with taxon_ids
  taxa <- stats::setNames(taxa, taxon_ids)
  
  # Validate `supertaxon_ids` -----------------------------------------------------------------------------
  # Check that `supertaxon_ids` is the same length as `taxa`
  if (length(taxa) != length(supertaxon_ids)) {
    stop("'supertaxon_ids' must be the same length as 'taxa'")
  }
  # Convert indexes to values of `taxon_ids`
  if (is.numeric(supertaxon_ids)) {
    supertaxon_ids <- taxon_ids[supertaxon_ids]
  } 
  supertaxon_ids[! supertaxon_ids %in% taxon_ids] <- NA
  
  # Validate `obs_taxon_ids` ---------------------------------------------------------------------------
  # Check that all `obs_taxon_ids` are in `taxon_ids`
  if (is.character(obs_taxon_ids)) {
    obs_taxon_ids[! obs_taxon_ids %in% taxon_ids] <- NA
    if (any(is.na(obs_taxon_ids))) {
      warning("Some `obs_taxon_ids` could not be found in `taxon_ids`. They will be `NA`. ")
    }
  } 
  # Convert indexes to values of `taxon_ids`
  if (is.numeric(obs_taxon_ids)) {
    obs_taxon_ids <- taxon_ids[obs_taxon_ids]
  }
  
  # Validate `taxon_data` and `obs_data` ----------------------------------------------------------
  # Check that the tables are structured correctly
  validate_table <- function(data, ids, ...) {
    reserved_col_names = c("taxon_ids", "supertaxon_ids", "obs_taxon_ids")
    result <- dplyr::tbl_df(as.data.frame(list(...), stringsAsFactors = FALSE))
    if (! is.null(data)) {
      data <- as.data.frame(data)
      if (nrow(data) != length(ids)) {
        stop(paste("'", match.call()[["data"]],
                   "' must have the same number of rows as '",
                   match.call()[["id"]], "'"))
      }
      if (any(colnames(data) %in% reserved_col_names)) {
        stop(paste("Column names cannot be one of the following:",
                   paste0(reserved_col_names, collapse = ", ")))
      }
      result <- dplyr::bind_cols(result, data)
    } 
    return(result)
  }
  taxon_data <- validate_table(taxon_data, taxon_ids, taxon_ids = taxon_ids, 
                              supertaxon_ids = supertaxon_ids)
  obs_data <- validate_table(obs_data, obs_taxon_ids, obs_taxon_ids = obs_taxon_ids)
  # Check that tables do not share column names
  all_col_names <- c(colnames(taxon_data), colnames(obs_data))
  if (length(unique(all_col_names)) != length(all_col_names)) {
    stop("'taxon_data' and 'obs_data' can not share column names.")
  }

  # Validate column-generating functions -----------------------------------------------------------
  validate_col_funcs <- function(col_funcs, col_funcs_name) {
    is_named <- function(x) { (! is.null(names(x))) && all(names(x) != '') }
    if ( length(col_funcs) > 1 && (! all(sapply(col_funcs, is.function)) || ! is_named(col_funcs)) ) {
      stop(paste("'", col_funcs_name, "' must all be named functions"))
    }
  }
  validate_col_funcs(taxon_funcs, "taxon_funcs") 
  validate_col_funcs(obs_funcs, "obs_funcs")

  # Make object
  output <- list(taxa = taxa,
                 taxon_data = taxon_data,
                 obs_data = obs_data,
                 taxon_funcs = taxon_funcs,
                 obs_funcs = obs_funcs)
  class(output) <- "taxmap"
  return(output)
}


#' Print a \code{\link{taxmap}} object
#' 
#' Print a \code{\link{taxmap}} object
#' 
#' @param x
#' object to print
#' @param max_rows (\code{integer} of length 1)
#' The maximum number of rows to print in tables. 
#' @param ... Not used
#' 
#' @export
print.taxmap <- function(x, max_rows = 7, ...) {
  loadNamespace("dplyr") # used for print methods
  max_chars <- getOption("width") - 12
  
   print_header <- function(var_name) {
     target_width <- max_chars
     spacer_count <- (target_width - nchar(var_name) - 2) / 2
     spacer <- paste0(rep("-", spacer_count), collapse = "")
     cat(paste0("\n", spacer, " ", var_name, " ", spacer, "\n"))
   }
   
   print_chars <- function(chars) {
     
     interleave <- function(v1,v2) { # https://stat.ethz.ch/pipermail/r-help/2006-March/101023.html
       ord1 <- 2*(1:length(v1))-1
       ord2 <- 2*(1:length(v2))
       c(v1,v2)[order(c(ord1,ord2))]
     }
     
     q = "'"
     interleaved <- interleave(chars[1:(length(chars) / 2)], 
                               rev(chars[(length(chars) / 2 + 1):length(chars)]))
     is_greater_than_max <- cumsum(nchar(interleaved) + 2) + 10 > max_chars
     if (all(! is_greater_than_max)) { 
       max_printed <- length(chars)
     } else {
       max_printed <- which.max(is_greater_than_max)
     }
     if (max_printed < length(chars)) {
       first_part <-  chars[1:as.integer(max_printed / 2 - 0.5)]
       second_part <- chars[as.integer(length(chars) - (max_printed / 2) + 1.5):length(chars)]
       output <- paste0(paste0(collapse = ", ", first_part),
                        " ... ",
                        paste0(collapse = ", ", second_part),
                        "\n")
     } else {
       output <- paste0(paste0(collapse = ", ", chars), "\n")
     }
     cat(output)
   }
  
  
  cat(paste0('`taxmap` object with data for ', nrow(x$taxon_data),
             ' taxa and ', nrow(x$obs_data), ' observations:\n'))
  print_header("taxa")
  print_chars(names(x$taxa))
  print_header("taxon_data")
  print(x$taxon_data, n = max_rows)
  print_header("obs_data")
  print(x$obs_data, n = max_rows)
  if (length(x$taxon_funcs) > 0) {
    print_header("taxon_funcs")
    print_chars(names(x$taxon_funcs))
  }
  if (length(x$obs_funcs) > 0) {
    print_header("obs_funcs")
    print_chars(names(x$obs_funcs))
  }
  invisible(x)
}


#' Return taxon data from \code{\link{taxmap}}
#'
#' Return a table of data associated with taxa of and object of type
#' \code{\link{taxmap}}.
#'
#' @param obj (\code{\link{taxmap}})
#' @param row_subset (\code{character})
#' The taxon_ids of a subset of \code{obj}.
#' Default: All rows.
#' @param col_subset (\code{character})
#' The names of columns, either user defined or generated using \code{taxon_funcs}.
#' Default: All columns.
#' @param calculated_cols (\code{logical} of length 1)
#' If \code{TRUE}, return calculated columns using  functions in \code{\link{taxmap}$taxon_funcs}.
#' These values are calculated each time \code{taxon_data} is called since their values can change if 
#' the data is subset.
#' @param sort_by (\code{character} of length 1)
#' The name of a column in \code{obj$taxon_data} or a function name in  \code{obj$taxon_funcs}.
#' This column will be used to sort the output rows.
#' If \code{NULL}, no sorting will be done.
#' @param decreasing (\code{logical} of length 1)
#' If \code{TRUE}, \code{sort_by} order is decreasing.
#' @param drop (\code{logical} of length 1)
#' If \code{TRUE}, if \code{subset} is a single column
#' name, then a \code{vector} is returned instead of a \code{data.frame}
#'
#' @return A \code{data.frame} or \code{vector} with rows corresponding to taxa in input
#'
#' @export
taxon_data <- function(obj,
                       row_subset = NULL,
                       col_subset = NULL,
                       calculated_cols = TRUE,
                       sort_by = NULL,
                       decreasing = FALSE,
                       drop = FALSE) {
  # Parse options
  if (is.null(row_subset)) {
    row_subset <- 1:nrow(obj$taxon_data)
  } else {
    row_subset <- format_taxon_subset(obj, row_subset)
  }
  if (is.null(col_subset)) {
    col_subset <- c(colnames(obj$taxon_data), names(obj$taxon_funcs))
  }
  
  # Check that the user is making sense
  if (calculated_cols == FALSE && any(col_subset %in% names(obj$taxon_funcs))) {
    stop("Cannot use a calculated column when `calculated_cols = FALSE`.")
  }
  # Make copy of taxon data
  data <- obj$taxon_data
  # Remove any user-defined rows not specified
  data <- dplyr::slice(data, row_subset)
  # Check if any of the column-generating functions are needed
  functions <- obj$taxon_funcs[names(obj$taxon_funcs) %in% col_subset]
  # Apply column-generating functions and append to output
  if (calculated_cols && length(functions) > 0) {
    calculated_data <- lapply(functions, function(f) f(obj, row_subset))
    names(calculated_data) <- names(functions)
    data <- dplyr::bind_cols(data, calculated_data)
  }
  # Remove any user-defined columns not specified
  data <- data[, colnames(data) %in% col_subset, drop = FALSE]
  # Reorder output columns to match order of col_subset
  data <- data[ , col_subset, drop = FALSE]
  # Reorder output rows according to `sort_by`
  if (! is.null(sort_by)) {
    if (is.character(sort_by) && sort_by %in% colnames(data)) {
      sort_by_col <- data[ , sort_by][[1]]
    } else if (is.function(sort_by) && all(c("obj", "subset") %in% names(formals(sort_by)))) {
      sort_by_col <- sort_by(obj, row_subset)
    } else if (sort_by %in% names(obj$taxon_funcs)) {
      sort_by_col <- obj$taxon_funcs[[sort_by]](obj, row_subset)
    } else {
      stop("Could not identify `sort_by` value. It must be a column displayed by `taxon_data`, a function to has `obj` and `subset` arguments, or a function in `obj$taxon_funcs`.")
    }
    data <- data[order(sort_by_col, decreasing = decreasing), ]
  }
  
  # Apply drop
  if (ncol(data) == 1 && drop) {
    data <- data[[1]]
  }
  
  return(data)
}


#' Return observation data from \code{\link{taxmap}}
#'
#' Return a table of data associated with taxa of and object of type
#' \code{\link{taxmap}}.
#'
#' @param obj (\code{\link{taxmap}})
#' @param row_subset (\code{character})
#' The obs_ids of a subset of \code{obj}.
#' Default: All rows.
#' @param col_subset (\code{character})
#' The names of columns, either user defined or generated using \code{obs_funcs}.
#' Default: All columns.
#' @param calculated_cols (\code{logical} of length 1)
#' If \code{TRUE}, return calculated columns using  functions in \code{\link{taxmap}$obs_funcs}.
#' These values are calculated each time \code{obs_data} is called since their values can change if 
#' the data is subset.
#' @param sort_by (\code{character} of length 1)
#' The name of a column in \code{obj$obs_data} or a function name in  \code{obj$obs_funcs}.
#' This column will be used to sort the output rows.
#' If \code{NULL}, no sorting will be done.
#' @param decreasing (\code{logical} of length 1)
#' If \code{TRUE}, \code{sort_by} order is decreasing.
#' @param drop (\code{logical} of length 1)
#' If \code{TRUE}, if \code{subset} is a single column
#' name, then a \code{vector} is returned instead of a \code{data.frame}
#'
#' @return A \code{data.frame} or \code{vector} with rows corresponding to taxa in input
#'
#' @export
obs_data <- function(obj,
                      row_subset = NULL,
                      col_subset = NULL,
                      calculated_cols = TRUE,
                      sort_by = NULL,
                      decreasing = FALSE,
                      drop = FALSE) {
  # Parse options
  if (is.null(row_subset)) {
    row_subset <- 1:nrow(obj$obs_data)
  }
  if (is.null(col_subset)) {
    col_subset <- c(colnames(obj$obs_data), names(obj$obs_funcs))
  }
  
  # Check that the user is making sense
  if (calculated_cols == FALSE && any(col_subset %in% names(obj$obs_funcs))) {
    stop("Cannot use a calculated column when `calculated_cols = FALSE`.")
  }
  # Make copy of observation data
  data <- obj$obs_data
  # Remove any user-defined rows not specified
  data <- data[row_subset, , drop = FALSE]
  # Check if any of the column-generating functions are needed
  functions <- obj$obs_funcs[names(obj$obs_funcs) %in% col_subset]
  # Apply column-generating functions and append to output
  if (calculated_cols && length(functions) > 0) {
    calculated_data <- lapply(functions, function(f) f(obj, row_subset))
    names(calculated_data) <- names(functions)
    data <- cbind(data, as.data.frame(calculated_data))
  }
  # Remove any user-defined columns not specified
  data <- data[, colnames(data) %in% col_subset, drop = FALSE]
  # Reorder output columns to match order of col_subset
  data <- data[ , col_subset, drop = FALSE]
  # Reorder output rows according to `sort_by`
  if (! is.null(sort_by)) {
    if (is.character(sort_by) && sort_by %in% colnames(data)) {
      sort_by_col <- data[ , sort_by]
    } else if (is.function(sort_by) && all(c("obj", "subset") %in% names(formals(sort_by)))) {
      sort_by_col <- sort_by(obj, row_subset)
    } else if (sort_by %in% names(obj$obs_funcs)) {
      sort_by_col <- obj$obs_funcs[[sort_by]](obj, row_subset)
    } else {
      stop("Could not identify `sort_by` value. It must be a column displayed by `obs_data`, a function to has `obj` and `subset` arguments, or a function in `obj$obs_funcs`.")
    }
    data <- data[order(sort_by_col, decreasing = decreasing), ]
  }
  
  # Apply drop
  if (ncol(data) == 1 && drop) {
    data <- data[[1]]
  } else {
    data <- dplyr::tbl_df(data)
  }
  
  return(data)
}
