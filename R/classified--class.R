#' Create an instance of \code{classified}
#'
#' Create an instance of \code{classified} containing items classified by a taxonomy
#'
#' @param taxa (\code{character})
#' Unique taxa.
#' Currently, these should be respresented by a \code{character} vector, but this might 
#' change in the future.  
#' @param parents (\code{character} OR (\code{numeric}))
#' Parent taxa (i.e. supertaxa) of \code{taxa}.
#' If a \code{character} vector, then these should be in the same format as \code{taxa}.
#' If a \code{numeric} vector, then it is interpreted as the indexes of \code{taxa}.
#' Taxa without parents should be \code{NA}.
#' @param item_taxa (\code{character} OR (\code{numeric}))
#' Taxon assignments of items.
#' Parent taxa (i.e. supertaxa) of \code{taxa}.
#' If a \code{character} vector, then these should be in the same format as \code{taxa}.
#' If a \code{numeric} vector, then it is interpreted as the indexes of \code{taxa}.
#' @param taxon_data (\code{data.frame})
#' A table with rows pretaining to \code{taxa}
#' @param item_data A (\code{data.frame})
#' A table with rows pretaining to \code{item_taxa}
#' @param taxon_funcs (\code{list} of named \code{function}s)
#' These the values produced by these functions will be accessible as a column in \code{taxon_data}.
#' The first parameter of each function should be a single \code{classified} object.
#' @param item_funcs (\code{list} of named \code{function}s)
#' These the values produced by these functions will be accessible as a column in \code{item_data}.
#' The first parameter of each function should be a single \code{classified} object.
#' 
#' @return An object of type \code{classified}
#'
#' @export
classified <- function(taxa, parents, item_taxa = numeric(0),
                       taxon_data = NULL, item_data = NULL,
                       taxon_funcs = list(item_counts = item_counts,
                                          taxon_ranks = taxon_ranks,
                                          classifications = classifications),
                       item_funcs = list()) {
  # Validate `taxa` --------------------------------------------------------------------------------
  # Coerce into a character vector
  taxa <- as.character(taxa)
  # Check that taxa are unique
  if (length(unique(taxa)) != length(taxa)) { stop("'taxa' must be unique") }
  # Make `taxon_ids` vector 
  taxon_ids <- seq_along(taxa)
  
  # Validate `parents` -----------------------------------------------------------------------------
  # Check that `parents` is the same length of `taxa`
  if (length(taxa) != length(parents)) {
    stop("'parents' must be the same length as 'taxa'")
  }
  # Make `parent_ids` vector
  if (is.character(parents)) {
    parent_ids <- vapply(parents, FUN.VALUE = numeric(1),
                         function(x) {
                           result = which(x == taxa)
                           if (length(result) == 0) {
                             result = as.numeric(NA)
                           }
                           return(result)
                         })
    
  } else if (is.numeric(parents)) {
    parent_ids <- parents
    parent_ids[! parent_ids %in% taxon_ids] <- NA
  } else {
    stop("'parents' is invalid.")
  }
  
  # Validate `item_taxa` ---------------------------------------------------------------------------
  if (is.character(item_taxa)) {
    if (any(! item_taxa %in% c(taxa, NA))) {
      stop("All 'item_taxa' must be in 'taxa'")
    }
    item_taxon_ids <- vapply(item_taxa, FUN.VALUE = numeric(1),
                         function(x) {
                           result = which(x == taxa)
                           if (length(result) == 0) {
                             result = as.numeric(NA)
                           }
                           return(result)
                         })
    
  } else if (is.numeric(item_taxa)) {
    if (any(! item_taxa %in% c(taxon_ids, NA))) {
      stop("All 'item_taxa' must be in 'taxa'")
    }
    item_taxon_ids <- item_taxa
  } else {
    stop("'item_taxa' is invalid.")
  }

  # Validate `taxon_data` and `item_data` ----------------------------------------------------------
  # Check that the tables are structured correctly
  validate_data <- function(data, data_var_name, ids, ids_var_name,
                            reserved_col_names = c("taxon_ids", "parent_ids", "item_taxon_ids")) {
    if (! is.null(data)) {
      if (! "data.frame" %in% class(data)) {
        stop(paste("'", data_var_name, "' must be convertable to a data.frame"))
      }
      if (nrow(data) != length(ids)) {
        stop(paste("'", data_var_name, "' must have the same number of rows as '", ids_var_name, "'"))
      }
      if (any(colnames(data) %in% reserved_col_names)) {
        stop(paste("Column names cannot be one of the following:",
                   paste0(reserved_col_names, collapse = ", ")))
      }
    } else {
      data <- data.frame(matrix(numeric(0), nrow = length(ids)))
    }
    return(data)
  }
  taxon_data <- validate_data(taxon_data, "taxon_data", taxon_ids, "taxon_ids")
  item_data <- validate_data(item_data, "item_data", item_taxon_ids, "item_taxa")
  # Check that tables do not share column names
  all_col_names <- c(colnames(taxon_data), colnames(item_data))
  if (length(unique(all_col_names)) != length(all_col_names)) {
    stop("'taxon_data' and 'item_data' can not share column names.")
  }
  # Add ID columns
  taxon_data <- dplyr::tbl_dt(cbind(data.frame(taxon_ids = taxon_ids, 
                                               parent_ids = parent_ids,
                                               row.names = NULL), 
                                    taxon_data))
  data.table::setkey(taxon_data, taxon_ids)
  item_data <- dplyr::tbl_dt(cbind(data.frame(item_taxon_ids = item_taxon_ids,
                                              row.names = NULL),
                                   item_data))

  # Validate column-generating functions -----------------------------------------------------------
  validate_col_funcs <- function(col_funcs, col_funcs_name) {
    is_named <- function(x) { (! is.null(names(x))) && all(names(x) != '') }
    if ( length(col_funcs) > 1 && (! all(sapply(col_funcs, is.function)) || ! is_named(col_funcs)) ) {
      stop(paste("'", col_funcs_name, "' must all be named functions"))
    }
  }
  validate_col_funcs(taxon_funcs, "taxon_funcs") 
  validate_col_funcs(item_funcs, "item_funcs")

  # Make object
  output <- list(taxa = taxa,
                 taxon_data = taxon_data,
                 item_data = item_data,
                 taxon_funcs = taxon_funcs,
                 item_funcs = item_funcs)
  class(output) <- "classified"
  return(output)
}


#' Print a \code{\link{classified}} object
#' 
#' Print a \code{\link{classified}} object
#' 
#' @param x object to print
#' @param ... Passed to other methods
#' 
#' @export
print.classified <- function(x, ...) {
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
     is_greater_than_max <- cumsum(nchar(interleaved) + 4) + 10 > max_chars
     if (all(! is_greater_than_max)) { 
       max_printed <- length(chars)
     } else {
       max_printed <- which.max(is_greater_than_max)
     }
     if (max_printed < length(chars)) {
       first_part <-  chars[1:(max_printed / 2)]
       second_part <- chars[(length(chars) - (max_printed / 2)):length(chars)]
       output <- paste0(q, paste0(collapse = paste0(q, ", ", q), first_part), q,
                        " ... ",
                        q, paste0(collapse = paste0(q, ", ", q), second_part), q,
                        "\n")
     } else {
       output <- paste0(q, paste0(collapse = paste0(q, ", ", q), chars), q, "\n")
     }
     cat(output)
   }
  
  
  cat(paste0('`classified` object with data for ', nrow(x$taxon_data),
             ' taxa and ', nrow(x$item_data), ' items/observations:\n'))
  print_header("taxa")
  print_chars(x$taxa)
  print_header("taxon_data")
  dplyr:::print.tbl_dt(x$taxon_data)
  print_header("item_data")
  dplyr:::print.tbl_dt(x$item_data)
  if (length(x$taxon_funcs) > 0) {
    print_header("taxon_funcs")
    print_chars(names(x$taxon_funcs))
  }
  if (length(x$item_funcs) > 0) {
    print_header("item_funcs")
    print_chars(names(x$item_funcs))
  }
  invisible(x)
}


#' Return taxon data column names
#'
#' Return taxon data column names of and object of type \code{\link{classified}}.
#' This includes "taxon_ids", "parent_ids", user-defined columns, and columns generated by \code{taxon_funcs}.
#'
#' @param obj (\code{\link{classified}})
#'
#' @export
taxon_data_colnames <- function(obj) {
  c("taxon_ids", "parent_ids", colnames(obj$taxon_data),  names(obj$taxon_funcs))
}


#' Return item data column names
#'
#' Return item data column names of and object of type \code{\link{classified}}.
#' This includes "taxon_ids", user-defined columns, and columns generated by \code{item_funcs}.
#'
#' @param obj (\code{\link{classified}})
#'
#' @export
item_data_colnames <- function(obj) {
  c("taxon_ids", colnames(obj$item_data), names(obj$item_funcs))
}


#' Return taxon data from \code{\link{classified}}
#'
#' Return a table of data associated with taxa of and object of type
#' \code{\link{classified}}.
#'
#' @param obj (\code{\link{classified}})
#' @param col_subset (\code{character})
#' The names of columns, either user defined or generated using \code{taxon_funcs}.
#' Default: All columns.
#' @param row_subset (\code{character})
#' The taxon_ids of a subset of \code{obj}.
#' Default: All rows.
#' @param calculated_cols (\code{logical} of length 1)
#' If \code{TRUE}, return calculated columns using  functions in \code{\link{classified}$taxon_funcs}.
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
                       col_subset = taxon_data_colnames(obj),
                       row_subset = taxon_ids(obj),
                       calculated_cols = TRUE,
                       sort_by = NULL,
                       decreasing = FALSE,
                       drop = FALSE) {
  row_subset <- format_taxon_subset(obj, row_subset)
  # Check that the user is making sense
  if (calculated_cols == FALSE && any(col_subset %in% names(obj$taxon_funcs))) {
    stop("Cannot use a calculated column when `calculated_cols = FALSE`.")
  }
  # Combine taxon id information and arbitrary user-defined data
  data <- cbind(data.frame(taxon_ids = obj$taxon_ids, parent_ids = obj$parent_ids, stringsAsFactors = FALSE),
                obj$taxon_data)
  # Remove any user-defined rows/columns not specified
  data <- data[row_subset, colnames(data) %in% col_subset, drop = FALSE]
  # Check if any of the column-generating functions are needed
  functions <- obj$taxon_funcs[names(obj$taxon_funcs) %in% col_subset]
  # Apply column-generating functions and append to output
  if (calculated_cols && length(functions) > 0) {
    calculated_data <- lapply(functions, function(f) f(obj, row_subset))
    names(calculated_data) <- names(functions)
    data <- cbind(data, as.data.frame(calculated_data))
  }
  
  # Reorder output columns to match order of col_subset
  data <- dplyr::tbl_df(data)[ , col_subset]
  
  # Reorder output rows according to `sort_by`
  if (! is.null(sort_by)) {
    if (is.character(sort_by) && sort_by %in% colnames(data)) {
      sort_by_col <- data[ , sort_by]
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
    data = data[[1]]
  }
  
  return(data)
}


#' Return item data from \code{\link{classified}}
#'
#' Return a table of data associated with items of an object of type
#' \code{\link{classified}}.
#'
#' @param obj (\code{\link{classified}})
#' @param col_subset (\code{character})
#' The names of columns, either user defined or generated using \code{taxon_funcs}.
#' Default: All columns.
#' @param row_subset (\code{character})
#' The taxon_ids of a subset of \code{obj}.
#' Default: All rows.
#' @param drop (\code{logical} of length 1) If \code{TRUE}, if \code{subset} is a single column
#' name, then a \code{vector} is returned instead of a \code{data.frame}
#'
#' @return \code{data.frame} with rows corresponding to items in input
#'
#' @export
item_data <- function(obj,
                      col_subset = item_data_colnames(obj),
                      row_subset = 1:nrow(obj$item_data),
                      drop = FALSE) {
  # Check that the user is making sense
#   if (calculated_cols == FALSE && any(col_subset %in% names(obj$item_funcs))) {
#     stop("Cannot use a calculated column when `calculated_cols = FALSE`.")
#   }
  # Combine taxon id information and arbitrary user-defined data
  data <- cbind(data.frame(taxon_ids = obj$item_taxon_ids, stringsAsFactors = FALSE),
                obj$item_data)
  # Remove any user-defined columns not specified
  data <- data[ , colnames(data) %in% col_subset, drop = FALSE]
  # Check if any of the column-generating functions are needed
  functions <- obj$item_funcs[names(obj$item_funcs) %in% col_subset]
  # Apply column-generating functions and append to output
#   if (calculated_cols && length(functions) > 0) {
#     calculated_data <- lapply(functions, item_apply, obj = obj, item = obj$item_taxon_ids)
#     names(calculated_data) <- names(functions)
#     data <- cbind(data, as.data.frame(calculated_data))
#   }
  # Reorder output to match order of col_subset
  data <- dplyr::tbl_df(data)[row_subset, col_subset]
  
  # Apply drop
  if (ncol(data) == 1 && drop) {
    data = data[[1]]
  }
  
  return(data)
}

#' Split \code{\link{classified}} into individual taxa
#'
#' Splits an object of type \code{\link{classified}} into a list  of
#' \code{\link{classified}} objects, one for each taxon in the input.
#'
#' @param obj (\code{\link{classified}}) The object to split.
#'
#' @return \code{list} of \code{\link{classified}}
#'
#' @export
split_by_taxon <- function(obj) {
  process_one <- function(sub_taxa_ids, super_taxa_ids, taxon_item_ids) {
    new_taxa_id <- c(sub_taxa_ids, super_taxa_ids)
    new_item_id <- obj$item_taxon_ids[taxon_item_ids]
    classified(taxon_ids = new_taxa_id,
               parent_ids =  obj$parent_ids[new_taxa_id],
               item_taxon_ids = new_item_id,
               taxon_data = obj$taxon_data[new_taxa_id, , drop = FALSE],
               item_data = obj$item_data[taxon_item_ids, , drop = FALSE],
               taxon_funcs = obj$taxon_funcs,
               item_funcs = obj$item_funcs)
  }

  mapply(SIMPLIFY = FALSE, process_one,
         subtaxa(obj), supertaxa(obj, include_input = TRUE), items(obj))
}


#' Get all supertaxa of a taxon
#'
#' Return the taxon IDs of all supertaxa (i.e. all taxa the target taxa are a part of) in an
#' object of type \code{classified}
#'
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which supertaxa will be returned.
#' Default: All taxon in \code{obj} will be used.
#' @param recursive (\code{logical})
#' If \code{FALSE}, only return the supertaxa one level above the
#'   target taxa. If \code{TRUE}, return all the supertaxa of every supertaxa, etc.
#' @param simplify (\code{logical})
#' If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs
#' @param include_input (\code{logical})
#' If \code{TRUE}, the input taxa are included in the output
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
#'   
#' @export
supertaxa <- function(obj, subset = obj$taxon_data$taxon_ids, recursive = TRUE,
                      simplify = FALSE, include_input = FALSE) {
  recursive_part <- function(taxon) {
    supertaxon <- obj$taxon_data$parent_ids[taxon]
    if (recursive) {
      if (is.na(supertaxon)) {
        output <- taxon
      } else {
        output <- c(taxon, recursive_part(supertaxon))
      }
    } else {
      output <- c(taxon, supertaxon)
    }
    return(unname(output))
  }

  subset <- format_taxon_subset(obj, subset)
  supertaxa <- stats::setNames(lapply(subset, recursive_part), subset)
  if (!include_input) {
    supertaxa <- lapply(supertaxa, `[`, -1)
    }
  # Reduce dimensionality if specified
  if (simplify) {
    supertaxa <- unlist(supertaxa)
  }
  return(supertaxa)
}


#' Get all subtaxa of a taxon
#'
#' Return the taxon IDs of all subtaxa in an object of type \code{classified}
#'
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which subtaxa will be returned.
#' @param recursive (\code{logical})
#' If \code{FALSE}, only return the subtaxa one level above the
#'   target taxa. If \code{TRUE}, return all the subtaxa of every subtaxa, etc.
#' @param simplify (\code{logical})
#' If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs
#' @param include_input (\code{logical})
#' If \code{TRUE}, the input taxa are included in the output
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
#'
#' @export
subtaxa <- function(obj, subset = obj$taxon_data$taxon_ids, recursive = TRUE,
                    simplify = FALSE, include_input = FALSE) {
  get_children <- function(taxon) {
    unname(obj$taxon_ids[obj$parent_ids == taxon & ! is.na(obj$parent_ids)])
  }
  
  recursive_part <- function(taxon) {
    # Get immediate children of current taxon
    children <- get_children(taxon)
    # Run this function on them to get their output
    child_output <- lapply(children, recursive_part)
    child_output <- stats::setNames(unlist(child_output, recursive = FALSE),
                             unlist(lapply(child_output, names)))
    # Get all subtaxa from the names of the child output
    if (include_input) {
      child_taxa <- c(taxon, names(child_output))
    } else {
      child_taxa <- names(child_output)
      if (is.null(child_taxa)) {
        child_taxa <- character(0)
      }
    }
    # Combine the child output with the subtaxa for the current taxon
    output <- stats::setNames(c(list(child_taxa), child_output),
                       c(taxon, names(child_output)))
    return(output)
  }
  
  
  subset <- format_taxon_subset(obj, subset)  # Get output content
  if (recursive) {
    starting_taxa <- roots(obj, subset)
    output <- unlist(lapply(starting_taxa, recursive_part), recursive = FALSE)[subset]
  } else {
    output <- stats::setNames(lapply(subset, get_children), subset)
    if (include_input) {
      output <- mapply(function(x, n) c(n, x), output, names(output), SIMPLIFY = FALSE)
    }
  }
  # Reduce dimensionality if specified
  if (simplify) {
    output <- unique(unlist(output))
  }
  return(output)
}


#' Get items associated with taxa
#'
#' Given one or more taxa IDs and a \code{\link{classified}} object, return the items
#' (e.g. sequence information) associated with each taxon.
#'
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which items will be returned.
#' @param recursive (\code{logical})
#' If \code{FALSE}, only return the item assigned to the specified input taxa, not subtaxa.
#' If \code{TRUE}, return all the items of every subtaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique item indexes.
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of item indexes are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the item indexes
#'   for all \code{target} taxa are returned in a single vector.
#'
#' @export
items <- function(obj, subset = obj$taxon_data$taxon_ids, recursive = TRUE, simplify = FALSE) {
  # Get output content
  my_subtaxa <- subtaxa(obj, subset, recursive = recursive, include_input = TRUE)
  unique_subtaxa <- unique(unlist(my_subtaxa))
  item_key <- stats::setNames(lapply(unique_subtaxa, function(x) which(x == obj$item_taxon_ids)), unique_subtaxa)
  output <- lapply(my_subtaxa, function(x) unname(unlist(item_key[x])))
  # Reduce dimensionality if specified
  if (simplify) {
    output <- unique(unlist(output))
  }
  return(output)
}


#' Get root taxa
#' 
#' Return the root taxa for a \code{\link{classified}} object.
#' Can also be used to get the roots of a subset of taxa.
#' 
#' @param obj (\code{classified}) The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which supertaxa will be returned.
#' Default: All taxon in \code{obj} will be used.
#' 
#' @return \code{character}
#'  
#' @export
roots <- function(obj, subset = taxon_ids(obj)) {
  parents <- supertaxa(obj, subset = subset, include_input = TRUE)
  is_global_root <- vapply(parents, function(x) length(x) == 1, logical(1))
  if (missing(subset)) {
    is_root <- is_global_root
  } else {
    subset <- format_taxon_subset(obj, subset)
    is_root <- is_global_root | vapply(parents, function(x) ! any(x[-1] %in% subset), logical(1))
  }
  subset[is_root]
}
