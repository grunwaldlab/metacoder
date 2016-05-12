#' Create an instance of \code{classified}
#'
#' Create an instance of \code{classified} containing items classified by a taxonomy
#'
#' @param taxon_ids Unique taxon ids
#' @param parent_ids Taxon ids of parent taxa of \code{taxon_ids}
#' @param item_taxon_ids Taxon ids of items
#' @param taxon_data A \code{data.frame} with rows pretaining to \code{taxon_ids}
#' @param item_data A \code{data.frame} with rows pretaining to \code{item_taxon_ids}
#' @param taxon_funcs A function that accepts a subset of this object containing information for a single taxa.
#' A single value must be returned derived from that information.
#' These values will be acessible as if it were a column in \code{taxon_data}.
#' @param item_funcs A function that accepts a subset of this object containing information for a single item.
#' A single value must be returned derived from that information.
#' These values will be acessible as if it were a column in \code{item_data}.
#'
#' @return An object of type \code{classified}
#'
#' @export
classified <- function(taxon_ids, parent_ids, item_taxon_ids,
                       taxon_data = NULL, item_data = NULL,
                       taxon_funcs = list(item_counts = item_counts, taxon_ranks = taxon_ranks, classifications = classifications),
                                          # child_counts = child_counts, subtaxon_counts = subtaxon_counts),
                       item_funcs = list()) {
  # Check that taxon ids are unique
  if (length(unique(taxon_ids)) != length(taxon_ids)) { stop("'taxon_ids' must be unique") }
  # Check that parent_ids is the same length of taxon_ids
  if (length(taxon_ids) != length(parent_ids)) { stop("'parent_ids' must be the same length as 'taxon_ids'") }
  # All parent_ids should be in taxon_ids
  parent_ids[! parent_ids %in% taxon_ids] <- NA
  # All item_taxon_ids should be in taxon_ids
  if (any(! item_taxon_ids %in% c(taxon_ids, NA))) { stop("All 'item_taxon_ids' must be in 'taxon_ids'") }

  # check taxon_data and item_data
  if (is.null(taxon_data)) {
    taxon_data <- as.data.frame(matrix(numeric(0), nrow = length(taxon_ids)))
  }
  if (is.null(item_data)) {
    item_data <- as.data.frame(matrix(numeric(0), nrow = length(item_taxon_ids)))
  }
  if (class(taxon_data) != "data.frame") {
    stop("'taxon_data' must a data.frame")
  }
  if (class(item_data) != "data.frame") {
    stop("'item_data' must a data.frame")
  }
  if (nrow(taxon_data) != length(taxon_ids)) {
    stop("'taxon_data' must have the same number of rows as 'taxon_ids'")
  }
  if (nrow(item_data) != length(item_taxon_ids)) {
    stop("'item_data' must have the same number of rows as 'item_taxon_ids'")
  }
  reserved_col_names <- c("taxon_ids", "parent_ids", "item_taxon_ids")
#   if ( any(colnames(taxon_data) %in% reserved_col_names ) {
#     stop(paste("Column names cannot be one of the following:",
#                paste0(reserved_col_names, collapse = ", ")))
#   }

  # Check column functions
  is_named <- function(x) {
    (! is.null(names(x))) && all(names(x) != '')
  }
  if ( length(taxon_funcs) > 1 && (! all(sapply(taxon_funcs, is.function)) || ! is_named(taxon_funcs)) ) {
    stop("'taxon_funcs' must all be named functions")
  }
  if ( length(item_funcs) > 1 && (! all(sapply(item_funcs, is.function)) || ! is_named(item_funcs)) ) {
    stop("'item_funcs' must all be named functions")
  }

  # Ensure correct type
  taxon_ids <- as.character(taxon_ids)
  parent_ids <- as.character(parent_ids)
  item_taxon_ids <- as.character(item_taxon_ids)


  # Make object
  rownames(taxon_data) <- taxon_ids
  output <- list(taxon_ids = setNames(taxon_ids, taxon_ids),
                 parent_ids = setNames(parent_ids, taxon_ids),
                 item_taxon_ids = item_taxon_ids,
                 taxon_data = taxon_data,
                 item_data = item_data,
                 taxon_funcs = taxon_funcs,
                 item_funcs = item_funcs)
  class(output) <- "classified"
  return(output)
}


#' Create a inclusive subset of \code{\link{classified}}
#'
#' Create a subset of items classified by a taxonomy.
#' Only unspecified taxa with no items or children with items are discarded.
#'
#' @param x \code{\link{classified}}
#' @param taxon A key to filter the taxon data.frame rows on
#' @param item A key to filter the item data.frame rows on
#' @param subtaxa (\code{logical} of length 1) If \code{TRUE}, return subtaxa of specified taxa
#' @param supertaxa (\code{logical} of length 1) If \code{TRUE}, return supertaxa of specified taxa
#' @param itemless (\code{logical} of length 1) If \code{TRUE}, return taxa even if they have no items assigned to them
#' @param ... not used
#'
#' @return \code{\link{classified}}
#'
#' @export
subset.classified <- function(x, taxon = taxon_ids(x), item = seq_along(x$item_taxon_ids),
                              subtaxa = TRUE, supertaxa = FALSE, itemless = TRUE, ...) {
  # non-standard argument evaluation
  parsed_taxon <- lazyeval::lazy_eval(lazyeval::lazy(taxon), data = taxon_data(x))
  parsed_item <- lazyeval::lazy_eval(lazyeval::lazy(item), data = item_data(x))
  
  # Get taxa of subset
  new_taxa <- unique(c(x$taxon_ids[parsed_taxon],
                       if (subtaxa) {
                         subtaxa(x, subset = parsed_taxon, simplify = TRUE)
                       },
                       if (supertaxa) {
                         supertaxa(x, subset = parsed_taxon, simplify = TRUE, include_input = FALSE)
                       }))

  # Get items of subset
  inluded_items <- x$item_taxon_ids[parsed_item] %in% new_taxa

  # Make output
  output <- classified(taxon_ids = new_taxa,
                       parent_ids =  x$parent_ids[new_taxa],
                       item_taxon_ids = x$item_taxon_ids[inluded_items],
                       taxon_data = x$taxon_data[new_taxa, , drop = FALSE],
                       item_data = x$item_data[inluded_items, , drop = FALSE],
                       taxon_funcs = x$taxon_funcs,
                       item_funcs = x$item_funcs)

  # Remove taxa with no items
  if (! itemless) {
    taxa_with_items <- item_counts(output) > 0
    output$taxon_ids <- output$taxon_ids[taxa_with_items]
    output$parent_ids <- output$parent_ids[taxa_with_items]
    output$taxon_data <- output$taxon_data[taxa_with_items, , drop = FALSE]
  }
  return(output)
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
                       sort_by = "classifications",
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
  data <- data[ , col_subset, drop = drop]
  
  # Reorder output rows according to `sort_by`
  if (! is.null(sort_by)) {
    if (sort_by %in% colnames(data)) {
      sort_by_col <- data[ , sort_by]
    } else if (is.function(sort_by)) {
      sort_by_col <- sort_by(obj, subset)
    } else if (sort_by %in% names(obj$taxon_funcs)) {
      sort_by_col <- obj$taxon_funcs[[sort_by]](obj, subset)
    } else {
      sort_by_col <- get(sort_by)(obj, subset)
    }
    data <- data[order(sort_by_col, decreasing = decreasing), ]
  }
  
  return(data)
}


#' Return item data from \code{\link{classified}}
#'
#' Return a table of data associated with items of an object of type
#' \code{\link{classified}}.
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character}) The names of columns, either user defined or generated using \code{item_funcs}.
#' Default: All columns.
#' @param drop (\code{logical} of length 1) If \code{TRUE}, if \code{subset} is a single column
#' name, then a \code{vector} is returned instead of a \code{data.frame}
#'
#' @return \code{data.frame} with rows corresponding to items in input
#'
#' @export
item_data <- function(obj,
                      subset = item_data_colnames(obj),
                      # calculated_cols = TRUE,
                      drop = FALSE) {
  # Check that the user is making sense
#   if (calculated_cols == FALSE && any(subset %in% names(obj$item_funcs))) {
#     stop("Cannot use a calculated column when `calculated_cols = FALSE`.")
#   }
  # Combine taxon id information and arbitrary user-defined data
  data <- cbind(data.frame(taxon_ids = obj$item_taxon_ids, stringsAsFactors = FALSE),
                obj$item_data)
  # Remove any user-defined columns not specified
  data <- data[ , colnames(data) %in% subset, drop = FALSE]
  # Check if any of the column-generating functions are needed
  functions <- obj$item_funcs[names(obj$item_funcs) %in% subset]
  # Apply column-generating functions and append to output
#   if (calculated_cols && length(functions) > 0) {
#     calculated_data <- lapply(functions, item_apply, obj = obj, item = obj$item_taxon_ids)
#     names(calculated_data) <- names(functions)
#     data <- cbind(data, as.data.frame(calculated_data))
#   }
  # Reorder output to match order of subset
  data <- data[ , subset, drop = drop]
  return(data)
}


#' Apply a function to every taxon
#'
#' Apply  a function to every taxon in an object of type \code{\link{classified}}.
#'
#' @param obj (\code{\link{classified}})
#' @param func (\code{function}) A function that accepts an object of type \code{\link{classified}}.
#' @param ... Passed to \code{mapply}
#'
#' @return \code{list} of length equal to the number of taxa in \code{obj}
#'
#' @export
taxon_apply <- function(obj, func, ...) {
  split_data <- split_by_taxon(obj)
  mapply_args <- c(list(func), list(split_data), list(...))
  unlist(do.call(mapply, mapply_args), recursive = FALSE)
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


#' Get rank of taxa 
#'
#' Get rank of taxa in an object of type \code{\link{classified}}
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character}) The \code{taxon_ids}s to get ranks for.
#'
#' @return \code{numeric}
#'
#' @export
taxon_ranks <- function(obj, subset = taxon_ids(obj)) {
  vapply(supertaxa(obj, subset, include_input = TRUE), length, numeric(1))
}


#' Count items in \code{\link{classified}}
#'
#' Count items in \code{\link{classified}}
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character}) The \code{taxon_ids}s to get counts for.
#'
#' @return \code{numeric}
#'
#' @export
item_counts <- function(obj, subset = taxon_ids(obj)) {
  vapply(items(obj, subset), length, numeric(1))
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
supertaxa <- function(obj, subset = taxon_ids(obj), recursive = TRUE,
                      simplify = FALSE, include_input = FALSE) {
  recursive_part <- function(taxon) {
    supertaxon <- obj$parent_ids[taxon]
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
  supertaxa <- setNames(lapply(subset, recursive_part), subset)
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
subtaxa <- function(obj, subset = taxon_ids(obj), recursive = TRUE,
                    simplify = FALSE, include_input = FALSE) {
  get_children <- function(taxon) {
    unname(obj$taxon_ids[obj$parent_ids == taxon & ! is.na(obj$parent_ids)])
  }
  
  recursive_part <- function(taxon) {
    # Get immediate children of current taxon
    children <- get_children(taxon)
    # Run this function on them to get their output
    child_output <- lapply(children, recursive_part)
    child_output <- setNames(unlist(child_output, recursive = FALSE),
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
    output <- setNames(c(list(child_taxa), child_output),
                       c(taxon, names(child_output)))
    return(output)
  }
  
  
  subset <- format_taxon_subset(obj, subset)  # Get output content
  if (recursive) {
    starting_taxa <- roots(obj, subset)
    output <- unlist(lapply(starting_taxa, recursive_part), recursive = FALSE)[subset]
  } else {
    output <- setNames(lapply(subset, get_children), subset)
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
items <- function(obj, subset = taxon_ids(obj), recursive = TRUE, simplify = FALSE) {
  # Get output content
  my_subtaxa <- subtaxa(obj, subset, recursive = recursive, include_input = TRUE)
  unique_subtaxa <- unique(unlist(my_subtaxa))
  item_key <- setNames(lapply(unique_subtaxa, function(x) which(x == obj$item_taxon_ids)), unique_subtaxa)
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


#' Get taxon IDs
#' 
#' Return the taxon IDs for a \code{\link{classified}} object.
#' 
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' 
#' @return \code{character}
#' 
#' @export
taxon_ids <- function(obj) {
  unname(obj$taxon_ids)
}


#' Get classification of taxa 
#'
#' Get classification strings of taxa in an object of type \code{\link{classified}}.
#' Each classification is constructed by concatenating the taxon ids of the given taxon and its supertaxa.
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character})
#' The \code{taxon_ids}s to get classifications for.
#' @param sep (\code{character} of length 1)
#' The character(s) to place between taxon IDs
#'
#' @return \code{character} of length equal to \code{subset}
#'
#' @export
classifications <- function(obj, subset = taxon_ids(obj), sep = ";") {
  vapply(supertaxa(obj, subset, include_input = TRUE), paste0, character(1), collapse = sep)
}



#' Format taxon subset value
#' 
#' Format an input to a \code{subset} option on functions like \code{subset}.
#' Converts logical and numeric vectors into taxon ids
#' 
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param index
#' 
#' @return \code{character}
#' 
#' @keywords internal
format_taxon_subset <- function(obj, index) {
  if (is.logical(index) || is.numeric(index)) {
    index <- taxon_ids(obj)[index]
  }
  index <- unname(index)
  return(index)
}
