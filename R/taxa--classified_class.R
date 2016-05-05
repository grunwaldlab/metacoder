#' Create an instance of \code{classified}
#'
#' Create an instance of \code{classified} containing items classified by a taxonomy
#'
#' @param taxon_id Unique taxon ids
#' @param parent_id Taxon ids of parent taxa of \code{taxon_id}
#' @param item_taxon_id Taxon ids of items
#' @param taxon_data A \code{data.frame} with rows pretaining to \code{taxon_id}
#' @param item_data A \code{data.frame} with rows pretaining to \code{item_taxon_id}
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
classified <- function(taxon_id, parent_id, item_taxon_id,
                       taxon_data = NULL, item_data = NULL,
                       taxon_funcs = list(item_count = count_items, rank = get_rank), item_funcs = list()) {
  #


  # Check that taxon ids are unique
  if (length(unique(taxon_id)) != length(taxon_id)) { stop("'taxon_id' must be unique") }
  # Check that parent_id is the same length of taxon_id
  if (length(taxon_id) != length(parent_id)) { stop("'parent_id' must be the same length as 'taxon_id'") }
  # All parent_id should be in taxon_id
  parent_id[! parent_id %in% taxon_id] <- NA
  # All item_taxon_id should be in taxon_id
  if (any(! item_taxon_id %in% c(taxon_id, NA))) { stop("All 'item_taxon_id' must be in 'taxon_id'") }

  # check taxon_data and item_data
  if (is.null(taxon_data)) {
    taxon_data <- as.data.frame(matrix(numeric(0), nrow = length(taxon_id)))
  }
  if (is.null(item_data)) {
    item_data <- as.data.frame(matrix(numeric(0), nrow = length(item_taxon_id)))
  }
  if (class(taxon_data) != "data.frame") {
    stop("'taxon_data' must a data.frame")
  }
  if (class(item_data) != "data.frame") {
    stop("'item_data' must a data.frame")
  }
  if (nrow(taxon_data) != length(taxon_id)) {
    stop("'taxon_data' must have the same number of rows as 'taxon_id'")
  }
  if (nrow(item_data) != length(item_taxon_id)) {
    stop("'item_data' must have the same number of rows as 'item_taxon_id'")
  }
#   if ( any(colnames(taxon_data) %in% colnames(item_data)) ) {
#     stop("'taxon_data' and 'item_data' may not share column names")
#   }

#   if ( any(colnames(taxon_data) %in% reserved_col_names()) ) {
#     stop(paste("Column names cannot be one of the following:",
#                paste0(reserved_col_names(), collapse = ", ")))
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
  taxon_id <- as.character(taxon_id)
  parent_id <- as.character(parent_id)
  item_taxon_id <- as.character(item_taxon_id)


  # Make object
  rownames(taxon_data) <- taxon_id
  output <- list(taxon_id = setNames(taxon_id, taxon_id),
                 parent_id = setNames(parent_id, taxon_id),
                 item_taxon_id = item_taxon_id,
                 taxon_data = taxon_data,
                 item_data = item_data,
                 taxon_funcs = taxon_funcs,
                 item_funcs = item_funcs)
  class(output) <- "classified"
  return(output)
}

#' invalid column names for user data
#'
#' invalid column names for user data
#' 
#' @keywords internal
reserved_col_names <- function() {
  c("taxon_id", "parent_id", "item_taxon_id")
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
subset.classified <- function(x, taxon, item, subtaxa = TRUE, supertaxa = FALSE, itemless = TRUE, ...) {
  # non-standard argument evaluation
  my_taxon_data <- taxon_data(x)
  column_var_name <- colnames(my_taxon_data)
  unused_result <- lapply(column_var_name, function(x) assign(x, my_taxon_data[[x]], envir = parent.frame(2)))

  # Defaults
#   if (missing(taxon) && missing(item)) {
#     return(x)
#   }
  if (missing(taxon)) {
    taxon_id <- x$taxon_id
  } else {
    parsed_taxon <- tryCatch(suppressWarnings(eval(taxon)), error = function(x) NULL)
    if (is.null(parsed_taxon)) { parsed_taxon <- eval(substitute(taxon)) }
    taxon_id <- x$taxon_id[parsed_taxon]
  }
  parent_id <- x$parent_id[taxon_id]
  if (missing(item)) {
    item_id <- seq_along(x$item_taxon_id)
  } else {
    parsed_item_id <- tryCatch(suppressWarnings(eval(item)), error = function(x) NULL)
    if (is.null(parsed_item_id)) { parsed_item_id <- eval(substitute(item)) }
    item_id <- parsed_item_id
  }

  new_taxa <- taxon_id

  if (subtaxa) {
    new_taxa <- c(new_taxa,
                  get_subtaxa(targets = taxon_id,
                              taxa = x$taxon_id,
                              parents = x$parent_id,
                              recursive = TRUE,
                              simplify = TRUE))
  }

  if (supertaxa) {
    new_taxa <- c(new_taxa,
                  get_supertaxa(targets = taxon_id,
                                taxa = x$taxon_id,
                                parents = x$parent_id,
                                recursive = TRUE,
                                simplify = TRUE,
                                include_target = FALSE))
  }
  new_taxa <- unique(new_taxa)

  inluded_items <- item_id[x$item_taxon_id[item_id] %in% new_taxa]



  # Make output
  output <- classified(taxon_id = new_taxa,
                       parent_id =  x$parent_id[new_taxa],
                       item_taxon_id = x$item_taxon_id[inluded_items],
                       taxon_data = x$taxon_data[new_taxa, , drop = FALSE],
                       item_data = x$item_data[inluded_items, , drop = FALSE],
                       taxon_funcs = x$taxon_funcs,
                       item_funcs = x$item_funcs)

  # Remove taxa with no items
  if (! itemless) {
    taxa_with_items <- taxon_data(output, "item_count") > 0
    output$taxon_id <- output$taxon_id[taxa_with_items]
    output$parent_id <- output$parent_id[taxa_with_items]
    output$taxon_data <- output$taxon_data[taxa_with_items, , drop = FALSE]
  }
  return(output)
}


#' Return taxon data column names
#'
#' Return taxon data column names of and object of type \code{\link{classified}}.
#' This includes "taxon_id", "parent_id", user-defined columns, and columns generated by \code{taxon_funcs}.
#'
#' @param obj (\code{\link{classified}})
#'
#' @export
taxon_data_colnames <- function(obj) {
  c("taxon_id", "parent_id", names(obj$taxon_funcs), colnames(obj$taxon_data))
}

#' Return item data column names
#'
#' Return item data column names of and object of type \code{\link{classified}}.
#' This includes "taxon_id", user-defined columns, and columns generated by \code{item_funcs}.
#'
#' @param obj (\code{\link{classified}})
#'
#' @export
item_data_colnames <- function(obj) {
  c("taxon_id", colnames(obj$item_data))
}




#' Return taxon data from \code{\link{classified}}
#'
#' Return a table of data associated with taxa of and object of type
#' \code{\link{classified}}.
#'
#' @param obj (\code{\link{classified}})
#' @param calculated_cols (\code{logical} of length 1) If \code{TRUE}, return calculated columns using
#' \code{\link{classified}$taxon_funcs}.
#' @param subset (\code{character}) The names of columns, either user defined or generated using \code{taxon_funcs}.
#' Default: All columns.
#' @param drop (\code{logical} of length 1) If \code{TRUE}, if \code{subset} is a single column
#' name, then a \code{vector} is returned instead of a \code{data.frame}
#' @param ... Passed to \code{taxon_data.classified}
#'
#' @return A \code{data.frame} or \code{vector} with rows corresponding to taxa in input
#'
#' @export
#' @rdname taxon_data
taxon_data <- function(...) {
  UseMethod("taxon_data")
}

#' @method taxon_data classified
#' @export
#' @rdname taxon_data
taxon_data.classified <- function(obj,
                                  subset = taxon_data_colnames(obj),
                                  calculated_cols = TRUE,
                                  drop = TRUE, ...) {
  # Check that the user is making sense
  if (calculated_cols == FALSE && any(subset %in% names(obj$taxon_funcs))) {
    stop("Cannot use a calculated column when `calculated_cols = FALSE`.")
  }
  # Combine taxon id information and arbitrary user-defined data
  data <- cbind(data.frame(taxon_id = obj$taxon_id, parent_id = obj$parent_id, stringsAsFactors = FALSE),
                obj$taxon_data)
  # Remove any user-defined columns not specified
  data <- data[ , colnames(data) %in% subset, drop = FALSE]
  # Check if any of the column-generating functions are needed
  functions <- obj$taxon_funcs[names(obj$taxon_funcs) %in% subset]
  # Apply column-generating functions and append to output
  if (calculated_cols && length(functions) > 0) {
    calculated_data <- lapply(functions, taxon_apply, obj = obj, taxon = obj$taxon_id)
    names(calculated_data) <- names(functions)
    data <- cbind(data, as.data.frame(calculated_data))
  }
  # Reorder output to match order of subset
  data <- data[ , subset, drop = drop]
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
#' @param ... Passed to \code{taxon_data.classified}
#'
#' @return \code{data.frame} with rows corresponding to items in input
#'
#' @export
#' @rdname item_data
item_data <- function(...) {
  UseMethod("item_data")
}

#' @method item_data classified
#' @export
#' @rdname item_data
item_data.classified <- function(obj,
                                 subset = item_data_colnames(obj),
                                 # calculated_cols = TRUE,
                                 drop = TRUE, ...) {
  # Check that the user is making sense
#   if (calculated_cols == FALSE && any(subset %in% names(obj$item_funcs))) {
#     stop("Cannot use a calculated column when `calculated_cols = FALSE`.")
#   }
  # Combine taxon id information and arbitrary user-defined data
  data <- cbind(data.frame(taxon_id = obj$item_taxon_id, stringsAsFactors = FALSE),
                obj$item_data)
  # Remove any user-defined columns not specified
  data <- data[ , colnames(data) %in% subset, drop = FALSE]
  # Check if any of the column-generating functions are needed
  functions <- obj$item_funcs[names(obj$item_funcs) %in% subset]
  # Apply column-generating functions and append to output
#   if (calculated_cols && length(functions) > 0) {
#     calculated_data <- lapply(functions, item_apply, obj = obj, item = obj$item_taxon_id)
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
#' @rdname taxon_apply
#' @keywords internal
taxon_apply <- function(obj, func, ...) {
  UseMethod("taxon_apply")
}

#' @method taxon_apply classified
#' @rdname taxon_apply
#' @export
taxon_apply.classified <- function(obj, func, ...) {
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
#' @keywords internal
#' @rdname split_by_taxon
split_by_taxon <- function(obj) {
  UseMethod("split_by_taxon")
}

#' @export
#' @rdname split_by_taxon
#' @method split_by_taxon classified
split_by_taxon.classified <- function(obj) {
  split_sub_taxa <- get_subtaxa(targets = obj$taxon_id,
                                taxa = obj$taxon_id,
                                parents = obj$parent_id,
                                recursive = TRUE,
                                simplify = FALSE)
  split_super_taxa <- get_supertaxa(targets = obj$taxon_id,
                                    taxa = obj$taxon_id,
                                    parents = obj$parent_id,
                                    recursive = TRUE,
                                    simplify = FALSE,
                                    include_target = TRUE)
  split_items <- get_taxon_items(targets = obj$taxon_id,
                                 taxa = obj$taxon_id,
                                 parents = obj$parent_id,
                                 items = obj$item_taxon_id,
                                 recursive = TRUE,
                                 simplify = FALSE)

  process_one <- function(sub_taxa_ids, super_taxa_ids, taxon_item_ids) {
    new_taxa_id <- c(sub_taxa_ids, super_taxa_ids)
    new_item_id <- obj$item_taxon_id[taxon_item_ids]
    classified(taxon_id = new_taxa_id,
               parent_id =  obj$parent_id[new_taxa_id],
               item_taxon_id = new_item_id,
               taxon_data = obj$taxon_data[new_taxa_id, , drop = FALSE],
               item_data = obj$item_data[taxon_item_ids, , drop = FALSE],
               taxon_funcs = obj$taxon_funcs,
               item_funcs = obj$item_funcs)
  }

  mapply(process_one, split_sub_taxa, split_super_taxa, split_items, SIMPLIFY = FALSE)
}

#' Count items in \code{\link{classified}}
#'
#' Count items in \code{\link{classified}}
#'
#' @param obj (\code{\link{classified}})
#' @param taxon (\code{character} of length 1) The taxon_id used
#'
#' @return \code{numeric} of length 1
#'
#' @keywords internal
count_items <- function(obj, taxon) {
  length(obj$item_taxon_id)
}

#' Get rank of taxon in \code{\link{classified}}
#'
#' Get rank of taxon in \code{\link{classified}}
#'
#' @param obj (\code{\link{classified}})
#' @param taxon (\code{character} of length 1) The taxon_id used
#'
#' @return \code{numeric} of length 1
#'
#' @keywords internal
get_rank <- function(obj, taxon) {
  super_taxa <- get_supertaxa(targets = taxon,
                              taxa = obj$taxon_id,
                              parents = obj$parent_id,
                              recursive = TRUE,
                              simplify = TRUE,
                              include_target = TRUE)
  length(super_taxa)
}


#' Combine classified data
#'
#' Combine the contents of multiple objects of type \code{\link{classified}}
#'
#' @param ... (\code{\link{classified}}) One or more objects of type \code{\link{classified}}
#'
#' @return An object of type \code{\link{classified}}
#'
#' @keywords internal
#' @method sum classified
#' @keywords internal
sum.classified <- function(...) {
  input <- list(...)
  input <- input[-length(input)] # remove `na.rm` from `sum` generic
  # Check that all inputs are of correct class
  if (! all(vapply(input, class, character(1)) == "classified")) {
    stop("All inputs must be of type `classified`")
  }

  # Make output object
  new_taxon_data <- do.call(rbind, lapply(input, taxon_data, drop = FALSE))
  new_taxon_data <- new_taxon_data[ , ! colnames(new_taxon_data) %in% reserved_col_names()]
  new_item_data <- do.call(rbind, lapply(input, item_data, drop = FALSE))
  new_item_data <- new_item_data[ , ! colnames(new_item_data) %in% reserved_col_names()]
  classified(taxon_id = unlist(lapply(input, taxon_data, subset = "taxon_id")),
             parent_id = unlist(lapply(input, taxon_data, subset = "parent_id")),
             item_taxon_id = unlist(lapply(input, item_data, subset = "taxon_id")),
             taxon_data = new_taxon_data,
             item_data = new_item_data)
}


#' Count items in \code{\link{classified}}
#'
#' Count items in \code{\link{classified}}
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character}) The \code{taxon_id}s to get counts for.
#'
#' @return \code{numeric} of length 1
#'
#' @keywords internal
item_count <- function(obj, subset = taxon_ids(obj)) {
  recursive_part <- function(taxon) {
    
  }
  
  # Find starting points for recursive search/count (trim the tree)
  subroots <- roots(obj, subset)
  
  # Recursivly count taxa in starting points
  unlist(lapply(subroots, recursive_part), recursive = FALSE)
  
  # ? Filter out non-target counts ?
  
}


#' Get all supertaxa of a taxon
#'
#' Return the taxon IDs of all supertaxa (i.e. all taxa the target taxa are a part of) in an
#' object of type \code{classified}
#'
#' @param obj (\code{classified}) The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which supertaxa will be returned.
#' Default: All taxon in \code{obj} will be used.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the supertaxa one level above the
#'   target taxa. If \code{TRUE}, return all the supertaxa of every supertaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs
#' @param include_target (\code{logical}) If \code{TRUE}, the target taxa are included in the output
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
#'   
#' @keywords internal
supertaxa <- function(obj, subset = taxon_ids(obj), recursive = TRUE, simplify = FALSE, include_target = FALSE) {
  get_supertaxa(targets = subset, taxa = obj$taxon_id, parents = obj$parent_id,
                recursive = recursive, simplify = simplify, include_target = include_target)
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
#' @keywords internal
roots <- function(obj, subset = taxon_ids(obj)) {
  parents <- supertaxa(obj, subset = subset, include_target = TRUE)
  is_global_root <- vapply(parents, function(x) length(x) == 1, logical(1))
  if (missing(subset)) {
    is_root <- is_global_root
  } else {
    is_root <- is_global_root | vapply(parents, function(x) ! any(x[-1] %in% subset), logical(1))
  }
  subset[is_root]
}


#' Get taxon IDs
#' 
#' Return the taxon IDs for a \code{\link{classified}} object.
#' 
#' @param obj (\code{classified}) The \code{classified} object containing taxon information to be queried.
#' 
#' @return \code{character}
#' 
#' @keywords internal
taxon_ids <- function(obj) {
  unname(obj$taxon_id)
}


#' Get parent taxon IDs
#' 
#' Return the taxon IDs for a parents of taxon IDs in an \code{\link{classified}} object.
#' 
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which parents will be returned.
#' Default: All taxon in \code{obj} will be used.
#' 
#' @return \code{character}
#' 
#' @keywords internal
parent_ids <- function(obj, subset = taxon_ids(obj)) {
  if (missing(subset)) {
    parents <- obj$parent_id
  } else {
    parents <- obj$parent_id[subset]
  }
  unname(parents)
}
