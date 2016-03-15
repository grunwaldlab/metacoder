f <- function(x) {
  sfsddsdeognmvdfs = 2
  eval(substitute(x), envir = environment())
}



#' Create an instance of \code{classified}
#' 
#' Create an instance of \code{classified} containing items classified by a taxonomy
#' 
#' @param taxon_id Unique taxon ids
#' @param parent_id Parent taxa of \code{taxon_id}
#' @param item_taxon_id Taxon ids of items
#' @param taxon_data A \code{data.frame} with rows pretaining to \code{taxon_id} 
#' @param item_data A \code{data.frame} with rows pretaining to \code{item_taxon_id}
#' @param taxon_funcs A function that accepts a subset of this object containing information for a single taxa.
#' A single value must be returned derived from that information. 
#' These values will be acessible as if it were a column in \code{taxon_data}.
#' @param item_funcs A function that accepts a subset of this object containing information for a single item.
#' A single value must be returned derived from that information. 
#' These values will be acessible as if it were a column in \code{item_data}.
classified <- function(taxon_id, parent_id, item_taxon_id,
                       taxon_data = NULL, item_data = NULL,
                       taxon_funcs = list(), item_funcs = list()) {
  # Check that taxon ids are unique
  if (unique(taxon_id) == taxon_id) { stop("'taxon_id' must be unique") }
  # Check that parent_id is the same length of taxon_id
  if (length(taxon_id) != length(parent_id)) { stop("'parent_id' must be the same length as 'taxon_id'") }
  # All parent_id should be in taxon_id
  parent_id[! parent_id %in% taxon_id] <- NA
  # All item_taxon_id should be in taxon_id
  if (any(! item_taxon_id %in% taxon_id)) { stop("All 'item_taxon_id' must be in 'taxon_id'") }
  
  # check taxon_data and item_data
  if (is.null(taxon_data)) {
    taxon_data <- as.data.frame(matrix(numeric(0), nrow = nrow(taxon_data)))
  }
  if (is.null(item_data)) {
    item_data <- as.data.frame(matrix(numeric(0), nrow = nrow(item_data)))
  }
  if (nrow(taxon_data) != length(taxon_id)) {
    stop("'taxon_data' must have the same number of rows as 'taxon_id'")
  }
  if (nrow(item_data) != length(item_taxon_id)) {
    stop("'item_data' must have the same number of rows as 'item_taxon_id'")
  }
  if ( any(colnames(taxon_data) %in% colnames(item_data)) ) {
    stop("'taxon_data' and 'item_data' may not share column names")
  }
  reserved_col_names <- c("taxon_id", "parent_id", "item_taxon_id")
  if ( any(colnames(taxon_data) %in% reserved_col_names) ) {
    stop(paste("Column names cannot be one of the following:",
               paste0(reserved_col_names, collapse = ", ")))
  }
  
  # Check column functions
  is_named <- function(x) {
    (! is.null(names(x))) && all(names(x) != '')
  }
  if ( (! all(sapply(taxon_funcs, is.function))) || (! is_named(taxon_funcs)) ) {
    stop("'taxon_funcs' must all be named functions")
  }
  if ( (! all(sapply(item_funcs, is.function))) || (! is_named(item_funcs)) ) {
    stop("'item_funcs' must all be named functions")
  }
  
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



#' Create a restrictive subset of \code{\link{classified}}
#' 
#' Create a subset of items classified by a taxonomy.
#' Only specified taxa and shared parent taxa of specified taxa with specified items are preserved.
#' Only specified items assigned to specified taxa are preserved.
#' 
#' @param taxon A key to filter the taxon data.frame rows on
#' @param item A key to filter the item data.frame rows on
#' 
#' @return \code{\link{classified}}
`[[.classified` <- function(taxon, item) {
  
}


#' Create a inclusive subset of \code{\link{classified}}
#' 
#' Create a subset of items classified by a taxonomy.
#' Only unspecified taxa with no items or children with items are discarded.
#' 
#' @param taxon A key to filter the taxon data.frame rows on
#' @param item A key to filter the item data.frame rows on
#' 
#' @return \code{\link{classified}}
`[.classified` <- function(obj, taxon, item) {
  if (is.missing(taxon) && is.missing(item)) {
    return(obj)
  }
  if (is.missing(taxon)) {
    taxon_id <- obj$taxon_id
  } else {
    taxon_id <- obj$taxon_id[taxon]
  }
  parent_id <- obj$parent_id[taxon_id]
  if (is.missing(item)) {
    item_id <- obj$taxon_item_id
  } else {
    item_id <- obj$taxon_item_id[item]
  }
  new_sub_taxa <- get_subtaxa(targets = taxon_id,
                              taxa = obj$taxon_id,
                              parents = obj$parent_id,
                              recursive = TRUE,
                              simplify = TRUE)
  new_super_taxa <- get_supertaxa(targets = taxon_id,
                                  taxa = obj$taxon_id,
                                  parents = obj$parent_id,
                                  recursive = TRUE,
                                  simpliTfy = TRUE,
                                  include_target = TRUE)
  new_items <- get_taxon_items(targets = taxon_id,
                               taxa = obj$taxon_id,
                               parents = obj$parent_id,
                               items = item_id,
                               recursive = TRUE,
                               simplify = TRUE)
  new_taxa <-  c(new_sub_taxa, new_super_taxa)
  
  classified(taxon_id = new_taxa,
             parent_id =  obj$parent_id[new_taxa],
             item_taxon_id = new_items,
             taxon_data = obj$taxon_data[new_taxa, ],
             item_data = obj$item_data[new_items, ],
             taxon_funcs = obj$taxon_funcs,
             item_funcs = obj$item_funcs)
  
}




taxon_data.classified <- function(obj, col_names = NULL, calculated_cols = TRUE) {
  data <- cbind(data.frame(taxon_id = taxon_id, parent_id = parent_id),
                obj$taxon_data)
  if (calculated_cols) {
    calculated_data <- lapply(obj$taxon_funcs, taxon_apply)
    names(calculated_data) <- names(obj$taxon_funcs)
    data <- cbind(data, as.data.frame(calculated_data))
  }
  if (! is.null(col_names)) {
    data <- data[ , col_names]
  }
  return(data)
}



taxon_apply.classfied <- function(obj, func) {
  unlist(lapply(split_by_taxon(obj), func), recursive = FALSE)
}




item_apply.classfied <- function(obj, func) {
  unlist(lapply(split_by_item(obj), func), recursive = FALSE)
}


split_by_taxon.classfied <- function(obj) {
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
    classified(taxon_id = new_taxa_id,
               parent_id =  obj$parent_id[new_taxa_id],
               item_taxon_id = taxon_item_ids,
               taxon_data = obj$taxon_data[new_taxa_id, ],
               item_data = obj$item_data[taxon_item_ids, ],
               taxon_funcs = obj$taxon_funcs,
               item_funcs = obj$item_funcs)
  }
  
  mapply(process_one, split_sub_taxa, split_super_taxa, split_items, SIMPLIFY = FALSE)
}




split_by_item.classfied <- function(obj) {
  
}