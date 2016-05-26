#' UNDER CONSTRUCTION
#' 
#' UNDER CONSTRUCTION
#' 
#' @export
filter <- function(.data, ...) {
  UseMethod("filter")
}


#' @export
filter.default <- function(.data, ...) {
  dplyr::filter(.data, ...)
}

#' Filter \code{\link{classified}} on a list of conditions
#'
#' Create a subset of a \code{\link{classified}} object.
#' Can filter based on columns in \code{obj$taxon_data} or \code{obj$taxon_data}.
#'
#' @param .data \code{\link{classified}}
#' @param ...
#' Filtering conditions.
#' Each condition of must contain the name of at least one column from \code{obj$taxon_data} or \code{obj$taxon_data}.
#' To filter by index or \code{TRUE}/\code{FALSE} vector, use \code{\link{filter_taxa}} or \code{\link{filter_items}}.
#' @param subtaxa (\code{logical} of length 1)
#' If \code{TRUE}, return subtaxa of specified taxa.
#' @param supertaxa (\code{logical} of length 1)
#' If \code{TRUE}, return supertaxa of specified taxa.
#' @param itemless (\code{logical} of length 1)
#' If \code{TRUE}, return taxa even if they have no items assigned to them.
#' @param taxonless (\code{logical} of length 1)
#' If \code{TRUE}, return items even if they are not assigned to taxa.
#'
#' @return \code{\link{classified}}
#'
#' @export
#' @rdname filter
filter.classified <- function(.data, ...,
                              subtaxa = TRUE, supertaxa = FALSE,
                              itemless = TRUE, taxonless = FALSE) {
  filter_.classified(.data, ..., .dots = lazyeval::lazy_dots(...),
                     subtaxa = subtaxa, supertaxa = supertaxa,
                     itemless = itemless, taxonless = taxonless)
}


#' Standard evaluation version of \code{\link{filter}}
#' 
#' @param .dots
#' Filtering conditions.
#' 
#' @inheritParams filter
#' 
#' @export
#' @rdname filter
filter_.classified <- function(.data, ..., .dots,
                               subtaxa = TRUE, supertaxa = FALSE,
                               itemless = TRUE, taxonless = FALSE) {
  var_names_in_calls <- function(.dots) {
    components <- lapply(.dots, function(x) as.character(x$expr))
    lapply(components,
           function(x) x[x %in% c(colnames(.data$taxon_data), colnames(.data$item_data))])
  }
  
  
  # non-standard argument evaluation
  parsed_data <- lazyeval::lazy_eval(.dots, data = taxon_data(obj))
  
  # Check that taxon and item columns are not mixed
  cols_used <- var_names_in_calls(.dots)
  in_taxon_data <- vapply(cols_used, function(x) any(x %in% colnames(obj$taxon_data)), logical(1))
  in_item_data <- vapply(cols_used, function(x) any(x %in% colnames(obj$item_data)), logical(1))
  invalid_mixtures <- in_taxon_data & in_item_data
  if (sum(invalid_mixtures) > 0) {
    char_conditions <- vapply(.dots, function(x) deparse(x$expr), character(1))
    invalid_list <- paste("   ", which(invalid_mixtures), ": ", char_conditions[invalid_mixtures], "\n")
    stop(paste0(collapse = "",
                c("The following filtering conditions are an invalid mixture of taxon_data and item_data column names:\n",
                  invalid_list)))
  }
  
  # Combine taxa filters
  taxa_selection <- Reduce(`&`, parsed_data[in_taxon_data])

  # Combine item filters
  item_selection <- Reduce(`&`, parsed_data[in_item_data])
  
  # Get taxa of subset
  new_taxa <- unique(c(obj$taxon_data$taxon_ids[taxa_selection],
                       if (subtaxa) {
                         subtaxa(obj, subset = taxa_selection, simplify = TRUE)
                       },
                       if (supertaxa) {
                         supertaxa(obj, subset = taxa_selection, simplify = TRUE, include_input = FALSE)
                       }))
  
  # Get items of subset
  inluded_items <- intersect(which(obj$item_data$item_taxon_ids %in% new_taxa),
                             item_selection)
  
  # Apply selection
  obj$taxa <- obj$taxa[new_taxa]
  obj$taxon_data <- obj$taxon_data[obj$taxon_data$taxon_ids %in% new_taxa, , drop = FALSE]
  obj$item_data <- obj$item_data[obj$item_data$item_taxon_ids %in% inluded_items, , drop = FALSE]
  
  # Remove taxa with no items
  if (! itemless) {
    taxa_with_items <- item_counts(obj) > 0
    obj$taxa <- obj$taxa[taxa_with_items]
    obj$taxon_data <- obj$taxon_data[obj$taxon_data$taxon_ids %in% taxa_with_items, , drop = FALSE]
  }
  
  # Remove items with no taxa
  if (! taxonless) {
    obj$item_data <- obj$taxon_data[obj$item_data$item_taxon_ids %in% obj$taxon_data$taxon_ids, , drop = FALSE]
  }
  
  return(obj)
}




#' UNDER CONSTRUCTION
#' 
#' UNDER CONSTRUCTION
#' 
#' @param .data \code{\link{classified}}
#' @param ... Filtering conditions.
#' 
#' @export
filter_taxa <- function(.data, ...) {
}



#' UNDER CONSTRUCTION
#' 
#' UNDER CONSTRUCTION
#' 
#' @param .data \code{\link{classified}}
#' @param ... Filtering conditions.
#' 
#' @export
filter_items <- function(.data, ...) {
}

