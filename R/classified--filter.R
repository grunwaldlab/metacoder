#' @export
#' @rdname filter
filter <- function(.data, ...) {
  UseMethod("filter")
}


#' @export
#' @rdname filter
filter.default <- function(.data, ...) {
  dplyr::filter(.data, ...)
}

#' Filter \code{\link{classified}} on a list of conditions
#'
#' Create a subset of a \code{\link{classified}} object.
#' Can filter based on columns in \code{.data$taxon_data} or \code{.data$taxon_data}.
#'
#' @param ... One or more filtering conditions.
#' This can be one of two things:
#' \describe{
#'    \item{\code{integer}}{One or more indexes of \code{item_data}}
#'    \item{\code{logical}}{A \code{TRUE}/\code{FALSE} vector of length equal to the number of rows in \code{item_data}}
#'  }
#' Any column name that appears in \code{item_data(.data)} or \code{taxon_data(.data)} can be used as if it was a vector on its own.
#' However, olumn names from \code{item_data(.data)} and \code{taxon_data(.data)} cannot be mixed in a single condition.
#' 
#' @inheritParams filter_taxa
#' @inheritParams filter_items
#' @inheritParams filter_.classified
#' 
#'
#' @return An object of type \code{\link{classified}}
#'
#' @export
#' @rdname filter
filter.classified <- function(.data, ...,
                              subtaxa = TRUE, supertaxa = FALSE,
                              itemless = TRUE, taxonless = FALSE, reassign = TRUE) {
  filter_.classified(.data, ..., .dots = lazyeval::lazy_dots(...),
                     subtaxa = subtaxa, supertaxa = supertaxa,
                     itemless = itemless, taxonless = taxonless,reassign = TRUE)
}


#' @param .dots
#' A list of values to filter on. Used for standard evaluation.
#' 
#' @export
#' @rdname filter
filter_.classified <- function(.data, ..., .dots,
                               subtaxa = TRUE, supertaxa = FALSE,
                               itemless = TRUE, taxonless = FALSE, reassign = TRUE) {
  var_names_in_calls <- function(.dots) {
    components <- lapply(.dots, function(x) as.character(x$expr))
    lapply(components,
           function(x) x[x %in% c(colnames(.data$taxon_data), colnames(.data$item_data))])
  }
  
  
  # non-standard argument evaluation
  parsed_data <- lazyeval::lazy_eval(.dots, data = taxon_data(.data))
  
  # Check that taxon and item columns are not mixed
  cols_used <- var_names_in_calls(.dots)
  in_taxon_data <- vapply(cols_used, function(x) any(x %in% colnames(.data$taxon_data)), logical(1))
  in_item_data <- vapply(cols_used, function(x) any(x %in% colnames(.data$item_data)), logical(1))
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
  
  # Filter by taxon 
  
  # Filter by item
  
  return(.data)
}




#' Filter taxa with a list of conditions
#' 
#' Filter taxa in a \code{\link{classified}} object with a list of conditions.
#' Any column name that appears in \code{taxon_data(.data)} can be used as if it was a vector on its own.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more filtering conditions.
#' This can be one of three things:
#' \describe{
#'    \item{\code{character}}{One or more \code{taxon_id}s}
#'    \item{\code{integer}}{One or more indexes of \code{taxon_data}}
#'    \item{\code{logical}}{A \code{TRUE}/\code{FALSE} vector of length equal to the number of rows in \code{taxon_data}}
#'  }
#' Any column name that appears in \code{taxon_data(.data)} can be used as if it was a vector on its own.
#' @param subtaxa (\code{logical} of length 1)
#' If \code{TRUE}, include subtaxa of taxa passing the filter.
#' @param supertaxa (\code{logical} of length 1)
#' If \code{TRUE}, include supertaxa of taxa passing the filter.
#' @param taxonless (\code{logical} of length 1)
#' If \code{TRUE}, include items even if the taxon they are assigned to is filtered out.
#' Item assigned to removed taxa will be assigned to \code{NA}.
#' See the \code{reassign} option below for further complications.
#' @param reassign (\code{logical} of length 1)
#' If \code{TRUE}, items assigned to removed taxa will be reassigned to the closest supertaxon that passed the filter.
#' If there are no supertaxa of such an item that passed the filter, they will be filtered out if \code{taxonless} is \code{TRUE}.
#' 
#' @return An object of type \code{\link{classified}}
#' 
#' @export
filter_taxa <- function(.data, ..., subtaxa = TRUE, supertaxa = FALSE,
                        taxonless = FALSE, reassign = TRUE) {
  
  # non-standard argument evaluation ---------------------------------------------------------------
  selection <- lazyeval::lazy_eval(lazyeval::lazy_dots(...), data = taxon_data(.data)) 
  
  # convert taxon_ids to logical -------------------------------------------------------------------
  is_char <- vapply(selection, is.character, logical(1))
  selection[is_char] <- lapply(selection[is_char], function(x) .data$taxon_data$taxon_ids %in% x)
  
  # convert indexes to logical ---------------------------------------------------------------------
  is_index <- vapply(selection, is.numeric, logical(1))
  selection[is_index] <- lapply(selection[is_index], function(x) 1:nrow(.data$taxon_data) %in% x)
  
  # combine filters --------------------------------------------------------------------------------
  selection <- Reduce(`&`, selection)
  
  # Get taxa of subset -----------------------------------------------------------------------------
  taxa_subset <- unique(c(.data$taxon_data$taxon_ids[selection],
                          if (subtaxa) {
                            subtaxa(.data, subset = selection, simplify = TRUE)
                          },
                          if (supertaxa) {
                            supertaxa(.data, subset = selection, simplify = TRUE, include_input = FALSE)
                          }))
  
  # Reassign taxonless items -----------------------------------------------------------------------
  if (reassign) {
    reassign_one <- function(x) {
      parents <- supertaxa(.data, subset = x, simplify = TRUE)
      included_parents <- parents[parents %in% taxa_subset]
      if (length(included_parents) > 0) {
        return(included_parents[1])
      } else {
        return(as.numeric(NA))
      }
    }
    
    to_reassign <- ! .data$item_data$item_taxon_ids %in% taxa_subset
    .data$item_data[to_reassign, "item_taxon_ids"] <- vapply(.data$item_data[to_reassign, "item_taxon_ids"], 
                                                             reassign_one, numeric(1))
  }
  
  # Remove taxonless items -------------------------------------------------------------------------
  if (! taxonless) {
    .data$item_data <- .data$item_data[.data$item_data$item_taxon_ids %in% taxa_subset, , drop = FALSE]
  }
  
  # Remove filtered taxa ---------------------------------------------------------------------------
  .data$taxa <- .data$taxa[taxa_subset]
  .data$taxon_data <- .data$taxon_data[.data$taxon_data$taxon_ids %in% taxa_subset, , drop = FALSE]
  
  # Rename taxon ids -------------------------------------------------------------------------------
  custom_which <- function(x, data) {
    if (is.na(x)) {
      return(as.numeric(NA))
    } else {
      out <- which(data$taxa == x)
      if (length(out) == 0) {
        return(as.numeric(NA))
      } else {
        return(out)
      }
    }
  }
  .data$taxon_data$taxon_ids <- vapply(.data$taxon_data$taxon_ids, function(x) custom_which(x, .data), 
                                       numeric(1))
  .data$taxon_data$parent_ids <- vapply(.data$taxon_data$parent_ids, function(x) custom_which(x, .data), 
                                        numeric(1))
  .data$item_data$item_taxon_ids <- vapply(.data$item_data$item_taxon_ids, function(x) custom_which(x, .data), 
                                           numeric(1))
  
  return(.data)
}



#' Filter items with a list of conditions
#' 
#' Filter items in a \code{\link{classified}} object with a list of conditions.
#' Any column name that appears in \code{item_data(.data)} can be used as if it was a vector on its own.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more filtering conditions.
#' This can be one of two things:
#' \describe{
#'    \item{\code{integer}}{One or more indexes of \code{item_data}}
#'    \item{\code{logical}}{A \code{TRUE}/\code{FALSE} vector of length equal to the number of rows in \code{item_data}}
#'  }
#' Any column name that appears in \code{item_data(.data)} can be used as if it was a vector on its own.
#' @param itemless (\code{logical} of length 1)
#' If \code{TRUE}, preserve taxa even if all of their items are filtered out.
#' 
#' @return An object of type \code{\link{classified}}
#' 
#' @export
filter_items <- function(.data, ..., itemless = TRUE) {
}

