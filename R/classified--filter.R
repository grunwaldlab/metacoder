#' Filter taxa with a list of conditions
#' 
#' Filter taxa in a \code{\link{classified}} object with a list of conditions. Any column name that
#' appears in \code{taxon_data(.data)} can be used as if it was a vector on its own. See
#' \code{\link[dplyr]{filter}} for inspiration and more information.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more filtering conditions. This can be one of three things: \describe{
#'   \item{\code{character}}{One or more \code{taxon_id}s} \item{\code{integer}}{One or more indexes
#'   of \code{taxon_data}} \item{\code{logical}}{A \code{TRUE}/\code{FALSE} vector of length equal
#'   to the number of rows in \code{taxon_data}} } Any column name that appears in 
#'   \code{taxon_data(.data)} can be used as if it was a vector on its own.
#' @param subtaxa (\code{logical} of length 1) If \code{TRUE}, include subtaxa of taxa passing the
#'   filter.
#' @param supertaxa (\code{logical} of length 1) If \code{TRUE}, include supertaxa of taxa passing
#'   the filter.
#' @param taxonless (\code{logical} of length 1) If \code{TRUE}, include items even if the taxon
#'   they are assigned to is filtered out. Item assigned to removed taxa will be assigned to
#'   \code{NA}. See the \code{reassign} option below for further complications.
#' @param reassign (\code{logical} of length 1) If \code{TRUE}, items assigned to removed taxa will
#'   be reassigned to the closest supertaxon that passed the filter. If there are no supertaxa of
#'   such an item that passed the filter, they will be filtered out if \code{taxonless} is
#'   \code{TRUE}.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
filter_taxa <- function(.data, ..., subtaxa = TRUE, supertaxa = FALSE,
                        taxonless = FALSE, reassign = TRUE) {
  
  # non-standard argument evaluation ---------------------------------------------------------------
  selection <- lazyeval::lazy_eval(lazyeval::lazy_dots(...),
                                   data = taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))) 
  
  # convert taxon_ids to logical -------------------------------------------------------------------
  is_char <- vapply(selection, is.character, logical(1))
  selection[is_char] <- lapply(selection[is_char], function(x) .data$taxon_data$taxon_ids %in% x)
  
  # convert indexes to logical ---------------------------------------------------------------------
  is_index <- vapply(selection, is.numeric, logical(1))
  selection[is_index] <- lapply(selection[is_index], function(x) 1:nrow(.data$taxon_data) %in% x)
  
  # combine filters --------------------------------------------------------------------------------
  selection <- Reduce(`&`, selection)
  
  # Get taxa of subset -----------------------------------------------------------------------------
  taxa_subset <- unique(c(which(selection),
                          if (subtaxa) {
                            subtaxa(.data, subset = selection, recursive = TRUE, index = TRUE,
                                    include_input = FALSE, simplify = TRUE)
                          },
                          if (supertaxa) {
                            supertaxa(.data, subset = selection, recursive = TRUE, index = TRUE,
                                      na = FALSE, simplify = TRUE, include_input = FALSE)
                          }))
  
  # Reassign taxonless items -----------------------------------------------------------------------
  if (reassign) {
    reassign_one <- function(parents) {
      included_parents <- parents[parents %in% taxa_subset]
      return(.data$taxon_data$taxon_ids[included_parents[1]])
    }
    
    to_reassign <- ! .data$item_data$item_taxon_ids %in% .data$taxon_data$taxon_ids[taxa_subset]
    supertaxa_key <- supertaxa(.data, subset = unique(.data$item_data$item_taxon_ids[to_reassign]),
                              recursive = TRUE, simplify = FALSE, include_input = FALSE, index = TRUE, na = FALSE)
    reassign_key <- vapply(supertaxa_key, reassign_one, character(1))
    .data$item_data[to_reassign, "item_taxon_ids"] <- reassign_key[.data$item_data$item_taxon_ids[to_reassign]]
  }
  
  # Remove taxonless items -------------------------------------------------------------------------
  item_subset <- .data$item_data$item_taxon_ids %in% .data$taxon_data$taxon_ids[taxa_subset]
  if (taxonless) {
    .data$item_data[! item_subset, "item_taxon_ids"] <- as.character(NA)
  } else {
    .data$item_data <- .data$item_data[item_subset, , drop = FALSE]
  }
  
  # Remove filtered taxa ---------------------------------------------------------------------------
  .data$taxa <- .data$taxa[.data$taxon_data$taxon_ids[taxa_subset]]
  .data$taxon_data <- .data$taxon_data[taxa_subset, , drop = FALSE]
  .data$taxon_data[! .data$taxon_data$parent_ids %in% .data$taxon_data$taxon_ids, "parent_ids"] <- as.character(NA)
  
  return(.data)
}



#' Filter items with a list of conditions
#' 
#' Filter items in a \code{\link{classified}} object with a list of conditions. Any column name that
#' appears in \code{item_data(.data)} can be used as if it was a vector on its own. See 
#' \code{\link[dplyr]{filter}} for inspiration and more information.
#' 
#' @param .data \code{\link{classified}}
#' @param ... One or more filtering conditions. This can be one of two things: \describe{ 
#'   \item{\code{integer}}{One or more indexes of \code{item_data}} \item{\code{logical}}{A 
#'   \code{TRUE}/\code{FALSE} vector of length equal to the number of rows in \code{item_data}} } 
#'   Any column name that appears in \code{item_data(.data)} can be used as if it was a vector on 
#'   its own.
#' @param itemless (\code{logical} of length 1) If \code{TRUE}, preserve taxa even if all of their 
#'   items are filtered out. If \code{FALSE}, remove taxa for which all items were filtered out. 
#'   Note that only taxa that are itemless due to this filtering will be removed; there might be 
#'   other taxa without items to begin with that will not be removed.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
filter_items <- function(.data, ..., itemless = TRUE) {
  # non-standard argument evaluation ---------------------------------------------------------------
  selection <- lazyeval::lazy_eval(lazyeval::lazy_dots(...), data = item_data(.data)) 
  
  # convert taxon_ids to logical -------------------------------------------------------------------
  is_char <- vapply(selection, is.character, logical(1))
  if (sum(is_char) > 0) {
    stop("Item filtering with taxon IDs or item IDs (which dont exist yet) is not currently supported. If you want to filter item by taxon IDs, use something like: `item_taxon_ids %in% my_subset`")
  }
  
  # convert indexes to logical ---------------------------------------------------------------------
  is_index <- vapply(selection, is.numeric, logical(1))
  selection[is_index] <- lapply(selection[is_index], function(x) 1:nrow(.data$item_data) %in% x)
  
  # combine filters --------------------------------------------------------------------------------
  selection <- Reduce(`&`, selection)
  
  # Remove items -----------------------------------------------------------------------------------
  itemless_taxa <- unique(.data$item_data$item_taxon_ids[! selection])
  .data$item_data <- .data$item_data[selection, , drop = FALSE]
  
  # Remove itemless taxa ---------------------------------------------------------------------------
  if (! itemless) {
    taxa_to_remove <- .data$taxon_data$taxon_ids %in% itemless_taxa & item_counts(.data) == 0
    .data$taxon_data <- .data$taxon_data[! taxa_to_remove, , drop = FALSE]
    .data$taxon_data[! .data$taxon_data$parent_ids %in% .data$taxon_data$taxon_ids, "parent_ids"] <- as.character(NA)
  }
  
  return(.data)
}

