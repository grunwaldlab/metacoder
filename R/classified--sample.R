#' Sample n items from \code{\link{classified}}
#' 
#' Randomly sample some number of items from a \code{\link{classified}} object. Weighs can be
#' specified for items or the taxa they are classified by.
#' 
#' See \link[dplyr]{sample_n} for the inspiration for this function.
#' 
#' @param .data (\code{\link{classified}})
#' @param size (\code{numeric} of length 1) The number of items to sample.
#' @param replace (\code{logical} of length 1) If \code{TRUE}, sample with replacement.
#' @param taxon_weight (\code{numeric}) Non-negative sampling weights of each taxon. If
#'   \code{supertaxa} is \code{TRUE}, the weights for each taxon in an item's classification are
#'   multiplied to get the item weight. The expression given is evaluated in the context of
#'   \code{\link{taxon_data}. In other words, any column name that appears in
#'   \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. If
#'   \code{item_weight} is also specified, the two weights are multiplied (after \code{taxon_weight}
#'   for each item is calculated).
#' @param item_weight (\code{numeric}) Sampling weights of each item. The expression given is
#'   evaluated in the context of \code{\link{taxon_data}. In other words, any column name that
#'   appears in \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. If
#'   \code{taxon_weight} is also specified, the two weights are multiplied (after
#'   \code{taxon_weight} for each item is calculated).
#' @param supertaxa (\code{logical} of length 1) Affects how the \code{taxon_weight} is used. If
#'   \code{TRUE}, the weights for each taxon in an item's classification are multiplied to get the
#'   item weight. Otherwise, just the taxonomic level the item is assign to it considered.
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
sample_n_items <- function(.data, size, replace = FALSE, taxon_weight = NULL,
                           item_weight = NULL, supertaxa = TRUE) {
  # Calculate taxon component of item weights ------------------------------------------------------
  
  
  # Calculate item component of item weights -------------------------------------------------------
  
  
  # Combine item and taxon weight components  ------------------------------------------------------
  
  
  # Sample items -----------------------------------------------------------------------------------
}