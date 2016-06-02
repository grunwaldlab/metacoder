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
#'   supplied to \code{collapse_func} to get the item weight. The expression given is evaluated in 
#'   the context of \code{\link{taxon_data}. In other words, any column name that appears in 
#'   \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. If 
#'   \code{item_weight} is also specified, the two weights are multiplied (after \code{taxon_weight}
#'   for each item is calculated).
#' @param item_weight (\code{numeric}) Sampling weights of each item. The expression given is 
#'   evaluated in the context of \code{\link{taxon_data}. In other words, any column name that 
#'   appears in \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. If 
#'   \code{taxon_weight} is also specified, the two weights are multiplied (after 
#'   \code{taxon_weight} for each item is calculated).
#' @param use_supertaxa (\code{logical} of length 1) Affects how the \code{taxon_weight} is used. If
#'   \code{TRUE}, the weights for each taxon in an item's classification are multiplied to get the 
#'   item weight. Otherwise, just the taxonomic level the item is assign to it considered.
#' @param itemless (\code{logical} of length 1) If \code{TRUE}, preserve taxa even if all of their 
#'   items are filtered out. If \code{FALSE}, remove taxa for which all items were filtered out. 
#'   Note that only taxa that are itemless due to this filtering will be removed; there might be 
#'   other taxa without items to begin with that will not be removed.
#' @param collapse_func (\code{function} of length 1) If \code{taxon_weight} is used and 
#'   \code{supertaxa} is \code{TRUE}, the weights for each taxon in an item's classification are 
#'   supplied to \code{collapse_func} to get the item weight. This function should take  numeric
#'   vector and return a single number.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
sample_n_items <- function(.data, size, replace = FALSE, taxon_weight = NULL, item_weight = NULL,
                           use_supertaxa = TRUE, itemless = TRUE, collapse_func = mean) {
  # Calculate taxon component of taxon weights -----------------------------------------------------
  my_taxon_data <- taxon_data(.data)
  taxon_weight <- lazyeval::lazy_eval(lazyeval::lazy(taxon_weight), data = my_taxon_data)
  if (is.null(taxon_weight)) {
    item_taxon_weight <- rep(1, nrow(.data$item_data))
  } else {
    item_index <- match(.data$item_data$item_taxon_ids, .data$taxon_data$taxon_ids)
    my_supertaxa <- supertaxa(.data, recursive = use_supertaxa, simplify = FALSE,
                              include_input = TRUE, index = TRUE, na = FALSE)
    taxon_weight_product <- vapply(my_supertaxa, function(x) collapse_func(taxon_weight[x]), numeric(1))
    item_taxon_weight <- taxon_weight_product[item_index]
  }
  item_taxon_weight <- item_taxon_weight / sum(item_taxon_weight)
  
  # Calculate item component of item weights -------------------------------------------------------
  my_item_data <- item_data(.data)
  item_weight <- lazyeval::lazy_eval(lazyeval::lazy(item_weight), data = my_item_data)
  if (is.null(item_weight)) {
    item_item_weight <- rep(1, nrow(.data$item_data)) 
  }
  item_item_weight <- item_item_weight / sum(item_item_weight)
  
  # Combine item and taxon weight components  ------------------------------------------------------
  item_weight <- item_taxon_weight * item_item_weight
  item_weight <- item_weight / sum(item_weight)
  
  # Sample items -----------------------------------------------------------------------------------
  sampled_rows <- sample.int(nrow(my_item_data), size = size, replace = replace, prob = item_weight)
  filter_items(.data, sampled_rows, itemless = itemless)
}