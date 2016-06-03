#' Sample n items from \code{\link{classified}}
#' 
#' Randomly sample some number of items from a \code{\link{classified}} object. Weights can be 
#' specified for items or the taxa they are classified by.
#' See \link[dplyr]{sample_n} for the inspiration for this function.
#' 
#' @param .data (\code{\link{classified}}) The object to sample from.
#' @param size (\code{numeric} of length 1) The number of items to sample.
#' @param replace (\code{logical} of length 1) If \code{TRUE}, sample with replacement.
#' @param taxon_weight (\code{numeric}) Non-negative sampling weights of each taxon. If 
#'   \code{use_supertaxa} is \code{TRUE}, the weights for each taxon in an item's classification are
#'   supplied to \code{collapse_func} to get the item weight. The expression given is evaluated in 
#'   the context of \code{\link{taxon_data}}. In other words, any column name that appears in 
#'   \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. If 
#'   \code{item_weight} is also specified, the two weights are multiplied (after \code{taxon_weight}
#'   for each item is calculated).
#' @param item_weight (\code{numeric}) Sampling weights of each item. The expression given is 
#'   evaluated in the context of \code{\link{item_data}}. In other words, any column name that 
#'   appears in \code{\link{item_data}(.data)} can be used as if it was a vector on its own. If 
#'   \code{taxon_weight} is also specified, the two weights are multiplied (after 
#'   \code{taxon_weight} for each item is calculated).
#' @param use_supertaxa (\code{logical} of length 1) Affects how the \code{taxon_weight} is used. If
#'   \code{TRUE}, the weights for each taxon in an item's classification are multiplied to get the 
#'   item weight. Otherwise, just the taxonomic level the item is assign to it considered.
#' @param collapse_func (\code{function} of length 1) If \code{taxon_weight} option is used and 
#'   \code{supertaxa} is \code{TRUE}, the weights for each taxon in an item's classification are 
#'   supplied to \code{collapse_func} to get the item weight. This function should take  numeric 
#'   vector and return a single number.
#' @param ... Additional options are passed to \code{\link{filter_items}}.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
sample_n_items <- function(.data, size, replace = FALSE, taxon_weight = NULL, item_weight = NULL,
                           use_supertaxa = TRUE, collapse_func = mean, ...) {
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
    item_weight <- rep(1, nrow(.data$item_data)) 
  }
  item_weight <- item_weight / sum(item_weight)
  
  # Combine item and taxon weight components  ------------------------------------------------------
  combine_func <- prod
  weight <- mapply(item_taxon_weight, item_weight, FUN = function(x, y) combine_func(c(x,y)))
  weight <- weight / sum(weight)
  
  # Sample items -----------------------------------------------------------------------------------
  sampled_rows <- sample.int(nrow(my_item_data), size = size, replace = replace, prob = weight)
  filter_items(.data, sampled_rows, ...)
}

#' Sample a proportion of items from \code{\link{classified}}
#' 
#' Randomly sample some propoortion of items from a \code{\link{classified}} object. Weights can be 
#' specified for items or the taxa they are classified by.
#' See \link[dplyr]{sample_frac} for the inspiration for this function.
#' 
#' @inheritParams sample_n_items
#' @param size (\code{numeric} of length 1) The proportion of items to sample.
#' 
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
sample_frac_items <- function(.data, size = 1, replace = FALSE, taxon_weight = NULL, item_weight = NULL,
                           use_supertaxa = TRUE, collapse_func = mean, ...) {
  sample_n_items(.data = .data, size = size * nrow(.data$item_data), replace = replace,
                 taxon_weight = taxon_weight, item_weight = item_weight,
                 use_supertaxa = use_supertaxa, collapse_func = collapse_func, ...)
}


#' Sample n taxa from \code{\link{classified}}
#' 
#' Randomly sample some number of taxa from a \code{\link{classified}} object. Weights can be 
#' specified for taxa or the items assigned to them.
#' See \link[dplyr]{sample_n} for the inspiration for this function.
#' 
#' @param .data (\code{\link{classified}}) The object to sample from.
#' @param size (\code{numeric} of length 1) The number of taxa to sample.
#' @param taxon_weight (\code{numeric}) Non-negative sampling weights of each taxon. The expression 
#'   given is evaluated in the context of \code{\link{taxon_data}}. In other words, any column name 
#'   that appears in \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. 
#'   If \code{item_weight} is also specified, the two weights are multiplied (after 
#'   \code{item_weight} for each taxon is calculated).
#' @param item_weight (\code{numeric}) Sampling weights of each item. The weights for each item 
#'   assigned to a given taxon are supplied to \code{collapse_func} to get the taxon weight. If 
#'   \code{use_subtaxa} is \code{TRUE} then the items assigned to every subtaxa are also used. The 
#'   expression given is evaluated in the context of \code{\link{item_data}}. In other words, any 
#'   column name that appears in \code{\link{item_data}(.data)} can be used as if it was a vector on
#'   its own. If \code{taxon_weight} is also specified, the two weights are multiplied (after 
#'   \code{item_weight} for each item is calculated).
#' @param use_subtaxa (\code{logical} of length 1) Affects how the \code{item_weight} option is
#'   used. If \code{TRUE}, the weights for each taxon in an item's classification are multiplied to
#'   get the item weight. Otherwise, just the taxonomic level the item is assign to it considered.
#' @param collapse_func (\code{function} of length 1) If \code{taxon_weight} is used and 
#'   \code{supertaxa} is \code{TRUE}, the weights for each taxon in an item's classification are 
#'   supplied to \code{collapse_func} to get the item weight. This function should take  numeric 
#'   vector and return a single number.
#' @param ... Additional options are passed to \code{\link{filter_taxa}}.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
sample_n_taxa <- function(.data, size, taxon_weight = NULL, item_weight = NULL,
                          use_subtaxa = TRUE, collapse_func = mean, ...) {
  # Calculate item component of taxon weights ------------------------------------------------------
  my_item_data <- item_data(.data)
  item_weight <- lazyeval::lazy_eval(lazyeval::lazy(item_weight), data = my_item_data)
  if (is.null(item_weight)) {
    taxon_item_weight <- rep(1, nrow(.data$taxon_data))
  } else {
    my_items <- items(.data, recursive = use_subtaxa, simplify = FALSE)
    taxon_item_weight <- vapply(my_items, function(x) collapse_func(item_weight[x]), numeric(1))
  }
  taxon_item_weight <- taxon_item_weight / sum(taxon_item_weight)
  
  # Calculate taxon component of taxon weights -------------------------------------------------------
  my_taxon_data <- taxon_data(.data)
  taxon_weight <- lazyeval::lazy_eval(lazyeval::lazy(taxon_weight), data = my_taxon_data)
  if (is.null(taxon_weight)) {
    taxon_weight <- rep(1, nrow(.data$taxon_data)) 
  }
  taxon_weight <- taxon_weight / sum(taxon_weight)
  
  # Combine item and taxon weight components  ------------------------------------------------------
  combine_func <- prod
  weight <- mapply(taxon_weight, taxon_item_weight, FUN = function(x, y) combine_func(c(x,y)))
  weight <- weight / sum(weight)
  
  # Sample items -----------------------------------------------------------------------------------
  sampled_rows <- sample.int(nrow(.data$taxon_data), size = size, replace = FALSE, prob = weight)
  filter_taxa(.data, sampled_rows, ...)
}


#' Sample a proportion of taxa from \code{\link{classified}}
#' 
#' Randomly sample some propoortion of taxa from a \code{\link{classified}} object. Weights can be 
#' specified for taxa or the items assigned to them. See \link[dplyr]{sample_frac} for the
#' inspiration for this function.
#' 
#' @inheritParams sample_n_taxa
#' @param size (\code{numeric} of length 1) The proportion of taxa to sample.
#'   
#' @return An object of type \code{\link{classified}}
#'   
#' @family dplyr-like functions
#'   
#' @export
sample_frac_taxa <- function(.data, size = 1, taxon_weight = NULL, item_weight = NULL,
                              use_subtaxa = TRUE, collapse_func = mean, ...) {
  sample_n_taxa(.data = .data, size = size * nrow(.data$item_data),
                 taxon_weight = taxon_weight, item_weight = item_weight,
                use_subtaxa = use_subtaxa, collapse_func = collapse_func, ...)
}

