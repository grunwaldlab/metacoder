

#' Filter \code{\link{classified}} on a list of conditions
#'
#' Create a subset of a \code{\link{classified}} object.
#' Can filter based on columns in \code{obj$taxon_data} or \code{obj$taxon_data}.
#'
#' @param obj \code{\link{classified}}
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
filter.classified <- function(obj, ..., subtaxa = TRUE, supertaxa = FALSE, itemless = TRUE, taxonless = FALSE) {
  filter_(obj, .dots = lazyeval::lazy_dots(...))
}


#' @export
#' @rdname filter
filter_.classified <- function(obj, taxon = taxon_ids(x), item = seq_along(x$item_taxon_ids),
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
  inluded_items <- intersect(which(x$item_taxon_ids %in% new_taxa),
                             parsed_item)
  
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
