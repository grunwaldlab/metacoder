#===================================================================================================
#' Recursivly sample a set of taxonomic assignments
#' 
#' Recursivly sample a set of items with taxonomic assignments and an associated taxonomy.
#' 
#' @param classified_data (An object of type \link{classified})
#' @param max_counts (\code{numeric}) A named vector that defines that maximum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. If more than the maximum number of items exist for a given taxon, it is randomly
#' subsampled to this number. 
#' @param min_counts (\code{numeric}) A named vector that defines that minimum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param max_children (\code{numeric}) A named vector that defines that maximum number of
#' subtaxa per taxon for each level specified. The names of the vector specifies that level each
#' number applies to. If more than the maximum number of subtaxa exist for a given taxon, they
#' are randomly subsampled to this number of subtaxa. 
#' @param min_children (\code{numeric}) A named vector that defines that minimum number of
#' subtaxa in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param item_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple items and a taxon id.
#' Returns a object of the same type with some of the items potentially removed.  
#' @param subtaxa_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple subtaxa IDs and the current taxon id.
#' Returns a object of the same type with some of the subtaxa potentially removed. If a function returns
#' \code{NULL}, then no items for the current taxon are returned.
#' @param stop_conditions (\code{list} of \code{function(id)}) A list of functions that take the
#' current taxon id. If any of the functions return \code{TRUE}, the items for the current taxon are 
#' returned rather than looking for items of subtaxa, stopping the recursion.
#' @param ... Additional parameters are passed to all of the function options.
#' 
#' @return Returns an object of type \code{classified}
#' 
#' @examples
#' 
#' \dontrun{
#' #Plot data before subsampling
#' plot(unite_ex_data_3,
#'      vertex_size = item_count,
#'      vertex_color = item_count,
#'      vertex_label = item_count)
#'      
#' # Subsampling
#' subsampled <- taxonomic_sample(unite_ex_data_3,
#'                                max_counts = c("4" = 20, "7" = 5),
#'                                min_counts = c("7" = 3))
#'      
#' # Remove itemless taxa and plot
#' plot(subset(subsampled, item_count > 0, itemless = FALSE),
#'      vertex_size = item_count,
#'      vertex_color = item_count,
#'      vertex_label = item_count)
#' }
#' 
#' @export
taxonomic_sample <- function(classified_data,
                             max_counts = c(), min_counts = c(), max_children = c(),
                             min_children = c(), item_filters = list(), subtaxa_filters = list(),
                             stop_conditions = list(), ...) {
  process_one_tree <- function(root_taxon) {
    # subset for just tree with root
    # tree <- subset(classified_data, root_taxon)
    # extract information from `classified` (This is a retrofit to use `classfied` objects)
    taxon_ids <- classified_data$taxon_data$taxon_ids
    parent_ids <- classified_data$taxon_data$parent_ids
    item_ids <- classified_data$item_data$item_taxon_ids
    ranks <- taxon_ranks(classified_data)
    # Define functions to interact with the taxonomic information ------------------------------------
    get_items_func <- function(id, ...) which(item_ids == id)
    get_subtaxa_func <- function(id, ...) taxon_ids[!is.na(parent_ids) & parent_ids == id]
    get_rank_func <- function(id, ...) ranks[taxon_ids == id]
    # recursive sampling -----------------------------------------------------------------------------
    recursive_sample(root_id = root_taxon, get_items = get_items_func, get_subtaxa = get_subtaxa_func,
                     get_rank = get_rank_func, cat_items = unlist, max_counts = max_counts, 
                     min_counts = min_counts, max_children = max_children, min_children = min_children, 
                     item_filters = item_filters, subtaxa_filters = subtaxa_filters, 
                     stop_conditions = stop_conditions)
  }
  
  root_taxa <- classified_data$taxon_data$taxon_ids[is.na(classified_data$taxon_data$parent_ids)]
  item_indexes <- unlist(lapply(root_taxa, process_one_tree))
  filter_items(classified_data, item = item_indexes)
}



#===================================================================================================
#' Recursivly sample items with a heirarchical classification
#' 
#' Recursivly sample a set of items with a heirarchical classification.
#' This function takes other functions as arguments and is intended to be used to make other more 
#' user-friendly functions.
#' 
#' @param root_id (\code{character} of length 1) The taxon to sample. By default, the root of the
#' taxonomy used.
#' @param get_items (\code{function(character)}) A function that returns the items assigned to the
#' a given taxon. The function's first argument should be the taxon id and it should return a data
#' structure possibly representing multiple items.
#' @param get_subtaxa (\code{function(character)}) A function that returns the sub taxa for a given
#' taxon. The function's first argument should be the taxon id and it should return a vector of
#' taxon IDs. 
#' @param get_rank (\code{function(character)}) A function that returns the rank of a given taxon
#' id. The function's first argument should be the taxon id and it should return the rank of that
#' taxon.
#' @param cat_items (\code{function(list)}) A function that takes a list of whatever is returned by
#' \code{get_items} and concatenates them into a single data structure of the type returned by
#' \code{get_items}.
#' @param max_counts (\code{numeric}) A named vector that defines that maximum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. If more than the maximum number of items exist for a given taxon, it is randomly
#' subsampled to this number. 
#' @param min_counts (\code{numeric}) A named vector that defines that minimum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param max_children (\code{numeric}) A named vector that defines that maximum number of
#' subtaxa per taxon for each level specified. The names of the vector specifies that level each
#' number applies to. If more than the maximum number of subtaxa exist for a given taxon, they
#' are randomly subsampled to this number of subtaxa. 
#' @param min_children (\code{numeric}) A named vector that defines that minimum number of
#' subtaxa in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param item_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple items and a taxon id.
#' Returns a object of the same type with some of the items potentially removed.  
#' @param subtaxa_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple subtaxa IDs and the current taxon id.
#' Returns a object of the same type with some of the subtaxa potentially removed. If a function returns
#' \code{NULL}, then no items for the current taxon are returned.
#' @param stop_conditions (\code{list} of \code{function(id)}) A list of functions that take the
#' current taxon id. If any of the functions return \code{TRUE}, the items for the current taxon are 
#' returned rather than looking for items of subtaxa, stopping the recursion.
#' @param ... Additional parameters are passed to all of the function options.
#' 
#' @seealso \code{\link{taxonomic_sample}}
#' @keywords internal
recursive_sample <- function(root_id, get_items, get_subtaxa, get_rank = NULL, cat_items = unlist,
                             max_counts = c(), min_counts = c(), max_children = c(),
                             min_children = c(), item_filters = list(), subtaxa_filters = list(),
                             stop_conditions = list(), ...) {
  # Parse options ----------------------------------------------------------------------------------
  validate_filter_options <- function(filter) {
    if (length(get(filter)) > 0 && is.null(names(get(filter)))) {
      if (!is.null(get_rank)) stop(paste0("`", filter, "` must be named if `get_rank` is defined."))
      return(stats::setNames(get(filter), as.character(seq_along(get(filter)))))
    }
    return(get(filter))
  }
  max_counts <- validate_filter_options("max_counts")
  min_counts <- validate_filter_options("min_counts")
  max_children <- validate_filter_options("max_children")
  min_children <- validate_filter_options("min_children")
  # Set default get_rank function ------------------------------------------------------------------
  if (is.null(get_rank)) {
    get_rank <- function(id, ...) {
      return(get("depth"))
    }
  }
  # Make max filter function factory ---------------------------------------------------------------
  max_filter_factory <- function(filter_key) {
    function(items, id, ...) {
      rank <- as.character(get_rank(id))
      if (rank %in% names(filter_key) && length(items) > filter_key[rank]) {
        items <- sample(items, filter_key[rank])
      }
      return(items)
    }
  }
  # Make min filter function factory ---------------------------------------------------------------
  min_filter_factory <- function(filter_key) {
    function(items, id, ...) {
      rank <- as.character(get_rank(id))
      if (rank %in% names(filter_key) && length(items) < filter_key[rank]) {
        items <- NULL
      }
      return(items)
    }
  }
  # Add standard filter functions ------------------------------------------------------------------
  subtaxa_filters <- c(subtaxa_filters, max_filter_factory(max_children), 
                       min_filter_factory(min_children))
  item_filters <- c(item_filters, max_filter_factory(max_counts), min_filter_factory(min_counts))
  # Recursivly sample taxon ------------------------------------------------------------------------
  recursive_part <- function(id, depth = 1, ...) {
    # Determine if to stop search  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    stop_recursion = any(vapply(stop_conditions, function(func) func(id, ...), logical(1)))
    # Get and filter subtaxa of current taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!stop_recursion) {
      sub_taxa <- get_subtaxa(id)
      for (func in subtaxa_filters) {
        sub_taxa <- func(sub_taxa, id, ...)
        if (is.null(sub_taxa)) return(NULL)
        if (length(sub_taxa) == 0) {
          stop_recursion = TRUE
          break
        }
      }
    }
    # Get items for current taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (stop_recursion) {
      items <- get_items(id, ...)
    } else {
      items <- cat_items(lapply(sub_taxa, recursive_part, depth = depth + 1))
    }
    # Filter items - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (func in item_filters) {
      items <- func(items, id, ...)
      if (is.null(items) || length(items) == 0) break
    }
    return(items)
  }
  
  recursive_part(root_id, ...)
}

