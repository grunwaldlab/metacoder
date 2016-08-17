#===================================================================================================
#' Recursivly sample a set of taxonomic assignments
#' 
#' Recursivly sample a set of observations with taxonomic assignments and an associated taxonomy.
#' 
#' @param taxmap_data (An object of type \link{taxmap})
#' @param max_counts (\code{numeric}) A named vector that defines that maximum number of
#' observations in for each level specified. The names of the vector specifies that level each number
#' applies to. If more than the maximum number of observations exist for a given taxon, it is randomly
#' subsampled to this number. 
#' @param min_counts (\code{numeric}) A named vector that defines that minimum number of
#' observations in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param max_children (\code{numeric}) A named vector that defines that maximum number of
#' subtaxa per taxon for each level specified. The names of the vector specifies that level each
#' number applies to. If more than the maximum number of subtaxa exist for a given taxon, they
#' are randomly subsampled to this number of subtaxa. 
#' @param min_children (\code{numeric}) A named vector that defines that minimum number of
#' subtaxa in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param obs_filters  (\code{list} of \code{function(observations, id)}) A list of functions that take a data
#' structure containing the information of multiple observations and a taxon id.
#' Returns a object of the same type with some of the observations potentially removed.  
#' @param subtaxa_filters  (\code{list} of \code{function(observations, id)}) A list of functions that take a data
#' structure containing the information of multiple subtaxa IDs and the current taxon id.
#' Returns a object of the same type with some of the subtaxa potentially removed. If a function returns
#' \code{NULL}, then no observations for the current taxon are returned.
#' @param stop_conditions (\code{list} of \code{function(id)}) A list of functions that take the
#' current taxon id. If any of the functions return \code{TRUE}, the observations for the current taxon are 
#' returned rather than looking for observations of subtaxa, stopping the recursion.
#' @param ... Additional parameters are passed to all of the function options.
#' 
#' @return Returns an object of type \code{taxmap}
#' 
#' @examples
#' 
#' \dontrun{
#' #Plot data before subsampling
#' heat_tree(unite_ex_data_3,
#'           node_size = n_obs,
#'           node_color = n_obs,
#'           node_label = n_obs)
#'      
#' # Subsampling
#' subsampled <- taxonomic_sample(unite_ex_data_3,
#'                                max_counts = c("4" = 20, "7" = 5),
#'                                min_counts = c("7" = 3))
#'      
#' # Remove unobserved taxa and plot
#' heat_tree(subset(subsampled, n_obs > 0, unobserved = FALSE),
#'           node_size = n_obs,
#'           node_color = n_obs,
#'           node_label = n_obs)
#' }
#' 
#' @export
taxonomic_sample <- function(taxmap_data,
                             max_counts = c(), min_counts = c(), max_children = c(),
                             min_children = c(), obs_filters = list(), subtaxa_filters = list(),
                             stop_conditions = list(), ...) {
  process_one_tree <- function(root_taxon) {
    # subset for just tree with root
    # tree <- subset(taxmap_data, root_taxon)
    # extract information from `taxmap` (This is a retrofit to use `classfied` objects)
    taxon_ids <- taxmap_data$taxon_data$taxon_ids
    supertaxon_ids <- taxmap_data$taxon_data$supertaxon_ids
    obs_ids <- taxmap_data$obs_data$obs_taxon_ids
    ranks <- n_supertaxa(taxmap_data)
    # Define functions to interact with the taxonomic information ------------------------------------
    get_obs_func <- function(id, ...) which(obs_ids == id)
    get_subtaxa_func <- function(id, ...) taxon_ids[!is.na(supertaxon_ids) & supertaxon_ids == id]
    get_rank_func <- function(id, ...) ranks[taxon_ids == id]
    # recursive sampling -----------------------------------------------------------------------------
    recursive_sample(root_id = root_taxon, get_obs = get_obs_func, get_subtaxa = get_subtaxa_func,
                     get_rank = get_rank_func, cat_obs = unlist, max_counts = max_counts, 
                     min_counts = min_counts, max_children = max_children, min_children = min_children, 
                     obs_filters = obs_filters, subtaxa_filters = subtaxa_filters, 
                     stop_conditions = stop_conditions)
  }
  
  root_taxa <- taxmap_data$taxon_data$taxon_ids[is.na(taxmap_data$taxon_data$supertaxon_ids)]
  obs_indexes <- unlist(lapply(root_taxa, process_one_tree))
  filter_obs(taxmap_data, obs = obs_indexes)
}



#===================================================================================================
#' Recursivly sample observations with a heirarchical classification
#' 
#' Recursivly sample a set of observations with a heirarchical classification.
#' This function takes other functions as arguments and is intended to be used to make other more 
#' user-friendly functions.
#' 
#' @param root_id (\code{character} of length 1) The taxon to sample. By default, the root of the
#' taxonomy used.
#' @param get_obs (\code{function(character)}) A function that returns the observations assigned to the
#' a given taxon. The function's first argument should be the taxon id and it should return a data
#' structure possibly representing multiple observations.
#' @param get_subtaxa (\code{function(character)}) A function that returns the sub taxa for a given
#' taxon. The function's first argument should be the taxon id and it should return a vector of
#' taxon IDs. 
#' @param get_rank (\code{function(character)}) A function that returns the rank of a given taxon
#' id. The function's first argument should be the taxon id and it should return the rank of that
#' taxon.
#' @param cat_obs (\code{function(list)}) A function that takes a list of whatever is returned by
#' \code{get_obs} and concatenates them into a single data structure of the type returned by
#' \code{get_obs}.
#' @param max_counts (\code{numeric}) A named vector that defines that maximum number of
#' observations in for each level specified. The names of the vector specifies that level each number
#' applies to. If more than the maximum number of observations exist for a given taxon, it is randomly
#' subsampled to this number. 
#' @param min_counts (\code{numeric}) A named vector that defines that minimum number of
#' observations in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param max_children (\code{numeric}) A named vector that defines that maximum number of
#' subtaxa per taxon for each level specified. The names of the vector specifies that level each
#' number applies to. If more than the maximum number of subtaxa exist for a given taxon, they
#' are randomly subsampled to this number of subtaxa. 
#' @param min_children (\code{numeric}) A named vector that defines that minimum number of
#' subtaxa in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param obs_filters  (\code{list} of \code{function(observations, id)}) A list of functions that take a data
#' structure containing the information of multiple observations and a taxon id.
#' Returns a object of the same type with some of the observations potentially removed.  
#' @param subtaxa_filters  (\code{list} of \code{function(observations, id)}) A list of functions that take a data
#' structure containing the information of multiple subtaxa IDs and the current taxon id.
#' Returns a object of the same type with some of the subtaxa potentially removed. If a function returns
#' \code{NULL}, then no observations for the current taxon are returned.
#' @param stop_conditions (\code{list} of \code{function(id)}) A list of functions that take the
#' current taxon id. If any of the functions return \code{TRUE}, the observations for the current taxon are 
#' returned rather than looking for observations of subtaxa, stopping the recursion.
#' @param ... Additional parameters are passed to all of the function options.
#' 
#' @seealso \code{\link{taxonomic_sample}}
#' @keywords internal
recursive_sample <- function(root_id, get_obs, get_subtaxa, get_rank = NULL, cat_obs = unlist,
                             max_counts = c(), min_counts = c(), max_children = c(),
                             min_children = c(), obs_filters = list(), subtaxa_filters = list(),
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
    function(observations, id, ...) {
      rank <- as.character(get_rank(id))
      if (rank %in% names(filter_key) && length(observations) > filter_key[rank]) {
        observations <- sample(observations, filter_key[rank])
      }
      return(observations)
    }
  }
  # Make min filter function factory ---------------------------------------------------------------
  min_filter_factory <- function(filter_key) {
    function(observations, id, ...) {
      rank <- as.character(get_rank(id))
      if (rank %in% names(filter_key) && length(observations) < filter_key[rank]) {
        observations <- NULL
      }
      return(observations)
    }
  }
  # Add standard filter functions ------------------------------------------------------------------
  subtaxa_filters <- c(subtaxa_filters, max_filter_factory(max_children), 
                       min_filter_factory(min_children))
  obs_filters <- c(obs_filters, max_filter_factory(max_counts), min_filter_factory(min_counts))
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
    # Get observations for current taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (stop_recursion) {
      observations <- get_obs(id, ...)
    } else {
      observations <- cat_obs(lapply(sub_taxa, recursive_part, depth = depth + 1))
    }
    # Filter observations - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (func in obs_filters) {
      observations <- func(observations, id, ...)
      if (is.null(observations) || length(observations) == 0) break
    }
    return(observations)
  }
  
  recursive_part(root_id, ...)
}

