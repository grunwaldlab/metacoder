#' Get all supertaxa of a taxon
#' 
#' Return the taxon IDs or \code{taxon_data} indexes of all supertaxa (i.e. all taxa the target taxa
#' are a part of) in an object of type \code{taxmap}.
#' 
#' @param obj (\code{taxmap}) The \code{taxmap} object containing taxon information to be 
#'   queried.
#' @param subset (\code{character}) \code{taxon_ids} or indexes of \code{taxon_data} for which
#'   supertaxa will be returned. Default: All taxa in \code{obj} will be used.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the supertaxa one level above the 
#'   target taxa. If \code{TRUE}, return all the supertaxa of every supertaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single 
#'   vector of unique values.
#' @param include_input (\code{logical}) If \code{TRUE}, the input taxa are included in the output
#' @param index (\code{logical}) If \code{TRUE}, return the indexes of supertaxa in 
#'   \code{taxon_data} instead of \code{taxon_ids}
#' @param na (\code{logical}) If \code{TRUE}, return \code{NA} where information is not available.
#'   
#' @return If \code{simplify = FALSE}, then a list of vectors are returned corresponding to the 
#'   \code{subset} argument. If \code{simplify = TRUE}, then unique values are returned in a single 
#'   vector.
#'   
#' @examples 
#' \dontrun{
#' supertaxa(contaminants, subset = 1:10)}
#'   
#' @family taxmap taxonomy functions
#'   
#' @export
supertaxa <- function(obj, subset = NULL, recursive = TRUE,
                      simplify = FALSE, include_input = FALSE, index = FALSE, na = FALSE) {
  # Parse arguments --------------------------------------------------------------------------------
  subset <- format_taxon_subset(obj, subset)
  
  # Get supertaxa ----------------------------------------------------------------------------------
  parent_index <- match(obj$taxon_data$supertaxon_ids, obj$taxon_data$taxon_ids) # precomputing makes it much faster
  recursive_part <- function(taxon) {
    supertaxon <- parent_index[taxon]
    if (recursive) {
      if (is.na(supertaxon)) {
        output <- c(taxon, supertaxon)
      } else {
        output <- c(taxon, recursive_part(supertaxon))
      }
    } else {
      output <- c(taxon, supertaxon)
    }
    return(unname(output))
  }
  output <- lapply(subset, recursive_part)
  
  # Remove query taxa from output ------------------------------------------------------------------
  if (! include_input) {
    output <- lapply(output, `[`, -1)
  }
  
  # Remove NAs from output -------------------------------------------------------------------------
  if (! na) {
    output <- lapply(output, function(x) x[!is.na(x)])
  }
  
  # Convert to taxon_ids ---------------------------------------------------------------------------
  if (! index) {
    output <- lapply(output, function(x) obj$taxon_data$taxon_ids[x])
  }
  
  # Reduce dimensionality --------------------------------------------------------------------------
  if (simplify) {
    output <- unique(unname(unlist(output)))
  }
  
  return(output)
}


#' Get root taxa
#' 
#' Return the root taxa for a \code{\link{taxmap}} object. Can also be used to get the roots of
#' a subset of taxa.
#' 
#' @param obj (\code{taxmap}) The \code{taxmap} object containing taxon information to be
#'   queried.
#' @param subset (\code{character}) Taxon IDs for which supertaxa will be returned. Default: All
#'   taxon in \code{obj} will be used.
#' @param index (\code{logical}) If \code{TRUE}, return the indexes of roots in 
#'   \code{taxon_data} instead of \code{taxon_ids}
#'   
#' @return \code{character}
#'   
#' @family taxmap taxonomy functions
#'   
#' @export
roots <- function(obj, subset = NULL, index = FALSE) {
  # Parse arguments --------------------------------------------------------------------------------
  subset <- format_taxon_subset(obj, subset)
  
  # Get roots --------------------------------------------------------------------------------------
  parents <- supertaxa(obj, subset = subset, recursive = TRUE, include_input = TRUE, index = TRUE, na = FALSE)
  is_global_root <- vapply(parents, length, numeric(1)) == 1
  if (missing(subset)) {
    is_root <- is_global_root
  } else {
    is_root <- is_global_root | vapply(parents, function(x) ! any(x[-1] %in% subset), logical(1))
  }
  output <- unname(subset[is_root])
  
  # Convert to taxon_ids ---------------------------------------------------------------------------
  if (! index) {
    output <- obj$taxon_data$taxon_ids[output]
  }
  
  return(output)
}


#' Get all subtaxa of a taxon
#' 
#' Return the taxon IDs or \code{taxon_data} indexes of all subtaxa in an object of type \code{taxmap}
#' 
#' @param obj (\code{taxmap}) The \code{taxmap} object containing taxon information to be
#'   queried.
#' @param subset (\code{character}) \code{taxon_ids} or indexes of \code{taxon_data} for which
#'   supertaxa will be returned. Default: All taxa in \code{obj} will be used.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the subtaxa one level below the 
#'   target taxa. If \code{TRUE}, return all the subtaxa of every subtaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single 
#'   vector of unique values.
#' @param include_input (\code{logical}) If \code{TRUE}, the input taxa are included in the output
#' @param index (\code{logical}) If \code{TRUE}, return the indexes of supertaxa in 
#'   \code{taxon_data} instead of \code{taxon_ids}
#'   
#' @return If \code{simplify = FALSE}, then a list of vectors are returned corresponding to the
#'   \code{target} argument. If \code{simplify = TRUE}, then the unique values are returned in a
#'   single vector.
#'
#' @examples 
#' \dontrun{
#' subtaxa(contaminants, subset = 1:10)}
#' 
#' @family taxmap taxonomy functions
#' 
#' @export
subtaxa <- function(obj, subset = NULL, recursive = TRUE,
                    simplify = FALSE, include_input = FALSE, index = FALSE) {
  # Parse arguments --------------------------------------------------------------------------------
  subset <- format_taxon_subset(obj, subset)
  
  # Return empty list if `subset` has no values
  if (length(subset) == 0) {
    if (simplify) {
      return(character(0))
    } else {
      return(list())
    }
  }
  
  # Get subtaxa ------------------------------------------------------------------------------------
  parent_index <- match(obj$taxon_data$supertaxon_ids, obj$taxon_data$taxon_ids) 

  get_children <- function(taxon) {
    which(parent_index == taxon)
  }
  
  recursive_part <- function(taxon) {
    # Get immediate children of current taxon
    children <- get_children(taxon)
    # Run this function on them to get their output
    child_output <- lapply(children, recursive_part)
    child_output <- stats::setNames(unlist(child_output, recursive = FALSE),
                                    unlist(lapply(child_output, names)))
    # Get all subtaxa from the names of the child output
    child_taxa <- c(taxon, as.numeric(names(child_output)))
    # Combine the child output with the subtaxa for the current taxon
    output <- stats::setNames(c(list(child_taxa), child_output),
                              c(taxon, names(child_output)))
    return(output)
  }
  
  if (recursive) {
    starting_taxa <- roots(obj, subset = subset, index = TRUE)
    output <- stats::setNames(unlist(lapply(starting_taxa, recursive_part), recursive = FALSE)[as.character(subset)],
                              names(subset))
  } else {
    output <- lapply(subset, function(x) c(x, get_children(x)))
  }
  
  # Remove query taxa from output ------------------------------------------------------------------
  if (! include_input) {
    output <- lapply(output, `[`, -1)
  }
  
  # Convert to taxon_ids ---------------------------------------------------------------------------
  if (! index) {
    output <- lapply(output, function(x) obj$taxon_data$taxon_ids[x])
  }
  
  # Reduce dimensionality --------------------------------------------------------------------------
  if (simplify) {
    output <- unique(unname(unlist(output)))
  }
  
  return(output)
}


#' Get observations associated with taxa
#'
#' Given one or more taxa IDs and a \code{\link{taxmap}} object, return the observation indexes
#' (e.g. sequence information) associated with each taxon.
#'
#' @param obj (\code{taxmap})
#' The \code{taxmap} object containing taxon information to be queried.
#' @param subset (\code{character}) \code{taxon_ids} or indexes of \code{taxon_data} for which
#'   supertaxa will be returned. Default: All taxa in \code{obj} will be used.
#' @param recursive (\code{logical})
#' If \code{FALSE}, only return the observation assigned to the specified input taxa, not subtaxa.
#' If \code{TRUE}, return all the observations of every subtaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique observation indexes.
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of observation indexes are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the observation indexes
#'   for all \code{target} taxa are returned in a single vector.
#'
#' @family taxmap taxonomy functions
#'
#' @export
obs <- function(obj, subset = NULL, recursive = TRUE, simplify = FALSE) {
  # Parse arguments --------------------------------------------------------------------------------
  subset <- format_taxon_subset(obj, subset)
  
  # Get observations of taxa ------------------------------------------------------------------------------
  my_subtaxa <- subtaxa(obj, subset = subset, recursive = recursive, include_input = TRUE, index = TRUE)
  unique_subtaxa <- unique(unlist(my_subtaxa))
  obs_taxon_index <- match(obj$obs_data$obs_taxon_ids, obj$taxon_data$taxon_ids)
  obs_key <- split(seq_along(obj$obs_data$obs_taxon_ids), obs_taxon_index)
  output <- stats::setNames(lapply(my_subtaxa, function(x) unname(unlist(obs_key[as.character(x)]))),
                            names(subset))
  is_null <- vapply(output, is.null, logical(1))
  output[is_null] <- lapply(1:sum(is_null), function(x) numeric(0))
  
  # Reduce dimensionality --------------------------------------------------------------------------
  if (simplify) {
    output <- unique(unname(unlist(output)))
  }
  
  return(output)
}
