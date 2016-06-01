#' Get all supertaxa of a taxon
#'
#' Return the taxon IDs of all supertaxa (i.e. all taxa the target taxa are a part of) in an
#' object of type \code{classified}
#'
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which supertaxa will be returned.
#' Default: All taxon in \code{obj} will be used.
#' @param recursive (\code{logical})
#' If \code{FALSE}, only return the supertaxa one level above the
#'   target taxa. If \code{TRUE}, return all the supertaxa of every supertaxa, etc.
#' @param simplify (\code{logical})
#' If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs
#' @param include_input (\code{logical})
#' If \code{TRUE}, the input taxa are included in the output
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
#'   
#' @export
supertaxa <- function(obj, subset = obj$taxon_data$taxon_ids, recursive = TRUE,
                      simplify = FALSE, include_input = FALSE) {
  recursive_part <- function(taxon) {
    supertaxon <- obj$taxon_data$parent_ids[taxon == obj$taxon_data$taxon_ids]
    if (recursive) {
      if (is.na(supertaxon)) {
        output <- taxon
      } else {
        output <- c(taxon, recursive_part(supertaxon))
      }
    } else {
      output <- c(taxon, supertaxon)
    }
    return(unname(output))
  }
  
  subset <- format_taxon_subset(obj, subset)
  supertaxa <- stats::setNames(lapply(subset, recursive_part), subset)
  if (!include_input) {
    supertaxa <- lapply(supertaxa, `[`, -1)
  }
  # Reduce dimensionality if specified
  if (simplify) {
    supertaxa <- unname(unlist(supertaxa))
  }
  return(supertaxa)
}


#' Get all subtaxa of a taxon
#'
#' Return the taxon IDs of all subtaxa in an object of type \code{classified}
#'
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which subtaxa will be returned.
#' @param recursive (\code{logical})
#' If \code{FALSE}, only return the subtaxa one level above the
#'   target taxa. If \code{TRUE}, return all the subtaxa of every subtaxa, etc.
#' @param simplify (\code{logical})
#' If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs
#' @param include_input (\code{logical})
#' If \code{TRUE}, the input taxa are included in the output
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
#'
#' @export
subtaxa <- function(obj, subset = obj$taxon_data$taxon_ids, recursive = TRUE,
                    simplify = FALSE, include_input = FALSE) {
  get_children <- function(taxon) {
    unname(obj$taxon_data$taxon_ids[obj$taxon_data$parent_ids == taxon & ! is.na(obj$taxon_data$parent_ids)])
  }
  
  recursive_part <- function(taxon) {
    # Get immediate children of current taxon
    children <- get_children(taxon)
    # Run this function on them to get their output
    child_output <- lapply(children, recursive_part)
    child_output <- stats::setNames(unlist(child_output, recursive = FALSE),
                                    unlist(lapply(child_output, names)))
    # Get all subtaxa from the names of the child output
    if (include_input) {
      child_taxa <- c(taxon, as.numeric(names(child_output)))
    } else {
      child_taxa <- as.numeric(names(child_output))
      if (is.null(child_taxa)) {
        child_taxa <- numeric(0)
      }
    }
    # Combine the child output with the subtaxa for the current taxon
    output <- stats::setNames(c(list(child_taxa), child_output),
                              c(taxon, names(child_output)))
    return(output)
  }
  
  
  subset <- format_taxon_subset(obj, subset)  # Get output content
  if (recursive) {
    starting_taxa <- roots(obj, subset)
    output <- unlist(lapply(starting_taxa, recursive_part), recursive = FALSE)[as.character(subset)]
  } else {
    output <- stats::setNames(lapply(subset, get_children), subset)
    if (include_input) {
      output <- mapply(function(x, n) c(n, x), output, names(output), SIMPLIFY = FALSE)
    }
  }
  # Reduce dimensionality if specified
  if (simplify) {
    output <- unname(unique(unlist(output)))
  }
  return(output)
}


#' Get items associated with taxa
#'
#' Given one or more taxa IDs and a \code{\link{classified}} object, return the items
#' (e.g. sequence information) associated with each taxon.
#'
#' @param obj (\code{classified})
#' The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which items will be returned.
#' @param recursive (\code{logical})
#' If \code{FALSE}, only return the item assigned to the specified input taxa, not subtaxa.
#' If \code{TRUE}, return all the items of every subtaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique item indexes.
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of item indexes are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the item indexes
#'   for all \code{target} taxa are returned in a single vector.
#'
#' @export
items <- function(obj, subset = obj$taxon_data$taxon_ids, recursive = TRUE, simplify = FALSE) {
  # Get output content
  my_subtaxa <- subtaxa(obj, subset, recursive = recursive, include_input = TRUE)
  unique_subtaxa <- unique(unlist(my_subtaxa))
  item_key <- stats::setNames(lapply(unique_subtaxa, function(x) which(x == obj$item_data$item_taxon_ids)),
                              unique_subtaxa)
  output <- lapply(my_subtaxa, function(x) unname(unlist(item_key[as.character(x)])))
  # Reduce dimensionality if specified
  if (simplify) {
    output <- unname(unique(unlist(output)))
  }
  return(output)
}


#' Get root taxa
#' 
#' Return the root taxa for a \code{\link{classified}} object.
#' Can also be used to get the roots of a subset of taxa.
#' 
#' @param obj (\code{classified}) The \code{classified} object containing taxon information to be queried.
#' @param subset (\code{character})
#' Taxon IDs for which supertaxa will be returned.
#' Default: All taxon in \code{obj} will be used.
#' 
#' @return \code{character}
#'  
#' @export
roots <- function(obj, subset = obj$taxon_data$taxon_ids) {
  parents <- supertaxa(obj, subset = subset, include_input = TRUE)
  is_global_root <- vapply(parents, function(x) length(x) == 1, logical(1))
  if (missing(subset)) {
    is_root <- is_global_root
  } else {
    subset <- format_taxon_subset(obj, subset)
    is_root <- is_global_root | vapply(parents, function(x) ! any(x[-1] %in% subset), logical(1))
  }
  subset[is_root]
}
