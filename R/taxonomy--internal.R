#===================================================================================================
#' Get all supertaxa of a taxon
#'
#' Given one or more taxa IDs and the edge list defining the taxonomy, return the taxon IDs of all
#' supertaxa (i.e. all taxa the target taxa are a part of).
#'
#' @param targets (\code{character}) Taxon IDs for which supertaxa will be returned.
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the supertaxa one level above the
#'   target taxa. If \code{TRUE}, return all the supertaxa of every supertaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs
#' @param include_target (\code{logical}) If \code{TRUE}, the target taxa are included in the output
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
get_supertaxa <- function(targets, taxa, parents, recursive = TRUE, simplify = FALSE,
                          include_target = FALSE) {
  # Argument validataion ---------------------------------------------------------------------------
  if (length(targets) == 0) stop("Argument 'targets' has 0 length")
  if (!all(targets %in% taxa)) stop("All 'targets' taxon IDs not found in 'taxa' taxon IDs")
  validate_edge_list(taxa, parents)
  parents[!(parents %in% taxa)] <- NA

  # Recursive function for one target --------------------------------------------------------------
  get_one <- function(target) {
    supertaxon <- parents[taxa == target]
    if (recursive) {
      if (is.na(supertaxon)) {
        return(target)
      } else {
        return(c(target, get_one(supertaxon)))
      }
    } else {
      return(c(target, supertaxon))
    }
  }

  # Apply function to all targets ------------------------------------------------------------------
  supertaxa <- lapply(targets, get_one)
  if (!include_target) supertaxa <- lapply(supertaxa, `[`, -1)
  if (simplify) supertaxa <- unlist(supertaxa)
  return(supertaxa)
}





#===================================================================================================
#' Get all subtaxa of a taxon
#'
#' Given one or more taxa IDs and the edge list defining the taxonomy, return the taxon IDs of all
#' subtaxa
#'
#' @param targets (\code{character}) Taxon IDs for which subtaxa will be returned.
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the subtaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the subtaxa one level above the
#'   target taxa. If \code{TRUE}, return all the subtaxa of every subtaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
#'
#' @export
get_subtaxa <- function(targets, taxa, parents, recursive = TRUE, simplify = FALSE) {
  # Argument validataion ---------------------------------------------------------------------------
  if (length(targets) == 0) stop("Argument 'targets' has 0 length")
  if (!all(targets %in% taxa)) stop("All 'targets' taxon IDs not found in 'taxa' taxon IDs")
  validate_edge_list(taxa, parents)
  parents[!(parents %in% taxa)] <- NA

  # Get level of each input to determin processing order -------------------------------------------
  levels <- vapply(get_supertaxa(targets, taxa, parents), length, numeric(1))
  original_targets <- targets
  targets <- targets[order(levels, decreasing = TRUE)]

  # Recursive function for one target --------------------------------------------------------------
  get_one <- function(target, output) {
    if (target %in% names(output)) { #if this taxa has already been processed
      return(output[[target]])
    } else {
      subtaxa <- taxa[parents == target]
      subtaxa <- subtaxa[!is.na(subtaxa)]
      if (recursive) {
        if (length(subtaxa) == 0) {
          return(subtaxa)
        } else {
          subsubtaxa <- unlist(lapply(subtaxa, get_one, output = output))
          return(c(subtaxa, subsubtaxa))
        }
      } else {
        return(subtaxa)
      }
    }
  }

  # Apply recursive function -----------------------------------------------------------------------
  output <- list()
  for (target in targets)  {
    result <- list(get_one(target, output))
    output <- c(output, result)
  }
  names(output) <- targets
  output <- output[original_targets]

  # Apply `simplify` option ------------------------------------------------------------------------
  if (simplify) output <- unique(unname(unlist(output)))

  return(output)
}



#===================================================================================================
#' Get items associated with taxa
#'
#' Given one or more taxa IDs and the edge list defining the taxonomy, return the items
#' (e.g. sequence information) associated with each taxon.
#'
#' @param targets (\code{character}) Taxon IDs for which items will be returned.
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the subtaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#' @param items (\code{character}) Taxon IDs for a set of items.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the item assigned to the specified
#'   target taxa, not its subtaxa. If \code{TRUE}, return all the items of every subtaxa, etc.
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique item indexes.
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of item indexes are returned
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the item indexes
#'   for all \code{target} taxa are returned in a single vector.
#'
#' @export
get_taxon_items <- function(targets, taxa, parents, items, recursive = TRUE, simplify = FALSE) {
  # Argument validataion ---------------------------------------------------------------------------
  if (length(targets) == 0) stop("Argument 'targets' has 0 length")
  if (!all(targets %in% taxa)) stop("All 'targets' taxon IDs not found in 'taxa' taxon IDs")
  validate_edge_list(taxa, parents)
  parents[!(parents %in% taxa)] <- NA

  # Get level of each input to determin processing order -------------------------------------------
  levels <- vapply(get_supertaxa(targets, taxa, parents), length, numeric(1))
  original_targets <- targets
  targets <- targets[order(levels, decreasing = TRUE)]

  # Recursive function for one target --------------------------------------------------------------
  get_one <- function(target, output) {
    if (target %in% names(output)) { #if this taxa has already been processed
      return(output[[target]])
    } else {
      subtaxa <- taxa[parents == target]
      subtaxa <- subtaxa[!is.na(subtaxa)]
      target_items <- which(items == target)
      if (recursive) {
        if (length(subtaxa) == 0) {
          return(target_items)
        } else {
          subtaxa_items <- unlist(lapply(subtaxa, get_one, output = output))
          return(c(target_items, subtaxa_items))
        }
      } else {
        return(target_items)
      }
    }
  }

  # Apply recursive function -----------------------------------------------------------------------
  output <- list()
  for (target in targets)  {
    result <- list(get_one(target, output))
    output <- c(output, result)
  }
  names(output) <- targets
  output <- output[original_targets]

  # Apply `simplify` option ------------------------------------------------------------------------
  if (simplify) output <- unique(unname(unlist(output)))

  return(output)
}



#===================================================================================================
#' Validate a taxon edge list
#'
#' Validate a taxon edge list for use in other functions
#'
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#'
#' @return (\code{logical}) Returns \code{TRUE} if the list is valid.
validate_edge_list <- function(taxa, parents) {
  if (length(taxa) == 0) stop("'taxa' has 0 length")
  if (length(parents) == 0) stop("'parents' has 0 length")
  if (length(taxa) != length(parents)) stop("'taxa' and 'parents' are of unequal length")
  if (!all(parents %in% c(taxa, NA))) stop("All 'parent' taxon IDs not found in 'taxa' taxon IDs")
  if (length(unique(taxa)) != length(taxa)) stop("All 'taxa' not unique.")
}

