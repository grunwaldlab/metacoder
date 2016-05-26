#' Get classification of taxa 
#'
#' Get classification strings of taxa in an object of type \code{\link{classified}}.
#' Each classification is constructed by concatenating the taxon ids of the given taxon and its supertaxa.
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character})
#' The \code{taxon_ids}s to get classifications for.
#' @param sep (\code{character} of length 1)
#' The character(s) to place between taxon IDs
#'
#' @return \code{character} of length equal to \code{subset}
#'
#' @export
classifications <- function(obj, subset = obj$taxon_data$taxon_ids, sep = ";") {
  vapply(supertaxa(obj, subset, include_input = TRUE), function(x) paste0(rev(x), collapse = sep), character(1))
}


#' Get rank of taxa 
#'
#' Get rank of taxa in an object of type \code{\link{classified}}
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character}) The \code{taxon_ids}s to get ranks for.
#'
#' @return \code{numeric}
#'
#' @export
taxon_ranks <- function(obj, subset = obj$taxon_data$taxon_ids) {
  vapply(supertaxa(obj, subset, include_input = TRUE), length, numeric(1))
}


#' Count items in \code{\link{classified}}
#'
#' Count items in \code{\link{classified}}
#'
#' @param obj (\code{\link{classified}})
#' @param subset (\code{character}) The \code{taxon_ids}s to get counts for.
#'
#' @return \code{numeric}
#'
#' @export
item_counts <- function(obj, subset = obj$taxon_data$taxon_ids) {
  vapply(items(obj, subset), length, numeric(1))
}