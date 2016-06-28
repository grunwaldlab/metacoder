#' Get classification of taxa 
#'
#' Get classification strings of taxa in an object of type \code{\link{taxmap}}.
#' Each classification is constructed by concatenating the taxon ids of the given taxon and its supertaxa.
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character})
#' The \code{taxon_ids}s to get classifications for.
#' @param sep (\code{character} of length 1)
#' The character(s) to place between taxon IDs
#'
#' @return \code{character} of length equal to \code{subset}
#'
#' @export
classifications <- function(obj, subset = obj$taxon_data$taxon_ids, sep = ";") {
  vapply(supertaxa(obj, subset = subset, recursive = TRUE, include_input = TRUE, index = FALSE, na = FALSE),
         function(x) paste0(rev(x), collapse = sep), character(1))
}


#' Get rank of taxa 
#'
#' Get rank of taxa in an object of type \code{\link{taxmap}}
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids}s to get ranks for.
#'
#' @return \code{numeric}
#'
#' @export
taxon_levels <- function(obj, subset = obj$taxon_data$taxon_ids) {
  vapply(supertaxa(obj, subset = subset, recursive = TRUE, include_input = TRUE, index = TRUE, na = FALSE),
         length, numeric(1))
}


#' Count observations in \code{\link{taxmap}}
#'
#' Count observations in \code{\link{taxmap}}
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids}s to get counts for.
#'
#' @return \code{numeric}
#'
#' @export
n_obs <- function(obj, subset = obj$taxon_data$taxon_ids) {
  vapply(obs(obj, subset = subset, recursive = TRUE, simplify = FALSE), length, numeric(1))
}


#' Count observation assigned in \code{\link{taxmap}}
#'
#' Count observations assigned to a specific taxon in \code{\link{taxmap}}
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids}s to get counts for.
#'
#' @return \code{numeric}
#'
#' @export
assigned <- function(obj, subset = obj$taxon_data$taxon_ids) {
  vapply(obs(obj, subset = subset, recursive = FALSE, simplify = FALSE), length, numeric(1))
}