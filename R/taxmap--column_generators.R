#' Get classification of taxa 
#'
#' Get classification strings of taxa in an object of type \code{\link{taxmap}}.
#' Each classification is constructed by concatenating the taxon ids of the given taxon and its supertaxa.
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character})
#' The \code{taxon_ids}s or indexes in \code{taxon_data} to get classifications for.
#' @param sep (\code{character} of length 1)
#' The character(s) to place between taxon IDs
#'
#' @return \code{character} of length equal to \code{subset}
#'
#' @export
hierarchies <- function(obj, subset = 1:nrow(obj$taxon_data), sep = ";") {
  vapply(supertaxa(obj, subset = subset, recursive = TRUE, include_input = TRUE, index = FALSE, na = FALSE),
         function(x) paste0(rev(x), collapse = sep), character(1))
}


#' Get rank of taxa 
#'
#' Get rank of taxa in an object of type \code{\link{taxmap}}
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids}s or indexes in \code{taxon_data} to get ranks for.
#'
#' @return \code{numeric}
#'
#' @export
taxon_levels <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  vapply(supertaxa(obj, subset = subset, recursive = TRUE, include_input = TRUE, index = TRUE, na = FALSE),
         length, numeric(1))
}


#' Count observations in \code{\link{taxmap}}
#'
#' Count observations for each taxon in \code{\link{taxmap}}.
#' This includes observations for the specific taxon and its subtaxa.
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids}s or indexes in \code{taxon_data} to get counts for.
#'
#' @return \code{numeric}
#'
#' @export
n_obs <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  vapply(obs(obj, subset = subset, recursive = TRUE, simplify = FALSE), length, numeric(1))
}


#' Count observation assigned in \code{\link{taxmap}}
#'
#' Count observations assigned to a specific taxon in \code{\link{taxmap}}.
#' This does not include observations assigned to subtaxa.
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids}s or indexes in \code{taxon_data} to get counts for.
#'
#' @return \code{numeric}
#'
#' @export
n_obs_1 <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  vapply(obs(obj, subset = subset, recursive = FALSE, simplify = FALSE), length, numeric(1))
}