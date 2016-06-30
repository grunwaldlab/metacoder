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
#' @family taxon_funcs
#'
#' @export
hierarchies <- function(obj, subset = 1:nrow(obj$taxon_data), sep = ";") {
  vapply(supertaxa(obj, subset = subset, recursive = TRUE, include_input = TRUE, index = FALSE, na = FALSE),
         function(x) paste0(rev(x), collapse = sep), character(1))
}


#' Get number of supertaxa 
#'
#' Get  number of supertaxa for each taxon in an object of type \code{\link{taxmap}}
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids} or indexes in \code{taxon_data}.
#'
#' @return \code{numeric}
#' 
#' @family taxon_funcs
#'
#' @export
n_supertaxa <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  vapply(supertaxa(obj, subset = subset, recursive = TRUE, include_input = FALSE, index = TRUE, na = FALSE),
         length, numeric(1))
}

#' Get number of subtaxa 
#'
#' Get number of subtaxa for each taxon in an object of type \code{\link{taxmap}}
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids} or indexes in \code{taxon_data}.
#'
#' @return \code{numeric}
#' 
#' @family taxon_funcs
#'
#' @export
n_subtaxa <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  vapply(subtaxa(obj, subset = subset, recursive = TRUE, include_input = FALSE, index = TRUE),
         length, numeric(1))
}

#' Get number of subtaxa 
#'
#' Get number of subtaxa for each taxon in an object of type \code{\link{taxmap}}, not including subtaxa of subtaxa etc. 
#' This does not include subtaxa assigned to subtaxa.
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids} or indexes in \code{taxon_data}.
#'
#' @return \code{numeric}
#'
#' @family taxon_funcs
#' 
#' @export
n_subtaxa_1 <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  vapply(subtaxa(obj, subset = subset, recursive = FALSE, include_input = FALSE, index = TRUE),
         length, numeric(1))
}


#' Count observations in \code{\link{taxmap}}
#'
#' Count observations for each taxon in \code{\link{taxmap}}.
#' This includes observations for the specific taxon and its subtaxa.
#'
#' @param obj (\code{\link{taxmap}})
#' @param subset (\code{character}) The \code{taxon_ids} or indexes in \code{taxon_data} to get counts for.
#'
#' @return \code{numeric}
#' 
#' @family taxon_funcs
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
#' @param subset (\code{character}) The \code{taxon_ids} or indexes in \code{taxon_data} to get counts for.
#'
#' @return \code{numeric}
#' 
#' @family taxon_funcs
#'
#' @export
n_obs_1 <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  vapply(obs(obj, subset = subset, recursive = FALSE, simplify = FALSE), length, numeric(1))
}