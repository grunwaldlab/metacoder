#' Filter taxa with a list of conditions
#' 
#' Filter taxa in a \code{\link{taxmap}} object with a list of conditions. Any column name that 
#' appears in \code{taxon_data(.data)} can be used as if it was a vector on its own. See 
#' \code{\link[dplyr]{filter}} for inspiration and more information.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more filtering conditions. This can be one of three things: \describe{ 
#'   \item{\code{character}}{One or more \code{taxon_id}s} \item{\code{integer}}{One or more indexes
#'   of \code{taxon_data}} \item{\code{logical}}{A \code{TRUE}/\code{FALSE} vector of length equal 
#'   to the number of rows in \code{taxon_data}} } Any column name that appears in 
#'   \code{taxon_data(.data)} can be used as if it was a vector on its own.
#' @param subtaxa (\code{logical} of length 1) If \code{TRUE}, include subtaxa of taxa passing the 
#'   filter.
#' @param supertaxa (\code{logical} of length 1) If \code{TRUE}, include supertaxa of taxa passing 
#'   the filter.
#' @param taxonless (\code{logical} of length 1) If \code{TRUE}, include observations even if the
#'   taxon they are assigned to is filtered out. observation assigned to removed taxa will be
#'   assigned to \code{NA}. See the \code{reassign} option below for further complications.
#' @param reassign_obs (\code{logical} of length 1) If \code{TRUE}, observations assigned to removed
#'   taxa will be reassigned to the closest supertaxon that passed the filter. If there are no
#'   supertaxa of such an observation that passed the filter, they will be filtered out if
#'   \code{taxonless} is \code{TRUE}.
#' @param reassign_taxa (\code{logical} of length 1) If \code{TRUE}, subtaxa of removed taxa will be
#'   reassigned to the closest supertaxon that passed the filter.
#' @param invert (\code{logical} of length 1) If \code{TRUE}, do NOT include the selection.
#'   This is different than just replacing a \code{==} with a \code{!=} because this option negates
#'   the selection after taking into account the \code{subtaxa} and \code{supertaxa} options.
#'   This is useful for removing a taxon and all its subtaxa for example.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#' 
#' @examples
#' # Remove singleton taxa, but reassign singletons to supertaxa that pass filter
#' filter_taxa(unite_ex_data_3, n_obs > 1)
#' # Remove singleton taxa and associated seqeuence data
#' filter_taxa(unite_ex_data_3, n_obs > 1, taxonless = FALSE, reassign_obs = FALSE)
#' # Subset to a single taxon and its subtaxa
#' filter_taxa(unite_ex_data_3, name == "Basidiomycota", subtaxa = TRUE)
#' # Remove a taxon and its subtaxa
#' filter_taxa(unite_ex_data_3, name == "Basidiomycota", subtaxa = TRUE, invert = TRUE)
#' # Remove taxa, reassigning supertaxa and subtaxa 
#' filter_taxa(unite_ex_data_3, unite_rank != "p")
#' 
#' @export
filter_taxa <- function(.data, ..., subtaxa = FALSE, supertaxa = FALSE, taxonless = FALSE,
                        reassign_obs = TRUE, reassign_taxa = TRUE, invert = FALSE) {
  
  # non-standard argument evaluation ---------------------------------------------------------------
  selection <- lazyeval::lazy_eval(lazyeval::lazy_dots(...),
                                   data = taxon_data(.data, col_subset = taxon_data_cols_used(.data, ...))) 
  
  # convert taxon_ids to logical -------------------------------------------------------------------
  is_char <- vapply(selection, is.character, logical(1))
  selection[is_char] <- lapply(selection[is_char], function(x) .data$taxon_data$taxon_ids %in% x)
  
  # convert indexes to logical ---------------------------------------------------------------------
  is_index <- vapply(selection, is.numeric, logical(1))
  selection[is_index] <- lapply(selection[is_index], function(x) 1:nrow(.data$taxon_data) %in% x)
  
  # combine filters --------------------------------------------------------------------------------
  selection <- Reduce(`&`, selection)
  
  # Get taxa of subset -----------------------------------------------------------------------------
  taxa_subset <- unique(c(which(selection),
                          if (subtaxa) {
                            subtaxa(.data, subset = selection, recursive = TRUE, index = TRUE,
                                    include_input = FALSE, simplify = TRUE)
                          },
                          if (supertaxa) {
                            supertaxa(.data, subset = selection, recursive = TRUE, index = TRUE,
                                      na = FALSE, simplify = TRUE, include_input = FALSE)
                          }))
  
  # Invert selection -------------------------------------------------------------------------------
  if (invert) {
    taxa_subset <- (1:nrow(.data$taxon_data))[-taxa_subset]
  }
  
  # Reassign taxonless observations ----------------------------------------------------------------
  if (reassign_obs) {
    reassign_one <- function(parents) {
      included_parents <- parents[parents %in% taxa_subset]
      return(.data$taxon_data$taxon_ids[included_parents[1]])
    }
    
    to_reassign <- ! .data$obs_data$obs_taxon_ids %in% .data$taxon_data$taxon_ids[taxa_subset]
    supertaxa_key <- supertaxa(.data, subset = unique(.data$obs_data$obs_taxon_ids[to_reassign]),
                              recursive = TRUE, simplify = FALSE, include_input = FALSE, index = TRUE, na = FALSE)
    reassign_key <- vapply(supertaxa_key, reassign_one, character(1))
    .data$obs_data[to_reassign, "obs_taxon_ids"] <- reassign_key[.data$obs_data$obs_taxon_ids[to_reassign]]
  }
  
  # Reassign subtaxa  ------------------------------------------------------------------------------
  if (reassign_taxa) {
    reassign_one <- function(parents) {
      included_parents <- parents[parents %in% taxa_subset]
      return(.data$taxon_data$taxon_ids[included_parents[1]])
    }
    
    to_reassign <- ! .data$taxon_data$supertaxon_ids %in% .data$taxon_data$taxon_ids[taxa_subset]
    supertaxa_key <- supertaxa(.data, subset = unique(.data$taxon_data$taxon_ids[to_reassign]),
                               recursive = TRUE, simplify = FALSE, include_input = FALSE, index = TRUE, na = FALSE)
    reassign_key <- vapply(supertaxa_key, reassign_one, character(1))
    .data$taxon_data[to_reassign, "supertaxon_ids"] <- reassign_key[.data$taxon_data$taxon_ids[to_reassign]]
  }
  
  
  # Remove taxonless observations -------------------------------------------------------------------------
  obs_subset <- .data$obs_data$obs_taxon_ids %in% .data$taxon_data$taxon_ids[taxa_subset]
  if (taxonless) {
    .data$obs_data[! obs_subset, "obs_taxon_ids"] <- as.character(NA)
  } else {
    .data$obs_data <- .data$obs_data[obs_subset, , drop = FALSE]
  }
  
  # Remove filtered taxa ---------------------------------------------------------------------------
  .data$taxa <- .data$taxa[.data$taxon_data$taxon_ids[taxa_subset]]
  .data$taxon_data <- .data$taxon_data[taxa_subset, , drop = FALSE]
  .data$taxon_data[! .data$taxon_data$supertaxon_ids %in% .data$taxon_data$taxon_ids, "supertaxon_ids"] <- as.character(NA)
  
  return(.data)
}



#' Filter observations with a list of conditions
#' 
#' Filter observations in a \code{\link{taxmap}} object with a list of conditions. Any column name that
#' appears in \code{obs_data(.data)} can be used as if it was a vector on its own. See 
#' \code{\link[dplyr]{filter}} for inspiration and more information.
#' 
#' @param .data \code{\link{taxmap}}
#' @param ... One or more filtering conditions. This can be one of two things: \describe{ 
#'   \item{\code{integer}}{One or more indexes of \code{obs_data}} \item{\code{logical}}{A 
#'   \code{TRUE}/\code{FALSE} vector of length equal to the number of rows in \code{obs_data}} } 
#'   Any column name that appears in \code{obs_data(.data)} can be used as if it was a vector on 
#'   its own.
#' @param unobserved (\code{logical} of length 1) If \code{TRUE}, preserve taxa even if all of their 
#'   observations are filtered out. If \code{FALSE}, remove taxa for which all observations were filtered out. 
#'   Note that only taxa that are unobserved due to this filtering will be removed; there might be 
#'   other taxa without observations to begin with that will not be removed.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'
#' @examples 
#' # Filter by sequence name, but preserve all taxa
#' filter_obs(unite_ex_data_3, grepl("Lachnum", seq_name))
#' # Filter by sequence name and only keep taxa with sequences that pass the filter
#' filter_obs(unite_ex_data_3, grepl("Lachnum", seq_name), unobserved = FALSE)
#' 
#' @export
filter_obs <- function(.data, ..., unobserved = TRUE) {
  # non-standard argument evaluation ---------------------------------------------------------------
  selection <- lazyeval::lazy_eval(lazyeval::lazy_dots(...), data = obs_data(.data)) 
  
  # convert taxon_ids to indexes -------------------------------------------------------------------
  is_char <- vapply(selection, is.character, logical(1))
  if (sum(is_char) > 0) {
    stop("observation filtering with taxon IDs or observation IDs (which dont exist yet) is not currently supported. If you want to filter observation by taxon IDs, use something like: `obs_taxon_ids %in% my_subset`")
  }
  # selection[is_char] <- lapply(selection[is_char], function(x) match(x, .data$taxon_data$taxon_ids))
  
  # convert logical to indexes ---------------------------------------------------------------------
  is_logical <- vapply(selection, is.logical, logical(1))
  selection[is_logical] <- lapply(selection[is_logical], which)
  
  # combine filters --------------------------------------------------------------------------------
  intersect_with_dups <-function(a, b) {
    #taken from http://r.789695.n4.nabble.com/intersect-without-discarding-duplicates-td2225377.html
    rep(sort(intersect(a, b)), pmin(table(a[a %in% b]), table(b[b %in% a])))
  }
  selection <- Reduce(intersect_with_dups, selection)
  
  # Remove observations -----------------------------------------------------------------------------------
  unobserved_taxa <- supertaxa(.data, unique(.data$obs_data$obs_taxon_ids[-selection]), na = FALSE,
                             recursive = TRUE, simplify = TRUE, include_input = TRUE, index = TRUE)
  .data$obs_data <- .data$obs_data[selection, , drop = FALSE]
  
  # Remove unobserved taxa ---------------------------------------------------------------------------
  if (! unobserved) {
    taxa_to_remove <- 1:nrow(.data$taxon_data) %in% unobserved_taxa & n_obs(.data) == 0
    .data$taxon_data <- .data$taxon_data[! taxa_to_remove, , drop = FALSE]
    .data$taxon_data[! .data$taxon_data$supertaxon_ids %in% .data$taxon_data$taxon_ids, "supertaxon_ids"] <- as.character(NA)
  }
  
  return(.data)
}

