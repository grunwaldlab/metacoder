#' Retrieve classifications from item IDs
#' 
#' Retrieve taxonomic classifications from item (e.g. sequence) IDs using a specified database.
#' 
#' @param item_id (\code{character})
#' An unique item (e.g. sequence) identifier for a particular \code{database}.
#' Requires an internet connection. 
#' @param database (\code{character} of length 1)
#' The name of the database that patterns given in  \code{parser} will apply to.
#' Valid databases include "ncbi", "itis", "eol", "col", "tropicos",
#' "nbn", and "none". \code{"none"} will cause no database to be quired; use this if you want to not use the
#' internet. NOTE: Only \code{"ncbi"} has been tested so far.
#' @param ... Not used
#' 
#' @return \code{list} of \code{data.frame}
#' 
#' @keywords internal
class_from_item_id <- function(item_id, database = c("ncbi", "none"), ...) {
  
  using_ncbi <- function(item_id) {
    taxize::classification(taxize::genbank2uid(item_id))
  }
  
  using_none <- function(item_id) {
    rep(NA, length(item_id))
  }
  
  database <- match.arg(database)
  map_unique(item_id, get(paste0("using_", database)))
}

# id_from_name_funcs <- list(ncbi = taxize::get_uid,
#                            itis = taxize::get_tsn,
#                            eol = taxize::get_eolid,
#                            col = taxize::get_colid,
#                            tropicos = taxize::get_tpsid,
#                            nbn = taxize::get_nbnid, 
#                            none = NA)
