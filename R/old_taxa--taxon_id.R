#' Taxon ID class
#'
#' Used to store taxon IDs, either arbitrary or from a taxonomy database. This
#' is typically used to store taxon IDs in [taxon()] objects.
#'
#' @export
#' @param id (character/integer/numeric) a taxonomic id, required
#' @param database (database) database class object, optional
#'
#' @return An `R6Class` object of class `TaxonId`
#' @family classes
#'
#' @examples
#' (x <- taxon_id(12345))
#' x$id
#' x$database
#'
#' (x <- taxon_id(
#'   12345,
#'   database_list$ncbi
#' ))
#' x$id
#' x$database
#'
#' # a null taxon_name object
#' taxon_name(NULL)
taxon_id <- function(id, database = NULL) {
  TaxonId$new(
    id = id,
    database = database
  )
}

TaxonId <- R6::R6Class(
  "TaxonId",
  public = list(
    id = NULL,
    database = NULL,

    initialize = function(id = NULL, database = NULL) {
      assert(id, c("character", "integer", "numeric"))
      assert(database, c("character", "TaxonDatabase"))

      # Convert characters to appropriate classes
      if (is.character(database)) {
        database <- taxon_database(database)
      }

      self$id <- id
      self$database <- database
    },

    print = function(indent = "") {
      cat(paste0(indent, sprintf("<TaxonId> %s\n", self$id %||% "none")))
      cat(paste0(indent, paste0("  database: ",
                                get_database_name(self$database) %||% "none",
                                "\n")))
      invisible(self)
    }
  )
)
