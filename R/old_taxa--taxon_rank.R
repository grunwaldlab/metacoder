#' Taxon rank class
#'
#' Stores the rank of a taxon. This is typically used to store where taxon
#' information came from in [taxon()] objects.
#'
#' @export
#' @param name (character) rank name. required
#' @param database (character) database class object, optional
#'
#' @return An `R6Class` object of class `TaxonRank`
#' @family classes
#'
#' @examples
#' taxon_rank("species")
#' taxon_rank("genus")
#' taxon_rank("kingdom")
#'
#' (x <- taxon_rank(
#'   "species",
#'   database_list$ncbi
#' ))
#' x$rank
#' x$database
#'
#' # a null taxon_name object
#' taxon_name(NULL)
taxon_rank <- function(name, database = NULL) {
  TaxonRank$new(
    name = name,
    database = database
  )
}

TaxonRank <- R6::R6Class(
  "TaxonRank",
  public = list(
    name = NULL,
    database = NULL,

    initialize = function(name = NULL, database = NULL) {
      assert(name, c("character", "TaxonName"))
      assert(database, c("character", "TaxonDatabase"))

      # Convert characters to appropriate classes
      if (is.character(database)) {
        database <- taxon_database(database)
      }

      self$name <- name
      self$database <- database
    },

    print = function(indent = "") {
      cat(paste0(indent, sprintf("<TaxonRank> %s\n", self$name %||% "none")))
      cat(paste0(indent, paste0("  database: ",
                                get_database_name(self$database) %||% "none",
                                "\n")))
      invisible(self)
    }
  )
)
