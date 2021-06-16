#' Taxon name class
#'
#' Used to store the name of taxa. This is typically used to
#' store where taxon names in [taxon()] objects.
#'
#' @export
#' @param name (character) a taxonomic name. required
#' @param database (character) database class object, optional
#'
#' @return An `R6Class` object of class `TaxonName`
#'
#' @family classes
#' @examples
#' (poa <- taxon_name("Poa"))
#' (undef <- taxon_name("undefined"))
#' (sp1 <- taxon_name("species 1"))
#' (poa_annua <- taxon_name("Poa annua"))
#' (x <- taxon_name("Poa annua L."))
#'
#' x$name
#' x$database
#'
#' (x <- taxon_name(
#'   "Poa annua",
#'   database_list$ncbi
#' ))
#' x$rank
#' x$database
#'
#' # a null taxon_name object
#' taxon_name(NULL)
taxon_name <- function(name, database = NULL) {
  TaxonName$new(
    name = name,
    database = database
  )
}

TaxonName <- R6::R6Class(
  "TaxonName",
  public = list(
    name = NULL,
    database = NULL,

    initialize = function(
      name = NULL, database = NULL
    ) {
      assert(name, "character")
      assert(database, c("character", "TaxonDatabase"))

      # Convert characters to appropriate classes
      if (is.character(database)) {
        database <- taxon_database(database)
      }

      self$name <- name
      self$database <- database
    },

    print = function(indent = "") {
      cat(paste0(indent, sprintf("<TaxonName> %s\n", self$name %||% "none")))
      cat(paste0(indent, paste0("  database: ",
                                get_database_name(self$database) %||% "none",
                                "\n")))
      invisible(self)
    }
  )
)
