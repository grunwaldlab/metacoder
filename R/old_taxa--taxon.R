#' Taxon class
#'
#' A class used to define a single taxon. Most other classes in the taxa package
#' include one or more objects of this class.
#'
#' @export
#' @param name a TaxonName object [taxon_name()] or character string. if character
#' passed in, we'll coerce to a TaxonName object internally, required
#' @param rank a TaxonRank object [taxon_rank()] or character string. if character
#' passed in, we'll coerce to a TaxonRank object internally, required
#' @param id a TaxonId object [taxon_id()], numeric/integer, or character string.
#' if numeric/integer/character passed in, we'll coerce to a TaxonId object
#' internally, required
#' @param authority (character) a character string, optional
#'
#' @details Note that there is a special use case of this function - you can
#' pass `NULL` as the first parameter to get an empty `taxon` object. It makes
#' sense to retain the original behavior where nothing passed in to the first
#' parameter leads to an error, and thus creating a `NULL` taxon is done very
#' explicitly.
#'
#' @return An `R6Class` object of class `Taxon`
#' @family classes
#'
#' @examples
#' (x <- taxon(
#'   name = taxon_name("Poa annua"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(93036)
#' ))
#' x$name
#' x$rank
#' x$id
#'
#' # a null taxon object
#' taxon(NULL)
#' ## with all NULL objects from the other classes
#' taxon(
#'   name = taxon_name(NULL),
#'   rank = taxon_rank(NULL),
#'   id = taxon_id(NULL)
#' )
taxon <- function(name, rank = NULL, id = NULL, authority = NULL) {
  Taxon$new(
    name = name,
    rank = rank,
    id = id,
    authority = authority
  )
}

Taxon <- R6::R6Class(
  "Taxon",
  public = list(
    name = NULL,
    rank = NULL,
    id = NULL,
    authority = NULL,

    initialize = function(
      name = NULL, rank = NULL, id = NULL, authority = NULL
    ) {
      assert(name, c('TaxonName', 'character'))
      assert(rank, c('TaxonRank', 'character'))
      assert(id, c('TaxonId', 'character', 'numeric', 'integer'))
      assert(authority, 'character')

      # Convert characters to appropriate classes
      if (is.character(name)) {
        name <- taxon_name(name)
      }
      if (is.character(rank)) {
        rank <- taxon_rank(rank)
      }
      if (is.character(id)) {
        id <- taxon_id(id)
      }

      self$name <- name
      self$rank <- rank
      self$id <- id
      self$authority <- authority
    },

    print = function(indent = "") {
      cat(paste0(indent, "<Taxon>\n"))
      cat(paste0(indent, paste0("  name: ",
                                self$get_name() %||% "none", "\n")))
      cat(paste0(indent, paste0("  rank: ",
                                self$get_rank() %||% "none", "\n")))
      cat(paste0(indent, paste0("  id: ",
                                self$get_id() %||% "none", "\n")))
      cat(paste0(indent, paste0("  authority: ",
                                self$authority %||% "none", "\n")))
      invisible(self)
    },

    is_empty = function(x) {
      is.null(self$name) && is.null(self$rank) && is.null(self$id)
    },

    get_name = function() {
      if ("TaxonName" %in% class(self$name)) {
        output <- self$name$name
      } else {
        output <- self$name
      }
      return(output)
    },

    get_rank = function() {
      if ("TaxonRank" %in% class(self$rank)) {
        output <- self$rank$name
      } else {
        output <- self$rank
      }
      return(output)
    },

    get_id = function() {
      if ("TaxonId" %in% class(self$id)) {
        output <- self$id$id
      } else {
        output <- self$id
      }
      return(output)
    }

  ),

  private = list(
  )
)
