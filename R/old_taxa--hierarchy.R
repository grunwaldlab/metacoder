#' The Hierarchy class
#'
#' A class containing an ordered list of [taxon()] objects that represent a
#' hierarchical classification.
#'
#' @export
#' @param ... Any number of object of class `Taxon` or taxonomic names as
#' character strings
#' @param .list An alternate to the `...` input. Any number of object of class
#'   [taxon()] or character vectors in a list. Cannot be used with `...`.
#' @return An `R6Class` object of class `Hierarchy`
#'
#' @details On initialization, taxa are sorted if they have ranks with a known
#'   order.
#'
#' **Methods**
#'   \describe{
#'     \item{`pop(rank_names)`}{
#'       Remove `Taxon` elements by rank name, taxon name or taxon ID. The
#'       change happens in place, so you don't need to assign output to a new
#'       object. returns self - rank_names (character) a vector of rank names
#'     }
#'     \item{`pick(rank_names)`}{
#'       Select `Taxon` elements by rank name, taxon name or taxon ID. The
#'       change happens in place, so you don't need to assign output to a new
#'       object. returns self - rank_names (character) a vector of rank names
#'     }
#'   }
#'
#' @family classes
#' @examples
#' (x <- taxon(
#'   name = taxon_name("Poaceae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(4479)
#' ))
#'
#' (y <- taxon(
#'   name = taxon_name("Poa"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(4544)
#' ))
#'
#' (z <- taxon(
#'   name = taxon_name("Poa annua"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(93036)
#' ))
#'
#' (res <- hierarchy(z, y, x))
#'
#' res$taxa
#' res$ranklist
#'
#' # pop off a rank
#' pop(res, ranks("family"))
#'
#' # pick a rank
#' (res <- hierarchy(z, y, x))
#' pick(res, ranks("family"))
#'
#'
#' # null taxa
#' x <- taxon(NULL)
#' (res <- hierarchy(x, x, x))
#' ## similar to hierarchy(), but `taxa` slot is not empty

hierarchy <- function(..., .list = NULL) {
  Hierarchy$new(..., .list = .list)
}

Hierarchy <- R6::R6Class(
  "Hierarchy",
  lock_objects = TRUE,
  public = list(
    taxa = NULL,
    ranklist = NULL,

    initialize = function(..., .list = NULL) {
      # Get intput
      input <- unlist(get_dots_or_list(..., .list = .list))

      if (!all(vapply(input, function(x)
        any(class(x) %in% c('character', 'Taxon')), logical(1)))
      ) {
        stop(
          "all inputs to 'hierarchy' must be of class 'Taxon' or 'character'",
          call. = FALSE)
      }

      # Convert factors to characters
      fact_input_index <- which(lapply(input, class) == "factor")
      input[fact_input_index] <- lapply(input[fact_input_index], as.character)

      # If character strings are supplied, convert to taxa
      char_input_index <- which(lapply(input, class) == "character")
      input[char_input_index] <- lapply(input[char_input_index], taxon)

      # Parse input
      all_have_ranks <- all(vapply(input,
                                   function(x) !is.null(x$rank$name),
                                   logical(1)))
      if (all_have_ranks) {
        self$taxa <- private$sort_hierarchy(input)
      } else {
        self$taxa <- input
      }
    },

    print = function(indent = "") {
      cat(paste0(indent, "<Hierarchy>\n"))
      if (length(self$taxa) > 0 &&
          !all(vapply(self$taxa, function(z) z$is_empty(), logical(1)))) {
        cat("  no. taxon's: ", length(self$taxa), "\n")
        for (i in seq_along(self$taxa[1:min(10, length(self$taxa))])) {
          cat(
            sprintf("  %s / %s / %s",
                    self$taxa[[i]]$get_name() %||% "",
                    self$taxa[[i]]$get_rank() %||% "",
                    self$taxa[[i]]$get_id() %||% ""
            ), "\n")
        }
        if (length(self$taxa) > 10) cat("  ...")
      } else {
        cat("  Empty hierarchy")
      }
      invisible(self)
    },

    pop = function(ranks = NULL, names = NULL, ids = NULL) {
      if (all_empty(self$taxa)) stop("no taxa found")
      alldat <- ct(unlist(c(ranks, names, ids), TRUE))
      if (is.null(alldat) || length(alldat) == 0) {
        stop("one of 'ranks', 'names', or 'ids' must be used")
      }
      taxa_rks <- vapply(self$taxa, function(x) x$rank$name, "")
      taxa_nms <- vapply(self$taxa, function(x) x$name$name, "")
      taxa_ids <- vapply(self$taxa, function(x) x$id$id, numeric(1))
      todrop <- which(taxa_rks %in% ranks |
                        taxa_nms %in% names |
                        taxa_ids %in% ids)
      self$taxa[todrop] <- NULL
      private$sort_hierarchy(self$taxa)
      return(self)
    },

    pick = function(ranks = NULL, names = NULL, ids = NULL) {
      if (all_empty(self$taxa)) stop("no taxa found")
      alldat <- ct(unlist(c(ranks, names, ids), TRUE))
      if (is.null(alldat) || length(alldat) == 0) {
        stop("one of 'ranks', 'names', or 'ids' must be used")
      }
      taxa_rks <- vapply(self$taxa, function(x) x$rank$name, "")
      taxa_nms <- vapply(self$taxa, function(x) x$name$name, "")
      taxa_ids <- vapply(self$taxa, function(x) x$id$id, numeric(1))
      todrop <- which(!(taxa_rks %in% ranks |
                        taxa_nms %in% names |
                        taxa_ids %in% ids))
      self$taxa[todrop] <- NULL
      private$sort_hierarchy(self$taxa)
      return(self)
    },

    span = function(ranks = NULL, names = NULL, ids = NULL) {
      if (all_empty(self$taxa)) stop("no taxa found")
      alldat <- ct(unlist(c(ranks, names, ids), TRUE))
      if (is.null(alldat) || length(alldat) == 0) {
        stop("one of 'ranks', 'names', or 'ids' must be used")
      }

      taxa_rks <- vapply(self$taxa, function(x) x$rank$name, "")

      if (length(ranks) != 0) {
        if (!is.null(attr(ranks[[1]], "operator"))) {
          ranks <- private$make_ranks2(ranks)
        } else {
          # if no operator, names must be length > 1
          if (length(ranks[[1]]) != 2) {
            stop("if no operator, must pass in 2 names")
          }
          ranks <- private$do_ranks(ranks[[1]])
        }
      }
      if (length(names) != 0) {
        if (!is.null(attr(names[[1]], "operator"))) {
          ranks <- private$make_ranks2(private$taxaswaprank(names))
        } else {
          # if no operator, names must be length > 1
          if (length(names[[1]]) != 2) {
            stop("if no operator, must pass in 2 names")
          }
          ranks <- private$do_ranks(private$taxa2rank(names[[1]]))
        }
      }
      if (length(ids) != 0) {
        if (!is.null(attr(ids[[1]], "operator"))) {
          ranks <- private$make_ranks2(private$idsswaprank(ids))
        } else {
          # if no operator, names must be length > 1
          if (length(ids[[1]]) != 2) {
            stop("if no operator, must pass in 2 names")
          }
          ranks <- private$do_ranks(private$ids2rank(ids[[1]]))
        }
      }

      todrop <- which(!taxa_rks %in% unique(ranks))
      self$taxa[todrop] <- NULL
      private$sort_hierarchy(self$taxa)
      return(self)
    }

  ),

  private = list(
    sort_hierarchy = function(x) {
      if (length(x) == 0) {
        return(x)
      }
      ranks <- tolower(vapply(x, function(z) z$rank$name, ""))
      # check that each rank is in the acceptable set
      invisible(lapply(ranks, function(z) {
        if (!z %in% private$poss_ranks()) {
          stop(z, " not in the acceptable set of rank names", call. = FALSE)
        }
      }))
      self$ranklist <- as.list(vapply(ranks, which_ranks, numeric(1)))
      x[order(unname(unlist(self$ranklist)))]
    },

    poss_ranks = function() {
      unique(
        do.call(
          c,
          sapply(ranks_ref$ranks, strsplit, split = ",", USE.NAMES = FALSE)
        )
      )
    },

    do_ranks = function(x) {
      idz <- vapply(x, which_ranks, numeric(1), USE.NAMES = FALSE)
      keep <- ranks_ref[ranks_ref$rankid >= idz[1] &
                          ranks_ref$rankid <= idz[2], ]
      csep2vec(keep$ranks)
    },

    make_ranks2 = function(x) {
      idz <- vapply(x, which_ranks, numeric(1), USE.NAMES = FALSE)
      op <- vapply(x, attr, "", which = "operator")
      funs <- lapply(op, function(z) {
        switch(z, `>` = `<`, `>=` = `<=`, `<` = `>`, `<=` = `>=`)
      })
      logs <- if (length(idz) > 1) {
        do.call(`&`, Map(function(a, b) eval(a)(as.numeric(ranks_ref$rankid), b), funs, idz))
      } else {
        eval(funs[[1]])(as.numeric(ranks_ref$rankid), idz[[1]])
      }
      keep <- ranks_ref[logs, ]
      csep2vec(keep$ranks)
    },

    taxa2rank = function(x) {
      tmp <- vapply(self$taxa, function(z) z$name$name, "")
      rcks <- vapply(self$taxa, function(z) z$rank$name, "")
      rcks[which(tmp %in% unlist(x))]
    },

    taxaswaprank = function(x) {
      tmp <- vapply(self$taxa, function(z) z$name$name, "")
      rcks <- vapply(self$taxa, function(z) z$rank$name, "")
      for (i in seq_along(x)) {
        x[[i]] <- structure(
          rcks[which(tmp %in% x[[i]][[1]])],
          type = "ranks",
          operator = attr(x[[i]], "operator"),
          class = "taxapicker"
        )
      }
      return(x)
    },

    ids2rank = function(w) {
      tmp <- vapply(self$taxa, function(z) z$id$id, 1)
      rcks <- vapply(self$taxa, function(z) z$rank$name, "")
      rcks[which(tmp %in% w)]
    },

    idsswaprank = function(x) {
      tmp <- vapply(self$taxa, function(z) z$id$id, 1)
      rcks <- vapply(self$taxa, function(z) z$rank$name, "")
      for (i in seq_along(x)) {
        x[[i]] <- structure(
          rcks[which(tmp %in% x[[i]][[1]])],
          type = "ranks",
          operator = attr(x[[i]], "operator"),
          class = "taxapicker"
        )
      }
      return(x)
    }

  )
)

which_ranks <- function(x) {
  as.numeric(ranks_ref[which(sapply(ranks_ref$ranks, function(z) {
    any(unlist(strsplit(z, split = ",")) == x)
  }, USE.NAMES = FALSE)), "rankid"])
}

all_empty <- function(x) {
  all(vapply(x, function(z) z$is_empty(), logical(1)))
}
