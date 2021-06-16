#' A class for multiple taxon objects
#'
#' Stores one or more [taxon()] objects. This is just a thin wrapper for a list
#' of [taxon()] objects.
#'
#' This is the documentation for the class called `taxa`. If you are looking for
#' the documentation for the package as a whole: [taxa-package].
#'
#' @export
#' @param ... Any number of object of class [taxon()]
#' @param .list An alternate to the `...` input. Any number of object of class
#'   [taxon()]. Cannot be used with `...`.
#' @return An `R6Class` object of class `Taxon`
#' @family classes
#' @examples
#' (a <- taxon(
#'   name = taxon_name("Poa annua"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(93036)
#' ))
#' taxa(a, a, a)
#'
#' # a null set
#' x <- taxon(NULL)
#' taxa(x, x, x)
#'
#' # combo non-null and null
#' taxa(a, x, a)
taxa <- function(..., .list = NULL) {
  tt <- get_dots_or_list(..., .list = .list)
  if (!all(vapply(tt, inherits, logical(1), what = "Taxon"))) {
    stop("all inputs to 'taxa' must be of class 'Taxon'",
         call. = FALSE)
  }
  structure(tt, class = "taxa")
}

#' @export
print.taxa <- function(x, ...) {
  cat("<taxa>", "\n")
  cat("  no. taxa: ", length(x), "\n")
  if (length(x)) {
    if (all(vapply(x, function(z) z$is_empty(), logical(1)))) {
      cat("   empty set", "\n")
    } else {
      for (i in seq_along(x[1:min(10, length(x))])) {
        if (x[[i]]$is_empty()) {
          cat("  empty", "\n")
        } else {
          cat(
            sprintf("  %s / %s / %s",
                    x[[i]]$name$name %||% "",
                    x[[i]]$rank$name %||% "",
                    x[[i]]$id$id %||% ""
            ), "\n")
        }
      }
    }
  }
  if (length(x) > 10) cat("  ...")
}
