#' Make a set of many [hierarchy()] class objects
#'
#' Make a set of many [hierarchy()] class objects.
#' This is just a thin wrapper over a standard list.
#'
#' @export
#' @param ... Any number of object of class [hierarchy()]
#' @param .list Any number of object of class [hierarchy()] in a list
#'
#' @family classes
#'
#' @return An `R6Class` object of class [hierarchy()]
#'
#' @examples
#' x <- taxon(
#'   name = taxon_name("Poaceae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(4479)
#' )
#' y <- taxon(
#'   name = taxon_name("Poa"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(4544)
#' )
#' z <- taxon(
#'   name = taxon_name("Poa annua"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(93036)
#' )
#' hier1 <- hierarchy(z, y, x)
#'
#' a <- taxon(
#'   name = taxon_name("Felidae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(9681)
#' )
#' b <- taxon(
#'   name = taxon_name("Puma"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(146712)
#' )
#' c <- taxon(
#'   name = taxon_name("Puma concolor"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(9696)
#' )
#' hier2 <- hierarchy(c, b, a)
#'
#' d <- taxon(
#'   name = taxon_name("Chordata"),
#'   rank = taxon_rank("phylum"),
#'   id = taxon_id(158852)
#' )
#' e <- taxon(
#'   name = taxon_name("Vertebrata"),
#'   rank = taxon_rank("subphylum"),
#'   id = taxon_id(331030)
#' )
#' f <- taxon(
#'   name = taxon_name("Teleostei"),
#'   rank = taxon_rank("class"),
#'   id = taxon_id(161105)
#' )
#' g <- taxon(
#'   name = taxon_name("Salmonidae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(161931)
#' )
#' h <- taxon(
#'   name = taxon_name("Salmo"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(161994)
#' )
#' i <- taxon(
#'   name = taxon_name("Salmo salar"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(161996)
#' )
#' hier3 <- hierarchy(d, e, f, g, h, i)
#'
#' hiers <- hierarchies(hier1, hier2, hier3)
#'
#' # pass into the .list parameter
#' hierarchies(.list = list(hier1, hier2, hier3))
hierarchies <- function(..., .list = NULL) {
  # Get intput
  input <- get_dots_or_list(..., .list = .list)

  if (!all(vapply(input, inherits, logical(1), what = "Hierarchy"))) {
    stop("all inputs to 'hierarchies' must be of class 'Hierarchy'",
         call. = FALSE)
  }
  structure(input, class = "hierarchies")
}

#' @export
print.hierarchies <- function(x, ...) {
  cat("<Hierarchies>", "\n")
  cat("  no. hierarchies: ", length(x), "\n")
  if (length(x)) {
    for (i in seq_along(x[1:min(10, length(x))])) {
      if (is.null(x[[i]]$taxa)) {
        cat("  Empty hierarchy", sep = "\n")
      } else {
        cat(
          paste0("  ", paste0(vapply(x[[i]]$taxa, function(x) x$name$name, ""),
                              collapse = " / ")),
          "\n"
        )
      }
    }
  }
  if (length(x) > 10) cat("  ...")
}
