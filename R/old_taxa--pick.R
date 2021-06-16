#' @title Pick taxa
#'
#' @description Pick out specific taxa, while others are dropped
#'
#' @export
#' @param .data Input, object of class `Hierarchy`, or `hierarchies`
#' @param ... quoted rank names (e.g., family) via [ranks()], taxon names
#' (e.g., Poa annua) via [nms()], or taxonomic IDs (e.g., 93036) via [ids()].
#' You can't pass in arbitrary strings or numbers.
#' @details supports `Hierarchy` and `hierarchies` objects
#' @return an object of the same class as passed in
#' @template pick_egs
#' @seealso See [filtering-helpers], including for more explanation
#' of how this function works.
pick <- function(.data, ...) {
  UseMethod("pick")
}

#' @export
pick.default <- function(.data, ...) {
  stop("no 'pick' method for ", class(.data), call. = FALSE)
}

#' @export
pick.Hierarchy <- function(.data, ...) {
  .data <- .data$clone(deep = TRUE)
  tmp <- Taxapickers$new(...)
  if (length(tmp$x) == 0) stop("no acceptable selectors passed in")
  .data$pick(ranks = unlist(tmp$ranks()),
             names = unlist(tmp$names()),
             ids = unlist(tmp$ids()))
}

#' @export
pick.hierarchies <- function(.data, ...) {
  hierarchies(.list = lapply(.data, pick, ...))
}
