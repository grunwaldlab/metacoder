#' @title Span taxa
#'
#' @description Select a range of taxa, either by two names, or relational
#' operators
#'
#' @export
#' @param .data Input, object of class `Hierarchy`, or `hierarchies`
#' @param ... quoted rank names (e.g., family) via [ranks()], taxon names
#' (e.g., Poa annua) via [nms()], or taxonomic IDs (e.g., 93036) via [ids()].
#' You can't pass in arbitrary strings or numbers.
#' @details supports `Hierarchy` and `hierarchies` objects
#' @return an object of the same class as passed in
#' @template span_egs
#' @seealso See [filtering-helpers], including for more explanation
#' of how this function works.
span <- function(.data, ...) {
  UseMethod("span")
}

#' @export
span.default <- function(.data, ...) {
  stop("no 'span' method for ", class(.data), call. = FALSE)
}

#' @export
span.Hierarchy <- function(.data, ...) {
  .data <- .data$clone(deep = TRUE)
  tmp <- Taxapickers$new(...)
  if (length(tmp$x) == 0) stop("no acceptable selectors passed in")
  .data$span(ranks = tmp$ranks(), names = tmp$names(), ids = tmp$ids())
}

#' @export
span.hierarchies <- function(.data, ...) {
  hierarchies(.list = lapply(.data, span, ...))
}
