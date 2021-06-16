#' @title Pop taxa out
#'
#' @description Pop out taxa, that is, drop them
#'
#' @export
#' @param .data Input, object of class `Hierarchy`, or `hierarchies`
#' @param ... quoted rank names (e.g., family) via [ranks()], taxon names
#' (e.g., Poa annua) via [nms()], or taxonomic IDs (e.g., 93036) via [ids()].
#' You can't pass in arbitrary strings or numbers.
#' @details supports `Hierarchy` and `hierarchies` objects
#' @return an object of the same class as passed in
#' @template pop_egs
#' @seealso See [filtering-helpers], including for more explanation
#' of how this function works.
pop <- function(.data, ...) {
  UseMethod("pop")
}

#' @export
pop.default <- function(.data, ...) {
  stop("no 'pop' method for ", class(.data), call. = FALSE)
}

#' @export
pop.Hierarchy <- function(.data, ...) {
  .data <- .data$clone(deep = TRUE)
  tmp <- Taxapickers$new(...)
  if (length(tmp$x) == 0) stop("no acceptable selectors passed in")
  .data$pop(ranks = unlist(tmp$ranks()),
             names = unlist(tmp$names()),
             ids = unlist(tmp$ids()))
}

#' @export
pop.hierarchies <- function(.data, ...){
  hierarchies(.list = lapply(.data, pop, ...))
}

# pop.Taxonomy <- function(.data, ...){
#   .data <- .data$clone(deep = TRUE)
#   tmp <- unlist(lapply(lazyeval::lazy_dots(...), lazyeval::lazy_eval), FALSE)
#   .data$pop(ranks = tmp$ranks, names = tmp$names, ids = tmp$ids)
# }
