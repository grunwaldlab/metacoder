#' Make a set of many [hierarchy()] class objects
#'
#' NOTE: This will soon be depreciated.
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
