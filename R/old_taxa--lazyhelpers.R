#' Taxonomic filtering helpers
#'
#' @name filtering-helpers
#' @param ... quoted rank names, taxonomic names, taxonomic ids, or
#' any of those with supported operators (See \strong{Supported Relational
#' Operators} below)
#'
#' @note NSE is not supported at the moment, but may be in the future
#'
#' @section How do these functions work?:
#' Each function assigns some metadata so we can more easily process
#' your query downstream. In addition, we check for whether you've
#' used any relational operators and pull those out to make downstream
#' processing easier
#'
#' The goal of these functions is to make it easy to combine queries
#' based on each of rank names, taxonomic names, and taxonomic ids.
#'
#' These are designed to be used inside of [pop()], [pick()], [span()]. Inside
#' of those functions, we figure out what rank names you want to filter
#' on, then check against a reference dataset ([ranks_ref]) to allow
#' ordered queries like \emph{I want all taxa between Class and Genus}. If you
#' provide rank names, we just use those, then do the filtering you requested.
#' If you provide taxonomic names or ids we figure out what rank names you are
#' referring to, then we can proceed as in the previous sentence.
#'
#' @section Supported Relational Operators:
#' \itemize{
#'  \item `>` all items above rank of x
#'  \item `>=` all items above rank of x, inclusive
#'  \item `<` all items below rank of x
#'  \item `<=` all items below rank of x, inclusive
#' }
#'
#' @section ranks:
#' Ranks can be any character string in the set of acceptable rank
#' names.
#'
#' @section nms:
#' `nms` is named to avoid using `names` which would collide with the
#' fxn [base::names()] in Base R. Can pass in any character taxonomic names.
#'
#' @section ids:
#' Ids are any alphanumeric taxonomic identifier. Some database providers
#' use all digits, but some use a combination of digits and characters.
#'
#' @examples
#' ranks("genus")
#' ranks("order", "genus")
#' ranks("> genus")
#'
#' nms("Poaceae")
#' nms("Poaceae", "Poa")
#' nms("< Poaceae")
#'
#' ids(4544)
#' ids(4544, 4479)
#' ids("< 4479")
NULL

#' @export
#' @rdname filtering-helpers
ranks <- function(...) {
  helpers_fxn_se("ranks", ...)
}

#' @export
#' @rdname filtering-helpers
nms <- function(...) {
  helpers_fxn_se("names", ...)
}

#' @export
#' @rdname filtering-helpers
ids <- function(...) {
  helpers_fxn_se("ids", ...)
}

###### ---------------------
helpers_fxn_se <- function(name, ...) {
  x <- unlist(list(...))
  if (is.null(x)) {
    stop("you must pass in one or more values, see ?filtering-helpers")
  }
  # check for operators
  op <- NULL
  if (any(grepl(">|>=|<|<=", x))) {
    op <- strex(x, ">|>=|<|<=")
    x <- gsub(">|>=|<|<=|\\s+", "", x)
  }
  structure(x, class = 'taxapicker', type = name, operator = op)
}

print.taxapicker <- function(x, ...) {
  cat("<taxapicker>", sep = "\n")
  cat(sprintf(
    "  (%s) (operator: `%s`): %s", attr(x, "type"), attr(x, "operator") %||% "",
    paste0(unclass(x), collapse = ", ")
  ), sep = "\n")
}

Taxapickers <- R6::R6Class(
  "Taxapickers",
  lock_objects = TRUE,
  public = list(
    x = NULL,

    initialize = function(...) {
      self$x <- Filter(function(x) inherits(x, "taxapicker"), list(...))
    },

    print = function() {
      cat(paste0("<Taxapickers> n=", length(self$x)), sep = "\n")
      for (i in seq_along(self$x)) {
        cat(sprintf(
          "  (%s): %s", attr(self$x[[i]], "type"),
          paste0(unclass(self$x[[i]]), collapse = ", ")
        ), sep = "\n")
      }
    },

    ranks = function() private$pluck("ranks"),
    names = function() private$pluck("names"),
    ids = function() private$pluck("ids")
  ),

  private = list(
    pluck = function(z) {
      self$x[vapply(self$x, attr, "", which = "type") == z]
    }
  )
)


## keeping around if/when we want to support NSE
# helpers_fxn <- function(name, ...) {
#   rks <- rlang::quos(...)
#   clzzs <- vapply(rks, function(z) class(rlang::quo_expr(z)), "",
#                   USE.NAMES = FALSE)
#   if (any(clzzs == "call")) {
#     rlang::quo_expr(rks[[1]])
#     tmp <- as.character(rlang::quo_expr(rks[[1]]))
#     ops <- paste0(rev(strex(tmp, "^>$|^<$|^\\.$|^:$|^::$")), collapse = " ")
#     structure(list(strex(tmp, "[A-Za-z0-9]+")),
#               operator = ops, names = name,
#               class = 'taxapicker', type = name)
#   } else if (any(clzzs == "name")) {
#     rlang::eval_tidy(rks[[1]])
#   } else {
#     tmp <- vapply(rks, function(z) rlang::quo_name(z), "", USE.NAMES = FALSE)
#     structure(list(tmp), names = name, class = 'taxapicker', type = name)
#   }
# }
