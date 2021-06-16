ct <- function(l) Filter(Negate(is.null), l)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

add_names <- function(...) {
  tt <- sapply(match.call(), deparse)[-1]
  nms <- unname(sapply(tt, function(x) strsplit(x, split = "\\$")[[1]][2]))
  stats::setNames(list(...), nms)
}

# assert <- function(x, y) {
#   if (!is.null(x)) {
#     if (!class(x) %in% y) {
#       stop(deparse(substitute(x)), " must be of class ",
#            paste0(y, collapse = ", "), call. = FALSE)
#     }
#   }
# }

assert <- function(x, y) {
  if (!is.null(x)) {
    if (!any(class(x) %in% y)) {
      stop(deparse(substitute(x)), " must be of class ",
           paste0(y, collapse = ", "), call. = FALSE)
    }
  }
}

csep2vec <- function(z, uniq = TRUE) {
  tmp <- unlist(lapply(z, function(w) strsplit(w, split = ",")[[1]]), FALSE)
  if (uniq) return(unique(tmp))
  return(tmp)
}

strex <- function(str, pattern) regmatches(str, regexpr(pattern, str))
