# -----------------------------------------------------------------------------
#' @export
obs <- function(obj, data, value = NULL, subset = NULL, recursive = TRUE,
                simplify = FALSE) {
  UseMethod("obs")
}

#' @export
obs.default <- function(obj, data, value = NULL, subset = NULL, recursive = TRUE,
                        simplify = FALSE) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
obs.Taxmap <- function(obj, data, value = NULL, subset = NULL, recursive = TRUE,
                       simplify = FALSE) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$obs(data, value = value, subset = subset, recursive = recursive,
                          simplify = simplify)))
}


# -----------------------------------------------------------------------------
#' @export
obs_apply <- function(obj, data, func, simplify = FALSE, value = NULL,
                      subset = NULL, recursive = TRUE, ...) {
  UseMethod("obs_apply")
}

#' @export
obs_apply.default <- function(obj, data, func, simplify = FALSE, value = NULL,
                              subset = NULL, recursive = TRUE, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
obs_apply.Taxmap <- function(obj, data, func, simplify = FALSE, value = NULL,
                             subset = NULL, recursive = TRUE, ...) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$obs_apply(data, func, simplify = simplify, value = value,
                                subset = subset, recursive = recursive, ...)))
}


# -----------------------------------------------------------------------------
#' @export
filter_obs <- function(obj, data, ..., drop_taxa = FALSE, drop_obs = TRUE,
                       subtaxa = FALSE, supertaxa = TRUE,
                       reassign_obs = FALSE, target = NULL) {
  UseMethod("filter_obs")
}

#' @export
filter_obs.default <- function(obj, data, ..., drop_taxa = FALSE, drop_obs = TRUE,
                               subtaxa = FALSE, supertaxa = TRUE,
                               reassign_obs = FALSE, target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
filter_obs.Taxmap <- function(obj, data, ..., drop_taxa = FALSE, drop_obs = TRUE,
                              subtaxa = FALSE, supertaxa = TRUE,
                              reassign_obs = FALSE, target = NULL) {
  obj <- obj$clone(deep = TRUE)
  obj$filter_obs(data, ..., drop_taxa = drop_taxa, drop_obs = drop_obs,
                 subtaxa = subtaxa, supertaxa = supertaxa,
                 reassign_obs = reassign_obs, target = target)
}


# -----------------------------------------------------------------------------
#' @export
select_obs <- function(obj, data, ..., target = NULL) {
  UseMethod("select_obs")
}

#' @export
select_obs.default <- function(obj, data, ..., target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
select_obs.Taxmap <- function(obj, data, ..., target = NULL) {
  obj <- obj$clone(deep = TRUE)
  obj$select_obs(data, ..., target = target)
}


# -----------------------------------------------------------------------------
#' @export
mutate_obs <- function(obj, data, ..., target = NULL) {
  UseMethod("mutate_obs")
}

#' @export
mutate_obs.default <- function(obj, data, ..., target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
mutate_obs.Taxmap <- function(obj, data, ..., target = NULL) {
  obj <- obj$clone(deep = TRUE)
  obj$mutate_obs(data, ..., target = target)
}


# -----------------------------------------------------------------------------
#' @export
transmute_obs <- function(obj, data, ..., target = NULL) {
  UseMethod("transmute_obs")
}

#' @export
transmute_obs.default <- function(obj, data, ..., target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
transmute_obs.Taxmap <- function(obj, data, ..., target = NULL) {
  obj <- obj$clone(deep = TRUE)
  obj$transmute_obs(data, ..., target = target)
}


# -----------------------------------------------------------------------------
#' @export
arrange_obs <- function(obj, data, ..., target = NULL) {
  UseMethod("arrange_obs")
}

#' @export
arrange_obs.default <- function(obj, data, ..., target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
arrange_obs.Taxmap <- function(obj, data, ..., target = NULL) {
  obj <- obj$clone(deep = TRUE)
  obj$arrange_obs(data, ..., target = target)
}


# -----------------------------------------------------------------------------
#' @export
sample_n_obs <- function(obj, data, size, replace = FALSE, taxon_weight = NULL,
                         obs_weight = NULL, use_supertaxa = TRUE,
                         collapse_func = mean, ..., target = NULL) {
  UseMethod("sample_n_obs")
}

#' @export
sample_n_obs.default <- function(obj, data, size, replace = FALSE, taxon_weight = NULL,
                                 obs_weight = NULL, use_supertaxa = TRUE,
                                 collapse_func = mean, ..., target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
sample_n_obs.Taxmap <- function(obj, data, size, replace = FALSE, taxon_weight = NULL,
                                obs_weight = NULL, use_supertaxa = TRUE,
                                collapse_func = mean, ..., target = NULL) {
  obj <- obj$clone(deep = TRUE)
  eval(substitute(obj$sample_n_obs(data, size, replace = replace, taxon_weight = taxon_weight,
                                   obs_weight = obs_weight, use_supertaxa = use_supertaxa,
                                   collapse_func = collapse_func, ..., target = target)))
}


# -----------------------------------------------------------------------------
#' @export
sample_frac_obs <- function(obj, data, size, replace = FALSE,
                            taxon_weight = NULL, obs_weight = NULL,
                            use_supertaxa = TRUE,
                            collapse_func = mean, ..., target = NULL) {
  UseMethod("sample_frac_obs")
}

#' @export
sample_frac_obs.default <- function(obj, data, size, replace = FALSE,
                                    taxon_weight = NULL, obs_weight = NULL,
                                    use_supertaxa = TRUE,
                                    collapse_func = mean, ..., target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
sample_frac_obs.Taxmap <- function(obj, data, size, replace = FALSE,
                                   taxon_weight = NULL, obs_weight = NULL,
                                   use_supertaxa = TRUE,
                                   collapse_func = mean, ..., target = NULL) {
  obj <- obj$clone(deep = TRUE)
  eval(substitute(obj$sample_frac_obs(data, size, replace = replace,
                                      taxon_weight = taxon_weight, obs_weight = obs_weight,
                                      use_supertaxa = use_supertaxa,
                                      collapse_func = collapse_func, ..., target = target)))
}


# -----------------------------------------------------------------------------
#' @export
n_obs <- function(obj, data = NULL, target = NULL) {
  UseMethod("n_obs")
}

#' @export
n_obs.default <- function(obj, data = NULL, target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_obs.Taxmap <- function(obj, data = NULL, target = NULL) {
  obj$n_obs(data = data, target = target)
}


# -----------------------------------------------------------------------------
#'@export
n_obs_1 <- function(obj, data = NULL, target = NULL) {
  UseMethod("n_obs_1")
}

#' @export
n_obs_1.default <- function(obj, data = NULL, target = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_obs_1.Taxmap <- function(obj, data = NULL, target = NULL) {
  obj$n_obs_1(data = data, target = target)
}


# -----------------------------------------------------------------------------
#'@export
get_dataset <- function(obj, data = NULL) {
  UseMethod("get_dataset")
}

#' @export
get_dataset.default <- function(obj, data = NULL) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
get_dataset.Taxmap <- function(obj, data = NULL) {
  obj$get_dataset(data = data)
}
