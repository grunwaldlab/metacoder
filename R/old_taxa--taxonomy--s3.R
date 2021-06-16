# -----------------------------------------------------------------------------
#' @export
taxon_ids <- function(obj) {
  UseMethod("taxon_ids")
}

#' @export
taxon_ids.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
taxon_ids.Taxonomy <- function(obj) {
  obj$taxon_ids()
}


# -----------------------------------------------------------------------------
#' @export
taxon_indexes <- function(obj, ...) {
  UseMethod("taxon_indexes")
}

#' @export
taxon_indexes.default <- function(obj, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
taxon_indexes.Taxonomy <- function(obj, ...) {
  obj$taxon_indexes(...)
}


# -----------------------------------------------------------------------------
#' @export
taxon_names <- function(obj, ...) {
  UseMethod("taxon_names")
}

#' @export
taxon_names.default <- function(obj, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
taxon_names.Taxonomy <- function(obj, ...) {
  obj$taxon_names(...)
}


# -----------------------------------------------------------------------------
#' @export
taxon_ranks <- function(obj) {
  UseMethod("taxon_ranks")
}

#' @export
taxon_ranks.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
taxon_ranks.Taxonomy <- function(obj) {
  obj$taxon_ranks()
}


# -----------------------------------------------------------------------------
#' @export
supertaxa <- function(obj, subset = NULL, recursive = TRUE, simplify = FALSE,
                      include_input = FALSE, value = "taxon_indexes", na = FALSE) {
  UseMethod("supertaxa")
}

#' @export
supertaxa.default <- function(obj, subset = NULL, recursive = TRUE, simplify = FALSE,
                              include_input = FALSE, value = "taxon_indexes", na = FALSE) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
supertaxa.Taxonomy <- function(obj, subset = NULL, recursive = TRUE, simplify = FALSE,
                               include_input = FALSE, value = "taxon_indexes", na = FALSE) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$supertaxa(subset = subset, recursive = recursive, simplify = simplify,
                                include_input = include_input, value = value, na = na)))
}


# -----------------------------------------------------------------------------
#' @export
supertaxa_apply <- function(obj, func, subset = NULL, recursive = TRUE,
                            simplify = FALSE, include_input = FALSE,
                            value = "taxon_indexes", na = FALSE, ...) {
  UseMethod("supertaxa_apply")
}

#' @export
supertaxa_apply.default <- function(obj, func, subset = NULL, recursive = TRUE,
                                    simplify = FALSE, include_input = FALSE,
                                    value = "taxon_indexes", na = FALSE, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
supertaxa_apply.Taxonomy <- function(obj, func, subset = NULL, recursive = TRUE,
                                     simplify = FALSE, include_input = FALSE,
                                     value = "taxon_indexes", na = FALSE, ...) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$supertaxa_apply(func, subset = subset, recursive = recursive,
                                      simplify = simplify, include_input = include_input,
                                      value = value, na = na, ...)))
}


# -----------------------------------------------------------------------------
#' @export
roots <- function(obj, subset = NULL, value = "taxon_indexes") {
  UseMethod("roots")
}

#' @export
roots.default <- function(obj, subset = NULL, value = "taxon_indexes") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
roots.Taxonomy <- function(obj, subset = NULL, value = "taxon_indexes") {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$roots(subset = subset, value = value)))
}


# -----------------------------------------------------------------------------
#' @export
branches <- function(obj, subset = NULL, value = "taxon_indexes") {
  UseMethod("branches")
}

#' @export
branches.default <- function(obj, subset = NULL, value = "taxon_indexes") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
branches.Taxonomy <- function(obj, subset = NULL, value = "taxon_indexes") {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$branches(subset = subset, value = value)))
}


# -----------------------------------------------------------------------------
#' @export
internodes <- function(obj, subset = NULL, value = "taxon_indexes") {
  UseMethod("internodes")
}

#' @export
internodes.default <- function(obj, subset = NULL, value = "taxon_indexes") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
internodes.Taxonomy <- function(obj, subset = NULL, value = "taxon_indexes") {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$internodes(subset = subset, value = value)))
}


# -----------------------------------------------------------------------------
#' @export
subtaxa <- function(obj, subset = NULL, recursive = TRUE,
                    simplify = FALSE, include_input = FALSE,
                    value = "taxon_indexes") {
  UseMethod("subtaxa")
}

#' @export
subtaxa.default <- function(obj, subset = NULL, recursive = TRUE,
                            simplify = FALSE, include_input = FALSE,
                            value = "taxon_indexes") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
subtaxa.Taxonomy <- function(obj, subset = NULL, recursive = TRUE,
                             simplify = FALSE, include_input = FALSE,
                             value = "taxon_indexes") {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$subtaxa(subset = subset, recursive = recursive,
                              simplify = simplify, include_input = include_input,
                              value = value)))
}


# -----------------------------------------------------------------------------
#' @export
subtaxa_apply <- function(obj, func, subset = NULL, recursive = TRUE,
                          simplify = FALSE, include_input = FALSE,
                          value = "taxon_indexes", ...) {
  UseMethod("subtaxa_apply")
}

#' @export
subtaxa_apply.default <- function(obj, func, subset = NULL, recursive = TRUE,
                                  simplify = FALSE, include_input = FALSE,
                                  value = "taxon_indexes", ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
subtaxa_apply.Taxonomy <- function(obj, func, subset = NULL, recursive = TRUE,
                                   simplify = FALSE, include_input = FALSE,
                                   value = "taxon_indexes", ...) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$subtaxa_apply(func, subset = subset, recursive = recursive,
                                    simplify = simplify, include_input = include_input,
                                    value = value, ...)))
}


# -----------------------------------------------------------------------------
#' @export
stems <- function(obj, subset = NULL, value = "taxon_indexes", simplify = FALSE,
                  exclude_leaves = FALSE) {
  UseMethod("stems")
}

#' @export
stems.default <- function(obj, subset = NULL, value = "taxon_indexes", simplify = FALSE,
                          exclude_leaves = FALSE) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
stems.Taxonomy <- function(obj, subset = NULL, value = "taxon_indexes", simplify = FALSE,
                           exclude_leaves = FALSE) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$stems(subset = subset, value = value, simplify = simplify,
                            exclude_leaves = exclude_leaves)))
}


# -----------------------------------------------------------------------------
#' @export
leaves <- function(obj, subset = NULL, recursive = TRUE, simplify = FALSE, value = "taxon_indexes") {
  UseMethod("leaves")
}

#' @export
leaves.default <- function(obj, subset = NULL, recursive = TRUE, simplify = FALSE, value = "taxon_indexes") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
leaves.Taxonomy <- function(obj, subset = NULL, recursive = TRUE, simplify = FALSE, value = "taxon_indexes") {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$leaves(subset = subset, recursive = recursive, simplify = simplify, value = value)))
}


# -----------------------------------------------------------------------------
#' @export
leaves_apply <- function(obj, func, subset = NULL, recursive = TRUE, simplify = FALSE,
                         value = "taxon_indexes", ...) {
  UseMethod("leaves_apply")
}

#' @export
leaves_apply.default <- function(obj, func, subset = NULL, recursive = TRUE, simplify = FALSE,
                                 value = "taxon_indexes", ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
leaves_apply.Taxonomy <- function(obj, func, subset = NULL, recursive = TRUE, simplify = FALSE,
                                  value = "taxon_indexes", ...) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$leaves_apply(func, subset = subset, recursive = recursive,
                                   simplify = simplify, value = value, ...)))
}


# -----------------------------------------------------------------------------
#' @export
id_classifications <- function(obj, sep = ";") {
  UseMethod("id_classifications")
}

#' @export
id_classifications.default <- function(obj, sep = ";") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
id_classifications.Taxonomy <- function(obj, sep = ";") {
  obj$id_classifications(sep = sep)
}


# -----------------------------------------------------------------------------
#' @export
classifications <- function(obj, value = "taxon_names", sep = ";") {
  UseMethod("classifications")
}

#' @export
classifications.default <- function(obj, value = "taxon_names", sep = ";") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
classifications.Taxonomy <- function(obj, value = "taxon_names", sep = ";") {
  obj$classifications(value = value, sep = sep)
}


# -----------------------------------------------------------------------------
#' @export
n_supertaxa <- function(obj) {
  UseMethod("n_supertaxa")
}

#' @export
n_supertaxa.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_supertaxa.Taxonomy <- function(obj) {
  obj$n_supertaxa()
}


# -----------------------------------------------------------------------------
#' @export
n_supertaxa_1 <- function(obj) {
  UseMethod("n_supertaxa_1")
}

#' @export
n_supertaxa_1.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_supertaxa_1.Taxonomy <- function(obj) {
  obj$n_supertaxa_1()
}


# -----------------------------------------------------------------------------
#' @export
n_subtaxa <- function(obj) {
  UseMethod("n_subtaxa")
}

#' @export
n_subtaxa.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_subtaxa.Taxonomy <- function(obj) {
  obj$n_subtaxa()
}


# -----------------------------------------------------------------------------
#' @export
n_subtaxa_1 <- function(obj) {
  UseMethod("n_subtaxa_1")
}

#' @export
n_subtaxa_1.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_subtaxa_1.Taxonomy <- function(obj) {
  obj$n_subtaxa_1()
}


# -----------------------------------------------------------------------------
#' @export
n_leaves <- function(obj) {
  UseMethod("n_leaves")
}

#' @export
n_leaves.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_leaves.Taxonomy <- function(obj) {
  obj$n_leaves()
}


# -----------------------------------------------------------------------------
#' @export
n_leaves_1 <- function(obj) {
  UseMethod("n_leaves_1")
}

#' @export
n_leaves_1.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
n_leaves_1.Taxonomy <- function(obj) {
  obj$n_leaves_1()
}


# -----------------------------------------------------------------------------
#' @export
all_names <- function(obj, tables = TRUE, funcs = TRUE, others = TRUE,
                      builtin_funcs = TRUE, warn = FALSE) {
  UseMethod("all_names")
}

#' @export
all_names.default <- function(obj, tables = TRUE, funcs = TRUE, others = TRUE,
                              builtin_funcs = TRUE, warn = FALSE) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
all_names.Taxonomy <- function(obj, tables = TRUE, funcs = TRUE, others = TRUE,
                               builtin_funcs = TRUE, warn = FALSE) {
  obj$all_names()
}

#' @export
all_names.Taxmap <- function(obj, tables = TRUE, funcs = TRUE, others = TRUE,
                             builtin_funcs = TRUE, warn = FALSE) {
  obj$all_names(tables = tables, funcs = funcs, others = others,
                builtin_funcs = builtin_funcs, warn = warn)
}


# -----------------------------------------------------------------------------
#' @export
get_data <- function(obj, name = NULL, ...) {
  UseMethod("get_data")
}

#' @export
get_data.default <- function(obj, name = NULL, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
get_data.Taxonomy <- function(obj, name = NULL, ...) {
  obj$get_data(name = name, ...)
}


# -----------------------------------------------------------------------------
#' @export
get_data_frame <- function(obj, ...) {
  UseMethod("get_data_frame")
}

#' @export
get_data_frame.default <- function(obj, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
get_data_frame.Taxonomy <- function(obj, ...) {
  obj$get_data_frame(...)
}


# -----------------------------------------------------------------------------
#' @export
filter_taxa <- function(obj, ..., subtaxa = FALSE, supertaxa = FALSE,
                        drop_obs = TRUE, reassign_obs = TRUE,
                        reassign_taxa = TRUE, invert = FALSE) {
  UseMethod("filter_taxa")
}

#' @export
filter_taxa.default <- function(obj, ..., subtaxa = FALSE, supertaxa = FALSE,
                                drop_obs = TRUE, reassign_obs = TRUE,
                                reassign_taxa = TRUE, invert = FALSE) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
filter_taxa.Taxonomy <- function(obj, ..., subtaxa = FALSE, supertaxa = FALSE,
                                 drop_obs = TRUE, reassign_obs = TRUE,
                                 reassign_taxa = TRUE, invert = FALSE,
                                 keep_order = TRUE) {
  # Check that a taxmap option is not used with a taxonomy object
  if (! "Taxmap" %in% class(obj)) {
    if (!missing(reassign_obs)) {
      warning(paste('The option "reassign_obs" can only be used with',
                    '`taxmap` objects. It will have no effect on a',
                    '`taxonomy` object.'))
    }
    if (!missing(drop_obs)) {
      warning(paste('The option "drop_obs" can only be used with',
                    '`taxmap` objects. It will have no effect on a',
                    '`taxonomy` object.'))
    }
  }

  # Execute R6 function
  obj <- obj$clone(deep = TRUE)
  obj$filter_taxa(..., subtaxa = subtaxa, supertaxa = supertaxa,
                  reassign_taxa = reassign_taxa, invert = invert,
                  keep_order = keep_order)
}

#' @export
filter_taxa.Taxmap <- function(obj, ..., subtaxa = FALSE, supertaxa = FALSE,
                               drop_obs = TRUE, reassign_obs = TRUE,
                               reassign_taxa = TRUE, invert = FALSE,
                               keep_order = TRUE) {
  obj <- obj$clone(deep = TRUE)
  obj$filter_taxa(..., subtaxa = subtaxa, supertaxa = supertaxa,
                  drop_obs = drop_obs, reassign_obs = reassign_obs,
                  reassign_taxa = reassign_taxa, invert = invert,
                  keep_order = keep_order)
}


# -----------------------------------------------------------------------------
#' @export
arrange_taxa <- function(obj, ...) {
  UseMethod("arrange_taxa")
}

#' @export
arrange_taxa.default <- function(obj, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
arrange_taxa.Taxonomy <- function(obj, ...) {
  obj <- obj$clone(deep = TRUE)
  obj$arrange_taxa(...)
}


# -----------------------------------------------------------------------------
#' @export
sample_n_taxa <- function(obj, size, taxon_weight = NULL, obs_weight = NULL,
                          obs_target = NULL, use_subtaxa = TRUE,
                          collapse_func = mean, ...) {
  UseMethod("sample_n_taxa")
}

#' @export
sample_n_taxa.default <- function(obj, size, taxon_weight = NULL, obs_weight = NULL,
                                  obs_target = NULL, use_subtaxa = TRUE,
                                  collapse_func = mean, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
sample_n_taxa.Taxonomy <- function(obj, size, taxon_weight = NULL, obs_weight = NULL,
                                   obs_target = NULL, use_subtaxa = TRUE,
                                   collapse_func = mean, ...) {
  # Check that a taxmap option is not used with a taxonomy object
  if (! "Taxmap" %in% class(obj)) {
    if (!missing(obs_weight)) {
      warning(paste('The option "obs_weight" can only be used with',
                    '`taxmap` objects. It will have no effect on a',
                    '`taxonomy` object.'))
    }
    if (!missing(obs_target)) {
      warning(paste('The option "obs_target" can only be used with',
                    '`taxmap` objects. It will have no effect on a',
                    '`taxonomy` object.'))
    }
  }

  # Execute R6 function
  obj <- obj$clone(deep = TRUE)
  eval(substitute(obj$sample_n_taxa(size, taxon_weight = taxon_weight,
                                    use_subtaxa = use_subtaxa,
                                    collapse_func = collapse_func, ...)))
}

#' @export
sample_n_taxa.Taxmap <- function(obj, size, taxon_weight = NULL, obs_weight = NULL,
                                 obs_target = NULL, use_subtaxa = TRUE,
                                 collapse_func = mean, ...) {
  obj <- obj$clone(deep = TRUE)
  eval(substitute(obj$sample_n_taxa(size, taxon_weight = taxon_weight, obs_weight = obs_weight,
                                    obs_target = obs_target, use_subtaxa = use_subtaxa,
                                    collapse_func = collapse_func, ...)))
}



# -----------------------------------------------------------------------------
#' @export
sample_frac_taxa <- function(obj, size = 1, taxon_weight = NULL,
                             obs_weight = NULL, obs_target = NULL,
                             use_subtaxa = TRUE, collapse_func = mean, ...) {
  UseMethod("sample_frac_taxa")
}

#' @export
sample_frac_taxa.default <- function(obj, size = 1, taxon_weight = NULL,
                                     obs_weight = NULL, obs_target = NULL,
                                     use_subtaxa = TRUE, collapse_func = mean, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
sample_frac_taxa.Taxonomy <- function(obj, size = 1, taxon_weight = NULL,
                                      obs_weight = NULL, obs_target = NULL,
                                      use_subtaxa = TRUE, collapse_func = mean, ...) {
  obj <- obj$clone(deep = TRUE)
  eval(substitute(obj$sample_frac_taxa(size = size, taxon_weight = taxon_weight,
                                       use_subtaxa = use_subtaxa, collapse_func = collapse_func, ...)))
}

#' @export
sample_frac_taxa.Taxmap <- function(obj, size = 1, taxon_weight = NULL,
                                    obs_weight = NULL, obs_target = NULL,
                                    use_subtaxa = TRUE, collapse_func = mean, ...) {
  obj <- obj$clone(deep = TRUE)
  eval(substitute(obj$sample_frac_taxa(size = size, taxon_weight = taxon_weight,
                                       obs_weight = obs_weight, obs_target = obs_target,
                                       use_subtaxa = use_subtaxa, collapse_func = collapse_func, ...)))
}



# -----------------------------------------------------------------------------
#' @export
is_root <- function(obj) {
  UseMethod("is_root")
}

#' @export
is_root.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
is_root.Taxonomy <- function(obj) {
  obj$is_root()
}


# -----------------------------------------------------------------------------
#' @export
is_internode <- function(obj) {
  UseMethod("is_internode")
}

#' @export
is_internode.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
is_internode.Taxonomy <- function(obj) {
  obj$is_internode()
}


# -----------------------------------------------------------------------------
#' @export
is_stem <- function(obj) {
  UseMethod("is_stem")
}

#' @export
is_stem.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
is_stem.Taxonomy <- function(obj) {
  obj$is_stem()
}


# -----------------------------------------------------------------------------
#' @export
is_branch <- function(obj) {
  UseMethod("is_branch")
}

#' @export
is_branch.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
is_branch.Taxonomy <- function(obj) {
  obj$is_branch()
}


# -----------------------------------------------------------------------------
#' @export
is_leaf <- function(obj) {
  UseMethod("is_leaf")
}

#' @export
is_leaf.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
is_leaf.Taxonomy <- function(obj) {
  obj$is_leaf()
}


# -----------------------------------------------------------------------------
#' @export
map_data <- function(obj, from, to, warn = TRUE) {
  UseMethod("map_data")
}

#' @export
map_data.default <- function(obj, from, to, warn = TRUE) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
map_data.Taxonomy <- function(obj, from, to, warn = TRUE) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$map_data(from = from, to = to, warn = warn)))
}


# -----------------------------------------------------------------------------
#' @export
map_data_ <- function(obj, from, to) {
  UseMethod("map_data_")
}

#' @export
map_data_.default <- function(obj, from, to) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
map_data_.Taxonomy <- function(obj, from, to) {
  obj$map_data_(from = from, to = to)
}


# -----------------------------------------------------------------------------
#' @export
replace_taxon_ids <- function(obj, ...) {
  UseMethod("replace_taxon_ids")
}

#' @export
replace_taxon_ids.default <- function(obj, ...) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
replace_taxon_ids.Taxonomy <- function(obj, ...) {
  obj <- obj$clone(deep = TRUE)
  obj$replace_taxon_ids(...)
}


# -----------------------------------------------------------------------------
#' @export
remove_redundant_names <- function(obj) {
  UseMethod("remove_redundant_names")
}

#' @export
remove_redundant_names.default <- function(obj) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
remove_redundant_names.Taxonomy <- function(obj) {
  obj <- obj$clone(deep = TRUE)
  obj$remove_redundant_names()
}

# -----------------------------------------------------------------------------
#' @export
taxonomy_table <- function(obj, subset = NULL, value = "taxon_names",
                           use_ranks = NULL, add_id_col = FALSE) {
  UseMethod("taxonomy_table")
}

#' @export
taxonomy_table.default <- function(obj, subset = NULL, value = "taxon_names",
                                   use_ranks = NULL, add_id_col = FALSE) {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
taxonomy_table.Taxonomy <- function(obj, subset = NULL, value = "taxon_names",
                                    use_ranks = NULL, add_id_col = FALSE) {
  obj <- eval(obj) # Needed by testthat for some reason
  eval(substitute(obj$taxonomy_table(subset = subset, value = value,
                                     use_ranks = use_ranks, add_id_col = add_id_col)))
}


# -----------------------------------------------------------------------------
#' @export
print_tree <- function(obj, value = "taxon_names") {
  UseMethod("print_tree")
}

#' @export
print_tree.default <- function(obj, value = "taxon_names") {
  stop("Unsupported class: ", class(obj)[[1L]], call. = FALSE, domain = NA)
}

#' @export
print_tree.Taxonomy <- function(obj, value = "taxon_names") {
  obj$print_tree(value = value)
}
