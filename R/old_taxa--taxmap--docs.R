#' Get data indexes associated with taxa
#'
#' Given a [taxmap()] object, return data associated with each taxon in a
#' given table included in that [taxmap()] object.
#' \preformatted{
#' obj$obs(data, value = NULL, subset = NULL,
#'   recursive = TRUE, simplify = FALSE)
#' obs(obj, data, value = NULL, subset = NULL,
#'   recursive = TRUE, simplify = FALSE)}
#'
#' @param obj ([taxmap()]) The [taxmap()] object containing taxon information to
#'   be queried.
#' @param data Either the name of something in `obj$data` that has taxon
#'   information or a an external object with taxon information. For tables,
#'   there must be a column named "taxon_id" and lists/vectors must be named by
#'   taxon ID.
#' @param value What data to return. This is usually the name of column in a
#'   table in `obj$data`. Any result of `all_names(obj)` can be used. If the
#'   value used has names, it is assumed that the names are taxon ids and the
#'   taxon ids are used to look up the correct values.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to find observations
#'   for. Default: All taxa in `obj` will be used. Any variable name that
#'   appears in [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the
#'   observation assigned to the specified input taxa, not subtaxa. If `TRUE`,
#'   return all the observations of every subtaxa, etc. Positive numbers
#'   indicate the number of ranks below the each taxon to get observations for
#'   `0` is equivalent to `FALSE`. Negative numbers are equivalent to `TRUE`.
#' @param simplify (`logical`) If `TRUE`, then combine all the results into a
#'   single vector of unique observation indexes.
#'
#' @return If `simplify = FALSE`, then a list of vectors of observation indexes
#'   are returned corresponding to the `data` argument. If `simplify = TRUE`,
#'   then the observation indexes for all `data` taxa are returned in a single
#'   vector.
#'
#' @name obs
#'
#' @examples
#' # Get indexes of rows corresponding to each taxon
#' obs(ex_taxmap, "info")
#'
#' # Get only a subset of taxon indexes
#' obs(ex_taxmap, "info", subset = 1:2)
#'
#' # Get only a subset of taxon IDs
#' obs(ex_taxmap, "info", subset = c("b", "c"))
#'
#' # Get only a subset of taxa using logical tests
#' obs(ex_taxmap, "info", subset = taxon_ranks == "genus")
#'
#' # Only return indexes of rows assinged to each taxon explicitly
#' obs(ex_taxmap, "info", recursive = FALSE)
#'
#' # Lump all row indexes in a single vector
#' obs(ex_taxmap, "info", simplify = TRUE)
#'
#' # Return values from a dataset instead of indexes
#' obs(ex_taxmap, "info", value = "name")
#'
NULL

#' Apply function to observations per taxon
#'
#' Apply a function to data for the observations for each taxon. This is similar
#' to using [obs()] with [lapply()] or [sapply()].
#' \preformatted{
#' obj$obs_apply(data, func, simplify = FALSE, value = NULL,
#'   subset = NULL, recursive = TRUE, ...)
#' obs_apply(obj, data, func, simplify = FALSE, value = NULL,
#'   subset = NULL, recursive = TRUE, ...)}
#'
#' @param obj The [taxmap()] object containing taxon information to
#'   be queried.
#' @param data Either the name of something in `obj$data` that has taxon
#'   information or a an external object with taxon information. For tables,
#'   there must be a column named "taxon_id" and lists/vectors must be named by
#'   taxon ID.
#' @param func (`function`) The function to apply.
#' @param simplify (`logical`) If `TRUE`, convert lists to vectors.
#' @param value What data to give to the function. This is usually the name of
#'   column in a table in `obj$data`. Any result of `all_names(obj)` can be
#'   used, but it usually only makes sense to use columns in the dataset
#'   specified by the `data` option. By default, the indexes of observation in
#'   `data` are returned.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to use.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the
#'   observation assigned to the specified input taxa, not subtaxa. If `TRUE`,
#'   return all the observations of every subtaxa, etc. Positive numbers
#'   indicate the number of ranks below the each taxon to get observations for
#'   `0` is equivalent to `FALSE`. Negative numbers are equivalent to `TRUE`.
#' @param ... Extra arguments are passed to the function.
#'
#' @name obs_apply
#'
#' @examples
#' # Find the average number of legs in each taxon
#' obs_apply(ex_taxmap, "info", mean, value = "n_legs", simplify = TRUE)
#'
#' # One way to implement `n_obs` and find the number of observations per taxon
#' obs_apply(ex_taxmap, "info", length, simplify = TRUE)
#'
NULL


#' Filter observations with a list of conditions
#'
#' Filter data in a [taxmap()] object (in `obj$data`) with a
#' set of conditions.  See
#' [dplyr::filter()] for the inspiration for this function and more
#' information. Calling the function using the `obj$filter_obs(...)` style
#' edits "obj" in place, unlike most R functions. However, calling the function
#' using the `filter_obs(obj, ...)` imitates R's traditional copy-on-modify
#' semantics, so "obj" would not be changed; instead a changed version would be
#' returned, like most R functions.
#' \preformatted{
#' obj$filter_obs(data, ..., drop_taxa = FALSE, drop_obs = TRUE,
#'                subtaxa = FALSE, supertaxa = TRUE, reassign_obs = FALSE)
#' filter_obs(obj, data, ..., drop_taxa = FALSE, drop_obs = TRUE,
#'            subtaxa = FALSE, supertaxa = TRUE, reassign_obs = FALSE)}
#'
#' @param obj An object of type [taxmap()]
#' @param data Dataset names, indexes, or a logical vector that indicates which datasets in
#'   `obj$data` to filter. If multiple datasets are filterd at once, then they must be the same
#'   length.
#' @param ... One or more filtering conditions. Any variable name that appears
#'   in [all_names()] can be used as if it was a vector on its own. Each
#'   filtering condition can be one of two things:
#'   * `integer`: One or more dataset indexes.
#'   * `logical`: A `TRUE`/`FALSE` vector of length equal to the number of
#'   items in the dataset.
#' @param drop_taxa (`logical` of length 1) If `FALSE`, preserve taxa
#'   even if all of their observations are filtered out. If `TRUE`, remove
#'   taxa for which all observations were filtered out. Note that only taxa that
#'   are unobserved due to this filtering will be removed; there might be other
#'   taxa without observations to begin with that will not be removed.
#' @param drop_obs (`logical`) This only has an effect when `drop_taxa` is
#'   `TRUE`. When `TRUE`, observations for other data sets (i.e. not `data`)
#'   assigned to taxa that are removed when filtering `data` are also removed.
#'   Otherwise, only data for taxa that are not present in all other data sets
#'   will be removed. This option can be either simply `TRUE`/`FALSE`, meaning
#'   that all data sets will be treated the same, or a logical vector can be
#'   supplied with names corresponding one or more data sets in `obj$data`. For
#'   example, `c(abundance = TRUE, stats = FALSE)` would remove observations in
#'   `obj$data$abundance`, but not in `obj$data$stats`.
#' @param subtaxa (`logical` or `numeric` of length 1) This only has an effect
#'   when `drop_taxa` is `TRUE`. If `TRUE`, include subtaxa of taxa passing the
#'   filter. Positive numbers indicate the number of ranks below the target taxa
#'   to return. `0` is equivalent to `FALSE`. Negative numbers are equivalent to
#'   `TRUE`.
#' @param supertaxa (`logical`  or `numeric` of length 1) This only has an
#'   effect when `drop_taxa` is `TRUE`. If `TRUE`, include supertaxa of taxa
#'   passing the filter. Positive numbers indicate the number of ranks above the
#'   target taxa to return. `0` is equivalent to `FALSE`. Negative numbers are
#'   equivalent to `TRUE`.
#' @param reassign_obs (`logical`) This only has an effect when `drop_taxa` is
#'   `TRUE`. If `TRUE`, observations assigned to removed taxa will be reassigned
#'   to the closest supertaxon that passed the filter. If there are no supertaxa
#'   of such an observation that passed the filter, they will be filtered out if
#'   `drop_obs` is `TRUE`. This option can be either simply `TRUE`/`FALSE`,
#'   meaning that all data sets will be treated the same, or a logical vector
#'   can be supplied with names corresponding one or more data sets in
#'   `obj$data`. For example, `c(abundance = TRUE, stats = FALSE)` would
#'   reassign observations in `obj$data$abundance`, but not in `obj$data$stats`.
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return An object of type [taxmap()]
#'
#' @examples
#' # Filter by row index
#' filter_obs(ex_taxmap, "info", 1:2)
#'
#' # Filter by TRUE/FALSE
#' filter_obs(ex_taxmap, "info", dangerous == FALSE)
#' filter_obs(ex_taxmap, "info", dangerous == FALSE, n_legs > 0)
#' filter_obs(ex_taxmap, "info", n_legs == 2)
#'
#' # Remove taxa whose obserservations were filtered out
#' filter_obs(ex_taxmap, "info", n_legs == 2, drop_taxa = TRUE)
#'
#' # Preserve other data sets while removing taxa
#' filter_obs(ex_taxmap, "info", n_legs == 2, drop_taxa = TRUE,
#'            drop_obs = c(abund = FALSE))
#'
#' # When filtering taxa, do not return supertaxa of taxa that are preserved
#' filter_obs(ex_taxmap, "info", n_legs == 2, drop_taxa = TRUE,
#'            supertaxa = FALSE)
#'
#' # Filter multiple datasets at once
#' filter_obs(ex_taxmap, c("info", "phylopic_ids", "foods"), n_legs == 2)
#'
#' @family taxmap manipulation functions
#'
#' @name filter_obs
NULL


#' Subset columns in a [taxmap()] object
#'
#' Subsets columns in a [taxmap()] object. Takes and returns a
#' [taxmap()] object. Any variable name that appears in
#' [all_names()] can be used as if it was a vector on its own. See
#' [dplyr::select()] for the inspiration for this function and more
#' information. Calling the function using the `obj$select_obs(...)` style
#' edits "obj" in place, unlike most R functions. However, calling the function
#' using the `select_obs(obj, ...)` imitates R's traditional copy-on-modify
#' semantics, so "obj" would not be changed; instead a changed version would be
#' returned, like most R functions.
#' \preformatted{
#' obj$select_obs(data, ...)
#' select_obs(obj, data, ...)}
#'
#' @param obj An object of type [taxmap()]
#' @param data Dataset names, indexes, or a logical vector that indicates which tables in
#'   `obj$data` to subset columns in. Multiple tables can be subset at once.
#' @param ... One or more column names to return in the new object. Each can be
#'   one of two things: \describe{ \item{expression with unquoted column
#'   name}{The name of a column in the dataset typed as if it was
#'   a variable on its own.} \item{`numeric`}{Indexes of columns in
#'   the dataset} } To match column names with a character vector,
#'   use `matches("my_col_name")`. To match a logical vector, convert it to
#'   a column index using `which`.
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return An object of type [taxmap()]
#'
#' @family taxmap manipulation functions
#'
#' @examples
#' # Selecting a column by name
#' select_obs(ex_taxmap, "info", dangerous)
#'
#' # Selecting a column by index
#' select_obs(ex_taxmap, "info", 3)
#'
#' # Selecting a column by regular expressions
#' select_obs(ex_taxmap, "info", matches("^n"))
#'
#' @name select_obs
NULL


#' Add columns to [taxmap()] objects
#'
#' Add columns to tables in `obj$data` in [taxmap()] objects.  See
#' [dplyr::mutate()] for the inspiration for this function and more information.
#' Calling the function using the `obj$mutate_obs(...)` style edits "obj" in
#' place, unlike most R functions. However, calling the function using the
#' `mutate_obs(obj, ...)` imitates R's traditional copy-on-modify semantics, so
#' "obj" would not be changed; instead a changed version would be returned, like
#' most R functions.
#' \preformatted{
#' obj$mutate_obs(data, ...)
#' mutate_obs(obj, data, ...)}
#'
#' @param obj An object of type [taxmap()]
#' @param data Dataset name, index, or a logical vector that indicates which dataset in
#'   `obj$data` to add columns to.
#' @param ... One or more named columns to add. Newly created columns can be
#'   referenced in the same function call. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return An object of type [taxmap()]
#'
#' @examples
#'
#' # Add column to existing tables
#' mutate_obs(ex_taxmap, "info",
#'            new_col = "Im new",
#'            newer_col = paste0(new_col, "er!"))
#'
#' # Create columns in a new table
#' mutate_obs(ex_taxmap, "new_table",
#'            nums = 1:10,
#'            squared = nums ^ 2)
#'
#' # Add a new vector
#' mutate_obs(ex_taxmap, "new_vector", 1:10)
#'
#' # Add a new list
#' mutate_obs(ex_taxmap, "new_list", list(1, 2))
#'
#' @family taxmap manipulation functions
#' @name mutate_obs
NULL


#' Replace columns in [taxmap()] objects
#'
#' Replace columns of tables in `obj$data` in [taxmap()] objects. See
#' [dplyr::transmute()] for the inspiration for this function and more
#' information. Calling the function using the `obj$transmute_obs(...)` style
#' edits "obj" in place, unlike most R functions. However, calling the function
#' using the `transmute_obs(obj, ...)` imitates R's traditional copy-on-modify
#' semantics, so "obj" would not be changed; instead a changed version would be
#' returned, like most R functions.
#' \preformatted{
#' obj$transmute_obs(data, ...)
#' transmute_obs(obj, data, ...)}
#'
#' @param obj An object of type [taxmap()]
#' @param data Dataset name, index, or a logical vector that indicates which dataset in
#'   `obj$data` to use.
#' @param ... One or more named columns to add. Newly created columns can be
#'   referenced in the same function call. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return An object of type [taxmap()]
#' @examples
#' # Replace columns in a table with new columns
#' transmute_obs(ex_taxmap, "info", new_col = paste0(name, "!!!"))
#'
#' @family taxmap manipulation functions
#'
#' @name transmute_obs
NULL


#' Sort user data in [taxmap()] objects
#'
#' Sort rows of tables  or the elements of lists/vectors in the `obj$data` list
#' in [taxmap()] objects. Any variable name that appears in [all_names()] can be
#' used as if it was a vector on its own. See [dplyr::arrange()] for the
#' inspiration for this function and more information. Calling the function
#' using the `obj$arrange_obs(...)` style edits "obj" in place, unlike most R
#' functions. However, calling the function using the `arrange_obs(obj, ...)`
#' imitates R's traditional copy-on-modify semantics, so "obj" would not be
#' changed; instead a changed version would be returned, like most R functions.
#' \preformatted{
#' obj$arrange_obs(data, ...)
#' arrange_obs(obj, data, ...)}
#'
#' @param obj An object of type [taxmap()].
#' @param data Dataset names, indexes, or a logical vector that indicates which datasets in
#'   `obj$data` to sort If multiple datasets are sorted at once, then they must be the same
#'   length.
#' @param ... One or more expressions (e.g. column names) to sort on.
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return An object of type [taxmap()]
#'
#' @examples
#' # Sort in ascending order
#' arrange_obs(ex_taxmap, "info", n_legs)
#' arrange_obs(ex_taxmap, "foods", name)
#'
#' # Sort in decending order
#' arrange_obs(ex_taxmap, "info", desc(n_legs))
#'
#' # Sort multiple datasets at once
#' arrange_obs(ex_taxmap, c("info", "phylopic_ids", "foods"), n_legs)
#'
#' @family taxmap manipulation functions
#'
#' @name arrange_obs
NULL


#' Sample n observations from [taxmap()]
#'
#' Randomly sample some number of observations from a [taxmap()] object. Weights
#' can be specified for observations or the taxa they are classified by. Any
#' variable name that appears in [all_names()] can be used as if it was a vector
#' on its own. See [dplyr::sample_n()] for the inspiration for this function.
#' Calling the function using the `obj$sample_n_obs(...)` style edits "obj" in
#' place, unlike most R functions. However, calling the function using the
#' `sample_n_obs(obj, ...)` imitates R's traditional copy-on-modify semantics,
#' so "obj" would not be changed; instead a changed version would be returned,
#' like most R functions.
#' \preformatted{
#' obj$sample_n_obs(data, size, replace = FALSE,
#'   taxon_weight = NULL, obs_weight = NULL,
#'   use_supertaxa = TRUE, collapse_func = mean, ...)
#' sample_n_obs(obj, data, size, replace = FALSE,
#'   taxon_weight = NULL, obs_weight = NULL,
#'   use_supertaxa = TRUE, collapse_func = mean, ...)}
#'
#' @param obj ([taxmap()]) The object to sample from.
#' @param data Dataset names, indexes, or a logical vector that indicates which datasets in
#'   `obj$data` to sample. If multiple datasets are sampled at once, then they must be the same
#'   length.
#' @param size (`numeric` of length 1) The number of observations to
#'   sample.
#' @param replace (`logical` of length 1) If `TRUE`, sample with
#'   replacement.
#' @param taxon_weight (`numeric`) Non-negative sampling weights of each
#'   taxon. If `use_supertaxa` is `TRUE`, the weights for each taxon
#'   in an observation's classification are supplied to `collapse_func` to
#'   get the observation weight. If `obs_weight` is also specified, the two
#'   weights are multiplied (after `taxon_weight` for each observation is
#' calculated).
#' @param obs_weight (`numeric`) Sampling weights of each observation.  If
#'   `taxon_weight` is also specified, the two weights are multiplied (after
#'   `taxon_weight` for each observation is calculated).
#' @param use_supertaxa (`logical` or `numeric` of length 1) Affects how the
#'   `taxon_weight` is used. If `TRUE`, the weights for each taxon in an
#'   observation's classification are multiplied to get the observation weight.
#'   Otherwise, just the taxonomic level the observation is assign to it
#'   considered. If `TRUE`, use all supertaxa. Positive numbers indicate the
#'   number of ranks above each taxon to use. `0` is equivalent to `FALSE`.
#'   Negative numbers are equivalent to `TRUE`.
#' @param collapse_func (`function` of length 1) If `taxon_weight` option is
#'   used and `supertaxa` is `TRUE`, the weights for each
#'   taxon in an observation's classification are supplied to
#'   `collapse_func` to get the observation weight. This function should
#'   take  numeric vector and return a single number.
#' @param ... Additional options are passed to [filter_obs()].
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return An object of type [taxmap()]
#'
#' @examples
#' # Sample 2 rows without replacement
#' sample_n_obs(ex_taxmap, "info", 2)
#' sample_n_obs(ex_taxmap, "foods", 2)
#'
#' # Sample with replacement
#' sample_n_obs(ex_taxmap, "info", 10, replace = TRUE)
#'
#' # Sample some rows for often then others
#' sample_n_obs(ex_taxmap, "info", 3, obs_weight = n_legs)
#'
#' # Sample multiple datasets at once
#' sample_n_obs(ex_taxmap, c("info", "phylopic_ids", "foods"), 3)
#'
#' @family taxmap manipulation functions
#'
#' @name sample_n_obs
NULL


#' Sample a proportion of observations from [taxmap()]
#'
#' Randomly sample some proportion of observations from a [taxmap()]
#' object. Weights can be specified for observations or their taxa. See
#' [dplyr::sample_frac()] for the inspiration for this function. Calling the
#' function using the `obj$sample_frac_obs(...)` style edits "obj" in place, unlike
#' most R functions. However, calling the function using the `sample_frac_obs(obj,
#' ...)` imitates R's traditional copy-on-modify semantics, so "obj" would not
#' be changed; instead a changed version would be returned, like most R
#' functions.
#' \preformatted{
#' obj$sample_frac_obs(data, size, replace = FALSE,
#'   taxon_weight = NULL, obs_weight = NULL,
#'   use_supertaxa = TRUE, collapse_func = mean, ...)
#' sample_frac_obs(obj, data, size, replace = FALSE,
#'   taxon_weight = NULL, obs_weight = NULL,
#'   use_supertaxa = TRUE, collapse_func = mean, ...)}
#'
#' @param obj ([taxmap()]) The object to sample from.
#' @param data Dataset names, indexes, or a logical vector that indicates which datasets in
#'   `obj$data` to sample. If multiple datasets are sample at once, then they must be the same
#'   length.
#' @param size (`numeric` of length 1) The proportion of observations to
#'   sample.
#' @param replace (`logical` of length 1) If `TRUE`, sample with
#'   replacement.
#' @param taxon_weight (`numeric`) Non-negative sampling weights of each
#'   taxon. If `use_supertaxa` is `TRUE`, the weights for each taxon
#'   in an observation's classification are supplied to `collapse_func` to
#'   get the observation weight. If `obs_weight` is also specified, the two
#'   weights are multiplied (after `taxon_weight` for each observation is
#'   calculated).
#' @param obs_weight (`numeric`) Sampling weights of each observation.  If
#'   `taxon_weight` is also specified, the two weights are multiplied
#'   (after `taxon_weight` for each observation is calculated).
#' @param use_supertaxa (`logical` or `numeric` of length 1) Affects how the
#'   `taxon_weight` is used. If `TRUE`, the weights for each taxon in
#'   an observation's classification are multiplied to get the observation
#'   weight. If `FALSE` just the taxonomic level the observation is assign to it
#'   considered. Positive numbers indicate the number of ranks above the
#'   each taxon to use. `0` is equivalent to `FALSE`. Negative numbers
#'   are equivalent to `TRUE`.
#' @param collapse_func (`function` of length 1) If `taxon_weight`
#'   option is used and `supertaxa` is `TRUE`, the weights for each
#'   taxon in an observation's classification are supplied to
#'   `collapse_func` to get the observation weight. This function should
#'   take  numeric vector and return a single number.
#' @param ... Additional options are passed to [filter_obs()].
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return An object of type [taxmap()]
#'
#' @examples
#' # Sample half of the rows fram a table
#' sample_frac_obs(ex_taxmap, "info", 0.5)
#'
#' # Sample multiple datasets at once
#' sample_frac_obs(ex_taxmap, c("info", "phylopic_ids", "foods"), 0.5)
#'
#' @family taxmap manipulation functions
#'
#' @name sample_frac_obs
NULL


#' Count observations in [taxmap()]
#'
#' Count observations for each taxon in a data set in a [taxmap()] object. This
#' includes observations for the specific taxon and the observations of its
#' subtaxa. "Observations" in this sense are the items (for list/vectors) or
#' rows (for tables) in a dataset. By default, observations in the first data
#' set in the [taxmap()] object is used.  For example, if the data set is a
#' table, then a value of 3 for a taxon means that their are 3 rows in that
#' table assigned to that taxon or one of its subtaxa.
#' \preformatted{
#' obj$n_obs(data)
#' n_obs(obj, data)}
#'
#' @param obj ([taxmap()])
#' @param data Dataset name, index, or a logical vector that indicates which dataset in
#'   `obj$data` to add columns to.
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return `numeric`
#'
#' @examples
#' # Get number of observations for each taxon in first dataset
#' n_obs(ex_taxmap)
#'
#' # Get number of observations in a specified data set
#' n_obs(ex_taxmap, "info")
#' n_obs(ex_taxmap, "abund")
#'
#' # Filter taxa using number of observations in the first table
#' filter_taxa(ex_taxmap, n_obs > 1)
#'
#' @family taxmap data functions
#'
#' @name n_obs
NULL


#' Count observation assigned in [taxmap()]
#'
#' Count observations for each taxon in a data set in a [taxmap()] object. This
#' includes observations for the specific taxon but NOT the observations of its
#' subtaxa. "Observations" in this sense are the items (for list/vectors) or
#' rows (for tables) in a dataset. By default, observations in the first data
#' set in the [taxmap()] object is used.  For example, if the data set is a
#' table, then a value of 3 for a taxon means that their are 3 rows in that
#' table assigned to that taxon.
#' \preformatted{
#' obj$n_obs_1(data)
#' n_obs_1(obj, data)}
#'
#' @param obj ([taxmap()])
#' @param data Dataset name, index, or a logical vector that indicates which dataset in
#'   `obj$data` to add columns to.
#' @param target DEPRECIATED. use "data" instead.
#'
#' @return `numeric`
#'
#' @examples
#' # Get number of observations for each taxon in first dataset
#' n_obs_1(ex_taxmap)
#'
#' # Get number of observations in a specified data set
#' n_obs_1(ex_taxmap, "info")
#' n_obs_1(ex_taxmap, "abund")
#'
#' # Filter taxa using number of observations in the first table
#' filter_taxa(ex_taxmap, n_obs_1 > 0)
#'
#' @family taxmap data functions
#'
#' @name n_obs_1
NULL


#' Get a data set from a taxmap object
#'
#' Get a data set from a taxmap object and complain if it does not
#' exist.
#'
#' @param obj A taxmap object
#' @param data Dataset name, index, or a logical vector that indicates which dataset in
#'   `obj$data` to add columns to.
#'
#' @examples
#' \dontrun{
#' # Get data set by name
#' get_dataset(ex_taxmap, "info")
#'
#' # Get data set by indeex_taxmap
#' get_dataset(ex_taxmap, 1)
#'
#' # Get data set by T/F vector
#' get_dataset(ex_taxmap, startsWith(names(ex_taxmap$data), "i"))
#'
#' }
#'
#' @name get_dataset
NULL
