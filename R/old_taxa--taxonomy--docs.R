#' Get taxon IDs
#'
#' Return the taxon IDs in a [taxonomy()] or [taxmap()] object.
#' They are in the order they appear in the edge list.
#' \preformatted{
#' obj$taxon_ids()
#' taxon_ids(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Return the taxon IDs for each taxon
#' taxon_ids(ex_taxmap)
#'
#' # Filter using taxon IDs
#' filter_taxa(ex_taxmap, ! taxon_ids %in% c("c", "d"))
#'
#' @name taxon_ids
NULL


#' Get taxon indexes
#'
#' Return the taxon indexes in a [taxonomy()] or [taxmap()] object.
#' They are the indexes of the edge list rows.
#' \preformatted{
#' obj$taxon_indexes()
#' taxon_indexes(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Return the indexes for each taxon
#' taxon_indexes(ex_taxmap)
#'
#' # Use in another function (stupid example; 1:5 would work too)
#' filter_taxa(ex_taxmap, taxon_indexes < 5)
#'
#' @name taxon_indexes
NULL


#' Get taxon names
#'
#' Return the taxon names in a [taxonomy()] or [taxmap()] object.
#' They are in the order they appear in the edge list.
#' \preformatted{
#' obj$taxon_names()
#' taxon_names(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Return the names for each taxon
#' taxon_names(ex_taxmap)
#'
#' # Filter by taxon name
#' filter_taxa(ex_taxmap, taxon_names == "Felidae", subtaxa = TRUE)
#'
#' @name taxon_names
NULL


#' Get taxon ranks
#'
#' Return the taxon ranks in a [taxonomy()] or [taxmap()] object.
#' They are in the order taxa appear in the edge list.
#' \preformatted{
#' obj$taxon_ranks()
#' taxon_ranks(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Get ranks for each taxon
#' taxon_ranks(ex_taxmap)
#'
#' # Filter by rank
#' filter_taxa(ex_taxmap, taxon_ranks == "family", supertaxa = TRUE)
#'
#' @name taxon_ranks
NULL


#' Get all supertaxa of a taxon
#'
#' Return data for supertaxa (i.e. all taxa the target
#' taxa are a part of) of each taxon in a [taxonomy()] or [taxmap()] object.
#' \preformatted{
#' obj$supertaxa(subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE,
#'   value = "taxon_indexes", na = FALSE)
#' supertaxa(obj, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE,
#'   value = "taxon_indexes", na = FALSE)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to find supertaxa for.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the
#'   supertaxa one rank above the target taxa. If `TRUE`, return all the
#'   supertaxa of every supertaxa, etc. Positive numbers indicate the number of
#'   recursions (i.e. number of ranks above the target taxon to return). `1` is
#'   equivalent to `FALSE`. Negative numbers are equivalent to `TRUE`.
#' @param simplify (`logical`) If `TRUE`, then combine all the results into a
#'   single vector of unique values.
#' @param include_input (`logical`) If `TRUE`, the input taxa are included in
#'   the output
#' @param value What data to return. Any result of [all_names()] can be used, but it
#'   usually only makes sense to use data that has an associated taxon id.
#' @param na (`logical`) If `TRUE`, return `NA` where information
#'   is not available.
#'
#' @return If `simplify = FALSE`, then a list of vectors are returned
#'   corresponding to the `subset` argument. If `simplify = TRUE`,
#'   then unique values are returned in a single vector.
#'
#' @family taxonomy indexing functions
#'
#' @name supertaxa
#'
#' @examples
#' # return the indexes for supertaxa for each taxon
#' supertaxa(ex_taxmap)
#'
#' # Only return data for some taxa using taxon indexes
#' supertaxa(ex_taxmap, subset = 1:3)
#'
#' # Only return data for some taxa using taxon ids
#' supertaxa(ex_taxmap, subset = c("d", "e"))
#'
#' # Only return data for some taxa using logical tests
#' supertaxa(ex_taxmap, subset = taxon_ranks == "species")
#'
#' # Only return supertaxa one level above
#' supertaxa(ex_taxmap, recursive = FALSE)
#'
#' # Only return supertaxa some number of ranks above
#' supertaxa(ex_taxmap, recursive = 2)
#'
#' # Return something besides taxon indexes
#' supertaxa(ex_taxmap, value = "taxon_names")
NULL


#' Apply function to supertaxa of each taxon
#'
#' Apply a function to the supertaxa for each taxon. This is similar
#' to using [supertaxa()] with [lapply()] or [sapply()].
#' \preformatted{
#' obj$supertaxa_apply(func, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE, value = "taxon_indexes",
#'   na = FALSE, ...)
#' supertaxa_apply(obj, func, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE, value = "taxon_indexes",
#'   na = FALSE, ....)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param func (`function`) The function to apply.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes of taxa to use.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the
#'   supertaxa one rank above the target taxa. If `TRUE`, return all the
#'   supertaxa of every supertaxa, etc. Positive numbers indicate the number of
#'   recursions (i.e. number of ranks above the target taxon to return). `1` is
#'   equivalent to `FALSE`. Negative numbers are equivalent to `TRUE`.
#' @param simplify (`logical`) If `TRUE`, then combine all the results into a
#'   single vector of unique values.
#' @param include_input (`logical`) If `TRUE`, the input taxa are included in
#'   the output
#' @param value What data to give to the function. Any result of
#'   `all_names(obj)` can be used, but it usually only makes sense to use data
#'   that has an associated taxon id.
#' @param na (`logical`) If `TRUE`, return `NA` where information
#'   is not available.
#' @param ... Extra arguments are passed to the function.
#'
#' @name supertaxa_apply
#'
#' @examples
#' # Get number of supertaxa that each taxon is contained in
#' supertaxa_apply(ex_taxmap, length)
#'
#' # Get classifications for each taxon
#' # Note; this can be done with `classifications()` easier
#' supertaxa_apply(ex_taxmap, paste, collapse = ";", include_input = TRUE,
#'                 value = "taxon_names")
#'
NULL


#' Get root taxa
#'
#' Return the root taxa for a [taxonomy()] or [taxmap()] object. Can also be used to
#' get the roots of a subset of taxa.
#' \preformatted{
#' obj$roots(subset = NULL, value = "taxon_indexes")
#' roots(obj, subset = NULL, value = "taxon_indexes")}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to find roots for.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param value What data to return. This is usually the name of column in a
#'   table in `obj$data`. Any result of `all_names(obj)` can be used, but it
#'   usually only makes sense to data that corresponds to taxa 1:1, such as
#'   [taxon_ranks()]. By default, taxon indexes are returned.
#'
#' @family taxonomy indexing functions
#'
#' @return `character`
#'
#' @examples
#' # Return indexes of root taxa
#' roots(ex_taxmap)
#'
#' # Return indexes for a subset of taxa
#' roots(ex_taxmap, subset = 2:17)
#'
#' # Return something besides taxon indexes
#' roots(ex_taxmap, value = "taxon_names")
#'
#' @name roots
NULL


#' Get "branch" taxa
#'
#' Return the "branch" taxa for a [taxonomy()] or [taxmap()] object. A branch is
#' anything that is not a root, stem, or leaf. Its the interior of the tree
#' after the first split starting from the roots. Can also be used to get the
#' branches of a subset of taxa.
#' \preformatted{
#' obj$branches(subset = NULL, value = "taxon_indexes")
#' branches(obj, subset = NULL, value = "taxon_indexes")}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes used to subset
#'   the tree prior to determining branches. Default: All taxa in `obj` will be
#'   used. Any variable name that appears in [all_names()] can be used as if it
#'   was a vector on its own. Note that branches are determined after the
#'   filtering, so a given taxon might be a branch on the unfiltered tree, but
#'   not a branch on the filtered tree.
#' @param value What data to return. This is usually the name of column in a
#'   table in `obj$data`. Any result of [all_names()] can be used, but it
#'   usually only makes sense to use data that corresponds to taxa 1:1, such as
#'   [taxon_ranks()]. By default, taxon indexes are returned.
#'
#' @family taxonomy indexing functions
#'
#' @return `character`
#'
#' @examples
#' # Return indexes of branch taxa
#' branches(ex_taxmap)
#'
#' # Return indexes for a subset of taxa
#' branches(ex_taxmap, subset = 2:17)
#' branches(ex_taxmap, subset = n_obs > 1)
#'
#' # Return something besides taxon indexes
#' branches(ex_taxmap, value = "taxon_names")
#'
#' @name branches
NULL


#' Get "internode" taxa
#'
#' Return the "internode" taxa for a [taxonomy()] or [taxmap()] object. An
#' internode is any taxon with a single immediate supertaxon and a single
#' immediate subtaxon. They can be removed from a tree without any loss of
#' information on the relative relationship between remaining taxa. Can also be
#' used to get the internodes of a subset of taxa.
#' \preformatted{
#' obj$internodes(subset = NULL, value = "taxon_indexes")
#' internodes(obj, subset = NULL, value = "taxon_indexes")}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes used to subset the tree prior to
#'   determining internodes. Default: All taxa in `obj` will be used. Any variable
#'   name that appears in [all_names()] can be used as if it was a vector on its
#'   own. Note that internodes are determined after the filtering, so a given
#'   taxon might be a internode on the unfiltered tree, but not a internode
#'   on the filtered tree.
#' @param value What data to return. This is usually the name of column in a
#'   table in `obj$data`. Any result of [all_names()] can be used, but it
#'   usually only makes sense to use data that corresponds to taxa 1:1, such as
#'   [taxon_ranks()]. By default, taxon indexes are returned.
#'
#' @family taxonomy indexing functions
#'
#' @return `character`
#'
#' @examples
#' \dontrun{
#'
#' # Return indexes of branch taxa
#' internodes(ex_taxmap)
#'
#' # Return indexes for a subset of taxa
#' internodes(ex_taxmap, subset = 2:17)
#' internodes(ex_taxmap, subset = n_obs > 1)
#'
#' # Return something besides taxon indexes
#' internodes(ex_taxmap, value = "taxon_names")
#'
#' }
#' @name internodes
NULL


#' Get subtaxa
#'
#' Return data for the subtaxa of each taxon in an [taxonomy()] or [taxmap()]
#' object.
#' \preformatted{
#' obj$subtaxa(subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE, value = "taxon_indexes")
#' subtaxa(obj, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE, value = "taxon_indexes")}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to find subtaxa for.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the subtaxa
#'   one rank below the target taxa. If `TRUE`, return all the subtaxa of every
#'   subtaxa, etc. Positive numbers indicate the number of ranks below the
#'   immediate subtaxa to return. `1` is equivalent to `FALSE`. Negative numbers
#'   are equivalent to `TRUE`. Since the algorithm is optimized for traversing
#'   all of large trees, `numeric` values greater than 0 for this option
#'   actually take slightly longer to compute than either TRUE or FALSE.
#' @param simplify (`logical`) If `TRUE`, then combine all the results
#'   into a single vector of unique values.
#' @param include_input (`logical`) If `TRUE`, the input taxa are
#'   included in the output
#' @param value What data to return. This is usually the name of column in a
#'   table in `obj$data`. Any result of [all_names()] can be used, but it
#'   usually only makes sense to data that corresponds to taxa 1:1, such as
#'   [taxon_ranks()]. By default, taxon indexes are returned.
#'
#' @return If `simplify = FALSE`, then a list of vectors are returned
#'   corresponding to the `target` argument. If `simplify = TRUE`,
#'   then the unique values are returned in a single vector.
#'
#' @family taxonomy indexing functions
#'
#' @name subtaxa
#'
#' @examples
#' # return the indexes for subtaxa for each taxon
#' subtaxa(ex_taxmap)
#'
#' # Only return data for some taxa using taxon indexes
#' subtaxa(ex_taxmap, subset = 1:3)
#'
#' # Only return data for some taxa using taxon ids
#' subtaxa(ex_taxmap, subset = c("d", "e"))
#'
#' # Only return data for some taxa using logical tests
#' subtaxa(ex_taxmap, subset = taxon_ranks == "genus")
#'
#' # Only return subtaxa one level below
#' subtaxa(ex_taxmap, recursive = FALSE)
#'
#' # Only return subtaxa some number of ranks below
#' subtaxa(ex_taxmap, recursive = 2)
#'
#' # Return something besides taxon indexes
#' subtaxa(ex_taxmap, value = "taxon_names")
NULL


#' Apply function to subtaxa of each taxon
#'
#' Apply a function to the subtaxa for each taxon. This is similar
#' to using [subtaxa()] with [lapply()] or [sapply()].
#' \preformatted{
#' obj$subtaxa_apply(func, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE, value = "taxon_indexes", ...)
#' subtaxa_apply(obj, func, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, include_input = FALSE, value = "taxon_indexes", ...)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param func (`function`) The function to apply.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to use.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the
#'   subtaxa one rank below the target taxa. If `TRUE`, return all the
#'   subtaxa of every subtaxa, etc. Positive numbers indicate the number of
#'   recursions (i.e. number of ranks below the target taxon to return). `1` is
#'   equivalent to `FALSE`. Negative numbers are equivalent to `TRUE`.
#' @param simplify (`logical`) If `TRUE`, then combine all the results into a
#'   single vector of unique values.
#' @param include_input (`logical`) If `TRUE`, the input taxa are included in
#'   the output
#' @param value What data to give to the function. Any result of
#'   `all_names(obj)` can be used, but it usually only makes sense to use data
#'   that has an associated taxon id.
#' @param ... Extra arguments are passed to the function.
#'
#' @name subtaxa_apply
#'
#' @examples
#' # Count number of subtaxa in each taxon
#' subtaxa_apply(ex_taxmap, length)
#'
#' # Paste all the subtaxon names for each taxon
#' subtaxa_apply(ex_taxmap, value = "taxon_names",
#'               recursive = FALSE, paste0, collapse = ", ")
NULL


#' Get stem taxa
#'
#' Return the stem taxa for a [taxonomy()] or a [taxmap()]
#' object. Stem taxa are all those from the roots to the first taxon with more
#' than one subtaxon.
#' \preformatted{
#' obj$stems(subset = NULL, simplify = FALSE,
#'   value = "taxon_indexes", exclude_leaves = FALSE)
#' stems(obj, subset = NULL, simplify = FALSE,
#'   value = "taxon_indexes", exclude_leaves = FALSE)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to find stems for.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param value What data to return. This is usually the name of column in a
#'   table in `obj$data`. Any result of `all_names(obj)` can be used, but it
#'   usually only makes sense to data that corresponds to taxa 1:1, such as
#'   [taxon_ranks()]. By default, taxon indexes are returned.
#' @param simplify (`logical`) If `TRUE`, then combine all the results
#'   into a single vector of unique values.
#' @param exclude_leaves (`logical`) If `TRUE`, the do not include
#'   taxa with no subtaxa.
#'
#' @return `character`
#'
#' @family taxonomy indexing functions
#'
#' @examples
#' # Return indexes of stem taxa
#' stems(ex_taxmap)
#'
#' # Return indexes for a subset of taxa
#' stems(ex_taxmap, subset = 2:17)
#'
#' # Return something besides taxon indexes
#' stems(ex_taxmap, value = "taxon_names")
#'
#' # Return a vector instead of a list
#' stems(ex_taxmap, value = "taxon_names", simplify = TRUE)
#'
#' @name stems
NULL


#' Get leaf taxa
#'
#' Return the leaf taxa for a [taxonomy()] or [taxmap()] object. Leaf taxa are taxa
#' with no subtaxa.
#' \preformatted{
#' obj$leaves(subset = NULL, recursive = TRUE, simplify = FALSE, value = "taxon_indexes")
#' leaves(obj, subset = NULL, recursive = TRUE, simplify = FALSE, value = "taxon_indexes")}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to find leaves for.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the
#'   leaves if they occur one rank below the target taxa. If `TRUE`, return all of the
#'   leaves for each taxon. Positive numbers indicate the number of
#'   recursions (i.e. number of ranks below the target taxon to return). `1` is
#'   equivalent to `FALSE`. Negative numbers are equivalent to `TRUE`.
#' @param simplify (`logical`) If `TRUE`, then combine all the results into a
#'   single vector of unique values.
#' @param value What data to return. This is usually the name of column in a
#'   table in `obj$data`. Any result of `all_names(obj)` can be used, but it
#'   usually only makes sense to data that corresponds to taxa 1:1, such as
#'   [taxon_ranks()]. By default, taxon indexes are returned.
#'
#' @return `character`
#'
#' @family taxonomy indexing functions
#'
#' @examples
#' # Return indexes of leaf taxa
#' leaves(ex_taxmap)
#'
#' # Return indexes for a subset of taxa
#' leaves(ex_taxmap, subset = 2:17)
#' leaves(ex_taxmap, subset = taxon_names == "Plantae")
#'
#' # Return something besides taxon indexes
#' leaves(ex_taxmap, value = "taxon_names")
#' leaves(ex_taxmap, subset = taxon_ranks == "genus", value = "taxon_names")
#'
#' # Return a vector of all unique values
#' leaves(ex_taxmap, value = "taxon_names", simplify = TRUE)
#'
#' # Only return leaves for their direct supertaxa
#' leaves(ex_taxmap, value = "taxon_names", recursive = FALSE)
#'
#' @name leaves
NULL


#' Apply function to leaves of each taxon
#'
#' Apply a function to the leaves of each taxon. This is similar
#' to using [leaves()] with [lapply()] or [sapply()].
#' \preformatted{
#' obj$leaves_apply(func, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, value = "taxon_indexes", ...)
#' leaves_apply(obj, func, subset = NULL, recursive = TRUE,
#'   simplify = FALSE, value = "taxon_indexes", ...)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object containing taxon
#'   information to be queried.
#' @param func (`function`) The function to apply.
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to use.
#'   Default: All taxa in `obj` will be used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own.
#' @param recursive (`logical` or `numeric`) If `FALSE`, only return the
#'   leaves if they occur one rank below the target taxa. If `TRUE`, return all of the
#'   leaves for each taxon. Positive numbers indicate the number of
#'   recursions (i.e. number of ranks below the target taxon to return). `1` is
#'   equivalent to `FALSE`. Negative numbers are equivalent to `TRUE`.
#' @param simplify (`logical`) If `TRUE`, then combine all the results into a
#'   single vector of unique values.
#' @param value What data to give to the function. Any result of
#'   `all_names(obj)` can be used, but it usually only makes sense to use data
#'   that has an associated taxon id.
#' @param ... Extra arguments are passed to the function `func`.
#'
#' @name leaves_apply
#'
#' @examples
#' # Count number of leaves under each taxon or its subtaxa
#' leaves_apply(ex_taxmap, length)
#'
#' # Count number of leaves under each taxon
#' leaves_apply(ex_taxmap, length, recursive = FALSE)
#'
#' # Converting output of leaves to upper case
#' leaves_apply(ex_taxmap, value = "taxon_names", toupper)
#'
#' # Passing arguments to the function
#' leaves_apply(ex_taxmap, value = "taxon_names", paste0, collapse = ", ")
#'
NULL


#' Get classifications of taxa
#'
#' Get character vector classifications of taxa in an object of type
#' [taxonomy()] or [taxmap()] composed of data associated with taxa. Each
#' classification is constructed by concatenating the data of the given taxon
#' and all of its supertaxa.
#' \preformatted{
#' obj$classifications(value = "taxon_names", sep = ";")
#' classifications(obj, value = "taxon_names", sep = ";")}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#' @param value What data to return. Any result of `all_names(obj)` can be used,
#'   but it usually only makes sense to data that corresponds to taxa 1:1, such
#'   as [taxon_ranks()]. By default, taxon indexes are returned.
#' @param sep (`character` of length 1) The character(s) to place between
#'   taxon IDs
#'
#' @return `character`
#'
#' @examples
#' # Defualt settings returns taxon names separated by ;
#' classifications(ex_taxmap)
#'
#' # Other values can be returned besides taxon names
#' classifications(ex_taxmap, value = "taxon_ids")
#'
#' # The separator can also be changed
#' classifications(ex_taxmap, value = "taxon_ranks", sep = "||")
#'
#' @family taxonomy data functions
#'
#' @name classifications
NULL


#' Get ID classifications of taxa
#'
#' Get classification strings of taxa in an object of type [taxonomy()] or [taxmap()]
#' composed of taxon IDs. Each classification is constructed by concatenating
#' the taxon ids of the given taxon and its supertaxa.
#' \preformatted{
#' obj$id_classifications(sep = ";")
#' id_classifications(obj, sep = ";")}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#' @param sep (`character` of length 1) The character(s) to place between
#'   taxon IDs
#'
#' @return `character`
#'
#' @examples
#' # Get classifications of IDs for each taxon
#' id_classifications(ex_taxmap)
#'
#' # Use a different seperator
#' id_classifications(ex_taxmap, sep = '|')
#'
#' @family taxonomy data functions
#'
#' @name id_classifications
NULL


#' Get number of supertaxa
#'
#' Get number of supertaxa for each taxon in an object of type
#' [taxonomy()] or [taxmap()].
#' \preformatted{
#' obj$n_supertaxa()
#' n_supertaxa(obj)}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#'
#' @return \code{numeric}
#'
#' @examples
#' # Count number of supertaxa that contain each taxon
#' n_supertaxa(ex_taxmap)
#'
#' # Filter taxa based on the number of supertaxa
#' #  (this command removes all root taxa)
#' filter_taxa(ex_taxmap, n_supertaxa > 0)
#'
#' @family taxonomy data functions
#'
#' @name n_supertaxa
NULL


#' Get number of supertaxa
#'
#' Get number of immediate supertaxa (i.e. not supertaxa of supertaxa, etc) for
#' each taxon in an object of type [taxonomy()] or [taxmap()]. This should
#' always be either 1 or 0.
#' \preformatted{
#' obj$n_supertaxa_1()
#' n_supertaxa_1(obj)}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#'
#' @return \code{numeric}
#'
#' @examples
#' # Test for the presence of supertaxa containing each taxon
#' n_supertaxa_1(ex_taxmap)
#'
#' # Filter taxa based on the presence of supertaxa
#' #  (this command removes all root taxa)
#' filter_taxa(ex_taxmap, n_supertaxa_1 > 0)
#'
#' @family taxonomy data functions
#'
#' @name n_supertaxa_1
NULL


#' Get number of subtaxa
#'
#' Get number of subtaxa for each taxon in an object of type
#' [taxonomy()] or [taxmap()]
#' \preformatted{
#' obj$n_subtaxa()
#' n_subtaxa(obj)}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#'
#' @return \code{numeric}
#'
#' @examples
#' # Count number of subtaxa within each taxon
#' n_subtaxa(ex_taxmap)
#'
#' # Filter taxa based on number of subtaxa
#' #  (this command removed all leaves or "tips" of the tree)
#' filter_taxa(ex_taxmap, n_subtaxa > 0)
#'
#' @family taxonomy data functions
#'
#' @name n_subtaxa
NULL


#' Get number of subtaxa
#'
#' Get number of subtaxa for each taxon in an object of type
#' [taxonomy()] or [taxmap()], not including subtaxa of subtaxa etc. This does not
#' include subtaxa assigned to subtaxa.
#' \preformatted{
#' obj$n_subtaxa_1()
#' n_subtaxa_1(obj)}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#'
#' @return \code{numeric}
#'
#' @examples
#' # Count number of immediate subtaxa in each taxon
#' n_subtaxa_1(ex_taxmap)
#'
#' # Filter taxa based on number of subtaxa
#' #  (this command removed all leaves or "tips" of the tree)
#' filter_taxa(ex_taxmap, n_subtaxa_1 > 0)
#'
#' @family taxonomy data functions
#'
#' @name n_subtaxa_1
NULL


#' Get number of leaves
#'
#' Get number of leaves for each taxon in an object of type
#' [taxonomy()] or [taxmap()]
#' \preformatted{
#' obj$n_leaves()
#' n_leaves(obj)}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#'
#' @return \code{numeric}
#'
#' @examples
#' # Get number of leaves for each taxon
#' n_leaves(ex_taxmap)
#'
#' # Filter taxa based on number of leaves
#' filter_taxa(ex_taxmap, n_leaves > 0)
#'
#' @family taxonomy data functions
#'
#' @name n_leaves
NULL


#' Get number of leaves
#'
#' Get number of leaves for each taxon in an object of type
#' [taxonomy()] or [taxmap()], not including leaves of subtaxa etc.
#' \preformatted{
#' obj$n_leaves_1()
#' n_leaves_1(obj)}
#'
#' @param obj ([taxonomy()] or [taxmap()])
#'
#' @return \code{numeric}
#'
#' @examples
#' # Get number of leaves for each taxon
#' n_leaves_1(ex_taxmap)
#'
#' # Filter taxa based on number of leaves
#' filter_taxa(ex_taxmap, n_leaves_1 > 0)
#'
#' @family taxonomy data functions
#'
#' @name n_leaves_1
NULL


#' Return names of data in [taxonomy()] or [taxmap()]
#'
#' Return the names of data that can be used with functions in the taxa
#' package that use [non-standard evaluation](http://adv-r.had.co.nz/Computing-on-the-language.html) (NSE),
#' like [filter_taxa()].
#' \preformatted{
#' obj$all_names(tables = TRUE, funcs = TRUE,
#'   others = TRUE, warn = FALSE)
#' all_names(obj, tables = TRUE, funcs = TRUE,
#'   others = TRUE, warn = FALSE)}
#'
#' @param obj ([taxonomy()] or [taxmap()]) The object containing
#'   taxon information to be queried.
#' @param tables This option only applies to [taxmap()] objects. If `TRUE`,
#'   include the names of columns of tables in `obj$data`
#' @param funcs This option only applies to [taxmap()] objects. If `TRUE`,
#'   include the names of user-definable functions in `obj$funcs`.
#' @param others This option only applies to [taxmap()] objects. If `TRUE`,
#'   include the names of data in `obj$data` besides tables.
#' @param builtin_funcs This option only applies to [taxmap()] objects. If
#'   `TRUE`, include functions like [n_supertaxa()] that provide information for
#'   each taxon.
#' @param warn option only applies to [taxmap()] objects. If `TRUE`, warn if
#'   there are duplicate names. Duplicate names make it unclear what data is
#'   being referred to.
#'
#' @return `character`
#'
#' @examples
#' # Get the names of all data accesible by non-standard evaluation
#' all_names(ex_taxmap)
#'
#' # Dont include the names of automatically included functions.
#' all_names(ex_taxmap, builtin_funcs = FALSE)
#'
#' @family NSE helpers
#'
#' @name all_names
NULL


#' Get names of data used in expressions
#'
#' Get names of available data used in expressions. This is used to find data
#' for use with [non-standard evaluation](http://adv-r.had.co.nz/Computing-on-the-language.html) (NSE) in
#' functions like [filter_taxa()]. Expressions are not evaluated and do not need
#' to make sense.
#' \preformatted{
#' obj$names_used(...)}
#'
#' @param obj a [taxonomy()] or [taxmap()] object
#' @param ... One or more expressions
#'
#' @return Named `character`
#'
#' @examples
#' ex_taxmap$names_used(n_legs + dangerous == invalid_expression)
#'
#' @family NSE helpers
#'
#' @name names_used
#' @keywords internal
NULL


#' Get data in a taxmap object by name
#'
#' Given a vector of names, return a list of data (usually lists/vectors)
#' contained in a [taxonomy()] or [taxmap()] object. Each item will be named by
#' taxon ids when possible.
#' \preformatted{
#' obj$get_data(name = NULL, ...)
#' get_data(obj, name = NULL, ...)}
#'
#' @param obj A [taxonomy()] or [taxmap()]  object
#' @param name (`character`) Names of data to return. If not supplied, return
#'   all data listed in [all_names()].
#' @param ... Passed to [all_names()]. Used to filter what kind of data is
#'   returned (e.g. columns in tables or function output?) if `name` is not
#'   supplied or what kinds are allowed if `name` is supplied.
#'
#' @return `list` of vectors or lists. Each vector or list will be named by
#'   associated taxon ids if possible.
#'
#' @examples
#' # Get specific values
#' get_data(ex_taxmap, c("reaction", "n_legs", "taxon_ranks"))
#'
#' # Get all values
#' get_data(ex_taxmap)
#'
#' @family NSE helpers
#'
#' @name get_data
NULL


#' Get data in a taxonomy or taxmap object by name
#'
#' Given a vector of names, return a  table of the indicated data
#' contained in a [taxonomy()] or [taxmap()] object.
#' \preformatted{
#' obj$get_data_frame(name = NULL, ...)
#' get_data_frame(obj, name = NULL, ...)}
#'
#' Note: This function will not work with variables in datasets in [taxmap()]
#' objects unless their rows correspond 1:1 with all taxa.
#'
#' @param obj A [taxonomy()] or [taxmap()]  object
#' @param name (`character`) Names of data to return. If not supplied, return
#'   all data listed in [all_names()].
#' @param ... Passed to [all_names()]. Used to filter what kind of data is
#'   returned (e.g. columns in tables or function output?) if `name` is not
#'   supplied or what kinds are allowed if `name` is supplied.
#'
#' @return `data.frame`
#'
#' @examples
#' # Get specific values
#' get_data_frame(ex_taxonomy, c("taxon_names", "taxon_indexes", "is_stem"))
#'
#' # Get all values
#' get_data_frame(ex_taxonomy)
#'
#' @family accessors
#'
#' @name get_data_frame
NULL


#' Get values of data used in expressions
#'
#' Get values available for
#' [non-standard evaluation](http://adv-r.had.co.nz/Computing-on-the-language.html)
#' in a [taxonomy()] or [taxmap()] object used in expressions. Expressions are
#' not evaluated and do not need to make sense.
#' \preformatted{
#' obj$data_used(...)}
#'
#' @param obj a [taxonomy()] or [taxmap()] object
#' @param ... One or more expressions
#'
#' @return `list`
#'
#' @examples
#' # Get values for variables names used in expressions
#' ex_taxmap$data_used(n_legs + dangerous == invalid_expression)
#' ex_taxmap$data_used(length(unique(taxon_names)))
#'
#' @family NSE helpers
#'
#' @name data_used
#' @keywords internal
NULL


#' Filter taxa with a list of conditions
#'
#' Filter taxa in a [taxonomy()] or [taxmap()] object with a series of
#' conditions. Any variable name that appears in [all_names()] can be used
#' as if it was a vector on its own. See [dplyr::filter()] for the inspiration
#' for this function and more information. Calling the function using the
#' `obj$filter_taxa(...)` style edits "obj" in place, unlike most R functions.
#' However, calling the function using the `filter_taxa(obj, ...)` imitates R's
#' traditional copy-on-modify semantics, so "obj" would not be changed; instead
#' a changed version would be returned, like most R functions.
#' \preformatted{
#' filter_taxa(obj, ..., subtaxa = FALSE, supertaxa = FALSE,
#'   drop_obs = TRUE, reassign_obs = TRUE, reassign_taxa = TRUE,
#'   invert = FALSE, keep_order = TRUE)
#' obj$filter_taxa(..., subtaxa = FALSE, supertaxa = FALSE,
#'   drop_obs = TRUE, reassign_obs = TRUE, reassign_taxa = TRUE,
#'   invert = FALSE, keep_order = TRUE)}
#'
#' @param obj An object of class [taxonomy()] or [taxmap()]
#' @param ... One or more filtering conditions. Any variable name that appears
#'   in [all_names()] can be used as if it was a vector on its own. Each
#'   filtering condition must resolve to one of three things:
#'   * `character`: One or more taxon IDs contained in `obj$edge_list$to`
#'   * `integer`: One or more row indexes of `obj$edge_list`
#'   * `logical`: A `TRUE`/`FALSE` vector of length equal to the number of rows
#'   in `obj$edge_list`
#'   * `NULL`: ignored
#' @param subtaxa (`logical` or `numeric` of length 1) If `TRUE`, include
#'   subtaxa of taxa passing the filter. Positive numbers indicate the number of
#'   ranks below the target taxa to return. `0` is equivalent to `FALSE`.
#'   Negative numbers are equivalent to `TRUE`.
#' @param supertaxa (`logical`  or `numeric` of length 1) If `TRUE`, include
#'   supertaxa of taxa passing the filter. Positive numbers indicate the number
#'   of ranks above the target taxa to return. `0` is equivalent to `FALSE`.
#'   Negative numbers are equivalent to `TRUE`.
#' @param drop_obs (`logical`)  This option only applies to [taxmap()] objects.
#'   If `FALSE`, include observations (i.e. user-defined data in `obj$data`)
#'   even if the taxon they are assigned to is filtered out. Observations
#'   assigned to removed taxa will be assigned to \code{NA}. This option can be
#'   either simply `TRUE`/`FALSE`, meaning that all data sets will be treated
#'   the same, or a logical vector can be supplied with names corresponding one
#'   or more data sets in `obj$data`. For example, `c(abundance = FALSE, stats =
#'   TRUE)` would include observations whose taxon was filtered out in
#'   `obj$data$abundance`, but not in `obj$data$stats`. See the `reassign_obs`
#'   option below for further complications.
#' @param reassign_obs (`logical` of length 1) This option only applies to
#'   [taxmap()] objects. If `TRUE`, observations (i.e. user-defined data in
#'   `obj$data`) assigned to removed taxa will be reassigned to the closest
#'   supertaxon that passed the filter. If there are no supertaxa of such an
#'   observation that passed the filter, they will be filtered out if `drop_obs`
#'   is `TRUE`. This option can be either simply `TRUE`/`FALSE`, meaning that
#'   all data sets will be treated the same, or a logical vector can be supplied
#'   with names corresponding one or more data sets in `obj$data`. For example,
#'   `c(abundance = TRUE, stats = FALSE)` would reassign observations in
#'   `obj$data$abundance`, but not in `obj$data$stats`.
#' @param reassign_taxa (`logical` of length 1) If `TRUE`, subtaxa of removed
#'   taxa will be reassigned to the closest supertaxon that passed the filter.
#'   This is useful for removing intermediate levels of a taxonomy.
#' @param invert (`logical` of length 1) If `TRUE`, do NOT include the
#'   selection. This is different than just replacing a `==` with a `!=` because
#'   this option negates the selection after taking into account the `subtaxa`
#'   and `supertaxa` options. This is useful for removing a taxon and all its
#'   subtaxa for example.
#' @param keep_order (`logical` of length 1) If `TRUE`, keep relative order of
#'   taxa not filtered out. For example, the result of `filter_taxa(ex_taxmap,
#'   1:3)` and `filter_taxa(ex_taxmap, 3:1)` would be the same. Does not affect
#'   dataset order, only taxon order. This is useful for maintaining order
#'   correspondence with a dataset that has one value per taxon.
#'
#' @return An object of type [taxonomy()] or [taxmap()]
#'
#' @examples
#' # Filter by index
#' filter_taxa(ex_taxmap, 1:3)
#'
#' # Filter by taxon ID
#' filter_taxa(ex_taxmap, c("b", "c", "d"))
#'
#' # Fiter by TRUE/FALSE
#' filter_taxa(ex_taxmap, taxon_names == "Plantae", subtaxa = TRUE)
#' filter_taxa(ex_taxmap, n_obs > 3)
#' filter_taxa(ex_taxmap, ! taxon_ranks %in% c("species", "genus"))
#' filter_taxa(ex_taxmap, taxon_ranks == "genus", n_obs > 1)
#'
#' # Filter by an observation characteristic
#' dangerous_taxa <- sapply(ex_taxmap$obs("info"),
#'                          function(i) any(ex_taxmap$data$info$dangerous[i]))
#' filter_taxa(ex_taxmap, dangerous_taxa)
#'
#' # Include supertaxa
#' filter_taxa(ex_taxmap, 12, supertaxa = TRUE)
#' filter_taxa(ex_taxmap, 12, supertaxa = 2)
#'
#' # Include subtaxa
#' filter_taxa(ex_taxmap, 1, subtaxa = TRUE)
#' filter_taxa(ex_taxmap, 1, subtaxa = 2)
#'
#' # Dont remove rows in user-defined data corresponding to removed taxa
#' filter_taxa(ex_taxmap, 2, drop_obs = FALSE)
#' filter_taxa(ex_taxmap, 2, drop_obs = c(info = FALSE))
#'
#' # Remove a taxon and it subtaxa
#' filter_taxa(ex_taxmap, taxon_names == "Mammalia",
#'             subtaxa = TRUE, invert = TRUE)
#'
#' @family taxmap manipulation functions
#'
#' @name filter_taxa
NULL


#' Sort the edge list of [taxmap()] objects
#'
#' Sort the edge list and taxon list in [taxonomy()] or [taxmap()] objects. See
#' [dplyr::arrange()] for the inspiration for this function and more
#' information. Calling the function using the `obj$arrange_taxa(...)` style
#' edits "obj" in place, unlike most R functions. However, calling the function
#' using the `arrange_taxa(obj, ...)` imitates R's traditional copy-on-modify
#' semantics, so "obj" would not be changed; instead a changed version would be
#' returned, like most R functions.
#' \preformatted{
#' obj$arrange_taxa(...)
#' arrange_taxa(obj, ...)}
#'
#' @param obj [taxonomy()] or [taxmap()]
#' @param ... One or more expressions (e.g. column names) to sort on. Any
#'   variable name that appears in [all_names()] can be used as if it was a
#'   vector on its own.
#'
#' @return An object of type [taxonomy()] or [taxmap()]
#'
#' @examples
#' # Sort taxa in ascending order
#' arrange_taxa(ex_taxmap, taxon_names)
#'
#' # Sort taxa in decending order
#' arrange_taxa(ex_taxmap, desc(taxon_names))
#'
#' # Sort using an expression. List genera first.
#' arrange_taxa(ex_taxmap, taxon_ranks != "genus")
#'
#' @family taxmap manipulation functions
#'
#' @name arrange_taxa
NULL

#' Sample n taxa from [taxonomy()] or [taxmap()]
#'
#' Randomly sample some number of taxa from a [taxonomy()] or [taxmap()] object.
#' Weights can be specified for taxa or the observations assigned to them. See
#' [dplyr::sample_n()] for the inspiration for this function.
#' \preformatted{
#' obj$sample_n_taxa(size, taxon_weight = NULL,
#'   obs_weight = NULL, obs_target = NULL,
#'   use_subtaxa = TRUE, collapse_func = mean, ...)
#' sample_n_taxa(obj, size, taxon_weight = NULL,
#'   obs_weight = NULL, obs_target = NULL,
#'   use_subtaxa = TRUE, collapse_func = mean, ...)}
#'
#' @param obj ([taxonomy()] or [taxmap()]) The object to sample from.
#' @param size (`numeric` of length 1) The number of taxa to sample.
#' @param taxon_weight (`numeric`) Non-negative sampling weights of each
#'   taxon. If `obs_weight` is also specified, the two weights are
#' multiplied (after `obs_weight` for each taxon is calculated).
#' @param obs_weight (`numeric`)  This option only applies to [taxmap()]
#'   objects. Sampling weights of each observation. The weights for each
#'   observation assigned to a given taxon are supplied to `collapse_func` to
#'   get the taxon weight. If `use_subtaxa` is `TRUE` then the observations
#'   assigned to every subtaxa are also used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own. If
#'   `taxon_weight` is also specified, the two weights are multiplied (after
#'   `obs_weight` for each observation is calculated). `obs_target` must be used
#'   with this option.
#' @param obs_target (`character` of length 1)  This option only applies to
#'   [taxmap()] objects. The name of the data set in `obj$data` that values in
#'   `obs_weight` corresponds to. Must be used when `obs_weight` is used.
#' @param use_subtaxa (`logical` or `numeric` of length 1) Affects how the
#'   `obs_weight` option is used. If `TRUE`, the weights for each taxon in an
#'   observation's classification are multiplied to get the observation weight.
#'   If `FALSE` just the taxonomic level the observation is assign to it
#'   considered. Positive numbers indicate the number of ranks below the each
#'   taxon to use. `0` is equivalent to `FALSE`. Negative numbers are equivalent
#'   to `TRUE`.
#' @param collapse_func (`function` of length 1) If `taxon_weight` is used and
#'   `supertaxa` is `TRUE`, the weights for each taxon in an observation's
#'   classification are supplied to `collapse_func` to get the observation
#'   weight. This function should take  numeric vector and return a single
#'   number.
#' @param ... Additional options are passed to [filter_taxa()].
#'
#' @return An object of type [taxonomy()] or [taxmap()]
#'
#' @examples
#' # Randomly sample three taxa
#' sample_n_taxa(ex_taxmap, 3)
#'
#' # Include supertaxa
#' sample_n_taxa(ex_taxmap, 3, supertaxa = TRUE)
#'
#' # Include subtaxa
#' sample_n_taxa(ex_taxmap, 1, subtaxa = TRUE)
#'
#' # Sample some taxa more often then others
#' sample_n_taxa(ex_taxmap, 3, supertaxa = TRUE,
#'               obs_weight = n_legs, obs_target = "info")
#'
#' @family taxmap manipulation functions
#'
#' @name sample_n_taxa
NULL


#' Sample a proportion of taxa from [taxonomy()] or [taxmap()]
#'
#' Randomly sample some proportion of taxa from a [taxonomy()] or [taxmap()]
#' object. Weights can be specified for taxa or the observations assigned to
#' them. See
#' [dplyr::sample_frac()] for the inspiration for this function.
#' \preformatted{
#' obj$sample_frac_taxa(size, taxon_weight = NULL,
#'   obs_weight = NULL, obs_target = NULL,
#'   use_subtaxa = TRUE, collapse_func = mean, ...)
#' sample_frac_taxa(obj, size, taxon_weight = NULL,
#'   obs_weight = NULL, obs_target = NULL,
#'   use_subtaxa = TRUE, collapse_func = mean, ...)}
#'
#' @param obj ([taxonomy()] or [taxmap()]) The object to sample from.
#' @param size (`numeric` of length 1) The proportion of taxa to sample.
#' @param taxon_weight (`numeric`) Non-negative sampling weights of each
#'   taxon. If `obs_weight` is also specified, the two weights are
#' multiplied (after `obs_weight` for each taxon is calculated).
#' @param obs_weight (`numeric`)  This option only applies to [taxmap()]
#'   objects. Sampling weights of each observation. The weights for each
#'   observation assigned to a given taxon are supplied to `collapse_func` to
#'   get the taxon weight. If `use_subtaxa` is `TRUE` then the observations
#'   assigned to every subtaxa are also used. Any variable name that appears in
#'   [all_names()] can be used as if it was a vector on its own. If
#'   `taxon_weight` is also specified, the two weights are multiplied (after
#'   `obs_weight` for each observation is calculated). `obs_target` must be used
#'   with this option.
#' @param obs_target (`character` of length 1)  This option only applies to
#'   [taxmap()] objects. The name of the data set in `obj$data` that values in
#'   `obs_weight` corresponds to. Must be used when `obs_weight` is used.
#' @param use_subtaxa (`logical` or `numeric` of length 1) Affects how the
#'   `obs_weight` option is used. If `TRUE`, the weights for each taxon in an
#'   observation's classification are multiplied to get the observation weight.
#'   If `TRUE` just the taxonomic level the observation is assign to it
#'   considered. Positive numbers indicate the number of ranks below the target
#'   taxa to return. `0` is equivalent to `FALSE`. Negative numbers are
#'   equivalent to `TRUE`.
#' @param collapse_func (`function` of length 1) If `taxon_weight` is
#'   used and `supertaxa` is `TRUE`, the weights for each taxon in an
#'   observation's classification are supplied to `collapse_func` to get
#'   the observation weight. This function should take  numeric vector and
#'   return a single number.
#' @param ... Additional options are passed to [filter_taxa()].
#'
#' @return An object of type [taxonomy()] or [taxmap()]
#'
#'
#' @examples
#' # sample half of the taxa
#' sample_frac_taxa(ex_taxmap, 0.5, supertaxa = TRUE)
#'
#' @family taxmap manipulation functions
#'
#' @name sample_frac_taxa
NULL


#' Test if taxa are roots
#'
#' Test if taxa are roots in a [taxonomy()] or [taxmap()] object. Roots are taxa
#' without supertaxa, typically things like "Bacteria", or "Life".
#' \preformatted{
#' obj$is_root()
#' is_root(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @return A `logical` of length equal to the number of taxa.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Test for which taxon IDs correspond to roots
#' is_root(ex_taxmap)
#'
#' # Filter out roots
#' filter_taxa(ex_taxmap, ! is_root)
#'
#' @name is_root
NULL


#' Test if taxa are "internodes"
#'
#' Test if taxa are "internodes" in a [taxonomy()] or [taxmap()] object.  An
#' internode is any taxon with a single immediate supertaxon and a single immediate
#' subtaxon. They can be removed from a tree without any loss of information on
#' the relative relationship between remaining taxa.
#' \preformatted{
#' obj$is_internode()
#' is_internode(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @return A `logical` of length equal to the number of taxa.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Test for which taxon IDs correspond to internodes
#' is_internode(ex_taxmap)
#'
#' # Filter out internodes
#' filter_taxa(ex_taxmap, ! is_internode)
#'
#' @name is_internode
NULL


#' Test if taxa are stems
#'
#' Test if taxa are stems in a [taxonomy()] or [taxmap()] object. Stems are taxa
#' from the [roots()] taxa to the first taxon with more than one subtaxon. These
#' can usually be filtered out of the taxonomy without removing any information
#' on how the remaining taxa are related.
#' \preformatted{
#' obj$is_stem()
#' is_stem(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @return A `logical` of length equal to the number of taxa.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Test which taxon IDs correspond to stems
#' is_stem(ex_taxmap)
#'
#' # Filter out stems
#' filter_taxa(ex_taxmap, ! is_stem)
#'
#' @name is_stem
NULL


#' Test if taxa are branches
#'
#' Test if taxa are branches in a [taxonomy()] or [taxmap()] object. Branches
#' are taxa in the interior of the tree that are not [roots()], [stems()], or
#' [leaves()].
#' \preformatted{
#' obj$is_branch()
#' is_branch(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @return A `logical` of length equal to the number of taxa.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Test which taxon IDs correspond to branches
#' is_branch(ex_taxmap)
#'
#' # Filter out branches
#' filter_taxa(ex_taxmap, ! is_branch)
#'
#' @name is_branch
NULL


#' Test if taxa are leaves
#'
#' Test if taxa are leaves in a [taxonomy()] or [taxmap()] object. Leaves are taxa
#' without subtaxa, typically species.
#' \preformatted{
#' obj$is_leaf()
#' is_leaf(obj)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#'
#' @return A `logical` of length equal to the number of taxa.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Test which taxon IDs correspond to leaves
#' is_leaf(ex_taxmap)
#'
#' # Filter out leaves
#' filter_taxa(ex_taxmap, ! is_leaf)
#'
#' @name is_leaf
NULL


#' Create a mapping between two variables
#'
#' Creates a named vector that maps the values of two variables associated with
#' taxa in a [taxonomy()] or [taxmap()] object. Both values must be named by
#' taxon ids.
#' \preformatted{
#' obj$map_data(from, to, warn = TRUE)
#' map_data(obj, from, to, warn = TRUE)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#' @param from The value used to name the output. There will be one output value
#'   for each value in `from`. Any variable that appears in [all_names()] can be
#'   used as if it was a variable on its own.
#' @param to The value returned in the output. Any variable that appears in
#'   [all_names()] can be used as if it was a variable on its own.
#' @param warn If `TRUE`, issue a warning if there are multiple unique values of
#'   `to` for each value of `from`.
#'
#' @return A vector of `to` values named by values in `from`.
#'
#' @family taxonomy data functions
#'
#' @examples
#' # Mapping between two variables in `all_names(ex_taxmap)`
#' map_data(ex_taxmap, from = taxon_names, to = n_legs > 0)
#'
#' # Mapping with external variables
#' x = c("d" = "looks like a cat", "h" = "big scary cats",
#'       "i" = "smaller cats", "m" = "might eat you", "n" = "Meow! (Feed me!)")
#' map_data(ex_taxmap, from = taxon_names, to = x)
#'
#' @name map_data
NULL


#' Create a mapping without NSE
#'
#' Creates a named vector that maps the values of two variables associated with
#' taxa in a [taxonomy()] or [taxmap()] object without using Non-Standard
#' Evaluation (NSE). Both values must be named by taxon ids. This is the same as
#' [map_data()] without NSE and can be useful in some odd cases where NSE fails
#' to work as expected.
#' \preformatted{
#' obj$map_data(from, to)
#' map_data(obj, from, to)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#' @param from The value used to name the output. There will be one output value
#'   for each value in `from`.
#' @param to The value returned in the output.
#'
#' @return A vector of `to` values named by values in `from`.
#'
#' @family taxonomy data functions
#'
#' @examples
#' x = c("d" = "looks like a cat", "h" = "big scary cats",
#'       "i" = "smaller cats", "m" = "might eat you", "n" = "Meow! (Feed me!)")
#' map_data_(ex_taxmap, from = ex_taxmap$taxon_names(), to = x)
#'
#' @name map_data_
NULL


#' Replace taxon ids
#'
#' Replace taxon ids in a [taxmap()] or [taxonomy()] object.
#' \preformatted{
#' obj$replace_taxon_ids(new_ids)
#' replace_taxon_ids(obj, new_ids)}
#'
#' @param obj The [taxonomy()] or [taxmap()] object.
#' @param new_ids A vector of new ids, one per taxon. They must be unique and in
#'   the same order as the corresponding ids in `obj$taxon_ids()`.
#'
#' @return A [taxonomy()] or [taxmap()] object with new taxon ids
#' @name replace_taxon_ids
#' @examples
#' # Replace taxon IDs with numbers
#' replace_taxon_ids(ex_taxmap, seq_len(length(ex_taxmap$taxa)))
#'
#' # Make taxon IDs capital letters
#' replace_taxon_ids(ex_taxmap, toupper(taxon_ids(ex_taxmap)))
#'
NULL


#' Remove redundant parts of taxon names
#'
#' Remove the names of parent taxa in the beginning of their children's names in a \code{taxonomy} or \code{taxmap} object.
#' This is useful for removing genus names in species binomials.
#' \preformatted{
#' obj$remove_redundant_names()
#' remove_redundant_names(obj)}
#'
#' @param obj A \code{taxonomy} or \code{taxmap} object
#'
#' @return A \code{taxonomy} or \code{taxmap} object
#' @name remove_redundant_names
#' @examples
#' # Remove genus named from species taxa
#' species_data <- c("Carnivora;Felidae;Panthera;Panthera leo",
#'                   "Carnivora;Felidae;Panthera;Panthera tigris",
#'                   "Carnivora;Ursidae;Ursus;Ursus americanus")
#' obj <-  parse_tax_data(species_data, class_sep = ";")
#' remove_redundant_names(obj)
NULL

#' Convert taxonomy info to a table
#'
#' Convert per-taxon information, like taxon names, to a table of taxa (rows) by
#' ranks (columns).
#'
#' @param obj A \code{taxonomy} or \code{taxmap} object
#' @param subset Taxon IDs, TRUE/FALSE vector, or taxon indexes to find
#'   supertaxa for. Default: All leaves will be used. Any variable name that
#'   appears in [all_names()] can be used as if it was a vector on its own.
#' @param value What data to return. Default is taxon names. Any result of
#'   [all_names()] can be used, but it usually only makes sense to use data with
#'   one value per taxon, like taxon names.
#' @param use_ranks Which ranks to use. Must be one of the following:
#' * `NULL` (the default): If there is rank information, use the ranks that
#' appear in the lineage with the most ranks. Otherwise, assume the number of
#' supertaxa corresponds to rank and use placeholders for the rank column names
#' in the output.
#' * `TRUE`: Use the ranks that appear in the lineage with the most ranks. An
#' error will occur if no rank information is available.
#' * `FALSE`: Assume the number of supertaxa corresponds to rank and use
#' placeholders for the rank column names in the output. Do not use included
#' rank information.
#' * `character`: The names of the ranks to use. Requires included rank information.
#' * `numeric`: The "depth" of the ranks to use. These are equal to `n_supertaxa` + 1.
#' @param add_id_col If `TRUE`, include a taxon ID column.
#'
#' @return A tibble of taxa (rows) by ranks (columns).
#'
#' @examples
#' # Make a table of taxon names
#' taxonomy_table(ex_taxmap)
#'
#' # Use a differnt value
#' taxonomy_table(ex_taxmap, value = "taxon_ids")
#'
#' # Return a subset of taxa
#' taxonomy_table(ex_taxmap, subset = taxon_ranks == "genus")
#'
#' # Use arbitrary ranks names based on depth
#' taxonomy_table(ex_taxmap, use_ranks = FALSE)
#'
#' @name taxonomy_table
NULL


#' Print a text tree
#'
#' Print a text-based tree of a [taxonomy()] or [taxmap()] object.
#'
#' @param obj A \code{taxonomy} or \code{taxmap} object
#' @param value What data to return. Default is taxon names. Any result of
#'   [all_names()] can be used, but it usually only makes sense to use data with
#'   one value per taxon, like taxon names.
#'
#' @examples
#' print_tree(ex_taxmap)
#'
#' @name print_tree
NULL

