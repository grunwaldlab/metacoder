#' An example taxmap object
#'
#' An example taxmap object built from the ground up. Typically, data stored in
#' taxmap would be parsed from an input file, but this data set is just for
#' demonstration purposes.
#'
#' @name ex_taxmap
#' @format A [taxmap()] object.
#' @source Created from the example code in the [taxmap()]
#'   documentation.
#' @family taxa-datasets
#' @keywords data
NULL

#' An example Taxonomy object
#'
#' An example Taxonomy object built from the ground up.
#'
#' @name ex_taxonomy
#' @format A [taxonomy()] object.
#' @source Created from the example code in the [taxonomy()]
#'   documentation.
#' @family taxa-datasets
#' @keywords data
NULL

#' An example Hierarchy object
#'
#' An example Hierarchy object built from the ground up.
#'
#' @name ex_hierarchy1
#' @format A [hierarchy()] object with
#' \itemize{
#'  \item name: Poaceae / rank: family / id: 4479
#'  \item name: Poa / rank: genus / id: 4544
#'  \item name: Poa annua / rank: species / id: 93036
#' }
#' Based on NCBI taxonomic classification
#' @source Created from the example code in the [hierarchy()]
#'   documentation.
#' @family taxa-datasets
#' @keywords data
NULL

#' An example Hierarchy object
#'
#' An example Hierarchy object built from the ground up.
#'
#' @name ex_hierarchy2
#' @format A [hierarchy()] object with
#' \itemize{
#'  \item name: Felidae / rank: family / id: 9681
#'  \item name: Puma / rank: genus / id: 146712
#'  \item name: Puma concolor / rank: species / id: 9696
#' }
#' Based on NCBI taxonomic classification
#' @source Created from the example code in the [hierarchy()]
#'   documentation.
#' @family taxa-datasets
#' @keywords data
NULL

#' An example Hierarchy object
#'
#' An example Hierarchy object built from the ground up.
#'
#' @name ex_hierarchy3
#' @format A [hierarchy()] object with
#' \itemize{
#'  \item name: Chordata / rank: phylum / id: 158852
#'  \item name: Vertebrata / rank: subphylum / id: 331030
#'  \item name: Teleostei / rank: class / id: 161105
#'  \item name: Salmonidae / rank: family / id: 161931
#'  \item name: Salmo / rank: genus / id: 161994
#'  \item name: Salmo salar / rank: species / id: 161996
#' }
#' Based on ITIS taxonomic classification
#' @source Created from the example code in the [hierarchy()]
#'   documentation.
#' @family taxa-datasets
#' @keywords data
NULL

#' An example hierarchies object
#'
#' An example hierarchies object built from the ground up.
#'
#' @name ex_hierarchies
#' @format A [hierarchies()] object.
#' @source Created from the example code in the [hierarchies()]
#'   documentation.
#' @family taxa-datasets
#' @keywords data
NULL

#' Lookup-table for IDs of taxonomic ranks
#'
#' Composed of two columns:
#' \itemize{
#'  \item rankid - the ordered identifier value. lower values mean higher rank
#'  \item ranks - all the rank names that belong to the same level, with
#'  different variants that mean essentially the same thing
#' }
#'
#' @name ranks_ref
#' @docType data
#' @keywords data
NULL
