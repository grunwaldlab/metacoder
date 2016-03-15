f <- function(x) {
  sfsddsdeognmvdfs = 2
  eval(substitute(x), envir = environment())
}



#' Create a restrictive subset of \code{\link{classified}}
#' 
#' Create a subset of items classified by a taxonomy.
#' Only specified taxa and shared parent taxa of specified taxa with specified items are preserved.
#' Only specified items assigned to specified taxa are preserved.
#' 
#' @param taxon A key to filter the taxon data.frame rows on
#' @param item A key to filter the item data.frame rows on
#' @param column The name of a column in either the item or the taxon data.frame
#' 
#' @return \code{\link{classified}}
`[[.classified` <- function(taxon, item, column) {
  
}


#' Create a inclusive subset of \code{\link{classified}}
#' 
#' Create a subset of items classified by a taxonomy.
#' Only unspecified taxa with no items or children with items are discarded.
#' 
#' @param taxon A key to filter the taxon data.frame rows on
#' @param item A key to filter the item data.frame rows on
#' @param column The name of a column in either the item or the taxon data.frame
#' 
#' @return \code{\link{classified}}
`[.classified` <- function(taxon, item, column) {
  
}