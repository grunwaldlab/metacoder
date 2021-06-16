#' @examples
#' # Hierarchy class
#' ex_hierarchy1
#'
#' ## ranks
#' ### keep all taxa between family and genus
#' span(ex_hierarchy1, ranks("family", "genus"))
#' span(ex_hierarchy1, nms("Poaceae", "Poa"))
#' span(ex_hierarchy1, ids(4479, 4544))
#'
#' ### keep all taxa between genus and species
#' span(ex_hierarchy1, ranks("genus", "species"))
#'
#' ### keep all taxa greater than genus
#' span(ex_hierarchy1, ranks("> genus"))
#'
#' ### keep all taxa greater than or equal to genus
#' span(ex_hierarchy1, ranks(">= genus"))
#'
#' ### keep all taxa less than genus
#' span(ex_hierarchy1, ranks("< genus"))
#'
#' ### keep all taxa less than or equal to genus
#' span(ex_hierarchy1, ranks("<= genus"))
#'
#' ### same as above, with different dataset
#' span(ex_hierarchy2, ranks("> genus"))
#' span(ex_hierarchy2, ranks(">= genus"))
#' span(ex_hierarchy2, ranks("< genus"))
#' span(ex_hierarchy2, ranks("<= genus"))
#'
#' # using taxonomic names
#' span(ex_hierarchy2, nms("< Felidae"))
#'
#' # using taxonomic ids
#' span(ex_hierarchy2, ids("< 9681"))
#'
#' ## Multiple operator statements - useful with larger classifications
#' ex_hierarchy3
#' span(ex_hierarchy3, ranks("> genus"), ranks("< phylum"))
#' span(ex_hierarchy3, ids("> 161994"), ids("< 158852"))
#'
#'
#' ## taxon names
#' ### keep all taxa between Poaceae and Poa
#' ### - matches to ranks first
#' ex_hierarchy1 %>% span(nms("Poaceae", "Poa"))
#'
#' ## taxon ids
#' ### keep all taxa between 4479 and 4544 taxonomic IDs
#' ### - matches to ranks first
#' ex_hierarchy1 %>% span(ids(4479, 4544))
#'
#'
#' # hierarchies class
#' invisible(lapply(ex_hierarchies, print))
#' ex_hierarchies %>% span(ranks("family", "genus")) %>% lapply(., print) %>%
#' invisible
