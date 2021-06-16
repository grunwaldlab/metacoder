#' @examples
#' # With Hierarchy class object
#' ex_hierarchy1
#' ## ranks
#' pop(ex_hierarchy1, ranks("family"))
#' ex_hierarchy1 %>% pop(ranks("family"))
#' ex_hierarchy1 %>% pop(ranks("family", "genus"))
#' ## taxon names
#' ex_hierarchy1 %>% pop(nms("Poa"))
#' ex_hierarchy1 %>% pop(nms("Poaceae", "Poa"))
#' ## taxon ids
#' ex_hierarchy1 %>% pop(ids(4479))
#' ex_hierarchy1 %>% pop(ids(4479, 4544))
#' ## mixed: ids and names
#' ex_hierarchy1 %>% pop(ranks("family"), ids(4544))
#'
#' # With hierarchies class object
#' # single taxonomic group
#' invisible(lapply(ex_hierarchies, print))
#' ex_hierarchies %>% pop(ranks("family")) %>% lapply(., print) %>% invisible
#' ## more than one taxonomic group
#' invisible(lapply(ex_hierarchies, print))
#' ex_hierarchies %>% pop(ranks("family", "genus")) %>% lapply(., print) %>%
#'   invisible
