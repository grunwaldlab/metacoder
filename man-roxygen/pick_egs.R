#' @examples
#' # ranks
#' ex_hierarchy1
#' ex_hierarchy1 %>% pick(ranks("family"))
#' ex_hierarchy1 %>% pick(ranks("family", "genus"))
#' # taxon names
#' ex_hierarchy1 %>% pick(nms('Poa'))
#' ex_hierarchy1 %>% pick(nms("Poaceae", "Poa"))
#' # taxon ids
#' ex_hierarchy1 %>% pick(ids(4479))
#' ex_hierarchy1 %>% pick(ids(4479, 4544))
#' # mixed: ids and names
#' ex_hierarchy1 %>% pick(ranks("family"), ids(4544))
#'
#' ## single taxonomic group
#' ex_hierarchy1 %>% pick(ranks("family"))
#' pick(ex_hierarchy1, ranks("family"))
#' ### more than 1 - remake res object above first
#' ex_hierarchy1 %>% pick(ranks("family", "genus"))
#'
#'
#' # hierarchies
#' # single taxonomic group
#' invisible(lapply(ex_hierarchies, print))
#' ex_hierarchies %>% pick(ranks("family")) %>% lapply(., print) %>% invisible
#'
#' ## more than one taxonomic group
#' invisible(lapply(ex_hierarchies, print))
#' ex_hierarchies %>% pick(ranks("family", "genus")) %>% lapply(., print) %>%
#'   invisible
