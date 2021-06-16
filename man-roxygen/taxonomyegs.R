#' @examples
#'
#' # Making a taxonomy object with vectors
#' taxonomy(c("mammalia", "felidae", "panthera", "tigris"),
#'          c("mammalia", "felidae", "panthera", "leo"),
#'          c("mammalia", "felidae", "felis", "catus"))
#'
#' # Making a taxonomy object from scratch
#' #   Note: This information would usually come from a parsing function.
#' #         This is just for demonstration.
#' x <- taxon(
#'   name = taxon_name("Notoryctidae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(4479)
#' )
#' y <- taxon(
#'   name = taxon_name("Notoryctes"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(4544)
#' )
#' z <- taxon(
#'   name = taxon_name("Notoryctes typhlops"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(93036)
#' )
#'
#' a <- taxon(
#'   name = taxon_name("Mammalia"),
#'   rank = taxon_rank("class"),
#'   id = taxon_id(9681)
#' )
#' b <- taxon(
#'   name = taxon_name("Felidae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(9681)
#' )
#'
#' cc <- taxon(
#'   name = taxon_name("Puma"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(146712)
#' )
#' d <- taxon(
#'   name = taxon_name("Puma concolor"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(9696)
#' )
#'
#' m <- taxon(
#'   name = taxon_name("Panthera"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(146712)
#' )
#' n <- taxon(
#'   name = taxon_name("Panthera tigris"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(9696)
#' )
#'
#' (hier1 <- hierarchy(z, y, x, a))
#' (hier2 <- hierarchy(cc, b, a, d))
#' (hier3 <- hierarchy(n, m, b, a))
#'
#' (hrs <- hierarchies(hier1, hier2, hier3))
#'
#' taxonomy(hier1, hier2, hier3)
