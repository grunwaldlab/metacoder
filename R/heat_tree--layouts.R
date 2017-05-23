#' Layout functions
#' 
#' Functions used to determine graph layout.
#' Calling the function with no parameters returns available function names.
#' Calling the function with only the name of a function returns that function.
#' Supplying a name and a \code{\link[igraph]{graph}} object to run the layout function on the graph.
#' 
#' @param name (\code{character} of length 1 OR NULL) name of algorithm. Leave \code{NULL} to 
#' see all options. 
#' @param graph (\code{igraph}) The graph to generate the layout for.
#' @param intitial_coords (\code{matrix}) Initial node layout to base new layout off of.
#' @param effort  (\code{numeric} of length 1) The amount of effort to put into layouts. Typically
#' determines the the number of iterations. 
#' @param ... (other arguments) Passed to igraph layout function used.
#' 
#' @return The name available functions, a layout functions,
#' or a two-column matrix depending on how arguments are provided.
#' 
#' @examples 
#' # List available function names:
#' layout_functions()
#' 
#' # Execute layout function on graph:
#' layout_functions("davidson-harel", igraph::make_ring(5))
#' 
#' @export
layout_functions <- function(name = NULL, graph = NULL, intitial_coords = NULL, effort = 1, ...) {
  funcs <- list("automatic" = igraph::nicely,
                "reingold-tilford" =  igraph::as_tree,
                "davidson-harel" = igraph::with_dh,
                "gem" = igraph::with_gem,
                "graphopt" = igraph::with_graphopt,
                "mds" = igraph::with_mds(),
                "fruchterman-reingold" = igraph::with_fr,
                "kamada-kawai" = igraph::with_kk,
                "large-graph" = igraph::with_lgl,
                "drl" = igraph::with_drl)
  return_names <- is.null(name) && is.null(graph) && is.null(intitial_coords)
  if (return_names) {
    return(names(funcs))
  } else {
    v_weight <- igraph::V(graph)$weight_factor
    e_weight <- igraph::E(graph)$weight_factor
    e_density <- igraph::edge_density(graph)
    defaults <- list("automatic" = list(),
                     "reingold-tilford" = list(circular = TRUE,
                                               mode = "out"),
                     "davidson-harel" = list(coords = intitial_coords,
                                             maxiter = 10 * effort,
                                             fineiter = max(10, log2(igraph::vcount(graph))) * effort,
                                             cool.fact = 0.75 - effort * 0.1,
                                             weight.node.dist = 13, #* ifelse(is.null(v_weight), 1, list(rescale(v_weight, c(.1, 10))))[[1]], #higher values spread out nodes 
                                             weight.border = 0,
                                             weight.edge.lengths = 0.5, #* ifelse(is.null(e_weight), 1, list(rescale(e_weight, c(10, 1))))[[1]], # higher number spread the graph out more
                                             weight.edge.crossings = 100,
                                             weight.node.edge.dist = 1), #* ifelse(is.null(v_weight), 1, list(rescale(v_weight, c(.1, 10))))[[1]]), 
                     "gem" = list(coords = intitial_coords,
                                  maxiter = 40 * igraph::vcount(graph)^2 * effort,
                                  temp.max = igraph::vcount(graph) * (1 + effort * 0.1),
                                  temp.min = 1/10,
                                  temp.init = sqrt(igraph::vcount(graph))),
                     "graphopt" = list(start = intitial_coords,
                                       niter = 500 * effort,
                                       charge = 0.0005,
                                       mass = 30,
                                       spring.length = 0,
                                       spring.constant = 1,
                                       max.sa.movement = 5),
                     "mds" = list(),
                     "fruchterman-reingold" = list(coords = intitial_coords,
                                                   niter = 500 * effort,
                                                   start.temp = sqrt(igraph::vcount(graph)) * (1 + effort * 0.1) ,
                                                   grid = "nogrid",
                                                   weights = e_weight), #ifelse(is.null(e_weight), 1, list(rescale(e_weight, c(1, 10))))[[1]]^4), #edge weights
                     "kamada-kawai" = list(coords = intitial_coords,
                                           maxiter = 100 * igraph::vcount(graph),
                                           epsilon = 0,
                                           kkconst = igraph::vcount(graph),
                                           weights = NULL),
                     "large-graph" = list(maxiter = 200,
                                          maxdelta = igraph::vcount(graph),
                                          area = igraph::vcount(graph)^2,
                                          coolexp = 1.5 - effort * 0.1,
                                          repulserad = igraph::vcount(graph)^2 * igraph::vcount(graph),
                                          cellsize = sqrt(sqrt(igraph::vcount(graph)^2)),
                                          root = 1),
                     "drl" = list(use.seed = ! is.null(intitial_coords),
                                  seed = ifelse(is.null(intitial_coords), 
                                                matrix(stats::runif(igraph::vcount(graph) * 2), ncol = 2),
                                                intitial_coords),
                                  options = igraph::drl_defaults$default,
                                  weights = NULL,
                                  fixed = NULL))
    arguments <- utils::modifyList(defaults[[name]], list(...))
    coords <- igraph::layout_(graph = graph, layout = do.call(funcs[[name]], arguments))
    return(coords)
  }
}
