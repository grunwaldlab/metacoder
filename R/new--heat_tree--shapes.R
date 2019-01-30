#' All built-in node shape functions
#' 
#' Returns a named list of all built-in node shape functions.
#' 
#' @param simple If TRUE, only return shape function for simple geometric shapes
#' 
#' @return named list of functions
#' 
#' @export
node_shape_functions <- function(simple = TRUE) {
  simple_shapes <- list(circle = node_shape_circle,
                        square = node_shape_circle,
                        triangle = node_shape_circle,
                        pentagon = node_shape_circle)
  
  complex_shapes <- list()
  
  if (simple) {
    return(simple_shapes)
  } else {
    return(c(simple_shapes, complex_shapes))
  }
}


#' All built-in node shape functions
#' 
#' Returns a named list of all built-in node shape functions.
#' 
#' @param simple If TRUE, only return shape function for simple geometric shapes
#' 
#' @return named list of functions
#' 
#' @export
edge_shape_functions <- function(simple = TRUE) {
  simple_shapes <- list(default = edge_shape_default)
  
  complex_shapes <- list()
  
  if (simple) {
    return(simple_shapes)
  } else {
    return(c(simple_shapes, complex_shapes))
  }
}


#' Node shape function: circle
#' 
#' A function used to plot nodes in a heat tree. This function plots a simple circle in a solid color. 
node_shape_circle <- function(node_color = "grey", n_verticies = 30) {
  
}


#' Edge shape function: circle
#' 
#' A function used to plot edges in a heat tree. This function plots the default lines in a solid color. 
edge_shape_default <- function(edge_color = "grey") {
  
}