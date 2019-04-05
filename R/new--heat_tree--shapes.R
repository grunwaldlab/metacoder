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
                        square = node_shape_square,
                        triangle = node_shape_triangle,
                        pentagon = node_shape_pentagon,
                        hexagon = node_shape_hexagon)
  
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
#' 
#' Inspired by (i.e. stolen from) https://gist.github.com/baptiste/2224724, which was
#' itself inspired from a post by William Dunlap on r-help (10/09/09)
#' 
#' @param node_color The color of the circle fill
#' @param n_verticies The number of verticies used to simulate a circle using a regular polygon.
#' @param angle (\code{numeric} of length 1) Angle to rotate points around the center of the circle.
#' 
#' @export
node_shape_circle <- function(node_color = "grey", n_verticies = 30, angle = 0) {
  if(n_verticies <= 2) stop("n must be more than 2!")
  
  # Calculate imaginary coords 
  coords <- exp(seq(0, n_verticies)*2i*pi/n_verticies)
  
  # Rotate around center
  coords <- coords * exp(1i * angle)
  
  # Translate to x, y coords 
  coords <- tibble::tibble(x = Re(coords), y = Im(coords), color = node_color)
  
  return(coords)
}


#' Node shape function: square
#' 
#' A function used to plot nodes in a heat tree. This function plots a simple square in a solid color. 
#' 
#' @param node_color The color of the square fill
#' 
#' @export
node_shape_square <- function(node_color = "grey") {
  tibble::tibble(x = c(-1, -1, 1, 1),
                 y = c(-1, 1, 1, -1),
                 color = node_color)
}


#' Node shape function: triangle
#' 
#' A function used to plot nodes in a heat tree. This function plots a simple triangle in a solid color. 
#' 
#' @param node_color The color of the triagle fill
#' 
#' @export
node_shape_triangle <- function(node_color = "grey") {
  h = sqrt(3)
  tibble::tibble(x = c(-1, 0, 1),
                 y = c(-h/3, h * 2 / 3, -h/3),
                 color = node_color)
}


#' Node shape function: pentagon
#' 
#' A function used to plot nodes in a heat tree. This function plots a simple pentagon in a solid color. 
#' 
#' @param node_color The color of the pentagon fill
#' 
#' @export
node_shape_pentagon <- function(node_color = "grey") {
  node_shape_circle(node_color = node_color, n_verticies = 5, angle = 0.314)
}


#' Node shape function: hexagon
#' 
#' A function used to plot nodes in a heat tree. This function plots a simple hexagon in a solid color. 
#' 
#' @param node_color The color of the hexagon fill
#' 
#' @export
node_shape_hexagon <- function(node_color = "grey") {
  node_shape_circle(node_color = node_color, n_verticies = 6, angle = 0.314)
}


#' Edge shape function: circle
#' 
#' A function used to plot edges in a heat tree. This function plots the default lines in a solid color. 
edge_shape_default <- function(edge_color = "grey") {
  
}