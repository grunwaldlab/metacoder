#===================================================================================================
#' Makes coordinates for a regular polygon
#' 
#' Generates an n x 2 matrix containing x and y coordinates between 1 and 0 for the points of a 
#' regular polygon. 
#' 
#' Inspired by (i.e. stolen from) https://gist.github.com/baptiste/2224724, which was
#' itself inspired from a post by William Dunlap on r-help (10/09/09)
#' 
#' @param n (\code{numeric} of length 1) The number of nodes in the polygon.
#' @param x (\code{numeric} of length 1) x coordinate of center
#' @param y (\code{numeric} of length 1) y coordinate of center
#' @param radius (\code{numeric} of length 1) The diameter of the circle.
#' @param angle (\code{numeric} of length 1) Angle to rotate points around the center of the circle.
#' 
#' @keywords internal
polygon_coords <- function(n = 5, x = 0, y = 0, radius = 1, angle = 0){
  # Define function to make points for a single polygon --------------------------------------------
  process_one <- function(n, x, y, r, a) {
    if(n<3) stop("n must be more than 3!")
    # Calculate imaginary coords - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    coords <- exp(seq(0, n)*2i*pi/n)
    # Rotate around center - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    coords <- coords * exp(1i*a)
    # Translate to x, y coords - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    coords <- data.frame(x = Re(coords), y = Im(coords))
    # Scale to radius  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    coords <- coords * r
    # Offset center to given x, y coords - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    coords$x <- coords$x + x
    coords$y <- coords$y + y
    return(coords)
  }
  # Compile of the points for multiple polygons ----------------------------------------------------
  output <- mapply(process_one, n, x, y, radius, angle, SIMPLIFY = FALSE)
  group <- rep(seq_along(output), vapply(output, nrow, numeric(1)))
  cbind(group = as.factor(group), do.call(rbind, output))
} 



#===================================================================================================
#' Makes coordinates for a line
#' 
#' Generates an n x 2 matrix containing x and y coordinates between 1 and 0 for the points of a 
#' line with a specified width in cartesian coordinates. 
#' 
#' @param x1 (\code{numeric} of length 1) x coordinate of the center of one end
#' @param y1 (\code{numeric} of length 1) y coordinate of the center of one end
#' @param x2 (\code{numeric} of length 1) x coordinate of the center of the other end
#' @param y2 (\code{numeric} of length 1) y coordinate of the center of the other end
#' @param width (\code{numeric} of length 1) The width of the line.
#' 
#' @keywords internal
line_coords <- function(x1, y1, x2, y2, width) {
  # Define function to make points for a single line rect ------------------------------------------
  process_one <- function(x1, y1, x2, y2, w) {
    slope <- (y2 - y1) / (x2 - x1)
    inv_slope <- -1/slope
    angle <- atan(inv_slope)
    off_x <- w / 2 * cos(angle)
    off_y <- w / 2 * sin(angle)
    data.frame(x = c(x1 + off_x, x1 - off_x, x2 - off_x, x2 + off_x),
               y = c(y1 + off_y, y1 - off_y, y2 - off_y, y2 + off_y))
  }
  # Compile of the points for multiple polygons ----------------------------------------------------
  output <- mapply(process_one, x1, y1, x2, y2, width, SIMPLIFY = FALSE)
  group <- rep(seq_along(output), vapply(output, nrow, numeric(1)))
  cbind(group = as.factor(group), do.call(rbind, output))
}
