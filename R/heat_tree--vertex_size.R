#===================================================================================================
#' Get all distances between points
#' 
#' Returns the distances between every possible combination of two points.
#' 
#' @param x (\code{numeric} of length 1) x coordinate
#' @param y (\code{numeric} of length 1) y coordinate
#' 
#' @return A \code{data.frame}
#' @keywords internal
molten_dist <- function(x, y) {
  data <- as.matrix(stats::dist(cbind(x, y)))
  data[!lower.tri(data)] <- NA
  molten_data <- data.frame(index_1 = rep(1:nrow(data), ncol(data)), index_2 = rep(1:nrow(data), each = ncol(data)))
  molten_data$distance <- mapply(function(i1, i2) data[i1, i2], molten_data$index_1, molten_data$index_2)
  molten_data[!is.na(molten_data$distance), ]    
}

#===================================================================================================
#' Finds the gap/overlap of circle coordinates
#' 
#' Given a set of x, y coordinates and corresponding radii return the gap between every possible 
#' combination.
#' 
#' @param x (\code{numeric} of length 1) x coordinate of center
#' @param y (\code{numeric} of length 1) y coordinate of center
#' @param r (\code{numeric} of length 1) The diameter of the circle.
#' 
#' @keywords internal
inter_circle_gap <- function(x, y, r) {
  # Force x, y, and r to same length ---------------------------------------------------------------
  temp <- as.data.frame(cbind(x, y, r))
  x <- temp$x
  y <- temp$y
  r <- temp$r
  # Get distance between all points ----------------------------------------------------------------
  data <- molten_dist(x, y)
  # Get distance between circles -------------------------------------------------------------------
  data$gap <- data$distance - r[data$index_1] - r[data$index_2]
  return(data)
}



#===================================================================================================
#' Find optimal range
#' 
#' Finds optimal max and min value using an optimality criterion.
#' 
#' @param max_range (\code{numeric} of length 2) The min and max boundaries to the search space for
#' the optimal maximum value.
#' @param min_range (\code{numeric} of length 2) The min and max boundaries to the search space for
#' the optimal minimum value.
#' @param resolution (\code{numeric} of length 2) The number of increments in each dimension.
#' @param opt_crit (\code{function}) A function that takes two arguments, the max and min, and
#' returns the optimality statistic.
#' @param choose_best (\code{function}) A function that takes a list of \code{opt_crit} outputs
#' and returns the index of the best one.
#' 
#' 
#' @keywords internal
get_optimal_range <- function(max_range, min_range, resolution, opt_crit, choose_best, minimize = TRUE) {
  # Validate arguments -----------------------------------------------------------------------------
  if (length(max_range) != 2 || max_range[1] > max_range[2]) stop('Invalid `max_range`')
  if (length(min_range) != 2 || min_range[1] > min_range[2]) stop('Invalid `min_range`')
  # List of ranges to test -------------------------------------------------------------------------
  max_increments <- seq(from = max_range[1], to = max_range[2], length.out = resolution[2])
  min_increments <- seq(from = min_range[1], to = min_range[2], length.out = resolution[1])
  search_space <- lapply(min_increments, function(x) lapply(max_increments, function(y) c(x, y)))
  search_space <- unlist(search_space, recursive = F)
  valid_range <- vapply(search_space, function(x) x[1] < x[2], logical(1))
  search_space <- search_space[valid_range]
  # Find optimal range -----------------------------------------------------------------------------
  scores <- lapply(search_space, function(x) opt_crit(x[1], x[2]))
  search_space[[choose_best(scores)]]
}

