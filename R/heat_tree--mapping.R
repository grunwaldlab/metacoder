#' Rescale numeric vector to have specified minimum and maximum.
#' 
#' Rescale numeric vector to have specified minimum and maximum, but allow for hard boundaries.
#' It is a slightly modified version of scales::rescale, incorporating scales::zero_range, both by Hadley Wickham used under the conditions of the MIT license.
#' 
#' @param x values to rescale
#' @param to range to scale to
#' @param from range of values the x could have been
#' @param hard_bounds If \code{TRUE}, all values will be forced into the range of \code{to}.
#' 
#' @keywords internal
rescale <- function (x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE), hard_bounds = TRUE) 
{
  # COPIED FROM scales::zero_range by Hadley Wickham
  zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
    if (length(x) == 1) return(TRUE)
    if (length(x) != 2) stop("x must be length 1 or 2")
    if (any(is.na(x))) return(NA)
    
    # Special case: if they are equal as determined by ==, then there
    # is zero range. Also handles (Inf, Inf) and (-Inf, -Inf)
    if (x[1] == x[2]) return(TRUE)
    
    # If we reach this, then x must be (-Inf, Inf) or (Inf, -Inf)
    if (all(is.infinite(x))) return(FALSE)
    
    # Take the smaller (in magnitude) value of x, and use it as the scaling
    # factor.
    m <- min(abs(x))
    
    # If we get here, then exactly one of the x's is 0. Return FALSE
    if (m == 0) return(FALSE)
    
    # If x[1] - x[2] (scaled to 1) is smaller than tol, then return
    # TRUE; otherwise return FALSE
    abs((x[1] - x[2]) / m) < tol
  }
  
  # COPIED FROM scales::rescale by Hadley Wickham
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  }
  result <- (x - from[1]) / diff(from) * diff(to) + to[1]
  
  # Hard bound implementations 
  if (hard_bounds) {
    result[result > max(to)] <- max(to)
    result[result < min(to)] <- min(to)
  }
  return(result)
}



#' Covert numbers to colors
#' 
#' Convert numbers to colors.
#' If colors are already supplied, return the input
#' 
#' @param values (\code{numeric}) The numbers to represent as colors
#' @param color_series (\code{character}) Hex values or a character in \code{colors}
#' @param no_color_in_palette (\code{numeric} of length 1) The number of distinct colors to use.
#' @param interval (\code{numeric} of length 2) The range \code{values} could have taken.
#' 
#' 
#' @return \code{character} Hex color codes. 
#' 
#' @keywords internal
apply_color_scale <- function(values, color_series, interval = NULL, no_color_in_palette = 1000) {
  numeric_values <- get_numerics(values)
  if (length(numeric_values) > 0) { ## Not factors, characters, or hex codes
    palette <- grDevices::colorRampPalette(color_series)(no_color_in_palette)
    if (is.null(interval)) {
      interval <- range(numeric_values, na.rm = TRUE, finite = TRUE)
    }
    color_index <- as.integer(rescale(numeric_values, to = c(1, no_color_in_palette), from = interval))
    output <- values
    output[can_be_num(values)] <- palette[color_index]
    return(output)
  } else {
    return(values)
  }
}



#' The default quantative color palette
#' 
#' Returns the default color palette for quantative data.
#' 
#' @return \code{character} of hex color codes
#' 
#' @examples
#' quantative_palette()
#' 
#' @export
quantative_palette <- function() {
  # produced with: c("#BBBBBB", rev(viridisLite::viridis(7, begin = .4, end = .9)))
  return(c("#BBBBBB", "#BBDF27FF", "#85D44AFF", "#54C568FF", "#2FB47CFF", 
           "#1FA188FF", "#228C8DFF", "#2A788EFF"))
}


#' The default qualitative color palette
#' 
#' Returns the default color palette for qualitative data
#' 
#' @return \code{character} of hex color codes
#' 
#' @examples
#' qualitative_palette()
#' 
#' @export
qualitative_palette <- function() {
  # produced with c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(9, "Pastel1"))
  return(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
           "#A65628", "#F781BF", "#999999", "#FBB4AE", "#B3CDE3", "#CCEBC5", 
           "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2"))
}

#' The default diverging color palette
#' 
#' Returns the default color palette for diverging data
#' 
#' @return \code{character} of hex color codes
#' 
#' @examples
#' diverging_palette()
#' 
#' @export
diverging_palette <- function() {
  return(c("#a6611a", "#DDDDDD", "#018571"))
}
