

#' Make color/size legend
#' 
#' Make color/size legend
#' 
#' @param length (\code{numeric} of length 1) the length of the scale bar
#' @param width_range (\code{numeric} of length 1 or 2) the width of the scale bar or the range
#' @param width_stat_range (\code{numeric} of length 1 or 2) The stat range to display in the size labels
#' @param width_stat_trans (\code{function}) The transformation used to convert the statistic to size
#' @param width_title (\code{character} of length 1) The title of the size labels.
#' @param width_sig_fig (\code{numeric} of length 1) The number of significant figures to use in size labels.
#' @param color_range (\code{character}) One ore more hex codes constituting a color scale.
#' @param color_stat_range (\code{numeric} of length 1 or 2) The stat range to display in the color labels
#' @param color_stat_trans (\code{function}) The transformation used to convert the statistic to size
#' @param color_title (\code{character} of length 1) The title of the color labels.
#' @param color_sig_fig (\code{numeric} of length 1) The number of significant figures to use in color labels.
#' @param divisions (\code{numeric} of length 1) The number of colors to display. 
#' @param label_count (\code{numeric} of length 1) The number of labels.
make_plot_legend <- function(length, width_range, width_stat_range, width_stat_trans = function(x) {x},
                             width_title = "Size", width_sig_fig = 3,
                             color_range, color_stat_range, color_stat_trans = function(x) {x},
                             color_title = "Color", color_sig_fig = 3,
                             divisions = 100, label_count = 5) {
  inverse = function (f, lower = -100, upper = 100) {
    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
  }
  
  
  # Generate scale bar coordinates
  prop_div <- seq(0, 1, length.out = divisions)
  point_data <- data.frame(x = max(width_range) - prop_div * (max(width_range) - min(width_range)),
                           y = prop_div * length)
  prop_seg <- vapply(1:(divisions - 1), function(i) mean(prop_div[c(i, i + 1)]), FUN.VALUE = numeric(1))
  seq_color <- apply_color_scale(prop_seg, color_range) 
  scale_data <- lapply(1:(divisions - 1), 
                       function(i) scale_bar_coords(x1 = point_data$x[i + 1],
                                                    x2 = point_data$x[i],
                                                    y1 = point_data$y[i + 1],
                                                    y2 = point_data$y[i],
                                                    color = seq_color[i],
                                                    group = i))
  scale_data <- do.call(rbind, scale_data)
  # Generate tick mark coordinates
  
  # Generate label coordinates
  
  # Generate title coordinates
  
  
  # Graph 
  ggplot2::ggplot(data = scale_data) +
    ggplot2::geom_polygon(data = scale_data, ggplot2::aes(x = x, y = y, group = group),
                          fill = scale_data$color) 
}


#' Make scale bar division
#' 
#' Make scale bar division
#' 
#' @param x1 (\code{numeric} of length 1) x of top right 
#' @param x2 (\code{numeric} of length 1) x of bottom right
#' @param y1 (\code{numeric} of length 1) y of top right 
#' @param y2 (\code{numeric} of length 1) y of bottom right
#' @param color 
#' @param group 
#' 
#' @return \code{data.frame}
scale_bar_coords <- function(x1, x2, y1, y2, color, group) {
  data.frame(x = c(x1, x2, 0, 0),
             y = c(y1, y2, y2, y1),
             color = color, 
             group = group)
}

