#' Make color/size legend
#' 
#' Make color/size legend
#' 
#' @param x bottom left
#' @param y bottom left
#' @param length (\code{numeric} of length 1) the length of the scale bar
#' @param tick_size (\code{numeric} of length 1) the thickness of tick marks
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
#' @keywords internal
make_plot_legend <- function(x, y, length, tick_size, width_range, width_stat_range, width_stat_trans = function(x) {x},
                             width_title = "Size", width_sig_fig = 3,
                             color_range, color_stat_range, color_stat_trans = function(x) {x},
                             color_title = "Color", color_sig_fig = 3,
                             divisions = 100, label_count = 5) {
  
  
  # Generate scale bar coordinates
  prop_div <- seq(0, 1, length.out = divisions)
  point_data <- data.frame(x = max(width_range) - prop_div * (max(width_range) - min(width_range)),
                           y = prop_div * length)
  prop_seg <- vapply(1:(divisions - 1), function(i) mean(prop_div[c(i, i + 1)]), FUN.VALUE = numeric(1))
  seq_color <- apply_color_scale(rev(prop_seg), color_range) 
  scale_data <- lapply(1:(divisions - 1), 
                       function(i) scale_bar_coords(x1 = point_data$x[i + 1],
                                                    x2 = point_data$x[i],
                                                    y1 = point_data$y[i + 1],
                                                    y2 = point_data$y[i],
                                                    color = seq_color[i],
                                                    group = paste0("scale-",i)))
  scale_data <- do.call(rbind, scale_data)
  
  # Generate tick mark coordinates
  tick_color <- "#555555"
  tick_div <- seq(0, 1, length.out = label_count)
  label_point_y <- tick_div * length
  tick_coords <- function(center_y) {
    max_y <- center_y + tick_size / 2
    min_y <- center_y - tick_size / 2
    data.frame(group = paste0("tick-", center_y),
               x = c(0, 0, rep(max(width_range), 2)),
               y = c(min_y, max_y, max_y, min_y),
               color = tick_color)
  }
  
  
  tick_data <- lapply(label_point_y, tick_coords)
  tick_data <- do.call(rbind, tick_data)
  
  # Generate label coordinates
  format_label <- function(n) {
    # format(n, scientific = FALSE, drop0trailing = TRUE, digits = 3)
    as.character(signif(n, digits = 3))
  }
  
  scale_undo_trans <- function(points, my_range, my_trans) {
    trans_points <- vapply(points, my_trans, FUN.VALUE = numeric(1))
    format_label(scales::rescale(trans_points, to = range(my_range)))
  }
  
  label_color = "#000000"
  label_size = max(width_range) / 5
  
  
  make_size_label_data <- function() {
    data.frame(stringsAsFactors = FALSE, 
               label = scale_undo_trans(label_point_y, width_stat_range, width_stat_trans),
               x = max(width_range) * 1.1,
               y = rev(label_point_y),
               size = label_size,
               color = label_color,
               rotation = 0,
               justification = "left")
  }
  
  make_color_label_data <- function() {
    data.frame(stringsAsFactors = FALSE, 
               label = scale_undo_trans(label_point_y, color_stat_range, color_stat_trans),
               x = 0 - max(width_range) * 0.1,
               y = rev(label_point_y),
               size = label_size,
               color = label_color,
               rotation = 0,
               justification = "right")
  } 
  
  
  if (!is.null(width_stat_range) && !is.null(color_stat_range)) {
    label_data <- rbind(make_size_label_data(), make_color_label_data())
  } else if (!is.null(width_stat_range)) {
    label_data <- make_size_label_data()
  } else if (!is.null(color_stat_range)) {
    label_data <- make_color_label_data()
  } else {
    label_data <- NULL
  }
  
  # Generate title coordinates
  
  
  shape_data <- rbind(tick_data, scale_data)
  shape_data$x <- shape_data$x + x
  shape_data$y <- shape_data$y + y
  
  if (!is.null(label_data)) {
    label_data$x <- label_data$x + x
    label_data$y <- label_data$y + y
    
  }
  
  output <- list(shapes = shape_data,
                 labels = label_data)
  return(output)
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
#' @keywords internal
scale_bar_coords <- function(x1, x2, y1, y2, color, group) {
  data.frame(group = group,
             x = c(x1, x2, 0, 0),
             y = c(y1, y2, y2, y1),
             color = color)
  
}


