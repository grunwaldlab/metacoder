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
#' @param group_prefix (\code{character} of length 1) The prefix of the group field in the shape data returned
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
#' @param title (\code{character} of length 1) The title of the legend
#' @param color_axis_label (\code{character} of length 1) The label for the color axis
#' @param size_axis_label (\code{character} of length 1) The label for the size axis
#' @param hide_size (\code{logical} of length 1) If \code{TRUE} hide size axis
#' @param hide_color (\code{logical} of length 1) If \code{TRUE} hide color axis
#' @param axis_label_size (\code{numeric} of length 1) 
#' @keywords internal 
make_plot_legend <- function(x, y, length, width_range, width_trans_range = NULL, width_stat_range, group_prefix,
                             tick_size = .008, width_stat_trans = function(x) {x},
                             width_title = "Size", width_sig_fig = 3,
                             color_range, color_trans_range = NULL, color_stat_range, color_stat_trans = function(x) {x},
                             color_title = "Color", color_sig_fig = 3,
                             divisions = 100, label_count = 7, title = NULL, label_size = 0.09,
                             title_size = 0.11, axis_label_size = 0.11,
                             color_axis_label = NULL, size_axis_label = NULL, hide_size = FALSE,
                             hide_color = FALSE) {
  # if the color is defined explicitly, do not print color scale labels
  explicit_color_scale <- any(is.character(color_stat_range)) || is.null(color_stat_range)
  if (explicit_color_scale) {
    hide_color <- TRUE
  }
  
  # If only color is used, set width range to average
  if (hide_size && !hide_color) {
    width_range <- rep(mean(width_range), 2)
  }

    # If size and color are the same, only show color scale
  if (hide_color == FALSE && all(width_stat_range == color_stat_range) && isTRUE(all.equal(width_stat_trans, color_stat_trans))) {
    hide_size = TRUE
  }
  
  # if both scales are hidden then return null
  if (hide_size && hide_color) {
    return(NULL)
  }
  
  
  # if (length(unique(width_range)) == 1) {
  #   size_axis_label = NULL
  # }
  # if (length(unique(color_range)) == 1) {
  #   color_axis_label = NULL
  # }
  
  # Scale bar shapes ===============================================================================
  # Scale bar coordinates --------------------------------------------------------------------------
  prop_div <- seq(0, 1, length.out = divisions)
  point_data <- data.frame(x = max(width_range) - prop_div * (max(width_range) - min(width_range)),
                           y = prop_div * length)
  prop_seg <- vapply(1:(divisions - 1), function(i) mean(prop_div[c(i, i + 1)]), FUN.VALUE = numeric(1))
  if (hide_color) {
    seq_color <- rep("#000000", length(prop_seg))
  } else {
    seq_color <- apply_color_scale(rev(prop_seg), color_range) 
  }
  scale_data <- lapply(1:(divisions - 1), 
                       function(i) scale_bar_coords(x1 = point_data$x[i + 1],
                                                    x2 = point_data$x[i],
                                                    y1 = point_data$y[i + 1] + length * 0.001,
                                                    y2 = point_data$y[i] - length * 0.001,
                                                    color = seq_color[i],
                                                    group = paste0("scale-", group_prefix, i)))
  scale_data <- do.call(rbind, scale_data)
  shape_data <- scale_data
  
  # Tick mark coordinates --------------------------------------------------------------------------
  tick_color <- "#555555"
  tick_div <- seq(0, 1, length.out = label_count)
  label_point_y <- tick_div * length
  if (!hide_size) {
    tick_coords <- function(center_y) {
      max_y <- center_y + tick_size * length / 2
      min_y <- center_y - tick_size * length / 2
      data.frame(group = paste0("tick-", group_prefix, center_y),
                 x = c(0, 0, rep(max(width_range), 2)),
                 y = c(min_y, max_y, max_y, min_y),
                 color = tick_color)
    }
    tick_data <- lapply(label_point_y, tick_coords)
    tick_data <- do.call(rbind, tick_data)
    shape_data <- rbind(tick_data, shape_data)
  }
  
  # Text output ====================================================================================
  format_label <- function(n) {
    as.character(format(signif(n, digits = 3), scientific = FALSE))
  }
  
  scale_undo_trans <- function(points, my_range, my_trans) {
    pre_scaled <- points#rescale(points, to = c(1, 2)) # needed to avoid giving log inverses big inputs
    trans_points <- my_trans(pre_scaled)
    format_label(rescale(trans_points, to = my_range))
  }
  
  label_color = "#000000"
  label_size = length * label_size 
  label_data <- NULL
  
  # Size axis labels -------------------------------------------------------------------------------
  if (!hide_size) {
    label_data <- rbind(label_data, 
                        data.frame(stringsAsFactors = FALSE, 
                                   label = scale_undo_trans(label_point_y, width_stat_range, width_stat_trans),
                                   x = max(width_range) * 1.1,
                                   y = rev(label_point_y),
                                   size = label_size,
                                   color = label_color,
                                   rotation = 0,
                                   justification = "left",
                                   group = "legend"))
    
  }
  # Color axis labels ------------------------------------------------------------------------------
  if (!hide_color) {
    label_data <- rbind(label_data, 
                        data.frame(stringsAsFactors = FALSE, 
                                   label = scale_undo_trans(rescale(label_point_y, to = color_trans_range), color_stat_range, color_stat_trans),
                                   x = 0 - max(width_range) * 0.1,
                                   y = rev(label_point_y),
                                   size = label_size,
                                   color = label_color,
                                   rotation = 0,
                                   justification = "right",
                                   group = "legend"))
  }
  
  # Add color axis title ---------------------------------------------------------------------------
  y_range <- range(point_data$y)
  axis_label_size <- length * axis_label_size
  axis_spacer <- axis_label_size * 0.5
  if (!hide_color && !is.null(color_axis_label)) {
    grob_lengths <- text_grob_length(label_data$label) * label_data$size
    label_x_min <- min(label_data$x - ifelse(label_data$justification == "right", grob_lengths, 0))
    label_data <- rbind(label_data, 
                        data.frame(stringsAsFactors = FALSE, 
                                   label = color_axis_label,
                                   x = label_x_min - axis_spacer,
                                   y = mean(y_range),
                                   size = axis_label_size,
                                   color = "#000000",
                                   rotation = pi / 2,
                                   justification = "center-bottom",
                                   group = "legend"))
  }
  
  # Add size axis title ----------------------------------------------------------------------------
  if (!hide_size && !is.null(size_axis_label)) {
    grob_lengths <- text_grob_length(label_data$label) * label_data$size
    label_x_max <- max(label_data$x + ifelse(label_data$justification == "left", grob_lengths, 0))
    label_data <- rbind(label_data, 
                        data.frame(stringsAsFactors = FALSE, 
                                   label = size_axis_label,
                                   x = label_x_max + axis_spacer,
                                   y = mean(y_range),
                                   size = axis_label_size,
                                   color = "#000000",
                                   rotation = pi / 2,
                                   justification = "center-top",
                                   group = "legend"))
  }
  
  # Add legend title -------------------------------------------------------------------------------
  if (!is.null(title)) {
    title_size <- length * title_size
    title_spacer <- title_size * 0.5
    label_y_max <- max(label_data$y + label_size / 2)
    label_data <- rbind(label_data, 
                        data.frame(stringsAsFactors = FALSE, 
                                   label = title,
                                   x = min(scale_data$x),
                                   y = label_y_max + title_spacer,
                                   size = title_size,
                                   color = "#000000",
                                   rotation = 0,
                                   justification = "left-bottom",
                                   group = "legend"))
  } 
  
  
  # Set origin to c(0, 0)
  min_x <- min(c(shape_data$x, label_data$x))
  if (min_x < 0) {
    x <- x + abs(min_x)
  }

  # Transpose ouput
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



#' Generate the inverse of a function
#' 
#' http://stackoverflow.com/questions/10081479/solving-for-the-inverse-of-a-function-in-r
#' 
#' @param f (\code{function} with one argument) A function to derive and inverse from
#' @param interval (\code{character} of length 2) The range of the value the inverse function can return.
#' 
#' @return (\code{function}) Return the inverse of the function given
#' 
#' @keywords internal
inverse = function (f, interval) {
  interval <- c(interval[1] - .001, interval[2] + .001) # ensures that min < max
  function (y) {
    process_one <- function(one_y) {
      stats::uniroot((function (x) f(x) - one_y), interval = interval, extendInt = "yes")[1]
    }
    unname(unlist(lapply(y, process_one)))
  }
}
