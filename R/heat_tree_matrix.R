heat_tree_matrix <- function(obj, seed = 1, ...) {
  # Make plot layout
  treatments <- unique(c(obj$data$diff_table$treatment_1, obj$data$diff_table$treatment_2))
  combinations <- t(combn(seq_along(treatments), 2))
  layout_matrix <- matrix(rep(NA, (length(treatments))^2), nrow = length(treatments))
  for (index in 1:nrow(combinations)) {
    layout_matrix[combinations[index, 1], combinations[index, 2]] <- index
  }
  
  # Make individual plots
  sub_plots <- lapply(seq_len(nrow(combinations)),
                      function(index) {
                        set.seed(seed)
                        obj %>%
                          taxa::filter_obs("diff_table",
                                           treatment_1 == treatments[combinations[index, 1]] & treatment_2 == treatments[combinations[index, 2]]) %>%
                          metacoder::heat_tree(...,
                                               make_legend = FALSE)
                      })
  
  # Make key plot
  set.seed(seed)
  key_plot <- metacoder::heat_tree(obj, ..., make_legend = TRUE)
  
  calc_subplot_coords <- function(a_matrix, x1 = 0, y1 = 0, x2 = 1, y2 = 1) {
    # lowerleft = c(x1, y1), upperright = c(x2, y2)
    x_coords <- seq(from = x1, to = x2, length.out = ncol(a_matrix) + 1)[- (ncol(a_matrix) + 1)]
    y_coords <- seq(from = y1, to = y2, length.out = nrow(a_matrix) + 1)[- (nrow(a_matrix) + 1)]
    do.call(rbind, lapply(1:ncol(a_matrix), function(x) data.frame(plot_index = a_matrix[, x],
                                                                   x = x_coords[x], 
                                                                   y = rev(y_coords))))
  }
  
  # remove empty column/row
  layout_matrix <- layout_matrix[! apply(layout_matrix, MARGIN = 1, function(x) all(is.na(x))), ]
  layout_matrix <- layout_matrix[ , ! apply(layout_matrix, MARGIN = 2, function(x) all(is.na(x)))]
  
  # Get subplot layout data
  matrix_data <- calc_subplot_coords(layout_matrix, x1 = 0.2, y1 = 0.2, x2 = 0.95, y2 = 0.95)
  matrix_data$treatment_1 <- treatments[combinations[matrix_data$plot_index, 1]]
  matrix_data$treatment_2 <- treatments[combinations[matrix_data$plot_index, 2]]
  matrix_data <- matrix_data[!is.na(matrix_data$plot_index), ]
  matrix_data <- matrix_data[order(matrix_data$plot_index), ]
  rownames(matrix_data) <- matrix_data$plot_index
  
  # Make label data
  named_row <- which(apply(layout_matrix, MARGIN = 1, function(x) all(!is.na(x))))
  named_col <- which(apply(layout_matrix, MARGIN = 2, function(x) all(!is.na(x))))
  horz_label_data <- matrix_data[match(layout_matrix[named_row, ], matrix_data$plot_index), ]
  vert_label_data <- matrix_data[match(layout_matrix[, named_col], matrix_data$plot_index), ]
  subgraph_width <- abs(horz_label_data$x[1] - horz_label_data$x[2])
  subgraph_height <- abs(vert_label_data$y[1] - vert_label_data$y[2])
  horz_label_data$label_x <- horz_label_data$x + subgraph_width / 2 # center of label
  horz_label_data$label_y <- 0.96 # bottom of label
  vert_label_data$label_x <- 0.96 # bottom of rotated label 
  vert_label_data$label_y <- vert_label_data$y + subgraph_height / 2 # center of rotated label 
  
  # Make plot
  label_size <- 12
  matrix_plot <- cowplot::ggdraw() + 
    cowplot::draw_plot(key_plot, x = 0, y = 0, width = 0.75, height = 0.6) +
    cowplot::draw_text(gsub("_", " ", horz_label_data$treatment_2), 
                       x = horz_label_data$label_x, y = horz_label_data$label_y, 
                       size = label_size, colour = diverging_palette()[1],
                       hjust = "center", vjust = "bottom") +
    cowplot::draw_text(gsub("_", " ", vert_label_data$treatment_1), 
                       x = vert_label_data$label_x, y = vert_label_data$label_y, 
                       size = label_size, colour = diverging_palette()[3],
                       hjust = "center", vjust = "bottom", angle = -90)
  for (i in seq_along(sub_plots)) {
    matrix_plot <- matrix_plot + cowplot::draw_plot(sub_plots[[i]], 
                                                    x = matrix_data[i, "x"],
                                                    y = matrix_data[i, "y"],
                                                    width = subgraph_width,
                                                    height = subgraph_height)
  }
  
  return(matrix_plot)
}
