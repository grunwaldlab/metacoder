#' Plot a matrix of heat trees
#' 
#' Plot a matrix of heat trees for showing pairwise comparisons. A larger,
#' labelled tree serves as a key for the matrix of smaller unlabelled trees. The
#' data for this function is typically created with \code{\link{compare_groups}},
#' 
#' @param obj A \code{\link[taxa]{taxmap}} object
#' @param data The name of a table in \code{obj$data} that is the output of 
#'   \code{\link{compare_groups}} or in the same format.
#' @param label_small_trees If \code{TRUE} add labels to small trees as well as 
#'   the key tree. Otherwise, only the key tree will be labeled.
#' @param key_size The size of the key tree relative to the whole graph. For
#'   example, 0.5 means half the width/height of the graph.
#' @param seed That random seed used to make the graphs.
#' @param output_file The path to one or more files to save the plot in using \code{\link[ggplot2]{ggsave}}. 
#' The type of the file will be determined by the extension given.
#' Default: Do not save plot.
#' @param row_label_color The color of the row labels on the right side of the matrix. Default:
#'   based on the node_color_range.
#' @param col_label_color The color of the columns labels along the top of the matrix. Default:
#'   based on the node_color_range.
#' @param row_label_size The size of the row labels on the right side of the matrix. Default: 12.
#' @param col_label_size The size of the columns labels along the top of the matrix. Default: 12.
#' @param ... Passed to \code{\link{heat_tree}}. Some options will be overwritten.
#' @param dataset DEPRECIATED. use "data" instead.
#' 
#' @examples
#' \dontrun{
#' # Parse dataset for plotting
#' x <- parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                     class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                     class_regex = "^(.+)__(.+)$")
#' 
#' # Convert counts to proportions
#' x$data$otu_table <- calc_obs_props(x, data = "tax_data", cols = hmp_samples$sample_id)
#' 
#' # Get per-taxon counts
#' x$data$tax_table <- calc_taxon_abund(x, data = "otu_table", cols = hmp_samples$sample_id)
#' 
#' # Calculate difference between treatments
#' x$data$diff_table <- compare_groups(x, data = "tax_table",
#'                                     cols = hmp_samples$sample_id,
#'                                     groups = hmp_samples$body_site)
#'
#' # Plot results (might take a few minutes)
#' heat_tree_matrix(x,
#'                  data = "diff_table",
#'                  node_size = n_obs,
#'                  node_label = taxon_names,
#'                  node_color = log2_median_ratio,
#'                  node_color_range = diverging_palette(),
#'                  node_color_trans = "linear",
#'                  node_color_interval = c(-3, 3),
#'                  edge_color_interval = c(-3, 3),
#'                  node_size_axis_label = "Number of OTUs",
#'                  node_color_axis_label = "Log2 ratio median proportions")
#' 
#' }
#' 
#' @export
heat_tree_matrix <- function(obj, data, label_small_trees =  FALSE,
                             key_size = 0.6, seed = 1, output_file = NULL,
                             row_label_color = diverging_palette()[3],
                             col_label_color = diverging_palette()[1], 
                             row_label_size = 12, col_label_size = 12,
                             ..., dataset = NULL) {
  
  # Check for use of "dataset"
  if (! is.null(dataset)) {
    warning(call. = FALSE,
            'Use of "dataset" is depreciated. Use "data" instead.')
    data <- dataset
  }
  
  # Get plot data table 
  diff_table <- get_taxmap_table(obj, data)
  
  # Make plot layout
  treat_1 <- as.character(diff_table$treatment_1)
  treat_2 <- as.character(diff_table$treatment_2)
  treatments <- unique(c(treat_1, treat_2))
  combinations <- t(utils::combn(seq_along(treatments), 2))
  layout_matrix <- matrix(rep(NA, (length(treatments))^2), nrow = length(treatments))
  for (index in 1:nrow(combinations)) {
    layout_matrix[combinations[index, 1], combinations[index, 2]] <- index
  }
  
  # Make individual plots
  plot_sub_plot <- ifelse(label_small_trees, # This odd thing is used to overwrite options without evaluation
    function(..., make_node_legend = FALSE, make_edge_legend = FALSE, output_file = NULL) {
      metacoder::heat_tree(..., make_node_legend = FALSE, make_edge_legend = FALSE, output_file = NULL)
    },
    function(..., node_label = NULL, make_node_legend = FALSE, make_edge_legend = FALSE, output_file = NULL) {
      metacoder::heat_tree(..., make_node_legend = FALSE, make_edge_legend = FALSE, output_file = NULL)
    }
  )
  
  sub_plots <- lapply(seq_len(nrow(combinations)),
                      function(index) {
                        set.seed(seed)
                        obj %>%
                          taxa::filter_obs(data,
                                           (treat_1 == treatments[combinations[index, 1]] &
                                              treat_2 == treatments[combinations[index, 2]]) |
                                             (treat_1 == treatments[combinations[index, 2]] &
                                                treat_2 == treatments[combinations[index, 1]])) %>%
                          plot_sub_plot(...)
                      })
  
  # Make key plot
  plot_key_plot <- function(..., node_color = NULL) {
    heat_tree(..., node_color = "grey")
  }
  
  set.seed(seed)
  key_plot <- plot_key_plot(obj, ...)
  
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
  matrix_plot <- cowplot::ggdraw() + 
    cowplot::draw_plot(key_plot, x = 0, y = 0, width = key_size, height = key_size) +
    cowplot::draw_text(gsub("_", " ", horz_label_data$treatment_2), 
                       x = horz_label_data$label_x, y = horz_label_data$label_y, 
                       size = col_label_size, colour = col_label_color,
                       hjust = "center", vjust = "bottom") +
    cowplot::draw_text(gsub("_", " ", vert_label_data$treatment_1), 
                       x = vert_label_data$label_x, y = vert_label_data$label_y, 
                       size = row_label_size, colour = row_label_color,
                       hjust = "center", vjust = "bottom", angle = -90) +
    ggplot2::theme(aspect.ratio = 1)
  for (i in seq_along(sub_plots)) {
    matrix_plot <- matrix_plot + cowplot::draw_plot(sub_plots[[i]], 
                                                    x = matrix_data[i, "x"],
                                                    y = matrix_data[i, "y"],
                                                    width = subgraph_width,
                                                    height = subgraph_height)
  }
  
  # Save plot
  if (!is.null(output_file)) {
    for (path in output_file) {
      ggplot2::ggsave(path, matrix_plot, bg = "transparent", width = 10, height = 10)
    }
  }
  
  
  return(matrix_plot)
}
