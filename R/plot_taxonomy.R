#===================================================================================================
#' Plot a taxonomic tree
#' 
#' Plots the distribution of values associated with a taxonomic classification.
#' Taxonomic classifications can have multiple roots, resulting in multiple trees on the same plot.
#' Sizes and colors of vertexes, edges, labels, and individual trees can be displayed relative to
#' numbers (e.g. taxon statistics, such as abundance).
#' The displayed range of colors and sizes can be explicitly defined or automatically genereated.
#' Various transforamtions can be applied to numbers sizes/colors are mapped to.
#' Several types of tree layout algorithms from \code{\link{igraph}} can be used. 
#' 
#' @param taxon_id The unique ids of taxa.
#' @param parent_id The unique id of supertaxon \code{taxon_id} is a part of.
#' 
#' @param vertex_label See details on labels.
#' Default: no labels.
#' @param edge_label See details on labels.
#' Default: no labels.
#' @param tree_label See details on labels.
#' The label to display above each graph.
#' The value of the root of each graph will be used.
#' Default: None.
#' 
#' @param vertex_size See details on size.
#' Default: constant size.
#' @param edge_size See details on size.
#' Default: relative to vertex size. 
# #' @param tree_size See details on size.
# #' The value of the root of each graph will be used.
# #' This scales the space used to display graphs, but does not effect vertex/edge size.
# #' Default: Not used. 
#' 
#' @param vertex_label_size See details on size.
#' Default: relative to veterx size.
#' @param edge_label_size See details on size.
#' Default: relative to edge size.
#' @param tree_label_size See details on size.
#' Default: relative to graph size.
#' 
#' @param vertex_color See details on colors.
#' Default: grey.
#' @param edge_color See details on colors.
#' Default: same as vertex color.
#' @param tree_color See details on colors.
#' The value of the root of each graph will be used.
#' Overwrites the vertex and edge color if specified.
#' Default: Not used.
#' 
#' @param vertex_label_color See details on colors.
#' Default: black.
#' @param edge_label_color See details on colors.
#' Default: black.
#' @param tree_label_color See details on colors.
#' Default: black.
#' 
#' @param vertex_size_trans See details on transformations.
#' Default: \code{"area"}.
#' @param edge_size_trans See details on transformations. 
#' Default: same as \code{vertex_size_trans}. 
# #' @param tree_size_trans See details on transformations.
# #' Default: \code{"area"}.
#' 
#' @param vertex_label_size_trans See details on transformations. 
#' Default: same as \code{vertex_size_trans}.
#' @param edge_label_size_trans See details on transformations. 
#' Default: same as \code{edge_size_trans}.
#' @param tree_label_size_trans See details on transformations.
#' Default: \code{"area"}.
#' 
#' @param vertex_color_trans See details on transformations. 
#' Default: \code{"area"}.
#' @param edge_color_trans See details on transformations.
#' Default: same as vertex color transformation.
#' @param tree_color_trans See details on transformations.
#' Default: \code{"area"}.
#' 
#' @param vertex_label_color_trans See details on transformations.
#' Default: \code{"area"}.
#' @param edge_label_color_trans See details on transformations.
#' Default: \code{"area"}.
#' @param tree_label_color_trans See details on transformations. 
#' Default: \code{"area"}.
#' 
#' @param vertex_size_range See details on ranges.
#' Defualt: Optimize to balance overlaps and range size.
#' @param edge_size_range See details on ranges.
#' Default: relative to vertex size range. 
# #' @param tree_size_range See details on ranges.
# #' Default: Not set.
#' 
#' @param vertex_label_size_range See details on ranges.
#' Default: relative to vertex size. 
#' @param edge_label_size_range See details on ranges.
#' Default: relative to edge size.
#' @param tree_label_size_range See details on ranges.
#' Default: relative to tree size.
#' 
#' @param vertex_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' @param edge_color_range See details on ranges.
#' Default: same as vertex color.
#' @param tree_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' 
#' @param vertex_label_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' @param edge_label_color_range See details on ranges.
#' Default: Color-blind friendly palette.
#' @param  tree_label_color_range See details on ranges.
#' Default: Color-blind friendly palette. 
#' 
#' @param vertex_size_interval See details on intervals.
#' Default: The range of values in \code{vertex_size}. 
#' @param vertex_color_interval See details on intervals.
#' Default: The range of values in \code{vertex_color}. 
#' @param edge_size_interval See details on intervals.
#' Default: The range of values in \code{edge_size}. 
#' @param edge_color_interval See details on intervals.
#' Default: The range of values in \code{edge_color}. 
#' 
#' 
#' @param vertex_label_max The maximum number of vertex labels.
#' Default: 20.
#' @param edge_label_max The maximum number of edge labels.
#' Default: 20.
#' @param tree_label_max The maximum number of tree labels.
#' Default: 20.
#' 
#' @param overlap_avoidance (\code{numeric})
#' The relative importance of avoiding overlaps vs maximizing size range.
#' Higher numbers will cause vertex size optimazation to avoid overlaps more.
#' Default: \code{1}.
#' 
#' @param margin_size (\code{numeric} of length 2)
#' The horizontal and vertical margins.
#' Default: \code{0, 0}.
#' 
#' @param layout The layout algorithm used to position vertexes.
#' See details on layouts.
#' Default: \code{"reingold-tilford"}.
#' @param initial_layout he layout algorithm used to set the initial position
#' of vertexes, passed as input to the \code{layout} algorithm.
#' See details on layouts.
#' Default: Not used.
#' @param make_legend if TRUE...
#' @param title Name to print above the graph.
#' @param title_size The size of the title realtive to the rest of the graph. 
#' 
#' @param vertex_color_axis_label The label on the scale axis corresponding to \code{vertex_color}.
#' Default: The expression given to \code{vertex_color}.
#' @param vertex_size_axis_label The label on the scale axis corresponding to \code{vertex_size}.
#' Default: The expression given to \code{vertex_size}.
#' @param edge_color_axis_label The label on the scale axis corresponding to \code{edge_color}.
#' Default: The expression given to \code{edge_color}.
#' @param edge_size_axis_label The label on the scale axis corresponding to \code{edge_size}.
#' Default: The expression given to \code{edge_size}.
#' 
#' @param background_color The background color of the plot.
#' Default: Transparent
#' @param output_file The path to a file to save the plot in using \code{\link[ggplot2]{ggsave}}. 
#' The type of the file will be determined by the extension given.
#' Default: Do not save plot.
#' 
#' @param ... (other named arguments)
#' Passed to the \code{\link{igraph}} layout function used.
#' 
#' 
#' @section size:
#' 
#' 
#' The size of vertexes, edges, labels, and trees can be mapped to arbitrary numbers.
#' This is useful for displaying statistics for taxa, such as abundance.
#' Only the relative size of numbers is used, not the values themeselves.
#' They can be transformed to make the mapping non-linear using the transformation options.
#' The range of actual sizes displayed on the graph can be set using the range options.
#' 
#' Accepts a \code{numeric} vector, the same length \code{taxon_id} or a
#' factor of its length.
#' 
#' @section colors:
#' 
#' The colors of vertexes, edges, labels, and trees can be mapped to arbitrary numbers.
#' This is useful for highlighting groups of taxa.
#' Only the relative size of numbers is used, not the values themeselves.
#' They can be transformed to make the mapping non-linear using the transformation options.
#' The range of actual colors displayed on the graph can be set using the range options.
#' 
#' Accepts a vector, the same length \code{taxon_id} or a factor of its length.
#' If a numeric vector is given, it is mapped to a color scale.
#' Hex values or color names can be used (e.g. \code{#000000} or \code{"black"}).
#' 
#' @section labels:
#' 
#' The labels of vertexes, edges, and trees can be added.
#' Vertex labels are centered over their vertex.
#' Edge labels are displayed over edges, in the same orientation.
#' Tree labels are displayed over their tree.
#' 
#' Accepts a vector, the same length \code{taxon_id} or a factor of its length.
#' 
#' @section transformations:
#' 
#' Before any numbers specified are mapped to color/size, they can be transformed to make
#' the mapping non-linear. 
#' Any of the transformations listed below can be used by specifying their name.
#' A customized function can also be supplied to do the transformation.
#' 
#' \describe{
#'   \item{"radius"}{Proprotional to radius/diameter of vertex}
#'   \item{"area"}{circular area; better perceptual accuracy than \code{"radius"}}
#'   \item{"log10 radius"}{Log base 10 of radius}
#'   \item{"log2 radius"}{Log base 2 of radius}
#'   \item{"ln radius"}{Log base e of radius}
#'   \item{"log10 area"}{Log base 10 of circular area}
#'   \item{"log2 area"}{Log base 2 of circular area}
#'   \item{"ln area"}{Log base e of circular area}
#' }
#' 
#' @section ranges:
#' 
#' The displayed range of colors and sizes can be explicitly defined or automatically genereated.
#' Size ranges are specified by supplying a \code{numeric} vector with two values: the minimum and maximum.
#' The units used should be between 0 and 1, representing the proportion of a dimension of the graph.
#' Since the dimensions of the graph are determined by layout, and not always square, the value
#' that \code{1} corresponds to is the square root of the graph area (i.e. the side of a square with 
#' the same area as the plotted space).
#' Color ranges can be any number of color values as either HEX codes (e.g. \code{#000000}) or
#' color names (e.g. \code{"black"}).
#' 
#' @section layout:
#' 
#' Layouts determine the position of vertexes on the graph.
#' The are implemented using the \code{\link{igraph}} package.
#' Any additional arguments passed to \code{plot_taxonomy} are passed to the  \code{\link{igraph}}
#' function used.
#' The following \code{character} values are understood:
#' 
#' \describe{
#'   \item{"automatic"}{Use \code{\link[igraph]{nicely}}. Let \code{\link{igraph}} choose the layout.}
#'   \item{"reingold-tilford"}{Use \code{\link[igraph]{as_tree}}. A circular tree-like layout.}
#'   \item{"davidson-harel"}{Use \code{\link[igraph]{with_dh}}. A type of simulated annealing.}
#'   \item{"gem"}{Use \code{\link[igraph]{with_gem}}. A force-directed layout.}
#'   \item{"graphopt"}{Use \code{\link[igraph]{with_graphopt}}. A force-directed layout.}
#'   \item{"mds"}{Use \code{\link[igraph]{with_mds}}. Multidimensional scaling.}
#'   \item{"fruchterman-reingold"}{Use \code{\link[igraph]{with_fr}}. A force-directed layout.}
#'   \item{"kamada-kawai"}{Use \code{\link[igraph]{with_kk}}. A layout based on a phyisical model of springs.}
#'   \item{"large-graph"}{Use \code{\link[igraph]{with_lgl}}. Meant for larger graphs.}
#'   \item{"drl"}{Use \code{\link[igraph]{with_drl}}. A force-directed layout.}
#' }
#' 
#' 
#' @section intervals:
#' 
#' This is the minimum and maximum of values displayed on the legend scales.
#' Intervals are specified by supplying a \code{numeric} vector with two values: the minimum and maximum.
#' These are defined in the same units as element size/color.
#' By default, the minimum and maximum equals the range of the values used to infer size/color.
#' Setting a custom interval is useful for making size/color in multiple graphs correspond to the same statistics,
#' or setting logical bounderies (such as \code{c(0,1)} for proportions.
#' Note that this is different from the "range" options, which determine the range of graphed sizes/colors.
#'  
#' @examples 
#' plot(contaminants,
#'      vertex_size = item_counts,
#'      vertex_color = item_counts,
#'      vertex_label = name,
#'      tree_label = name,
#'      layout = "fruchterman-reingold")
#' 
#' @keywords internal
#' @rdname plot_taxonomy
plot_taxonomy <- function(taxon_id, parent_id, 
                          vertex_label = NA,
                          edge_label = NA,
                          tree_label = NA,
                          
                          vertex_size = 1,
                          edge_size = vertex_size,
                          # tree_size = 1,
                          
                          vertex_label_size = vertex_size,
                          edge_label_size = edge_size,
                          tree_label_size = as.numeric(NA), 
                          
                          vertex_color = "#999999",
                          edge_color = vertex_color,
                          tree_color = NA,
                          
                          vertex_label_color = "#000000",
                          edge_label_color = "#000000",
                          tree_label_color = "#000000",
                          
                          vertex_size_trans = "area",
                          edge_size_trans = vertex_size_trans,
                          # tree_size_trans = "area",
                          
                          vertex_label_size_trans = vertex_size_trans,
                          edge_label_size_trans = edge_size_trans,
                          tree_label_size_trans = "area",
                          
                          vertex_color_trans = "area",
                          edge_color_trans = vertex_color_trans,
                          tree_color_trans = "area",
                          
                          vertex_label_color_trans = "area",
                          edge_label_color_trans = "area",
                          tree_label_color_trans = "area",
                          
                          vertex_size_range = c(NA, NA),
                          edge_size_range = c(NA, NA),
                          # tree_size_range = c(NA, NA),
                          
                          vertex_label_size_range = c(NA, NA),
                          edge_label_size_range = c(NA, NA),
                          tree_label_size_range = c(NA, NA),
                          
                          vertex_color_range = quantative_palette(),
                          edge_color_range = vertex_color_range,
                          tree_color_range = quantative_palette(),
                          
                          vertex_label_color_range = quantative_palette(),
                          edge_label_color_range = quantative_palette(),
                          tree_label_color_range = quantative_palette(),
                          
                          vertex_size_interval = range(vertex_size, na.rm = TRUE, finite = TRUE),
                          vertex_color_interval = NULL,
                          edge_size_interval = range(edge_size, na.rm = TRUE, finite = TRUE),
                          edge_color_interval = NULL,
                          
                          vertex_label_max = 40,
                          edge_label_max = 40,
                          tree_label_max = 40,
                          
                          overlap_avoidance = 1,
                          margin_size = c(0, 0),
                          layout = "reingold-tilford",
                          initial_layout = "fruchterman-reingold",
                          make_legend = TRUE,
                          title = NULL,
                          title_size = 0.08,
                          
                          vertex_color_axis_label = NULL, 
                          vertex_size_axis_label = NULL,
                          edge_color_axis_label = NULL, 
                          edge_size_axis_label = NULL,
                          
                          background_color = "#FFFFFF00",
                          output_file = NULL,
                          
                          ...) {
  #| ### Verify arguments =========================================================================
  if (length(taxon_id) != length(parent_id)) {
    stop("'taxon_id' and 'parent_id' must be of equal length.")
  }
  if (length(taxon_id) == 0) {
    stop("'taxon_id' and 'parent_id' are empty.")
  }
  if (length(unique(taxon_id)) != length(taxon_id)) {
    stop("All values of 'taxon_id' are not unique.")
  }
  check_element_length(c("vertex_size", "edge_size",# "tree_size",
                         "vertex_label_size", "edge_label_size",  "tree_label_size",
                         "vertex_color", "edge_color", "tree_color",
                         "vertex_label_color", "edge_label_color", "tree_label_color",
                         "vertex_label", "edge_label", "tree_label"))
  verify_size(c("vertex_size", "edge_size", #"tree_size",
                "vertex_label_size", "edge_label_size", "tree_label_size"))
  verify_size_range(c("vertex_size_range",  "edge_size_range", # "tree_size_range",
                      "vertex_label_size_range", "edge_label_size_range", "tree_label_size_range",
                      "vertex_size_interval", "edge_size_interval"))
  verify_trans(c("vertex_size_trans", "edge_size_trans", #"tree_size_trans",
                 "vertex_color_trans", "edge_color_trans", "tree_color_trans",
                 "vertex_label_size_trans", "edge_label_size_trans", "tree_label_size_trans", 
                 "vertex_label_color_trans", "edge_label_color_trans", "tree_label_color_trans"))
  verify_color_range(c("vertex_color_range", "edge_color_range", "tree_color_range",
                       "vertex_label_color_range", "edge_label_color_range", "tree_label_color_range"))
  verify_label_count(c("vertex_label_max", "edge_label_max", "tree_label_max"))
  if (length(overlap_avoidance) == 0 || ! is.numeric(overlap_avoidance)) {
    stop("Argument 'overlap_avoidance' must be a numeric of length 1.")
  }
  if (length(margin_size) != 2 || ! is.numeric(margin_size)) {
    stop("Argument 'margin_size' must be a numeric of length 2.")
  }
  layout <- match.arg(layout, layout_functions())
  if (!is.null(initial_layout)) {
    initial_layout <- match.arg(initial_layout, layout_functions())
  }
  
  #| ### Parse arguments
  
  if (is.null(vertex_color_interval)) {
    vertex_color_interval <- range(vertex_color, na.rm = TRUE, finite = TRUE)
  }
  if (is.null(edge_color_interval)) {
    edge_color_interval <- range(edge_color, na.rm = TRUE, finite = TRUE)
  }
  
  #| ### Standardize source data ==================================================================
  data <- data.frame(stringsAsFactors = FALSE,
                     tid_user = as.character(taxon_id),
                     pid_user = as.character(parent_id),
                     
                     vl_user = as.character(vertex_label),
                     el_user = as.character(edge_label),
                     tl_user = as.character(tree_label),
                     
                     vs_user = as.numeric(vertex_size),
                     es_user = as.numeric(edge_size),
                     # ts_user = as.numeric(tree_size),
                     
                     vls_user = as.numeric(vertex_label_size),
                     els_user = as.numeric(edge_label_size),
                     tls_user = as.numeric(tree_label_size),
                     
                     vc_user = vertex_color,
                     ec_user = edge_color,
                     tc_user = tree_color,
                     
                     vlc_user = vertex_label_color,
                     elc_user = edge_label_color,
                     tlc_user = tree_label_color)
  row.names(data) <- data$tid_user
  
  #| #### Apply statistic transformations =========================================================
  trans_key <- c(vs_user = vertex_size_trans, es_user = edge_size_trans, #ts_user = tree_size_trans,
                 vls_user = vertex_label_size_trans, els_user = edge_label_size_trans,  tls_user = tree_label_size_trans,
                 vc_user = vertex_color_trans, ec_user = edge_color_trans, tc_user = edge_color_trans,
                 vlc_user = vertex_label_color_trans, elc_user = edge_label_color_trans, tlc_user = tree_label_color_trans)
  transformed_names <- gsub(pattern = "_user$", x = names(trans_key), replacement = "_trans")
  apply_trans <- function(col_name) {
    if (is.numeric(data[ , col_name])) { 
      transform_data(trans_key[col_name], data[ , col_name]) # if numbers are supplied
    } else {
      data[ , col_name] # if colors are defined explicitly, then no transformation is done
    }
  }
  data[, transformed_names] <- lapply(names(trans_key), apply_trans)
  # transform intervals
  vertex_size_interval_trans <- transform_data(vertex_size_trans, vertex_size_interval)
  edge_size_interval_trans <- transform_data(edge_size_trans, edge_size_interval)
  vertex_color_interval_trans <- transform_data(vertex_color_trans, vertex_color_interval)
  edge_color_interval_trans <- transform_data(edge_color_trans, edge_color_interval)
  
  
  #| ### Make layout ==============================================================================
  #| The layout is used to generate a list of coordinates to places graph verticies
  #| First the edge list consituted by the `taxon_id` and `parent_id` columns is used to construct 
  #| an `igraph` graph object and then the layout is generated for that object. 
  #|
  #| #### Make a graph for each root in the graph -------------------------------------------------
  get_sub_graphs <- function(taxa) {
    if (length(taxa) == 1) {
      # Make a graph with only a single vertex
      adj_matrix <- matrix(c(0), ncol = 1, dimnames =  list(taxa, taxa))
      sub_graph <- igraph::graph.adjacency(adj_matrix)
    } else {
      # Make edge list from taxon_id and parent_id
      edgelist <- as.matrix(data[taxa, c("pid_user", "tid_user")])
      # Remove edges to taxa that dont exist in this subset of the dataset
      edgelist <- edgelist[! is.na(edgelist[, "pid_user"]), , drop = FALSE]
#       # Randomly resort if layout is "reingold-tilford". NOTE: This is kinda hackish and should be replaced
#       if (layout == "reingold-tilford") { 
#         grouped_index <- split(rownames(edgelist), f = edgelist[, "pid_user"])
#         grouped_index <- unlist(grouped_index[sample(seq_along(grouped_index))])
#         edgelist <- edgelist[grouped_index, , drop = FALSE] 
#       }
      sub_graph <- igraph::graph_from_edgelist(edgelist)
    }
    igraph::V(sub_graph)$weight_factor <- data[taxa, c("vs_trans")]
    edge_end_vertex <- gsub("^[0-9]+\\|", "", attr(igraph::E(sub_graph), "vnames"))
    igraph::E(sub_graph)$weight_factor <- data[edge_end_vertex, c("vs_trans")]
    return(sub_graph)
  }
  data$is_root <- !(data$pid_user %in% data$tid_user)
  data[data$is_root, "pid_user"] <- NA # Needed by split_by_level
  sub_graph_taxa <- split_by_level(data$tid_user, data$pid_user, level =  1)
  sub_graphs <- lapply(sub_graph_taxa, get_sub_graphs)
  #|
  #| #### Generate a layout for each graph --------------------------------------------------------
  #|
  get_sub_layouts <- function(graph, backup_layout = 'fruchterman-reingold') {
    # Calculate an initial layout if specified
    if (! is.null(initial_layout) && layout != initial_layout) {
      intitial_coords <- layout_functions(initial_layout, graph)
      intitial_coords <- rescale(intitial_coords, to = c(-100, 100))
    } else {
      intitial_coords <- NULL
    }
    # Calculate the primary layout 
    coords <- layout_functions(layout, graph, intitial_coords = intitial_coords, ...)
    # Calculate backup layout if primary one does not work
    if (any(is.na(coords) | is.nan(unlist(coords)))) {
      coords <- layout_functions(backup_layout, graph)
      warning(paste0("Could not apply layout '", layout,
                     "' to subgraph. Using 'fruchterman-reingold' instead."))
    }
    return(coords)
  }
  
  sub_coords <- lapply(sub_graphs, get_sub_layouts)
  subgraph_key <- stats::setNames(rep(names(sub_graph_taxa), vapply(sub_graph_taxa, length, numeric(1))),
                           unlist(sub_graph_taxa))
  data$subgraph_root <- subgraph_key[data$tid_user]
  #|
  #| #### Merge layout coordinates into an overall graph ------------------------------------------
  #|
  coords <- igraph::merge_coords(sub_graphs, sub_coords) # merge vertex coordinates for each tree
  graph <- igraph::disjoint_union(sub_graphs) # merge graphs of each tree
  row.names(coords) <- names(igraph::V(graph))
  data$vx_plot <- coords[data$tid_user, 1]
  data$vy_plot <- coords[data$tid_user, 2]
  
  #| ### Core plot data ===========================================================================
  #|
  #| #### Optimize vertex size range --------------------------------------------------------------
  #|
  # Get range of potential vertex size ranges - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (nrow(data) > 1) {
    all_pairwise <- molten_dist(x = data$vx_plot, y = data$vy_plot) # get distance between all vertexes
    x_diff <- max(data$vx_plot) - min(data$vx_plot)
    y_diff <- max(data$vy_plot) - min(data$vy_plot)
    square_side_length <- sqrt(x_diff * y_diff)
    if (is.na(vertex_size_range[1])) { # if minimum vertex size not set
      min_range <- c(square_side_length / 1000, min(all_pairwise$distance))
    } else {
      min_range <- rep(vertex_size_range[1], 2) * square_side_length
    }
    if (is.na(vertex_size_range[2])) { # if maximum vertex size not set
      max_range <- c(min_range[1], square_side_length / 4)
    } else {
      max_range <- c(vertex_size_range[2], 2) * square_side_length
    }
    # Subset pairwise pairs to increase speed - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    max_important_pairs <- 2000 # Takes into account both size and distance
    max_biggest_pairs <- 2000
    max_closest_pairs <- 2000
    if (nrow(all_pairwise) > sum(max_important_pairs, max_biggest_pairs, max_closest_pairs)) {
      all_pairwise$size_sum <- data$vs_trans[all_pairwise$index_1] + data$vs_trans[all_pairwise$index_2]
      all_pairwise$importance <- all_pairwise$size_sum / all_pairwise$distance
      # all_pairwise <- all_pairwise[order(all_pairwise$importance, decreasing = TRUE), ]
      pair_subset <- c(order(all_pairwise$importance, decreasing = TRUE)[1:max_important_pairs],
                       order(all_pairwise$size_sum, decreasing = TRUE)[1:max_biggest_pairs],
                       order(all_pairwise$distance)[1:max_closest_pairs])
      all_pairwise <- all_pairwise[pair_subset, ]
    }
    # Define search space for potential vertex size ranges  - - - - - - - - - - - - - - - - - - - - -
    get_search_space <- function(min_range, max_range, breaks_per_dim) {
      min_breaks <- seq(from = min_range[1], to = min_range[2], length.out = breaks_per_dim)
      max_breaks <- seq(from = max_range[1], to = max_range[2], length.out = breaks_per_dim)
      max_breaks <- rescale(max_breaks^4, to = max_range)
      space <- data.frame(min = rep(min_breaks, each = length(max_breaks)),
                          max = rep(max_breaks, length(min_breaks)))
      space[space$min <= space$max, ]
    }
    search_space <- get_search_space(min_range, max_range, breaks_per_dim = 35)
    search_space$range_size <- search_space$max - search_space$min
    # Calculate vertex overlap resulting from possible ranges - - - - - - - - - - - - - - - - - - - -
    find_overlap <- function(a_min, a_max, distance) {
      scaled_vs <- rescale(data$vs_t, to = c(a_min, a_max), from = vertex_size_interval_trans)
      names(scaled_vs) <- data$tid_user
      gap <- distance$distance - scaled_vs[distance$index_1] - scaled_vs[distance$index_2]
      gap <- ifelse(gap < 0, abs(gap), 0)
      gap <- gap^2
      sqrt(sum(gap))
    } 
    search_space$overlap <- apply(search_space, MARGIN = 1,
                                  function(x) find_overlap(x["min"], x["max"], all_pairwise))
    
    # Choose base range based on optimality criteria  - - - - - - - - - - - - - - - - - - - - - - - -
    optimality_stat <- function(overlap, range_size, minimum) {
      overlap_weight <- 1 / nrow(data)^(1/3)
      minimum_weight <- 75 / sqrt(nrow(data))
      (1 + range_size + minimum * minimum_weight) / (1 + overlap * overlap_avoidance * overlap_weight)
    }
    search_space$opt_stat <- apply(search_space, MARGIN = 1,
                                   function(x) optimality_stat(x["overlap"], x["range_size"], x["min"]))
    vsr_plot <- unlist(search_space[which.max(search_space$opt_stat), c("min", "max")])
  } else {
    square_side_length = 1
    vsr_plot <- rep(square_side_length / 4, 2)
  }
  data$vs_plot <- rescale(data$vs_t, to = vsr_plot, from = vertex_size_interval_trans)
  #|
  #| #### Infer edge size range -------------------------------------------------------------------
  #|
  infer_size_range <- function(specified_range, reference_range, defualt_scale) {
    result <- specified_range * square_side_length
    if (is.na(result[1]) && is.na(result[2])) { # If the user has not set range
      result <- reference_range * defualt_scale
    } else if (is.na(result[1])) { # If the user has set a maximum but not a minimum
      result[1] <- result[2] * reference_range[1] / reference_range[2]
    } else if (is.na(result[2])) { # If the user has set a minimum but not a maximum
      result[2] <- result[1] * reference_range[2] / reference_range[1]
    }
    return(result)
  }
  
  esr_plot <- infer_size_range(edge_size_range, vsr_plot, defualt_scale = 0.5)
  data$es_plot <- rescale(data$es_t, to = esr_plot, from = edge_size_interval_trans)
  #|
  #| #### Infer tree size range -------------------------------------------------------------------
  #|
  get_tree_area <- function(a_root) {
    size <- data[data$subgraph_root == a_root, "vs_plot"]
    x <- data[data$subgraph_root == a_root, "vx_plot"]
    x <- c(x + size, x - size)
    y <- data[data$subgraph_root == a_root, "vy_plot"]
    y <- c(y + size, y - size)
    (max(x) - min(x)) * (max(y) - min(y)) 
  }
  tree_area <- vapply(unique(data$subgraph_root), get_tree_area, FUN.VALUE = numeric(1))
  data$tree_area <- tree_area[data$subgraph_root]
  tsr_plot <- range(sqrt(tree_area))
  #|
  #| #### Infer label size ranges -----------------------------------------------------------------
  #|
  if (all(is.na(data$tls_user))) {
    data$tls_user <- sqrt(data$tree_area)
    data$tls_trans <- apply_trans("tls_user") 
  }
  vlsr_plot <- infer_size_range(vertex_label_size_range, vsr_plot, defualt_scale = 0.8)
  elsr_plot <- infer_size_range(edge_label_size_range, esr_plot, defualt_scale = 0.8)
  tlsr_plot <- infer_size_range(tree_label_size_range, tsr_plot, defualt_scale = 0.1)
  data$vls_plot <- rescale(data$vls_trans, to = vlsr_plot)
  data$els_plot <- rescale(data$els_trans, to = elsr_plot)
  data$tls_plot <- rescale(data$tls_trans, to = tlsr_plot)
  #|
  #| #### Assign color scales ---------------------------------------------------------------------
  #|
  
  color_colume_key <- list("ec_trans" = edge_color_range, "vc_trans" = vertex_color_range, 
                           "tc_trans" = tree_color_range, "vlc_trans" = vertex_label_color_range,
                           "elc_trans" = edge_label_color_range, "tlc_trans" = tree_label_color_range)
  color_interval_key <- list("ec_trans" = edge_color_interval_trans, "vc_trans" = vertex_color_interval_trans)
  plot_value_names <- gsub(pattern = "_trans$", x = names(color_colume_key), replacement = "_plot")
  data[, plot_value_names] <- lapply(names(color_colume_key),
                                     function(x) apply_color_scale(data[ , x],
                                                                   color_colume_key[[x]],
                                                                   interval = color_interval_key[[x]]))
  # If tree_color is used, overwrite other colors - - - - - - - - - - - - - - - - - - - - - - - - -
  data$tc_plot <- data[data$subgraph_root, "tc_plot"]
  to_replace <- ! is.na(data$tc_plot)
  data[to_replace, "vc_plot"] <- data[to_replace, "tc_plot"]
  data[to_replace, "ec_plot"] <- data[to_replace, "tc_plot"]
  
  #| ### Secondary plot data ======================================================================
  #|
  #| #### Calculate coordinants of graph elements -------------------------------------------------
  #| The vertexes and edges must be specified by a dataframe of coordinates, with a colume 
  #| grouping the coordinates of each shape.
  #| These shapes must be added to the graph in a specific order.
  #| A list of vertexes is sorted by first vertex depth in the heirarchy and then by vertex size.
  taxon_elements <- function(tid) {
    circle_resolution <- 25
    edge_data <- line_coords(x1 = data[tid, 'vx_plot'],
                             y1 = data[tid, 'vy_plot'],
                             x2 = data[data[tid, 'pid_user'], "vx_plot"],
                             y2 = data[data[tid, 'pid_user'], "vy_plot"],
                             width = data[tid, 'es_plot'] * 2)
    edge_data$group <- paste0(tid, "_edge")
    edge_data$color <- rep(data[tid, 'ec_plot'], each = 4)
    vertex_data <- polygon_coords(n = circle_resolution,
                                  x = data[tid, 'vx_plot'],
                                  y = data[tid, 'vy_plot'],
                                  radius = data[tid, 'vs_plot'])
    vertex_data$group <- paste0(tid, "_vertex")
    vertex_data$color <- rep(data[tid, 'vc_plot'], each = circle_resolution + 1)
    output <- rbind(edge_data, vertex_data)
    # output$tid_user <- tid
    return(output[stats::complete.cases(output),])
  }
  data$level = edge_list_depth(data$tid_user, data$pid_user)
  element_order <- data$tid_user[order(data$level, 1 / data$vs_plot, decreasing = TRUE)]
  element_data <- do.call(rbind, lapply(element_order, taxon_elements))
  element_data$group <- factor(element_data$group, levels = unique(element_data$group))
  #|
  #| #### Make text data ------------------------------------------------------------------
  #|
  # Get vertex label data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$vl_is_shown <- select_labels(data, vertex_label_max,
                                    sort_by_column = c("vls_plot", "vs_plot"),
                                    label_column = "vl_user")
  if (any(data$vl_is_shown)) {
    vl_data <- data[data$vl_is_shown, , drop = FALSE]
    text_data <- data.frame(stringsAsFactors = FALSE,
                            label = vl_data$vl_user,
                            x = vl_data$vx_plot,
                            y = vl_data$vy_plot,
                            size = vl_data$vls_plot,
                            color = vl_data$vlc_plot,
                            rotation = 0,
                            justification = "center")
  } else {
    text_data <- NULL
  }
  # Get edge label data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$el_is_shown <- select_labels(data, edge_label_max,
                                    sort_by_column = c("els_plot", "es_plot"),
                                    label_column = "el_user")
  data[is.na(data$pid_user), "el_is_shown"] <- FALSE # taxa with no parents get no line label
  if (any(data$el_is_shown)) {
    el_data <- data[data$el_is_shown, ]
    # edge label rotation 
    el_data$el_slope <- (el_data$vy_plot - data[el_data$pid_user, "vy_plot"]) / (el_data$vx_plot - data[el_data$pid_user, "vx_plot"])
    el_data$el_slope[is.na(el_data$el_slope)] <- 0
    el_data$el_rotation <- atan(el_data$el_slope)
    # edge label coordinate 
    line_label_offset = 1
    justify <- data[el_data$pid_user, "vx_plot"] > el_data$vx_plot
    justify[is.na(justify)] <- TRUE
    justification <- ifelse(justify, "left-center", "right-center")
    line_label_x_offset <- line_label_offset * el_data$vs_plot * cos(el_data$el_rotation)
    line_label_y_offset <- line_label_offset * el_data$vs_plot * sin(el_data$el_rotation)
    el_data$elx_plot <- el_data$vx_plot + ifelse(justify, 1, -1) * line_label_x_offset
    el_data$ely_plot <- el_data$vy_plot + ifelse(justify, 1, -1) * line_label_y_offset
    # create text data   
    text_data <- rbind(text_data,
                       data.frame(stringsAsFactors = FALSE, 
                                  label = el_data$el_user,
                                  x = el_data$elx_plot,
                                  y = el_data$ely_plot,
                                  size = el_data$els_plot,
                                  color = el_data$elc_plot,
                                  rotation = el_data$el_rotation,
                                  justification = justification))
  }
  # Get tree label data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$tl_is_shown <- FALSE
  data[data$is_root, "tl_is_shown"] <- select_labels(data[data$is_root, ], tree_label_max,
                                                     sort_by_column = c("tls_plot", "vs_plot"), label_column = "tl_user")
  if (any(data$tl_is_shown)) {
    title_data <- data[data$tl_is_shown, , drop = FALSE]
    tx_plot <- vapply(split(data$vx_plot, data$subgraph_root), FUN.VALUE = numeric(1),
                      function(x) mean(range(x)))
    title_data$tx_plot <- tx_plot[title_data$subgraph_root]
    ty_plot <- vapply(split(data$vy_plot, data$subgraph_root), FUN.VALUE = numeric(1),
                      function(y) mean(range(y)))
    title_data$ty_plot <- ty_plot[title_data$subgraph_root]
    title_data$tlx_plot <- title_data$tx_plot 
    tly_plot <- mapply(function(y, size) max(y + size),
                       y = split(data$vy_plot, data$subgraph_root),
                       size = split(data$vs_plot, data$subgraph_root))
    title_data$tly_plot <- tly_plot[title_data$subgraph_root] + title_data$tls_plot * 1.1
    text_data <- rbind(text_data,
                       data.frame(stringsAsFactors = FALSE, 
                                  label = title_data$tl_user,
                                  x = title_data$tlx_plot,
                                  y = title_data$tly_plot,
                                  size = title_data$tls_plot,
                                  color = title_data$tlc_plot,
                                  rotation = 0,
                                  justification = "center"))
  }
  # Add tree title data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (! is.null(title)) {
    max_label_size <- ifelse(is.null(text_data), 0, max(text_data$size))
    x_range <- range(element_data$x)
    
    text_data <- rbind(text_data,
                       data.frame(stringsAsFactors = FALSE, 
                                  label = title,
                                  x = mean(x_range),
                                  y = max(element_data$y) * 1.05 + diff(x_range) * title_size + max_label_size,
                                  size = diff(x_range) * title_size,
                                  color = "#000000",
                                  rotation = 0,
                                  justification = "center"))
  }
  #|
  #| #### Make vertex legend -----------------------------------------------------------------------
  #|
  y_range <- max(element_data$y) - min(element_data$y)
  legend_length <- square_side_length * 0.3
  right_plot_boundry <- max(element_data$x) * 1.1
  if (make_legend) {
    vertex_legend <- make_plot_legend(x = right_plot_boundry,
                                    y = min(element_data$y), 
                                    length = legend_length, 
                                    width_range = range(data$vs_plot) * 2, 
                                    width_stat_range =  vertex_size_interval,
                                    group_prefix = "vertex_legend",
                                    width_stat_trans = transform_data(func = vertex_size_trans),
                                    color_range = vertex_color_range,
                                    color_stat_range = vertex_color_interval, 
                                    color_stat_trans =  transform_data(func = vertex_color_trans),
                                    title = "Vertices",
                                    color_axis_label = vertex_color_axis_label,
                                    size_axis_label = vertex_size_axis_label,
                                    hide_size = missing(vertex_size),
                                    hide_color = missing(vertex_color))
  #|
  #| #### Make edge legend -----------------------------------------------------------------------
  #|
    edge_legend <- make_plot_legend(x = right_plot_boundry,
                                    y = max(element_data$y) - legend_length - 0.1 * legend_length, 
                                    length = legend_length, 
                                    width_range = range(data$es_plot) * 2, 
                                    width_stat_range =  edge_size_interval,
                                    group_prefix = "edge_legend",
                                    width_stat_trans = transform_data(func = edge_size_trans),
                                    color_range = edge_color_range,
                                    color_stat_range = edge_color_interval, 
                                    color_stat_trans =  transform_data(func = edge_color_trans),
                                    title = "Edges",
                                    color_axis_label = edge_color_axis_label,
                                    size_axis_label = edge_size_axis_label,
                                    hide_size = missing(edge_size),
                                    hide_color = missing(edge_color))
    element_data <- rbind(element_data, vertex_legend$shapes, edge_legend$shapes)
    text_data <- rbind(text_data, vertex_legend$labels, edge_legend$labels)
  } else {
    legend_data <- NULL
  }
  #| ### Draw plot ================================================================================
  
  label_x_bounds <- function(x, size, label) {
    spread <- size  * text_grob_length(label)
    c(x + spread, x - spread)
  }
  label_y_bounds <- function(y, size, label) {
    spread <- size / 2 * 1.1
    c(y + spread, y - spread)
  }
  
  
  
  label_x <- unlist(do.call(mapply, args = c(text_data[ , c("x", "size", "label")],
                                             FUN = label_x_bounds)))
  label_y <- unlist(do.call(mapply, args = c(text_data[ , c("y", "size", "label")],
                                             FUN = label_y_bounds)))
  
  x_points <- c(element_data$x, label_x)
  y_points <- c(element_data$y, label_y)
  margin_size_plot <- margin_size * square_side_length
  x_range <- c(min(x_points) - margin_size_plot[1], max(x_points) + margin_size_plot[1])
  y_range <- c(min(y_points) - margin_size_plot[2], max(y_points) + margin_size_plot[2])
  
  # theme(panel.background = element_rect(fill = '#00000000', colour = '#00000000'))
  
  the_plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_polygon(data = element_data, ggplot2::aes_string(x = "x", y = "y", group = "group"),
                          fill = element_data$color) +
    ggplot2::guides(fill = "none") +
    ggplot2::coord_fixed(xlim = x_range, ylim = y_range) +
    ggplot2::scale_y_continuous(expand = c(0,0), limits = y_range) +
    ggplot2::scale_x_continuous(expand = c(0,0), limits = x_range) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   panel.background = ggplot2::element_rect(fill = background_color, colour = background_color),
                   plot.background = ggplot2::element_rect(fill = background_color, colour = background_color),
                   axis.title = ggplot2::element_blank(),
                   axis.text  = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), 
                   axis.line  = ggplot2::element_blank(),
                   plot.margin = grid::unit(c(0,0,0,0) , "in"))
  if (!is.null(text_data)) {
    text_grobs <- do.call(make_text_grobs, c(text_data, list(x_range = x_range, y_range = y_range)))
    for (a_grob in text_grobs) {
      the_plot <- the_plot + ggplot2::annotation_custom(grob = a_grob)
    }
  }
  
  #| ### Save output file
  if (!is.null(output_file)) {
    img_width <- diff(x_range)
    img_height <- diff(y_range)
    ggplot2::ggsave(output_file, the_plot, bg = "transparent", width = 10, height = 10 * (img_height / img_width))
  }
  return(the_plot)
}


#' @param x An object of type \code{\link{classified}}
#' 
#' @method plot classified
#' @export
#' @rdname plot_taxonomy
plot.classified <- function(x, ...) {
  # Non-standard argument evaluation
  data <- taxon_data(x, sort_by = classifications, 
                     col_subset = unique(c(taxon_data_cols_used(x, ...), "taxon_ids", "parent_ids")))
  arguments <- c(list(taxon_id = data$taxon_ids, parent_id = data$parent_ids),
                 lazyeval::lazy_eval(lazyeval::lazy_dots(...), data = data))
  
  # Use variable name for scale axis labels
  if (! "vertex_color_axis_label" %in% names(arguments)) {
    arguments$vertex_color_axis_label <- deparse(as.list(match.call())$vertex_color)
  }
  if (! "vertex_size_axis_label" %in% names(arguments)) {
    arguments$vertex_size_axis_label <- deparse(as.list(match.call())$vertex_size)
  }
  if (! "edge_color_axis_label" %in% names(arguments)) {
    arguments$edge_color_axis_label <- deparse(as.list(match.call())$edge_color)
  }
  if (! "edge_size_axis_label" %in% names(arguments)) {
    arguments$edge_size_axis_label <- deparse(as.list(match.call())$edge_size)
  }
  
  # Call plot_taxonomy
  do.call(plot_taxonomy, arguments)
}
