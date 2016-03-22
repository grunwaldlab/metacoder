#' Plot a taxonomic tree
#' 
#' Plot a taxonomic tree
#' 
#' @param classified_data An object of class \link{classified}.
#' @param ... Passed to \link{plot_taxonomy}
#' 
#' @inheritParams plot_taxonomy
#' 
#' @export
plot.classified <- function(classified_data, ...) {
  my_taxon_data <- taxon_data(classified_data)
  column_var_name <- colnames(my_taxon_data)
  unused_result <- lapply(column_var_name, function(x) assign(x, my_taxon_data[[x]], envir = parent.frame(2)))
  arguments <- c(list(taxon_id = classified_data$taxon_id, parent_id = classified_data$parent_id),
                 eval(substitute(list(...))))
  do.call(plot_taxonomy, arguments)
}







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
#' @param tree_size_trans See details on transformations.
#' Default: \code{"area"}.
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
#' @export
plot_taxonomy <- function(taxon_id, parent_id, 
                          vertex_label = NA,
                          edge_label = NA,
                          tree_label = NA,
                          
                          vertex_size = 1,
                          edge_size = vertex_size,
                          tree_size = 1,
                          
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
                          
                          vertex_label_max = 50,
                          edge_label_max = 50,
                          tree_label_max = 50,
                          
                          overlap_avoidance = 1,
                          margin_size = c(0.1, 0.05),
                          layout = "reingold-tilford",
                          initial_layout = "fruchterman-reingold",
                          make_legend = TRUE,
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
                      "vertex_label_size_range", "edge_label_size_range", "tree_label_size_range"))
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
  if (! layout %in% layout_functions()) {
    stop("Argument 'layout' must be an output of layout_functions().")
  }
  if (! initial_layout %in% layout_functions()) {
    stop("Argument 'initial_layout' must be an output of layout_functions().")
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
      transform_data(data[ , col_name], trans_key[col_name]) # if numbers are supplied
    } else {
      data[ , col_name] # if colors are defined explicitly, then no transformation is done
    }
  }
  data[, transformed_names] <- lapply(names(trans_key), apply_trans)
  
  
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
      return(igraph::graph.adjacency(adj_matrix))
    } else {
      # Make edge list from taxon_id and parent_id
      edgelist <- as.matrix(data[taxa, c("pid_user", "tid_user")])
      # Remove edges to taxa that dont exist in this subset of the dataset
      edgelist <- edgelist[! is.na(edgelist[, "pid_user"]), , drop = FALSE]
      return(igraph::graph_from_edgelist(edgelist))
    }
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
      intitial_coords <- igraph::layout_(graph, layout_functions(initial_layout, graph))
    } else {
      intitial_coords <- NULL
    }
    # Calculate the primary layout 
    coords <- igraph::layout_(graph, layout_functions(layout, graph = graph,
                                                      intitial_coords = intitial_coords, 
                                                      ...))
    # Calculate backup layout if primary one does not work
    if (any(is.na(coords) | is.nan(unlist(coords)))) {
      coords <- igraph::layout_(graph, layout_functions(backup_layout, graph = graph))
      warning(paste0("Could not apply layout '", layout,
                     "' to subgraph. Using 'fruchterman-reingold' instead."))
    }
    return(coords)
  }
  
  sub_coords <- lapply(sub_graphs, get_sub_layouts)
  data$subgraph_root <- rep(names(sub_coords), vapply(sub_coords, nrow, numeric(1)))
  # scaled_ts_trans <- scales::rescale(data$ts_trans, to = c(1, 2)) # make sure numbers are reasonable
  # sub_coords <- mapply(`*`, sub_coords, scaled_ts_trans[data$is_root]) # Scale coordinates by tree_size
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
  # Define search space for potential vertex size ranges  - - - - - - - - - - - - - - - - - - - - -
  get_search_space <- function(min_range, max_range, breaks_per_dim = 10) {
    min_breaks <- seq(from = min_range[1], to = min_range[2], length.out = breaks_per_dim)
    max_breaks <- seq(from = max_range[1], to = max_range[2], length.out = breaks_per_dim)
    max_breaks <- scales::rescale(max_breaks^4, to = max_range)
    space <- data.frame(min = rep(min_breaks, each = length(max_breaks)),
                        max = rep(max_breaks, length(min_breaks)))
    space[space$min <= space$max, ]
  }
  search_space <- get_search_space(min_range, max_range, breaks_per_dim = 10)
  search_space$range_size <- search_space$max - search_space$min
  # Calculate vertex overlap resulting from possible ranges - - - - - - - - - - - - - - - - - - - -
  find_overlap <- function(a_min, a_max, distance) {
    scaled_vs <- scales::rescale(data$vs_t, to = c(a_min, a_max))
    names(scaled_vs) <- data$tid_user
    gap <- distance$distance - scaled_vs[distance$index_1] - scaled_vs[distance$index_2]
    gap <- ifelse(gap < 0, abs(gap), 0)
    sum(gap)
  } 
  search_space$overlap <- apply(search_space, MARGIN = 1,
                                function(x) find_overlap(x["min"], x["max"], all_pairwise))
  
  # Choose base range based on optimality criteria  - - - - - - - - - - - - - - - - - - - - - - - -
  optimality_stat <- function(overlap, range_size, minimum) {
    overlap_weight <- 0.10
    minimum_weight <- 3
    (1 + range_size + minimum * minimum_weight) / (1 + overlap * overlap_avoidance * overlap_weight)
  }
  search_space$opt_stat <- apply(search_space, MARGIN = 1,
                                 function(x) optimality_stat(x["overlap"], x["range_size"], x["min"]))
  vsr_plot <- unlist(search_space[which.max(search_space$opt_stat), c("min", "max")])
  data$vs_plot <- scales::rescale(data$vs_t, to = vsr_plot)
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
  data$es_plot <- scales::rescale(data$es_t, to = esr_plot)
  #|
  #| #### Infer tree size range -------------------------------------------------------------------
  #|
  get_tree_area <- function(a_root) {
    x <- data[data$subgraph_root == a_root, "vx_plot"]
    y <- data[data$subgraph_root == a_root, "vy_plot"]
    (max(x) - min(x)) * (max(y) - min(y)) 
  }
  tree_area <- vapply(unique(data$subgraph_root), get_tree_area, FUN.VALUE = numeric(1))
  tree_vertex_counts <-  as.numeric(table(data$subgraph_root)[unique(data$subgraph_root)])
  data$tree_area <- rep(tree_area, tree_vertex_counts)
  tsr_plot <- range(sqrt(tree_area))
  #|
  #| #### Infer label size ranges -----------------------------------------------------------------
  #|
  if (all(is.na(data$tls_user))) {
    data$tls_user <- sqrt(data$tree_area)
    data$tls_trans <- apply_trans("tls_user") 
  }
  vlsr_plot <- infer_size_range(vertex_label_size_range, vsr_plot, defualt_scale = 0.5)
  elsr_plot <- infer_size_range(edge_label_size_range, esr_plot, defualt_scale = 0.7)
  tlsr_plot <- infer_size_range(tree_label_size_range, tsr_plot, defualt_scale = 0.1)
  data$vls_plot <- scales::rescale(data$vls_trans, to = vlsr_plot)
  data$els_plot <- scales::rescale(data$els_trans, to = elsr_plot)
  data$tls_plot <- scales::rescale(data$tls_trans, to = tlsr_plot)
  #|
  #| #### Assign color scales ---------------------------------------------------------------------
  #|
  
  color_colume_key <- list("ec_trans" = edge_color_range, "vc_trans" = vertex_color_range, 
                           "tc_trans" = tree_color_range, "vlc_trans" = vertex_label_color_range,
                           "elc_trans" = edge_label_color_range, "tlc_trans" = tree_label_color_range)
  plot_value_names <- gsub(pattern = "_trans$", x = names(color_colume_key), replacement = "_plot")
  data[, plot_value_names] <- lapply(names(color_colume_key),
                                     function(x) apply_color_scale(data[ , x], color_colume_key[[x]]))
  # If tree_color is used, overwrite other colors - - - - - - - - - - - - - - - - - - - - - - - - -
  data$tc_plot <- rep(data[data$is_root, "tc_plot"], tree_vertex_counts) # apply root color to trees
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
    circle_resolution <- 50
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
    output$tid_user <- tid
    return(output[complete.cases(output),])
  }
  data$level = edge_list_depth(data$tid_user, data$pid_user)
  element_order <- data$tid_user[order(data$level, 1 / data$vs_plot, decreasing = TRUE)]
  element_data <- do.call(rbind, lapply(element_order, taxon_elements))
  element_data$group <- factor(element_data$group, levels = unique(element_data$group))
  #|
  #| #### Make legend 
  #|
  if (make_legend && (!missing(vertex_size) || !missing(vertex_color))) {
    if (missing(vertex_size)) {
      width_stat_range <- NULL
    } else {
      width_stat_range <- range(vertex_size)
    }
    if (missing(vertex_color)) {
      color_stat_range <- NULL
    } else {
      color_stat_range <- range(vertex_color)
    }
    y_range <- max(element_data$y) - min(element_data$y)
    legend_data <- make_plot_legend(x = max(element_data$x) * 1.1,
                                    y = min(element_data$y), 
                                    length = y_range * 0.3, 
                                    tick_size = y_range * 0.003, 
                                    width_range = range(data$vs_plot) * 2, 
                                    width_stat_range = width_stat_range,
                                    width_stat_trans = transform_data(func = vertex_size_trans, inverse = TRUE),
                                    color_range = vertex_color_range,
                                    color_stat_range = color_stat_range, 
                                    color_stat_trans =  transform_data(func = vertex_color_trans, inverse = TRUE),
                                    divisions = 100, label_count = 6)
    legend_data$shapes$tid_user = NA
    element_data <- rbind(element_data, legend_data$shapes)
  } else {
    legend_data <- NULL
  }
  
  #|
  #| #### Make vertex text grobs ------------------------------------------------------------------
  #|
  # Determine which labels will be shown -
  select_labels <- function(my_data, label_max, sort_by_colume, label_colume, subset = TRUE) {
    subset_data <- data[subset, ]
    if (is.null(label_max) || is.na(label_max) || nrow(subset_data) <= label_max) {
      labels_shown <- subset_data[,  "tid_user"]
    } else {
      index_order <- rev(do.call(order, subset_data[ , sort_by_colume]))
      top_indexes <- index_order[1:label_max]
      labels_shown <- subset_data[top_indexes,  "tid_user"]
    }
    labels_shown <- labels_shown[!is.na(subset_data[labels_shown, label_colume])]  # Do not make grobs for NA
    return(my_data[ ,  "tid_user"] %in% labels_shown)
  }
  
  data$vl_is_shown <- select_labels(data, vertex_label_max,
                                    sort_by_colume = c("vls_plot", "vs_plot"), label_colume = "vl_user")
  data$el_is_shown <- select_labels(data, edge_label_max,
                                    sort_by_colume = c("els_plot", "es_plot"), label_colume = "el_user")
  data$tl_is_shown <- select_labels(data, tree_label_max, subset = data$is_root,
                                    sort_by_colume = c("tls_plot", "vs_plot"), label_colume = "tl_user")
  # Estimate plotted radius of vertex and tree labels
  tx_plot <- vapply(split(data$vx_plot, data$subgraph_root), FUN.VALUE = numeric(1),
                    function(x) mean(range(x)))
  data$tx_plot <- rep(tx_plot[unique(data$subgraph_root)], tree_vertex_counts)
  ty_plot <- vapply(split(data$vy_plot, data$subgraph_root), FUN.VALUE = numeric(1),
                    function(x) mean(range(x)))
  data$ty_plot <- rep(ty_plot[unique(data$subgraph_root)], tree_vertex_counts)
  data$tlx_plot <- data$tx_plot 
  tly_plot <- vapply(split(data$vy_plot, data$subgraph_root), FUN.VALUE = numeric(1), max)
  data$tly_plot <- rep(tly_plot[unique(data$subgraph_root)], tree_vertex_counts) + data$tls_plot * 1.1
  data$tlx_plot <- data$tx_plot 
  vl_radius_plot <- data$vls_plot * nchar(data$vl_user) / 2
  tl_radius_plot <- data$tls_plot * nchar(data$tl_user) / 2
  vl_data_shown <- data[data$vl_is_shown, ]
  tl_data_shown <- data[data$tl_is_shown, ]
  
  vlx_min <- (data$vx_plot - vl_radius_plot)[data$vl_is_shown]
  vlx_max <- (data$vx_plot + vl_radius_plot)[data$vl_is_shown]
  tlx_min <- (data$tx_plot - tl_radius_plot)[data$tl_is_shown]
  tlx_max <- (data$tlx_plot + tl_radius_plot)[data$tl_is_shown]
  vly_min <- (data$vy_plot - data$vls_plot)[data$vl_is_shown]
  vly_max <- (data$vy_plot + data$vls_plot)[data$vl_is_shown]
  tly_min <- (data$tly_plot - data$tls_plot)[data$tl_is_shown]
  tly_max <- (data$tly_plot + data$tls_plot)[data$tl_is_shown]
  
  x_points <- c(element_data$x, vlx_min, vlx_max, tlx_min, tlx_max)
  y_points <- c(element_data$y, vly_min, vly_max, tly_min, tly_max)
  
  
  margin_size_plot <- margin_size * square_side_length
  x_min <- min(x_points) - margin_size_plot[1]
  x_max <- max(x_points) + margin_size_plot[1]
  y_min <- min(y_points) - margin_size_plot[2]
  y_max <- max(y_points) + margin_size_plot[2]
  data$vls_plot <- scales::rescale(data$vls_plot, to = c(0, 1), from = c(0, square_side_length))
  data$vlx_plot <- scales::rescale(data$vx_plot, to = c(0, 1), from = c(x_min, x_max))
  data$vly_plot <- scales::rescale(data$vy_plot, to = c(0, 1), from = c(y_min, y_max))
  
  vertex_label_grobs <- plyr::dlply(data[data$vl_is_shown, ], "tid_user",
                                    function(row) resizingTextGrob(label = row['vl_user'],
                                                                   y = row['vly_plot'],
                                                                   x = row['vlx_plot'],
                                                                   gp = grid::gpar(text_prop = row['vls_plot'],
                                                                                   col = as.character(row['vlc_plot']))))    
  #|
  #| #### Make legend text grobs -------------------------------------------------------------------
  #|
  if (!is.null(legend_data)) {
    legend_label_data <- legend_data$labels
    legend_label_data$x_prop <- scales::rescale(legend_label_data$x, to = c(0, 1), from = c(x_min, x_max))
    legend_label_data$y_prop <- scales::rescale(legend_label_data$y, to = c(0, 1), from = c(y_min, y_max))
    legend_label_data$s_prop <- scales::rescale(legend_label_data$size, to = c(0, 1), from = c(0, square_side_length))
    legend_label_data$id <- 1:nrow(legend_label_data)
    # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    legend_label_grobs <- plyr::dlply(legend_label_data, "id",
                                      function(row) resizingTextGrob(label = row['label'],
                                                                     x = row['x_prop'],
                                                                     y = row['y_prop'],
                                                                     gp = grid::gpar(text_prop = row['s_prop'],
                                                                                     col = as.character(row['color'])),
                                                                     rot = row['rot'],
                                                                     just =  row[['just']]))
    
  } else {
    legend_label_grobs <- list()
  }
  
  #|
  #| #### Make edge text grobs -------------------------------------------------------------------
  #|
  data$el_user[is.na(data$pid_user)] <- ""
  # line label rotation  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$el_slope <- (data$vy_plot - data[data$pid_user, "vy_plot"]) / (data$vx_plot - data[data$pid_user, "vx_plot"])
  data$el_slope[is.na(data$el_slope)] <- 0
  data$el_rotation <- atan(data$el_slope)
  # line label coordinate  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  line_label_offset = 1
  justify <- data[data$pid_user, "vx_plot"] > data$vx_plot
  justify[is.na(justify)] <- TRUE
  justification <- lapply(1:nrow(data), function(i) if (justify[i]) c("left", "center") else c("right", "center"))
  names(justification) <- data$tid_user
  line_label_x_offset <- line_label_offset * data$vs_plot * cos(data$el_rotation)
  line_label_y_offset <- line_label_offset * data$vs_plot * sin(data$el_rotation)
  data$elx_plot <- data$vx_plot + ifelse(justify, 1, -1) * line_label_x_offset
  data$ely_plot <- data$vy_plot + ifelse(justify, 1, -1) * line_label_y_offset
  data$elx_plot <- scales::rescale(data$elx_plot, to = c(0, 1), from = c(x_min, x_max))
  data$ely_plot <- scales::rescale(data$ely_plot, to = c(0, 1), from = c(y_min, y_max))
  # line label text size - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$els_plot <-  scales::rescale(data$els_plot, to = c(0, 1), from = c(0, square_side_length))
  # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  edge_label_grobs <- plyr::dlply(data[data$el_is_shown, ], "tid_user",
                                  function(row) resizingTextGrob(label = row['el_user'],
                                                                 y = row['ely_plot'],
                                                                 x = row['elx_plot'],
                                                                 gp = grid::gpar(text_prop = row['els_plot'],
                                                                                 col = as.character(row['elc_plot'])),
                                                                 rot = row['el_rotation'] * 180 / pi,
                                                                 just = justification[[row[['tid_user']]]]))
  #|
  #| #### Make subplot title text grobs -----------------------------------------------------------
  #|
  title_data <- data[data$tl_is_shown, ]
  title_data$tlx_prop <- scales::rescale(title_data$tlx_plot, to = c(0, 1), from = c(x_min, x_max))
  title_data$tly_prop <- scales::rescale(title_data$tly_plot, to = c(0, 1), from = c(y_min, y_max))
  title_data$tls_prop <- scales::rescale(title_data$tls_plot, to = c(0, 1), from = c(0, square_side_length))
  # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tree_label_grobs <- plyr::dlply(title_data, "tid_user",
                                  function(row) resizingTextGrob(label = row['tl_user'],
                                                                 x = row['tlx_prop'],
                                                                 y = row['tly_prop'],
                                                                 gp = grid::gpar(text_prop = row['tls_prop'],
                                                                                 col = as.character(row['tlc_plot'])),
                                                                 rot = 0,
                                                                 just = "center"))
  #| ### Draw plot ================================================================================
  the_plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_polygon(data = element_data, ggplot2::aes(x = x, y = y, group = group),
                          fill = element_data$color) +
    ggplot2::guides(fill = "none") +
    ggplot2::coord_fixed(xlim = c(x_max, x_min), ylim = c(y_max, y_min)) +
    ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   panel.background = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text  = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), 
                   axis.line  = ggplot2::element_blank())
  for (a_grob in edge_label_grobs) {
    the_plot <- the_plot + ggplot2::annotation_custom(grob = a_grob)
  }    
  for (a_grob in vertex_label_grobs) {
    the_plot <- the_plot + ggplot2::annotation_custom(grob = a_grob)
  }    
  for (a_grob in tree_label_grobs) {
    the_plot <- the_plot + ggplot2::annotation_custom(grob = a_grob)
  }    
  for (a_grob in legend_label_grobs) {
    the_plot <- the_plot + ggplot2::annotation_custom(grob = a_grob)
  }    
  the_plot
}



#' Verify size range parameters
#' 
#' Verify size range parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
verify_size_range <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (length(value) != 2) {
      stop(paste0("Argument ", arg, " must be of length 2."))
    }
    if (all(!is.na(value)) && value[2] < value[1]) {
      stop(paste0("The min value of ", arg, " is greater than its max."))
    }
  }
}


#' Verify transformation function parameters
#' 
#' Verify transformation function parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
verify_trans <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (! is.function(value) && ! value %in% transform_data()) {
      stop(paste0("Argument '", arg,
                  "' must be a function or the name of a built-in transformation function."))
    }
  }
}

#' Verify size parameters
#' 
#' Verify size parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
verify_size <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (! is.numeric(value)) {
      stop(paste0("Argument '", arg, "' is not numeric."))
    }
  }
}


#' Verify color range parameters
#' 
#' Verify color range parameters
#' 
#' @param args (\code{character}) The names of arguments to verify.
verify_color_range <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (any(! grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", value) & ! value %in% colors())) {
      stop(paste0("Argument '", arg, "' must be hex color codes or a name returned by 'colors()'"))
    }
  }
}



#' Verify label count
#' 
#' Verify label count
#' 
#' @param args (\code{character}) The names of arguments to verify.
verify_label_count <- function(args) {
  for (arg in args) {
    value <- get(arg, pos = parent.frame())
    if (value %% 1 != 0 || value < 0) {
      stop(paste0("Argument '", arg, "' must be a positive integer."))
    }
  }
}


#' Check length of graph attributes
#' 
#' Length should divind evenly into the number of taxon/parent IDs
check_element_length <- function(args) {
  for (arg in args) {
    observed_length <- length(get(arg, pos = parent.frame()))
    correct_length <- length(get("taxon_id", pos = parent.frame()))
    if (observed_length < 1) {
      stop(paste0("Argument '", arg, "' is empty."))
    }
    if (correct_length %% observed_length != 0) {
      stop(paste0("Length of argument'", arg, "' must be a factor of the length of 'taxon_id'"))
    }
  }
}


#' Transformation functions
#' 
#' Functions used by plotting funtions to transform data.
#' Calling the function with no parameters returns available function names.
#' 
#' @param data (\code{numeric}) Data to transform
#' @param func (\code{character}) Name of transformation to apply.
#' @param inverse If TRUE, return inverse
transform_data <- function(data = NULL, func = NULL, inverse = FALSE) {
  sign <- function(x) {
    ifelse(x < 0, -1, 1)
  }
  
  if (inverse) {
    funcs <- list("radius" = function(x) {x},
                  "area" = function(x) {sign(x) * pi*x^2})
    
  } else {
    funcs <- list("radius" = function(x) {x},
                  "area" = function(x) {sign(x) * (abs(x)/pi)^(1/2)},
                  "log10 radius" = function(x) {log(x, base = 10)},
                  "log2 radius" = function(x) {log(x, base = 2)},
                  "ln radius" = function(x) {log(x)},
                  "log10 area" = function(x) {log((x/pi)^(1/2), base = 10)},
                  "log2 area" = function(x) {log((x/pi)^(1/2), base = 2)},
                  "ln area" =  function(x) {log((x/pi)^(1/2))})
    
  }
  
  if (is.null(data) & is.null(func)) {
    return(names(funcs))
  } else if (is.null(data)) {
    return(funcs[[func]])
  } else {
    return(vapply(X = data, FUN = funcs[[func]], FUN.VALUE = numeric(1)))
  }
}


#' Layout functions
#' 
#' Functions used to determine graph layout.
#' Calling the function with no parameters returns available function names.
#' 
#' @param name (\code{character} of length 1 OR NULL) name of algorithm. Leave \code{NULL} to 
#' see all options. 
#' @param graph (\code{igraph}) The graph to generate the layout for.
#' @param  intitial_coords (\code{matrix}) Initial vertex layout to bawse new layout off of. 
#' @param ... (other arguments) Passed to igraph layout function used.
#' 
#' @export
layout_functions <- function(name = NULL, graph = NULL, intitial_coords = NULL, ...) {
  return_names <- is.null(name) && is.null(graph) && is.null(intitial_coords)
  if (return_names) {
    graph <- igraph::make_ring(1) # Dummy graph so that the list can be defined, but only names used
  }
  defaults <- list("automatic" = list(),
                   "reingold-tilford" = list(circular = TRUE,
                                             mode = "out"),
                   "davidson-harel" = list(coords = intitial_coords,
                                           maxiter = 15,
                                           fineiter = max(15, log2(igraph::vcount(graph))),
                                           cool.fact = 0.75,
                                           weight.node.dist = 1,
                                           weight.border = 0,
                                           weight.edge.lengths = 1,
                                           weight.edge.crossings = 100,
                                           weight.node.edge.dist = 1),
                   "gem" = list(coords = intitial_coords,
                                maxiter = 40 * igraph::vcount(graph)^2,
                                temp.max = igraph::vcount(graph),
                                temp.min = 1/10,
                                temp.init = sqrt(igraph::vcount(graph))),
                   "graphopt" =list(start = intitial_coords,
                                    niter = 500,
                                    charge = 0.0005,
                                    mass = 30,
                                    spring.length = 0,
                                    spring.constant = 1,
                                    max.sa.movement = 5),
                   "mds" = list(),
                   "fruchterman-reingold" = list(coords = intitial_coords,
                                                 niter = 500,
                                                 start.temp = sqrt(igraph::vcount(graph)),
                                                 grid = "nogrid",
                                                 weights = NULL),
                   "kamada-kawai" = list(coords = intitial_coords,
                                         maxiter = 50 * igraph::vcount(graph),
                                         epsilon = 0,
                                         kkconst = igraph::vcount(graph),
                                         weights = NULL),
                   "large-graph" = list(maxiter = 150,
                                        maxdelta = igraph::vcount(graph),
                                        area = igraph::vcount(graph)^2,
                                        coolexp = 1.5,
                                        repulserad = igraph::vcount(graph)^2 * igraph::vcount(graph),
                                        cellsize = sqrt(sqrt(igraph::vcount(graph)^2)),
                                        root = 1),
                   "drl" = list(use.seed = ! is.null(intitial_coords),
                                seed = ifelse(is.null(intitial_coords), 
                                              matrix(runif(igraph::vcount(graph) * 2), ncol = 2),
                                              intitial_coords),
                                options = igraph::drl_defaults$default,
                                weights = igraph::E(graph)$weight,
                                fixed = NULL))
  funcs <- list("automatic" = igraph::nicely,
                "reingold-tilford" = igraph::as_tree,
                "davidson-harel" = igraph::with_dh,
                "gem" = igraph::with_gem,
                "graphopt" = igraph::with_graphopt,
                "mds" = igraph::with_mds(),
                "fruchterman-reingold" = igraph::with_fr,
                "kamada-kawai" = igraph::with_kk,
                "large-graph" = igraph::with_lgl,
                "drl" = igraph::with_drl)
  if (return_names) {
    return(names(funcs)) 
  } else {
    arguments <- modifyList(defaults[[name]], list(...))
    return(do.call(funcs[[name]], arguments))
  }
}


#' The defualt quantative color palette
#' 
#' The default color palette for quantative data
quantative_palette <- function() {
  return(c("grey", "#018571", "#80cdc1", "#dfc27d", "#a6611a"))
}


#' The defualt qualitative color palette
#' 
#' The default color palette for qualitative data
qualitative_palette <- function() {
  return(c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(9, "Pastel1")))
}




#' Covert numbers to colors
#' 
#' Convert numbers to colors.
#' If colors are already supplied, return the input
#' 
#' @param values (\code{numeric}) The numbers to represent as colors
#' @param color_series (\code{character}) Hex values or a character in \code{colors}
#' @param no_color_in_palette (\code{numeric} of length 1) The number of distinct colors to use.
#' 
#' 
#' @return \code{character} Hex color codes. 
apply_color_scale <- function(values, color_series, no_color_in_palette = 1000) {
  if (is.numeric(values)) { ## Not factors, characters, or hex codes
    palette <- colorRampPalette(color_series)(no_color_in_palette)
    color_index <- as.integer(scales::rescale(values, to = c(1, no_color_in_palette)))
    return(palette[color_index])
  } else {
    return(values)
  }
}


#===================================================================================================
#' Makes coordinates for a regualr polygon
#' 
#' Generates an n x 2 matrix containing x and y coordinates between 1 and 0 for the points of a 
#' regular polygon. 
#' 
#' Inspired by (i.e. stolen from) https://gist.github.com/baptiste/2224724, which was
#' itself inspired from a post by William Dunlap on r-help (10/09/09)
#' 
#' @param n (\code{numeric} of length 1) The number of vertices in the polygon.
#' @param x (\code{numeric} of length 1) x coordinate of center
#' @param y (\code{numeric} of length 1) y coordinate of center
#' @param radius (\code{numeric} of length 1) The diameter of the circle.
#' @param angle (\code{numeric} of length 1) Angle to rotate points around the center of the circle.
#' 
#' @examples
#' ggplot(data = polygon_coords(n = 4:13, x = rnorm(10), y = rnorm(10), radius = .5)) + 
#'   geom_polygon(aes(x = x, y = y, fill = group))
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
#' @examples
#' ggplot(data = line_coords(x1 = 1, y1 = 1, x2 = 2, y2 = 2, width = .1)) + 
#'   geom_polygon(aes(x = x, y = y, fill = group))
#' ggplot(data = line_coords(x1 = rnorm(10), y1 = rnorm(10), x2 = rnorm(10),
#'                           y2 = rnorm(10), width = rnorm(10)/5)) + 
#'   geom_polygon(aes(x = x, y = y, fill = group))
#'   
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


#===================================================================================================
#' Get all distances beween points
#' 
#' Returns the distances bewteen every possible combination of two points.
#' 
#' @param x (\code{numeric} of length 1) x coordinate
#' @param y (\code{numeric} of length 1) y coordinate
#' 
#' @return A \code{data.frame}
#' 
#' @examples
#' molten_dist(x = 1:5, y = 1:5)
#' 
molten_dist <- function(x, y) {
  data <- as.matrix(dist(cbind(x, y)))
  data[!lower.tri(data)] <- NA
  data <- reshape2::melt(data)
  names(data) <- c("index_1", "index_2", "distance")
  data[!is.na(data$distance), ]    
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
#' @examples
#' inter_circle_gap(x = 1:5, y = 1:5, r = 1:5)
#' 
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
#' Finds optimal max and min value useing an optimality criterion.
#' 
#' @param max_range (\code{numeric} of length 2) The min and max boundries to the search space for
#' the optimal maximum value.
#' @param min_range (\code{numeric} of length 2) The min and max boundries to the search space for
#' the optimal minimum value.
#' @param resolution (\code{numeric} of length 2) The number of increments in each dimension.
#' @param opt_crit (\code{function}) A function that takes two arguments, the max and min, and
#' returns the optimality statistic.
#' @param choose_best (\code{function}) A function that takes a list of \code{opt_crit} outputs
#' and returns the index of the best one.
#' 
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
    signif(n, digits = 3)
  }
  
  scale_undo_trans <- function(points, my_range, my_trans) {
    trans_points <- vapply(points, my_trans, FUN.VALUE = numeric(1))
    format_label(scales::rescale(trans_points, to = range(my_range)))
  }
  
  label_color = "#000000"
  label_size = max(width_range) / 5
  
  
  make_size_label_data <- function() {
    size_label_data <- data.frame(stringsAsFactors = FALSE, 
                                  label = scale_undo_trans(label_point_y, width_stat_range, width_stat_trans),
                                  x = max(width_range) * 1.1,
                                  y = rev(label_point_y),
                                  color = label_color,
                                  size = label_size,
                                  rot = 0,
                                  just = "left")
  }
  
  make_color_label_data <- function() {
    
    color_label_data <- data.frame(stringsAsFactors = FALSE, 
                                   label = scale_undo_trans(label_point_y, color_stat_range, color_stat_trans),
                                   x = 0 - max(width_range) * 0.1,
                                   y = rev(label_point_y),
                                   color = label_color,
                                   size = label_size,
                                   rot = 0,
                                   just = "right")
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
scale_bar_coords <- function(x1, x2, y1, y2, color, group) {
  data.frame(group = group,
             x = c(x1, x2, 0, 0),
             y = c(y1, y2, y2, y1),
             color = color)
  
}


