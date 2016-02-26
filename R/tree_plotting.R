#| ## Tree plotting
#|
#|
#|
#|
#|
#|
#| 
#|
#| 
#| 
#|  
#|
#| 
#|  
#|
#|
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
#' @param tree_size See details on size.
#' The value of the root of each graph will be used.
#' This scales the space used to display graphs, but does not effect vertex/edge size.
#' Default: Not used. 
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
#' @param vertex_color_trans See details on transformations. 
#' Default: \code{"area"}.
#' @param edge_color_trans See details on transformations.
#' Default: same as vertex color transformation.
#' @param tree_color_trans See details on transformations.
#' Default: \code{"area"}.
#' 
#' @param vertex_label_size_trans See details on transformations. 
#' Default: same as \code{vertex_size_trans}.
#' @param edge_label_size_trans See details on transformations. 
#' Default: same as \code{edge_size_trans}.
#' @param tree_label_size_trans See details on transformations.
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
#' @param tree_size_range See details on ranges.
#' Default: Not set.
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
new_plot_taxonomy <- function(taxon_id, parent_id, 
                              vertex_size = 1,
                              edge_size = vertex_size,
                              tree_size = 1,

                              vertex_color = "#999999",
                              edge_color = vertex_color,
                              tree_color = NULL,
                              
                              vertex_label = NA,
                              edge_label = NA,
                              tree_label = NA,

                              vertex_label_size = vertex_size,
                              edge_label_size = edge_size,
                              tree_label_size = NULL, 
                              
                              vertex_label_color = "#000000",
                              edge_label_color = "#000000",
                              tree_label_color = "#000000",
                              
                              
                              vertex_size_range = c(NA, NA),
                              vertex_size_trans = "area",
                              vertex_color_range = quantative_palette(),
                              vertex_color_trans = vertex_size_trans,
                              edge_size_range = c(NA, NA),
                              edge_size_trans = vertex_size_trans,
                              edge_color_range = vertex_color_range,
                              edge_color_trans = vertex_color_trans,
                              vertex_label_size_range = c(NA, NA),
                              vertex_label_size_trans = vertex_size_trans,
                              vertex_label_color_range = quantative_palette(),
                              vertex_label_color_trans = "area",
                              edge_label_size_range = c(NA, NA),
                              edge_label_size_trans = edge_size_trans,
                              edge_label_color_range = quantative_palette(),
                              edge_label_color_trans = "area",
                              tree_label_size_range = c(NA, NA),
                              tree_label_color_range = quantative_palette(),
                              tree_label_color_trans = "area",
                              vertex_label_max = 20,
                              edge_label_max = 20,
                              tree_label_max = 20,
                              overlap_avoidance = 1,
                              margin_size = c(0, 0),
                              # aspect_ratio = NULL,
                              layout = "reingold-tilford",
                              initial_layout = "fruchterman-reingold",
                              ...) {
  #| ### Verify arguments =========================================================================
  if (length(taxon_id) != length(parent_id)) {
    stop("'taxon_id' and 'parent_id' must be of equal length.")
  }
  if (length(taxon_id) == 0) {
    stop("'taxon_id' and 'parent_id' are empty.")
  }
  check_element_length(c("vertex_size", "edge_size", "vertex_label_size", "edge_label_size",
                         "vertex_color", "edge_color", "vertex_label_color", "edge_label_color",
                         "vertex_label", "edge_label", "tree_label", "tree_label_color"))
  verify_size(c("vertex_size", "edge_size", "vertex_label_size", "edge_label_size"))
  verify_size_range(c("vertex_size_range", "edge_size_range",
                      "vertex_label_size_range", "edge_label_size_range"))
  verify_trans(c("vertex_size_trans", "vertex_color_trans", 
                 "edge_size_trans", "edge_color_trans",
                 "vertex_label_size_trans", "vertex_label_color_trans",
                 "edge_label_size_trans", "edge_label_color_trans", "tree_label_color_trans"))
  verify_color_range(c("vertex_color_range", "edge_color_range",
                       "vertex_label_color_range", "edge_label_color_range",
                       "tree_label_color_range"))
  verify_label_count(c("vertex_label_max", "edge_label_max", "tree_label_max"))
  if (length(overlap_avoidance) == 0 || ! is.numeric(overlap_avoidance)) {
    stop("Argument 'overlap_avoidance' must be a numeric of length 1.")
  }
  if (length(margin_size) != 2 || ! is.numeric(margin_size)) {
    stop("Argument 'margin_size' must be a numeric of length 2.")
  }
  #   if (! is.null(aspect_ratio) && ! is.numeric(aspect_ratio)) {
  #     stop("Argument 'aspect_ratio' must be a numeric of length 1.")
  #   }
  if (! layout %in% layout_functions()) {
    stop("Argument 'layout' must be an output of layout_functions().")
  }
  if (! initial_layout %in% layout_functions()) {
    stop("Argument 'initial_layout' must be an output of layout_functions().")
  }
  #| ### Standardize source data ==================================================================
  data <- data.frame(stringsAsFactors = FALSE,
                     tid = as.character(taxon_id),
                     pid = as.character(parent_id),
                     vs = as.numeric(vertex_size),
                     vc = vertex_color,
                     es = as.numeric(edge_size),
                     ec = edge_color,
                     vls = as.numeric(vertex_label_size),
                     vlc = vertex_label_color,
                     els = as.numeric(edge_label_size),
                     elc = edge_label_color,
                     glc = tree_label_color,
                     vl = as.character(vertex_label),
                     el = as.character(edge_label),
                     gl = as.character(tree_label))
  row.names(data) <- data$tid
  
  
  #| #### Apply statistic transformations ---------------------------------------------------------
  #|
  data$subgraph_prop <- vapply(data$subgraph_root, function(x) sum(data$subgraph_root == x) / nrow(data), numeric(1))
  data$gls <- data$subgraph_prop
  trans_key <- c(vs = vertex_size_trans, vc = vertex_color_trans,
                 vls = vertex_label_size_trans, vlc = vertex_label_color_trans,
                 es = edge_size_trans, ec = edge_color_trans,
                 els = edge_label_size_trans, elc = edge_label_color_trans,
                 glc = graph_label_color_trans, gls = "area")
  new_names <- paste0(names(trans_key), "_t")
  apply_trans <- function(col_name) {
    if (is.numeric(data[ , col_name])) { 
      transform_data(data[ , col_name], trans_key[col_name]) # if numbers are supplied
    } else {
      data[ , col_name] # if colors are defined explicitly, then no transformation is done
    }
  }
  data[, new_names] <- lapply(names(trans_key), apply_trans)
  
  
  #| ### Make layout ==============================================================================
  #| The layout is used to generate a list of coordinates to places graph verticies
  #| First the edge list consituted by the `taxon_id` and `parent_id` columns is used to construct 
  #| an `igraph` graph object and then the layout is generated for that object. 
  #|
  #| #### Make a graph for each root in the graph -------------------------------------------------
  get_sub_graphs <- function(taxa) {
    if (length(taxa) == 1) {
      # Make a graph with only a single vertex
      return(igraph::graph.adjacency(matrix(c(taxa), ncol = 1)))
    } else {
      # Make edge list from taxon_id and parent_id
      edgelist <- as.matrix(data[taxa, c("pid", "tid")])
      # Remove edges to taxa that dont exist in this subset of the dataset
      edgelist <- edgelist[! is.na(edgelist[, "pid"]), ]
      return(igraph::graph_from_edgelist(edgelist))
    }
  }
  data$pid[!(data$pid %in% data$tid)] <- NA # Needed by split_by_level
  sub_graph_taxa <- split_by_level(data$tid, data$pid, level =  1)
  sub_graphs <- lapply(sub_graph_taxa, get_sub_graphs)
  #|
  #| #### Generate a layout for each graph --------------------------------------------------------
  #|
  get_sub_layouts <- function(graph, backup_layout = 'fruchterman-reingold') {
    # Calculate an initial layout if specified
    if (! is.null(initial_layout) && layout != initial_layout) {
      intitial_coords <- igraph::layout_(graph, layout_functions(initial_layout, graph = graph))
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
      warning(paste0("Could not apply layout ", layout,
                     " to subgraph with. Using 'fruchterman-reingold' instead."))
    }
    return(coords)
  }
  sub_coords <- lapply(sub_graphs, get_sub_layouts)
  #|
  #| #### Merge layout coordinates into an overall graph ------------------------------------------
  #|
  coords <- merge_coords(sub_graphs, sub_coords)
  graph <- disjoint_union(sub_graphs)
  row.names(coords) <- names(V(graph))
  data$vx <- coords[data$tid, 1]
  data$vy <- coords[data$tid, 2]
  data$subgraph_root <- rep(names(sub_coords), vapply(sub_coords, nrow, numeric(1)))
  
  
  
  #| ### Core plot data ===========================================================================
  #|
  
  #|
  #| #### Optimize vertex size --------------------------------------------------------------------
  #|
  
  all_pairwise <- molten_dist(x = data$vx, y = data$vy)
  x_diff <- max(data$vx) - min(data$vx)
  y_diff <- max(data$vy) - min(data$vy)
  square_side_length <- sqrt(x_diff * y_diff)
  if (is.na(vertex_size_range[1])) {
    min_range <- c(square_side_length / 1000, min(all_pairwise$distance))
  } else {
    min_range <- rep(vertex_size_range[1], 2) * square_side_length
  }
  if (is.na(vertex_size_range[2])) {
    max_range <- c(min_range[1], square_side_length / 4)
  } else {
    max_range <- c(vertex_size_range[2], 2) * square_side_length
  }
  
  get_search_space <- function(min_range, max_range, breaks_per_dim = 10) {
    min_breaks <- seq(from = min_range[1], to = min_range[2], length.out = breaks_per_dim)
    max_breaks <- seq(from = max_range[1], to = max_range[2], length.out = breaks_per_dim)
    space <- data.frame(min = rep(min_breaks, each = length(max_breaks)),
                        max = rep(max_breaks, length(min_breaks)))
    space[space$min <= space$max, ]
  }
  
  
  search_space <- get_search_space(min_range, max_range, breaks_per_dim = 10)
  search_space$range_size <- search_space$max - search_space$min
  
  find_overlap <- function(a_min, a_max, distance) {
    scaled_vs <- scales::rescale(data$vs_t, to = c(a_min, a_max))
    names(scaled_vs) <- data$tid
    gap <- distance$distance - scaled_vs[distance$index_1] - scaled_vs[distance$index_2]
    gap <- ifelse(gap < 0, abs(gap), 0)
    sum(gap)
  } 
  
  search_space$overlap <- apply(search_space, MARGIN = 1,
                                function(x) find_overlap(x["min"], x["max"], all_pairwise))
  
  optimality_stat <- function(overlap, range_size, minimum) {
    overlap_weight <- 0.05
    minimum_weight <- 2
    (1 + range_size + minimum * minimum_weight) / (1 + overlap * overlap_avoidance * overlap_weight)
  }
  
  search_space$opt_stat <- apply(search_space, MARGIN = 1,
                                 function(x) optimality_stat(x["overlap"], x["range_size"], x["min"]))
  vertex_size_range_g <- unlist(search_space[which.max(search_space$opt_stat), c("min", "max")])
  data$vs_g <- scales::rescale(data$vs_t, to = vertex_size_range_g)
  # ggplot(search_space) + geom_point(aes(x = max, y = min, size = opt_stat))
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
  
  edge_size_range_g <- infer_size_range(edge_size_range, vertex_size_range_g, defualt_scale = 0.5)
  data$es_g <- scales::rescale(data$es_t, to = edge_size_range_g)
  #|
  #| #### Infer label size ranges -----------------------------------------------------------------
  #|
  vertex_label_size_range_g <- infer_size_range(vertex_label_size_range, vertex_size_range_g, 
                                                defualt_scale = 0.5)
  edge_label_size_range_g <- infer_size_range(edge_label_size_range, edge_size_range_g, 
                                              defualt_scale = 0.7)
  tree_label_size_range_g <- infer_size_range(tree_label_size_range, range(data$subgraph_prop), 
                                              defualt_scale = 0.1)
  data$vls_g <- scales::rescale(data$vls_t, to = vertex_label_size_range_g)
  data$els_g <- scales::rescale(data$els_t, to = edge_label_size_range_g)
  data$gls_g <- scales::rescale(data$gls_t, to = tree_label_size_range_g)
  #|
  #| #### Assign color scales ---------------------------------------------------------------------
  #|
  apply_color_scale <- function(values, color_series, no_color_in_palette = 1000) {
    if (is.numeric(values)) { ## Not factors, characters, or hex codes
      palette <- colorRampPalette(color_series)(no_color_in_palette)
      color_index <- as.integer(scales::rescale(values, to = c(1, no_color_in_palette)))
      return(palette[color_index])
    } else {
      return(values)
    }
  }
  
  color_colume_key <- list("ec_t" = edge_color_range, "vc_t" = vertex_color_range,
                           "glc_t" = tree_label_color_range,
                           "elc_t" = edge_label_color_range, "vlc_t" = vertex_label_color_range)
  new_names <- gsub(pattern = "_t$", x = names(color_colume_key), replacement = "_g")
  data[, new_names] <- lapply(names(color_colume_key),
                              function(x) apply_color_scale(data[ , x], color_colume_key[[x]]))
  
  #| ### Secondary plot data ======================================================================
  #|
  #| #### Calculate coordinants of graph elements -------------------------------------------------
  #| The vertexes and edges must be specified by a dataframe of coordinates, with a colume 
  #| grouping the coordinates of each shape.
  #| These shapes must be added to the graph in a specific order.
  #| A list of vertexes is sorted by first vertex depth in the heirarchy and then by vertex size.
  taxon_elements <- function(taxon_id) {
    circle_resolution <- 50
    edge_data <- line_coords(x1 = data[taxon_id, 'vx'],
                             y1 = data[taxon_id, 'vy'],
                             x2 = data[data[taxon_id, 'pid'], "vx"],
                             y2 = data[data[taxon_id, 'pid'], "vy"],
                             width = data[taxon_id, 'es_g'] * 2)
    edge_data$group <- paste0(taxon_id, "_edge")
    edge_data$color <- rep(data[taxon_id, 'ec_g'], each = 4)
    vertex_data <- polygon_coords(n = circle_resolution,
                                  x = data[taxon_id, 'vx'],
                                  y = data[taxon_id, 'vy'],
                                  radius = data[taxon_id, 'vs_g'])
    vertex_data$group <- paste0(taxon_id, "_vertex")
    vertex_data$color <- rep(data[taxon_id, 'vc_g'], each = circle_resolution + 1)
    output <- rbind(edge_data, vertex_data)
    output$tid <- taxon_id
    return(output[complete.cases(output),])
  }
  data$level = edge_list_depth(data$tid, data$pid)
  element_order <- data$tid[order(data$level, 1 / data$vs_g, decreasing = TRUE)]
  element_data <- do.call(rbind, lapply(element_order, taxon_elements))
  element_data$group <- factor(element_data$group, levels = unique(element_data$group))
  #|
  #| #### Make vertex text grobs ------------------------------------------------------------------
  #|
<<<<<<< HEAD
  margin_size <- margin_size * square_side_length
  x_min <-  min(vertex_data$x) - margin_size[1]
  x_max <- max(vertex_data$x) + margin_size[1]
  y_min <- min(vertex_data$y) - margin_size[2]
  y_max <- max(vertex_data$y) + margin_size[2]
  data$vls_g <- scales::rescale(data$vls_t, to = c(0, 1), from = c(0, square_side_length))
=======
  x_min <- min(element_data$x) - margin_size[1] * square_side_length
  x_max <- max(element_data$x) + margin_size[1] * square_side_length
  y_min <- min(element_data$y) - margin_size[2] * square_side_length
  y_max <- max(element_data$y) + margin_size[2] * square_side_length
  data$vls_g <- scales::rescale(data$vls_g, to = c(0, 1), from = c(0, square_side_length))
>>>>>>> 636f95e9952b9b0baea31f57245c8ae759fac105
  data$vlx_g <- scales::rescale(data$vx, to = c(0, 1), from = c(x_min, x_max))
  data$vly_g <- scales::rescale(data$vy, to = c(0, 1), from = c(y_min, y_max))
  
  select_labels <- function(my_data, label_max, sort_by_colume, label_colume) {
    if (is.null(label_max) || is.na(label_max) || nrow(my_data) <= label_max) {
      labels_shown <- my_data[ ,  "tid"]
    } else {
      top_indexes <- order(my_data[ , sort_by_colume], decreasing = TRUE)[1:label_max]
      labels_shown <- my_data[top_indexes,  "tid"]
    }
    return( labels_shown[!is.na(my_data[labels_shown, label_colume])] ) # Do not make grobs for NA
  }
  
  vertex_label_shown <- select_labels(data, vertex_label_max,
                                      sort_by_colume = "vls_g", label_colume = "vl")
  vertex_label_grobs <- plyr::dlply(data[vertex_label_shown, ], "tid",
                               function(row) resizingTextGrob(label = row['vl'],
                                                            y = row['vly_g'],
                                                            x = row['vlx_g'],
                                                            gp = grid::gpar(text_prop = row['vls_g'])))    
  #|
  #| #### Make edge text grobs -------------------------------------------------------------------
  #|
  data$el[is.na(data$pid)] <- ""
  # line label rotation  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$el_slope <- (data$vy - data[data$pid, "vy"]) / (data$vx - data[data$pid, "vx"])
  data$el_slope[is.na(data$el_slope)] <- 0
  data$el_rotation <- atan(data$el_slope)
  # line label coordinate  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  line_label_offset = 1
  justify <- data[data$pid, "vx"] > data$vx
  justify[is.na(justify)] <- TRUE
  justification <- lapply(1:nrow(data), function(i) if (justify[i]) c("left", "center") else c("right", "center"))
  names(justification) <- data$tid
  line_label_x_offset <- line_label_offset * data$vs_g * cos(data$el_rotation)
  line_label_y_offset <- line_label_offset * data$vs_g * sin(data$el_rotation)
  data$elx_g <- data$vx + ifelse(justify, 1, -1) * line_label_x_offset
  data$ely_g <- data$vy + ifelse(justify, 1, -1) * line_label_y_offset
  data$elx_g <- scales::rescale(data$elx_g, to = c(0, 1), from = c(x_min, x_max))
  data$ely_g <- scales::rescale(data$ely_g, to = c(0, 1), from = c(y_min, y_max))
  # line label text size - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data$els_g <-  scales::rescale(data$els_g, to = c(0, 1), from = c(0, square_side_length))
  # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  edge_label_shown <- select_labels(data, edge_label_max, 
                                    sort_by_colume = "els_g", label_colume = "el")
  edge_label_grobs <- plyr::dlply(data[edge_label_shown, ], "tid",
                                    function(row) resizingTextGrob(label = row['el'],
                                                                   y = row['ely_g'],
                                                                   x = row['elx_g'],
                                                                   gp = grid::gpar(text_prop = row['els_g']),
                                                                   rot = row['el_rotation'] * 180 / pi,
                                                                   just = justification[[row[['tid']]]]))
  #|
  #| #### Make subplot title text grobs -----------------------------------------------------------
  #|
  root_indexes <- unique(data$subgraph_root)
  title_data <- data[root_indexes, ]
  title_data$glx <- vapply(title_data$subgraph_root, FUN.VALUE = numeric(1),
                           function(x) mean(range(data[data$subgraph_root == x, "vx"])))
  title_data$gly <- vapply(title_data$subgraph_root, FUN.VALUE = numeric(1),
                           function(x) max(data[data$subgraph_root == x, "vy"]))
  title_data$glx <- scales::rescale(title_data$glx, to = c(0, 1), from = c(x_min, x_max))
  title_data$gly <- scales::rescale(title_data$gly, to = c(0, 1), from = c(y_min, y_max)) + title_data$gls_g * 1.1
  # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tree_label_shown <- select_labels(title_data, tree_label_max, 
                                     sort_by_colume = "gls_g", label_colume = "gl")
  tree_label_grobs <- plyr::dlply(title_data[tree_label_shown, ], "tid",
                                  function(row) resizingTextGrob(label = row['gl'],
                                                                 x = row['glx'],
                                                                 y = row['gly'],
                                                                 gp = grid::gpar(text_prop = row['gls_g']),
                                                                 rot = 0,
                                                                 just = "center"))
  #| ### Draw plot ================================================================================
  the_plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_polygon(data = element_data, ggplot2::aes(x = x, y = y, group = group),
                          fill = element_data$color) +
    ggplot2::guides(fill = "none") +
    ggplot2::coord_fixed(xlim = c(x_max, x_min), ylim = c(y_max, y_min)) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
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
    if (! all(grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", value) | value %in% colors())) {
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

# 
# #' Standardize color input 
# #' 
# #' This forces user input to be the correct type.
# standardize_color <- function(data) {
#   if (all(grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", data) | data %in% colors())) {
#     return(as.character(data))
#   } else 
#   
#   if (is.factor(data)) {
#     if (length)
#   }
# }
# 

#' Transformation functions
#' 
#' Functions used by plotting funtions to transform data.
#' Calling the function with no parameters returns available function names.
#' 
#' @param data (\code{numeric}) Data to transform
#' @param func (\code{character}) Name of transformation to apply.
transform_data <- function(data = NULL, func = NULL) {
  funcs <- list("radius" = function(x) {x},
                "area" = function(x) {(x/pi)^(1/2)},
                "log10 radius" = function(x) {log(x, base = 10)},
                "log2 radius" = function(x) {log(x, base = 2)},
                "ln radius" = function(x) {log(x)},
                "log10 area" = function(x) {log((x/pi)^(1/2), base = 10)},
                "log2 area" = function(x) {log((x/pi)^(1/2), base = 2)},
                "ln area" =  function(x) {log((x/pi)^(1/2))})
  if (is.null(data) & is.null(func)) {
    return(names(funcs))
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
    graph <- make_ring(1) # Dummy graph so that the list can be defined, but only names used
  }
  defaults <- list("automatic" = list(),
                   "reingold-tilford" = list(circular = TRUE,
                                             mode = "out"),
                   "davidson-harel" = list(coords = intitial_coords,
                                           maxiter = 15,
                                           fineiter = max(15, log2(vcount(graph))),
                                           cool.fact = 0.75,
                                           weight.node.dist = 1,
                                           weight.border = 0,
                                           weight.edge.lengths = 1,
                                           weight.edge.crossings = 100,
                                           weight.node.edge.dist = 1),
                   "gem" = list(coords = intitial_coords,
                                maxiter = 40 * vcount(graph)^2,
                                temp.max = vcount(graph),
                                temp.min = 1/10,
                                temp.init = sqrt(vcount(graph))),
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
                                                 start.temp = sqrt(vcount(graph)),
                                                 grid = "nogrid",
                                                 weights = NULL),
                   "kamada-kawai" = list(coords = intitial_coords,
                                         maxiter = 50 * vcount(graph),
                                         epsilon = 0,
                                         kkconst = vcount(graph),
                                         weights = NULL),
                   "large-graph" = list(maxiter = 150,
                                        maxdelta = vcount(graph),
                                        area = vcount(graph)^2,
                                        coolexp = 1.5,
                                        repulserad = vcount(graph)^2 * vcount(graph),
                                        cellsize = sqrt(sqrt(vcount(graph)^2)),
                                        root = 1),
                   "drl" = list(use.seed = ! is.null(intitial_coords),
                                seed = ifelse(is.null(intitial_coords), 
                                              matrix(runif(vcount(graph) * 2), ncol = 2),
                                              intitial_coords),
                                options = drl_defaults$default,
                                weights = E(graph)$weight,
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