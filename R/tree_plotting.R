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
#' Plot items on taxonomy
#' 
#' Plots the distribution of values associated with an taxonomic classification.
#' Colors and sizes of vertexes and edges can be mapped to variables.
#' Uses \code{igraph} to make layout and \code{ggplot2} to make plots.
#' 
#' @param taxon_id (\code{character}) The unique ids of the taxon for each row.
#' @param parent_id (\code{character}) The unique id of supertaxon \code{taxon_id} is a part of.
#' @param vertex_size (\code{numeric}) The value to base vertex size on. Default: inverse of depth 
#' of taxa in hierarchy. 
#' @param vertex_size_range (\code{numeric} of length 2) The minimum and maximum size of vertexes.
#' If either value is \code{NA}, the missing value(s) will be optimized to maximize range and
#' minimizing overlaps. Defualt: \code{c(NA, NA)}.
#' @param vertex_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_size}. Alternativly one of the following \code{character} values. Default:
#' \code{"area"}.
#' \describe{
#'   \item{"radius"}{Proprotional to radius/diameter of vertex}
#'   \item{"area"}{circular area; better perceptual accuracy than \code{"radius"}}
#'   \item{"log10 radius"}{Log base 10 of radius}
#'   \item{"log2 radius"}{Log base 2 of radius}
#'   \item{"ln radius"}{Log base e of radius}
#'   \item{"log10 area"}{Log base 10 of area}
#'   \item{"log2 area"}{Log base 2 of area}
#'   \item{"ln area"}{Log base e of area}
#' }
#' @param vertex_color (\code{numeric} OR \code{character}) The value to base vertex color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: grey.
#' @param vertex_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to \code{vertex_color}. Default: Color-blind friendly palette. 
#' @param vertex_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_size}. Alternativly, one of the \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: \code{"area"}.
#' @param edge_size (\code{numeric)} The value to base edge size on. Default: relative to vertex
#' size. 
#' @param edge_size_range (\code{numeric} of length 2) The minimum and maximum size of edges.
#' If either value is \code{NA}, the missing value(s) will be set relative to vertex size. 
#' Default: relative to vertex size. 
#' @param edge_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as \code{vertex_size_trans}. 
#' @param edge_color (\code{numeric} OR \code{character}) The value to base edge color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: same as vertex color.
#' @param edge_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to low-high statistics supplied to \code{edge_color}. Default: same as vertex color.
#' @param edge_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as vertex color transformation.
#' @param vertex_label (\code{character}) The values of labels over vertcies. Use \code{NA} to exclude labels.
#' Default: no labels.
#' @param vertex_label_size (\code{numeric)} The value to base veterx label size on. Default: 
#' relative to veterx size.
#' @param vertex_label_size_range (\code{numeric} of length 2) The minimum and maximum size of labels.
#' If either value is \code{NA}, the missing value(s) will be optimized relative to vertex size. 
#' Default: relative to vertex size. 
#' @param vertex_label_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_label_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as \code{vertex_size_trans}.
#' @param vertex_label_color (\code{numeric} OR \code{character}) The value to base label color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: black.
#' @param vertex_label_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to low-high statistics supplied to \code{vertex_label_color}. Default: Color-blind friendly palette. 
#' @param vertex_label_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{vertex_label_size}. Alternativly one of the \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: \code{"area"}.
#' @param vertex_label_max (\code{numeric}) The maximum number of vertex labels. Default: 20.
#' @param edge_label (\code{character}) The values of labels over edges. Default: no labels.
#' @param edge_label_size (\code{numeric)} The value to base edge label size on. Default: relative to
#' edge size.
#' @param edge_label_size_range (\code{numeric} of length 2) The minimum and maximum size of labels.
#' If either value is \code{NA}, the missing value(s) will be relative to edge size. Default: relative to
#' edge size.
#' @param edge_label_size_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_label_size}. Alternativly one of the \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: same as \code{edge_size_trans}.
#' @param edge_label_color (\code{numeric} OR \code{character}) The value to base label color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: black.
#' @param edge_label_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to low-high statistics supplied to \code{edge_label_color}. Default: Color-blind friendly palette.
#' @param edge_label_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{edge_label_size}. Alternativly one of the following \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: \code{"area"}.
#' @param edge_label_max (\code{numeric}) The maximum number of edge labels. Default: 20.
#' @param graph_label (\code{character}) The label to display above each graph. The value of the root
#' taxon of each graph will be used. Default: None.
#' @param graph_label_size_range (\code{numeric} of length 2) The minimum and maximum size of labels.
#' If either value is \code{NA}, the missing value(s) will be relative to proportion of total vertexes
#' in graph Default: relative to proportion of total vertexes in graph.
#' @param graph_label_color (\code{numeric} OR \code{character}) The value to base graph label color on. 
#' If a numeric vector is given, it is used to  construct a color scale. Hex values or color 
#' names can be used (e.g. \code{#000000} or \code{"black"}). Default: black.
#' @param  graph_label_color_range (valid argument of \code{col2rgb}) A series of colors corresponding 
#' to low-high statistics supplied to \code{graph_label_color}. Default: Color-blind friendly palette. 
#' @param graph_label_color_trans (\code{function(value)} OR \code{character}) A function to transform the
#' value of \code{graph_label_color}. Alternativly one of the \code{character} values displayed
#' under the \code{vertex_size_trans}. Default: \code{"area"}.
#' @param graph_label_max (\code{numeric}) The maximum number of graph labels. Default: 20.
#' @param overlap_bias (\code{numeric}) The relative importance of avoiding overlaps vs maximizing 
#' size range. Default: \code{1}.
#' @param margin_size (\code{numeric} of length 2) The horizontal and vertical margins.
#' Default: \code{0, 0}.
# #' @param aspect_ratio (\code{numeric}) The height / width of the plot. Default: Whatever the layout
# #' function produces.
#' @param layout (\code{character} of length 1) The layout function to use. 
#' Type \code{\link{layout_functions}()} for available layout names.
#' @param initial_layout (\code{character} of length 1) Optional starting layout to use to initialize
#' the final layout function.
#' Type \code{\link{layout_functions}()} for available layout names.
#'  
#' @export
new_plot_taxonomy <- function(taxon_id, parent_id, 
                              vertex_size = 1,
                              vertex_size_range = c(NA, NA),
                              vertex_size_trans = "area",
                              vertex_color = "#999999",
                              vertex_color_range = quantative_palette(),
                              vertex_color_trans = vertex_size_trans,
                              edge_size = vertex_size,
                              edge_size_range = c(NA, NA),
                              edge_size_trans = vertex_size_trans,
                              edge_color = vertex_color,
                              edge_color_range = vertex_color_range,
                              edge_color_trans = vertex_color_trans,
                              vertex_label = NA,
                              vertex_label_size = vertex_size,
                              vertex_label_size_range = c(NA, NA),
                              vertex_label_size_trans = vertex_size_trans,
                              vertex_label_color = "#222222",
                              vertex_label_color_range = quantative_palette(),
                              vertex_label_color_trans = "area",
                              edge_label = NA,
                              edge_label_size = edge_size,
                              edge_label_size_range = c(NA, NA),
                              edge_label_size_trans = edge_size_trans,
                              edge_label_color = "#555555",
                              edge_label_color_range = quantative_palette(),
                              edge_label_color_trans = "area",
                              graph_label = NA,
                              graph_label_size_range = c(NA, NA),
                              graph_label_color = "#000000",
                              graph_label_color_range = quantative_palette(),
                              graph_label_color_trans = "area",
                              vertex_label_max = 20,
                              edge_label_max = 20,
                              graph_label_max = 20,
                              overlap_bias = 1,
                              margin_size = c(0, 0),
                              # aspect_ratio = NULL,
                              layout = "reingold-tilford",
                              initial_layout = "fruchterman-reingold") {
  #| ### Verify arguments =========================================================================
  if (length(taxon_id) != length(parent_id)) {
    stop("'taxon_id' and 'parent_id' must be of equal length.")
  }
  if (length(taxon_id) == 0) {
    stop("'taxon_id' and 'parent_id' are empty.")
  }
  check_element_length(c("vertex_size", "edge_size", "vertex_label_size", "edge_label_size",
                         "vertex_color", "edge_color", "vertex_label_color", "edge_label_color",
                         "vertex_label", "edge_label", "graph_label", "graph_label_color"))
  verify_size(c("vertex_size", "edge_size", "vertex_label_size", "edge_label_size"))
  verify_size_range(c("vertex_size_range", "edge_size_range",
                      "vertex_label_size_range", "edge_label_size_range"))
  verify_trans(c("vertex_size_trans", "vertex_color_trans", 
                 "edge_size_trans", "edge_color_trans",
                 "vertex_label_size_trans", "vertex_label_color_trans",
                 "edge_label_size_trans", "edge_label_color_trans", "graph_label_color_trans"))
  verify_color_range(c("vertex_color_range", "edge_color_range",
                       "vertex_label_color_range", "edge_label_color_range",
                       "graph_label_color_range"))
  verify_label_count(c("vertex_label_max", "edge_label_max", "graph_label_max"))
  if (length(overlap_bias) == 0 || ! is.numeric(overlap_bias)) {
    stop("Argument 'overlap_bias' must be a numeric of length 1.")
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
                     glc = graph_label_color,
                     vl = as.character(vertex_label),
                     el = as.character(edge_label),
                     gl = as.character(graph_label))
  row.names(data) <- data$tid
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
                                                      intitial_coords = intitial_coords))
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
    (1 + range_size + minimum * minimum_weight) / (1 + overlap * overlap_bias * overlap_weight)
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
  graph_label_size_range_g <- infer_size_range(graph_label_size_range, range(data$subgraph_prop), 
                                              defualt_scale = 0.1)
  data$vls_g <- scales::rescale(data$vls_t, to = vertex_label_size_range_g)
  data$els_g <- scales::rescale(data$els_t, to = edge_label_size_range_g)
  data$gls_g <- scales::rescale(data$gls_t, to = graph_label_size_range_g)
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
                           "glc_t" = graph_label_color_range,
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
  x_min <- min(element_data$x) - margin_size[1] * square_side_length
  x_max <- max(element_data$x) + margin_size[1] * square_side_length
  y_min <- min(element_data$y) - margin_size[2] * square_side_length
  y_max <- max(element_data$y) + margin_size[2] * square_side_length
  data$vls_g <- scales::rescale(data$vls_g, to = c(0, 1), from = c(0, square_side_length))
  data$vlx_g <- scales::rescale(data$vx, to = c(0, 1), from = c(x_min, x_max))
  data$vly_g <- scales::rescale(data$vy, to = c(0, 1), from = c(y_min, y_max))
  if (is.null(vertex_label_max) || is.na(vertex_label_max)) {
    vertex_label_shown <- data$tid
  } else {
    vertex_label_shown <- data[order(data$vs_g, decreasing = TRUE)[1:vertex_label_max], "tid"]
  }
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
  if (is.null(edge_label_max) || is.na(edge_label_max)) {
    edge_label_shown <- data$tid
  } else {
    edge_label_shown <- data[order(data$es_g, decreasing = TRUE)[1:edge_label_max], "tid"]
  }
  # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  if (is.null(graph_label_max) || is.na(graph_label_max) || length(root_indexes) <= graph_label_max) {
    graph_label_shown <- root_indexes
  } else {
    graph_label_shown <- title_data[order(title_data$gls_g, decreasing = TRUE)[1:graph_label_max], "tid"]
  }
  graph_label_grobs <- plyr::dlply(title_data[graph_label_shown, ], "tid",
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
  for (a_grob in graph_label_grobs) {
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
#' 
#' @export
layout_functions <- function(name = NULL, graph = NULL, intitial_coords = NULL) {
  return_names <- is.null(name) && is.null(graph) && is.null(intitial_coords)
  if (return_names) {
    graph <- make_ring(1) # Dummy graph so that the list can be defined, but only names used
  }
  funcs <- list("automatic" = igraph::nicely(),
                "reingold-tilford" = igraph::as_tree(circular = TRUE,
                                                     mode = "out"),
                "davidson-harel" = igraph::with_dh(coords = intitial_coords,
                                                   maxiter = 15,
                                                   fineiter = max(15, log2(vcount(graph))),
                                                   cool.fact = 0.75,
                                                   weight.node.dist = 1,
                                                   weight.border = 0,
                                                   weight.edge.lengths = 1,
                                                   weight.edge.crossings = 100,
                                                   weight.node.edge.dist = 1),
                "gem" = igraph::with_gem(coords = intitial_coords,
                                         maxiter = 40 * vcount(graph)^2,
                                         temp.max = vcount(graph),
                                         temp.min = 1/10,
                                         temp.init = sqrt(vcount(graph))),
                "graphopt" = igraph::with_graphopt(start = intitial_coords,
                                                   niter = 500,
                                                   charge = 0.0005,
                                                   mass = 30,
                                                   spring.length = 0,
                                                   spring.constant = 1,
                                                   max.sa.movement = 5),
                "mds" = igraph::with_mds(),
                "fruchterman-reingold" = igraph::with_fr(coords = intitial_coords,
                                                         niter = 500,
                                                         start.temp = sqrt(vcount(graph)),
                                                         grid = "nogrid",
                                                         weights = NULL),
                "kamada-kawai" = igraph::with_kk(coords = intitial_coords,
                                                 maxiter = 50 * vcount(graph),
                                                 epsilon = 0,
                                                 kkconst = vcount(graph),
                                                 weights = NULL),
                "large-graph" = igraph::with_lgl(maxiter = 150,
                                                 maxdelta = vcount(graph),
                                                 area = vcount(graph)^2,
                                                 coolexp = 1.5,
                                                 repulserad = vcount(graph)^2 * vcount(graph),
                                                 cellsize = sqrt(sqrt(vcount(graph)^2)),
                                                 root = 1),
                "drl" = igraph::with_drl(use.seed = ! is.null(intitial_coords),
                                         seed = ifelse(is.null(intitial_coords), 
                                                       matrix(runif(vcount(graph) * 2), ncol = 2),
                                                       intitial_coords),
                                         options = drl_defaults$default,
                                         weights = E(graph)$weight,
                                         fixed = NULL))
  if (return_names) {
    return(names(funcs)) # Dummy graph so that the list can be defined, but only names used
  } else {
    return(funcs[[name]])
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