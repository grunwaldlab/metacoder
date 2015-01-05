#===================================================================================================
#' plot_threshold_optimization
#' 
#' @import reshape2
#' @import ggplot2 
plot_threshold_optimization <- function(input, title=NULL, save_png=NULL, display=FALSE, background="transparent") {
  if (class(input) == "character" || class(input) == "factor" ) {
    if (file.exists(as.character(input))) {
      data <- read.csv(as.character(input), sep="\t")
    } else {
      stop("Cannot read input file.")
    }
  } else if (class(input) == "data.frame") {
    data <- input
  } else {
    stop("Invalid input class. data.frame required.")
  }
  comparisons <- sum(data[1, c(2:5)])
  data[,c(2:6)] <- data[,c(2:6)] / comparisons
  error_at_max_x <- ((max(data$cumulative_error) - min(data$cumulative_error)) / 4) +  min(data$cumulative_error)
  max_index <- which(data$false_negative > error_at_max_x)[1]
  if (!is.na(max_index)) {
    data <- data[1:max_index, ]    
  }
  data <- melt(data, measure.vars=4:6, id.vars = 1,  na.rm=TRUE)
  my_plot <- ggplot(data[data$variable != "cumulative_error", ], aes(x=threshold, y=value)) + 
    geom_area(aes(fill = variable), alpha = .3, position='identity') +
    geom_line(data=data[data$variable == "cumulative_error", ], position='identity') +
    labs(title=title) +
    #     scale_x_continuous(limits = c(min_x_display, max_x_display)) +
    #     scale_y_continuous(limits = c(min_y_display, max_y_display)) +
    theme(title=element_text(size=30),
          axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())  
  if (!is.null(save_png)) {
    png(file = save_png, bg = background)
    print(my_plot)
    dev.off()
  }
  if (display) {
    print(my_plot)
  }
  return(my_plot)
}


#===================================================================================================
#' plot_distance_distribution
#' 
#' @importFrom zoo rollmean
#' @import reshape2
#' @import ggplot2 
plot_distance_distribution <- function(input, title=NULL, save_png=NULL, display=FALSE, background="transparent", smoothness=3, bin_width=NULL) {
  if (class(input) == "character" || class(input) == "factor" ) {
    if (file.exists(as.character(input))) {
      data <- read.csv(as.character(input), sep="\t")
    } else {
      stop("Cannot read input file.")
    }
  } else if (class(input) == "data.frame") {
    data <- input
  } else {
    stop("Invalid input class. data.frame required.")
  }
  #infer bin width from adjacent count values
  if (is.null(bin_width)) {
    bin_width <- data$count_middle[2] - data$count_middle[1]     
  }
  #scale same and diffent measurments to equal area relative to total (invalidates y axis)
  if (smoothness > nrow(data)) {
    smoothness <-  nrow(data)
  }
  if (sum(data$same) != 0 && sum(data$different) != 0) {    
    data$same <- data$same * (sum(data$total) / sum(data$same))
    data$different <- data$different * (sum(data$total) / sum(data$different))
  }
  if (nrow(data) > smoothness && sum(data$same) != 0 && sum(data$different) != 0) {    
    #Apply a moving average relative to differences in standard deviation
    different_window <- (1 + as.integer(sd(data$different) / sd(data$same))) * smoothness
    same_window <- (1 + as.integer(sd(data$same) / sd(data$different))) * smoothness
  } else {
    same_window <- smoothness
    different_window <- smoothness
  }
  if (different_window >= nrow(data)) {
    different_window <- nrow(data) - 1
  }
  if (same_window >= nrow(data)) {
    same_window <- nrow(data) - 1
  }
  if (different_window > 1 && sum(data$different) != 0) {
    data$different <- rollmean(data$different, different_window, fill="extend")
  }
  if (same_window > 1 && sum(data$same) != 0) { 
    data$same <- rollmean(data$same, same_window, fill="extend")
  }
  #plot data
  data <- melt(data, measure.vars=2:ncol(data), id.vars = 1,  na.rm=TRUE)
  if (all(c("different", "same") %in% data$variable)) { #if the same and different columns are present (typical)
    data <- data[data$variable != "total", ]
    my_plot <- ggplot(data, aes(x=count_middle, y=value)) + 
      geom_bar(aes(fill = variable), alpha = .5, position='identity', stat="identity", width=bin_width)
  } else {
    my_plot <- ggplot(data, aes(x=count_middle, y=value)) + 
      geom_bar(fill="grey", alpha = .7, position='identity', stat="identity", width=bin_width)
  }
  my_plot <- my_plot + 
    labs(title=title) +
    theme(title=element_text(size=30),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())  
  if (!is.null(save_png)) {
    png(file = save_png, bg = background)
    print(my_plot)
    dev.off()
  }
  if (display) {
    print(my_plot)
  }
  return(my_plot)
}


#===================================================================================================
#' Makes igraph from images
#' 
#' @import igraph 
plot_image_tree <- function(graph, image_file_paths, labels=NA, scaling=1, exclude=c(), root_index=1, label_color = "black") {
  #store the distance of all verticies and edges from the root
  root <- V(graph)[root_index]
  vertex_depth <- sapply(get.shortest.paths(graph, from=root)$vpath, length)
  edge_depth <- vertex_depth[get_edge_parents(graph)]
  
  #set vertex graphing parameters
  V(graph)$size <- (log(scaling + .5) / max(log(scaling) + .5)) * 10
  if (is.na(labels)) {
    V(graph)$label.cex <- 0
  } else {
    V(graph)$label <- labels
    V(graph)$label.cex <- V(graph)$size * .05 + .15
    V(graph)$label.color <- label_color
  }
  V(graph)$alpha <- (max(vertex_depth)*1.5 - vertex_depth) / (max(vertex_depth)*1.5)
  V(graph)$raster_file <- image_file_paths #not used in disaply, but should be subset below
  
  #set edge graphing parameters
  E(graph)$width <- V(graph)$size[get_edge_children(graph)] * 5
  E(graph)$color <- sapply(((max(edge_depth)*4 - edge_depth) / (max(edge_depth)*4)) * .3,
                           function(x) rgb(red=.3,green=.3,blue=.3,alpha=x))
  
  #exclude specific verticies and their decendents from display
  graph <- delete_vetices_and_children(graph, exclude)
  
  #Calculate vertex layout
  graph_layout <- layout.reingold.tilford(graph, root = root_index, circular = TRUE)
  
  #Load vertex images 
  V(graph)$raster <- lapply(as.character(V(graph)$raster_file), readPNG)
  
  #plot graph
  my_plot <- plot(graph,
                  layout=graph_layout,
                  margin=0, 
                  vertex.label.dist=0,
                  vertex.label.degree=0,
                  vertex.label=labels,
                  edge.arrow.size =0,
                  vertex.shape="raster", 
                  vertex.size=V(graph)$size*1.5,
                  vertex.size2=V(graph)$size*1.5)
  if (display) {
    print(my_plot)
  }
  return(plot)
}


#===================================================================================================
#' Makes igraph from values
#' 
#' @param graph An igraph graph object or an adjacency list (i.e. edge list) represented by
#' something that can be coerced to a two column \code{data.frame} where the first column is 
#' every unique taxon id and the second is its parent/super taxon id. 
#' 
#' @import igraph grid
#' @importFrom plotrix color.scale
#' @export
plot_value_tree <- function(graph, values, labels=NA, scaling=1, exclude=c(), root_index=1,
                            label_color = "black", display=FALSE, fade=FALSE, legend_text="",
                            value_range=c(0,1), highlight_outliers=TRUE, background="#00000000",
                            save=NULL) {
  #   init_igraph()
  #make graph if nessesary
  if (class(graph) != "igraph") {
    graph <- graph[complete.cases(graph), ]
    graph <- graph.edgelist(as.matrix(as.data.frame(graph)))
  }   
  #store the distance of all verticies and edges from the root
  root <- V(graph)[root_index]
  vertex_depth <- sapply(get.shortest.paths(graph, from=root)$vpath, length)
  edge_depth <- vertex_depth[get_edge_parents(graph)]
  
  #set vertex graphing parameters
  V(graph)$size <- (log(scaling + .5) / max(log(scaling) + .5)) * 10
  if (is.na(labels)) {
    V(graph)$label.cex <- 0
    V(graph)$label <- NA
  } else if (labels == TRUE) {
    V(graph)$label <- as.character(signif(values, 2))
    V(graph)$label.cex <- V(graph)$size * .45 + .15
    V(graph)$label.color <- label_color
  } else {
    V(graph)$label <- labels
    V(graph)$label.cex <- V(graph)$size * .45 + .15
    V(graph)$label.color <- label_color
  }
  if (fade == TRUE) {
    V(graph)$alpha <- (max(vertex_depth)*1.5 - vertex_depth) / (max(vertex_depth)*1.5)
  } else if (fade == FALSE) {
    V(graph)$alpha <- 1
  } else {
    V(graph)$alpha <- fade
  }
  V(graph)$values <- values
  
  #set edge graphing parameters
  E(graph)$width <- V(graph)$size[get_edge_children(graph)] * 5
  E(graph)$color <- sapply(((max(edge_depth)*4 - edge_depth) / (max(edge_depth)*4)) * .3,
                           function(x) rgb(red=.3,green=.3,blue=.3,alpha=x))
  
  #exclude specific verticies and their decendents from display
  graph <- delete_vetices_and_children(graph, exclude)
  
  #set vertex color
  color_values <- V(graph)$values
  value_range_quantile <- quantile(color_values, value_range, na.rm=TRUE)
  #   if (highlight_outliers) {
  #     outliers <- color_values < value_range_quantile[1] | color_values > value_range_quantile[2]
  #     V(graph)$frame.width <- ifelse(outliers, 25, .05)    
  #   }
  outliers <- color_values < value_range_quantile[1] | color_values > value_range_quantile[2]
  #   V(graph)$frame.width <- ifelse(outliers, 25, .05)
  color_values[color_values < value_range_quantile[1]] <- value_range_quantile[1]
  color_values[color_values > value_range_quantile[2]] <- value_range_quantile[2]
  no_color_in_palette <- 10000
  palette <- colorRampPalette(c("red","green", "blue"))(no_color_in_palette)
  displayed_colors <- palette[as.integer((no_color_in_palette - 1) * color_values/max(color_values)) + 1]
  V(graph)$color=mapply(add_alpha, displayed_colors, alpha=V(graph)$alpha)
  
  #Calculate vertex layout
  graph_layout <- layout.reingold.tilford(graph, root = root_index, circular = TRUE)
  
  #Load vertex images 
  #   V(graph)$raster <- lapply(as.character(V(graph)$raster_file), readPNG)
  
  if (!is.null(save)) {
    png(file = save, bg = background, width=5000, height=5000)
  }
  
  #Make legend (http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2)
  legend <- continuous_color_legend(color_values,
                                    low=V(graph)$color[which.min(color_values)], 
                                    mid=V(graph)$color[which_middle(color_values)], 
                                    high=V(graph)$color[which.max(color_values)],
                                    name=legend_text,
                                    background=background)
  
  main_vp <- viewport(x = .5, y = .5, width = 1, height = 1)
  legend_vp <- viewport(x = 1, y = 0, width = 0.1, height = 1, just = c("right", "bottom"))
  grid.newpage()
  pushViewport(main_vp)
  plot(graph, 
       layout=graph_layout,
       margin=0, 
       vertex.label.dist=0,
       vertex.label.degree=0,
       edge.arrow.size =0,
       vertex.shape="circle", 
       vertex.frame.color=NA)
  upViewport()
  pushViewport(legend_vp)
  grid.draw(legend[[1]][[1]])
  popViewport(1)
  if (!is.null(save)) {
    dev.off()
  }
}

#===================================================================================================
#' plot_value_distribution_by_level
#' 
#' @import ggplot2 
plot_value_distribution_by_level <- function(taxon_data, value_column, level_column = "level", ...) {
  ggplot(taxon_data, aes_string(x=level_column, y=value_column)) + 
    geom_boxplot(width=.5, outlier.colour="transparent") +
    geom_violin(alpha=.4, aes_string(fill=level_column,  colour = level_column)) + 
    geom_point(position = position_jitter(w = 0.1, h = 0), alpha = .2) + 
    facet_grid(~ clustering_level, scales = "free_x", space="free_x")  +
    labs(...) +
    theme(title=element_text(size=17),
          axis.text.y=element_text(size=12),
          axis.text.x=element_text(size=12, angle = 60, hjust = 1),
          axis.title=element_text(size=20),
          legend.position="none",
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          strip.text.x = element_text(size = 12))
}

#===================================================================================================
#' Plot items on taxonomy
#' 
#' Plots the distribution of values associated with item with an associated taxonomic classification.
#' Uses \code{igraph} to make layout and \code{ggplot2} to make plots.
#' 
#' @param data (\code{data.frame}) Must have columns with names specified by options \code{taxon_id}
#' and \code{parent_id}.
#' @param taxon_id  The unique ids of the taxon for each row.
#' @param parent_id The unique id of supertaxon \code{taxon_id} is a part of.
#' @param size The value to base vertex and line size on. If a numeric vector is given, it is rescaled
#' to an optimum size for display based on minimizing overlaps and spaces beteween adjacent verticies.
#' See option \code{overlap_bias} for how this optimization is done.
#' @param vertex_color The value to base vertex color on. If a numeric vector is given, it is used to 
#' construct a color scale. Hex values or color names can be used (e.g. \code{#000000} or
#' \code{"black"}). 
#' @param vertex_label The values of labels over vertcies. 
#' @param line_color The value to base line color on. If a numeric vector is given, it is used to 
#' construct a color scale. Hex values or color names can be used (e.g. \code{#000000} or
#' \code{"black"}). 
#' @param line_label The values of labels over lines. 
#' @param overlap_bias (\code{numeric} > 0) The factor by which overlaps are punished relative to
#' spaces when optimizing displayed vertex size range. Only nearby pairs of verticies are considered. 
#' For example, a value of \code{10} would mean that overlaps are considered 10 time more important than
#' gaps to avoid amoung nearby verticies when choosing an optimum vertex size range. 
#' @param min_label_size (\code{numeric} of length 1) The minimum label size that will be shown. A
#' proportion of viewport size. Labels that would be smaller are not added to save time. The type
#' of text grob used to display the labels that scale with output dimensions is time consuming to
#' draw, so printing all labels can take a long time. Giving a value of \code{0} will cause all labels
#' to be printed.
#' @param line_label_offset (\code{numeric}) The length the line label is offset from its associated
#' vertex centerin multiples of vertex radius. For example, a value of \code{1} will place the label on the
#' edge of the vertex circle and \code{0} would place it in the center.
#' @param margin_size (\code{numeric} < 0.5) The amount of space around the plot in terms of the proportion to
#' the range of space occupied by verticies. For example, a value of \code{.5} would make margins as wide
#' as the plotted data.
#' @param aspect_ratio (\code{numeric}) The height / width of the plot.
#' 
#' @import grid
#' @import ggplot2 
#' @import scales 
#' @export
plot_taxonomy <- function(taxon_id, parent_id, size = NULL, vertex_color = NULL, vertex_label = NULL, 
                          line_color = NULL, line_label = NULL, overlap_bias = 15, min_label_size = .015,
                          line_label_offset = 1, margin_size = 0.1, aspect_ratio = NULL, data_only = FALSE,
                          layout_func = NULL, layout_args = NULL, titles = NULL) {
  # Validate arguments -----------------------------------------------------------------------------
  if (length(taxon_id) != length(parent_id)) stop("unequal argument lengths")
  parent_id[!(parent_id %in% taxon_id)] <- NA
  if (is.null(layout_func)) {
    layout_func <- layout.reingold.tilford
    if (is.null(layout_args)) layout_args <- list(circular = TRUE)
  }
  if (is.null(layout_args)) layout_args <- list()
  # Get vertex coordinants  ------------------------------------------------------------------------
  get_vertex_coords <- function(index) {
    if (length(index) == 1) return(data.frame(x = 0, y = 0))
    part <- data[index, ]
    graph <- graph.edgelist(as.matrix(part[complete.cases(part), c("parent_id", "taxon_id")]))
    layout <- do.call(layout_func, c(list(graph), layout_args))
    if (!is.null(aspect_ratio)) layout[, 1] <- layout[, 1] * aspect_ratio
    if (any(is.na(layout) | is.nan(unlist(layout)))) {
      layout <- layout.fruchterman.reingold(graph)
      warning(paste('Could not apply layout_func to subgraph with root', part$taxon_id[1]))
    }
    return(list(graph, layout))
  }
  data <- data.frame(taxon_id = taxon_id, parent_id = parent_id)
  rownames(data) <- taxon_id
  subgraphs <- split_by_level(taxon_id, parent_id, level =  1)
  layouts <- lapply(subgraphs, get_vertex_coords)
  layout <- layout.merge(graphs = lapply(layouts, `[[`, 1), layouts = lapply(layouts, `[[`, 2))
  coords <- setNames(as.data.frame(layout), c('x', 'y'))
  data <- cbind(data, coords)
  data$group <- rep(seq_along(subgraphs), vapply(subgraphs, length, numeric(1)))
  if (!is.null(aspect_ratio)) {
    data$x <- data$x / aspect_ratio
    data$y <- rescale(data$y, to = range(data$x, na.rm = TRUE) * aspect_ratio)
  }
  # Get vertex size --------------------------------------------------------------------------------
  if (is.null(size)) {
    data$depth <- edge_list_depth(data$taxon_id, data$parent_id)
    data$size <- (max(data$depth) - data$depth + 1)^5
  } else {
    data$size <- size
    data$size <- sqrt(data$size / pi)
  }
  all_pairwise <- molten_dist(x = data$x, y = data$y)
  smallest_side <- min(c(max(data$x) - min(data$x), max(data$y) - min(data$y)))
  max_range <- c(min(all_pairwise$distance), smallest_side / 5)
  min_range <- c(min(all_pairwise$distance) / 5, min(all_pairwise$distance))
  pairwise <- all_pairwise[all_pairwise$distance <= max_range[2], ]
  size_opt_func <- function(a_max, a_min) {
    # Get pairwise distance metrics  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pairs <- pairwise
    size <- rescale(data$size, to = c(a_min, a_max))
    pairs$gap <- pairs$distance - size[pairs$index_1] - size[pairs$index_2]
    # Calculate optimality metric  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    overlap <- sum(abs(pairs$gap[pairs$gap < 0]))
    space <- sum(pairs$gap[pairs$gap >= 0])
    return(c(overlap, space))
  }
  choose_best <- function(options) {
    data <- as.data.frame(do.call(rbind, options))
    names(data) <- c("overlap", "space")
    data$overlap <- rescale(data$overlap, to = c(0, overlap_bias))
    data$space <- rescale(data$space, to = c(0, 1))
    data$score <- data$overlap  + data$space
    which.min(data$score)
  }
  opt_size_range <- get_optimal_range(max_range = max_range,
                                      min_range = min_range,
                                      resolution = c(15, 20),
                                      opt_crit = size_opt_func, 
                                      choose_best = choose_best)
  data$size <- rescale(data$size, to = opt_size_range)
  vertex_data <- polygon_coords(n = 50, x = data$x, y = data$y, radius = data$size)
  vertex_data <- cbind(vertex_data, data[as.numeric(vertex_data$group), c("size"), drop = F])
  # Get graph range data ---------------------------------------------------------------------------
  x_range <- max(vertex_data$x) - min(vertex_data$x)
  y_range <- max(vertex_data$y) - min(vertex_data$y)
  x_margin <- x_range * margin_size
  y_margin <- y_range * margin_size
  x_display <- x_range + 2 * x_margin
  y_display <- y_range + 2 * y_margin
  x_min <-  min(vertex_data$x) - x_margin
  x_max <- max(vertex_data$x) + x_margin
  y_min <- min(vertex_data$y) - y_margin
  y_max <- max(vertex_data$y) + y_margin
  ideal_diameter <- sqrt(x_display * y_display)
  # Get edge coordinants ---------------------------------------------------------------------------
  data$parent_x <- data$x[match(data$parent_id, data$taxon_id)]  
  data$parent_y <- data$y[match(data$parent_id, data$taxon_id)]
  line_data <- line_coords(x1 = data$x, y1 = data$y, x2 = data$parent_x, y2 = data$parent_y,
                           width = data$size)
  line_data <- cbind(line_data, data[as.numeric(line_data$group), c("size"), drop = F])
  # Get vertex color -------------------------------------------------------------------------------
  if (is.null(vertex_color)) {
    data$vertex_color <- "grey"
  } else {
    data$vertex_color <- vertex_color
  }
  if (is.numeric(data$vertex_color)) { ## Not factors or hex codes
    no_color_in_palette <- 10000
    palette <- colorRampPalette(c("red","green", "blue"))(no_color_in_palette)
    color_index <- (no_color_in_palette - 1) * data$vertex_color/max(data$vertex_color) + 1
    data$vertex_color <- palette[as.integer(color_index)]    
  }
  vertex_data$vertex_color <- data$vertex_color[as.numeric(vertex_data$group)]    
  # Get edge color ---------------------------------------------------------------------------------
  if (is.null(line_color)) {
    data$line_color <- data$vertex_color
  } else {
    data$line_color <- line_color
  }
  if (is.numeric(data$line_color)) { ## Not factors or hex codes
    no_color_in_palette <- 10000
    palette <- colorRampPalette(c("red","green", "blue"))(no_color_in_palette)
    color_index <- (no_color_in_palette - 1) * data$line_color/max(data$line_color) + 1
    data$line_color <- palette[as.integer(color_index)]    
  }
  line_data$line_color <- data$line_color[as.numeric(line_data$group)]
  # Get vertex labels ------------------------------------------------------------------------------
  if (!is.null(vertex_label)) {
    data$vertex_label <- as.character(vertex_label)
    data$vertex_label_x <- rescale(data$x, to = c(0, 1), from = c(x_min, x_max))
    data$vertex_label_y <- rescale(data$y, to = c(0, 1), from = c(y_min, y_max))
    data$vertex_label_size <-  rescale(data$size, to = c(0, 1), from = c(0, ideal_diameter))
    vertex_label_grobs <- lapply(which(data$vertex_label_size > min_label_size), 
                                 function(i) resizingTextGrob(label = data$vertex_label[i],
                                                              y = data$vertex_label_y[i],
                                                              x = data$vertex_label_x[i],
                                                              gp = gpar(text_prop = data$vertex_label_size[i])))    
  }
  # Get line labels --------------------------------------------------------------------------------
  if (!is.null(line_label)) {
    data$line_label <- as.character(line_label)
    data$line_label[is.na(data$parent_id)] <- ""
    # line label rotation  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    data$slope <- (data$y - data$parent_y) / (data$x - data$parent_x)
    data$slope[is.na(data$slope)] <- 0
    data$line_label_rot <- atan(data$slope)
    # line label coordinate  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    justify <- data$parent_x > data$x
    justify[is.na(justify)] <- TRUE
    justification <- lapply(1:nrow(data), function(i) if (justify[i]) c("left", "center") else c("right", "center"))
    line_label_x_offset <- line_label_offset * data$size * cos(data$line_label_rot)
    line_label_y_offset <- line_label_offset * data$size * sin(data$line_label_rot)
    data$line_label_x <-  data$x + ifelse(justify, 1, -1) * line_label_x_offset
    data$line_label_y <- data$y + ifelse(justify, 1, -1) * line_label_y_offset
    data$line_label_x <- rescale(data$line_label_x, to = c(0, 1), from = c(x_min, x_max))
    data$line_label_y <- rescale(data$line_label_y, to = c(0, 1), from = c(y_min, y_max))
    # line label text size - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    mean_inter_pair <- mean(sqrt((data$x - data$parent_x)^2 + (data$y - data$parent_y)^2), na.rm = TRUE)
    mean_inter_pair <- rescale(mean_inter_pair, to = c(0, 1), from = c(0, mean(c(x_display, y_display)))) 
    data$line_label_size <-  rescale(data$size / 2, to = c(0, 1), from = c(0, mean(c(x_display, y_display)))) 
    max_line_label_size <- mean_inter_pair / 10
    data$line_label_size[data$line_label_size > max_line_label_size] <-  max_line_label_size
    # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    line_label_grobs <- lapply(which(data$line_label_size > min_label_size / 3), 
                               function(i) resizingTextGrob(label = data$line_label[i],
                                                            y = data$line_label_y[i],
                                                            x = data$line_label_x[i],
                                                            rot = data$line_label_rot[i] * 180 / pi,
                                                            just = justification[[i]],
                                                            gp = gpar(text_prop = data$line_label_size[i])))
  }
  # Subgraph titles --------------------------------------------------------------------------------
  if (!is.null(titles)) {
    data$title <- titles
    root_index <- vapply(unique(data$group), function(x) which(data$group == x)[1], numeric(1))
    title_data <- data.frame(label = data$title[root_index])
    title_data$x <- vapply(unique(data$group), function(x) mean(range(data$x[data$group == x])), numeric(1))
    title_data$size <- vapply(unique(data$group), function(x) (mean(c(max(data$x[data$group == x]) -  min(data$x[data$group == x]), max(data$y[data$group == x]) -  min(data$y[data$group == x])))) * 0.075, numeric(1))
    title_data$y <- vapply(unique(data$group), function(x) max(data$y[data$group == x]), numeric(1)) + title_data$size * 1.1
    title_data$x <- rescale(title_data$x, to = c(0, 1), from = c(x_min, x_max))
    title_data$y <- rescale(title_data$y, to = c(0, 1), from = c(y_min, y_max))
    title_data$size <- rescale(title_data$size, to = c(0.001, 1), from = c(0, mean(c(x_display, y_display))))
    # create text grobs  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    title_grobs <- lapply(which(title_data$size > min_label_size / 4), 
                          function(i) resizingTextGrob(label = title_data$label[i],
                                                       y = title_data$y[i],
                                                       x = title_data$x[i],
                                                       rot = 0,
                                                       just = "center",
                                                       gp = gpar(text_prop = title_data$size[i])))
    
  }
  
  # Get graph range data ---------------------------------------------------------------------------
  if (!is.null(aspect_ratio)) {
    line_data$y <- rescale(line_data$y, to = range(line_data$x, na.rm = TRUE) * aspect_ratio)
    vertex_data$y <- rescale(vertex_data$y, to = range(vertex_data$x, na.rm = TRUE) * aspect_ratio)
    data$y <- rescale(data$y, to = range(data$x, na.rm = TRUE) * aspect_ratio)    
  }
  x_range <- max(vertex_data$x) - min(vertex_data$x)
  y_range <- max(vertex_data$y) - min(vertex_data$y)
  x_margin <- x_range * margin_size
  y_margin <- y_range * margin_size
  x_min <-  min(vertex_data$x) - x_margin
  x_max <- max(vertex_data$x) + x_margin
  y_min <- min(vertex_data$y) - y_margin
  y_max <- max(vertex_data$y) + y_margin
  
  # Plot it! ---------------------------------------------------------------------------------------
  if (data_only) {
    return(list(taxon_data = data, vertex_data = vertex_data, line_data = line_data))
  } else {
    the_plot <- ggplot(data = data) +
      geom_polygon(data = line_data, aes(x = x, y = y, group = group), fill = line_data$line_color) +
      geom_polygon(data = vertex_data, aes(x = x, y = y, group = group), fill = vertex_data$vertex_color) +
      guides(fill = "none") +
      coord_fixed(xlim = c(x_max, x_min), ylim = c(y_max, y_min)) +
      theme(panel.grid = element_blank(), 
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text  =  element_blank(),
            axis.ticks = element_blank(), 
            axis.line  = element_blank())
    if (!is.null(data$vertex_label)) {
      for (a_grob in vertex_label_grobs) {
        the_plot <- the_plot + annotation_custom(grob = a_grob)
      }    
    }
    if (!is.null(data$line_label)) {
      for (a_grob in line_label_grobs) {
        the_plot <- the_plot + annotation_custom(grob = a_grob)
      }    
    }
    if (!is.null(titles)) {
      for (a_grob in title_grobs) {
        the_plot <- the_plot + annotation_custom(grob = a_grob)
      }    
    }
    return(the_plot)    
  }
 }