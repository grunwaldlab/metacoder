#===================================================================================================
#' plot_threshold_optimization
#' 
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes geom_area geom_line labs theme element_text element_blank
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
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes geom_area geom_line labs theme element_text element_blank geom_bar
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
#' @importFrom igraph graph.edgelist V get.shortest.paths E
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
#' @importFrom igraph graph.edgelist V get.shortest.paths E
#' @importFrom plotrix color.scale
#' @export
plot_value_tree <- function(graph, values, labels=NA, scaling=1, exclude=c(), root_index=1,
                            label_color = "black", display=FALSE, fade=FALSE, legend_text="",
                            value_range=c(0,1), highlight_outliers=TRUE, background="#00000000",
                            save=NULL) {
  init_igraph()
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
  if (highlight_outliers) {
    outliers <- color_values < value_range_quantile[1] | color_values > value_range_quantile[2]
    V(graph)$frame.width <- ifelse(outliers, 25, .05)    
  }
  outliers <- color_values < value_range_quantile[1] | color_values > value_range_quantile[2]
  V(graph)$frame.width <- ifelse(outliers, 25, .05)
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
  #plot graph
  my_plot <- plot(graph,
                  layout=graph_layout,
                  margin=0, 
                  vertex.label.dist=0,
                  vertex.label.degree=0,
                  edge.arrow.size =0,
                  vertex.shape="fcircle", 
                  vertex.frame.color='black')
  
  #Make legend (http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2)
  legend <- continuous_color_legend(color_values,
                                    low=V(graph)$color[which.min(color_values)], 
                                    mid=V(graph)$color[which_middle(color_values)], 
                                    high=V(graph)$color[which.max(color_values)],
                                    name=legend_text,
                                    background=background)  
  pushViewport(viewport(x=0.9, y=0.15))
  grid.draw(legend)
  
  if (!is.null(save)) {
    dev.off()
  }
  return(my_plot)
}

#===================================================================================================
#' plot_value_distribution_by_level
#' 
#' @importFrom ggplot2 aes_string geom_boxplot geom_violin geom_point facet_grid position_jitter labs theme element_text element_blank
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
