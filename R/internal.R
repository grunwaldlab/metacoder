### Generic internal functions

#===================================================================================================
#' under development 
#' 
#' @export
offset_ordered_factor <- function(ordered_factor, offset) { 
  my_levels <-  levels(ordered_factor)
  new_level <- my_levels[which(my_levels == ordered_factor) + offset]
  ordered(new_level, my_levels)
}

#===================================================================================================
#' under development 
#' 
#' @export
fapply <- function(iterable, functions, 
                   .preprocessor={function(x) x},
                   .preprocessor_args=list(),
                   .allow_complex=TRUE,
                   ...) {
  apply_functions <- function(input, functions, ...) {
    if (!is.list(input)) {
      input <- list(input)
    }
    input <- append(input, list(...))
    results <- lapply(functions, function(f) do.call(f, input))
    atomics <- which(!sapply(results, is.recursive))
    if (length(atomics) > 0) {
      results[atomics] <- lapply(1:length(atomics), function(i) {y <- list(results[[atomics[i]]]); 
                                                                 names(y) <- functions[i];
                                                                 y})
    }
    results <- unlist(results, recursive=FALSE)
    if (!.allow_complex) {
      results <- results[!sapply(results, is.recursive)]
    }
    return(results)
  }
  if (length(iterable) < 1) {
    return(NULL)
  }
  if (is.data.frame(iterable) | is.matrix(iterable)) {
    iterable_length <- length(iterable[[1]])    
    row_names <- row.names(iterable)
    call_preprocessor <- function(i) {do.call(.preprocessor, append(list(iterable[i,]), .preprocessor_args))}
  } else if (is.list(iterable)) {
    iterable_length <- length(iterable)
    row_names <- unlist(iterable)
    call_preprocessor <- function(i) {do.call(.preprocessor, append(list(iterable[[i]]), .preprocessor_args))}    
  } else {
    iterable_length <- length(iterable)
    row_names <- iterable
    call_preprocessor <- function(i) {do.call(.preprocessor, append(list(iterable[i]), .preprocessor_args))}        
  }
  output <- lapply(1:iterable_length, function(i) apply_functions(call_preprocessor(i), functions, ...))
  column_names <- names(output[[1]])
  output <- lapply(1:length(output[[1]]), function(i) lapply(output, function(row) row[[i]]))
  output <- lapply(output, function(x) if (!is.recursive(x[[1]])) {unlist(x, recursive=FALSE)} else {I(x)})
  output <- do.call(data.frame, output)
  colnames(output) <- column_names
  row.names(output) <- row_names
  return(output)
}


#===================================================================================================
#' under development 
#' 
#' @export
remove_all_na_rows <- function(input) {
  na_rows <- sapply(1:nrow(input), function(x) sum(!is.na(input[x,])) != 0)
  input[na_rows, ]
}


#===================================================================================================
#' under development 
#' 
#' @export
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.01, .99), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


#===================================================================================================
#' under development 
#' 
#' @export
which_median <- function(x) which.min(abs(x - median(x)))


#===================================================================================================
#' under development 
#' 
#' @export
which_middle <- function(x) {
  middle <- (max(x) + min(x)) / 2
  which.min(abs(x - middle))
}


### File system functions

#===================================================================================================
#' under development 
#' 
#' @export
rm_ext <- function(file) {
  sub("[.][^.]*$", "", file, perl=TRUE)
}


#===================================================================================================
#' under development 
#' 
#' @export
next_incremental_file_number <-function(directory) {
  current_numbers <- as.integer(rm_ext(list.files(directory, no..=TRUE)))
  if (length(current_numbers) == 0) {
    current_numbers = 0
  }
  max(current_numbers) + 1
}

### iGraph-associated functions


#===================================================================================================
#' under development
#' 
#' @export
taxon_edge_list <- function(taxonomy, separator) {
  get_taxon_edge_list <- function(taxon) {
    apply(matrix(c(1:(length(taxon)-1),2:length(taxon)), ncol = 2), 1, function(x) c(taxon[x[1]], taxon[x[2]]))
  }
  taxons <- unique(taxonomy)
  taxons <- strsplit(taxons, separator, fixed=TRUE)
  taxons <- taxons[sapply(taxons, length) > 1]
  taxons <- lapply(taxons, function(x) sapply(seq(1, length(x)), function(y) paste(x[1:y], collapse=separator)))
  edge_list <- t(do.call(cbind,lapply(taxons, FUN=get_taxon_edge_list)))
  edge_list[!duplicated(edge_list),]
}


#===================================================================================================
#' get_edge_parents
#' 
#' @export
#' @importFrom igraph get.edges ecount
get_edge_parents <-function(graph) {
  get.edges(graph, 1:ecount(graph))[,1]
}


#===================================================================================================
#' get_edge_children
#' 
#' @export
#' @importFrom igraph get.edges ecount
get_edge_children <- function(graph) {
  get.edges(graph, 1:ecount(graph))[,2]
}


#===================================================================================================
#' get_vertex_children
#' 
#' @export
#' @importFrom igraph shortest.paths V
get_vertex_children <- function(graph, vertex) {
  which(shortest.paths(graph, V(graph)[vertex], mode="out") != Inf)
}

#===================================================================================================
#' delete_vetices_and_children
#'
#' @export
#' @importFrom igraph delete.vertices
delete_vetices_and_children <- function(graph, vertices) {
  #delete children
  vertices <- unlist(sapply(vertices, function(x) get_vertex_children(graph, x)))
  graph <- delete.vertices(graph, vertices)
  return(graph)
}


### Generic ploting functions


#===================================================================================================
#' add_alpha
#' 
#' @export
add_alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}


### iGraph-associated plotting functions


#===================================================================================================
#' mycircle generator
#' 
#' @export
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

#===================================================================================================
init_igraph <- function() {
  add.vertex.shape("fcircle", 
                   plot=mycircle, 
                   parameters=list(vertex.frame.color=1, vertex.frame.width=1))
}

#===================================================================================================
#' make a color legend
#'
#' @importFrom ggplot2 ggplot_gtable ggplot_build qplot scale_colour_gradient2 theme element_text element_rect
#' @importFrom grid unit
#' @export
continuous_color_legend <- function(values, background="#00000000", ...) {
  #Extract Legend (http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2)
  g_legend<-function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)} 
  mid_point = (max(values) + min(values)) / 2
  label_points <- seq(min(values), max(values), length.out=7)
  labels <- as.character(signif(label_points, 2))
  full_plot <- qplot(x,y, colour=value, data=data.frame(x=1, y=1, value=values)) + 
    scale_colour_gradient2(breaks = label_points, labels = labels, midpoint=mid_point, ...) + 
    theme(legend.key.size = unit(5, "cm"), 
          legend.text = element_text(size=85),
          legend.title = element_text(size=85),
          legend.background = element_rect(fill = background))
  g_legend(full_plot)
}


#===================================================================================================
#' get nodes from leafs
#' 
#' @export
get_nodes_from_leafs <- function(leafs, sep = ";") {
  if (is.atomic(leafs)) leafs <- strsplit(leafs, sep, fixed=TRUE)
  taxons <- lapply(leafs, function(x) sapply(seq(1, length(x)), function(y) paste(x[1:y], collapse=sep)))
  unique(unlist(taxons))
}


#===================================================================================================
#' count nodes in leafs
#' 
#' @export
count_nodes_in_leafs <- function(leafs, nodes = NULL, sep = ";") {
  if (is.null(nodes)) nodes <- get_nodes_from_leafs(leafs, sep = sep)
  vapply(nodes, function(x) sum(grepl(x, leafs, fixed = TRUE)), numeric(1))
}

get_tips <- function(nodes, sep = "__") {
  split_nodes <- strsplit(nodes, "__", fixed=TRUE)
  sapply(split_nodes, function(x) x[length(x)])
}
