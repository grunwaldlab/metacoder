### Generic internal functions

#===================================================================================================
#' under development 
#' 
#' @keywords internal
offset_ordered_factor <- function(ordered_factor, offset) { 
  my_levels <-  levels(ordered_factor)
  new_level <- my_levels[which(my_levels == ordered_factor) + offset]
  ordered(new_level, my_levels)
}

#===================================================================================================
#' under development 
#' 
#' 
#' @keywords internal
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
#' @keywords internal
remove_all_na_rows <- function(input) {
  na_rows <- sapply(1:nrow(input), function(x) sum(!is.na(input[x,])) != 0)
  input[na_rows, ]
}


#===================================================================================================
#' under development 
#' 
#' 
#' @keywords internal
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
#' 
#' @keywords internal
which_median <- function(x) which.min(abs(x - median(x)))


#===================================================================================================
#' under development 
#' 
#' 
#' @keywords internal
which_middle <- function(x) {
  middle <- (max(x) + min(x)) / 2
  which.min(abs(x - middle))
}


### File system functions

#===================================================================================================
#' under development 
#' 
#' 
#' @keywords internal
rm_ext <- function(file) {
  sub("[.][^.]*$", "", file, perl=TRUE)
}


#===================================================================================================
#' under development 
#' 
#' 
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
get_edge_parents <-function(graph) {
  igraph::get.edges(graph, 1:igraph::ecount(graph))[,1]
}


#===================================================================================================
#' get_edge_children
#' 
#' @keywords internal
get_edge_children <- function(graph) {
  igraph::get.edges(graph, 1:igraph::ecount(graph))[,2]
}


#===================================================================================================
#' get_vertex_children
#' 
#' @keywords internal
get_vertex_children <- function(graph, vertex) {
  which(igraph::shortest.paths(graph, igraph::V(graph)[vertex], mode="out") != Inf)
}

#===================================================================================================
#' delete_vetices_and_children
#' 
#' @keywords internal
delete_vetices_and_children <- function(graph, vertices) {
  vertices <- unlist(sapply(vertices, function(x) get_vertex_children(graph, x)))
  graph <- igraph::delete.vertices(graph, vertices)
  return(graph)
}


### Generic ploting functions


#===================================================================================================
#' add_alpha
#' 
#' @keywords internal
add_alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}


### iGraph-associated plotting functions



#===================================================================================================
#' get nodes from leafs
#' 
#' @keywords internal
get_nodes_from_leafs <- function(leafs, sep = ";") {
  if (is.atomic(leafs)) leafs <- strsplit(leafs, sep, fixed=TRUE)
  taxons <- lapply(leafs, function(x) sapply(seq(1, length(x)), function(y) paste(x[1:y], collapse=sep)))
  unique(unlist(taxons))
}


#===================================================================================================
#' count nodes in leafs
#' 
#' @keywords internal
count_nodes_in_leafs <- function(leafs, nodes = NULL, sep = ";") {
  if (is.null(nodes)) nodes <- get_nodes_from_leafs(leafs, sep = sep)
  vapply(nodes, function(x) sum(grepl(x, leafs, fixed = TRUE)), numeric(1))
}

get_tips <- function(nodes, sep = "__") {
  split_nodes <- strsplit(nodes, "__", fixed=TRUE)
  sapply(split_nodes, function(x) x[length(x)])
}

#===================================================================================================
#' get indexes of a unique set of the input
#' 
#' @keywords internal
unique_mapping <- function(input) {
  unique_input <- unique(input)
  vapply(input, function(x) {if (is.na(x)) which(is.na(unique_input)) else which(x == unique_input)}, numeric(1))
}


#===================================================================================================
#' run a function on unique values of a iterable 
#' 
#' @keywords internal
map_unique <- function(input, func, ...) {
  input_class <- class(input)
  unique_input = unique(input)
  class(unique_input) <- input_class
  func(unique_input, ...)[unique_mapping(input)]
}





#===================================================================================================
#' Converts DNAbin to a named character vector
#' 
#' Converts an object of class DNAbin (as produced by ape) to a named character vector.
#' 
#' @param dna_bin (\code{DNAbin} of length 1) the input.
#' 
#' @keywords internal
DNAbin_to_char <- function(dna_bin) {
  vapply(as.character(dna_bin), paste, character(1), collapse="")
}
