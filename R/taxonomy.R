#===================================================================================================
#' Generate new unique IDs
#' 
#' Makes a vector of unique IDs that differ from a previously defined set of unique IDs. 
#' @param count (\code{numeric} of length 1) The number of new unique IDs to generate
#' @param existing (\code{character}) Existing unique IDs. These will not appear in the output.
#' 
#' @keywords internal
make_new_ids <- function(count, existing) {
  output <- rep(NA, count)
  current <- 1
  while (any(is.na(output))) {
    if (!(as.character(current) %in% existing)) {
      output[which(is.na(output))[1]] <- as.character(current)
    }
    current <-  1 + current
  }
  return(output)
}


#===================================================================================================
#' Get classification for taxa in edge list
#' 
#' Extracts the classification of every taxon in a list of unique taxa and their supertaxa.
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon.
#' Root taxa should have \code{NA} in this column.
#' 
#' @return A list of vectors of taxa IDs. Each list entry corresponds to the \code{taxa} supplied.
#' 
#' @keywords internal
get_class_from_el <- function(taxa, parents) {
  process_one <- function(x) {
    output <- character(0)
    my_next <- x
    while (length(my_next) != 0 && !is.null(my_next) && !is.na(my_next)) {
      output <- c(my_next, output)
      my_next <- parents[taxa == my_next]
    }
    return(output)
  }
  if (!is.character(taxa)) taxa <- as.character(taxa)
  if (!is.character(parents)) parents <- as.character(parents)  
  setNames(lapply(taxa, process_one), taxa)
}


#===================================================================================================
#' Get distance from root of edgelist items
#' 
#' Gets the number of ancestors/supergroups for items of an edge/adjacency list
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon.
#' Root taxa should have \code{NA} in this column.
#' 
#' @keywords internal
edge_list_depth <-  function(taxa, parents) {
  vapply(get_class_from_el(taxa, parents), length, numeric(1))
}


#===================================================================================================
#' Get ordered ranks from taxonomy
#' 
#' Returns an ordered factor of the ranks of a given taxonomy.
#' It also checks if rank corresponds consistently with level.
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon.
#' Root taxa should have \code{NA} in this column.
#' @param rank (\code{character}) The rank designation (e.g. "genus") corresponding to each item in
#' \code{taxa}.
#' @param strict If \code{FALSE}, ranks with inconsistent levels will be allowed. Otherwise ranks 
#' with overlapping level ranges will cause an error. 
#' 
#' @keywords internal
taxonomy_ranks <- function(taxa, parents, rank, strict = TRUE) {
  # Get rank data ----------------------------------------------------------------------------------
  level_by_rank <- split(edge_list_depth(taxa, parents), rank)
  rank_data <- data.frame(rank = names(level_by_rank),
                          max = vapply(level_by_rank, max, numeric(1)),
                          min = vapply(level_by_rank, min, numeric(1)),
                          count = vapply(level_by_rank, length, numeric(1)),
                          mean = vapply(level_by_rank, mean, numeric(1)))
  # Identify overlapping levels bewteen ranks ------------------------------------------------------
  overlapping <- vapply(1:nrow(rank_data),
                        function(i) any(rank_data$min[i] <= rank_data$max[-i] & rank_data$min[-i] <= rank_data$max[i]),
                        logical(1))
  if (strict && any(overlapping)) stop("Inconsistent rank levels. Use strict == FALSE or check input data.")
  if (length(unique(rank_data$mean)) != length(rank_data$mean)) warning("Order of ranks might be partially arbitrary")
  # construct ordered factor for ranks -------------------------------------------------------------
  rank_data <- rank_data[order(rank_data$mean), ]
  factor(rank_data$rank, levels = rank_data$rank, ordered = TRUE)
}


#===================================================================================================
#' Splits a taxonomy at a specific level or rank
#' 
#' Breaks one taxonomy into multiple, each with a root of a specified distance from the root.
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon.
#' @param level (\code{character} or \code{numeric} of length 1)
#' @param rank (\code{character}) The rank designation (e.g. "genus") corresponding to each item in
#' 
#' @return a \code{list} of taxon id \code{character} vectors. 
#' \code{taxa}.
#' 
#' @keywords internal
split_by_level <- function(taxa, parents, level, rank = NULL) {
  class_data <- get_class_from_el(taxa, parents)
  data <- data.frame(taxa = taxa, parents = parents, 
                     level = vapply(class_data, length, numeric(1)))
  if (is.null(rank)) {
    new_roots <- data$taxa[data$level == level]
  } else {
    new_roots <- data$taxa[rank == level]
  }
  get_children <- function(id) {
    index <- vapply(class_data, function(x) id %in% x, logical(1))
    names(class_data[index])
  }
  setNames(lapply(new_roots, get_children), new_roots)
}


#' Get taxonomy levels
#' 
#' Return An ordered factor of taxonomy levels, such as "Subkingdom" and "Order", in order of the
#'   hierarchy.
#'   
#' @keywords internal
get_taxonomy_levels <- function() {
  unique_levels <- unique(sapply(strsplit(taxize::rank_ref$ranks, ","), `[`, 1))
  unique_levels <- tolower(unique_levels)
  output <- sort(factor(unique_levels, labels = unique_levels, ordered = TRUE))
  return(output)
}