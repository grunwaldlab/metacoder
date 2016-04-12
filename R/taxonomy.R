#===================================================================================================
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

#===================================================================================================
#' Format taxize output to string
#' 
#' @param taxa Output from taxize::classification
#' @return A character vecotr of taxonomy strings pasted together with delimiters.
#' 
#' @keywords internal
format_taxize <- function(taxa) {
  output <- vapply(taxa, function(x) paste(paste(x$rank, gsub(' ', '_', x$name), sep = "__"), collapse = ";"), character(1))
  paste("Superkingdom__Life;", output, sep="")
}


#===================================================================================================
#' filter_taxonomy_string
#' 
#' @keywords internal
filter_taxonomy_string <- function(taxon, min_level, max_level, taxon_levels) {
  parsed_taxonomy <- sapply(unlist(strsplit(taxon, split=';', fixed=T)),
                            strsplit, split='__', fixed=T)
  filter <- sapply(parsed_taxonomy, function(x) ordered(x[1], taxon_levels) >= min_level & ordered(x[1], taxon_levels) <= max_level)
  parsed_taxonomy <- parsed_taxonomy[filter]
  paste(sapply(parsed_taxonomy, paste, collapse='__'), collapse=';')
}

#===================================================================================================
#' subsample_by_taxonomy
#' 
#' subsample_by_taxonomy
#' 
#' @keywords internal
subsample_by_taxonomy <- function(distance_matrix, taxon, taxon_level, level_order, triangular=TRUE, level_to_analyze = 'subtaxon', max_subset=NA) {
  base_level <- offset_ordered_factor(taxon_level, 1)
  if (level_to_analyze == 'subtaxon') {
    level_to_analyze <- base_level
  }
  
  #if the level is not applicable return NA
  if (is.na(level_to_analyze) || taxon_level >= level_to_analyze) {
    return(NA)
  }
  
  #get indexes where taxon is present 
  indexes <- grep(taxon, row.names(distance_matrix), value = FALSE, fixed = TRUE)
  if (!is.na(max_subset) && length(indexes) > max_subset) {  #if their are too many instances, randomly subsample
    indexes = sample(indexes, max_subset)
  }
  
  #subsample matrix
  submatrix <- distance_matrix[indexes, indexes, drop = FALSE]
  names <- row.names(submatrix)
  names <- mapply(FUN=filter_taxonomy_string, names, MoreArgs=list(base_level, level_to_analyze, level_order))
  row.names(submatrix) <- names
  colnames(submatrix) <- names
  if (triangular) {
    submatrix[upper.tri(submatrix, diag=TRUE)] <- NA
  }
  return(submatrix)
}

#===================================================================================================
#' taxon_info
#' 
#' taxon_info
#' 
#' @keywords internal
taxon_info <- function(identifications, level_order, separator=';') {
  split_taxonomy <- strsplit(identifications, separator, fixed=TRUE)
  taxonomy <- unlist(lapply(split_taxonomy, function(x) sapply(seq(1, length(x)), function(y) paste(x[1:y], collapse=separator))))
  counts <- table(taxonomy)
  taxon_names <- names(counts)
  counts <- as.vector(counts)
  taxon_level <- sapply(strsplit(taxon_names, separator, fixed=TRUE), 
                        function(x) level_order[max(match(sub("__.*$", "", x), level_order))])
  taxon_level <- ordered(taxon_level, level_order)
  taxon_short_names <- sapply(1:length(taxon_names), function(i) 
    filter_taxonomy_string(taxon_names[i], taxon_level[i], taxon_level[i], level_order))
  taxon_short_names <- sub("^.*__", "", taxon_short_names)
  data.frame(row.names=taxon_names, 
             level=taxon_level, 
             name=taxon_short_names, 
             count=counts)
  
}

#===================================================================================================
#' Parse taxonomic classification strings
#' 
#' Parses taxonomic classification strings into a list of \code{data.frame}
#' Can extract taxonomic ranks if present.
#' 
#' @param lineage (\code{character}) Taxonomic classifications in the form of delimited strings.
#' @param taxon_sep (\code{character; length = 1}) Character that deliminates individual taxa 
#' and potentially their ranks
#' @param rank_sep (\code{character; length = 1}) Character that deliminates a taxon and its rank.
#' @param rev_taxon If \code{TRUE}, the rank order of taxa read in a lineage is reversed to be
#' specific to broad.
#' @param rev_rank If \code{TRUE}, the rank information come after the taxon information.
#' @param taxon_col_name (\code{character}) The name of the taxon column in the \code{data.frame}s 
#' representing each classification. 
#' @param rank_col_name (\code{character}) The name of the rank column in the \code{data.frame}s 
#' representing each classification.
#' 
#' 
#' @keywords internal
parse_lineage <- function(lineage, taxon_sep, rank_sep, rev_taxon = FALSE, rev_rank = FALSE,
                          taxon_col_name = "taxon", rank_col_name = "rank_name") {
  taxa <- strsplit(lineage, split = taxon_sep, fixed = TRUE)
  if (taxon_sep == rank_sep) {
    odd_taxa <- which(vapply(taxa, function(x) length(x) %% 2 == 1, logical(1)))
    taxa[odd_taxa] <- lapply(taxa[odd_taxa], function(x) c(x, NA))
    group_by_2 <- function(x) lapply(seq(1, length(x), 2), function(i) c(x[[i]], x[[i + 1]]))
    lineage <- lapply(taxa, group_by_2)
  } else {
    lineage <- lapply(taxa, strsplit, split = rank_sep, fixed = TRUE)
  }
  if (rev_taxon) lineage <- rev(lineage)
  #   if (taxon_sep == rank_sep) {
  #     lapply(lineage, function(x) lapply(seq(2, length(x), 2), function(i) unlist(x[i:(i+1)])))
  #   }
  if (rev_rank) lineage <- lapply(lineage, function(x) lapply(x, rev))
  lineage <- lapply(lineage, function(x) plyr::ldply(x))
  if (length(unique(vapply(lineage, ncol, numeric(1)))) > 1) stop("Inconsistent lineage.")
  if (ncol(lineage[[1]]) == 1) col_names <- taxon_col_name else if
  (ncol(lineage[[1]]) == 2) col_names <- c(rank_col_name, taxon_col_name) else 
    stop("Error parsing lineage.")
  lineage <- lapply(lineage, setNames, nm = col_names)
  return(lineage)
}


#===================================================================================================
#' @keywords internal
extract_last <- function(classifications, column) {
  sapply(classifications, function(x) x[nrow(x), column])
}

#===================================================================================================
#' Indentify unique taxa in classifications
#' 
#' Indentifies unique taxa in a list of taxonomic classifications.
#' 
#' @param classifications (\code{list} of \code{data.frame})
#' Each classification must be a \code{data.frame} with a column named "taxon".
#' @param id_column (\code{character}) The column name in each classification \code{data.frame}
#' that will be used to name the elements of the returned list.
#' 
#' @return (\code{list} of \code{data.frame}) Taxonomic classifications of every unique taxon in 
#' the list of classifications given. 
#' 
#' @keywords internal
unique_taxa <- function(classifications, id_column = NULL) {
  split_classification <- function(a_classification) {
    if (!is.data.frame(a_classification) && is.na(a_classification)) return(NA)
    lapply(1:nrow(a_classification), function(i) a_classification[1:i, ])
  }
  output <- unlist(lapply(classifications[!is.na(classifications)], split_classification), recursive = FALSE)
  output <- unique(output)
  if (is.null(id_column)) {
    names(output) <- seq_along(output)
  } else {
    names(output) <- extract_last(output, column = id_column)
  }
  return(output)
}

#===================================================================================================
#' Make an adjacency list from classifications
#' 
#' Makes an adjacency list from a list of taxonomic classifications meant to represent the tree 
#' structure shared by the classifications given. 
#' 
#' @param classifications (\code{list} of \code{data.frame}) Taxnomic classifications used to build
#' adjacency list. Each classification must be a \code{data.frame} with a column named "taxon".
#' 
#' @keywords internal
taxonomy_to_adj_list  <- function(classifications) {
  process_one <- function(taxon) {
    t(mapply(function(x, y) taxon$taxon[c(x,y)],
             1:(nrow(taxon) - 1),
             2:nrow(taxon)))
  }
  output <- do.call(rbind, lapply(classifications, process_one))
  output <- data.frame(unique(output))
  names(output) <- c("taxon", "subtaxon")
  return(output)
}


#===================================================================================================
#' Append each data.frame in list
#' 
#' Add a row to each \code{data.frame} in a list.
#' 
#' @param my_list  (\code{list} of \code{data.frame})
#' @param data (named \code{list}) The column content to add. The name of each element 
#' should match a column in the \code{data.frame}s in \code{my_list}. Each element should consist
#' of as many values as their are \code{data.frame}s in \code{my_list}.
#' 
#' @keywords internal
append_to_each <- function(my_list, data) {
  process_one <- function(element, index) {
    element <- rbind(element, rep(NA, length(element)))
    for (key in names(data)) {
      element[nrow(element), key] <- data[[key]][index]
    }
    return(element)
  }
  lapply(seq_along(my_list), function(i) process_one(my_list[[i]], i))
}


#===================================================================================================
#' Add IDs for classifications
#' 
#' @param classifications (\code{list} of \code{data.frame}) Taxnomic classifications. 
#' @param id_key (\code{list} of named \code{data.frame}) A list defining what classification every 
#' unique id is. If not provided, unique IDs will be generated from \code{classifications}.
#' @param id_col_name (\code{character}) The name of the column of unique IDs that will be added
#' to each classification.
#' 
#' @seealso unique_taxa
#' 
#' @keywords internal
add_taxon_ids <- function(classifications, id_key = NULL, id_col_name = "taxon_id") {
  if (is.null(id_key)) id_key <- unique_taxa(classifications)
  add_ids_to_one <- function(a_classification) {
    taxa_in_class <- lapply(1:nrow(a_classification), function(i) a_classification[1:i, ])
    ids <- vapply(taxa_in_class,
                  function(x) names(id_key)[sapply(id_key, identical, x)], character(1))
    a_classification[id_col_name] <- ids
    return(a_classification)
  }
  output <- lapply(classifications, add_ids_to_one)
  return(output)
}

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


#===================================================================================================
#' Make minimum taxonomy for a subset of taxa
#' 
#' Subsets a table representing every unique taxon based on a subset of taxa. 
#' This preserves supertaxa, not specified in the subset, but needed to define the internal tree 
#' structure of the taxonomy. 
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon.
#' Root taxa should have \code{NA} in this column.
#' @param subset (\code{character}) Taxon IDs to restrict the new taxonomy to.
#' 
#' @return A vector of unique taxon IDs.
#' 
#' @keywords internal
restrict_taxonomy <- function(taxa, parents, subset) {
  subset <- unique(subset)
  tax_class <- get_class_from_el(taxa, parents)
  unique(unlist(tax_class[subset]))
}


#===================================================================================================
#' Count occurrences throughout a taxonomy
#' 
#' Get counts of items at each unique taxon in a taxonomy, including internal taxa.
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon.
#' @param items (\code{character}) Taxon IDs for a set of items to count.
#' 
#' @return A named vector of counts corresponding to the taxa in \code{taxa}
#' 
#' @keywords internal
get_taxon_count <- function(taxa, parents, items) {
  tax_class <- get_class_from_el(taxa, parents)
  counts <- table(unlist(tax_class[items]))
  counts <- counts[taxa]
  counts[is.na(counts)] <- 0
  setNames(as.numeric(counts), taxa)
}


#===================================================================================================
#' Identify shared taxonomy root
#' 
#' Given an adjacency/edge list, return the indexes of taxa shared by all subtaxa.
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon. A
#'   lack of a parent should be coded as \code{NA}.
#' @param include_root (\code{logical}) If \code{TRUE}, the index of the new root (taxon with the 
#'   rank that still contains all other taxa). 
#' 
#' @return A interger vector of indexes corresponding to taxa shared by all subtaxa.
#' 
#' @keywords internal
get_stem_taxa <- function(taxa, parents, include_root = FALSE) {
  parents[!(parents %in% taxa)] <- NA
  stem_indexes <- c()
  current_index <- which(is.na(parents))
  for (index in current_index) {
    while (length(index) == 1) {
      stem_indexes <- c(stem_indexes, index)
      index <- which(parents == taxa[index])
    } 
    if (!include_root && length(stem_indexes) != 0) stem_indexes <- stem_indexes[-length(stem_indexes)]
  }
  return(stem_indexes)
}


#===================================================================================================
#' Validate a taxon edge list
#' 
#' Validate a taxon edge list for use in other functions
#' 
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#' 
#' @return (\code{logical}) Returns \code{TRUE} if the list is valid.
#' 
#' @keywords internal
validate_edge_list <- function(taxa, parents) {
  if (length(taxa) == 0) stop("'taxa' has 0 length")
  if (length(parents) == 0) stop("'parents' has 0 length")
  if (length(taxa) != length(parents)) stop("'taxa' and 'parents' are of unequal length")
  if (!all(parents %in% c(taxa, NA))) stop("All 'parent' taxon IDs not found in 'taxa' taxon IDs")
  if (length(unique(taxa)) != length(taxa)) stop("All 'taxa' not unique.")
}