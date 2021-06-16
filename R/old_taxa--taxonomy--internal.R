#' Infer edge list from hierarchies
#'
#' Infer edge list and unique taxa from hierarchies.
#'
#' @param A list of [hierarchy()] objects.
#'
#' @keywords internal
parse_heirarchies_to_taxonomy <- function(heirarchies) {

  # Look for input edge cases
  total_taxa_count <- sum(vapply(heirarchies, function(x) length(x$taxa), numeric(1)))
  if (length(heirarchies) == 0 || total_taxa_count == 0) {
    return(list(taxa = list(), edge_list = data.frame(from = character(), to = character(),
                                                      stringsAsFactors=FALSE)))
  }

  # Convert to character vecotor and use parse_raw_heirarchies_to_taxonomy
  taxon_names <- lapply(heirarchies,
                        function(heirarchy) vapply(heirarchy$taxa,
                                                   function(taxon) paste0(taxon$get_name(),
                                                                          taxon$get_rank(),
                                                                          taxon$get_id(),
                                                                          taxon$authority),
                                                   character(1)))
  input_taxa <- unlist(lapply(heirarchies, function(x) x$taxa))
  names(input_taxa) <- unlist(taxon_names)
  output <- parse_raw_heirarchies_to_taxonomy(taxon_names)

  # Replace taxa list with objects from input
  output$taxa <- stats::setNames(input_taxa[vapply(output$taxa, function(x) x$name$name, character(1))],
                                 names(output$taxa))

  return(output)
}


#' Infer edge list from hierarchies composed of character vectors
#'
#' Infer edge list and unique taxa from hierarchies.
#'
#' @param A list of character vectors.
#' @param named_by_rank (`TRUE`/`FALSE`) If  `TRUE` and the input is a list of
#'   vectors with each vector named by ranks, include that rank info in the
#'   output object, so it can be accessed by `out$taxon_ranks()`. If `TRUE`,
#'   taxa with different ranks, but the same name and location in the taxonomy,
#'   will be considered different taxa.
#'
#' @keywords internal
parse_raw_heirarchies_to_taxonomy <- function(heirarchies, named_by_rank = FALSE) {

  # Look for input edge cases
  total_taxa_count <- sum(vapply(heirarchies, length, numeric(1)))
  if (length(heirarchies) == 0 || total_taxa_count == 0) {
    return(list(taxa = list(), edge_list = data.frame(from = character(), to = character(),
                                                      stringsAsFactors=FALSE)))
  }

  # This is used to store both taxon names and their IDs once assigned. The IDs will be added as
  # names to the taxon character vectors
  taxon_names <- lapply(heirarchies,
                        function(heirarchy) stats::setNames(c(NA, heirarchy), NA))

  # initialize output lists
  unique_taxa <- list()
  edge_list <- list() # matrix?
  max_id <- 0 # used to keep track of what IDs have already been used

  # Find the maximum depth of the classifications
  heirarchies_depths <- vapply(taxon_names, length, numeric(1))
  max_depth <- max(heirarchies_depths)

  # For each level in the classifications, identify unique taxa and add them to the taxon list and
  # the edge list. NOTE: This function modifies variables outside of it and is not independent.
  process_one_level <- function(depth) {

    # Identify unique pairs of taxon ids and unclassified taxon names and make new IDs
    if (named_by_rank) { # If names are ranks, then consider rank part of the taxon name
      all_pairs <- lapply(taxon_names[heirarchies_depths >= depth],
                          function(x) c(names(x)[depth - 1], paste0(names(x)[depth], x[depth])))
    } else {
      all_pairs <- lapply(taxon_names[heirarchies_depths >= depth],
                          function(x) c(names(x)[depth - 1], x[depth]))
    }
    unique_encoding <- match(all_pairs, unique(all_pairs))
    new_ids <- unique_encoding + max_id
    max_id <<- max(new_ids)

    # Add new IDs to `taxon_names` as vector names
    unused_output <- mapply(seq_along(taxon_names)[heirarchies_depths >= depth], new_ids,
                            FUN = function(index, id) {names(taxon_names[[index]])[depth] <<- id})

    # Get representative taxa objects to add to the taxon list.
    # The `depth - 1` is because `NA` was added to each hierachy in `taxon_names` so the indexes
    # are one off.
    taxon_objects <- lapply(heirarchies[heirarchies_depths > (depth - 1)],
                            function(heirarchy) heirarchy[depth - 1])
    taxon_objects <- taxon_objects[match(unique(unique_encoding), unique_encoding)] # make unique
    if (named_by_rank) {
      taxon_objects <- mapply(taxon, taxon_objects, vapply(taxon_objects, names, character(1)))
    } else {
      taxon_objects <- mapply(taxon, taxon_objects)
    }
    new_taxa <- stats::setNames(taxon_objects, unique(new_ids))


    # TODO: either check that unique taxa identified by name are actually unique when considering
    # all fields (e.g. `id`) or rework function to optionally consider all fields with identifiying
    # unique taxa in the first place

    # Get edge list additions
    new_edges <- unique(lapply(taxon_names[heirarchies_depths >= depth],
                               function(x) c(names(x)[c(depth - 1, depth)])))

    # Append classified taxa and edges to ouptut of parent function
    unique_taxa <<- c(unique_taxa, new_taxa)
    edge_list <<- c(edge_list, new_edges)
  }
  no_ouput <- lapply(2:max_depth, process_one_level) # starts at 2 because of NA being prepended

  # Convert edge list to matrix
  edge_list <- stats::setNames(as.data.frame(do.call(rbind, edge_list), stringsAsFactors = FALSE),
                               c("from", "to"))

  # Get input taxon ids
  input_ids <- vapply(taxon_names, function(x) names(x)[length(x)], character(1))

  return(list(taxa = unique_taxa, edge_list = edge_list, input_ids = input_ids))
}


#' List to vector of unique elements
#'
#' Implements the `simplify` option in many functions like [supertaxa()].
#' Returns unique name-value pairs if all vectors are named.
#'
#' @param input A list of vectors
#'
#' @keywords internal
simplify <- function(input) {
  if (any(vapply(input, function(x) is.null(names(x)), logical(1)))) { # Any unnamed vectors
    output <- unique(unlist(unname(lapply(input, unname))))
  } else {
    key_value <- data.frame(stringsAsFactors = FALSE,
                            my_name = unlist(lapply(input, names)), my_value = unlist(input))
    key_value <- unique(key_value)
    output <- stats::setNames(key_value$my_value, key_value$my_name)
  }
  return(output)
}
