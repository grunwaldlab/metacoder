#===================================================================================================
#' Get taxonomy levels
#' 
#' @return An ordered factor of taxonomy levels, such as "Subkingdom" and "Order", in order of the
#'   hierarchy.
#' @export
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
#' @export
format_taxize <- function(taxa) {
  output <- vapply(taxa, function(x) paste(paste(x$rank, gsub(' ', '_', x$name), sep = "__"), collapse = ";"), character(1))
  paste("Superkingdom__Life;", output, sep="")
}


#===================================================================================================
#' filter_taxonomy_string
filter_taxonomy_string <- function(taxon, min_level, max_level, taxon_levels) {
  parsed_taxonomy <- sapply(unlist(strsplit(taxon, split=';', fixed=T)),
                            strsplit, split='__', fixed=T)
  filter <- sapply(parsed_taxonomy, function(x) ordered(x[1], taxon_levels) >= min_level & ordered(x[1], taxon_levels) <= max_level)
  parsed_taxonomy <- parsed_taxonomy[filter]
  paste(sapply(parsed_taxonomy, paste, collapse='__'), collapse=';')
}

#===================================================================================================
#' subsample_by_taxonomy
#' @export
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
#' @export
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
#' Extract taxonomy information from sequence headers
#' 
#' Extracts the taxonomy used by a set of sequences based on their header information. A data 
#' structure representing the heirerarchical nature of the taxonomy as well as a vector
#' identifing the taxon of each sequence is returned. Taxa are translated into unique codes if they
#' are not already encoded this way.
#' 
#' @param input (\code{character}) A vector from which to extract taxonomy information. 
#' @param parser (\code{character; length == 1}) A regular expression with custom replacement
#'  patterns indicating the locations of relevant information. Supported patterns include:
#'  \describe{
#'    \item{\code{\%\%tax_id\%\%}}{A unique numeric id for a taxon.}
#'    \item{\code{\%\%name\%\%}}{The name of a taxon. Not necessarily unique.}
#'    \item{\code{\%\%lineage<sep1><sep2>\%\%}}{A list of taxa names that consitute the full taxonomic classification
#'  from broad to specific. Individual names are not necessarily unique. Individual names are separated
#'  by the \code{<sep>} terms. \code{<sep1>}
#'  is required, as it defines the characters the deliminate taxa in the lineage. \code{<sep2>} is optional, 
#'  and further splits each taxon into its name and its level. For example: 
#'  \href{http://unite.ut.ee/repository.php}{UNITE} uses a format for taxonomy like
#'  "k__Fungi;p__Ascomycota;c__Dothideomycetes...". For this format use \code{\%\%lineage<;><__>\%\%}}
#'    \item{\code{\%\%lineage_id<sep1><sep2>\%\%}}{A list of taxa unique ids that consitute the full taxonomic
#'  classification from broad to specific. Same usage as \code{\%\%lineage<sep1><sep2>\%\%}}
#'    \item{\code{\%\%seq_id\%\%}}{An unique sequence identifier. The taxonomy information will be
#'  looked up if available. Requires an internet connection.}
#'  }
#' @return Returns a list of two elements:
#'  \describe{
#'    \item{\code{taxonomy}}{An type of \href{http://en.wikipedia.org/wiki/Adjacency_list}{adjacency list}
#'    represented as a data.frame with 3  columns. The first column is the taxon id and the second two
#'    are ids of a taxon and its subtaxon respectivly. If a taxon has more than subtaxon, it will be
#'    present multiple times in the 2nd "parent taxon" column. However, a taxon can only be 
#'    present once in the 3rd "sub taxon" column; otherwise, it would imply the taxon has more than
#'    one "parent". This data structure can be converted to more intuitive ones with ... TO BE CONTINUED.}
#'    \item{\code{taxon_info}}{A data.frame with one row per taxon, containing available information
#'    on each taxon. The number and nature of columns depend on the input data. Typically, a column
#'    of taxon names is present. }
#'    \item{\code{item_info}}{A data.frame corresponding to input items (typically sequences)
#'    containing their taxon ids. In addition, the content of any
#'    \href{http://en.wikipedia.org/wiki/Regular_expression}{regex} capturing groups (i.e. parts bound
#'    by parenthesis) in the \code{parser} argument will be present as columns of type
#'    \code{character}.}
#'  }
#' @param database (\code{character}): The name of the database that patterns given in 
#'  \code{parser} will apply to. Currently, only \code{ncbi} is being supported.
#' @export
#===================================================================================================
extract_taxonomy <- function(input, parser, database = 'ncbi') {
  replacement_delim <- "%%"
  # Parse parsing regex ----------------------------------------------------------------------------
  replacments <- strsplit("%%lineage<;><__>%%|%%seq_id%%", replacement_delim)[[1]]
  replacments <- replacments[seq(2, length(replacments), by = 2)]
  # Get taxon id -----------------------------------------------------------------------------------
  
  # Get taxon lineage ------------------------------------------------------------------------------
}