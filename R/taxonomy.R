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
#' @param regex (\code{character; length == 1}) A regular expression with capturing groups
#'  indicating the locations of relevant information. The identity of the information must
#'  be specified using the \code{key} argument.
#' @param key (\code{character}) The identity of the caturing groups defined using \code{regex}.
#'  The length of \code{key} must be equal to the number of capturing groups specified in \code{regex}.
#'  Any names added to the terms will be used as column names in the output.
#'  Each term must be one of those decribed below:
#'  \describe{
#'    \item{\code{tax_name}}{The name of a taxon. Not necessarily unique.}
#'    \item{\code{tax_id}}{A unique numeric id for a taxon.}
#'    \item{\code{tax_rank}}{A taxonomic rank name (e.g. "genus").}
#'    \item{\code{tax_info}}{Arbitrary taxon info you want included in the output. Can be used more than once.}
#'    \item{\code{lineage}}{A list of taxa names that consitute the full taxonomic classification
#'  from broad to specific (see \code{lineage.tax.rev}). Individual names are not necessarily unique.
#'  Individual taxa are separated by the \code{lineage.tax.sep} argument and the taxon-rank group is separated
#'  by the \code{lineage.rank.sep} argument.}
#'    \item{\code{lineage_id}}{A list of taxa unique ids that consitute the full taxonomic
#'  classification from broad to specific. Same usage as \code{lineage}}.
#'    \item{\code{item_name}}{An item (e.g. sequence) name. Not necessarily unique.}
#'    \item{\code{item_id}}{An unique item (e.g. sequence) identifier. The taxonomy information will be
#'  looked up if available. Requires an internet connection.}
#'    \item{\code{item_info}}{Arbitrary item info you want included in the output. Can be used more than once.}
#'  }
#' @param lineage.tax.sep (\code{character; length == 1}) Used with the \code{lineage} term in the \code{key}
#' argument. The characters used to separate individual taxa within a lineage.
#' @param lineage.rank.sep (\code{character; length == 1}) Used with the \code{lineage} term in the \code{key}
#' argument when a lineage contiains both taxon and rank information. This is the characters used to separate
#' th rank and the taxon name within an individual taxa in a lineage.
#' @param lineage.tax.rev Used with the \code{lineage} term in the \code{key} argument.
#' If TRUE, the rank order of taxa read in a lineage is reversed to be specific to broad.
#' @param lineage.rank.rev Used with the \code{lineage} term in the \code{key} argument  when a lineage
#' contiains both taxon and rank information. If TRUE, the rank information come after the taxon information.
#' @param database (\code{character; length == 1}): The name of the database that patterns given in 
#'  \code{parser} will apply to. Currently, only \code{ncbi} is being supported.
#' @return Returns a list of two elements:
#'  \describe{
#'    \item{\code{taxonomy}}{A data.frame with 2 columns representing A type of
#'    \href{http://en.wikipedia.org/wiki/Adjacency_list}{adjacency list}.
#'    The columns are the ids of a taxon and its subtaxon
#'    respectivly. If a taxon has more than subtaxon, it will be
#'    present multiple times in the first "parent taxon" column. However, a taxon can only be 
#'    present once in the second "sub taxon" column; otherwise, it would imply the taxon has more than
#'    one "parent". This data structure can be converted to more intuitive ones with ... TO BE CONTINUED.}
#'    \item{\code{taxa}}{A data.frame with one row per taxon, containing available information
#'    on each taxon. The number and nature of columns depend on the input data. Typically, a column
#'    of taxon names is present.}
#'    \item{\code{items}}{A data.frame with one row per input item (typically sequences)
#'    containing their taxon ids and other informtion, depending on input.}
#'    }
#' @export
#===================================================================================================
extract_taxonomy <- function(input, regex, key, lineage.tax.sep = ";", lineage.rank.sep = "__", 
                             lineage.tax.rev = FALSE, lineage.rank.rev = FALSE, database = 'ncbi') {
  replacement_delim <- "%%"
  # Parse parsing regex ----------------------------------------------------------------------------
  replacments <- strsplit("%%lineage<;><__>%%|%%seq_id%%", replacement_delim)[[1]]
  replacments <- replacments[seq(2, length(replacments), by = 2)]
  # Get taxon id -----------------------------------------------------------------------------------
  
  # Get taxon lineage ------------------------------------------------------------------------------
}