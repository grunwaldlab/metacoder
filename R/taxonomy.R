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
#' 
#' @export
parse_lineage <- function(lineage, taxon_sep, rank_sep, rev_taxon, rev_rank) {
  taxa <- strsplit(lineage, split = taxon_sep)
  if (rev_taxon) taxa <- rev(taxa)
  lineage <- lapply(taxa, strsplit, split = rank_sep, fixed = TRUE)
  if (rev_rank) lineage <- lapply(lineage, function(x) lapply(x, rev))
  lineage <- lapply(lineage, function(x) plyr::ldply(x))
  if (length(unique(vapply(lineage, ncol, numeric(1)))) > 1) stop("Inconsistent lineage.")
  if (ncol(lineage[[1]]) == 1) col_names <- "taxon" else if
  (ncol(lineage[[1]]) == 2) col_names <- c("rank", "taxon") else 
    stop("Error parsing lineage.")
  lineage <- lapply(lineage, setNames, nm = col_names)
  return(lineage)
}


#===================================================================================================
extract_most_specific <- function(classifications) {
  vapply(classifications, function(x) x[nrow(x), "taxon"], character(1))
}

#===================================================================================================
#' Indentify unique taxa in classifications
#' 
#' Indentifies unique taxa in a list of taxonomic classifications.
#' 
#' @param classifications (\code{list} of \code{data.frame})
#' Each classification must be a \code{data.frame} with a column named "taxon".
#' 
#' @return (\code{list} of \code{data.frame}) Taxonomic classifications of every unique taxon in 
#' the list of classifications given. 
#' 
#' @export
unique_taxa <- function(classifications) {
  split_classification <- function(a_classification) {
    lapply(1:nrow(a_classification), function(i) a_classification[1:i, ])
  }
  output <- unlist(lapply(classifications, split_classification), recursive = FALSE)
  output <- unique(output)
  if ("taxon" %in% names(output[[1]])) names(output) <- extract_most_specific(output)
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
#' @export
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
add_to_classification <- function(classifications, taxa, ranks = NULL) {
  process_one <- function(x, taxon, rank) {
    x <- rbind(x, rep(NA, ncol(x)))
    x[nrow(x), "taxon"] <- taxon
    if (!is.null(rank)) x[nrow(x), "rank"] <- rank
    return(x)
  }
  lapply(classifications, process_one)
}

#===================================================================================================
#' Generate unique ids for classifications
#' 
#' Generates unique ids for every taxon in a list of classifications.
#' 
#' @param classifications (\code{list} of \code{data.frame}) Taxnomic classifications. 
#' 
#' @export
make_unique_ids <- function(classification) {
  
}

#===================================================================================================
#' Add ids for classifications
#' 
#' @param classifications (\code{list} of \code{data.frame}) Taxnomic classifications. 
#' @param ids (\code{list} of named \code{data.frame}) A list defining what classification every 
#' unique id is. If not provided, unique ids will be generated from \code{classifications}.
#' 
#' @export
add_taxon_ids <- function(classification, ids = NULL) {
  
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
#'    \item{\code{taxon_name}}{The name of a taxon. Not necessarily unique.}
#'    \item{\code{taxon_id}}{A unique numeric id for a taxon.}
#'    \item{\code{taxon_rank}}{A taxonomic rank name (e.g. "genus").}
#'    \item{\code{taxon_info}}{Arbitrary taxon info you want included in the output. Can be used more than once.}
#'    \item{\code{lineage}}{A list of taxa names that consitute the full taxonomic classification
#'  from broad to specific (see \code{lineage_tax_rev}). Individual names are not necessarily unique.
#'  Individual taxa are separated by the \code{lineage_tax_sep} argument and the taxon-rank group is separated
#'  by the \code{lineage_rank_sep} argument.}
#'    \item{\code{lineage_id}}{A list of taxa unique ids that consitute the full taxonomic
#'  classification from broad to specific. Same usage as \code{lineage}}.
#'    \item{\code{item_name}}{An item (e.g. sequence) name. Not necessarily unique.}
#'    \item{\code{item_id}}{An unique item (e.g. sequence) identifier. The taxonomy information will be
#'  looked up if available. Requires an internet connection.}
#'    \item{\code{item_info}}{Arbitrary item info you want included in the output. Can be used more than once.}
#'  }
#' @param lineage_tax_sep (\code{character; length == 1}) Used with the \code{lineage} term in the \code{key}
#' argument. The characters used to separate individual taxa within a lineage.
#' @param lineage_rank_sep (\code{character; length == 1}) Used with the \code{lineage} term in the \code{key}
#' argument when a lineage contiains both taxon and rank information. This is the characters used to separate
#' th rank and the taxon name within an individual taxa in a lineage.
#' @param lineage_tax_rev Used with the \code{lineage} term in the \code{key} argument.
#' If TRUE, the rank order of taxa read in a lineage is reversed to be specific to broad.
#' @param lineage_rank_rev Used with the \code{lineage} term in the \code{key} argument  when a lineage
#' contiains both taxon and rank information. If TRUE, the rank information come after the taxon information.
#' @param taxon_in_lineage If \code{TRUE}, the lineage string included the taxon itself as its most specific
#' classification.
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
extract_taxonomy <- function(input, regex, key, lineage_tax_sep = ";", lineage_rank_sep = "__", 
                             lineage_tax_rev = FALSE, lineage_rank_rev = FALSE,
                             taxon_in_lineage = TRUE, database = 'ncbi') {
  unique_mapping <- function(input) {
    unique_input <- unique(input)
    vapply(input, function(x) which(x == unique_input), numeric(1))
  }
  map_unique <- function(input, func) {
    func(unique(input))[unique_mapping(input)]
  }
  valid_databases <- c("ncbi", "itis", "eol", "col", "tropicos", "nbn")
  valid_keys <- c("taxon_name", "taxon_id", "taxon_rank", "taxon_info", "lineage", "lineage_id",
                  "item_name", "item_id", "item_info")
  database_id_classes <- c(ncbi = "uid", itis = "tsn", eol = "eolid", col = "colid",
                           tropicos = "tpsid", nbn = "nbnid")
  get_functions <- c(ncbi = taxize::get_uid, itis = taxize::get_tsn, eol = taxize::get_eolid,
                     col = taxize::get_colid, tropicos = taxize::get_tpsid, nbn = taxize::get_nbnid)
  taxon_in_lineage = TRUE
  # Argument validation ----------------------------------------------------------------------------
  if (!all(key %in% valid_keys)) stop("Invalid key term. Look at documentation for valid terms.")
  database <- match.arg(database, choices = valid_databases)
  # Argument parseing ------------------------------------------------------------------------------
  id_class <- database_id_classes[database]
  get_id_func <- get_functions[database]
  # Parse input using regex ------------------------------------------------------------------------
  item_data <- data.frame(stringr::str_match(input, regex))
  names(item_data) <- c("input", key)
  if (ncol(item_data) != length(key)) stop("The number of capture groups and keys do not match.")
  # Get taxon id -----------------------------------------------------------------------------------
  report_found <- function(get_id_result) {
    not_found <- sum(attr(get_id_result, "match") ==  "not found")
    if (not_found > 0) warning(paste0("Could not find taxon ids for ", not_found, " items.")    
  }
  arbitrary_taxon_ids <- FALSE
  if ("taxon_id" %in% names(item_data)) {
    class(item_data$taxon_id) <- id_class
  } else {
    if ("lineage_id" %in% names(item_data) && taxon_in_lineage) {
      item_classification <- parse_lineage(item_data$lineage_id, taxon_sep = lineage_tax_sep,
                                           rank_sep = lineage_rank_sep, rev_taxon = lineage_tax_rev,
                                           rev_rank = lineage_rank_rev)
      item_data$taxon_id <- extract_most_specific(item_classification)
      class(item_data$taxon_id) <- id_class
    } else if ("taxon_name" %in% names(item_data)) {
      item_data$taxon_id <- map_unique(item_data$taxon_name, get_id_func)
      report_found(item_data$taxon_id)
    } else if ("lineage" %in% names(item_data) && taxon_in_lineage) {
      item_classification <- parse_lineage(item_data$lineage, taxon_sep = lineage_tax_sep,
                                           rank_sep = lineage_rank_sep, rev_taxon = lineage_tax_rev,
                                           rev_rank = lineage_rank_rev)
      item_data$taxon_id <- map_unique(extract_most_specific(item_classification), get_id_func)
      report_found(item_data$taxon_id)
    } else {
      warning("Insufficient information supplied to infer taxon ids. Assigning arbitrary ids.")
      arbitrary_taxon_ids <-TRUE
    }
  }
  # Get taxon id lineage ---------------------------------------------------------------------------
  if ("lineage_id" %in% names(item_data) && taxon_in_lineage) {
    item_classification <- parse_lineage(item_data$lineage_id, taxon_sep = lineage_tax_sep,
                                         rank_sep = lineage_rank_sep, rev_taxon = lineage_tax_rev,
                                         rev_rank = lineage_rank_rev)
  } else if ("lineage_id" %in% names(item_data) && !taxon_in_lineage && "taxon_id" %in% names(item_data) && !arbitrary_taxon_ids) {
    item_classification <- parse_lineage(item_data$lineage_id, taxon_sep = lineage_tax_sep,
                                         rank_sep = lineage_rank_sep, rev_taxon = lineage_tax_rev,
                                         rev_rank = lineage_rank_rev)
    item_classification <- add_to_classification(item_classification, item_data$taxon_id,
                                                 item_data$taxon_rank)
  } else if ("taxon_id" %in% names(item_data) && !arbitrary_taxon_ids) {
    item_classification <- map_unique(item_data$taxon_id, classification)
  } else if (arbitrary_taxon_ids && "lineage" %in% names(item_data) && taxon_in_lineage) {
    item_classification <- parse_lineage(item_data$lineage, taxon_sep = lineage_tax_sep,
                                         rank_sep = lineage_rank_sep, rev_taxon = lineage_tax_rev,
                                         rev_rank = lineage_rank_rev)
    item_classification <- add_unique_ids(item_classification)
  } else if (arbitrary_taxon_ids && "lineage" %in% names(item_data) && !taxon_in_lineage && "taxon_name" %in% names(item_data)) {
    item_classification <- parse_lineage(item_data$lineage, taxon_sep = lineage_tax_sep,
                                         rank_sep = lineage_rank_sep, rev_taxon = lineage_tax_rev,
                                         rev_rank = lineage_rank_rev)
    item_classification <- add_to_classification(item_classification, item_data$taxon_name,
                                                 item_data$taxon_rank)
    item_classification <- add_unique_ids(item_classification)
  } else {
    warning("Insufficient information supplied to infer lineage taxon ids. Taxonomy structure cannot be determined.")
    item_classification <- NULL
  } 
  
  # Get taxon name ---------------------------------------------------------------------------------
  if (!("taxon_name" %in% names(item_data))) {
    if ("lineage" %in% names(item_data) && taxon_in_lineage) {
      name_taxonomy <- parse_lineage(item_data$lineage)
      taxon_data$taxon_name <- vapply(name_taxonomy, function(x) x$taxon[nrow(x)], character(1))
    } else if ("taxon_id" %in% names(item_data) && !arbitrary_taxon_ids && !("lineage_id" %in% key)) {
      taxon_data$taxon_name <- vapply(taxonomy, function(x) x$name[nrow(x)], character(1))
    }
    item_data$taxon_name <-  taxon_data$taxon_name[unique_mapping(item_data$taxon_id)]
  }
  # Get taxon rank ---------------------------------------------------------------------------------
  if (!("taxon_rank" %in% names(item_data))) {
    rank_in_taxonomy <- !is.null(taxonomy) && ("rank" %in% unlist(lapply(taxonomy, names)))
    if (rank_in_taxonomy && taxon_in_lineage) {
      taxon_data$taxon_rank <- vapply(taxonomy, function(x) x$rank[nrow(x)], character(1))
    } else if ("lineage" %in% names(item_data) && taxon_in_lineage && rank_innam_taxonomy) {}
  }
  # Get item id ------------------------------------------------------------------------------------
  
  # Get item name ----------------------------------------------------------------------------------
  
  # Rename columns ---------------------------------------------------------------------------------
  
  item_data$taxon_id <- 1:nrow(data)
}