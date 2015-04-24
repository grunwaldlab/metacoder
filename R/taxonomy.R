#===================================================================================================
#' Get taxonomy levels
#' 
#' Return An ordered factor of taxonomy levels, such as "Subkingdom" and "Order", in order of the
#'   hierarchy.
#'   
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
#' @export
parse_lineage <- function(lineage, taxon_sep, rank_sep, rev_taxon = FALSE, rev_rank = FALSE,
                          taxon_col_name = "taxon", rank_col_name = "rank") {
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
extract_last <- function(classifications, column) {
  vapply(classifications, function(x) x[nrow(x), column], character(1))
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
#' @export
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
#' Append each data.frame in list
#' 
#' Add a row to each \code{data.frame} in a list.
#' 
#' @param my_list  (\code{list} of \code{data.frame})
#' @param data (named \code{list}) The column content to add. The name of each element 
#' should match a column in the \code{data.frame}s in \code{my_list}. Each element should consist
#' of as many values as their are \code{data.frame}s in \code{my_list}.
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
#' @param key (\code{character}) The identity of the capturing groups defined using \code{regex}.
#'  The length of \code{key} must be equal to the number of capturing groups specified in \code{regex}.
#'  Any names added to the terms will be used as column names in the output.
#'  Each term must be one of those decribed below:
#'  \describe{
#'    \item{\code{taxon_id}}{A unique numeric id for a taxon for a particular \code{database}}
#'    \item{\code{taxon_name}}{The name of a taxon. Not necessarily unique, but are specific (i.e. interperable)
#'    to a particular \code{database}.}
#'    \item{\code{taxon_info}}{Arbitrary taxon info you want included in the output. Can be used more than once.}
#'    \item{\code{class_id}}{A list of taxa unique IDs that consitute the full taxonomic classification
#'  from broad to specific (see \code{class_tax_rev}) for a particular \code{database}. Individual taxa
#'  are separated by the \code{class_tax_sep} argument and the taxon-rank group is separated by the
#'  \code{class_rank_sep} argument.}
#'    \item{\code{class_name}}{A list of taxa names that consitute the full taxonomic
#'  classification from broad to specific. Same usage as \code{class_id}}.
#'  Individual names are not necessarily unique, but are specific (i.e. interperable)
#'  to a particular \code{database}.
#'    \item{\code{item_id}}{An unique item (e.g. sequence) identifier. The taxonomy information will be
#'  looked up if available. Requires an internet connection.}
#'    \item{\code{item_name}}{An item (e.g. sequence) name. Not necessarily unique.}
#'    \item{\code{item_info}}{Arbitrary item info you want included in the output. Can be used more than once.}
#'  }
#' @param class_tax_sep (\code{character; length == 1}) Used with the \code{class_name} term in the \code{key}
#' argument. The characters used to separate individual taxa within a lineage.
#' @param class_rank_sep (\code{character; length == 1}) Used with the \code{class_name} term in the \code{key}
#' argument when a lineage contiains both taxon and rank information. This is the characters used to separate
#' th rank and the taxon name within an individual taxa in a lineage.
#' @param class_tax_rev Used with the \code{class_name} term in the \code{key} argument.
#' If TRUE, the rank order of taxa read in a lineage is reversed to be specific to broad.
#' @param class_rank_rev Used with the \code{class_name} term in the \code{key} argument  when a lineage
#' contiains both taxon and rank information. If TRUE, the rank information come after the taxon information.
#' @param taxon_in_lineage If \code{TRUE}, the lineage string included the taxon itself as its most specific
#' classification.
#' @param database (\code{character; length == 1}): The name of the database that patterns given in 
#'  \code{parser} will apply to. Valid databases include "ncbi", "itis", "eol", "col", "tropicos",
#'  and "nbn".
#' @param arbitrary_ids (\code{character} of length 1) Determines how the generation of arbitrary IDs is
#'  handled. Possible options are:
#'  \describe{
#'    \item{\code{"allow"}}{Arbitrary IDs are automatically generated if needed. These can occur intermixed
#'    with offical database IDs in the case of failed database lookups.}
#'    \item{\code{"warn"}}{Like \code{"allow"} but issue a warning when arbitrary IDs are used.}
#'    \item{\code{"error"}}{Cause an error if arbitrary IDs are needed.}
#'    \item{\code{"na"}}{Put \code{NA}s where arbitrary are needed.}
#'    \item{\code{"none"}}{Do not use a database to look up information.}
#'  } 
#' @return Returns a list of two elements:
#'  \describe{
#'    \item{\code{taxonomy}}{A list of \code{data.frame}s containing the classification of each
#'    unique taxon. The order of the elements corresponds to the rows in the "taxa" 
#'    \code{data.frame} described below.}
#'    \item{\code{taxa}}{A data.frame with one row per taxon, containing available information
#'    on each taxon. The number and nature of columns depend on the input data. Typically, a column
#'    of taxon names is present.}
#'    \item{\code{items}}{A data.frame with one row per input item (typically sequences)
#'    containing their taxon IDs and other informtion, depending on input.}
#'    }
#'    
#' @export
#===================================================================================================
extract_taxonomy <- function(input, regex, key, class_tax_sep = ";", class_rank_sep = "__", 
                             class_tax_rev = FALSE, class_rank_rev = FALSE,
                             taxon_in_lineage = TRUE, database = 'ncbi', arbitrary_ids = "warn") {
  # Constants --------------------------------------------------------------------------------------
  valid_databases <- c("ncbi", "itis", "eol", "col", "tropicos", "nbn", "none")
  valid_keys <- c("taxon_id", "taxon_name", "taxon_info", "class_id", "class_name", 
                  "item_id", "item_name", "item_info")
  valid_arb_id_opts <- c("allow", "warn", "error", "na")
  database_id_classes <- c(ncbi = "uid", itis = "tsn", eol = "eolid", col = "colid",
                           tropicos = "tpsid", nbn = "nbnid")
  id_from_name_funcs <- list(ncbi = taxize::get_uid, itis = taxize::get_tsn, eol = taxize::get_eolid,
                             col = taxize::get_colid, tropicos = taxize::get_tpsid, nbn = taxize::get_nbnid, 
                             none = NA)
  taxid_from_seqid_funcs <- list(ncbi = taxize::genbank2uid, none = NA)
  taxon_in_lineage = TRUE
  # Argument validation ----------------------------------------------------------------------------
  if (!all(key %in% valid_keys)) stop("Invalid key term. Look at documentation for valid terms.")
  database <- match.arg(tolower(database), choices = valid_databases)
  arbitrary_ids <- match.arg(tolower(arbitrary_ids), choices = valid_arb_id_opts)
  # Apply option default ---------------------------------------------------------------------------
  if (is.null(names(key))) names(key) <- key
  # Parse arguments --------------------------------------------------------------------------------
  id_class <- database_id_classes[database]
  id_from_name <- id_from_name_funcs[[database]]
  taxid_from_seqid <- taxid_from_seqid_funcs[[database]]
  # Parse input using regex ------------------------------------------------------------------------
  item_data <- data.frame(stringr::str_match(input, regex), stringsAsFactors = FALSE)
  if (ncol(item_data) != length(key) + 1) stop("The number of capture groups and keys do not match.")
  if (any(is.na(item_data))) stop("Could not parse one or more entries. Check that `regex` matches all of `input`.")
  names(item_data) <- c("input", key)
  # Get taxon id -----------------------------------------------------------------------------------
  report_found <- function(get_id_result) {
    not_found <- sum(attr(get_id_result, "match") ==  "not found")
    if (not_found > 0) warning(paste0("Could not find taxon IDs for ", not_found, " items."))    
  }
  arbitrary_taxon_ids <- FALSE
  if ("taxon_id" %in% names(item_data)) {
    if (database != "none") class(item_data$taxon_id) <- id_class
  } else if ("class_id" %in% names(item_data) && taxon_in_lineage) {
    item_classification <- parse_lineage(item_data$class_id, taxon_sep = class_tax_sep,
                                         rank_sep = class_rank_sep, rev_taxon = class_tax_rev,
                                         rev_rank = class_rank_rev, taxon_col_name = "taxon_id")
    item_data$taxon_id <- extract_last(item_classification, "taxon_id")
    if (database != "none") class(item_data$taxon_id) <- id_class
  } else if ("item_id" %in% names(item_data) && database != "none") {
    if (is.null(taxid_from_seqid)) stop("Cannot look up taxonomy from sequence id using current database.")
    item_data$taxon_id <- taxid_from_seqid(item_data$item_id)
  } else if ("taxon_name" %in% names(item_data) && database != "none") {
    item_data$taxon_id <- map_unique(item_data$taxon_name, id_from_name)
    report_found(item_data$taxon_id)
  } else {
    warning("Insufficient information supplied to infer taxon IDs. Assigning arbitrary IDs.")
    arbitrary_taxon_ids <- TRUE
  }
  # Get taxon id lineage ---------------------------------------------------------------------------
  get_id_from_name <- FALSE 
  if ("taxon_id" %in% names(item_data) && database != "none") {
    item_classification <- map_unique(item_data$taxon_id, taxize::classification,
                                      db = database, return_id = TRUE)
    item_classification[!is.na(item_classification)] <- lapply(item_classification[!is.na(item_classification)],
                                                               setNames, nm = c("name", "rank", "taxon_id"))
  } else if ("class_id" %in% names(item_data) && taxon_in_lineage) {
    item_classification <- parse_lineage(item_data$class_id, taxon_sep = class_tax_sep,
                                         rank_sep = class_rank_sep, rev_taxon = class_tax_rev,
                                         rev_rank = class_rank_rev, taxon_col_name = "taxon_id")
  } else if ("class_id" %in% names(item_data) && !taxon_in_lineage && "taxon_id" %in% names(item_data) && !arbitrary_taxon_ids) {
    item_classification <- parse_lineage(item_data$class_id, taxon_sep = class_tax_sep,
                                         rank_sep = class_rank_sep, rev_taxon = class_tax_rev,
                                         rev_rank = class_rank_rev, taxon_col_name = "taxon_id")
    item_classification <- append_to_each(item_classification, 
                                          item_data[ , c("taxon_id", "taxon_rank")])
  } else if ("class_name" %in% names(item_data) && taxon_in_lineage) {
    item_classification <- parse_lineage(item_data$class_name, taxon_sep = class_tax_sep,
                                         rank_sep = class_rank_sep, rev_taxon = class_tax_rev,
                                         rev_rank = class_rank_rev, taxon_col_name = "name")
    item_classification <- add_taxon_ids(item_classification)
    get_id_from_name <- TRUE 
  } else if ("class_name" %in% names(item_data) && !taxon_in_lineage && "taxon_name" %in% names(item_data)) {
    item_classification <- parse_lineage(item_data$class_name, taxon_sep = class_tax_sep,
                                         rank_sep = class_rank_sep, rev_taxon = class_tax_rev,
                                         rev_rank = class_rank_rev, taxon_col_name = "name")
    item_classification <- append_to_each(item_classification, 
                                          item_data[ , c("taxon_name", "taxon_rank")])
    item_classification <- add_taxon_ids(item_classification)
    get_id_from_name <- TRUE 
  } else {
    stop("Insufficient information supplied to infer lineage taxon IDs. Taxonomy structure cannot be determined.")
    item_classification <- NULL
  } 
  # Add arbitrary taxon IDs to item data if necessary ----------------------------------------------
  if (arbitrary_taxon_ids) item_data$taxon_id <- extract_last(item_classification, "taxon_id")
  # Get taxon data ---------------------------------------------------------------------------------
  taxon_id_key <- unique_taxa(item_classification, id_column = "taxon_id")
  taxon_data <- do.call(rbind, lapply(taxon_id_key, function(x) x[nrow(x), ]))
  if (database != "none") class(taxon_data$taxon_id) <- id_class
  # Get taxon id of classifications if necessary ---------------------------------------------------
  if (get_id_from_name && database != "none") {
    taxon_ids <- map_unique(taxon_data$name, id_from_name, rows = 1)
    not_unique <- names(which(table(taxon_ids) != 1))
    taxon_ids[taxon_ids %in% not_unique] <- NA
    ids_found <- !is.na(taxon_ids)
    if (arbitrary_ids == "error" && any(is.na(taxon_ids))) stop("Could not look up all taxon names. Use option `arbitrary_ids` to allow arbitrary IDs.")
    if (arbitrary_ids == "warn" && any(is.na(taxon_ids))) warning("Could not look up all taxon names. Arbitrary IDs will be applied.")
    taxon_data$taxon_id[ids_found] <- taxon_ids[ids_found]
    if (arbitrary_ids == "na") {
      taxon_data$taxon_id[!ids_found] <- NA
    } else {
      taxon_data$taxon_id[!ids_found] <- make_new_ids(sum(!ids_found), existing = taxon_data$taxon_id[ids_found]) 
    }
    taxon_data$id_type <- "arbitrary"
    taxon_data$id_type[ids_found] <- database
    names(taxon_id_key) <- taxon_data$taxon_id
    item_classification <- add_taxon_ids(item_classification, id_key = taxon_id_key)
    taxon_id_key <- add_taxon_ids(taxon_id_key, id_key = taxon_id_key)
    names(item_classification) <- extract_last(item_classification, "taxon_id")
    item_data$taxon_id <- names(item_classification)
  }
  # Get taxon name and rank if necessary -----------------------------------------------------------
  if (!arbitrary_taxon_ids && database != "none") {
    if (!("name" %in% names(taxon_data)) || !("rank" %in% names(taxon_data))) {
      result <- taxize::classification(taxon_data$taxon_id)
      result <- lapply(result, setNames, nm = c("name", "rank", "taxon_id"))
    }
    if (!("name" %in% names(taxon_data))) {
      taxon_data$name <- extract_last(result, column = "name")
      taxon_id_key <- lapply(seq_along(taxon_id_key),
                             function(i) cbind(taxon_id_key[[i]], result[[i]][ , "name", drop = FALSE]))
    }    
    if (!("rank" %in% names(taxon_data))) {
      taxon_data$rank <- extract_last(result, column = "rank")
      taxon_id_key <- lapply(seq_along(taxon_id_key),
                             function(i) cbind(taxon_id_key[[i]], result[[i]][ , "rank", drop = FALSE]))
    }
  }
  # Add taxon info ---------------------------------------------------------------------------------
  map_item_to_taxon <- function(col_index) {
    data <- lapply(taxon_data$taxon_id,
                   function(x) unique(item_data[x == item_data$taxon_id & !is.na(item_data$taxon_id), col_index]))
    data[lapply(data, length) == 0] <- NA
    if (max(vapply(data, length, numeric(1))) > 1)
      stop(paste0("taxon_info field ", col_index - 1, " content does not correspond to taxa IDs.", 
                  " Perhaps an item_info field would be more appropriate."))
    return(unlist(data))
  }
  taxon_info_cols <- 1 + which(key == "taxon_info")
  if (length(taxon_info_cols) > 0) {
    taxon_data <- cbind(taxon_data, data.frame(lapply(taxon_info_cols, map_item_to_taxon)))    
  }
  # Add arbitrary item IDs to item data if necessary -----------------------------------------------
  if (!("item_id" %in% names(item_data))) item_data$item_id <- 1:nrow(item_data)
  # Add taxon parent column ------------------------------------------------------------------------
  get_parent <- function(taxa, classifications) {
    taxon_classifications <- lapply(taxon_data$taxon_id, function(x) classifications[[x]])
    parents <- lapply(taxon_classifications, function(x) x$taxon_id[nrow(x) - 1])
    parents[vapply(parents, length, numeric(1)) == 0] <- NA
    unlist(parents)
  }
  taxon_data$parent_id <- get_parent(taxon_data$taxon_id, taxon_id_key)
  # Add counts to taxon_data -----------------------------------------------------------------------
  taxon_data$item_count <- get_taxon_count(taxon_data$taxon_id, taxon_data$parent_id, item_data$taxon_id)
  # Add taxon level column -------------------------------------------------------------------------
  get_taxon_level <- function(taxa, classifications) {
    taxon_classifications <- lapply(taxon_data$taxon_id, function(x) classifications[[x]])
    vapply(taxon_classifications, nrow, numeric(1))
  }
  taxon_data$level <- get_taxon_level(taxon_data$taxon_id, taxon_id_key)
  # Format and return output -----------------------------------------------------------------------
  names(key)[names(key) == ""] <- names(item_data)[seq_along(key) + 1][names(key) == ""] 
  names(item_data)[seq_along(key) + 1] <- names(key)
  row.names(item_data) <- item_data$item_id
  return(list(taxonomy = taxon_id_key, taxa = taxon_data, items = item_data))
}




#===================================================================================================
#' Recursivly sample items with a heirarchical classification
#' 
#' Recursivly sample a set of items with a heirarchical classification.
#' This function takes other functions as arguments and is intended to be used to make other more 
#' user-friendly functions.
#' 
#' @param root_id (\code{character} of length 1) The taxon to sample. By default, the root of the
#' taxonomy used.
#' @param get_items (\code{function(character)}) A function that returns the items assigned to the
#' a given taxon. The function's first argument should be the taxon id and it should return a data
#' structure possibly representing multiple items.
#' @param get_subtaxa (\code{function(character)}) A function that returns the sub taxa for a given
#' taxon. The function's first argument should be the taxon id and it should return a vector of
#' taxon IDs. 
#' @param get_rank (\code{function(character)}) A function that returns the rank of a given taxon
#' id. The function's first argument should be the taxon id and it should return the rank of that
#' taxon.
#' @param cat_items (\code{function(list)}) A function that takes a list of whatever is returned by
#' \code{get_items} and concatenates them into a single data structure of the type returned by
#' \code{get_items}.
#' @param max_counts (\code{numeric}) A named vector that defines that maximum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. If more than the maximum number of items exist for a given taxon, it is randomly
#' subsampled to this number. 
#' @param min_counts (\code{numeric}) A named vector that defines that minimum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param max_children (\code{numeric}) A named vector that defines that maximum number of
#' subtaxa per taxon for each level specified. The names of the vector specifies that level each
#' number applies to. If more than the maximum number of subtaxa exist for a given taxon, they
#' are randomly subsampled to this number of subtaxa. 
#' @param min_children (\code{numeric}) A named vector that defines that minimum number of
#' subtaxa in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param item_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple items and a taxon id.
#' Returns a object of the same type with some of the items potentially removed.  
#' @param subtaxa_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple subtaxa IDs and the current taxon id.
#' Returns a object of the same type with some of the subtaxa potentially removed. If a function returns
#' \code{NULL}, then no items for the current taxon are returned.
#' @param stop_conditions (\code{list} of \code{function(id)}) A list of functions that take the
#' current taxon id. If any of the functions return \code{TRUE}, the items for the current taxon are 
#' returned rather than looking for items of subtaxa, stopping the recursion.
#' @param ... Additional parameters are passed to all of the function options.
#' 
#' @seealso \code{\link{taxonomic_sample}}
#' 
#' @export
recursive_sample <- function(root_id, get_items, get_subtaxa, get_rank = NULL, cat_items = unlist,
                             max_counts = c(), min_counts = c(), max_children = c(),
                             min_children = c(), item_filters = list(), subtaxa_filters = list(),
                             stop_conditions = list(), ...) {
  # Parse options ----------------------------------------------------------------------------------
  validate_filter_options <- function(filter) {
    if (length(get(filter)) > 0 && is.null(names(get(filter)))) {
      if (!is.null(get_rank)) stop(paste0("`", filter, "` must be named if `get_rank` is defined."))
      return(setNames(get(filter), as.character(seq_along(get(filter)))))
    }
    return(get(filter))
  }
  max_counts <- validate_filter_options("max_counts")
  min_counts <- validate_filter_options("min_counts")
  max_children <- validate_filter_options("max_children")
  min_children <- validate_filter_options("min_children")
  # Set default get_rank function ------------------------------------------------------------------
  if (is.null(get_rank)) {
    get_rank <- function(id, ...) {
      return(depth)
    }
  }
  # Make max filter function factory ---------------------------------------------------------------
  max_filter_factory <- function(filter_key) {
    function(items, id, ...) {
      rank <- as.character(get_rank(id))
      if (rank %in% names(filter_key) && length(items) > filter_key[rank]) {
        items <- sample(items, filter_key[rank])
      }
      return(items)
    }
  }
  # Make min filter function factory ---------------------------------------------------------------
  min_filter_factory <- function(filter_key) {
    function(items, id, ...) {
      rank <- as.character(get_rank(id))
      if (rank %in% names(filter_key) && length(items) < filter_key[rank]) {
        items <- NULL
      }
      return(items)
    }
  }
  # Add standard filter functions ------------------------------------------------------------------
  subtaxa_filters <- c(subtaxa_filters, max_filter_factory(max_children), 
                       min_filter_factory(min_children))
  item_filters <- c(item_filters, max_filter_factory(max_counts), min_filter_factory(min_counts))
  # Recursivly sample taxon ------------------------------------------------------------------------
  recursive_part <- function(id, depth = 1, ...) {
    # Determine if to stop search  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    stop_recursion = any(vapply(stop_conditions, function(func) func(id, ...), logical(1)))
    # Get and filter subtaxa of current taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!stop_recursion) {
      sub_taxa <- get_subtaxa(id)
      for (func in subtaxa_filters) {
        sub_taxa <- func(sub_taxa, id, ...)
        if (is.null(sub_taxa)) return(NULL)
        if (length(sub_taxa) == 0) {
          stop_recursion = TRUE
          break
        }
      }
    }
    # Get items for current taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (stop_recursion) {
      items <- get_items(id, ...)
    } else {
      items <- cat_items(lapply(sub_taxa, recursive_part, depth = depth + 1))
    }
    # Filter items - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (func in item_filters) {
      items <- func(items, id, ...)
      if (is.null(items) || length(items) == 0) break
    }
    return(items)
  }
  
  recursive_part(root_id, ...)
}

#===================================================================================================
#' Recursivly sample a set of taxonomic assignments
#' 
#' Recursivly sample a set of items with taxonomic assignments and an associated taxonomy.
#' This function takes other functions as arguments that define how the taxonomy is iterpreted.
#' 
#' @param item_ids (\code{character}) Unique taxon IDs of a set of items. 
#' @param taxon_ids (\code{character}) All possible unique taxon IDs. Must have same order and
#' length as \code{parent_id}.
#' @param parent_ids (\code{character}) Unique taxon IDs of supertaxa of \code{taxon_ids}. Must have
#' same order and length as \code{taxon_ids}.
#' @param ranks (\code{character}) Ranks of \code{taxon_ids}. Must have
#' same order and length as \code{taxon_ids}.
#' Together these two arguments form an adjacency list, defining the taxonomy used to sample.
#' @param root_id (\code{character} of length 1) The taxon to sample. By default, the root of the
#' taxonomy used.
#' @param max_counts (\code{numeric}) A named vector that defines that maximum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. If more than the maximum number of items exist for a given taxon, it is randomly
#' subsampled to this number. 
#' @param min_counts (\code{numeric}) A named vector that defines that minimum number of
#' items in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param max_children (\code{numeric}) A named vector that defines that maximum number of
#' subtaxa per taxon for each level specified. The names of the vector specifies that level each
#' number applies to. If more than the maximum number of subtaxa exist for a given taxon, they
#' are randomly subsampled to this number of subtaxa. 
#' @param min_children (\code{numeric}) A named vector that defines that minimum number of
#' subtaxa in for each level specified. The names of the vector specifies that level each number
#' applies to. 
#' @param item_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple items and a taxon id.
#' Returns a object of the same type with some of the items potentially removed.  
#' @param subtaxa_filters  (\code{list} of \code{function(items, id)}) A list of functions that take a data
#' structure containing the information of multiple subtaxa IDs and the current taxon id.
#' Returns a object of the same type with some of the subtaxa potentially removed. If a function returns
#' \code{NULL}, then no items for the current taxon are returned.
#' @param stop_conditions (\code{list} of \code{function(id)}) A list of functions that take the
#' current taxon id. If any of the functions return \code{TRUE}, the items for the current taxon are 
#' returned rather than looking for items of subtaxa, stopping the recursion.
#' @param ... Additional parameters are passed to all of the function options.
#' 
#' @export
taxonomic_sample <- function(root_id, item_ids, taxon_ids, parent_ids, ranks = NULL,
                             max_counts = c(), min_counts = c(), max_children = c(),
                             min_children = c(), item_filters = list(), subtaxa_filters = list(),
                             stop_conditions = list(), ...) {
  # Define functions to interact with the taxonomic information ------------------------------------
  get_items_func <- function(id, ...) which(item_ids == id)
  get_subtaxa_func <- function(id, ...) taxon_ids[!is.na(parent_ids) & parent_ids == id]
  if (is.null(ranks)) get_rank_func <- NULL else get_rank_func <- function(id, ...) ranks[taxon_ids == id]
  # recursive sampling -----------------------------------------------------------------------------
  recursive_sample(root_id = root_id, get_items = get_items_func, get_subtaxa = get_subtaxa_func,
                   get_rank = get_rank_func, cat_items = unlist, max_counts = max_counts, 
                   min_counts = min_counts, max_children = max_children, min_children = min_children, 
                   item_filters = item_filters, subtaxa_filters = subtaxa_filters, 
                   stop_conditions = stop_conditions)
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
validate_edge_list <- function(taxa, parents) {
  if (length(taxa) == 0) stop("'taxa' has 0 length")
  if (length(parents) == 0) stop("'parents' has 0 length")
  if (length(taxa) != length(parents)) stop("'taxa' and 'parents' are of unequal length")
  if (!all(parents %in% c(taxa, NA))) stop("All 'parent' taxon IDs not found in 'taxa' taxon IDs")
  if (length(unique(taxa)) != length(taxa)) stop("All 'taxa' not unique.")
}



#===================================================================================================
#' Get all supertaxa of a taxon
#' 
#' Given one or more taxa IDs and the edge list defining the taxonomy, return the taxon IDs of all 
#' supertaxa (i.e. all taxa the target taxa are a part of).
#' 
#' @param targets (\code{character}) Taxon IDs for which supertaxa will be returned.
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the supertaxa one level above the
#'   target taxa. If \code{TRUE}, return all the supertaxa of every supertaxa, etc. 
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs 
#' @param include_target (\code{logical}) If \code{TRUE}, the target taxa are included in the output
#'  
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned 
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
get_supertaxa <- function(targets, taxa, parents, recursive = TRUE, simplify = FALSE,
                          include_target = FALSE) {
  # Argument validataion ---------------------------------------------------------------------------
  if (length(targets) == 0) stop("Argument 'targets' has 0 length")
  if (!all(targets %in% taxa)) stop("All 'targets' taxon IDs not found in 'taxa' taxon IDs")
  validate_edge_list(taxa, parents)
  parents[!(parents %in% taxa)] <- NA
  
  # Recursive function for one target --------------------------------------------------------------
  get_one <- function(target) {
    supertaxon <- parents[taxa == target]
    if (recursive) {
      if (is.na(supertaxon)) {
        return(target)
      } else {
        return(c(target, get_one(supertaxon)))
      }
    } else {
      return(c(target, supertaxon))
    }
  }
  
  # Apply function to all targets ------------------------------------------------------------------
  supertaxa <- lapply(targets, get_one)
  if (!include_target) supertaxa <- lapply(supertaxa, `[`, -1)
  if (simplify) supertaxa <- unlist(supertaxa)
  return(supertaxa)
}





#===================================================================================================
#' Get all subtaxa of a taxon
#' 
#' Given one or more taxa IDs and the edge list defining the taxonomy, return the taxon IDs of all 
#' subtaxa
#' 
#' @param targets (\code{character}) Taxon IDs for which subtaxa will be returned.
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the subtaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the subtaxa one level above the
#'   target taxa. If \code{TRUE}, return all the subtaxa of every subtaxa, etc. 
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique taxon IDs 
#' @param include_target (\code{logical}) If \code{TRUE}, the target taxa are included in the output
#'  
#' @return If \code{simplify = FALSE}, then a list of vectors of taxon IDs are returned 
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the unique taxon
#'   IDs for all \code{target} taxa are returned in a single vector.
#' 
#' @export
get_subtaxa <- function(targets, taxa, parents, recursive = TRUE, simplify = FALSE) {
  # Argument validataion ---------------------------------------------------------------------------
  if (length(targets) == 0) stop("Argument 'targets' has 0 length")
  if (!all(targets %in% taxa)) stop("All 'targets' taxon IDs not found in 'taxa' taxon IDs")
  validate_edge_list(taxa, parents)
  parents[!(parents %in% taxa)] <- NA
  
  # Get level of each input to determin processing order -------------------------------------------
  levels <- vapply(get_supertaxa(targets, taxa, parents), length, numeric(1))
  original_targets <- targets
  targets <- targets[order(levels, decreasing = TRUE)]
  
  # Recursive function for one target --------------------------------------------------------------
  get_one <- function(target, output) {
    if (target %in% names(output)) { #if this taxa has already been processed
      return(output[[target]])
    } else {
      subtaxa <- taxa[parents == target]
      subtaxa <- subtaxa[!is.na(subtaxa)]
      if (recursive) {
        if (length(subtaxa) == 0) {
          return(subtaxa)
        } else {
          subsubtaxa <- unlist(lapply(subtaxa, get_one, output = output))
          return(c(subtaxa, subsubtaxa))
        }
      } else {
        return(subtaxa)
      }
    }
  }
  
  # Apply recursive function -----------------------------------------------------------------------
  output <- list()
  for (target in targets)  {
    result <- list(get_one(target, output))
    output <- c(output, result)
  }
  names(output) <- targets
  output <- output[original_targets]
  
  # Apply `simplify` option ------------------------------------------------------------------------
  if (simplify) output <- unique(unname(unlist(output)))
  
  return(output)
}



#===================================================================================================
#' Get items associated with taxa
#' 
#' Given one or more taxa IDs and the edge list defining the taxonomy, return the items
#' (e.g. sequence information) associated with each taxon.
#' 
#' @param targets (\code{character}) Taxon IDs for which items will be returned.
#' @param taxa (\code{character}) Unique taxon IDs for every possible taxon.
#' @param parents (\code{character}) Unique taxon IDs for the subtaxa of every possible taxon. The
#'   lack of a parent should be coded as \code{NA}.
#' @param items (\code{character}) Taxon IDs for a set of items.
#' @param recursive (\code{logical}) If \code{FALSE}, only return the item assigned to the specified
#'   target taxa, not its subtaxa. If \code{TRUE}, return all the items of every subtaxa, etc. 
#' @param simplify (\code{logical}) If \code{TRUE}, then combine all the results into a single
#'   vector of unique item indexes. 
#'
#' @return If \code{simplify = FALSE}, then a list of vectors of item indexes are returned 
#'   corresponding to the \code{target} argument. If \code{simplify = TRUE}, then the item indexes
#'   for all \code{target} taxa are returned in a single vector.
#' 
#' @export
get_taxon_items <- function(targets, taxa, parents, items, recursive = TRUE, simplify = FALSE) {
  # Argument validataion ---------------------------------------------------------------------------
  if (length(targets) == 0) stop("Argument 'targets' has 0 length")
  if (!all(targets %in% taxa)) stop("All 'targets' taxon IDs not found in 'taxa' taxon IDs")
  validate_edge_list(taxa, parents)
  parents[!(parents %in% taxa)] <- NA
  
  # Get level of each input to determin processing order -------------------------------------------
  levels <- vapply(get_supertaxa(targets, taxa, parents), length, numeric(1))
  original_targets <- targets
  targets <- targets[order(levels, decreasing = TRUE)]
  
  # Recursive function for one target --------------------------------------------------------------
  get_one <- function(target, output) {
    if (target %in% names(output)) { #if this taxa has already been processed
      return(output[[target]])
    } else {
      subtaxa <- taxa[parents == target]
      subtaxa <- subtaxa[!is.na(subtaxa)]
      target_items <- which(items == target)
      if (recursive) {
        if (length(subtaxa) == 0) {
          return(target_items)
        } else {
          subtaxa_items <- unlist(lapply(subtaxa, get_one, output = output))
          return(c(target_items, subtaxa_items))
        }
      } else {
        return(target_items)
      }
    }
  }
  
  # Apply recursive function -----------------------------------------------------------------------
  output <- list()
  for (target in targets)  {
    result <- list(get_one(target, output))
    output <- c(output, result)
  }
  names(output) <- targets
  output <- output[original_targets]
  
  # Apply `simplify` option ------------------------------------------------------------------------
  if (simplify) output <- unique(unname(unlist(output)))
  
  return(output)
}




#===================================================================================================
#' Apply a function to items of each taxon
#' 
#' Calculates a value for each taxon using a function that accepts a vector of item data/IDs for
#' each taxon. Returns the a vector of results of the same length as order as the taxa used. 
#' 
#' @param taxa (\code{character}) Unique taxon IDs.
#' @param parents (\code{character}) Unique taxon IDs for the supertaxa of \code{taxa}. 
#' @param items (\code{character}) Unique taxon IDs of a set of items.
#' @param stat (\code{vector}) The input value given to the function corresponding to the
#'   \code{items} argument. If not supplied, the \code{items} taxon IDs are given to the function.
#' @param fun (\code{function}) The function to be applied. 
#' 
