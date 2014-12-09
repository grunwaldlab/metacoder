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
#' @param taxon_col_name (\code{character}) The name of the taxon column in the \code{data.frame}s 
#' representing each classification. 
#' @param rank_col_name (\code{character}) The name of the rank column in the \code{data.frame}s 
#' representing each classification. 
#' 
#' @export
parse_lineage <- function(lineage, taxon_sep, rank_sep, rev_taxon = FALSE, rev_rank = FALSE,
                          taxon_col_name = "taxon", rank_col_name = "rank") {
  taxa <- strsplit(lineage, split = taxon_sep)
  if (taxon_sep == rank_sep) {
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
    lapply(1:nrow(a_classification), function(i) a_classification[1:i, ])
  }
  output <- unlist(lapply(classifications, split_classification), recursive = FALSE)
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
#' Add ids for classifications
#' 
#' @param classifications (\code{list} of \code{data.frame}) Taxnomic classifications. 
#' @param id_key (\code{list} of named \code{data.frame}) A list defining what classification every 
#' unique id is. If not provided, unique ids will be generated from \code{classifications}.
#' @param id_col_name (\code{character}) The name of the column of unique ids that will be added
#' to each classification.
#' 
#' @seealso unique_taxa
#' 
#' @export
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
#' Generate new unique ids
#' 
#' Makes a vector of unique ids that differ from a previously defined set of unique ids. 
#' @param count (\code{numeric} of length 1) The number of new unique ids to generate
#' @param existing (\code{character}) Existing unique ids. These will not appear in the output.
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
#'    \item{\code{class_id}}{A list of taxa unique ids that consitute the full taxonomic classification
#'  from broad to specific (see \code{class_tax_rev}) for a particular \code{database}. Individual taxa
#'  are separated by the \code{class_tax_sep} argument and the taxon-rank group is separated by the
#'  \code{class_rank_sep} argument.}
#'    \item{\code{class_name}}{A list of taxa names that consitute the full taxonomic
#'  classification from broad to specific. Same usage as \code{class_id}}.
#'  Individual names are not necessarily unique, but are specific (i.e. interperable)
#'  to a particular \code{database}.
#'    \item{\code{item_id}}{An unique item (e.g. sequence) identifier. The taxonomy information will be
#'    \item{\code{item_name}}{An item (e.g. sequence) name. Not necessarily unique.}
#'  looked up if available. Requires an internet connection.}
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
#' @param arbitrary_ids (\code{character} of length 1) Determines how the generation of arbitrary ids is
#'  handled. Possible options are:
#'  \describe{
#'    \item{\code{"allow"}}{Arbitrary ids are automatically generated if needed. These can occur intermixed
#'    with offical database ids in the case of failed database lookups.}
#'    \item{\code{"warn"}}{Like \code{"allow"} but issue a warning when arbitrary ids are used.}
#'    \item{\code{"error"}}{Cause an error if arbitrary ids are needed.}
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
#'    containing their taxon ids and other informtion, depending on input.}
#'    }
#' @export
#===================================================================================================
extract_taxonomy <- function(input, regex, key, class_tax_sep = ";", class_rank_sep = "__", 
                             class_tax_rev = FALSE, class_rank_rev = FALSE,
                             taxon_in_lineage = TRUE, database = 'ncbi', arbitrary_ids = "warn") {
  browser()
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
  names(item_data) <- c("input", key)
  if (ncol(item_data) != length(key) + 1) stop("The number of capture groups and keys do not match.")
  if (any(is.na(item_data))) stop("Could not parse one or more entries. Check that `regex` matches all of `input`.")
  # Get taxon id -----------------------------------------------------------------------------------
  report_found <- function(get_id_result) {
    not_found <- sum(attr(get_id_result, "match") ==  "not found")
    if (not_found > 0) warning(paste0("Could not find taxon ids for ", not_found, " items."))    
  }
  arbitrary_taxon_ids <- FALSE
  if ("taxon_id" %in% names(item_data)) {
    class(item_data$taxon_id) <- id_class
  } else if ("class_id" %in% names(item_data) && taxon_in_lineage) {
    item_classification <- parse_lineage(item_data$class_id, taxon_sep = class_tax_sep,
                                         rank_sep = class_rank_sep, rev_taxon = class_tax_rev,
                                         rev_rank = class_rank_rev, taxon_col_name = "taxon_id")
    item_data$taxon_id <- extract_last(item_classification, "taxon_id")
    class(item_data$taxon_id) <- id_class
  } else if ("item_id" %in% names(item_data) && database != "none") {
    if (is.null(taxid_from_seqid)) stop("Cannot look up taxonomy from sequence id using current database.")
    item_data$taxon_id <- taxid_from_seqid(item_data$item_id)
  } else if ("taxon_name" %in% names(item_data) && database != "none") {
    item_data$taxon_id <- map_unique(item_data$taxon_name, id_from_name)
    report_found(item_data$taxon_id)
  } else {
    warning("Insufficient information supplied to infer taxon ids. Assigning arbitrary ids.")
    arbitrary_taxon_ids <- TRUE
  }
  # Get taxon id lineage ---------------------------------------------------------------------------
  get_id_from_name <- FALSE 
  if ("taxon_id" %in% names(item_data) && database != "none") {
    item_classification <- map_unique(item_data$taxon_id, taxize::classification,
                                      db = database, return_id = TRUE)
    item_classification <- lapply(item_classification, setNames, nm = c("name", "rank", "taxon_id"))
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
    stop("Insufficient information supplied to infer lineage taxon ids. Taxonomy structure cannot be determined.")
    item_classification <- NULL
  } 
  # Add arbitrary taxon ids to item data if necessary ----------------------------------------------
  if (arbitrary_taxon_ids) item_data$taxon_id <- extract_last(item_classification, "taxon_id")
  # Get taxon data ---------------------------------------------------------------------------------
  taxon_id_key <- unique_taxa(item_classification, id_column = "taxon_id")
  taxon_data <- do.call(rbind, lapply(taxon_id_key, function(x) x[nrow(x), ]))
  class(taxon_data$taxon_id) <- id_class
  # Get taxon id of classifications if necessary ---------------------------------------------------
  if (get_id_from_name && database != "none") {
    taxon_ids <- id_from_name(taxon_data$name)
    ids_found <- !is.na(taxon_ids)
    if (arbitrary_ids == "error" && any(is.na(taxon_ids))) stop("Coud not look up all taxon names. Use option `arbitrary_ids` to allow arbitrary ids.")
    if (arbitrary_ids == "warn" && any(is.na(taxon_ids))) warning("Coud not look up all taxon names. Arbitrary ids will be applied.")
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
      taxon_data$name <- extract_last(result, column = "rank")
      taxon_id_key <- lapply(seq_along(taxon_id_key),
                             function(i) cbind(taxon_id_key[[i]], result[[i]][ , "rank", drop = FALSE]))
    }
  }
  # Add taxon info ---------------------------------------------------------------------------------
  map_item_to_taxon <- function(col_index) {
    data <- lapply(taxon_data$taxon_id,
                   function(x) unique(item_data[x == item_data$taxon_id, col_index]))
    data[lapply(data, length) == 0] <- NA
    if (max(vapply(data, length, numeric(1))) > 1)
      stop(paste0("taxon_info field ", col_index - 1, " content does not correspond to taxa ids.", 
                  " Perhaps an item_info field would be more appropriate."))
    return(unlist(data))
  }
  taxon_info_cols <- 1 + which(key == "taxon_info")
  if (length(taxon_info_cols) > 0) {
    taxon_data <- cbind(taxon_data, data.frame(lapply(taxon_info_cols, map_item_to_taxon)))    
  }
  # Add arbitrary item ids to item data if necessary -----------------------------------------------
  if (!("item_id" %in% names(item_data))) item_data$item_id <- 1:nrow(item_data)
  # Add counts to taxon_data -----------------------------------------------------------------------
  taxon_data$item_count <- table(unlist(lapply(item_classification, `[[`, "taxon_id")))[taxon_data$taxon_id]
  # Format and return output -----------------------------------------------------------------------
  names(key)[names(key) == ""] <- names(item_data)[seq_along(key) + 1][names(key) == ""] 
  names(item_data)[seq_along(key) + 1] <- names(key)
  row.names(item_data) <- item_data$item_id
  return(list(taxonomy = taxon_id_key, taxa = taxon_data, items = item_data))
}