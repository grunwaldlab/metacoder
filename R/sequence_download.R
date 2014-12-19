#===================================================================================================
#' Queries ACNUC for taxa
#'
#' Produces a \pkg{seqinr} query object that can be used to download sequences with other
#' \pkg{seqinr} functions. It can also be parsed for result sequence IDs and lengths. Note that 
#' this function requires an open connection to a database provided by
#' \code{\link[seqinr]{choosebank}} (see examples). 
#' 
#' This function uses \code{\link[seqinr]{query}} to search for sequences.
#' \code{\link[seqinr]{query}} assignes results to a varaible in the global environment as well as
#' returning the variable. By default, this function avoids this behavior by using a temporary
#' variable and deleting it on function exit, even if an error is encountered. If you want to let
#' \code{\link[seqinr]{query}} assign to the global environment, use the \code{query_id} option.
#' 
#' @param taxon A character vector of taxa to search for.  
#' @param key A character vector of keywords. All keywords must be present for a given
#'   sequence to be returned.
#' @param type The type of the sequence to return (e.g. RRNA). Use \code{\link[seqinr]{getType}}
#' to see options. 
#' @param query_id A character vector of length 1. The name of the query passed to
#'   \code{\link[seqinr]{query}}. The result will be saved in a global variable with this name. If
#'   NULL, a temporary variable is made and deleted upon completion of the function. 
#' @param execute If FALSE, the query is not executed. Instead the query string will be returned.
#' @param parent If TRUE, results are replaced with their parent sequences (i.e. subsequences are
#' replaced with full sequences). See ACNUC documentation for more info.
#' @param all_taxa If TRUE, require all taxa to be present for a given seuquence. Typically, this
#'   would be used to specify additional taxonomic levels.
#' @return If \code{execute == TRUE} A instance of class \code{\link[seqinr]{qaw}}, as returned by
#' \code{\link[seqinr]{query}}. Otherwise, a query string will be returned as a character vector.
#' @examples
#' \dontrun{
#' choosebank("genbank")
#' x <- query_taxon(c("Phytophthora", "Pythium"), key=c("@@18S@@"), type="RRNA")
#' closebank()}
#' @seealso
#'   \code{\link[taxize]{ncbi_getbyname}}, 
#'   \code{\link[seqinr]{getSequence}}, 
#'   \code{\link[seqinr]{query}}, 
#'   \code{\link[seqinr]{getType}}
#' @export
query_taxon <- function(taxon, key = NULL, type = NULL, query_id = NULL, execute = TRUE,
                        parent = FALSE, all_taxa = FALSE, ...)   {
  
  query_paste <- function(x, key, sep) {
    if (length(x) == 0) return(NULL)
    sep <- paste0(" ", sep, " ")
    paste0("(", paste(paste0('"', key, "=", x, '"'), collapse = sep), ")")
  }
  # Input validation -------------------------------------------------------------------------------
  if (length(taxon) == 0 || is.null(taxon)) return(NULL)
  if (length(taxon) <= 1 && is.na(taxon)) return(NA)
  # Construct query text ---------------------------------------------------------------------------
  if (all_taxa) and_or <- "AND" else and_or <- "OR"
  search <- c(query_paste(taxon, "SP", and_or),
              query_paste(key, "K", "AND"),
              query_paste(type, "T", "AND"))
  query_text <- paste(search, collapse = " AND ")
  if (execute) {
    # Clean up after seqinr::query -----------------------------------------------------------------
    if (is.null(query_id)) {
      query_id <- R.utils::tempvar(prefix = "query_temp_", envir = .GlobalEnv)
      on.exit(force(rm(list = query_id, envir = .GlobalEnv)))
    }
    # Execute query --------------------------------------------------------------------------------
    result <- seqinr::query(query_id, query_text, ...)
    if (parent) {
      query_text <- paste("PAR", query_id)
      result <- seqinr::query(query_id, query_text, ...)     
    }
    return(result)
  } else {
    return(query_text)
  }
}



#===================================================================================================
#' Extract the binomial organism name from genbank annotations.
#' 
#' Takes a genebank annotation as a character vector with one line per element and returns the 
#' organism name
extract_organism <- function(annotation) {
  index <- vapply(annotation, grep, FUN.VALUE=numeric(1), pattern="^SOURCE")
  value <- mapply(`[`, annotation, index)
  stringr::str_match(value, "^SOURCE[ \t]+(.+)$")[, 2]
}


#===================================================================================================
#' Extract the description from genbank annotations.
#' 
#' Takes a genebank annotation as a character vector with one line per element and returns the 
#' description. 
extract_description <- function(annotation) {
  annotation <- vapply(annotation, paste, character(1), collapse="")
  result <- stringr::str_match(annotation, "DEFINITION[ \t]+(.+)ACCESSION")[, 2]
  gsub("[ \t]+", " ", result)
}


#===================================================================================================
#' Extract the gi from genbank annotations.
#' 
#' Takes a genebank annotation as a character vector with one line per element and returns the 
#' gi.
extract_gi <- function(annotation) {
  index <- vapply(annotation, grep, FUN.VALUE=numeric(1), pattern="^VERSION")
  value <- mapply(`[`, annotation, index)
  stringr::str_match(value, "^VERSION[ \t]+.+GI:([0-9]+)$")[, 2]
}


#===================================================================================================
#' Download the sequences from an seqinr query object and formats them with their annotations
#' 
#' Formats the output of \code{\link[seqinr]{getSequence}} and \code{\link[seqinr]{getAnnot}}
#' to the output of \code{\link[taxize]{ncbi_getbyid}}.
#' 
#' @param query_req A list of class  \code{\link[seqinr]{SeqAcnucWeb}} (e.g. what
#' \code{\link[seqinr]{query}} produces in the \code{x$req} element of the output.
#' @return A data.frame in the format of the output of \code{\link[taxize]{ncbi_getbyid}}.
#' @import seqinr
#' @export
download_gb_query <- function(query_req) {
  choosebank("genbank")
  on.exit(closebank())
  sequence <- getSequence(query_req)
  annotation <- getAnnot(query_req)
  sequence <- vapply(sequence, paste, character(1), collapse="")
  sequence <- toupper(sequence)
  data.frame(taxon = extract_organism(annotation),
             gene_desc = extract_description(annotation),
             gi_no = extract_gi(annotation),
             acc_no = as.character(query_req),
             length = vapply(sequence, nchar, numeric(1)),
             sequence = sequence)
}


#===================================================================================================
#' Converts seqinr::qaw to data.frame
#' 
#' Converts a list of class \code{\link[seqinr]{SeqAcnucWeb}} or \code{\link[seqinr]{qaw}} to a
#' data.frame.
#' 
#' @param query_req A list of class \code{\link[seqinr]{SeqAcnucWeb}} or a single
#'   \code{\link[seqinr]{qaw}} object.
#' @return A data frame with columns corresponding to \code{\link[seqinr]{SeqAcnucWeb}} attributes.
#' @export 
query_req_to_dataframe <- function(query_req) {
  if (class(query_req) == "qaw") query_req <- query_req$req
  data.frame(name = as.character(query_req),
             length = vapply(query_req, attr, numeric(1), "length"),
             frame = vapply(query_req, attr, numeric(1), "frame"),
             stringsAsFactors = FALSE)
}


#===================================================================================================
#' Coerce seqinr::qaw to a Data Frame
#' 
#' Coerce an object of class \code{\link[seqinr]{qaw}} to a Data Frame.
#' 
#' @param x An object of class \code{\link[seqinr]{qaw}}
#' @return A data frame with columns corresponding to \code{\link[seqinr]{SeqAcnucWeb}} attributes.
#' @export
as.data.frame.qaw <- function(x, ...) {
  query_req_to_dataframe(x)
}

#===================================================================================================
#' Downloads sequences that result from a taxon and keyword search
#' 
#' This function is meant to download sequences associated with a taxon and a set of keywords from
#' Genbank using a simple interface. More complicated queries should be done using the tools of the
#' `sequnir` and `taxize` packages, as this function does. 
#' 
#' This function first searches genbank via ACNUC and finds all sequences of a given taxon
#' with given keywords (e.g. genes). It then subsamples these search results if there are too many
#' and downloads the subsample sequences. The results are compiled into a data frame and
#' returned. 
#'   
#' @param taxon A character vector of taxa.
#' @param key A character vector of keywords. All keywords must be present for a given
#'    sequence to be returned.
#' @param type The type of the sequence to return (e.g. RRNA). Use `getType()` to see options. 
#' @param seq_length A numeric vector of length 2. The range of sequence lengths to allow.
#' @param max_count The maximum number of sequences to download. See note below for additional
#'   considerations.
#' @param subsample A character vector of length 1. Specifies how to subsample if more search 
#'   results are found than `max_count`. "subsample": randomly select using `sample`; "head": use
#'   first results; "tail": use last results.
#' @param standardize If TRUE, validate binomial taxon names using `taxize`.
#' @param separate If TRUE, each taxon supplied is searched for separately and the `max_count`
#'   option is applied to each separately, therefore the resulting max count of all results would be  
#'   \code{max_count * length(taxon)}.
#' @param use_acnuc If TRUE, sequences are downloaded using tools from `sequinr`. This is typically
#'   slower. 
#' @note When `use_acnuc = FALSE`, some sequences that are found during searching with `seqinr` are
#'   not found when downloading using `taxize`. This is because `seqinr::query` returns the ACNUC
#'   id, which appears to be the genbank "locus" field, whereas `taxize::ncbi_getbyid` uses the
#'   "accession" or "GI" field. Usually the "locus" and "accession" field are the same, but when 
#'   they are differnt results of the search are not downloaded. Therefore the maximum count of 
#'   sequences will not always be returned even if more than the maximum are found.
#' @seealso \code{\link{download_gb_query}} \code{\link{query_taxon}}
#' @examples
#' \dontrun{
#' x <- download_gb_taxon(c("phytophthora", "pythium"),
#'                        c("@@18S@@", "@@28S@@"),
#'                        "RRNA")}
#' @export
download_gb_taxon <- function(taxon, key, type,
                              seq_length = c(1,10000),
                              max_count = 100,
                              subsample = c("random", "head", "tail"),
                              standardize = TRUE,
                              use_acnuc = FALSE) {
  # Verify arguments -------------------------------------------------------------------------------
  subsample <- match.arg(subsample)
  stopifnot(length(seq_length) == 2, seq_length[1] <= seq_length[2])
  # Prepare connection -----------------------------------------------------------------------------
  seqinr::choosebank("genbank")
  on.exit(seqinr::closebank())
  
  run_once <- function(taxon) {
    cat("Searching: ", taxon)
    # Search for potential sequences ---------------------------------------------------------------
    query_taxon("taxon_query", taxon, key, type)
    results <- query_req_to_dataframe(taxon_query$req)
    results$index <- 1:nrow(results)
    # Filter by sequence length --------------------------------------------------------------------
    results <- results[results$length >= seq_length[1] & results$length <= seq_length[2], ]
    # Subsample if necessary -----------------------------------------------------------------------
    sub_count <- min(nrow(results), max_count)
    results <- switch(subsample,
                      "random" = results[sample(1:nrow(results), sub_count), ],
                      "head"   = head(results, sub_count), 
                      "tail"   = tail(results, sub_count))
    if (nrow(results) == 0) return(NULL)
    # Download sequences ---------------------------------------------------------------------------
    if (use_acnuc) {
      sequences <- download_gb_query(taxon_query$req[results$index]) #seqinr used
    } else {
      sequences <- taxize::ncbi_getbyid(rownames(results), verbose=FALSE) #taxize used
    }
    if (nrow(sequences) == 0) return(NULL)
    # Standardize binomial names -------------------------------------------------------------------
    if (standardize) {
      gnr_result <- taxize::gnr_resolve(sequences$taxon,
                                data_source_ids = 4, #4 is the code for NCBI
                                stripauthority = TRUE,
                                best_match_only = TRUE)
      gnr_result <- gnr_result$result
      sequences$taxon <- gnr_result$matched_name2[order(as.integer(rownames(gnr_result)))]
    }
    cat("complete")
    return(sequences)
  }
  
  if (separate) {
    return(do.call(rbind, lapply(taxon, run_once)))
  } else {
    return(run_once(taxon))
  }
}


#===================================================================================================
#' Download representative sequences for a taxon
#' 
#' Downloads a sample of sequences meant to evenly capture the diversity of a given taxon.
#' 
#' @param taxon A character vector of length 1. The taxon to download a sample of sequences for.
#' @param target_level A character vector of length 1. The finest taxonomic level at which
#'   to sample. The finest level at which replication occurs. Must be a finer level than 
#'   \code{taxon}.
#' @param max_counts A named numeric vector. The maximum number of sequences to download for each
#'   taxonomic level. The names correspond to taxonomic levels. See
#'   \code{\link{get_taxonomy_levels}} or \code{\link[taxize]{rank_ref}} for available taxonomic
#'   levels. 
#' @examples
#' get_taxon_sample(name = "oomycetes", target_level = "genus")
#' fungi <- get_taxon_sample(name = "fungi", target_level = "family", max_counts = c(family = 5, order = 30), entrez_query = "18S[All Fields] AND 28S[All Fields]", min_length = 600, max_length = 10000)
#' @export
get_taxon_sample <- function(name = NULL, id = NULL, target_level, max_counts = NULL,
                             interpolate_max = TRUE, min_counts = NULL, interpolate_min = TRUE,
                             verbose = TRUE, max_length = 10000, min_length = 1,
                             max_children = NULL, min_children = NULL, ...) {
  
  default_target_max <- 20
  default_target_min <- 5
  taxonomy_levels <- get_taxonomy_levels()
  
  # Argument validation ----------------------------------------------------------------------------
  if (sum(c(is.null(name), is.null(id))) != 1) {
    stop("Either name or id must be speficied, but not both")
  }
  if (!(target_level %in% levels(taxonomy_levels))) {
    stop("'target_level' is not a valid taxonomic level.")
  }
  
  # Argument parsing -------------------------------------------------------------------------------
  if (!is.null(name)) {
    result <- taxize::get_uid(name, verbose = verbose)
    if (is.na(result)) stop(cat("Could not find taxon ", name))
    id <- result
  }  else {
    id <- as.character(id)
    attr(id, "class") <- "uid"
  }
  taxon_classification <- taxize::classification(id, db = 'ncbi')[[1]]
  name <- taxon_classification[nrow(taxon_classification), "name"]
  taxon_level <- factor(taxon_classification[nrow(taxon_classification), "rank"],
                        levels = levels(taxonomy_levels),
                        ordered = TRUE)
  target_level <- factor(target_level,
                         levels = levels(taxonomy_levels),
                         ordered = TRUE)
  length_range <- paste(min_length, max_length, sep = ":")
  
  # Generate taxonomic level filtering limits ------------------------------------------------------
  get_level_limit <- function(user_limits, default_value, default_level, interpolate, 
                              extend_max = FALSE, extend_min = FALSE) {
    # Provide defaults if NULL - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(user_limits)) {
      user_limits <- c(default_value)
      names(user_limits) <- default_level
    } else if (length(user_limits) == 1 && is.null(names(user_limits))) {
      names(user_limits) <- default_level
    }
    # Order by taxonomic level - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    limit_levels <- factor(names(user_limits),
                           levels = levels(taxonomy_levels),
                           ordered = TRUE)
    user_limits <- user_limits[order(limit_levels)]
    # place input values in vector with all levels - - - - - - - - - - - - - - - - - - - - - - - - -
    all_user_limits <- rep(as.integer(NA), length(taxonomy_levels))
    names(all_user_limits) <- levels(taxonomy_levels)
    all_user_limits[names(user_limits)] <- user_limits
    # Interpolate limits for undefined intermediate levels - - - - - - - - - - - - - - - - - - - - -
    if (interpolate && length(user_limits) >= 2) {
      set_default_counts <- function(range) {
        between <- which(taxonomy_levels >= range[1] & taxonomy_levels <= range[2])
        all_user_limits[between] <<- as.integer(seq(user_limits[range[1]],
                                                    user_limits[range[2]],
                                                    along.with = between))
        return(NULL)
      }
      zoo::rollapply(names(user_limits), width = 2, set_default_counts)    
    }
    
    # Extend boundry values to adjacent undefined values - - - - - - - - - - - - - - - - - - - - - -
    defined <- which(!is.na(all_user_limits))
    if (length(defined) > 0) {
      if (extend_max) {
        highest_defined <- max(defined)
        all_user_limits[highest_defined:length(all_user_limits)] = all_user_limits[highest_defined]      
      }
      if (extend_min) {
        lowest_defined <- min(defined)
        all_user_limits[1:lowest_defined] = all_user_limits[lowest_defined]      
      }      
    }
    return(all_user_limits)
  }
  
  level_max_count <- get_level_limit(max_counts, default_target_max, target_level, interpolate_max,
                                     extend_max = TRUE)
  level_min_count <- get_level_limit(min_counts, default_target_min, target_level, interpolate_min,
                                     extend_min = TRUE)
  level_max_children <- get_level_limit(max_children, NA, target_level,
                                        interpolate_max, extend_max = TRUE)
  level_min_children <- get_level_limit(min_children, 0, target_level, interpolate_min,
                                        extend_min = TRUE)
  
  # Recursivly sample taxon ------------------------------------------------------------------------
  recursive_sample <- function(id, level, name) {
    cat("Processing '", name, "' (uid: ", id, ", level: ", as.character(level), ")", "\n",
        sep = "")
    # Get children of taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!(level %in% taxonomy_levels) || level < target_level) {
      sub_taxa <- taxize::ncbi_children(id = id)[[1]]
      print(sub_taxa)
    }
    # Filter by subtaxon count - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (level %in% taxonomy_levels && level < target_level) {
      if (!is.na(level_max_children[level]) && nrow(sub_taxa) > level_max_children[level]) {
        sub_taxa <- sub_taxa[sample(1:nrow(sub_taxa), level_max_children[level]), ]
      } else if (!is.na(level_min_children[level]) && nrow(sub_taxa) < level_min_children[level]) {
        return(NULL)
      }
    }
    # Search for sequences - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    if ((level %in% taxonomy_levels && level >= target_level) || (!is.null(sub_taxa) && nrow(sub_taxa) == 0)) {
      cat("Getting sequences for", name, "\n")
      result <- taxize::ncbi_search(id = id, limit = 1000, seqrange = length_range,
                                    hypothetical = TRUE, ...)
    } else {
      child_ranks <- factor(sub_taxa$childtaxa_rank,
                            levels = levels(taxonomy_levels), ordered = TRUE) 
      result <- Map(recursive_sample, sub_taxa$childtaxa_id, child_ranks, sub_taxa$childtaxa_name)
      result <- do.call(rbind, result)
    }
    # Filter by count limits - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (level %in% taxonomy_levels && !is.null(result)) {
      if (!is.na(level_max_count[level]) && nrow(result) > level_max_count[level]) {
        result <- result[sample(1:nrow(result), level_max_count[level]), ]
      } else if (!is.na(level_min_count[level]) && nrow(result) < level_min_count[level]) {
        return(NULL)
      }
    }
    return(result)
  }
  
  recursive_sample(id, taxon_level, name)
}