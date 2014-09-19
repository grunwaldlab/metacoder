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
#' @importFrom seqinr choosebank closebank 
#' @importFrom taxize ncbi_getbyid gnr_resolve
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
  choosebank("genbank")
  on.exit(closebank())
  
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
      sequences <- ncbi_getbyid(rownames(results), verbose=FALSE) #taxize used
    }
    if (nrow(sequences) == 0) return(NULL)
    # Standardize binomial names -------------------------------------------------------------------
    if (standardize) {
      gnr_result <- gnr_resolve(sequences$taxon,
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
#' @param target_level A character vector of length 1. The level finest taxonomic level at which
#'   to sample. The finest level at which replication occurs. Must be a finer level than 
#'   \code{taxon}.
#' @param max_counts A named numeric vector. The maximum number of sequences to download for each
#'   taxonomic level. The names correspond to taxonomic levels. See
#'   \code{\link{get_taxonomy_levels}} or \code{\link[taxize]{rank_ref}} for available taxonomic
#'   levels. 
#' @examples
#' get_taxon_sample(name = "oomycetes", target_level = "Genus")
get_taxon_sample <- function(name = NULL, id = NULL, target_level, max_counts = NULL,
                             interpolate_max = TRUE, min_counts = NULL, interpolate_min = TRUE,
                             verbose = TRUE, ...) {
  
  default_target_max <- 20
  default_target_min <- 5
  taxonomy_levels <- get_taxonomy_levels()
  
  # Argument validation ----------------------------------------------------------------------------
  if (sum(c(is.null(name), is.null(id))) != 1) {
    stop("Either name or id must be speficied, but not both")
  }
  
  # Argument parsing -------------------------------------------------------------------------------
  if (!is.null(name)) {
    result <- get_uid(name, verbose = verbose)
    if (is.na(result)) stop(cat("Could not find taxon ", name))
    id <- result
  }  else {
    id <- as.character(id)
    attr(id, "class") <- "uid"
  }
  taxon_classification <- classification(id, db = 'ncbi')[[1]]
  taxon_level <- factor(taxon_classification[nrow(taxon_classification), "rank"],
                        levels = levels(taxonomy_levels),
                        ordered = TRUE)
  target_level <- factor(target_level,
                         levels = levels(taxonomy_levels),
                         ordered = TRUE)
  
  # Generate taxonomic level sequences count limits ------------------------------------------------
  get_level_limit <- function(user_limits, default_value, default_level, interpolate) {
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
    return(all_user_limits)
  }
  level_max_count <- get_level_limit(max_counts, default_target_max, target_level, interpolate_max)
  level_min_count <- get_level_limit(min_counts, default_target_min, target_level, interpolate_min)
  
  # Recursivly sample taxon ------------------------------------------------------------------------
  recursive_sample <- function(id, level) {
    if (level >= target_level) {
      # Search for sequences - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      taxonomy <- classification(id = id)[[1]]
      taxon_name <- taxonomy$name[nrow(taxonomy)]
      result <- ncbi_search(taxon_name, limit = 1000)
      # Filter by count limits - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (level %in% levels(taxonomy_levels)) {
        if (nrow(results) > level_max_count[level]) {
          result <- result[sample(seq_along(result), level_max_count[level])]
        } else if (nrow(results) < level_min_count[level]) {
          return(NULL)
        }
      }
      return(results)
    } else {
      # Get children of taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      children <- ncbi_children(id)[[1]]
      results <- Map(recursive_sample, children$childtaxa_id, children$childtaxa_rank)
      results <- do.call(rbind, results)
      return(results)
    }
  }
  
}

#' Search NCBI for children of a taxon
#' 
#' Search the NCBI Taxonomy databse for uids of children of taxa. Taxa can be referenced by name
#' or uid. Referencing by name is faster.
#' 
#' In a few cases, different taxa have the same name (e.g. Satyrium; see examples). If one of these
#' are searched for then the children of both taxa will be returned. This can be avoided by
#' using a uid instead of the name or specifying an ancestor. If an ancestor is provided, only 
#' children of both the taxon and its ancestor are returned. This will only fail if there are two
#' taxa with the same name and the same specified ancestor. 
#' 
#' @param name (\code{character}) The string to search for. Only exact matches found the name given
#' will be returned. Not compatible with \code{id}.
#' @param id (\code{character}) The uid to search for. Not compatible with \code{name}.
#' @param start The first record to return. If omitted, the results are returned from the first
#'   record (start=0). 
#' @param max_return (\code{numeric; length=1}) The maximum number of children to return.
#' @param ancestor (\code{character}) The ancestor of the taxon being searched for. This is useful
#'   if there could be more than one taxon with the same name. Has no effect if \code{id} is used.
#' @return A list of character vectors of children uids. The names of the list elements are the
#'   input \code{name} or \code{id} values.
#' @examples
#' ncbi_children_uid(name="Satyrium") #Satyrium is the name of two different genera
#' ncbi_children_uid(name="Satyrium", ancestor="Eumaeini") # A genus of butterflies
#' ncbi_children_uid(name="Satyrium", ancestor="Orchidaceae") # A genus of orchids
#' ncbi_children_uid(id="266948") #"266948" is the uid for the butterfly genus
#' ncbi_children_uid(id="62858") #"62858" is the uid for the orchid genus
#' @export
ncbi_children_uid <- function(name = NULL, id = NULL, start = 0, max_return = 1000,
                              ancestor = NULL) {
  # Argument validation ----------------------------------------------------------------------------
  if (sum(c(is.null(name), is.null(id))) != 1) {
    stop("Either name or id must be speficied, but not both")
  }
  # Get name from id -------------------------------------------------------------------------------
  if (is.null(name)) {
    if (class(id) != 'uid') attr(id, 'class') <- 'uid'
    id_taxonomy <- classification(id, db = 'ncbi')
    name <- vapply(id_taxonomy, function(x) x$name[nrow(x)], character(1))
    ancestor <- vapply(id_taxonomy,
                     function(x) ifelse(nrow(x) > 1, x$name[nrow(x) - 1], NA),
                     character(1)) 
  } else if (is.null(ancestor)) {
    ancestor <- rep(NA, length(name))
  }
  single_search <- function(name, ancestor) {
    # Make eutils esearch query --------------------------------------------------------------------
    base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy"
    if (is.na(ancestor)) {
      ancestor_query <- NULL
    } else {
      ancestor_query <- paste0("+AND+", ancestor, "[subtree]")
    }
    taxon_query <- paste0("term=", name, "[Next+Level]", ancestor_query)
    max_return_query <- paste0("RetMax=", max_return)
    start_query <- paste0("RetStart=", start)
    query <- paste(base_url, taxon_query, max_return_query, start_query, sep="&")
    # Search ncbi for children ---------------------------------------------------------------------
    raw_results <- RCurl::getURL(query)
    # Parse results --------------------------------------------------------------------------------
    results <- XML::xmlTreeParse(raw_results, useInternalNodes = TRUE)
    children_uid <- XML::xpathSApply(results, "//eSearchResult/IdList/Id", XML::xmlValue)
    Sys.sleep(0.34) # NCBI limits requests to three per second
    return(children_uid)
  }
  #Combine the result of multiple searches ----------------------------------------------------------
  output <- Map(single_search, name, ancestor)
  if (is.null(id)) names(output) <- name else names(output) <- id
  return(output)
}



ncbi_get_taxon_summary <- function(id) {
  # Argument validation ----------------------------------------------------------------------------
  if (length(id) <= 1 && is.na(id)) return(NA)
  if (is.null(id)) return(NULL)
  id <- as.character(id)
  # Make eutils esummary query ---------------------------------------------------------------------
  base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy"
  query <- paste0(base_url, "&id=", paste(id, collapse = "+"))
  # Search ncbi taxonomy for uid -------------------------------------------------------------------
  raw_results <- RCurl::getURL(query)
  # Parse results ----------------------------------------------------------------------------------
  results <- XML::xmlTreeParse(raw_results, useInternalNodes = TRUE)
  output <- data.frame(stringsAsFactors = FALSE,
    uid = XML::xpathSApply(results, "/eSummaryResult//DocSum/Id", XML::xmlValue),
    name = XML::xpathSApply(results, "/eSummaryResult//DocSum/Item[@Name='ScientificName']",
                            XML::xmlValue),
    rank = XML::xpathSApply(results, "/eSummaryResult//DocSum/Item[@Name='Rank']", XML::xmlValue)
    )
  return(output)
}