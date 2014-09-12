#===================================================================================================
#' Queries sequences for a given taxon
#'
#' Produces a sequinr query object that can be used to download sequences. 
#' 
#' @param query_id A character vector of length 1 specifying the variable name the query will be
#'   bound to.
#' @param taxon A character vector of taxa. Specifing multiple taxa is equivalent to 
#'    running the command multiple times with each taxa and concatenating the results. 
#' @param key A character vector of keywords. All keywords must be present for a given
#'    sequence to be returned.
#' @param type The type of the sequence to return (e.g. RRNA). Use `getType()` to see options. 
#' @return A instance of class seqinr::qaw, as returned by seqinr::query. 
#' @details Based off a exercise in the seqinr vignette.
#' @examples
#' \dontrun{
#' choosebank("genbank")
#' query_taxon("test",
#'             c("phytophthora", "pythium"),
#'             c("@@18S@@", "@@28S@@"),
#'             "RRNA")
#' str(test)
#' closebank()}
#' @seealso \code{\link[taxize]{ncbi_getbyname}} \code{\link[seqinr]{query}}
#' @importFrom seqinr query
#' @export
query_taxon <- function(query_id, taxon, key = character(0), type)   {
  
  single_query <- function(query_id, taxon, key) {
    query("single_query", paste('sp=', taxon, sep=''), virtual=TRUE)
    if (length(key) >= 1) names(key) <- paste("key_query", seq_along(key), sep="_")
    for (item in names(key)) {
      query(item, paste('single_query AND T=', type, ' AND K=', key[item], sep=""), virtual=TRUE)
      query(item, paste("PAR", item), virtual=TRUE) # Replace by parent sequences
    }
    return(query(query_id, paste(names(key), collapse=" AND "), virtual=TRUE))
  }
  
  queries <- mapply(single_query, paste("taxon_query", seq_along(taxon), sep="_"), taxon, key)
  query_names <- apply(queries, MARGIN=2, function(x) x$name)
  query(query_id, paste(query_names, collapse=" OR "))
}


#===================================================================================================
#' Extract the binomial organism name from genbank annotations.
#' 
#' @importFrom stringr str_match
extract_organism <- function(annotation) {
  index <- vapply(annotation, grep, FUN.VALUE=numeric(1), pattern="^SOURCE")
  value <- mapply(`[`, annotation, index)
  str_match(value, "^SOURCE[ \t]+(.+)$")[, 2]
}


#===================================================================================================
#' Extract the description from genbank annotations.
#' 
#' @importFrom stringr str_match
extract_description <- function(annotation) {
  annotation <- vapply(annotation, paste, character(1), collapse="")
  result <- str_match(annotation, "DEFINITION[ \t]+(.+)ACCESSION")[, 2]
  gsub("[ \t]+", " ", result)
}


#===================================================================================================
#' Extract the gi from genbank annotations.
#' 
#' @importFrom stringr str_match
extract_gi <- function(annotation) {
  index <- vapply(annotation, grep, FUN.VALUE=numeric(1), pattern="^VERSION")
  value <- mapply(`[`, annotation, index)
  str_match(value, "^VERSION[ \t]+.+GI:([0-9]+)$")[, 2]
}


#===================================================================================================
#' Download the sequences from an seqinr query object and formats them with their annotations
#' 
#' Formats the output of `sequinr::getSequence` and `sequinr::getAnnot` to the output of 
#' `taxize::ncbi_getbyid`.
#' @param query_req A list of class `SeqAcnucWeb` (e.g. what `sequinr::query` produces in the
#'    `x$req` element of the output.
#' @return A data.frame in the format of the output of `taxize::ncbi_getbyid`.
#' @importFrom seqinr choosebank closebank getSequence getAnnot
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
#' Converts a list of class seqinr::SeqAcnucWeb to a data.frame
#' 
#' @param query_req A list of class seqinr::SeqAcnucWeb
#' @return A data frame with rows named after the sequence names and columns 'length' and 'frame'. 
query_req_to_dataframe <- function(query_req) {
  data.frame(length = vapply(query_req, attr, numeric(1), "length"),
             frame = vapply(query_req, attr, numeric(1), "frame"),
             row.names = as.character(query_req),
             stringsAsFactors = FALSE)
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
                              separate = FALSE,
                              use_acnuc = FALSE) {
  # Verify arguments -------------------------------------------------------------------------------
  subsample <- match.arg(subsample)
  stopifnot(length(seq_length) == 2, seq_length[1] <= seq_length[2])
  # Prepare connection -----------------------------------------------------------------------------
  choosebank("genbank")
  on.exit(closebank())
  
  run_once <- function(taxon) {
    if (taxon == "rhizoctonia") browser()
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
