#' Download representative sequences for a taxon
#'
#' Downloads a sample of sequences meant to evenly capture the diversity of a
#' given taxon. Can be used to get a shallow sampling of vast groups.
#' \strong{CAUTION:} This function can make MANY queries to Genbank depending on
#' arguments given and can take a very long time. Choose your arguments
#' carefully to avoid long waits and needlessly stressing NCBI's servers. Use a
#' downloaded database and a parser from the \code{taxa} package when possible.
#'
#' @param name (\code{character} of length 1) The taxon to download a sample of
#'   sequences for.
#' @param id (\code{character} of length 1) The taxon id to download a sample of
#'   sequences for.
#' @param target_rank (\code{character} of length 1) The finest taxonomic rank
#'   at which to sample. The finest rank at which replication occurs. Must be a
#'   finer rank than \code{taxon}.
#' @param min_counts (named \code{numeric}) The minimum number of sequences to
#'   download for each taxonomic rank. The names correspond to taxonomic ranks.
#' @param max_counts (named \code{numeric}) The maximum number of sequences to
#'   download for each taxonomic rank. The names correspond to taxonomic ranks.
#' @param interpolate_min (\code{logical}) If \code{TRUE}, values supplied to
#'   \code{min_counts} and \code{min_children} will be used to infer the values
#'   of intermediate ranks not specified. Linear interpolation between values of
#'   specified ranks will be used to determine values of unspecified ranks.
#' @param interpolate_max (\code{logical}) If \code{TRUE}, values supplied to
#'   \code{max_counts} and \code{max_children} will be used to infer the values
#'   of intermediate ranks not specified. Linear interpolation between values of
#'   specified ranks will be used to determine values of unspecified ranks.
#' @param min_children (named \code{numeric}) The minimum number sub-taxa of
#'   taxa for a given rank must have for its sequences to be searched. The names
#'   correspond to taxonomic ranks.
#' @param max_children (named \code{numeric}) The maximum number sub-taxa of
#'   taxa for a given rank must have for its sequences to be searched. The names
#'   correspond to taxonomic ranks.
#' @param seqrange (character) Sequence range, as e.g., "1:1000". This is the
#'   range of sequence lengths to search for. So "1:1000" means search for
#'   sequences from 1 to 1000 characters in length.
#' @param getrelated (logical) If TRUE, gets the longest sequences of a species
#'   in the same genus as the one searched for. If FALSE, returns nothing if no
#'   match found.
#' @param fuzzy (logical) Whether to do fuzzy taxonomic ID search or exact
#'   search. If \code{TRUE}, we use \code{xXarbitraryXx[porgn:__txid<ID>]}, but
#'   if \code{FALSE}, we use \code{txid<ID>}. Default: \code{FALSE}
#' @param limit (\code{numeric}) Number of sequences to search for and return.
#'   Max of 10,000. If you search for 6000 records, and only 5000 are found, you
#'   will of course only get 5000 back.
#' @param entrez_query (\code{character}; length 1) An Entrez-format query to
#'   filter results with. This is useful to search for sequences with specific
#'   characteristics. The format is the same as the one used to seach genbank.
#'   (\url{https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options})
#'
#'
#' @param hypothetical (\code{logical}; length 1) If \code{FALSE}, an attempt
#'   will be made to not return hypothetical or predicted sequences judging from
#'   accession number prefixs (XM and XR). This can result in less than the
#'   \code{limit} being returned even if there are more sequences available,
#'   since this filtering is done after searching NCBI.
#' @param verbose (\code{logical}) If \code{TRUE}, progress messages will be
#'   printed.
#'   
#' @examples
#' 
#' \donttest{
#' # Look up 5 ITS sequences from each fungal class
#' data <- ncbi_taxon_sample(name = "Fungi", target_rank = "class", limit = 5, 
#'                           entrez_query = '"internal transcribed spacer"[All Fields]')
#' 
#' # Look up taxonomic information for sequences
#' obj <- lookup_tax_data(data, type = "seq_id", column = "gi_no")
#' 
#' # Plot information
#' metacoder::filter_taxa(obj, taxon_names == "Fungi", subtaxa = TRUE) %>% 
#'   heat_tree(node_label = taxon_names, node_color = n_obs, node_size = n_obs)
#' }
#' @export
ncbi_taxon_sample <- function(name = NULL, id = NULL, target_rank,
                              min_counts = NULL, max_counts = NULL,
                              interpolate_min = TRUE, interpolate_max = TRUE,
                              min_children = NULL, max_children = NULL, 
                              seqrange = "1:3000", getrelated = FALSE,
                              fuzzy = TRUE, limit = 10, entrez_query = NULL,
                              hypothetical = FALSE, verbose = TRUE) {
  # Check that the "traits" package has been installed
  check_for_pkg("traits")
 
  run_once <- function(name, id) {
    default_target_max <- 20
    default_target_min <- 5
    
    taxonomy_levels <- get_taxonomy_levels()
    
    # Argument validation ----------------------------------------------------------------------------
    if (sum(c(is.null(name), is.null(id))) != 1) {
      stop("Either name or id must be speficied, but not both")
    }
    if (!(target_rank %in% levels(taxonomy_levels))) {
      stop("'target_rank' is not a valid taxonomic rank.")
    }
    
    # Argument parsing -------------------------------------------------------------------------------
    if (!is.null(name)) {
      result <- taxize::get_uid(name, verbose = verbose, rows = 1) #This needs attention
      if (is.na(result)) stop(paste0("Could not find taxon ", name))
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
    target_rank <- factor(target_rank,
                          levels = levels(taxonomy_levels),
                          ordered = TRUE)
    
    # Generate taxonomic rank filtering limits ------------------------------------------------------
    get_level_limit <- function(user_limits, default_value, default_level, interpolate, 
                                extend_max = FALSE, extend_min = FALSE) {
      # Provide defaults if NULL - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (is.null(user_limits)) {
        user_limits <- c(default_value)
        names(user_limits) <- default_level
      } else if (length(user_limits) == 1 && is.null(names(user_limits))) {
        names(user_limits) <- default_level
      }
      # Order by taxonomic rank - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
        lapply(seq_len(length(names(user_limits)) - 1), function(i) set_default_counts(names(user_limits)[i:(i+1)]))
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
    
    level_max_count <- get_level_limit(max_counts, default_target_max, target_rank, interpolate_max,
                                       extend_max = TRUE)
    level_min_count <- get_level_limit(min_counts, default_target_min, target_rank, interpolate_min,
                                       extend_min = TRUE)
    level_max_children <- get_level_limit(max_children, NA, target_rank,
                                          interpolate_max, extend_max = TRUE)
    level_min_children <- get_level_limit(min_children, 0, target_rank, interpolate_min,
                                          extend_min = TRUE)
    
    # Recursively sample taxon ------------------------------------------------------------------------
    recursive_sample <- function(id, rank, name) {
      message("Processing '", name, "' (uid: ", id, ", rank: ", as.character(rank), ")", "\n",
              sep = "")
      # Get children of taxon  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!(rank %in% taxonomy_levels) || rank < target_rank) {
        sub_taxa <- taxize::ncbi_children(id = id)[[1]]
      }
      # Filter by subtaxon count - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (rank %in% taxonomy_levels && rank < target_rank) {
        if (!is.na(level_max_children[rank]) && nrow(sub_taxa) > level_max_children[rank]) {
          sub_taxa <- sub_taxa[sample(1:nrow(sub_taxa), level_max_children[rank]), ]
        } else if (!is.na(level_min_children[rank]) && nrow(sub_taxa) < level_min_children[rank]) {
          return(NULL)
        }
      }
      # Search for sequences - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
      if ((rank %in% taxonomy_levels && rank >= target_rank) || (!is.null(sub_taxa) && nrow(sub_taxa) == 0)) {
        message("Getting sequences for ", name, "\n")
        result <- traits::ncbi_searcher(id = id, seqrange = seqrange, getrelated = getrelated,
                                        fuzzy = fuzzy, limit = limit, entrez_query = entrez_query,
                                        hypothetical = hypothetical, verbose = verbose)
      } else {
        child_ranks <- factor(sub_taxa$childtaxa_rank,
                              levels = levels(taxonomy_levels), ordered = TRUE) 
        result <- Map(recursive_sample, sub_taxa$childtaxa_id, child_ranks, sub_taxa$childtaxa_name)
        result <- do.call(rbind, result)
      }
      # Filter by count limits - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (rank %in% taxonomy_levels && !is.null(result)) {
        if (!is.na(level_max_count[rank]) && nrow(result) > level_max_count[rank]) {
          result <- result[sample(1:nrow(result), level_max_count[rank]), ]
        } else if (!is.na(level_min_count[rank]) && nrow(result) < level_min_count[rank]) {
          return(NULL)
        }
      }
      return(result)
    }
    
    recursive_sample(id, taxon_level, name)
  }
  
  
  if (is.null(name)) name = list(NULL)
  if (is.null(id)) id = list(NULL) 
  output <- mapply(run_once, id = id, name = name, SIMPLIFY = FALSE)
  output <- do.call(rbind, output)
  
  # Reformat output
  output <- dplyr::as_tibble(output)
  colnames(output) <- c("ncbi_name", "seq_length", "gene_desc", "acc_no", "gi_no")  
  
  return(output)
}


#' Downloads sequences from ids
#' 
#' Downloads the sequences associated with GenBank accession ids.
#' 
#' @param ids (\code{character}) One or more accession numbers to get sequences for
#' @param batch_size (\code{numeric} of length 1) The number of sequences to request in each query.
#' To large of values might case failures and too small will increase time to completion.
#' 
#' @return (\code{list} of \code{character})
#' 
#' @keywords internal
ncbi_sequence <- function (ids, batch_size = 100) {
  base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
  output <- list()
  while(length(ids) > 0) {
    id_subset <- 1:min(length(ids), batch_size)
    query <- paste(sep = "&",
                   "db=nucleotide",
                   "rettype=fasta",
                   "retmode=text",
                   paste0("id=", paste(ids[id_subset], collapse = ",")))
    raw_result <- RCurl::getURL(paste0(base_url, query))
    temp_path <- tempfile()
    writeChar(raw_result, temp_path)
    result <- ape::read.dna(temp_path, format = "fasta", as.character = TRUE, as.matrix = FALSE)
    if (length(id_subset) == 1)
      result <- stats::setNames(list(result[1,]), dimnames(result)[[1]])
    if (length(id_subset) != length(result))
      stop("Length of query and result do not match. Somthing went wrong.")
    output <- c(output, result)
    ids <- ids[-id_subset]
    Sys.sleep(time = 0.34)
  }
  return(output)
}


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
