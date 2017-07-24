#' Parse HMP QIIME results
#' 
#' NOTE: Not extensively tested
#' Parses the results of the Human Microbiome Project QIIME analysis of the 16S metagenomic data.
#' This is mostly a wrapper for \code{\link{extract_taxonomy}}.
#' 
#' @param otu_file (\code{character} of length 1) A file path or URL to the OTU table. 
#' @param mapping_file (\code{character} of length 1) A file path or URL to the sample mapping file.
#' @param min_abundance (\code{numeric} of length 1) Do not return rows less abundance. This can make the output object much smaller. 
#' @param max_otus (\code{numeric} of length 1) Only parse some number of OTUs. Useful for making small datasets for testing.
#' 
#' @return A \code{\link{taxmap}} object
#' 
#' @export
parse_hmp_qiime <- function(otu_file, mapping_file, min_abundance = 1, max_otus = -1) {
  
  # Download file if url
  download_if_url <- function(path) {
    if (grepl(pattern = "^http[s]?://", path)) {
      temp_file_path <- file.path(tempdir(), basename(path))
      utils::download.file(url = path, destfile = temp_file_path, quiet = TRUE)
      path <- temp_file_path
    }
    return(path)
  }
  otu_file <- download_if_url(otu_file)
  mapping_file <- download_if_url(mapping_file)
  
  # Parse OTU table
  otu_raw <- utils::read.table(otu_file, header = TRUE, skip = 1, comment.char = "",
                        stringsAsFactors = FALSE, sep = "\t", nrows = max_otus)
  # Fix header problems
  colnames(otu_raw) <- gsub(colnames(otu_raw), pattern = "^X\\.?", replacement = "")
  colnames(otu_raw)[1] <- "otu_id"
  # Parse taxonomy information
  otu_data <- extract_taxonomy(otu_raw[ , ncol(otu_raw)],
                               class_regex = "([a-z]{0,1})_{0,2}(.*)$",
                               class_key = c(rank = "taxon_info", "name"),
                               class_sep = ";")
  # Add OTU counts to observation data
  otu_data$obs_data <- dplyr::bind_cols(otu_data$obs_data, convert_numeric_cols(otu_raw[ , -ncol(otu_raw)]))
  # # Convert to long format
  # gather_one <- function(index) {
  #   part <- tidyr::gather(otu_data$obs_data[index, ], key = sample_id, value = abundance, 3:ncol(otu_data$obs_data[index, ]))
  #   part[part$abundance >= min_abundance, ] 
  # }
  # 
  # # otu_data$obs_data <- dplyr::bind_rows(lapply(1:nrow(otu_data$obs_data[1:10, ]), gather_one))
  # 
  # 
  # 
  # # Add taxon column generators
  # taxon_abundance <- function(obj, subset = 1:nrow(obj$taxon_data)) {
  #   vapply(obs(obj, subset = subset, recursive = TRUE, simplify = FALSE, index = TRUE), FUN.VALUE = numeric(1), 
  #          function(x) sum(obj$obs_data[x, ]))
  # }
  # otu_data$taxon_funcs <- c(otu_data$taxon_funcs, list(taxon_abundance = taxon_abundance))
  
  
  # Parse mapping table 
  mapping_data <- utils::read.table(mapping_file, header = TRUE, comment.char = "",
                             stringsAsFactors = FALSE, sep = "\t")[1:7]
  colnames(mapping_data) <- c("sample_id", "rsid", "visit_no", "sex", "run_center", "body_site", "description")
  otu_data$mapping <- dplyr::tbl_df(mapping_data)
  
  return(otu_data)
}



#' @export
parse_phyloseq <- function(obj) {
  datasets <- list()
  
  # Parse taxonomic data
  possible_ranks <- unique(unlist(strsplit(taxa::ranks_ref$ranks, split = ",")))
  tax_data <- as.data.frame(obj@tax_table)
  
  # Parse OTU tables
  if (! is.null(obj@otu_table)) {
    otu_table <- obj@otu_table
    if (! otu_table@taxa_are_rows) {
      otu_table <- t(otu_table)
    }
    otu_table <- as.data.frame(otu_table)
    datasets <- c(datasets, list(otu_table = otu_table))
  }
  
  # Parse sample data
  if (! is.null(obj@sam_data)) {
    datasets <- c(datasets, list(sam_data = obj@sam_data))
  }
  
  # Parse phylogenetic tree
  if (! is.null(obj@phy_tree)) {
    datasets <- c(datasets, list(phy_tree = obj@phy_tree))
  }
  
  # Parse reference sequences
  if (! is.null(obj@refseq)) {
    refseq <- as.character(obj@refseq)
    datasets <- c(datasets, list(refseq = refseq))
  }
  
  # Construct output
  parse_tax_data(tax_data = tax_data, 
                 datasets = datasets,
                 class_cols = which(tolower(colnames(tax_data)) %in% possible_ranks), 
                 mappings = c("{{name}}" = "{{name}}",
                              NA, 
                              NA, 
                              "{{name}}" = "{{name}}"))
}