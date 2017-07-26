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


#' Parse mothur *.tax.summary Classify.seqs output
#' 
#' Parse the `*.tax.summary` file that is returned by the `Classify.seqs` command
#' in mothur.
#' 
#' The input file has a format like:
#' 
#' \preformatted{
#' taxlevel	 rankID	 taxon	 daughterlevels	 total	
#' 0	0	Root	2	242	
#' 1	0.1	Bacteria	50	242	
#' 2	0.1.2	Actinobacteria	38	13	
#' 3	0.1.2.3	Actinomycetaceae-Bifidobacteriaceae	10	13	
#' 4	0.1.2.3.7	Bifidobacteriaceae	6	13	
#' 5	0.1.2.3.7.2	Bifidobacterium_choerinum_et_rel.	8	13	
#' 6	0.1.2.3.7.2.1	Bifidobacterium_angulatum_et_rel.	1	11	
#' 7	0.1.2.3.7.2.1.1	unclassified	1	11	
#' 8	0.1.2.3.7.2.1.1.1	unclassified	1	11	
#' 9	0.1.2.3.7.2.1.1.1.1	unclassified	1	11	
#' 10	0.1.2.3.7.2.1.1.1.1.1	unclassified	1	11
#' 11	0.1.2.3.7.2.1.1.1.1.1.1	unclassified	1	11	
#' 12	0.1.2.3.7.2.1.1.1.1.1.1.1	unclassified	1	11	
#' 6	0.1.2.3.7.2.5	Bifidobacterium_longum_et_rel.	1	2		
#' }
#' 
#' @param file_path (\code{character} of length 1)
#' The file path to the input file.
#' @param unclassified (\code{logical} of length 1)
#' If \code{FALSE}, remove any unclassified rows.
#' 
#' @return \code{\link{taxmap}}
#' 
#' @family parsers
#' 
#' @export
parse_mothur_tax_summary <- function(file_path = NULL, text = NULL,
                                     unclassified = FALSE) {
  
  # Read file
  content <- readLines(file_path)
  
  # Parse header to make key
  header <- strsplit(content[[1]], split = "\t")[[1]]
  key <- c("taxon_info", "class", rep("taxon_info", length(header) - 2))
  key_names <- header
  key_names[2] <- ""
  names(key) <- key_names
  
  # Make regex
  regex <- paste0("^", paste0(collapse = "\t", rep("(.*?)", length(header))), "$")
  
  # Remove 'unclassified' rows
  if (! unclassified) {
    unclassified_rows <- grepl(content, pattern = "^(.*?)\\t(.*?)\\tunclassified\\t")
    content <- content[! unclassified_rows]
    message(paste0("Removed ", sum(unclassified_rows), " unclassified rows."))
  }
  
  # Extract taxonomic data
  result <- extract_taxonomy(content[-1],
                             key = key,
                             regex = regex,
                             class_key = "name",
                             class_sep = "\\.",
                             return_input = FALSE,
                             return_match = FALSE)
  
  # Add 'all' calculated column
  result$taxon_funcs <- c(result$taxon_funcs,
                          list(all = function(obj, subset = obj$taxon_data$taxon_ids) {
                            sample_cols <- header[6:length(header)]
                            sample_cols <- sample_cols[sample_cols %in% colnames(obj$taxon_data)]
                            apply(obj$taxon_data[subset, sample_cols], MARGIN = 1, sum)
                          }))
  
  return(result)
}


#' Parse mothur Classify.seqs *.taxonomy output
#' 
#' Parse the `*.taxonomy` file that is returned by the `Classify.seqs` command
#' in mothur. If confidence scores are present, they are included in the output.
#' 
#' The input file has a format like:
#' 
#' \preformatted{
#' AY457915	Bacteria(100);Firmicutes(99);Clostridiales(99);Johnsone...
#' AY457914	Bacteria(100);Firmicutes(100);Clostridiales(100);Johnso...
#' AY457913	Bacteria(100);Firmicutes(100);Clostridiales(100);Johnso...
#' AY457912	Bacteria(100);Firmicutes(99);Clostridiales(99);Johnsone...
#' AY457911	Bacteria(100);Firmicutes(99);Clostridiales(98);Ruminoco...
#' }
#' 
#' or...
#' 
#' \preformatted{
#' AY457915	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457914	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457913	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457912	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;J...
#' AY457911	Bacteria;Firmicutes;Clostridiales;Ruminococcus_et_rel.;...
#' }
#' 
#' @param file (\code{character} of length 1) The file path to the input file.
#' Either "file" or "text" must be used, but not both.
#' @param text (\code{character}) An alternate input to "file". The contents of
#' the file as a character. Either "file" or "text" must be used, but not both.
#' 
#' @return \code{\link{taxmap}}
#' 
#' @family parsers
#' 
#' @export
parse_mothur_taxonomy <- function(file = NULL, text = NULL) {
  # Check that both `file` and `text` are not used together
  if ((!missing(file) && !missing(text)) || (missing(file) && missing(text))) {
    stop(paste0('Either "file" or "text" must be supplied, but not both.'))
  }
  
  # Convert file or text to a char vector of lines
  if (! missing(file)) {
    raw_lines <- readLines(file)
  } else { # "text" must have been supplied
    raw_lines <- unlist(strsplit(text, "\r\n?|\n"))
  }
  
  # Determine if there are scores associated with each taxon
  parts <- strsplit(raw_lines[1], ";")[[1]]
  has_scores <- all(grepl(parts, pattern = "^(.+)\\(([0-9]+)\\)$"))
  
  # Parse raw lines
  if (has_scores) {
    output <- taxa::extract_tax_data(tax_data = raw_lines,
                                     class_sep = ";",
                                     key = c("sequence_id" = "info",
                                             raw_tax = "class"),
                                     regex = "^(.+)\\t(.+);$",
                                     class_key = c(name = "taxon_name",
                                                   score = "info"),
                                     class_regex = "^(.+)\\(([0-9]+)\\)$")
  } else {
    output <- taxa::extract_tax_data(tax_data = raw_lines,
                                     class_sep = ";",
                                     key = c("sequence_id" = "info",
                                             raw_tax = "class"),
                                     regex = "^(.+)\\t(.+);$")
  }
  
  # report results of parsing
  message(paste0('Parsed ', length(raw_lines), ' lines as ',
                 length(output$taxa), ' unique taxa.'))
  
  return(output)
}