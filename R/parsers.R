#' Convert a phyloseq to taxmap
#'
#' Converts a phyloseq object to a taxmap object.
#'
#' @param obj A phyloseq object
#' @param class_regex A regular expression used to parse data in the taxon
#'   names. There must be a capture group (a pair of parentheses) for each item
#'   in \code{class_key}. See \code{\link[taxa]{parse_tax_data}} for examples of
#'   how this works.
#' @inheritParams taxa::parse_tax_data
#'
#' @return A taxmap object
#'
#' @family parsers
#'
#' @examples \dontrun{
#'
#' # Install phyloseq to get example data
#' # source('http://bioconductor.org/biocLite.R')
#' # biocLite('phyloseq')
#'
#' # Parse example dataset
#' library(phyloseq)
#' data(GlobalPatterns)
#' x <- parse_phyloseq(GlobalPatterns)
#'
#' # Plot data
#' heat_tree(x,
#'           node_size = n_obs,
#'           node_color = n_obs,
#'           node_label = taxon_names,
#'           tree_label = taxon_names)
#'
#' }
#'
#' @import taxa
#'
#' @export
parse_phyloseq <- function(obj, class_regex = "(.*)",
                           class_key = "taxon_name") {
  datasets <- list()
  mappings <- c()
  
  # Parse taxonomic data
  possible_ranks <- unique(unlist(strsplit(taxa::ranks_ref$ranks, split = ",")))
  tax_data <- as.data.frame(obj@tax_table, stringsAsFactors = FALSE)
  tax_cols <- colnames(tax_data)
  tax_data <- cbind(data.frame(otu_id = rownames(tax_data), stringsAsFactors = FALSE), tax_data)
  
  # Parse OTU tables
  if (! is.null(obj@otu_table)) {
    otu_table <- obj@otu_table
    if (! otu_table@taxa_are_rows) {
      otu_table <- t(otu_table)
    }
    otu_table <- as.data.frame(otu_table, stringsAsFactors = FALSE)
    otu_table <- cbind(data.frame(otu_id = rownames(otu_table), stringsAsFactors = FALSE), otu_table)
    datasets <- c(datasets, list(otu_table = otu_table))
    mappings <- c(mappings, c("{{name}}" = "{{name}}"))
  }
  
  # Parse sample data
  if (! is.null(obj@sam_data)) {
    sam_data <- as.data.frame(as.list(obj@sam_data), stringsAsFactors = FALSE)
    if (! is.null(rownames(obj@sam_data))) {
      sam_data <- cbind(sample_id = rownames(obj@sam_data), sam_data)
    }
    sam_data[] <- lapply(sam_data, as.character)
    datasets <- c(datasets, list(sample_data = sam_data))
    mappings <- c(mappings, NA)
  }
  
  # Parse phylogenetic tree
  if (! is.null(obj@phy_tree)) {
    datasets <- c(datasets, list(phy_tree = obj@phy_tree))
    mappings <- c(mappings, NA)
  }
  
  # Parse reference sequences
  if (! is.null(obj@refseq)) {
    refseq <- as.character(obj@refseq)
    datasets <- c(datasets, list(ref_seq = refseq))
    mappings <- c(mappings, c("{{name}}" = "{{name}}"))
  }
  
  # Construct output
  output <- taxa::parse_tax_data(tax_data = tax_data, 
                                 datasets = datasets,
                                 class_cols = tax_cols, 
                                 mappings = mappings,
                                 named_by_rank = TRUE,
                                 class_regex = class_regex,
                                 class_key = class_key)
  
  # Remove NA taxa
  withCallingHandlers({
    output$filter_taxa(output$taxon_names() != "NA")
  }, warning=function(w) {
    if (conditionMessage(w) %in% c(
      'There is no "taxon_id" column in the data set "3", so there are no taxon IDs.',
      'The data set "4" is named, but not named by taxon ids.'
    ))
      invokeRestart("muffleWarning")
  })
  
  # Move OTU table to front of data if it is there
  if ("otu_table" %in% names(output$data)) {
    otu_tab_index <- which(names(output$data) == "otu_table")
    output$data <- c(output$data[otu_tab_index], output$data[-otu_tab_index])
  }
  
  return(output)
}


#' Parse mothur *.tax.summary Classify.seqs output
#' 
#' Parse the `*.tax.summary` file that is returned by the `Classify.seqs` command
#' in mothur.
#' 
#' The input file has a format like:
#' 
#' \preformatted{
#' taxlevel	 rankID	 taxon	 daughterlevels	 total	A	B	C	
#' 0	0	Root	2	242	84	84	74	
#' 1	0.1	Bacteria	50	242	84	84	74	
#' 2	0.1.2	Actinobacteria	38	13	0	13	0	
#' 3	0.1.2.3	Actinomycetaceae-Bifidobacteriaceae	10	13	0	13	0	
#' 4	0.1.2.3.7	Bifidobacteriaceae	6	13	0	13	0	
#' 5	0.1.2.3.7.2	Bifidobacterium_choerinum_et_rel.	8	13	0	13	0	
#' 6	0.1.2.3.7.2.1	Bifidobacterium_angulatum_et_rel.	1	11	0	11	0	
#' 7	0.1.2.3.7.2.1.1	unclassified	1	11	0	11	0	
#' 8	0.1.2.3.7.2.1.1.1	unclassified	1	11	0	11	0	
#' 9	0.1.2.3.7.2.1.1.1.1	unclassified	1	11	0	11	0	
#' 10	0.1.2.3.7.2.1.1.1.1.1	unclassified	1	11	0	11	0	
#' 11	0.1.2.3.7.2.1.1.1.1.1.1	unclassified	1	11	0	11	0	
#' 12	0.1.2.3.7.2.1.1.1.1.1.1.1	unclassified	1	11	0	11	0	
#' 6	0.1.2.3.7.2.5	Bifidobacterium_longum_et_rel.	1	2	0	2	0	
#' 7	0.1.2.3.7.2.5.1	unclassified	1	2	0	2	0	
#' 8	0.1.2.3.7.2.5.1.1	unclassified	1	2	0	2	0	
#' 9	0.1.2.3.7.2.5.1.1.1	unclassified	1	2	0	2	0
#' }
#' 
#' or 
#' 
#' \preformatted{
#' taxon	total	A	B	C
#' "k__Bacteria";"p__Actinobacteria";"c__Actinobacteria";...	1	0	1	0
#' "k__Bacteria";"p__Actinobacteria";"c__Actinobacteria";...	1	0	1	0
#' "k__Bacteria";"p__Actinobacteria";"c__Actinobacteria";...	1	0	1	0
#' }
#' 
#' @param file (\code{character} of length 1) The file path to the input file. 
#'   Either "file", "text", or "table" must be used, but only one.
#' @param text (\code{character}) An alternate input to "file". The contents of 
#'   the file as a character. Either "file", "text", or "table" must be used,
#'   but only one.
#' @param table (\code{character} of length 1) An already parsed data.frame or
#'   tibble. Either "file", "text", or "table" must be used, but only one.
#' 
#' @return \code{\link{taxmap}}
#' 
#' @family parsers
#' 
#' @export
parse_mothur_tax_summary <- function(file = NULL, text = NULL, table = NULL) {
  
  # Check that `file` and `text` and `table` are not used together
  are_missing <- c(file = is.null(file),
                   text = is.null(text),
                   table = is.null(table))
  if (sum(are_missing) != 2) {
    stop(paste0('Either "file", "text", or "table" must be supplied, but only one.'))
  }
  
  # Read raw data
  if (! are_missing["file"]) {
    raw_data <- utils::read.csv(file = file, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE)
  } else if (! are_missing["text"]) {
    raw_data <- utils::read.csv(text = text, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE)
  } else {
    if (!is.data.frame(table)) {
      stop('The "table" input requires a data.frame or tibble.')
    }
    raw_data <- table
  }
  
  # Check that it is an accepted format
  detailed_cols <- c("taxlevel", "rankID", "taxon", "daughterlevels", "total")
  simple_cols <- c("taxon",	"total")
  if (all(detailed_cols %in% colnames(raw_data))) {
    is_detailed <- TRUE
  } else if (all(simple_cols %in% colnames(raw_data))) {
    is_detailed <- FALSE
  } else {
    stop("Format not recognized.")
  }
  
  
  if (is_detailed) {
    # parse raw table
    output <- taxa::parse_tax_data(tax_data = raw_data,
                                   class_cols = "rankID",
                                   class_sep = ".")
    # replace taxon names
    my_taxon_names <- output$map_data_(output$taxon_ids(),
                                       output$get_data("taxon")[[1]])
    output$taxa <- stats::setNames(lapply(seq_len(length(output$taxa)),
                                          function(i) {
                                            my_taxon <- output$taxa[[i]]
                                            my_taxon$name$name <- my_taxon_names[i]
                                            return(my_taxon)
                                          }),
                                   names(output$taxa))
  } else { # is simple format
    output <- taxa::parse_tax_data(tax_data = raw_data,
                                   class_cols = "taxon",
                                   class_sep = ";")
  }
  
  return(output)
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
#' @import taxa
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


#' Parse a BIOM output from QIIME
#' 
#' Parses a file in BIOM format from QIIME into a taxmap object.
#' This also seems to work with files from MEGAN.
#' I have not tested if it works with other BIOM files. 
#' 
#' This function was inspired by the tutorial created by Geoffrey Zahn at 
#' http://geoffreyzahn.com/getting-your-otu-table-into-r/.
#' 
#' @param file (\code{character} of length 1) The file path to the input file.
#' @param class_regex A regular expression used to parse data in the taxon
#'   names. There must be a capture group (a pair of parentheses) for each item
#'   in \code{class_key}. See \code{\link[taxa]{parse_tax_data}} for examples of
#'   how this works.
#' @inheritParams taxa::parse_tax_data
#' 
#' @return A taxmap object
#' 
#' @family parsers
#' 
#' @import taxa
#' 
#' @export
parse_qiime_biom <- function(file, class_regex = "(.*)",
                             class_key = "taxon_name") {
  # Check that the "biomformat" package has been installed
  check_for_pkg("biomformat")
  
  # Read biom file
  my_biom <- biomformat::read_biom(file)
  
  # Get taxonomy
  taxonomy <- biomformat::observation_metadata(my_biom)
  if (is.null(taxonomy)) {
    stop(call. = FALSE,
         'Could not find taxonomy data in "', file,
         '". If it does have taxonomy data, then it is not where this function expects it to be. ',
         'If you think this is a bug, let us know at "', repo_url(), '/issues"')
  }
  tax_cols <- colnames(taxonomy)
  
  # Get OTU IDs
  if (is.data.frame(taxonomy)) {
    otu_ids <- rownames(taxonomy)
  } else {
    otu_ids <- names(taxonomy)
  }
  
  # Coerce into a matrix
  otu_table <- as.data.frame(as.matrix(biomformat::biom_data(my_biom)),
                             stringsAsFactors = FALSE)
  otu_table <- cbind(list(otu_id = otu_ids), otu_table,
                     stringsAsFactors = FALSE)
  
  # Get sample metadata (not used yet)
  metadata <- biomformat::sample_metadata(my_biom)
  
  # Create a taxmap object
  if (is.data.frame(taxonomy)) {
    output <- taxa::parse_tax_data(tax_data = taxonomy,
                                   class_cols = colnames(taxonomy),
                                   datasets = list(otu_table = otu_table),
                                   mappings = c("{{name}}" = "{{name}}"),
                                   include_tax_data = FALSE,
                                   class_regex = class_regex,
                                   class_key = class_key)
  } else {
    output <- taxa::parse_tax_data(tax_data = taxonomy,
                                   datasets = list(otu_table = otu_table),
                                   mappings = c("{{name}}" = "{{name}}"),
                                   include_tax_data = FALSE,
                                   class_regex = class_regex,
                                   class_key = class_key)
  }
  
  return(output)
}


#' Parse a Newick file
#' 
#' Parse a Newick file into a taxmap object.
#' 
#' The input file has a format like:
#' 
#' \preformatted{
#' (ant:17, (bat:31, cow:22):7, dog:22, (elk:33, fox:12):40);
#' (dog:20, (elephant:30, horse:60):20):50;
#' }
#' 
#' @param file (\code{character} of length 1) The file path to the input file. Either \code{file} or \code{text} must be supplied but not both.
#' @param text (\code{character} of length 1) The raw text to parse. Either \code{file} or \code{text} must be supplied but not both.
#' 
#' @return \code{\link{taxmap}}
#' 
#' @family parsers
#' 
#' @import taxa
#' 
#' @export
parse_newick <- function(file = NULL, text = NULL) {
  # Check that `file` and `text` and `table` are not used together
 if (sum(c(is.null(file), is.null(text))) != 1) {
    stop(paste0('Either "file" or "text" must be supplied, but not both.'))
  }
  
  # Read raw data
  if (is.null(file)) {
    file <- tempfile()
    readr::write_lines(text, file)
  }
  
  # Read input
  raw_data <- phylotate::read_annotated(file, format = "newick")
  
  # Parse edge list
  output <- parse_phylo(raw_data)
  
  return(output)         
}



#' Parse UNITE general release FASTA
#' 
#' Parse the UNITE general release FASTA file
#' 
#' The input file has a format like:
#' 
#' \preformatted{
#' >Glomeromycota_sp|KJ484724|SH523877.07FU|reps|k__Fungi;p__Glomeromycota;c__unid...
#' ATAATTTGCCGAACCTAGCGTTAGCGCGAGGTTCTGCGATCAACACTTATATTTAAAACCCAACTCTTAAATTTTGTAT...
#' }
#' 
#' @inheritParams parse_seq_input
#' @param include_seqs (\code{logical} of length 1) If \code{TRUE}, include
#'   sequences in the output object.
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family parsers
#'   
#' @import taxa
#' 
#' @export
parse_unite_general <- function(input = NULL, file = NULL, include_seqs = TRUE) {
  
  # Read sequence info
  seqs <- parse_seq_input(input = input, file = file)
  headers <- names(seqs)

  # Create taxmap object
  output <- taxa::extract_tax_data(tax_data = headers,
                                   regex = "^(.*)\\|(.*)\\|(.*)\\|(.*)\\|(.*)$",
                                   key = c(organism = "info",
                                           acc_num = "info",
                                           unite_id = "info",
                                           unite_type = "info",
                                           tax_string = "class"),
                                   class_regex = "^(.*)__(.*)$",
                                   class_key = c(unite_rank = "info",
                                                 name = "taxon_name"),
                                   class_sep = ";")
  
  # Remove unneeded columns
  output$data$tax_data$input <- NULL
  output$data$tax_data$tax_string <- NULL
  
  # Add sequences 
  if (include_seqs) {
    output$data$tax_data$unite_seq <- toupper(seqs)
  }
  
  return(output)
}



#' Parse RDP FASTA release
#' 
#' Parses an RDP reference FASTA file.
#' 
#' The input file has a format like:
#' 
#' \preformatted{
#' >S000448483 Sparassis crispa; MBUH-PIRJO&ILKKA94-1587/ss5	Lineage=Root;rootrank;Fun...
#' ggattcccctagtaactgcgagtgaagcgggaagagctcaaatttaaaatctggcggcgtcctcgtcgtccgagttgtaa
#' tctggagaagcgacatccgcgctggaccgtgtacaagtctcttggaaaagagcgtcgtagagggtgacaatcccgtcttt
#' ...
#' }
#' 
#' @inheritParams parse_seq_input
#' @param include_seqs (\code{logical} of length 1) If \code{TRUE}, include 
#'   sequences in the output object.
#' @param add_species (\code{logical} of length 1) If \code{TRUE}, add the
#'   species information to the taxonomy. In this database, the species name
#'   often contains other information as well.
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family parsers
#'   
#' @import taxa
#' 
#' @export
parse_rdp <- function(input = NULL, file = NULL, include_seqs = TRUE, add_species = FALSE) {
  
  # Read sequence info
  raw_data <- parse_seq_input(input = input, file = file)
  headers <- names(raw_data)
  
  # Add species to classification if present
  if (add_species) {
    org_name <- stringr::str_match(headers, "^.*? (.*)\\tLineage=.*$")[, 2]
    genus <- stringr::str_match(headers, ";([a-zA-Z]+);genus$")[, 2]
    species <- vapply(seq_len(length(org_name)), FUN.VALUE = character(1),  function (i) {
      sub(org_name[i], pattern = paste0("^", genus[i], " ?"), replacement = "") %>%
        sub(pattern = ";.*$", replacement = "")
    })
    headers <- paste0(headers, ";", species, ";",
                      ifelse(endsWith(headers, ";"), "", "species"))
  }
  
  # Create taxmap object
  output <- taxa::extract_tax_data(tax_data = headers,
                                   regex = "^(.*?) (.*)\\tLineage=(.*)$",
                                   key = c(rdp_id = "info", seq_name = "info",
                                           tax_string = "class"),
                                   class_regex = "(.+?);(.*?)(?:;|$)",
                                   class_key = c(name = "taxon_name",
                                                 rdp_rank = "info"))
  
  # Add sequences 
  if (include_seqs) {
    seqs <- raw_data[output$data$tax_data$input]
    output$data$tax_data$rdp_seq <- tolower(seqs)
  }
  
  # Remove unneeded columns
  output$data$tax_data$input <- NULL
  output$data$tax_data$tax_string <- NULL
  
  return(output)
}



#' Parse SILVA FASTA release
#'
#' Parses an SILVA FASTA file that can be found at
#' \url{https://www.arb-silva.de/no_cache/download/archive/release_128/Exports/}.
#'
#' The input file has a format like:
#'
#' \preformatted{ >GCVF01000431.1.2369
#' Bacteria;Proteobacteria;Gammaproteobacteria;Oceanospiril...
#' CGUGCACGGUGGAUGCCUUGGCAGCCAGAGGCGAUGAAGGACGUUGUAGCCUGCGAUAAGCUCCGGUUAGGUGGCAAACA
#' ACCGUUUGACCCGGAGAUCUCCGAAUGGGGCAACCCACCCGUUGUAAGGCGGGUAUCACCGACUGAAUCCAUAGGUCGGU
#' ... }
#'
#' @inheritParams parse_seq_input
#' @param include_seqs (\code{logical} of length 1) If \code{TRUE}, include
#'   sequences in the output object.
#'
#' @return \code{\link{taxmap}}
#'
#' @family parsers
#'
#' @import taxa
#'
#' @export
parse_silva_fasta <- function(file = NULL, input = NULL, include_seqs = TRUE) {
  
  # Read sequence info
  raw_data <- parse_seq_input(input = input, file = file)
  raw_headers <- names(raw_data)
  
  # Make classifications easier to parse
  name_chars <- "A-Za-z0-9._+ ='\"\\-"
  parts <- stringr::str_match(raw_headers,
                              paste0("^(.+;)([", name_chars, "]+)(\\(?.*\\)?)$"))
  parts <- as.data.frame(parts[, -1], stringsAsFactors = FALSE)
  colnames(parts) <- c("tax", "binom", "common")
  parts$binom <- sub(parts$binom, pattern = "sp\\. ", replacement = "sp\\._")
  parts$binom <- sub(parts$binom, pattern = "uncultured ", replacement = "uncultured_")
  # parts$binom <- gsub(pattern = " ", replacement = ";", parts$binom) 
  parts$binom <- sub(pattern = ";$", replacement = " ", parts$binom)
  headers <- apply(parts, MARGIN = 1, paste0, collapse = "")
  headers <- gsub(headers, pattern = "\\[|\\]", replacement = "")
  
  # Create taxmap object
  output <- taxa::extract_tax_data(tax_data = headers,
                                   regex = "^(.*)\\.([0-9]*)\\.([0-9]*) (.*)$",
                                   key = c(ncbi_id = "info",
                                           start_pos = "info",
                                           end_pos = "info",
                                           tax_string = "class"),
                                   class_regex = paste0("^([", name_chars, "]+) ?\\(?([", name_chars, "]*)\\)?$"),
                                   class_key = c("name" = "taxon_name", other_name = "info"),
                                   class_sep = ";")
  
  # Clean up taxon names
  output$taxa <- lapply(output$taxa, function(x) {
    x$name$name <- sub(pattern = " $", replacement = "", x$name$name)
    x$name$name <- sub(pattern = "_", replacement = " ", x$name$name) # undo the sp._xxx hack above
    return(x)
  })
  
  # Add sequences 
  if (include_seqs) {
    output$data$tax_data$silva_seq <- raw_data
  }
  
  # Remove unneeded columns
  # output$data$tax_data$input <- NULL
  output$data$tax_data$tax_string <- NULL
  
  # Filter uninformative rows in class_data
  output$data$class_data <- output$data$class_data[output$data$class_data$other_name != "", ]
  output$data$class_data$name <- trimws(output$data$class_data$name)
  
  return(output)
}


#' Parse Greengenes release
#' 
#' Parses the greengenes database.
#' 
#' The taxonomy input file has a format like:
#' 
#' \preformatted{
#' 228054  k__Bacteria; p__Cyanobacteria; c__Synechococcophycideae; o__Synech...
#' 844608  k__Bacteria; p__Cyanobacteria; c__Synechococcophycideae; o__Synech...
#' ...
#' }
#' 
#' The optional sequence file has a format like:
#' 
#' \preformatted{
#' >1111886
#' AACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTA...
#' >1111885
#' AGAGTTTGATCCTGGCTCAGAATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGTACGAGAAATCCCGAGC...
#' ...
#' }
#' 
#' @param tax_file (\code{character} of length 1) The file path to the
#'   greengenes taxonomy file.
#' @param seq_file (\code{character} of length 1) The file path to the
#'   greengenes sequence fasta file. This is optional.
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family parsers
#'   
#' @import taxa
#' 
#' @export
parse_greengenes <- function(tax_file, seq_file = NULL) {
  # Parse taxonomy file
  tax_data <- utils::read.table(tax_file, sep = "\t")
  colnames(tax_data) <- c("gg_id", "classification")
  result <- taxa::parse_tax_data(tax_data, class_cols = "classification", 
                                 class_sep = "; ",
                                 class_regex = "^([a-z]{1})__(.*)$", 
                                 class_key = c("gg_rank" = "info", "name" = "taxon_name"))
  result$data$tax_data$gg_id <- as.character(result$data$tax_data$gg_id)
  
  # Remove data for ranks with no information
  result <- result$filter_taxa(result$taxon_names() != "", drop_obs = TRUE, 
                               reassign_obs = c(tax_data = TRUE, class_data = FALSE))
  
  # Integrating sequence and taxonomy
  if (! is.null(seq_file)) {
    gg_sequences <- seqinr::read.fasta(seq_file, as.string = TRUE)
    result <- result$mutate_obs("tax_data", gg_seq = toupper(unlist(gg_sequences)[result$data$tax_data$gg_id]))
  }
  
  return(result)
}



#' Parse a phylo object 
#' 
#' Parses a phylo object from the ape package.
#' 
#' @param obj A phylo object from the ape package.
#'   
#' @return \code{\link{taxmap}}
#'   
#' @family parsers
#'   
#' @import taxa
#' 
#' @export
parse_phylo  <- function(obj) {
  # Parse edge list
  edge_list <- as.data.frame(obj$edge)
  colnames(edge_list) <- c("from", "to")
  edge_list$from <- as.character(edge_list$from)
  edge_list$to <- as.character(edge_list$to)
  
  # Add roots to edge list
  is_root <- vapply(seq_len(nrow(edge_list)),
                    function(i) ! edge_list$from[i] %in% edge_list$to,
                    FUN.VALUE = logical(1))
  roots <- unique(edge_list$from[is_root])
  edge_list <- rbind(data.frame(from = NA, to = roots), edge_list)
  
  # Parse edge length
  edge_length <- rep(NA, nrow(edge_list))
  edge_length[!is.na(edge_list$from)] <- obj$edge.length
  
  # Parse tip labels
  is_tip <-  vapply(seq_len(nrow(edge_list)),
                    function(i) ! edge_list$to[i] %in% edge_list$from,
                    FUN.VALUE = logical(1))
  
  tip_label <- rep(NA, nrow(edge_list))
  tip_label[is_tip] <- obj$tip.label
  
  # Build taxmap object
  output <- taxa::taxmap()
  output$edge_list <- edge_list
  taxon_ids <- unique(unlist(output$edge_list))
  taxon_ids <- taxon_ids[!is.na(taxon_ids)]
  output$taxa <- stats::setNames(lapply(paste0("node_", taxon_ids), taxa::taxon),
                                 taxon_ids)
  tax_data <- dplyr::as.tbl(data.frame(stringsAsFactors = FALSE,
                                       taxon_id = edge_list$to,
                                       edge_length = edge_length,
                                       tip_label = tip_label))
  output$data <- c(output$data, list(tax_data = tax_data))
  output$replace_taxon_ids(convert_base(as.integer(output$taxon_ids())))
}




#' Converts the uBiome file format to taxmap
#' 
#' Converts the uBiome file format to taxmap. NOTE: This is experimental and might not work if
#' uBiome changes their format. Contact the maintainers if you encounter problems/
#' 
#' The input file has a format like:
#' 
#' \preformatted{
#'  tax_name,tax_rank,count,count_norm,taxon,parent
#'  root,root,29393,1011911,1,
#'  Bacteria,superkingdom,29047,1000000,2,131567
#'  Campylobacter,genus,23,791,194,72294
#'  Flavobacterium,genus,264,9088,237,49546
#' }
#' 
#' @param file (\code{character} of length 1) The file path to the input file. 
#'   Either "file", or "table" must be used, but only one.
#' @param table (\code{character} of length 1) An already parsed data.frame or
#'   tibble. Either "file", or "table" must be used, but only one.
#' 
#' @return \code{\link{taxmap}}
#' 
#' @family parsers
#' 
#' @export
parse_ubiome <- function(file = NULL, table = NULL) {
  
  # Check that `file` and `text` and `table` are not used together
  are_missing <- c(file = is.null(file),
                   text = is.null(table))
  if (sum(are_missing) != 1) {
    stop(paste0('Either "file" or "table" must be supplied, but not both.'))
  }
  
  # Read raw data
  if (is.null(file)) {
    raw_data <- table
  } else {
    raw_data <- readr::read_csv(file)
  }
  
  # Make taxmap object
  output <- parse_edge_list(input = raw_data,
                            taxon_id = "taxon",
                            supertaxon_id = "parent",
                            taxon_name = "tax_name",
                            taxon_rank = "tax_rank")
  
  return(output)
}


#' Convert a table with an edge list to taxmap
#'
#' Converts a table containing an edge list into a [taxa::taxmap()] object.
#' An "edge list" is two columns in a table, where each row defines a taxon-supertaxon relationship.
#' The contents of the edge list will be used as taxon IDs.
#' The whole table will be included as a data set in the output object.
#'
#' @param input A table containing an edge list encoded by two columns.
#' @param taxon_id The name/index of the column containing the taxon IDs.
#' @param supertaxon_id The name/index of the column containing the taxon IDs for the supertaxon of the IDs in `taxon_col`.
#'
#' @family parsers
#'
#' @keywords internal
parse_edge_list <- function(input, taxon_id, supertaxon_id, taxon_name, taxon_rank = NULL) {
  
  # Create empty taxmap object
  output <- taxmap()
  
  # Make taxon ID characters
  input[taxon_id] <- as.character(input[[taxon_id]])
  input[supertaxon_id] <- as.character(input[[supertaxon_id]])
  
  # Add edge list
  output$edge_list <- data.frame(from = input[[supertaxon_id]],
                                 to = input[[taxon_id]],
                                 stringsAsFactors = FALSE)
  
  # Add taxa
  output$taxa <- lapply(seq_len(nrow(input)), function(i) {
    my_name <- input[[taxon_name]][i]
    if (is.null(taxon_rank)) {
      my_rank <- NULL
    } else {
      my_rank <- input[[taxon_rank]][i]
    }
    my_id <- input[[taxon_id]][i]
    taxon(name = my_name, rank = my_rank, id = my_id)
  })
  names(output$taxa) <- input[[taxon_id]]
  
  # Add data
  input <- dplyr::mutate(input, taxon_id = taxon_ids(output))
  input <- dplyr::select(input, taxon_id, everything())
  output$data <- list(input = input)
  
  return(output)
}


#' Convert the output of dada2 to a taxmap object
#'
#' Convert the ASV table and taxonomy table returned by dada2 into a taxmap object. An example of
#' the input format can be found by following the dada2 tutorial here:
#' shttps://benjjneb.github.io/dada2/tutorial.html
#'
#' @param seq_table The ASV abundance matrix, with rows as samples and columns as ASV ids or
#'   sequences
#' @param tax_table The table with taxonomic classifications for ASVs, with ASVs in rows and
#'   taxonomic ranks as columns.
#' @inheritParams taxa::parse_tax_data
#'
#' @return \code{\link{taxmap}}
#'
#' @family parsers
#'
#' @export
parse_dada2 <- function(seq_table, tax_table, class_key = "taxon_name", class_regex = "(.*)", include_match = TRUE) {
  # Convert sequence table to tibble
  seq_table <- t(seq_table)
  seq_table <- dplyr::bind_cols(asv_id = rownames(seq_table), dplyr::as_tibble(seq_table))
  
  # Convert taxonomy table to tibble
  tax_table <- dplyr::as_tibble(cbind(asv_id = rownames(tax_table), tax_table))
  
  # Convert to taxmap format
  output <- taxa::parse_tax_data(tax_table,
                                 class_cols = -1,
                                 named_by_rank = TRUE,
                                 class_key = class_key,
                                 class_regex = class_regex,
                                 include_match = include_match,
                                 include_tax_data = FALSE,
                                 datasets = list(asv_table = seq_table),
                                 mappings = c("asv_id" = "asv_id"))
  
  # Remove NA taxa
  output$filter_taxa(!is.na(taxon_names))
  
  return(output)
}
