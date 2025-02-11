#' Convert taxmap to phyloseq
#' 
#' Convert a taxmap object to a phyloseq object.
#' 
#' @param obj The taxmap object.
#' @param otu_table The table in `obj$data` with OTU counts. Must be one of the following:
#' \describe{
#'   \item{\code{NULL}}{Look for a table named "otu_table" in `obj$data` with taxon IDs, OTU IDs, and OTU counts. If it exists, use it.}
#'   \item{\code{character}}{The name of the table stored in `obj$data` with taxon IDs, OTU IDs, and OTU counts}
#'   \item{\code{data.frame}}{A table with taxon IDs, OTU IDs, and OTU counts}
#'   \item{\code{FALSE}}{Do not include an OTU table, even if "otu_table" exists in `obj$data`}
#' }
#' @param otu_id_col The name of the column storing OTU IDs in the OTU table.
#' @param sample_data A table containing sample data with sample IDs matching
#'   column names in the OTU table. Must be one of the following:
#' \describe{
#'   \item{\code{NULL}}{Look for a table named "sample_data" in `obj$data`. If it exists, use it.}
#'   \item{\code{character}}{The name of the table stored in `obj$data` with sample IDs}
#'   \item{\code{data.frame}}{A table with sample IDs}
#'   \item{\code{FALSE}}{Do not include a sample data table, even if "sample_data" exists in `obj$data`}
#' }
#' @param sample_id_col The name of the column storing sample IDs in the sample data table.
#' @param phy_tree A phylogenetic tree of class \code{ape:phylo} from
#'   the \code{ape} package with tip labels matching OTU ids. Must be one of the following:
#' \describe{
#'   \item{\code{NULL}}{Look for a tree named "phy_tree" in `obj$data` with tip labels matching OTU ids. If it exists, use it.}
#'   \item{\code{character}}{The name of the tree stored in `obj$data` with tip labels matching OTU ids.}
#'   \item{\code{ape::phylo}}{A tree with tip labels matching OTU ids.}
#'   \item{\code{FALSE}}{Do not include a tree, even if "phy_tree" exists in `obj$data`}
#' }
#' 
#' @examples 
#' \donttest{
#' # Parse example dataset
#' library(phyloseq)
#' data(GlobalPatterns)
#' x <- parse_phyloseq(GlobalPatterns)
#' 
#' # Convert back to a phylseq object
#' as_phyloseq(x)
#' }
#' @export
as_phyloseq <- function(obj,
                        otu_table = NULL, otu_id_col = "otu_id",
                        sample_data = NULL, sample_id_col = "sample_id",
                        phy_tree = NULL) {
  
  # Check that phyloseq is intalled
  if (! requireNamespace("phyloseq", quietly = TRUE)) {
    stop('The "phyloseq" package needs to be installed for this function to work.',
         call. = FALSE)
  }
  
  # Get and check OTU table
  otu_table <- get_expected_data(obj, input = otu_table, default = "otu_table",
                                 expected_class = "data.frame")
  if (! is.null(otu_table)) {
    # Get OTU IDs
    if (otu_id_col %in% colnames(otu_table)) {
      otu_ids <- otu_table[[otu_id_col]]
    } else  {
      stop('OTU table does not have an OTU ID column named "', otu_id_col,
           '". Use the "otu_id_col" option if it is named something else.')
    }
    
    # Get taxon IDs
    if ("taxon_id" %in% colnames(otu_table)) {
      otu_taxon_ids <- otu_table$taxon_id
    } else {
      stop('OTU table is not named by taxon IDs.')
    }
    
    # Get OTU abundance matrix
    is_num_col <- vapply(otu_table, is.numeric, logical(1))
    is_invalid_col <- ! is_num_col & ! colnames(otu_table) %in% c(otu_id_col, 'taxon_id')
    if (any(is_invalid_col)) {
      warning(call. = FALSE,
              'Discarding non-numeric columns in OTU table:', 
              limited_print(colnames(otu_table)[is_invalid_col],
                            type = "silent", prefix = "  "))
    }
    parsed_otu_table <- as.matrix(otu_table[, is_num_col])
    rownames(parsed_otu_table) <- otu_ids
  }
  
  # Get taxonomy table
  if (is.null(otu_table)) {
    tax_table <- as.matrix(obj$taxonomy_table())
    rownames(tax_table) <- obj$taxon_ids()
  } else {
    tax_table <- as.matrix(obj$taxonomy_table(subset = otu_taxon_ids))
    rownames(tax_table) <- otu_ids
  }
  ps_tax_table <- phyloseq::tax_table(tax_table)
  rownames(ps_tax_table) <- rownames(tax_table)
  
  # Get sample data
  sample_table <- get_expected_data(obj, input = sample_data,
                                    default = "sample_data",
                                    expected_class = "data.frame")
  if (! is.null(sample_table)) {
    # Get sample IDs
    if (sample_id_col %in% colnames(sample_table)) {
      sample_ids <- sample_table[[sample_id_col]]
    } else {
      stop('Sample data table does not have an sample ID column named "', sample_id_col,
           '". Use the "sample_id_col" option if it is named something else.')
    }
    
    # Check for sample information not in the sample table
    if (! is.null(otu_table)) {
      unnkown_cols <- colnames(parsed_otu_table)[! colnames(parsed_otu_table) %in% sample_ids]
      if (length(unnkown_cols) > 0) {
        warning('The OTU table contains the following ', length(unnkown_cols),
                ' of ', ncol(parsed_otu_table), ' samples that do not appear in the sample data table:\n',
                limited_print(prefix = "  ", type = "silent", unnkown_cols))
      }
    }
    
    # reformt
    ps_sample_table <- as.data.frame(sample_table)
    rownames(ps_sample_table) <- ps_sample_table[[sample_id_col]]
    ps_sample_table <- ps_sample_table[colnames(ps_sample_table) != sample_id_col]
    ps_sample_table <- phyloseq::sample_data(ps_sample_table)
  } else {
    ps_sample_table <- NULL
  }
  
  # Get phylogenetic tree
  my_phy_tree <- get_expected_data(obj, input = phy_tree, default = "phy_tree",
                                   expected_class = "phylo")
  
  
  # Make phyloseq object
  phyloseq::phyloseq(phyloseq::otu_table(parsed_otu_table, taxa_are_rows = TRUE),
                     ps_tax_table,
                     ps_sample_table,
                     my_phy_tree)
}


#' Get a data set in as_phyloseq
#' 
#' Get a data set in as_phyloseq
#' 
#' @param obj The taxmap object
#' @param input The input to as_phyloseq options.
#' @param default The default name of the data set.
#' @param expected_class What the dataset is expected to be.
#' 
#' @keywords internal
get_expected_data <- function(obj, input, default, expected_class) {
  
  get_dataset_by_name <- function(data_name) {
    output <- obj$data[[data_name]]
    if (any(expected_class %in% class(output))) {
      return(output)
    } else { # named dataset found, but wrong class
      stop(call. = FALSE,
           'Data set named "', data_name, '" found by it is not one of the accepted classes:',
           limited_print(expected_class, type = "silent", prefix = "  "))
    }
  }
  
  # If FALSE, return nothing
  if (is.logical(input) && input == FALSE) {
    return(NULL)
  }
  
  # If NULL, look for expected table. return nothing if not found
  if (is.null(input)) {
    if (default %in% names(obj$data)) {
      return(get_dataset_by_name(default))
    } else {
      return(NULL)
    }
  }
  
  # If a character vector, then look for named dataset
  if (is.character(input)) {
    return(get_dataset_by_name(input))
  } 
  
  # Otherwise make sure input is right class
  if (any(expected_class %in% class(input))) {
    return(input)
  } else {
    stop(call. = FALSE,
         'Data set given is not one of the accepted classes:',
         limited_print(expected_class, type = "silent", prefix = "  "))
  }
}
