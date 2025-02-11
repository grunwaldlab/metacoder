#' Differential abundance with DESeq2
#'
#' EXPERIMENTAL: This function is still being tested and developed; use with caution. Uses the
#' \code{\link[DESeq2]{DESeq2-package}} package to conduct differential abundance analysis of count data. Counts can
#' be of OTUs/ASVs or taxa. The plotting function \code{\link{heat_tree_matrix}} is useful for
#' visualizing these results. See details section below for considerations on preparing data for
#' this analysis.
#'
#' Data should be raw read counts, not rarefied, converted to proportions, or modified with any
#' other technique designed to correct for sample size since \code{\link[DESeq2]{DESeq2-package}} is designed to be
#' used with count data and takes into account unequal sample size when determining differential
#' abundance. Warnings will be given if the data is not integers or all sample sizes are equal.
#'
#' @param obj A \code{\link[metacoder]{taxmap}} object
#' @param data The name of a table in \code{obj} that contains data for each sample in columns.
#' @param cols The names/indexes of columns in \code{data} to use. By default, all numeric columns
#'   are used. Takes one of the following inputs: \describe{ \item{TRUE/FALSE:}{All/No columns will
#'   used.} \item{Character vector:}{The names of columns to use} \item{Numeric vector:}{The indexes
#'   of columns to use} \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Use
#'   the columns corresponding to \code{TRUE} values.} }
#' @param groups A vector defining how samples are grouped into "treatments". Must be the same order
#'   and length as \code{cols}.
#' @param other_cols If \code{TRUE}, preserve all columns not in \code{cols} in the output. If
#'   \code{FALSE}, dont keep other columns. If a column names or indexes are supplied, only preserve
#'   those columns.
#' @param lfc_shrinkage What technique to use to adjust the log fold change results for low counts.
#'   Useful for ranking and visualizing log fold changes. Must be one of the following:
#' \describe{
#'   \item{'none'}{No log fold change adjustments.}
#'   \item{'normal'}{The original DESeq2 shrinkage estimator}
#'   \item{'ashr'}{Adaptive shrinkage estimator from the \code{ashr} package, using a fitted mixture of normals prior.}
#' }
#' @param ... Passed to \code{\link[DESeq2]{results}} if the \code{lfc_shrinkage} option is "none"
#'   and to \code{\link[DESeq2]{lfcShrink}} otherwise.
#' 
#' @return A tibble with at least the taxon ID of the thing tested, the groups compared, and the
#'   DESeq2 results. The \code{log2FoldChange} values will be positive if \code{treatment_1} is more
#'   abundant and \code{treatment_2}.
#'
#' @family calculations
#'
#' @examples
#' \donttest{
#' # Parse data for plotting
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#'
#' # Get per-taxon counts
#' x$data$tax_table <- calc_taxon_abund(x, data = "tax_data", cols = hmp_samples$sample_id)
#'
#' # Calculate difference between groups
#' x$data$diff_table <- calc_diff_abund_deseq2(x, data = "tax_table",
#'                                     cols = hmp_samples$sample_id,
#'                                     groups = hmp_samples$body_site)
#'                                     
#' # Plot results (might take a few minutes)
#' heat_tree_matrix(x,
#'                  data = "diff_table",
#'                  node_size = n_obs,
#'                  node_label = taxon_names,
#'                  node_color = ifelse(is.na(padj) | padj > 0.05, 0, log2FoldChange),
#'                  node_color_range = diverging_palette(),
#'                  node_color_trans = "linear",
#'                  node_color_interval = c(-3, 3),
#'                  edge_color_interval = c(-3, 3),
#'                  node_size_axis_label = "Number of OTUs",
#'                  node_color_axis_label = "Log2 fold change")
#'
#' }
#'
#' @export
calc_diff_abund_deseq2 <- function(obj, data, cols, groups, other_cols = FALSE, 
                                   lfc_shrinkage = c('none', 'normal', 'ashr'), ...) {
  
  # Check that DESeq2 is intalled
  if (! requireNamespace("DESeq2", quietly = TRUE)) {
    stop('The "DESeq2" package needs to be installed for this function to work.',
         call. = FALSE)
  }
  
  # Parse options
  lfc_shrinkage <- match.arg(lfc_shrinkage)

  # Get abundance by sample data
  abund_data <- get_taxmap_table(obj, data)
  
  # Parse columns to use
  cols <- get_numeric_cols(obj, data, cols)
  
  # Check that columns contain integers only
  col_is_int <- vapply(cols, FUN.VALUE = logical(1), function(c) {
    all(abund_data[[c]] %% 1 == 0)
  })
  non_int_cols <- cols[! col_is_int]
  if (length(non_int_cols) > 0) {
    stop(call. = FALSE,
         'All data given to DESeq2 should be untransformed counts. The following columns contain non-integers:\n',
         limited_print(type = 'silent', non_int_cols, prefix = '  '))
  }
  
  # Check for equal sample sizes
  if (length(unique(colSums(abund_data[cols]))) == 1) {
    warning(call. = FALSE,
            'All columns have equal counts, suggesting counts were normalized (e.g. rarefied) to correct for sample size variation.',
            ' DESeq2 is designed to be used with unnormalized data.',
            ' Use untransformed counts if available.')
  }
  
  # Check groups option
  groups <- check_option_groups(groups, cols)
  
  # Find other columns
  #   These might be added back to the output later
  other_cols <- get_taxmap_other_cols(obj, data, cols, other_cols)
  
  # Get every combination of groups to compare
  combinations <- t(utils::combn(unique(as.character(groups)), 2))
  combinations <- lapply(seq_len(nrow(combinations)), function(i) combinations[i, ])
  
  # Format data for DESeq2
  metadata <- data.frame(group = groups)
  deseq_data <- DESeq2::DESeqDataSetFromMatrix(countData = abund_data[cols],
                                               colData = metadata,
                                               design = ~ group)
  
  # Run DESeq2
  deseq_result <- DESeq2::DESeq(deseq_data)
  
  # Make function to compare one pair of groups
  one_comparison <- function(treat_1, treat_2) {
    
    # Extract results for a single pair
    if (lfc_shrinkage == 'none') {
      pair_result <- DESeq2::results(deseq_result, contrast = c('group', treat_1, treat_2), ...)
    } else {
      pair_result <- DESeq2::lfcShrink(deseq_result, contrast = c('group', treat_1, treat_2),
                                       type = lfc_shrinkage, ...)
    }

    # Add treatments compared
    output <- cbind(data.frame(stringsAsFactors = FALSE,
                               treatment_1 = treat_1,
                               treatment_2 = treat_2),
                    as.data.frame(pair_result))
    
    # Add in other columns if specified
    output <- cbind(abund_data[, other_cols], output)
    
    return(output)
  }
  
  # Compare all pairs of treatments and combine
  output <- lapply(seq_len(length(combinations)), function(i) {
    one_comparison(combinations[[i]][1], combinations[[i]][2])
  })
  output <- do.call(rbind, output)
  
  # Convert to tibble and return
  dplyr::as_tibble(output)
}
