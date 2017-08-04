#' Calculate proportions from observation counts
#'
#' For a given table in a taxmap object, convert one or more columns containing
#' counts to proportions. This is meant to be used with counts associated with
#' observations (e.g. OTUs), as opposed to counts that have already been summed
#' per taxon. If counts have already been summed per taxon, use
#' [calc_tax_props()] instead.
#'
#' @param obj A taxmap object
#' @param data The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} that have counts. By Default,
#'   all numeric columns in \code{data} are used.
#' @param keep_other_cols If \code{TRUE}, keep non-count cols in the input data.
#' The "taxon_id" column will always be preserved. 
#'
#' @return A tibble
#'
#' @export
calc_obs_props <- function(obj, dataset, cols = NULL, keep_other_cols = TRUE) {
  # Check that dataset exists and is a table
  if (! dataset %in% names(obj$data)) {
    stop(paste0('The dataset "', dataset,
                '" is not in the object supplied. Datasets found include:\n  ',
                limited_print(names(obj$data), type = "silent")))
  }
  if (! is.data.frame(obj$data[[dataset]])) {
    stop(paste0('The dataset "', dataset,  '" is not a table.'))
  }

  # Get count table
  count_table <- obj$data[[dataset]]

  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  }

  # Check that count columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }

  # Check that counts are not summed per taxon
  # NOT DONE YET
  # if ("taxon_id" %in% colnames(count_table)) {
  #   lapply(obj$subtaxa(value = "taxon_ids", recursive = FALSE), 
  #          function(ids) {
  #            test_col <- count_table[[cols[1]]]
  #            sum(as.numeric(test_col[match(ids, count_table$taxon_id)]))
  #          })
  # }

  # Calculate proportions
  prop_table <- count_table
  prop_table[cols] <- lapply(prop_table[cols], function(x) x / sum(x))
  
  # Remove other columns if specified
  if (! keep_other_cols) {
    cols_to_keep <- c(colnames(prop_table[cols]), "taxon_id")
    prop_table <- prop_table[colnames(prop_table) %in% cols_to_keep]
  }
  
  return(prop_table)
}