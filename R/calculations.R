#' Calculate proportions from observation counts
#'
#' For a given table in a taxmap object, convert one or more columns containing
#' counts to proportions. This is meant to be used with counts associated with
#' observations (e.g. OTUs), as opposed to counts that have already been summed
#' per taxon. If counts have already been summed per taxon, use
#' [calc_tax_props()] instead.
#'
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} that have counts. By Default,
#'   all numeric columns in \code{data} are used.
#' @param keep_other_cols If \code{TRUE}, keep non-count cols in the input data.
#' The "taxon_id" column will always be preserved. 
#'
#' @return A tibble
#'
#' @family calculations
#' 
#' @export
calc_obs_props <- function(obj, dataset, cols = NULL, keep_other_cols = TRUE) {
  
  # Get count table
  count_table <- get_taxmap_table(obj, dataset, cols)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  } else { # remove any columns that do not exist
    cols <- cols[! cols %in% get_invalid_cols(count_table, cols)]
  }
  
  # Check that count columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  # Calculate proportions
  count_table[cols] <- lapply(count_table[cols], function(x) x / sum(x))
  
  # Remove other columns if specified
  if (! keep_other_cols) {
    cols_to_keep <- c(colnames(count_table[cols]), "taxon_id")
    count_table <- count_table[colnames(count_table) %in% cols_to_keep]
  }
  
  return(count_table)
}


#' Calculate proportions from observation counts
#' 
#' For a given table in a taxmap object, rarefy counts to a constant total. This
#' is a wrapper around \code{\link[vegan]{rrarefy}} that automatically detects
#' which columns are numeric and handles the reformatting needed to use tibbles.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} that have counts. By 
#'   Default, all numeric columns in \code{data} are used.
#' @param sample_size The sample size counts will be rarefied to. This can be 
#'   either a single integer or a vector of integers of equal length to the 
#'   number of columns.
#' @param keep_other_cols If \code{TRUE}, keep non-count cols in the input data.
#'   The "taxon_id" column will always be preserved.
#'   
#' @return A tibble
#' 
#' @family calculations
#'   
#' @export
rarefy_obs <- function(obj, dataset, cols = NULL, sample_size = NULL, keep_other_cols = TRUE) {
  
  # Get count table
  count_table <- get_taxmap_table(obj, dataset, cols)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  } else { # remove any columns that do not exist
    cols <- cols[! cols %in% get_invalid_cols(count_table, cols)]
  }
  
  # Check that columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  # Calculate minimum count if no sample size is given
  if (is.null(sample_size)) {
    sample_size <- min(colSums(count_table[cols]))
  }
  
  # Rarefy
  count_table[cols] <- dplyr::as.tbl(as.data.frame(t(vegan::rrarefy(t(count_table[cols]), sample = sample_size))))

  # Remove other columns if specified
  if (! keep_other_cols) {
    cols_to_keep <- c(colnames(count_table[cols]), "taxon_id")
    count_table <- count_table[colnames(count_table) %in% cols_to_keep]
  }
  
  return(count_table)
}


#' Compare treatments
#' 
#' Apply a function to compare data, usually abundance, from pairs of 
#' treatments. By default, every pairwise combination of treatments are 
#' compared. A custom function can be supplied to perform the comparison.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains data for each 
#'   sample in columns.
#' @param sample_ids The names of sample columns in \code{dataset}
#' @param treatments The treatment associated with each sample. Must be the same
#'   order and length as \code{sample_ids}.
#' @param func The function to apply for each comparison. For each row in 
#'   \code{dataset}, for each combination of treatments, this function will 
#'   recieve the data for each treatment, passed a two character vecotors.
#'   Therefor the function must take at least 2 arguments corresponding to the
#'   two treatments compared. The function should return a vector or list or
#'   results of a fixed length. If named, the names will be used in the output.
#'   The names should be consistent as well. A simple example is
#'   \code{function(x, y) mean(x) - mean(y)}. By default, the following function
#'   is used:
#'   \preformatted{
#'   function(abund_1, abund_2) {
#'     log_ratio <- log2(median(abund_1) / median(abund_2))
#'     if (is.nan(log_ratio)) {
#'       log_ratio <- 0
#'     }
#'     list(log2_median_ratio = log_ratio,
#'          median_diff = median(abund_1) - median(abund_2),
#'          mean_diff = mean(abund_1) - mean(abund_2),
#'          wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
#'   }
#'   }
#' @param combinations Which combinations of treatments to use. Must be a list 
#'   of vectors, each containing the names of 2 treatments to compare. By 
#'   default, all pairwise combinations of treatments are compared.
#' @param keep_cols If \code{TRUE}, preserve all columns not in 
#'   \code{sample_ids} in the output. If \code{FALSE}, dont keep other columns. 
#'   If a column names or indexes are supplied, only preserve those columns.
#'   
#' @return A tibble
#'   
#' @family calculations
#' 
#' @export
compare_treatments <- function(obj, dataset, sample_ids, treatments,
                               func = NULL, combinations = NULL,
                               keep_cols = TRUE) {
  # Get abundance by sample data
  abund_data <- get_taxmap_table(obj, dataset, sample_ids)
  
  # Define defualt function
  if (is.null(func)) {
    func <- function(abund_1, abund_2) {
      log_ratio <- log2(stats::median(abund_1) / stats::median(abund_2))
      if (is.nan(log_ratio)) {
        log_ratio <- 0
      }
      list(log2_median_ratio = log_ratio,
           median_diff = stats::median(abund_1) - stats::median(abund_2),
           mean_diff = mean(abund_1) - mean(abund_2),
           wilcox_p_value = stats::wilcox.test(abund_1, abund_2)$p.value)
    }
  }
  
  # Account for columns that do not exist
  invalid_cols <- get_invalid_cols(abund_data, sample_ids)
  treatments <- treatments[! sample_ids %in% invalid_cols]
  sample_ids <- sample_ids[! sample_ids %in% invalid_cols]
  
  # Parse "keep_cols" option
  if (is.logical(keep_cols)) {
    if (length(keep_cols) != 1) {
      stop('The "keep_cols" option must either be TRUE/FALSE or a vector of valid column names/indexes.', call. = FALSE)
    } else if (keep_cols) {
      keep_cols <- colnames(abund_data)[! colnames(abund_data) %in% sample_ids]
    } else {
      keep_cols <- c()
    }
  } else {
    not_used <- get_taxmap_table(obj, dataset, keep_cols) # Checks that columns specified are valid
  }
  
  # Get every combination of treatments to compare
  if (is.null(combinations)) {
    combinations <- t(utils::combn(unique(treatments), 2))
    combinations <- lapply(seq_len(nrow(combinations)), function(i) combinations[i, ])
  }
  
  # Make function to compare one pair of treatments
  one_comparison <- function(treat_1, treat_2) {
    output <- lapply(seq_len(nrow(abund_data)), function(i) {
      # get samples ids for each treatment
      samples_1 <- sample_ids[treatments == treat_1]
      samples_2 <- sample_ids[treatments == treat_2]
      
      # get abundance data for each treatment
      abund_1 <- unlist(abund_data[i, samples_1])
      abund_2 <- unlist(abund_data[i, samples_2])
      
      # Run user-supplied function on abundance info
      result <- as.list(func(abund_1, abund_2))
      
      # Complain if nothing is returned
      if (length(result) == 0) {
        stop(paste0("The function supplied returned nothing when given the following data:\n",
                    '  Treatments compared: "', treat_1, '" vs "', treat_2, '"\n',
                    '  Row index: ', i, '\n',
                    limited_print(abund_1, prefix = paste0('  "', treat_1, '" data: '), type = "silent"),
                    limited_print(abund_2, prefix = paste0('  "', treat_2, '" data: '), type = "silent")))
      }
      
      # If the output does not have names, add defaults
      if (is.null(names(result))) {
        if (length(result) == 1) {
          names(result) <- "value"
        } else {
          names(result) <- paste("value_", seq_along(result)) 
        }
      }
      
      do.call(data.frame, result)
    })
    
    # Combine the results for each row into a data frame
    output <- as.data.frame(do.call(rbind, output))
    
    # Add treatments compared
    output <- cbind(data.frame(stringsAsFactors = FALSE,
                               treatment_1 = treat_1,
                               treatment_2 = treat_2),
                    output)
    
    # Add in other columns if specified
    output <- cbind(abund_data[, keep_cols], output)
    
    return(output)
  }
  
  # Compare all pairs of treatments and combine
  output <- lapply(seq_len(length(combinations)), function(i) {
    one_comparison(combinations[[i]][1], combinations[[i]][2])
  })
  output <- do.call(rbind, output)
  
  # Convert to tibble and return
  dplyr::as.tbl(output)
}


#' Sum observation values for each taxon
#' 
#' For a given table in a taxmap object, sum the values in each column for each 
#' taxon. This is useful to convert per-observation counts (e.g. OTU counts) to 
#' per-taxon counts.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} that have counts. By 
#'   Default, all numeric columns in \code{data} are used.
#' @param col_names The names of count columns in the output. Must be the same
#'   length as \code{cols}.
#'   
#' @return A tibble
#'   
#' @family calculations
#'   
#' @export
calc_taxon_abund <- function(obj, dataset, cols = NULL, col_names = NULL) {
  
  # Get count table
  count_table <- get_taxmap_table(obj, dataset, cols)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  }  else { # remove any columns that do not exist
    cols <- cols[! cols %in% get_invalid_cols(count_table, cols)]
  }
  
  # Get output column names
  if (is.null(col_names)) {
    col_names <- names(count_table[cols])
  } else if (length(col_names) != length(cols)) { 
    stop("The length of 'cols' and 'col_names' are not equal")
  }  

  # Check that count columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  # Sum counts per taxon for each sample
  obs_indexes <- obj$obs(dataset)
  output <- lapply(cols, function(col_index) {
    vapply(obs_indexes, function(i) sum(count_table[[col_index]][i]), numeric(1))
  })
  output <- as.data.frame(output, stringsAsFactors = FALSE)
  colnames(output) <- colnames(count_table[cols])
  
  # Add taxon_id column
  output <- cbind(data.frame(taxon_id = obj$taxon_ids(), stringsAsFactors = FALSE),
                  output)
  
  # Rename cols
  colnames(output)[match(cols, colnames(output))] <- col_names
  
  # Convert to tibble and return
  dplyr::as.tbl(output)
}

