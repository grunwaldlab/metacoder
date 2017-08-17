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
#' @export
compare_treatments <- function(obj, dataset, sample_ids, treatments,
                               func = NULL, combinations = NULL,
                               keep_cols = TRUE) {
  # Get abundance by sample data
  abund_data <- obj$data[[dataset]]
  
  # Define defualt function
  if (is.null(func)) {
    func <- function(abund_1, abund_2) {
      log_ratio <- log2(median(abund_1) / median(abund_2))
      if (is.nan(log_ratio)) {
        log_ratio <- 0
      }
      list(log2_median_ratio = log_ratio,
           median_diff = median(abund_1) - median(abund_2),
           mean_diff = mean(abund_1) - mean(abund_2),
           wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
    }
  }
  
  # Parse "keep_cols" option
  kc_error_msg <- 'The "keep_cols" option must either be TRUE/FALSE or a vector of valid column names/indexes.'
  if (is.logical(keep_cols)) {
    if (length(keep_cols) != 1) {
      stop(kc_error_msg)
    } else if (keep_cols) {
      keep_cols <- colnames(abund_data)[! colnames(abund_data) %in% sample_ids]
    } else {
      keep_cols <- c()
    }
  } else {
    if (is.numeric(keep_cols)) {
      invalid_cols <- keep_cols[! keep_cols %in% seq_len(ncol(abund_data))]
    } else {
      invalid_cols <- keep_cols[! keep_cols %in% colnames(abund_data)]
    }
    if (length(invalid_cols) > 0) {
      stop(paste0(kc_error_msg, 
                  " The following column names/indexes are not valid:\n",
                  "  ", limited_print(invalid_cols, type = "silent")))
    }
  }
  
  # Get every combination of treatments to compare
  if (is.null(combinations)) {
    combinations <- t(combn(unique(treatments), 2))
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
#'   
#' @return A tibble
#'   
#' @export
calc_taxon_abund <-function(obj, dataset, cols = NULL) {
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
  
  # Sum counts per taxon for each sample
  output <- lapply(cols, function(col_index) {
    count_table[col_index]
    taxa::obs_apply(obj, dataset, function(i) sum(count_table[i, col_index]), simplify = TRUE)
  })
  output <- as.data.frame(output, stringsAsFactors = FALSE)
  colnames(output) <- colnames(count_table[cols])
  
  # Add taxon_id column
  output <- cbind(data.frame(taxon_id = obj$taxon_ids(), stringsAsFactors = FALSE),
                  output)
  
  # Convert to tibble and return
  dplyr::as.tbl(output)
}

